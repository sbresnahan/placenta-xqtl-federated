#!/usr/bin/env python3
"""
xqtl_gxe_core.py -- torch core for SNP x exposure (GxE) cis-interaction scans.


Statistical model (Objective 2.1):
    Y = b_g * G + b_e * E + b_gi * (G x E) + Z_cov * gamma + eps

The regression math mirrors tensorQTL 1.0.10
(tensorqtl.core.calculate_interaction_nominal) exactly:
  * g, e, g:e and the phenotype are centered, then residualized on the
    covariates with the tensorQTL Residualizer (QR of the centered covariate
    matrix; residualizer.dof = n - 2 - k),
  * per variant the design is X = [g0, e0, ge0]  (n_samples x 3),
  * b = (X'X)^-1 X' p0,  rss = ||X b - p0||^2,
  * b_se = sqrt(diag((X'X)^-1) * rss / dof),
  * dof = residualizer.dof - 2 * n_interactions = n - k - 4  (one exposure).

rss is evaluated in closed form, rss = p0'p0 - b'(X'X)b = p0'p0 - b'w with
w = X'p0 (algebraically identical to the explicit residual sum of squares;
float32 relative error ~1e-6, verified against tensorQTL
calculate_interaction_nominal). The closed form is what makes the permutation scan cheap.

Permutation scheme (tier 1): PHENOTYPE permutation, the GTEx/tensorQTL
convention. Each permutation reorders the RAW phenotype (breaking its
alignment to genotypes AND covariates — that is the null), and the permuted
phenotype is re-centered and re-residualized against the fixed covariate
projection Q Q', exactly as tensorQTL's map_cis does inside
calculate_corr. Because genotypes and the exposure are never permuted, the
design X and (X'X)^-1 stay fixed; each permutation costs one [block x k]
covariate projection plus one batched matvec for the interaction fit.

  Caveat (documented in docs/INTERACTIONS.md): phenotype permutation is slightly
  conservative for the interaction null when a main G effect exists, because
  the permuted null carries neither main nor interaction signal while
  rss_perm is correspondingly larger than the observed rss. This is the same
  permutation scheme tensorQTL's own interaction mode uses.

Multiple testing (tier 1): per phenotype, r2 = t^2 / (t^2 + dof) of the top
cis variant (max over the window); the permutation max-statistic gives the
empirical pval_perm and the Beta-approximation pval_beta via
tensorqtl.core.calculate_beta_approx_pval — identical to cis.map_cis.
Storey q-values are added at merge time via bin/compute_qvalues.R
(GTEx convention; Benjamini-Hochberg fallback recorded in qvalue_method).

Adaptive permutation schedule: cumulative blocks (100, 500, 1000, 10000) by default;
a phenotype stops early when phat = (exceedances + 1) / (nperm + 1) > 0.10
(clearly null phenotypes do not waste the full 10k permutations; their
pval_beta is fit on fewer permutations and is correspondingly crude, which
does not matter — they are nowhere near the FDR threshold).

Also home to ivw_meta (inverse-variance-weighted meta-analysis + Q/I^2), used by the central
workflow for the tier-2 synthesis.
"""

import numpy as np
import torch
from scipy import stats


# ---------------------------------------------------------------------------
# Residualizer (mirrors tensorqtl.core.Residualizer)
# ---------------------------------------------------------------------------

def make_residualizer(covariates, device=None):
    """Build the tensorQTL Residualizer from a covariate matrix.

    Parameters
    ----------
    covariates : array-like [n_samples x k]
    device     : torch device (default: cuda if available)

    Returns
    -------
    Q_t  : torch.float32 [n_samples x k] — Q of the QR of centered covariates
    dof  : int — n_samples - 2 - k (tensorQTL convention)
    """
    if device is None:
        device = get_device()
    C_t = torch.as_tensor(np.asarray(covariates), dtype=torch.float32).to(device)
    Q_t, _ = torch.linalg.qr(C_t - C_t.mean(0))
    dof = C_t.shape[0] - 2 - C_t.shape[1]
    return Q_t, dof


def residualizer_transform(M_t, Q_t, center=False):
    """Residualize rows of M_t wrt the covariate space spanned by Q.

    Faithful port of tensorqtl.core.Residualizer.transform: the projection
    is computed from the centered rows; with center=False the row means are
    added back (the inputs of the interaction path are pre-centered, so the
    two are numerically identical there).
    """
    M0_t = M_t - M_t.mean(1, keepdim=True)
    proj = torch.mm(torch.mm(M0_t, Q_t), Q_t.t())
    if center:
        return M0_t - proj
    return M_t - proj


def get_device(prefer_gpu=True):
    if prefer_gpu and torch.cuda.is_available():
        return torch.device("cuda")
    return torch.device("cpu")


# ---------------------------------------------------------------------------
# Interaction scan
# ---------------------------------------------------------------------------

def prepare_interaction_terms(genotypes_t, exposure_t, Q_t):
    """Center + residualize g, e and g:e (mirrors calculate_interaction_nominal).

    Parameters
    ----------
    genotypes_t : torch.float32 [n_variants x n_samples] — hardcalls/dosages,
                  already mean-imputed (no missing values).
    exposure_t  : torch.float32 [n_samples] — raw exposure (NOT pre-centered;
                  centering happens here, as in tensorQTL).
    Q_t         : torch.float32 [n_samples x k] or None.

    Returns
    -------
    g0_t  : [n_variants x n_samples]
    e0_t  : [n_samples]
    ge0_t : [n_variants x n_samples]
    """
    g0_t = genotypes_t - genotypes_t.mean(1, keepdim=True)
    ge_t = genotypes_t * exposure_t.unsqueeze(0)          # raw product
    ge0_t = ge_t - ge_t.mean(1, keepdim=True)
    e0_t = exposure_t - exposure_t.mean(0)                # [n_samples]
    if Q_t is not None:
        g0_t = residualizer_transform(g0_t, Q_t, center=False)
        e0_t = residualizer_transform(e0_t.unsqueeze(0), Q_t, center=False).squeeze(0)
        ge0_t = residualizer_transform(ge0_t, Q_t, center=False)
    return g0_t, e0_t, ge0_t


def prepare_phenotype(phenotype_t, Q_t):
    """Center + residualize one phenotype vector [n_samples] -> [n_samples]."""
    p0_t = phenotype_t - phenotype_t.mean(0)
    if Q_t is not None:
        p0_t = residualizer_transform(p0_t.unsqueeze(0), Q_t, center=False).squeeze(0)
    return p0_t


def build_window_design(g0_w, e0_t, ge0_w, min_ge_var=1e-10):
    """Build X, (X'X)^-1 and X' for one cis window.

    Parameters
    ----------
    g0_w, ge0_w : torch.float32 [n_window x n_samples]
    e0_t        : torch.float32 [n_samples]

    Returns
    -------
    Xt_w   : [n_window x 3 x n_samples]
    Xinv_w : [n_window x 3 x 3]
    ok_w   : bool tensor [n_window] — False where ge0 is ~constant (X'X
             singular); those variants get NaN stats downstream.
    """
    nw = g0_w.shape[0]
    X_t = torch.stack([g0_w, e0_t.expand(nw, -1), ge0_w], dim=2)  # nw x ns x 3
    Xt_w = X_t.transpose(1, 2).contiguous()                       # nw x 3 x ns
    XtX = torch.matmul(Xt_w, X_t)                                 # nw x 3 x 3
    ge_var = ge0_w.var(1, unbiased=False)
    ok_w = ge_var > min_ge_var
    XtX_safe = XtX.clone()
    # Replace singular slices with identity so the batch inverse cannot fail;
    # the ok_w mask zeroes those variants out of every downstream statistic.
    if not bool(ok_w.all()):
        eye = torch.eye(3, dtype=XtX.dtype, device=XtX.device).expand_as(XtX)
        XtX_safe = torch.where(ok_w.unsqueeze(-1).unsqueeze(-1), XtX, eye)
    try:
        Xinv_w = torch.linalg.inv(XtX_safe)
    except Exception:
        # Extremely rare numerical failure: fall back to per-variant pinv.
        Xinv_w = torch.stack([torch.linalg.pinv(XtX_safe[i]) for i in range(nw)])
    Xinv_w = torch.where(ok_w.unsqueeze(-1).unsqueeze(-1), Xinv_w,
                         torch.full_like(Xinv_w, float('nan')))
    return Xt_w, Xinv_w, ok_w


def window_fit(Xt_w, Xinv_w, p0_t, dof):
    """Fit b, b_se, t for one (residualized) phenotype across a cis window.

    rss in closed form: rss = p0'p0 - b'w with w = X'p0 (see module docstring).

    Returns
    -------
    b     : [n_window x 3]  (columns: g, e, g:e)
    b_se  : [n_window x 3]
    tstat : [n_window x 3]  (float32; NaN where the variant was masked)
    """
    w = torch.matmul(Xt_w, p0_t)                                  # nw x 3
    b = torch.matmul(Xinv_w, w.unsqueeze(-1)).squeeze(-1)         # nw x 3
    rss = torch.clamp(p0_t.dot(p0_t) - (b * w).sum(1), min=1e-12)  # nw
    b_se = torch.sqrt(Xinv_w.diagonal(dim1=1, dim2=2) * rss.unsqueeze(1) / dof)
    tstat = (b.double() / b_se.double()).float()
    return b, b_se, tstat


def t_to_r2(t, dof):
    """Squared correlation scale of a t statistic: r2 = t^2 / (t^2 + dof)."""
    t2 = np.asarray(t, dtype=np.float64) ** 2
    return t2 / (t2 + dof)


def r2_to_pval(r2, dof):
    """Two-sided Student-t p-value from r2 (tensorqtl.core.pval_from_corr)."""
    r2 = np.clip(np.asarray(r2, dtype=np.float64), 0.0, 1.0 - 1e-16)
    tstat2 = dof * r2 / (1.0 - r2)
    return 2.0 * stats.t.cdf(-np.sqrt(tstat2), dof)


def adaptive_permutation_scan(Xt_w, Xinv_w, p_raw_t, Q_t, dof, r2_nominal_max,
                              rng, blocks=(100, 400, 500, 9000), stop_p=0.10):
    """Phenotype-permutation max-statistic scan for one phenotype.

    Each permutation reorders the raw phenotype; the permuted phenotype is
    re-centered and re-residualized against the fixed covariate projection
    (mirroring tensorQTL's map_cis permutation path). X and (X'X)^-1 are
    fixed, so each permutation costs one [block x k] projection plus one
    batched matvec. Permutation indices are drawn per phenotype from `rng`
    (seeded per phenotype by the caller — deterministic and chunk-position
    invariant).

    Parameters
    ----------
    Xt_w, Xinv_w : window design from build_window_design
    p_raw_t      : RAW phenotype [n_samples] (not centered/residualized)
    Q_t          : covariate projection Q [n_samples x k] (or None)
    dof          : interaction dof (n - k - 4)
    r2_nominal_max : max nominal interaction r2 over the window (scalar)
    rng          : np.random.Generator (per-phenotype)
    blocks       : permutation block sizes (100, 400, 500, 9000 ->
                   100/500/1000/10000 cumulative)
    stop_p       : early-stop threshold on the running empirical p-value

    Returns
    -------
    r2_perm   : np.ndarray [nperm] — per-permutation max interaction r2
    nperm     : int — permutations actually run
    pval_perm : float — (exceedances + 1) / (nperm + 1)
    """
    ns = p_raw_t.shape[0]
    c_t = p_raw_t - p_raw_t.mean(0)          # centered once (mean is perm-invariant)
    diag_gi = Xinv_w[:, 2, 2].unsqueeze(1)                      # nw x 1
    r2_perm_parts = []
    n_done = 0
    exceed = 0
    for block in blocks:
        perm_ix = np.stack([rng.permutation(ns) for _ in range(block)])
        ix_t = torch.from_numpy(perm_ix).to(p_raw_t.device)
        P_t = c_t[ix_t]                                         # block x ns
        if Q_t is not None:
            PQ = torch.mm(P_t, Q_t)                             # block x k
            P_t = P_t - torch.mm(PQ, Q_t.t())                   # re-residualize
        ptp = (P_t * P_t).sum(1)                                # block
        W = torch.einsum('vks,ps->vkp', Xt_w, P_t)              # nw x 3 x block
        B = torch.einsum('vkj,vjp->vkp', Xinv_w, W)             # nw x 3 x block
        rss = torch.clamp(ptp.unsqueeze(0) - (B * W).sum(1),
                          min=1e-12)                            # nw x block
        t2 = (B[:, 2, :] ** 2) / (diag_gi * rss / dof)          # nw x block
        r2 = t2 / (t2 + dof)
        r2 = torch.where(torch.isfinite(r2), r2, torch.full_like(r2, -1.0))
        r2_perm_parts.append(r2.max(dim=0).values.cpu().numpy())
        n_done += block
        r2_perm = np.concatenate(r2_perm_parts)
        exceed = int((r2_perm >= r2_nominal_max).sum())
        if (exceed + 1) / (n_done + 1) > stop_p:
            break
    r2_perm = np.concatenate(r2_perm_parts)
    pval_perm = (exceed + 1) / (n_done + 1)
    return r2_perm, n_done, pval_perm


def beta_approx_pval(r2_perm, r2_nominal, dof):
    """Beta-approximation tail p-value (tensorqtl.core.calculate_beta_approx_pval).

    Returns dict with pval_beta, beta_shape1, beta_shape2, true_df,
    pval_true_df; all NaN when the fit fails (degenerate permutation
    distributions — e.g. too few permutations after early stop — fall back
    to NaN and the merge step keeps pval_perm for those).
    """
    from tensorqtl.core import calculate_beta_approx_pval
    out = dict(pval_beta=np.nan, beta_shape1=np.nan, beta_shape2=np.nan,
               true_df=np.nan, pval_true_df=np.nan)
    try:
        with np.errstate(all='ignore'):
            pval_beta, b1, b2, true_dof, pval_true = calculate_beta_approx_pval(
                np.asarray(r2_perm, dtype=np.float64), float(r2_nominal), int(dof))
        if np.isfinite(pval_beta):
            out.update(pval_beta=float(pval_beta), beta_shape1=float(b1),
                       beta_shape2=float(b2), true_df=float(true_dof),
                       pval_true_df=float(pval_true))
    except Exception:
        pass
    return out


def allele_stats(genotypes_t):
    """af, ma_samples, ma_count per variant (tensorqtl.core.get_allele_stats)."""
    n2 = 2 * genotypes_t.shape[1]
    af_t = genotypes_t.sum(1) / n2
    ix_t = af_t <= 0.5
    m = genotypes_t > 0.5
    a = m.sum(1).int()
    b = (genotypes_t < 1.5).sum(1).int()
    ma_samples_t = torch.where(ix_t, a, b)
    a = (genotypes_t * m.float()).sum(1).int()
    ma_count_t = torch.where(ix_t, a, n2 - a)
    return af_t, ma_samples_t, ma_count_t


# ---------------------------------------------------------------------------
# Inverse-variance-weighted meta-analysis (tier 2)
# ---------------------------------------------------------------------------

def ivw_meta(betas, ses):
    """Fixed-effect IVW meta-analysis of per-ancestry interaction betas.

    Returns dict(beta, se, z, p, q_het, p_het, i2, k). NaN inputs are
    dropped; k < 2 returns NaNs for the heterogeneity stats; k == 0 returns
    all NaN.
    """
    betas = np.asarray(betas, dtype=np.float64)
    ses = np.asarray(ses, dtype=np.float64)
    ok = np.isfinite(betas) & np.isfinite(ses) & (ses > 0)
    betas, ses = betas[ok], ses[ok]
    k = len(betas)
    if k == 0:
        return dict(beta=np.nan, se=np.nan, z=np.nan, p=np.nan,
                    q_het=np.nan, p_het=np.nan, i2=np.nan, k=0)
    w = 1.0 / ses ** 2
    beta = float((w * betas).sum() / w.sum())
    se = float(np.sqrt(1.0 / w.sum()))
    z = beta / se
    p = float(2.0 * stats.norm.sf(abs(z)))
    if k >= 2:
        q_het = float((w * (betas - beta) ** 2).sum())
        p_het = float(stats.chi2.sf(q_het, k - 1))
        i2 = float(max(0.0, (q_het - (k - 1)) / q_het)) if q_het > 0 else 0.0
    else:
        q_het, p_het, i2 = np.nan, np.nan, np.nan
    return dict(beta=beta, se=se, z=z, p=p, q_het=q_het, p_het=p_het, i2=i2, k=k)
