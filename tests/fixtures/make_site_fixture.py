#!/usr/bin/env python3
"""Helper for make_site_fixture.sh (see that script for what is real and what is synthetic).

  select  --ped PED --tmp DIR          choose cohort and panel samples from 1kGP
  build   --tmp DIR --out DIR          write the TOPMed-style VCF, metadata, samplesheet,
                                       synthetic phenotypes (planted cis effects) and TRUTH.tsv
"""
import argparse
import gzip
import os

import numpy as np
import pandas as pd

SEED = 20261029
COHORT = {"EUR": 60, "EAS": 60, "AFR": 10}
PANEL_PER_GROUP = 40
# AMR is left out of the fixture panel: on ~1000 chr22 SNPs the admixed AMR samples overlap the
# EUR cluster completely, so no European cohort sample would pass the posterior threshold. The
# real panel (genome-wide, all 3202 1kGP samples) keeps all five superpopulations.
PANEL_GROUPS = ["AFR", "EAS", "EUR", "SAS"]
CELL_TYPES = ["Syncytiotrophoblast", "Cytotrophoblast", "Extravillous_trophoblast", "Stromal",
              "Endothelial", "Hofbauer", "Maternal"]
N_GENES = 40
GENE_FROM, GENE_TO = 19_300_000, 21_200_000
# planted cis effects: (modality, gene index, effect per ALT allele on the log2 / logit scale)
PLANTED = [("expression", 5, 1.2), ("expression", 15, -1.0), ("splicing", 25, 1.5)]


def select(a):
    rng = np.random.default_rng(SEED)
    ped = pd.read_csv(a.ped, sep=r"\s+")
    founders = ped[(ped.FatherID.astype(str) == "0") & (ped.MotherID.astype(str) == "0")]
    cohort, panel, panel_lab = [], [], []
    for sp in ["AFR", "AMR", "EAS", "EUR", "SAS"]:
        ids = sorted(founders.loc[founders.Superpopulation == sp, "SampleID"])
        ids = list(rng.permutation(ids))
        k = COHORT.get(sp, 0)
        cohort += ids[:k]
        if sp in PANEL_GROUPS:
            panel += ids[k:k + PANEL_PER_GROUP]
            panel_lab += [sp] * PANEL_PER_GROUP
    sex = dict(zip(ped.SampleID, ped.Sex))
    with open(os.path.join(a.tmp, "cohort_samples.txt"), "w") as f:
        f.write("\n".join(cohort) + "\n")
    with open(os.path.join(a.tmp, "panel_samples.txt"), "w") as f:
        f.write("\n".join(panel) + "\n")
    with open(os.path.join(a.tmp, "all_samples.txt"), "w") as f:
        f.write("\n".join(cohort + panel) + "\n")
    pd.DataFrame({"sample": panel, "superpop": panel_lab}).to_csv(
        os.path.join(a.tmp, "panel_superpop.tsv"), sep="\t", index=False)
    pd.DataFrame({"sample": cohort, "sex": [sex[s] for s in cohort],
                  "superpop": [ped.set_index("SampleID").Superpopulation[s] for s in cohort]}).to_csv(
        os.path.join(a.tmp, "cohort_info.tsv"), sep="\t", index=False)


def write_bed(df, samples, path):
    with gzip.open(path, "wt") as f:
        f.write("#chr\tstart\tend\tphenotype_id\t" + "\t".join(samples) + "\n")
        for _, r in df.iterrows():
            vals = "\t".join("NA" if (isinstance(v, float) and np.isnan(v)) else
                             (f"{v:.6g}" if isinstance(v, float) else str(v)) for v in r["values"])
            f.write(f"{r['chr']}\t{r['start']}\t{r['end']}\t{r['phenotype_id']}\t{vals}\n")


def build(a):
    rng = np.random.default_rng(SEED + 1)
    tmp, out = a.tmp, a.out
    info = pd.read_csv(os.path.join(tmp, "cohort_info.tsv"), sep="\t")
    cohort = list(info["sample"])
    vcf_samples = open(os.path.join(tmp, "gt_samples.txt")).read().split()
    gt = pd.read_csv(os.path.join(tmp, "gt.tsv"), sep="\t", header=None, dtype=str)
    meta_cols = ["chrom", "pos", "id", "ref", "alt"]
    gt.columns = meta_cols + vcf_samples
    gt["pos"] = gt["pos"].astype(int)
    gtc = gt[cohort]
    hap1 = gtc.apply(lambda c: c.str[0].astype(int)).to_numpy()
    hap2 = gtc.apply(lambda c: c.str[2].astype(int)).to_numpy()
    G = hap1 + hap2                                          # variants x samples, ALT dosage
    af = G.mean(axis=1) / 2
    keep = (af > 0) & (af < 1)
    gt, G, hap1, hap2, af = gt[keep].reset_index(drop=True), G[keep], hap1[keep], hap2[keep], af[keep]
    nv = len(gt)

    # synthetic imputation quality: ~8% of variants below the 0.4 filter
    r2 = np.where(rng.random(nv) < 0.08, rng.uniform(0.05, 0.39, nv), rng.uniform(0.6, 1.0, nv))
    typed = rng.random(nv) < 0.1
    r2[typed] = 1.0
    noise = rng.normal(0, 1, G.shape) * (0.4 * (1 - r2))[:, None]
    DS = np.clip(G + noise, 0, 2).round(3)
    # TOPMed-style VCF
    vcf_path = os.path.join(out, "genotypes", "chr22.dose.vcf")
    with open(vcf_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n##filedate=fixture\n##source=make_site_fixture.py "
                "(1kGP high-coverage genotypes; synthetic R2/DS)\n##contig=<ID=chr22,length=50818468>\n")
        f.write('##INFO=<ID=AF,Number=1,Type=Float,Description="Estimated Alternate Allele Frequency">\n')
        f.write('##INFO=<ID=MAF,Number=1,Type=Float,Description="Estimated Minor Allele Frequency">\n')
        f.write('##INFO=<ID=R2,Number=1,Type=Float,Description="Estimated Imputation Accuracy (R-square)">\n')
        f.write('##INFO=<ID=IMPUTED,Number=0,Type=Flag,Description="Marker was imputed">\n')
        f.write('##INFO=<ID=TYPED,Number=0,Type=Flag,Description="Marker was genotyped">\n')
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('##FORMAT=<ID=DS,Number=1,Type=Float,Description="Estimated Alternate Allele Dosage">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(cohort) + "\n")
        for i in range(nv):
            r = gt.iloc[i]
            flag = "TYPED" if typed[i] else "IMPUTED"
            inf = f"AF={af[i]:.5f};MAF={min(af[i], 1 - af[i]):.5f};R2={r2[i]:.5f};{flag}"
            cells = "\t".join(f"{hap1[i, j]}|{hap2[i, j]}:{DS[i, j]:.3f}" for j in range(len(cohort)))
            f.write(f"{r.chrom}\t{r.pos}\t{r.id}\t{r.ref}\t{r.alt}\t.\tPASS\t{inf}\tGT:DS\t{cells}\n")

    # ---- samples: RNA IDs, one technical replicate, one RNA sample without genotypes -----
    order = rng.permutation(len(cohort))
    rna_ids = [f"S{k + 1:03d}" for k in range(len(cohort))]
    rna_of = {cohort[order[k]]: rna_ids[k] for k in range(len(cohort))}
    meta = pd.DataFrame({"rnaseq_id": [rna_of[s] for s in cohort], "genotype_id": cohort})
    meta = meta.sort_values("rnaseq_id").reset_index(drop=True)
    eur = info.loc[info.superpop == "EUR", "sample"].iloc[0]
    meta = pd.concat([meta, pd.DataFrame({"rnaseq_id": ["S131", "S132"],
                                          "genotype_id": [eur, "NOGENO01"]})], ignore_index=True)
    sexmap = dict(zip(info["sample"], info["sex"].map({1: "M", 2: "F"})))
    meta["sex"] = [sexmap.get(g, "F") for g in meta.genotype_id]
    n = len(meta)
    meta["GA"] = np.clip(rng.normal(39, 1.2, n), 34, 42).round(1)
    meta["ppBMI"] = np.clip(rng.normal(26, 5, n), 16, 50).round(1)
    meta["gdm"] = (rng.random(n) < 0.15).astype(int)
    meta.to_csv(os.path.join(out, "metadata.tsv"), sep="\t", index=False)
    with open(os.path.join(out, "samplesheet.csv"), "w") as f:
        f.write("sample_id,fastq_1,fastq_2\n")
        for s in meta.rnaseq_id:
            f.write(f"{s},fastq/{s}_R1.fastq.gz,fastq/{s}_R2.fastq.gz\n")

    samples = list(meta.rnaseq_id)
    col = {g: i for i, g in enumerate(cohort)}
    # dosage matrix aligned to RNA samples (the RNA-only sample gets the cohort mean; it is
    # dropped by the pipeline because it has no genotypes)
    Gs = np.stack([G[:, col[g]] if g in col else G.mean(axis=1) for g in meta.genotype_id], axis=1)

    # ---- genes and planted variants --------------------------------------------------------
    tss = np.linspace(GENE_FROM, GENE_TO, N_GENES).astype(int)
    genes = [f"ENSGFIX{k + 1:07d}" for k in range(N_GENES)]
    dense = (gt.pos >= 19_000_000).to_numpy() & (gt.pos <= 21_500_000).to_numpy()
    sup = dict(zip(info["sample"], info["superpop"]))
    eur_idx = [col[g] for g in cohort if sup[g] == "EUR"]
    eas_idx = [col[g] for g in cohort if sup[g] == "EAS"]
    maf_eur = G[:, eur_idx].mean(1) / 2
    maf_eas = G[:, eas_idx].mean(1) / 2
    ok = dense & (r2 >= 0.8) & (np.minimum(maf_eur, 1 - maf_eur) > 0.2) & (np.minimum(maf_eas, 1 - maf_eas) > 0.2)
    truth = []
    causal = {}
    for mod, gi, eff in PLANTED:
        cand = np.where(ok & (np.abs(gt.pos.to_numpy() - tss[gi]) < 150_000))[0]
        v = int(cand[rng.integers(len(cand))])
        causal[(mod, gi)] = (v, eff)
        truth.append((mod, gi, eff, gt.id[v]))

    # ---- expression: counts (negative binomial) and TPM, with 3 hidden factors -------------
    F = rng.normal(0, 1, (3, n))
    L = rng.normal(0, 0.35, (N_GENES, 3))
    mu = rng.normal(6.5, 1.0, N_GENES)
    mu[[37, 38, 39]] = -4.0                                   # undetected genes (filtered)
    lib = rng.uniform(0.6, 1.6, n)
    logr = mu[:, None] + L @ F + rng.normal(0, 0.25, (N_GENES, n))
    for (mod, gi), (v, eff) in causal.items():
        if mod == "expression":
            logr[gi] += eff * np.log(2) * Gs[v]
    lam = np.exp(logr) * lib[None, :]
    disp = 0.05
    counts = rng.poisson(rng.gamma(1 / disp, lam * disp))
    # one sample with scrambled expression (a connectivity outlier)
    out_s = samples.index("S017")
    counts[:, out_s] = rng.permutation(counts[:, out_s])
    glen = rng.uniform(1, 40, N_GENES)
    rpk = counts / glen[:, None]
    tpm = rpk / rpk.sum(0, keepdims=True) * 1e6
    strand = ["+" if k % 2 == 0 else "-" for k in range(N_GENES)]
    gdf = pd.DataFrame({"chr": "chr22", "start": tss - 1, "end": tss, "phenotype_id": genes})
    ph = os.path.join(out, "phenotypes")
    write_bed(gdf.assign(values=[list(map(float, r)) for r in tpm]), samples,
              os.path.join(ph, "expression.bed.gz"))
    write_bed(gdf.assign(values=[list(map(int, r)) for r in counts]), samples,
              os.path.join(ph, "expression.counts.bed.gz"))

    # ---- splicing: one 3-intron cluster per gene, ratios from multinomial counts ----------
    rows, groups = [], []
    for k, g in enumerate(genes):
        if k in (37, 38, 39):
            continue
        base = rng.dirichlet([6, 3, 2])
        logit = np.log(base)[:, None] + rng.normal(0, 0.15, (3, n)) + 0.15 * (L[k, :1, None] * F[:1])
        if ("splicing", k) in causal:
            v, eff = causal[("splicing", k)]
            logit[0] += eff * Gs[v]
        p = np.exp(logit) / np.exp(logit).sum(0, keepdims=True)
        tot = rng.poisson(120, n)
        c = np.stack([rng.multinomial(tot[j], p[:, j]) for j in range(n)], axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            ratio = np.where(tot[None, :] > 0, c / tot[None, :], np.nan)
        s0 = int(tss[k]) + (500 if strand[k] == "+" else -9000)
        for i in range(3):
            st, en = s0 + 1000 * i, s0 + 1000 * i + 700 + 300 * i
            pid = f"{g}__chr22:{st}:{en}:clu_{k + 1}_{strand[k]}"
            # atlas convention: every phenotype BED row carries its gene's TSS (start = TSS - 1)
            rows.append({"chr": "chr22", "start": int(tss[k]) - 1, "end": int(tss[k]), "phenotype_id": pid,
                         "values": list(map(float, ratio[i]))})
            groups.append((pid, g))
            if ("splicing", k) in causal and i == 0:
                truth = [t if not (t[0] == "splicing" and t[1] == k) else (t[0], t[1], t[2], t[3], pid)
                         for t in truth]
    write_bed(pd.DataFrame(rows), samples, os.path.join(ph, "splicing.bed.gz"))
    pd.DataFrame(groups).to_csv(os.path.join(ph, "splicing.phenotype_groups.txt"), sep="\t",
                                header=False, index=False)

    # ---- sequencing QC metrics and cell-type proportions ----------------------------------
    intr = np.clip(25 + 4 * F[0] + rng.normal(0, 1.5, n), 5, 60)
    qc = pd.DataFrame({"sample": samples,
                       "total_reads": (lib * 30e6).round().astype(int),
                       "pct_uniquely_mapped": np.clip(88 + 2 * F[1] + rng.normal(0, 1, n), 60, 99).round(2),
                       "pct_intronic": intr.round(2),
                       "pct_exonic": (100 - intr - 5 + rng.normal(0, 0.3, n)).round(2),
                       "pct_rRNA": np.clip(1 + 0.5 * F[2] + rng.normal(0, 0.3, n), 0, 10).round(3),
                       "median_3prime_bias": np.clip(0.6 + 0.05 * rng.normal(0, 1, n), 0.2, 1.2).round(3)})
    qc.to_csv(os.path.join(ph, "qc_metrics.tsv"), sep="\t", index=False)
    props = rng.dirichlet([30, 15, 4, 12, 6, 5, 3], n)
    cp = pd.DataFrame(props.round(5), columns=CELL_TYPES)
    cp.insert(0, "sample", samples)
    cp.to_csv(os.path.join(ph, "cell_proportions.tsv"), sep="\t", index=False)

    tr = []
    for t in truth:
        mod, gi, eff, vid = t[:4]
        pid = t[4] if len(t) > 4 else genes[gi]
        tr.append({"modality": mod, "phenotype_id": pid, "gene_id": genes[gi], "variant_id": vid,
                   "effect_per_alt": eff})
    pd.DataFrame(tr).to_csv(os.path.join(out, "TRUTH.tsv"), sep="\t", index=False)
    with open(os.path.join(out, "README.txt"), "w") as f:
        f.write("Site-run test fixture written by tests/fixtures/make_site_fixture.sh.\n"
                "Genotypes: 130 unrelated 1000 Genomes high-coverage samples (60 EUR, 60 EAS, 10 AFR),\n"
                "chr22 subset, with SYNTHETIC imputation R2 and dosage noise.\n"
                "Phenotypes, QC metrics, cell proportions and covariates: SYNTHETIC SOFTWARE FIXTURES\n"
                "with planted cis effects listed in TRUTH.tsv. Not biological data.\n"
                "S131 is a technical replicate of another sample; S132 has no genotypes; S017 has\n"
                "scrambled expression (connectivity outlier).\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["select", "build"])
    ap.add_argument("--ped")
    ap.add_argument("--tmp", required=True)
    ap.add_argument("--out")
    a = ap.parse_args()
    select(a) if a.cmd == "select" else build(a)


if __name__ == "__main__":
    main()
