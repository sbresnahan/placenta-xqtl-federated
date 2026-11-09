"""Documentation coverage: every cohort.yml setting and every pipeline parameter is documented."""
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def schema_keys(node, prefix=""):
    for key, sub in (node.get("properties") or {}).items():
        path = f"{prefix}{key}"
        yield path
        yield from schema_keys(sub, path + ".")
        # map-like blocks (gxe.exposures, gxe.sensitivity): document the per-entry settings
        extra = sub.get("additionalProperties")
        if isinstance(extra, dict):
            for item in (extra.get("properties") or {}):
                yield f"{path}.<name>.{item}"


def test_every_schema_setting_in_configuration():
    schema = json.loads((ROOT / "assets/schema/cohort.schema.json").read_text())
    doc = (ROOT / "docs/CONFIGURATION.md").read_text()
    missing = []
    for key in schema_keys(schema):
        if "<name>." in key:
            base, item = key.split(".<name>.")
            if f"`{base}.<name>.{item}`" not in doc and f"`.{item}`" not in doc:
                missing.append(key)
        elif f"`{key}`" not in doc and f"`{key}." not in doc:   # a block counts when its keys are listed
            missing.append(key)
    assert not missing, f"undocumented settings in docs/CONFIGURATION.md: {missing}"


def test_every_param_documented():
    cfg = (ROOT / "nextflow.config").read_text()
    block = re.search(r"^params\s*\{(.*?)^\}", cfg, re.S | re.M).group(1)
    params = re.findall(r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=", block, re.M)
    docs = (ROOT / "docs/CONFIGURATION.md").read_text() + (ROOT / "docs/HPC.md").read_text()
    missing = [p for p in params if f"--{p}" not in docs and f"`{p}`" not in docs and p != "help"]
    assert not missing, f"undocumented parameters (docs/CONFIGURATION.md or docs/HPC.md): {missing}"
