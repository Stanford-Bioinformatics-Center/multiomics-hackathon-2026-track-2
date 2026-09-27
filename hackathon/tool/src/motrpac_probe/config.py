"""Config files (YAML or JSON) that mirror every CLI flag.

Method details: docs/METHODS.md#configpy
"""
import json
from pathlib import Path

LISTS = {"context_groups", "tissues", "pool_sets", "pool_groups", "sections", "layers"}


def _scalar(v):
    v = v.strip().strip('"').strip("'")
    if v.lower() in ("true", "false"):
        return v.lower() == "true"
    for f in (int, float):
        try:
            return f(v)
        except ValueError:
            pass
    return v


def _flat_yaml(text):
    out, block = {}, None
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].rstrip()
        if not line.strip():
            continue
        if not line.startswith(" ") and line.endswith(":"):
            block = line[:-1].strip()
            out[block] = {}
            continue
        k, _, v = line.strip().partition(":")
        v = v.strip()
        if v.startswith("[") and v.endswith("]"):
            val = [_scalar(x) for x in v[1:-1].split(",") if x.strip()]
        else:
            val = _scalar(v)
        (out[block] if block and line.startswith(" ") else out)[k.strip()] = val
    return out


def load_config(path, section):
    text = Path(path).read_text(encoding="utf-8")
    if str(path).endswith(".json"):
        cfg = json.loads(text)
    else:
        try:
            import yaml
            cfg = yaml.safe_load(text)
        except ImportError:
            cfg = _flat_yaml(text)
    cfg = cfg.get(section, cfg) if isinstance(cfg, dict) else {}
    out = {}
    for k, v in cfg.items():
        k = k.replace("-", "_")
        if k in LISTS and isinstance(v, str):
            v = [x.strip() for x in v.split(",") if x.strip()]
        out[k] = tuple(v) if isinstance(v, list) else v
    return out
