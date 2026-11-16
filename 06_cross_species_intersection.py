"""Intersect selected human DEGs, audited mouse orthologues, and MR genes.

Run from any directory: python scripts/06_cross_species_intersection.py
Supply real MGI-derived mappings and the HGNC snapshot; no online lookup occurs.
"""
from pathlib import Path
from collections import defaultdict
import hashlib
import json
import math
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/intersection"
PATHS = {
    "mouse": ROOT / "results/differential_expression/mouse_selected.csv",
    "human": ROOT / "results/differential_expression/human_selected.csv",
    "mr": ROOT / "data/mr/prioritized_genes.txt",
    "mapping": ROOT / "data/orthology/mouse_human_mapping.csv",
    "hgnc": ROOT / "data/orthology/hgnc_complete_set.txt",
}
def read_table(path, required, sep=","):
    frame = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False)
    missing = set(required) - set(frame.columns)
    if missing:
        raise ValueError(f"{path.name}: missing columns {sorted(missing)}")
    for column in required:
        frame[column] = frame[column].str.strip()
    return frame

hg = read_table(PATHS["hgnc"], ["symbol", "hgnc_id", "prev_symbol", "alias_symbol"], "\t")
if not hg.symbol.is_unique or not hg.hgnc_id.is_unique:
    raise ValueError("HGNC snapshot must contain unique approved symbols and HGNC IDs.")
exact = dict(zip(hg.symbol, hg.hgnc_id))
by_id = dict(zip(hg.hgnc_id, hg.symbol))
previous, aliases = defaultdict(set), defaultdict(set)
for row in hg.itertuples():
    for column, lookup in [("prev_symbol", previous), ("alias_symbol", aliases)]:
        for symbol in getattr(row, column).split("|"):
            if symbol:
                lookup[symbol].add(row.hgnc_id)
def normalize(symbol):
    if symbol in exact:
        return exact[symbol], by_id[exact[symbol]], "approved"
    for lookup, label in [(previous, "previous_symbol"), (aliases, "unique_alias")]:
        hits = lookup.get(symbol, set())
        if len(hits) == 1:
            key = next(iter(hits))
            return key, by_id[key], label
        if len(hits) > 1:
            return "UNRESOLVED:" + symbol, symbol, "ambiguous_" + label
    return "UNRESOLVED:" + symbol, symbol, "unresolved"
def annotate(frame, column):
    values = [normalize(s) for s in frame[column]]
    frame[["human_key", "human_symbol", "symbol_status"]] = pd.DataFrame(values, index=frame.index, columns=["human_key", "human_symbol", "symbol_status"])
    return frame
def read_de(name):
    frame = read_table(PATHS[name], ["gene", "log2FC", "p_value", "p_adj"])
    if not frame.gene.is_unique or frame.gene.eq("").any():
        raise ValueError(f"{name}: selected gene symbols must be nonempty and unique.")
    for column in ["log2FC", "p_value", "p_adj"]:
        frame[column] = pd.to_numeric(frame[column], errors="raise")
        if not frame[column].map(math.isfinite).all():
            raise ValueError(f"{name}: nonfinite {column}.")
    if not frame[["p_value", "p_adj"]].apply(lambda x: x.between(0, 1)).all().all():
        raise ValueError(f"{name}: P values must be between zero and one.")
    return frame

mouse = read_de("mouse").rename(columns={c: "mouse_" + c for c in ["gene", "log2FC", "p_value", "p_adj"]})
human = annotate(read_de("human"), "gene")
mr = annotate(pd.DataFrame({"gene": sorted(set(s.strip() for s in PATHS["mr"].read_text().splitlines() if s.strip()))}), "gene")
# The supplied audit must retain all MGI homologues. Its producer should match
# mouse NCBI Gene ID first, then mouse symbol, and retain HGNC provenance.
mapping = annotate(read_table(PATHS["mapping"], ["mouse_symbol", "human_symbol"]).rename(columns={"human_symbol": "human_source_symbol"}), "human_source_symbol")
if mapping.mouse_symbol.eq("").any() or mapping.human_source_symbol.eq("").any():
    raise ValueError("The audited mapping contains an empty source or target symbol.")
all_rows = mouse.merge(mapping, left_on="mouse_gene", right_on="mouse_symbol", how="left", indicator=True)
mapped = all_rows[all_rows._merge.eq("both")].drop(columns="_merge")
unmapped = all_rows[all_rows._merge.eq("left_only")][mouse.columns]
pairs = mapped.drop_duplicates(["mouse_gene", "human_key"])
sets = [set(human.human_key), set(pairs.human_key), set(mr.human_key)]
core = set.intersection(*sets)
membership = pd.DataFrame([{"human_key": k, "human_symbol": by_id.get(k, k.removeprefix("UNRESOLVED:")), "region": "".join(str(int(k in s)) for s in sets)} for k in sorted(set.union(*sets))], columns=["human_key", "human_symbol", "region"])
details = human[human.human_key.isin(core)].merge(pairs, on="human_key", suffixes=("_human", "_mouse"))
OUT.mkdir(parents=True, exist_ok=True)
for name, frame in {"human_symbol_audit": human, "MR_symbol_audit": mr, "orthology_input_audit": mapping, "mouse_mapping_all_rows": mapped, "mouse_human_mapping": pairs, "mouse_unmapped": unmapped, "exclusive_regions": membership, "core_gene_details": details}.items():
    frame.to_csv(OUT / (name + ".csv"), index=False)
regions = ["100", "010", "001", "110", "101", "011", "111"]
membership.region.value_counts().reindex(regions, fill_value=0).rename_axis("region").rename("gene_count").to_csv(OUT / "region_counts.csv")
names = sorted(by_id.get(k, k.removeprefix("UNRESOLVED:")) for k in core)
(OUT / "core_genes.txt").write_text("".join(g + "\n" for g in names))
summary = dict(zip(["human_DEGs", "mouse_mapped_human_genes", "MR_genes"], map(len, sets)))
summary.update(core_genes=names, unmapped_mouse_genes=len(unmapped), direction_filter="none", region_order="human,mouse_orthologues,MR", unresolved_core_keys=sorted(k for k in core if k.startswith("UNRESOLVED:")))
summary["input_sha256"] = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in PATHS.values()}
(OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
