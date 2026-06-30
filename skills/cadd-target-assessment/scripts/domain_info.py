#!/usr/bin/env python
"""Decompose a target into its domain / architectural features from UniProt.

Usage:
  python domain_info.py --acc Q8TDX7 [--out domains.json]
  python domain_info.py --query "NEK7" [--organism 9606] [--out domains.json]

Answers three questions for the report: how many domains the protein has, what
each one does, and the residue range (start-end) of each. All numbers come from
the UniProt feature table of the reviewed entry, so they're reproducible. When
UniProt's annotation is sparse (few/no architectural features), the `note` field
flags that the report writer should supplement from the literature / InterPro and
attribute the source.

Feature types kept are the architectural ones (Domain, Repeat, Zinc finger,
DNA binding, Transmembrane, Intramembrane, Topological domain, Region, Coiled
coil, Motif). Residue-level annotations (Active site, Binding site, Modified
residue, secondary-structure Helix/Strand, etc.) and the whole-chain `Chain`
feature are intentionally excluded — this is an architecture view, not a site list.
"""

import argparse

import cadd_common as C

# UniProt feature `type` strings that describe protein architecture (ordered by
# how prominently they define a fold). `domain_count` counts only true "Domain".
ARCH_TYPES = [
    "Domain",
    "Repeat",
    "Zinc finger",
    "DNA binding",
    "Transmembrane",
    "Intramembrane",
    "Topological domain",
    "Region",
    "Coiled coil",
    "Motif",
]

# Cross-reference databases that name the domain family behind a feature.
FAMILY_DBS = ("Pfam", "InterPro", "SMART")


def _pos(loc, end=False):
    """Pull a residue position out of a UniProt feature location, or None."""
    node = loc.get("end" if end else "start", {})
    return node.get("value")


def architectural_features(entry):
    feats = []
    for f in entry.get("features", []):
        ftype = f.get("type")
        if ftype not in ARCH_TYPES:
            continue
        loc = f.get("location", {})
        start, end = _pos(loc), _pos(loc, end=True)
        length = (end - start + 1) if (start is not None and end is not None) else None
        feats.append(
            {
                "type": ftype,
                "name": f.get("description") or ftype,
                "start": start,
                "end": end,
                "length": length,
                "evidence_count": len(f.get("evidences", [])),
            }
        )
    # Read N-terminus -> C-terminus; unplaced features sort last.
    feats.sort(key=lambda d: (d["start"] is None, d["start"] or 0))
    return feats


def domain_notes(entry):
    out = []
    for cc in entry.get("comments", []):
        if cc.get("commentType") == "DOMAIN":
            for t in cc.get("texts", []):
                if t.get("value"):
                    out.append(t["value"])
    return out


def domain_families(entry):
    out = []
    for x in entry.get("uniProtKBCrossReferences", []):
        db = x.get("database")
        if db not in FAMILY_DBS:
            continue
        name = None
        for p in x.get("properties", []):
            if p.get("key") == "EntryName":
                name = p.get("value")
        out.append({"database": db, "id": x.get("id"), "name": name})
    return out


def build(entry):
    genes = entry.get("genes", [{}])
    gene = genes[0].get("geneName", {}).get("value") if genes else None
    feats = architectural_features(entry)
    domain_count = sum(1 for f in feats if f["type"] == "Domain")

    note = (
        "UniProt feature table lists no annotated architectural features; describe "
        "the domain organization from the literature / InterPro and attribute the "
        "source in the report."
        if not feats
        else (
            "Only a single architectural feature in UniProt — if the protein is "
            "multi-domain, supplement from the literature / InterPro and attribute it."
            if len(feats) == 1 and domain_count <= 1
            else None
        )
    )

    return {
        "accession": entry.get("primaryAccession"),
        "gene": gene,
        "length": entry.get("sequence", {}).get("length"),
        "reviewed": C.is_reviewed(entry),
        "source": "UniProt feature table"
        + (" (reviewed entry)" if C.is_reviewed(entry) else " (unreviewed entry)"),
        "domain_count": domain_count,
        "feature_count": len(feats),
        "features": feats,
        "domain_notes": domain_notes(entry),
        "domain_families": domain_families(entry),
        "note": note,
    }


def resolve_acc(query, organism):
    """Resolve a name/gene to a single accession, mirroring protein_info.py."""
    org = f" AND organism_id:{organism}" if organism else ""
    fields = "accession,id,gene_names,protein_name,length,organism_name"
    for q in (f"(gene:{query}){org} AND reviewed:true", f"({query}){org} AND reviewed:true", f"({query}){org}"):
        rows = C.uniprot_search(q, fields, size=10)
        if rows:
            best = C.best_entry(rows)
            return best.get("primaryAccession")
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--acc", help="UniProt accession (skips resolution)")
    ap.add_argument("--query", help="gene symbol or protein name")
    ap.add_argument("--organism", default="9606", help="taxon id (default human 9606)")
    ap.add_argument("--out")
    args = ap.parse_args()

    acc = args.acc
    if not acc:
        if not args.query:
            ap.error("provide --acc or --query")
        acc = resolve_acc(args.query, args.organism)
        if not acc:
            C.emit({"error": f"no UniProt match for '{args.query}'"}, args.out)
            return

    C.emit(build(C.uniprot_entry(acc)), args.out)


if __name__ == "__main__":
    main()
