#!/usr/bin/env python
"""Resolve a target to a single UniProt protein and dump its core facts.

Usage:
  python protein_info.py --query "HRH1" [--organism 9606] [--out protein.json]
  python protein_info.py --acc P35367 [--out protein.json]
  python protein_info.py --query "H1 receptor" --list      # disambiguation only

--list prints the top candidate proteins (accession, gene, name, organism) WITHOUT
committing to one. Use it when a target name is ambiguous so the human can confirm
which protein is meant before the full assessment runs.

Default organism is human (9606); the cross-species table is built separately by
species_identity.py.
"""

import argparse

import cadd_common as C

FIELDS = (
    "accession,id,gene_names,protein_name,length,mass,organism_name,"
    "cc_subcellular_location,cc_function,cc_similarity,keyword,xref_pdb"
)


def candidates(query, organism):
    org = f" AND organism_id:{organism}" if organism else ""
    # Try gene-symbol match first (most precise), then a general name search.
    rows = C.uniprot_search(
        f"(gene:{query}){org} AND reviewed:true", FIELDS, size=10
    )
    if not rows:
        rows = C.uniprot_search(
            f"({query}){org} AND reviewed:true", FIELDS, size=10
        )
    if not rows:
        rows = C.uniprot_search(f"({query}){org}", FIELDS, size=10)
    return rows


def summarize_candidate(e):
    genes = e.get("genes", [{}])
    gene = genes[0].get("geneName", {}).get("value") if genes else None
    name = (
        e.get("proteinDescription", {})
        .get("recommendedName", {})
        .get("fullName", {})
        .get("value")
    )
    return {
        "accession": e.get("primaryAccession"),
        "gene": gene,
        "protein_name": name,
        "organism": e.get("organism", {}).get("scientificName"),
        "length": e.get("sequence", {}).get("length"),
        "reviewed": C.is_reviewed(e),
    }


def subcellular(e):
    locs = []
    for cc in e.get("comments", []):
        if cc.get("commentType") == "SUBCELLULAR LOCATION":
            for sl in cc.get("subcellularLocations", []):
                v = sl.get("location", {}).get("value")
                if v:
                    locs.append(v)
    return locs


def cc_text(e, ctype):
    out = []
    for cc in e.get("comments", []):
        if cc.get("commentType") == ctype:
            for t in cc.get("texts", []):
                if t.get("value"):
                    out.append(t["value"])
    return out


def alt_names(e):
    pd = e.get("proteinDescription", {})
    names = []
    for a in pd.get("alternativeNames", []):
        v = a.get("fullName", {}).get("value")
        if v:
            names.append(v)
    return names


def full_info(e):
    genes = e.get("genes", [{}])
    gene = genes[0].get("geneName", {}).get("value") if genes else None
    syn = [s.get("value") for s in (genes[0].get("synonyms", []) if genes else [])]
    pd = e.get("proteinDescription", {})
    name = pd.get("recommendedName", {}).get("fullName", {}).get("value")
    pdb_ids = [
        x["id"] for x in e.get("uniProtKBCrossReferences", []) if x.get("database") == "PDB"
    ]
    return {
        "accession": e.get("primaryAccession"),
        "uniprot_id": e.get("uniProtkbId") or e.get("uniProtKBId"),
        "gene": gene,
        "gene_synonyms": syn,
        "protein_name": name,
        "alt_names": alt_names(e),
        "organism": e.get("organism", {}).get("scientificName"),
        "length": e.get("sequence", {}).get("length"),
        "mass_da": e.get("sequence", {}).get("molWeight"),
        "mass_kda": round(e.get("sequence", {}).get("molWeight", 0) / 1000.0, 1),
        "subcellular_location": subcellular(e),
        "family": cc_text(e, "SIMILARITY"),
        "function": cc_text(e, "FUNCTION"),
        "keywords": [k.get("name") for k in e.get("keywords", [])],
        "pdb_count_uniprot": len(pdb_ids),
        "reviewed": C.is_reviewed(e),
        "sequence": e.get("sequence", {}).get("value"),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--query", help="gene symbol or protein name")
    ap.add_argument("--acc", help="UniProt accession (skips resolution)")
    ap.add_argument("--organism", default="9606", help="taxon id (default human 9606)")
    ap.add_argument("--list", action="store_true", help="list candidates only")
    ap.add_argument("--out")
    args = ap.parse_args()

    if args.acc:
        C.emit(full_info(C.uniprot_entry(args.acc)), args.out)
        return

    if not args.query:
        ap.error("provide --query or --acc")

    cands = candidates(args.query, args.organism)
    if not cands:
        C.emit({"error": f"no UniProt match for '{args.query}'", "candidates": []})
        return

    if args.list:
        C.emit({"query": args.query, "candidates": [summarize_candidate(e) for e in cands]})
        return

    chosen = C.best_entry(cands)
    # The search response omits sequence and most comment blocks (only the
    # requested `fields` come back), so re-fetch the full entry by accession.
    full = C.uniprot_entry(chosen["primaryAccession"])
    info = full_info(full)
    if len(cands) > 1:
        info["other_candidates"] = [summarize_candidate(e) for e in cands if e is not chosen][:5]
    C.emit(info, args.out)


if __name__ == "__main__":
    main()
