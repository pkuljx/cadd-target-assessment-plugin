#!/usr/bin/env python
"""Within-family (paralog) identity/similarity table for one human target.

You supply the family/related human gene symbols (Claude knows GPCR, kinase, ion
channel, etc. families and can web-search to confirm membership); this aligns each
paralog's human sequence to the target and reports % identity and % similarity.
This is the "same family, different members" table — it shows how much of a
selectivity problem the closest relatives pose for a drug-design program.

Usage:
  python family_identity.py --ref-acc P35367 \
      --members "HRH2,HRH3,HRH4,CHRM1,HTR2A" [--out family.json]
  python family_identity.py --ref-gene HRH1 --members "HRH2,HRH3,HRH4"
"""

import argparse

import cadd_common as C


def resolve_human(gene):
    return C.best_entry(
        C.uniprot_search(
            f"gene:{gene} AND organism_id:9606 AND reviewed:true",
            "accession,length,reviewed,protein_name",
        )
    ) or C.best_entry(
        C.uniprot_search(
            f"gene:{gene} AND organism_id:9606", "accession,length,reviewed,protein_name"
        )
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref-acc")
    ap.add_argument("--ref-gene")
    ap.add_argument("--members", required=True, help="comma-separated human gene symbols")
    ap.add_argument("--out")
    args = ap.parse_args()

    if args.ref_acc:
        ref = C.uniprot_entry(args.ref_acc)
    elif args.ref_gene:
        ref = resolve_human(args.ref_gene)
    else:
        ap.error("provide --ref-acc or --ref-gene")

    ref_acc = ref["primaryAccession"]
    ref_seq = ref["sequence"]["value"]
    ref_gene = (ref.get("genes", [{}])[0].get("geneName", {}) or {}).get("value")

    rows = [{
        "gene": ref_gene, "accession": ref_acc, "identity_pct": 100.0,
        "similarity_pct": 100.0, "is_target": True,
    }]
    for gene in [m.strip() for m in args.members.split(",") if m.strip()]:
        m = resolve_human(gene)
        if not m:
            rows.append({"gene": gene, "accession": None, "note": "not resolved"})
            continue
        macc = m["primaryAccession"]
        mseq = C.seq_of(macc)
        idn, sim = C.identity(ref_seq, mseq)
        rows.append({
            "gene": gene,
            "accession": macc,
            "protein_name": (m.get("proteinDescription", {})
                             .get("recommendedName", {})
                             .get("fullName", {}).get("value")),
            "length": m.get("sequence", {}).get("length"),
            "identity_pct": idn,
            "similarity_pct": sim,
        })

    C.emit({
        "target_gene": ref_gene,
        "target_accession": ref_acc,
        "method": "global BLOSUM62 align; identity/similarity over aligned non-gap columns",
        "members": rows,
    }, args.out)


if __name__ == "__main__":
    main()
