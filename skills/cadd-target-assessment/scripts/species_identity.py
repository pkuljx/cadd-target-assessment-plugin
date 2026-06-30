#!/usr/bin/env python
"""Cross-species sequence-identity table for one target.

For the human protein, find the ortholog in mouse / rat / dog / monkey (by gene
symbol), align each to the human sequence, and report % identity. This is the
"species homology" line in the report — it tells you whether an animal model
shares the binding site well enough for the program's in-vivo work to translate.

Usage:
  python species_identity.py --gene HRH1 [--ref-acc P35367] [--out species.json]
  python species_identity.py --ref-acc P35367            # gene read from the entry

Override the species set with --species "mouse:10090,rat:10116,...".
"""

import argparse

import cadd_common as C


def parse_species(s):
    out = []
    for tok in s.split(","):
        name, _, tid = tok.partition(":")
        out.append((name.strip(), int(tid)))
    return out


def ortholog(gene, organism):
    rows = C.uniprot_search(
        f"gene:{gene} AND organism_id:{organism}",
        "accession,length,reviewed,protein_name",
        size=15,
    )
    return C.best_entry(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", help="human gene symbol (e.g. HRH1)")
    ap.add_argument("--ref-acc", help="human UniProt accession")
    ap.add_argument("--species", help='override, e.g. "mouse:10090,rat:10116"')
    ap.add_argument("--out")
    args = ap.parse_args()

    if not args.gene and not args.ref_acc:
        ap.error("provide --gene and/or --ref-acc")

    ref_entry = None
    if args.ref_acc:
        ref_entry = C.uniprot_entry(args.ref_acc)
        if not args.gene:
            genes = ref_entry.get("genes", [{}])
            args.gene = genes[0].get("geneName", {}).get("value")
    if ref_entry is None:
        ref_entry = C.best_entry(
            C.uniprot_search(
                f"gene:{args.gene} AND organism_id:9606 AND reviewed:true",
                "accession,length,reviewed",
            )
        )
    ref_acc = ref_entry["primaryAccession"]
    ref_seq = ref_entry["sequence"]["value"]

    species = parse_species(args.species) if args.species else C.SPECIES
    rows = []
    for name, tid in species:
        o = ortholog(args.gene, tid)
        if not o:
            rows.append({"species": name, "organism_id": tid, "accession": None,
                         "note": "no ortholog found"})
            continue
        oacc = o["primaryAccession"]
        oseq = C.seq_of(oacc)
        idn, sim = C.identity(ref_seq, oseq)
        rows.append({
            "species": name,
            "organism_id": tid,
            "accession": oacc,
            "reviewed": C.is_reviewed(o),
            "length": o.get("sequence", {}).get("length"),
            "identity_pct": idn,
            "similarity_pct": sim,
        })

    C.emit({
        "gene": args.gene,
        "reference_accession": ref_acc,
        "reference_length": len(ref_seq),
        "method": "global BLOSUM62 align; identity over aligned non-gap columns",
        "species": rows,
    }, args.out)


if __name__ == "__main__":
    main()
