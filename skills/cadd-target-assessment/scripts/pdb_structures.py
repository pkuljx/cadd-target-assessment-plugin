#!/usr/bin/env python
"""Enumerate every experimental PDB structure of a target and summarize each.

Queries RCSB for all entries whose polymer entity references the target's UniProt
accession (this catches truncated constructs that a plain sequence search misses),
then for each entry pulls method, resolution, the modeled residue range in UniProt
numbering, the bound ligand(s) with buffer/cryo additives filtered out, and the
primary citation (title, DOI, PubMed, year).

The ligand + citation fields are what let the assessment answer "is this compound
an agonist or an antagonist, and what does the paper say about the mechanism" —
so the per-structure citation is deliberately included for downstream lookup.

Usage:
  python pdb_structures.py --acc P35367 [--out pdb.json] [--max 200]
"""

import argparse

import cadd_common as C

DATA = "https://data.rcsb.org/rest/v1/core"
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"


def entry_ids_for_uniprot(acc, limit):
    payload = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": (
                    "rcsb_polymer_entity_container_identifiers."
                    "reference_sequence_identifiers.database_accession"
                ),
                "operator": "exact_match",
                "value": acc,
            },
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": limit},
            "results_content_type": ["experimental"],
            "sort": [{"sort_by": "rcsb_accession_info.initial_release_date", "direction": "asc"}],
        },
    }
    res = C.http_json_post(f"{SEARCH}?", payload)
    ids = [r["identifier"] for r in res.get("result_set", [])]
    return ids, res.get("total_count", len(ids))


def ligands(entry):
    ids = entry.get("rcsb_entry_container_identifiers", {}).get(
        "non_polymer_entity_ids", []
    ) or []
    out = []
    for nid in ids:
        try:
            ne = C.http_json(f"{DATA}/nonpolymer_entity/{entry['_pdb']}/{nid}")
        except Exception:  # noqa: BLE001
            continue
        comp = ne.get("pdbx_entity_nonpoly", {})
        cid = comp.get("comp_id")
        if not cid or cid in C.ADDITIVES:
            continue
        out.append({"comp_id": cid, "name": comp.get("name")})
    return out


def target_positions(pdb_id, acc):
    """Modeled residue range(s) in UniProt numbering for the entity matching acc."""
    try:
        eids = C.http_json(f"{DATA}/entry/{pdb_id}")[
            "rcsb_entry_container_identifiers"
        ].get("polymer_entity_ids", [])
    except Exception:  # noqa: BLE001
        return None
    for eid in eids:
        try:
            pe = C.http_json(f"{DATA}/polymer_entity/{pdb_id}/{eid}")
        except Exception:  # noqa: BLE001
            continue
        refs = pe.get("rcsb_polymer_entity_container_identifiers", {}).get(
            "reference_sequence_identifiers", []
        ) or []
        if not any(r.get("database_accession") == acc for r in refs):
            continue
        regions = []
        for al in pe.get("rcsb_polymer_entity_align", []) or []:
            if al.get("reference_database_accession") != acc:
                continue
            for reg in al.get("aligned_regions", []) or []:
                beg = reg.get("ref_beg_seq_id")
                length = reg.get("length")
                if beg and length:
                    regions.append((beg, beg + length - 1))
        auth_chains = pe.get("rcsb_polymer_entity_container_identifiers", {}).get(
            "auth_asym_ids", []
        )
        if regions:
            regions.sort()
            # Construct-engineering artifacts (fusion-partner residues, single
            # modeled point mutations) show up as 1-2 residue "aligned regions"
            # that clutter the range. Drop fragments shorter than 5 residues and
            # merge spans separated by gaps <5, so the output reads like the
            # domain boundaries a crystallographer would quote (e.g. 20-217/408-483).
            merged = []
            for b, e in regions:
                if merged and b - merged[-1][1] <= 5:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], e))
                else:
                    merged.append((b, e))
            merged = [(b, e) for b, e in merged if e - b + 1 >= 5]
            pos = "/".join(f"{b}-{e}" for b, e in merged) or None
            return {"positions": pos, "chains": auth_chains}
        return {"positions": None, "chains": auth_chains}
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--acc", required=True)
    ap.add_argument("--max", type=int, default=200)
    ap.add_argument("--out")
    args = ap.parse_args()

    ids, total = entry_ids_for_uniprot(args.acc, args.max)
    rows = []
    for pdb_id in ids:
        try:
            entry = C.http_json(f"{DATA}/entry/{pdb_id}")
        except Exception:  # noqa: BLE001
            rows.append({"pdb_id": pdb_id, "note": "fetch failed"})
            continue
        entry["_pdb"] = pdb_id
        info = entry.get("rcsb_entry_info", {})
        res = info.get("resolution_combined")
        cit = entry.get("rcsb_primary_citation", {}) or {}
        pos = target_positions(pdb_id, args.acc) or {}
        rows.append({
            "pdb_id": pdb_id,
            "method": (entry.get("exptl", [{}])[0].get("method")),
            "resolution_A": (res[0] if isinstance(res, list) and res else None),
            "chains": pos.get("chains"),
            "positions": pos.get("positions"),
            "ligands": ligands(entry),
            "title": entry.get("struct", {}).get("title"),
            "citation": {
                "title": cit.get("title"),
                "doi": cit.get("pdbx_database_id_doi"),
                "pubmed": cit.get("pdbx_database_id_pub_med"),
                "year": cit.get("year"),
                "journal": cit.get("rcsb_journal_abbrev"),
            },
        })

    C.emit({
        "accession": args.acc,
        "total_structures": total,
        "returned": len(rows),
        "structures": rows,
    }, args.out)


if __name__ == "__main__":
    main()
