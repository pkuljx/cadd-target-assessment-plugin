"""Shared helpers for the CADD target-assessment skill.

Stdlib-only HTTP (urllib) + Biopython pairwise alignment. The only third-party
dependency is Biopython (>=1.80); everything else is the standard library, so any
Python 3.9+ with `pip install biopython` can run these (no `requests` needed).

All network calls hit public REST endpoints:
  - UniProt   https://rest.uniprot.org
  - RCSB data https://data.rcsb.org
  - RCSB search https://search.rcsb.org
"""

import json
import sys
import time
import urllib.parse
import urllib.request
import warnings

warnings.filterwarnings("ignore")  # silence Biopython deprecation chatter

UA = {"User-Agent": "cadd-target-assessment/1.0 (research script)"}

# Species used for the cross-species conservation table. Override on the CLI if a
# target's relevant model organism differs.
SPECIES = [
    ("mouse", 10090),
    ("rat", 10116),
    ("dog", 9615),
    ("monkey", 9544),  # Macaca mulatta (rhesus); 9541 = crab-eating macaque
]

# Ligands that are almost always crystallization / cryo / buffer additives rather
# than the biologically interesting bound molecule. Filtered out of PDB ligand
# lists so the "key compound" per structure is easy to see. Not exhaustive — when
# in doubt the structure title usually names the real ligand.
ADDITIVES = {
    "HOH", "DOD",                                   # water
    "NA", "K", "MG", "CA", "ZN", "MN", "CL", "BR", "IOD", "CD", "NI", "CU", "FE",
    "SO4", "PO4", "PEG", "EDO", "GOL", "MPD", "DMS", "ACT", "ACE", "FMT",
    "TRS", "EPE", "MES", "BTB", "IMD", "CIT", "TLA", "MLI",
    "OLC", "OLA", "OLB", "PLM", "MYR", "STE",        # monoolein / lipidic-cubic-phase
    "Y01", "CLR", "CHS", "PGE", "P6G", "1PE", "PG4", "2PE", "12P", "PE4",
    "BOG", "LMT", "LMN", "DMU", "9MA", "C8E", "BNG", "OGA",  # detergents
    "NAG", "BMA", "MAN", "FUC", "GAL",               # common N-glycan sugars
    "EDT", "BME", "DTT", "GSH",
}


def http_json(url, retries=3, timeout=45):
    """GET a URL and parse JSON, with simple backoff retry."""
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.load(r)
        except Exception as e:  # noqa: BLE001 - network flakiness, retry
            last = e
            time.sleep(1.5 * (i + 1))
    raise RuntimeError(f"GET failed after {retries} tries: {url}\n{last!r}")


def http_json_post(url, payload, retries=3, timeout=45):
    """POST JSON and parse the JSON response."""
    data = json.dumps(payload).encode()
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(
                url, data=data, headers={**UA, "Content-Type": "application/json"}
            )
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.load(r)
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(1.5 * (i + 1))
    raise RuntimeError(f"POST failed after {retries} tries: {url}\n{last!r}")


# ----------------------------- UniProt -----------------------------

UNIPROT = "https://rest.uniprot.org/uniprotkb"


def uniprot_search(query, fields, size=10):
    url = (
        f"{UNIPROT}/search?query={urllib.parse.quote(query)}"
        f"&fields={fields}&format=json&size={size}"
    )
    return http_json(url).get("results", [])


def uniprot_entry(acc):
    return http_json(f"{UNIPROT}/{acc}.json")


def is_reviewed(entry):
    return entry.get("entryType", "").startswith("UniProtKB reviewed")


def best_entry(entries):
    """Prefer a reviewed (Swiss-Prot) entry; else the longest TrEMBL entry.

    Reviewed entries are curated and isoform-canonical, so they give the most
    comparable sequence. Many non-human orthologs (dog, macaque) only exist as
    TrEMBL predictions — there we fall back to the longest, which is usually the
    full-length model rather than a fragment.
    """
    if not entries:
        return None
    reviewed = [e for e in entries if is_reviewed(e)]
    pool = reviewed or entries
    return max(pool, key=lambda e: e.get("sequence", {}).get("length", 0))


def seq_of(acc):
    return uniprot_entry(acc)["sequence"]["value"]


# ----------------------------- alignment -----------------------------

_ALIGNER = None


def _aligner():
    global _ALIGNER
    if _ALIGNER is None:
        from Bio.Align import PairwiseAligner, substitution_matrices

        a = PairwiseAligner()
        a.substitution_matrix = substitution_matrices.load("BLOSUM62")
        a.open_gap_score = -11
        a.extend_gap_score = -1
        a.mode = "global"
        _ALIGNER = a
    return _ALIGNER


def identity(seq_a, seq_b):
    """Return (percent_identity, percent_similarity) for two protein sequences.

    Global Needleman-Wunsch alignment with BLOSUM62 (gap open -11, extend -1) —
    the standard general-purpose protein setting. Identity and similarity are
    counted over aligned columns where BOTH sequences have a residue (gaps
    excluded from the denominator), which is the most common way "% identity" is
    reported for orthologs/paralogs. Similarity additionally counts conservative
    substitutions (BLOSUM62 score > 0).
    """
    from Bio.Align import substitution_matrices

    mat = substitution_matrices.load("BLOSUM62")
    aln = _aligner().align(seq_a, seq_b)[0]
    top, bot = aln[0], aln[1]
    ident = sim = cols = 0
    for x, y in zip(top, bot):
        if x == "-" or y == "-":
            continue
        cols += 1
        if x == y:
            ident += 1
            sim += 1
        else:
            try:
                if mat[x, y] > 0:
                    sim += 1
            except KeyError:
                pass
    if cols == 0:
        return 0.0, 0.0
    return round(100 * ident / cols, 1), round(100 * sim / cols, 1)


def emit(obj, out=None):
    """Print JSON to stdout and optionally also write it to a file."""
    text = json.dumps(obj, indent=2, ensure_ascii=False)
    if out:
        with open(out, "w") as f:
            f.write(text)
        print(f"[written] {out}", file=sys.stderr)
    print(text)
