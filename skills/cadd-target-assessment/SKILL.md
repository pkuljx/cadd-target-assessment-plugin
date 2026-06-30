---
name: cadd-target-assessment
description: >-
  Produce a computational drug-discovery (CADD) target-assessment dossier for a
  protein target, output as a single self-contained HTML report. Covers protein
  basics (gene, length, molecular weight, subcellular location, biological
  function), cross-species conservation (mouse/rat/dog/monkey identity), within-
  family/paralog identity, the catalog of PDB structures (with bound ligands and
  agonist/antagonist annotation), the structural mechanism of small-molecule
  modulation, and a final computational pros/cons + feasibility verdict. Use this
  whenever the user names a protein/drug target and asks to "assess", "evaluate",
  "research", "profile", "do a target review", "druggability", "feasibility", or
  "background" on it — and especially for GPCRs, ion channels, kinases, and other
  receptors where agonist vs antagonist matters, even if they don't say "CADD". If
  the target name is ambiguous, this skill runs a short intake Q&A first.
---

# CADD Target Assessment

Generate a structured, evidence-backed assessment of a protein target for a
computer-aided drug design program, delivered as one self-contained HTML file.

The job has three parts: **(A)** confirm exactly which protein the user means,
**(B)** gather hard data with the bundled scripts and read the relevant structural
papers, **(C)** synthesize it into the HTML report from the template. Hard numbers
come from scripts so they're reproducible; the judgment (mechanism, pros/cons,
feasibility) is yours, grounded in the papers.

## Environment

The scripts need **Python 3.9+ with Biopython** installed, plus network access to
the UniProt and RCSB REST APIs. Everything else is the Python standard library —
no `requests`, no other third-party packages. Install the one dependency once:

```
pip install biopython          # or: pip install -r requirements.txt
```

Then run each script with that interpreter:

```
python <script> ...
```

> Author's local setup: the `PDBAnalysis` conda env already has Biopython, so the
> author runs `/data/jinxin.liu/miniconda3/envs/PDBAnalysis/bin/python <script>`.
> Substitute whatever Python has Biopython on your machine.

All scripts print JSON to stdout and accept `--out <file>` to also save it. Save
the JSON into a working directory so the data is available while you write the
report. Scripts are in `scripts/`, the report template in `assets/`.

## Step A — Intake Q&A (do this before researching)

Target names are frequently ambiguous: "H1" could be histamine H1 receptor, a
histone, or a HERG splice form; "NHE", "PR", "CB1", and many gene aliases collide
across families. **Before committing to a protein, ask the user 2–4 short
questions** so you assess the right one. Cover:

- **Indication / therapeutic area** (e.g. "allergy", "oncology") — disambiguates
  fast and tells you the relevant tissue/biology.
- **Any other names / the gene symbol** if they know it.
- **Species of interest** if not human (default is human).
- **Program intent** if useful — agonist vs antagonist goal, modality.

Then resolve candidates and confirm before the full run:

```
python scripts/protein_info.py --query "<name>" --list
```

This lists matching reviewed proteins (accession, gene, organism). If more than
one plausibly matches, show them and confirm with the user. If it's unambiguous
(a clean gene symbol like `HRH1`), you can proceed and just state which UniProt
accession you locked onto.

Skip the questions only when the user already gave an unambiguous accession/gene
plus context. When in doubt, ask — a wrong protein wastes the whole report.

## Step B — Gather data

Once the protein is fixed (you have its UniProt accession, e.g. `P35367`):

**1. Protein facts**
```
python scripts/protein_info.py --acc P35367 --out work/protein.json
```
Gene, synonyms, length, mass, subcellular location, family, function text, PDB
count, and the sequence.

**2. Cross-species conservation** (mouse, rat, dog, monkey)
```
python scripts/species_identity.py --ref-acc P35367 --out work/species.json
```
Override the species set with `--species "mouse:10090,rat:10116,dog:9615,monkey:9544"`
if a different model organism matters.

**3. Within-family identity** — you supply the family/paralog human gene symbols.
You know the major families (GPCR subfamilies, kinase groups, channel families);
when unsure, web-search "<target> family members" or check the UniProt family note
first, then pass the genes:
```
python scripts/family_identity.py --ref-acc P35367 --members "HRH2,HRH3,HRH4" --out work/family.json
```
Include the closest relatives that pose selectivity risk. For a receptor with
cross-family look-alikes (e.g. H1 vs muscarinic/serotonin receptors), add those
genes too — the report's selectivity discussion depends on it.

**4. PDB structures**
```
python scripts/pdb_structures.py --acc P35367 --out work/pdb.json
```
Every experimental entry: method, resolution, modeled positions (UniProt
numbering), bound ligand(s) with buffer/cryo additives filtered out, and the
primary-citation title/journal/year (DOI/PMID when RCSB has them). The structure
title usually names the real drug (e.g. comp_id `Y5E` ↔ "mepyramine").

**5. Read the key papers.** The structure citations are your route into the
mechanism. Use `WebSearch`/`WebFetch` on the citation titles (and DOIs) to learn:
- For each ligand-bound structure, **is the ligand an agonist, antagonist, inverse
  agonist, or biased agonist?** This is the single most important pharmacology
  fact for GPCRs, ion channels, and nuclear receptors — capture it per structure.
- **How does the small molecule regulate the protein?** The activation/inactivation
  switch, key anchoring residues, conserved vs target-specific subpockets,
  conformational changes (TM movements, allosteric coupling). Describe the
  mechanism the paper actually proposes, with the paper as the source.
Prioritize the highest-resolution and most recent structures, and any review that
covers the activation mechanism. You don't need to read all papers for a target
with 40 structures — cover the apo/inactive, agonist-bound, and antagonist-bound
representatives.

## Step C — Write the HTML report

Copy `assets/report_template.html` to the output path (default
`<TARGET>_assessment.html` in the working directory) and fill it in. The template
has inline comments at every section explaining what goes where. Rules that keep
the report trustworthy:

- **Never hand-edit the computed numbers.** Identity %, lengths, masses, resolutions,
  and positions come from the JSON — transcribe them faithfully. If a script
  returned `null`/`note`, say "not available", don't invent.
- **Attribute the narrative.** Mechanism, agonist/antagonist calls, and pocket
  observations must trace to a paper you read; link or cite it. Distinguish what
  the structure/data shows from your inference.
- **The agonist/antagonist column is the point** for receptors and channels. Use the
  modulator pills (`agonist` / `antagonist` / `inverse` / `neutral`). If the
  literature doesn't specify, mark it `neutral`/unknown rather than guessing.
- **The Mechanism section (template §5)** is required for targets with directional
  pharmacology (GPCRs, ion channels, nuclear receptors). For a plain enzyme with no
  agonist/antagonist switch, retitle it "Mechanism of small-molecule regulation"
  and describe active-site/allosteric inhibition instead, or remove it and note why.
- **The final assessment (§6)** must follow the requested shape: a protein recap, a
  homology recap, a structure recap, a two-column **pros/cons** from the
  computational standpoint, and an explicit **feasibility** verdict (High / Moderate
  / Low) with a one-paragraph rationale and concrete suggested next steps. Anchor
  pros/cons in the data: e.g. abundant apo+holo structures and clear pockets are
  pros; high within-family identity (selectivity risk), only low-resolution EM, or
  no apo structure are cons.

Then tell the user the report path and give a 3–4 sentence verbal summary
(what the target is, structural readiness, the feasibility call).

## What "good" looks like

The report should let a computational chemist decide, in five minutes, whether the
target is worth a structure-based campaign: what it is and where it lives, whether
animal models translate, how hard selectivity will be, what structures and ligands
exist, how agonists differ from antagonists mechanistically, and a defensible
feasibility verdict with next steps — every hard number reproducible from the
scripts, every claim about mechanism traceable to a paper.
