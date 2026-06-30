# cadd-target-assessment

A [Claude Code](https://code.claude.com) skill that turns a protein **target name** into a
single self-contained **HTML dossier** for computer-aided drug design (CADD):

- **Protein basics** — gene, length, molecular weight, subcellular location, biological function
- **Cross-species conservation** — % identity to mouse / rat / dog / monkey orthologs
- **Within-family identity** — % identity/similarity to the closest paralogs (selectivity landscape)
- **PDB structure catalog** — every experimental structure with bound ligands, modeled positions, and (for receptors/channels) **agonist vs antagonist** annotation
- **Mechanism** — how small molecules / binders modulate or inhibit the target, from the primary literature
- **Computational verdict** — pros/cons and an explicit **feasibility** call with suggested next steps

Hard numbers (identity, masses, structure tables) are computed live from **UniProt** and **RCSB PDB** so they're reproducible; the mechanism and feasibility narrative is synthesized from cited papers. If a target name is ambiguous, the skill first runs a short intake Q&A.

---

## Requirements

- **Claude Code** (the skill runs inside it).
- **Python 3.9+ with [Biopython](https://biopython.org/)** on the machine — the scripts use it for sequence alignment. Everything else is the Python standard library (no `requests`), and the scripts need outbound network access to the UniProt and RCSB REST APIs.

Install the one dependency:

```bash
pip install -r skills/cadd-target-assessment/requirements.txt   # just biopython
```

(Conda users: `conda create -n cadd python=3.11 biopython -c conda-forge` works too — then point the skill at that interpreter.)

---

## Install

### Option A — as a plugin (one command)

```text
/plugin install https://github.com/<you>/cadd-target-assessment-plugin
```

That's it — Claude Code discovers the skill under `skills/`. (You can also publish via a
[plugin marketplace](https://code.claude.com/docs/en/plugin-marketplaces) and have users
`/plugin marketplace add <you>/<repo>` then `/plugin install cadd-target-assessment@<marketplace>`.)

### Option B — manual copy (no plugin)

Claude Code discovers skills from `~/.claude/skills/` (personal) and `<project>/.claude/skills/` (per-project). Just copy the skill folder into either:

```bash
git clone https://github.com/<you>/cadd-target-assessment-plugin
cp -r cadd-target-assessment-plugin/skills/cadd-target-assessment ~/.claude/skills/
```

Then restart Claude Code so it picks up the new skill.

---

## Usage

In Claude Code, name a target and ask for an assessment, e.g.:

> *"Do a CADD target assessment for KRAS G12C — I'm scoping a covalent inhibitor program."*

or invoke it explicitly:

> `/cadd-target-assessment TL1A — indication RA, want an inhibitor`

The skill resolves the protein (asking to disambiguate if needed), runs the data scripts, reads the key structural papers, and writes `<TARGET>_assessment.html`.

---

## What's inside

```
skills/cadd-target-assessment/
├── SKILL.md                 # workflow Claude follows
├── requirements.txt         # biopython
├── scripts/                 # data collection (UniProt + RCSB via stdlib urllib + Biopython)
│   ├── cadd_common.py        # HTTP + alignment helpers
│   ├── protein_info.py       # gene / length / mass / location / family / function
│   ├── species_identity.py   # mouse/rat/dog/monkey ortholog identity
│   ├── family_identity.py    # paralog identity / similarity
│   └── pdb_structures.py     # PDB structures + ligands + citations
└── assets/
    └── report_template.html  # single-page HTML report template
```

## Notes & caveats

- **Identity method:** global Needleman–Wunsch alignment (BLOSUM62), identity computed over aligned non-gap columns. Numbers are reproducible but may differ from other tools that use local alignment or different gap penalties.
- The scripts hit public REST APIs and have no API key requirement; they retry on transient failures.
- Mechanism and feasibility text is model-synthesized from cited literature — **verify key claims against the linked papers** before acting on them.

## License

[MIT](LICENSE) © 2026 Jinxin
