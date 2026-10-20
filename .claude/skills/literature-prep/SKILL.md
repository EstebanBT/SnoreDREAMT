---
name: literature-prep
description: Prepare scientific paper PDFs for later use. Renames each PDF to FirstAuthor_Year, converts the whole paper with the docling MCP server, and writes a structured markdown summary (methods parameters, results, limitations, relevance) plus a JSON index. Use this whenever the user adds, drops or mentions PDFs in a literature/papers/references folder, asks to "process", "prepare", "ingest", "read in" or "summarize" papers, wants papers renamed by author and year, or wants Claude to have digestible notes on the literature before planning or writing. Also use it before relying on a paper's methods details when only a raw PDF exists.
---

# Literature prep

The goal is that a paper is read carefully **once** and then exists in the project as:

```
literature/
├── Schreiner_2023.pdf              # renamed original
├── index.json                      # catalog of all processed papers (machine-readable)
└── notes/
    ├── Schreiner_2023.fulltext.md  # complete docling conversion (searchable, citable)
    └── Schreiner_2023.md           # structured summary written by you
```

Later sessions read the short summary first and grep the fulltext for specifics, instead of re-parsing PDFs. So the summary has to be trustworthy and specific enough to replicate analyses from, and the full text has to be complete.

## Workflow (per PDF)

Work out which PDFs need processing: every `*.pdf` in the literature folder (default `literature/` in the project root; use the folder the user names) that is not yet listed under `pdf` in `index.json`. If the user names specific files, process just those. Already-processed papers are skipped unless the user asks to redo them; for a redo, overwrite the existing notes and keep the existing key.

### 1. Convert with docling

```
mcp__docling__convert_document_into_docling_document(source=<absolute pdf path>)  -> document_key
mcp__docling__save_docling_document(document_key)                                  -> md_file, json_file
```

`save_docling_document` writes to docling's own cache directory; you only need the `md_file` path from its result. If the docling tools are not loaded, find them with ToolSearch (`docling`) — the server can take a moment to connect. If conversion fails, tell the user which file failed and why; don't fall back to guessing the content.

### 2. Read the whole paper

Read the full `md_file` with the Read tool — all of it, in chunks if it is long. Methods sections usually sit at the end (Nature-family journals) or in supplementary sections, and the details that matter for replication (filter settings, thresholds, exclusion criteria, statistical tests) are almost never in the abstract. Skimming the first part produces summaries that look fine but are useless for planning an analysis.

From the text, identify the first author's surname and the publication year:
- Use the surname exactly as printed (the script strips accents and spaces: "Bullón Tarrasó" → `BullonTarraso`).
- Year = year of publication (online publication date for journals; for preprints/manuscripts the year on the document; for conference proceedings the conference year). If the year is truly absent, ask the user rather than inventing one.

### 3. Rename and register

Run the bundled script (standard-library Python, any `python3` works):

```bash
python3 <skill-dir>/scripts/register_paper.py \
  --pdf "<pdf path>" --author "<Surname>" --year <YYYY> \
  --docling-md "<md_file>" \
  --meta-json '{"title": "...", "authors": ["A. Surname", "..."], "journal": "...", "doi": "...", "type": "empirical|review|methods|dataset|preprint|thesis", "tags": ["..."], "one_line": "..."}'
```

It handles the deterministic parts: the `Surname_Year` key (adding `b`, `c`, … when a different paper already has that key), the rename, copying a cleaned full text to `notes/<key>.fulltext.md`, and updating `index.json`. Use `--dry-run` first if you want to see the key without touching files. It prints the final key and paths — use those, since the key may carry a suffix.

### 4. Write the structured summary

Write `notes/<key>.md` using the template below. Guidelines that make the summary worth having:

- **Copy numbers exactly** (sample sizes, frequency bands, filter orders, thresholds, window lengths, p-values, effect sizes). These notes will be used to set analysis parameters, so a paraphrased "about 1 Hz" where the paper says "0.3–1.25 Hz" is a real error. If something is ambiguous in the text, say so instead of resolving it silently.
- **Point back to the source**: tag claims with the section they come from, e.g. `(Methods › Event detection)`, so they can be checked in the fulltext quickly.
- **Separate what the paper says from your own assessment.** Your own critical observations go under "Limitations → Not stated by authors" and "Relevance".
- **Relevance**: if the project has a stated research question (CLAUDE.md, memory, or the conversation), say concretely how this paper bears on it — which methods to borrow, which results are the comparison point, which confounds it raises. If there is no project context, keep this section short and general.
- Leave a section out (or write "Not reported") rather than padding it; e.g. reviews have no sample.

```markdown
---
key: Schreiner_2023
title: "..."
authors: [Thomas Schreiner, Marit Petzka, ...]
year: 2023
journal: Nature Communications
doi: 10.1038/...
type: empirical
tags: [respiration, slow oscillations, spindles, memory reactivation]
fulltext: Schreiner_2023.fulltext.md
---

# Schreiner et al. (2023) — <short title>

## TL;DR
2–3 sentences: question, approach, main finding.

## Research question & hypotheses

## Sample & design
Participants (N, age, sex, population, exclusions), protocol, recording setup (devices, channels, sampling rate, references).

## Methods (replication-level detail)
Subsections per analysis step — preprocessing, event detection, phase/coupling estimation, statistics — with every parameter as reported. A small table works well for detection parameters. Name toolboxes and versions.

## Key results
Bullet per finding with the statistic (test, value, p, effect size, direction) and the figure it appears in.

## Figures
One line per main figure: what it shows.

## Limitations
- Stated by authors:
- Not stated by authors:

## Relevance to this project

## Follow-up references
Cited works worth reading next, with one line on why.
```

### 5. Report back

After all papers are processed, give the user a short table: original filename → new filename, title, one-line takeaway. Mention any paper where something needed a judgment call (ambiguous year, duplicate key suffix, garbled conversion such as missing equations or tables).

## Notes

- Many PDFs at once: the per-paper work is independent, so for more than ~3 papers you can hand each paper (steps 1–4) to a subagent in parallel, then write the report yourself. Run `register_paper.py` calls one at a time though, since they all update the same `index.json`.
- Docling writes `<!-- image -->` where figures were; figure content itself is not extracted, so rely on figure captions and results text for the Figures section.
- Supplementary PDFs: process them like any paper; the script gives them a suffixed key (e.g. `Schreiner_2023b`). Say in the title that it is the supplement of the main paper and cross-link the two summaries.
