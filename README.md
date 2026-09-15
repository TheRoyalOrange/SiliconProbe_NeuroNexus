# [NeuroNexus Data Pipeline] — Neuroscience Data Processing Pipeline

MATLAB/Python/R data-processing pipeline for
neural recordings from NeuroNexus silicon probes in mice, covering
stimulus-data preprocessing, functional connectivity analysis, and
cross-animal summarization of basic analysis data (LFP, MUA, TF-Power).

## Repository structure

- **`PreProcessing_StimData/`** — Preprocesses raw stimulus-aligned recording
  data (e.g. trial alignment, epoching, filtering) into a format usable by
  downstream analysis.
- **`Functional_Connectivity/`** — Computes functional connectivity measures
  (e.g. from LFP epochs / spike times) on the preprocessed data.
- **`Summarize_MultiAnimal/`** — Aggregates and summarizes results across
  animals/sessions for group-level analysis and plotting.

## Pipeline execution order

1. Run scripts in `PreProcessing_StimData/` first, to go from raw data to
   trial-aligned / epoched data structures.
2. From there, the pipeline forks into two independent analysis branches,
   run depending on the question being asked — either or both can be run
   on the preprocessed output, in any order relative to each other:
   - `Functional_Connectivity/` — computes functional connectivity
     measures per animal/session.
   - `Summarize_MultiAnimal/` — aggregates and summarizes results across
     animals/sessions.

## Data

See `data/README.md` for a description of the data structures (variable
names, shapes, units) used across the pipeline.

## Conventions

Scripts are documented with header comments describing inputs, outputs,
variable shapes, units, and struct fields — this is the primary way script
compatibility is tracked, since MATLAB has no static type system. See
`CLAUDE.md` for the header comment template used across the repo.

## Working with Claude Code

This repo includes a `CLAUDE.md` file with pipeline context that Claude
Code loads automatically. Run `claude` from the repo root to start a
session with that context in place.
