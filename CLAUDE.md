# Pipeline Context for Claude Code

This repo is a MATLAB/Python/R data-processing pipeline for
neural recordings from NeuroNexus silicon probes in mice. This file is loaded automatically at the
start of every Claude Code session in this repo.

## Domain terminology
- **LFP** - Local field potential
- **TF** - Time-frequency (from Morlet Wavelet transform)
- **CSD** - Current source density
- **spike** - neural action potentials detected by threshold crossings in highpassed raw data
- **MUA** - multi-unit activity (another term for spikes)
- **stim_lfp_stimchunks** — peri-stimulus lfp data from each channel. Usually 1kHz sampling rate
- **stim_spike_stimchunks** — peri-stimulus MUA data from each channel. Usually 30kHz sampling rate
  to an event
- **ITPC** - Inter-trial phase clustering. Measure of frequency band phase alignment between trials at each time point 
- **Spike times** — timestamps of detected neural spiking events
- **ProbeInfo** — MATLAB struct of relevent experiment information used in downstream analyses


## Repository structure & pipeline order

1. `PreProcessing_StimData/` — raw data → trial-aligned/epoched data and general experiment info, plus plots for visualization
2. The pipeline then forks into two independent analysis branches, each
   answering a different follow-up question — either or both may be run
   on the preprocessed output, and neither depends on the other:
   - `Functional_Connectivity/` — connectivity measures such as ITPC per animal/session. 
   - `Summarize_MultiAnimal/` — cross-animal/condition aggregation and summary. This data is output as csv to plot with R

## Data contracts: header comment convention

MATLAB has no type system, so header comments are the primary way script
compatibility is tracked across the pipeline. Some scripts already follow
this convention; others don't yet — when editing or writing a script,
bring its header up to this standard:

```matlab
% SCRIPT_NAME.m
%
% Description: <what this script does>
%
% Inputs:
%   var_name (type, shape) - description, units if applicable
%     e.g. lfp_epochs (struct array, 1xN trials) - fields:
%       .data (double, channels x samples) - LFP signal, in microvolts
%       .trial_id (int) - trial identifier
%
% Outputs:
%   var_name (type, shape) - description, units if applicable
%
% Dependencies: <scripts/functions this expects to have run first, if any>
```

When asked to edit a script, check whether its header matches this
template; if it's missing or incomplete, update it as part of the edit
rather than leaving it stale.

## Working conventions

- Preserve existing variable/struct field naming across scripts unless
  explicitly asked to refactor — downstream scripts depend on exact field
  names.
- Flag any inferred data contract (a shape/field assumption not stated in
  a header comment) rather than silently assuming it's correct.
