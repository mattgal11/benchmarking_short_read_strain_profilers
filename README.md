# Benchmarking strain-level profiling of _Escherichia coli_ in short-read gut metagenomes (Galbraith _et al_., 2026)

## Summary

This repository contains the scripts from the Data Summary used for the benchmarking of short-read strain-level profiling tools (PanTax, PathoScope, StrainGE, Strainify, StrainR2 and StrainScan) on _Escherichia coli_. Reproducibility requirements for each tool are also provided. Supplementary data is available on the associated FigShare link from the paper ().

## Structure

- `envs/` - human-readable conda environment files used to define each tool’s dependencies. This includes InsilicoSeq, used for simulating metagenomes.
- `locks/` - fully reproducible conda-lock files for exact environment recreation.
  - Example usage:
    - `conda create -n <env_name> --file locks/<env_name>-linux-64.lock`
    - `conda activate <env_name>`
  - Notes:
    -  Strainify (v1.1.0) was installed from the upstream repository (https://github.com/treangenlab/Strainify). A `conda-lock` file could not be generated due to dependency resolution issues
    -  PanTax requires `Gurobi` (v11), which is not included in the lockfile because it requires a user-specific license.
- `scripts/`
  - `InSilicoSeq_simulated_metagenomes/` - contains scripts and associated files (including abundance files) for simulating metagenomes using InSilicoSeq.
  - `R_visualizations_and_stats/` - contains scripts to recreate R visualizations and any associated statistical analysis used for figures in the manuscript. R scripts are grouped by results narrative - e.g., Fig. 1 and Fig. S1 are presented together. See FigShare repository for the required input data files.
  - `SLURM_jobs_simulated_metagenomes/` - contains per-tool scripts from the simulated metagenomes for job submissions on SLURM clusters.
    - PanTax has a script included but this is not a SLURM job script, rather a file-driven batch processing loop.
    - PathoScope has two scripts, one per core module used (MAP & ID), which should be run sequentially.
    - StrainGE has an additional extraction/normalisation script to process the results.
