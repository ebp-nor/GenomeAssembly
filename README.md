## EBP-Nor Genome Assembly pipeline

This repository contains the EBP-Nor genome assembly pipeline. This pipeline is implemented in snakemake.
The current version of this pipeline is designed to work on a cluster with slurm.

## Requirements

- A working version of Conda or Mamba
- Snakemake (can be installed through conda/mamba)
- GCC v. 11.2.0* 
- CMake v. 3.21.1 (for GCCcore-11.2.0)*
- Python v. 3.7.2 (for GCCcore-8.2.0)*
- KMC v3.1.1 (https://github.com/tbenavi1/KMC)
- HiFiAdapterFilt (https://github.com/sheinasim/HiFiAdapterFilt)
- Yahs (https://github.com/c-zhou/yahs)

*These are loaded as linux modules

## Setup

### Software

Most software installation is automatically handled by conda.
However, some software and packages need to be manually installed (see requirements).
To set these up, you need to point the pipeline to the install directories or the name of the modules.
This can be done by editing the EBPNor_GenomeAssembly_config.yaml file in the  EBPNor_GenomeAssembly_profile folder.

For modules:
- Change the parameter "apterfilt_modules" to the name of the GCC and Cmake modules. You can use a python list to add both: ["GCC-module", "CMAKE-module"]
- Change the parameter "smudge_module" to the name of the Python module 

For non-conda installed software:
- Set the "adapterfilt_install_dir" parameter to the installation path of HiFiAdapterFilt
- Set the "smudge_path" parameter to the installation path of KMC
- Set the "yahs_install_path" parameter to the installation path of Yahs

### Data

In its current setup, the pipeline requires both PacBio HiFi data and paired-end Hi-C data.
You need to create a directory containing the sequencing data of one sample/species with following structure:
The directory contains a subdirectory named "genomic_data".
The "genomic_data" directory needs to contain two subdirectories:
1) A directory called "pacbio" containing one or multiple PacBio HiFi read files (each ending in .fastq.gz)
2) A directory called "hic" containing one pair of Hi-C Illumina reads (ending in _1.fastq.gz and _2.fastq.gz)

The absolute path of this directory needs to be specified in the EBPNor_GenomeAssembly_config file (see next section)


### Parameters

General snakemake and cluster submission parameters are defined in EBPNor_GenomeAssembly_profile/config.yaml
Software-specific parameters are defined in EBPNor_GenomeAssembly_profile/EBPNor_GenomeAssembly_config.yaml

The location of the input data needs to be specified:

Define the "species_dir" parameter in the EBPNor_GenomeAssembly_profile/EBPNor_GenomeAssembly_config.yaml file as the absolute path to the input directory containing the "genomic_data" directory (see previous section).


Some other important/useful parameters to verify:
- Billing account for cluster submission ("--account" in config.yaml) (if applicable)
- Email for job updates on the cluster ("--mail_user" in config.yaml) (the line needs to be uncommented to activate email notifications)
- Assumed ploidy level of the species (for smudgeplot) ("smudge_ploidy" in EBP_GenomeAssembly_config.yaml)
- BUSCO lineage to be used when running BUSCO ("busco_lineage" in EBP_GenomeAssembly_config.yaml)

## Usage and run modes

To run the pipeline run the following command from the directory containing the "Snakefile":

```
snakemake --profile EBPNor_GenomeAssembly_profile {run_mode}
```

The pipeline contains different run_modes:
- "all" (default): will run the full workflow including pre-assembly (genomescope & smudgeplot), assembly, and scaffolding
- "pre_assembly": will run only the pre-assembly steps (genomescope & smudgeplot)
- "assembly": will filter the HiFi reads and assemble them using hifiasm (also using the Hi-C reads), and run busco
- "scaffolding": will run all steps necessary for scaffolding (filtering, assembly, HiC filtering, scaffolding, busco), but without pre-assembly

## Output

All generated output will be present in the "results" directory, which will be created in your input directory.
This results directory contains different subdirectories related to the different steps in the assembly:
- results/pre_assembly: genomescope and smudgeplot output (each in its own subfolder)
- resulst/assembly: Hifiasm assembly output and corresponding busco results
- results/scaffolding: scaffolding output, separated in two folders:
  - meryl: meryl databases used for filtering HiC reads
  - yahs: scaffolding output, including final scaffolds and their corresponding busco results
