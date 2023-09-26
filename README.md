## EBP-Nor Genome Assembly pipeline

This repository contains the EBP-Nor genome assembly pipeline. This pipeline is implemented in snakemake.
The current version of this pipeline is designed to work on a linux cluster with slurm as workload manager.

## Requirements & Setup

Some software need to be configured/installed before the pipeline can be run

### Conda setup

Most required software, including snakemake itself, can be installed using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

Once conda is installed, you can create a new environment containing most necessary software from the provided asm_pipeline.yaml file as follows:

```shell
conda create -n asm_pipeline --file=asm_pipeline.yaml
```

### Other software setup

The following software need to be installed manually:

- KMC v3.1.1 (https://github.com/tbenavi1/KMC)
- HiFiAdapterFilt (https://github.com/sheinasim/HiFiAdapterFilt)
- NCBI FCS-Adaptor (https://github.com/ncbi/fcs/wiki/FCS-adaptor)
- NCBI FCS-GX (https://github.com/ncbi/fcs/wiki/FCS-GX)

### BUSCO database setup

As in general, computing nodes are not connected to the internet, BUSCO lineage datasets need to be downloaded manually before running the pipeline.
This can easily be done by running

```shell
busco --download eukaryota
```

You will need to specify the folder where you downloaded the busco lineages in the config file (see Parameter section).

### Data

In its current setup, the pipeline requires both PacBio HiFi data and paired-end Hi-C data.
You need to create a directory containing the sequencing data of one sample/species with following structure:
The directory contains a subdirectory named "genomic_data", which contains two subdirectories:
1) A directory called "pacbio" containing one or multiple PacBio HiFi read files (*each ending in .fastq.gz*)
2) A directory called "hic" containing one pair of Hi-C Illumina reads (*ending in _1.fastq.gz and _2.fastq.gz*)

The absolute path of this directory needs to be specified in the asm_params file (see Parameter section)

### Parameters

The necessary config files for running the pipeline can be found in the Asm_profile folder.

General snakemake and cluster submission parameters are defined in Asm_profile/config.yaml, 
software-specific parameters are defined in Asm_profile/asm_params.yaml.

For software not installed by conda, the installation path needs to be provided to the Snakemake pipeline by editing following parameters in the asm_params.yaml:

- Set the "adapterfilt_install_dir" parameter to the installation path of HiFiAdapterFilt
- Set the "KMC_path" parameter to the installation path of KMC
- Set the "fcs_path" parameter to the location of the ```run_fcsadaptor.sh``` and ```fcs.py``` scripts
- Set the "fcs_adaptor_image" and "fcs_gx_image" parameters to the paths to the ```fcs-adaptor.sif``` and ```fcs-gx.sif``` files respectively
- Set the "fcs_gx_db" parameter to the path of the fcs-gx database

A couple of other parameters need to be verified as well in the asm_params.yaml before running the pipeline:

- The location of the input data (```species_dir```) should be set to the folder containing the input data.
- The location of the busco lineages (```busco_db_dir```) should be set to the folder containing the busco lineages files downloaded earlier
- The required BUSCO lineage for running the BUSCO analysis needs to set (```busco_lineage``` parameter)
- The NCBI taxid of your species, required for the decontamination step (```taxid``` parameter)

## Usage and run modes

Before running, make sure to activate the conda environment containing the necessary software.
To run the pipeline run the following command from the directory containing the "Snakefile":

```
snakemake --profile Asm_profile {run_mode}
```

The ```--profile``` option must link to the folder containing the main config.yaml file (default ```Asm_profile```).
If the parameter yaml file (```asm_params.yaml```) is not in the profile folder, you can link to it by adding ```--configfile {path/to/asm_params.yaml}```
All parameters in the ```asm_params.yaml``` can also be set/changed on the command line by adding ```--config {parameter}=new_value```

The pipeline contains different run_modes, and the run mode should always be the last argument on the command line:

- "all" (default): will run the full workflow including pre-assembly (genomescope & smudgeplot), assembly, scaffolding, and decontamination
- "pre_assembly": will run only the pre-assembly steps (genomescope & smudgeplot)
- "assembly": will filter the HiFi reads and assemble them using hifiasm (also using the Hi-C reads), and run busco
- "scaffolding": will run all steps necessary for scaffolding (filtering, assembly, HiC filtering, scaffolding, busco), but without pre-assembly
- "decontamination": will run assembly, scaffolding, and decontamination, but without pre-assembly and busco analyses

## Output

All generated output will be present in the "results" directory, which will be created in the specified input directory containing the raw data.
This results directory contains different subdirectories related to the different steps in the assembly:
- results/pre_assembly: genomescope and smudgeplot output (each in its own subfolder)
- resulst/assembly: Hifiasm assembly output and corresponding busco results
- results/scaffolding: scaffolding output, separated in two folders:
  - meryl: meryl databases used for filtering HiC reads
  - yahs: scaffolding output, including final scaffolds and their corresponding busco results
- results/decontamination: decontamination output of the final scaffolded assembly
