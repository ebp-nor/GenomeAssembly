## EBP-Nor Genome Assembly pipeline

This repository contains the EBP-Nor genome assembly pipeline. This pipeline is implemented in snakemake.
This pipeline is developed to create haplotype-resolved genome assemblies from long reads (PacBio HiFi or Oxford Nanopore) and optionally using HiC data.
The pipeline is primarly designed for diploid eukaryotic organisms and is implemented to work on a linux cluster with slurm as workload manager.

## Requirements & Setup

Some software need to be configured/installed before the pipeline can be run.
This workflow is set up to work on a computing cluster with SLURM.
The Snakemake workflow will handle job submission of the different steps of the pipeline.

### Conda setup

Most required software, including Snakemake itself, can be installed using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

You can create an environment containing most necessary software from the provided asm_pipeline.yaml file as follows:

```shell
conda create -n asm_pipeline --file=worfklow/envs/asm_pipeline.yaml
```

### Other software setup

One software (smudgeplot) has to be installed in a separate environment. To install this environment, run:

```shell
conda create -c bioconda -c conda-forge -n smudge0.4 smudgeplot fastk
```

The pipeline will automatically detect the environment, and activate it when needed. If you decide to rename this enviroment, please renamen the "smudgeplot_conda_env" paremeter in the config file.

The following software need to be installed manually:

- HiFiAdapterFilt (https://github.com/sheinasim/HiFiAdapterFilt)
- Oatk (https://github.com/c-zhou/oatk)
- OatkDB (https://github.com/c-zhou/OatkDB)
- NCBI FCS-Adaptor (https://github.com/ncbi/fcs/wiki/FCS-adaptor)
- NCBI FCS-GX (https://github.com/ncbi/fcs/wiki/FCS-GX)

Optional for the first step of assembly curation:
- Functional singularity installation
- Rapid Curation Singularity images (```.sif```) (https://gitlab.com/wtsi-grit/rapid-curation/-/tree/65fa15125dbcce025d02471f30fc093170f48bf7) (Versions after 2024-03-04 do not have the required singularity images)

Please refer to their respective installation instructions to properly install them. You will need to privide the installation paths of these software to the config file (see Parameter section).

### BUSCO database setup

As in general, computing nodes are not connected to the internet, BUSCO lineage datasets need to be downloaded manually before running the pipeline.
This can easily be done by running

```shell
busco --download eukaryota
```

You will need to specify the folder where you downloaded the busco lineages in the config file (see Parameter section).

### Data

This pipeline is created for using PacBio HiFi or Oxford Nanopore reads, preferably together with paired-end HiC data (but HiC data can be omitted).
At least one set of long reads should be provided to run the pipeline. All input data should be in ```fastq.gz``` format. 

Which data will be used for the assembly of the genome depends on the input data:

- If only HiFi or ONT are provided, only those will be used.
- If both are provided, ONT will be used if the ONT reads are derived from the r10 technology
- Otherwise, HiFi will be used with the ONT reads as support. 

You will need to specify the absolute paths of all input files in the config file (see Parameters section).

### Parameters

The necessary template config files for running the pipeline can be found in the config folder.

General snakemake and cluster submission parameters are defined in ```config/config.yaml```, 
data- and software-specfic parameters are defined in ```config/asm_params.yaml```.

First, define the paths of the input files you want to use:
- pacbio: path to the location of the PacBio HiFi reads (```.fastq.gz```)
- hicF and hicR: path to the forward and reverse HiC reads respectively (both in ```.fastq.gz```)
- ont: path to the location of ONT reads (```.fastq.gz```)
- r10: Set to True if the ONT reads are r10, otherwise set to False

For software not installed by conda, the installation path needs to be provided to the Snakemake pipeline by editing following parameters in the ```config/asm_params.yaml```:

- Set the "adapterfilt_install_dir" parameter to the installation path of HiFiAdapterFilt
- Set the "oatk_dir" parameter to the installation path of oatk
- Set the "oatk_db" parameter to the directory where you downloaded the oatk_db files
- Set the "fcs_path" parameter to the location of the ```run_fcsadaptor.sh``` and ```fcs.py``` scripts
- Set the "fcs_adaptor_image" and "fcs_gx_image" parameters to the paths to the ```fcs-adaptor.sif``` and ```fcs-gx.sif``` files respectively
- Set the "fcs_gx_db" parameter to the path of the fcs-gx database
- Set the "curation_install_dir" parameter to the path containing the rapid-curation singularity images (```.sif```) if performing curation

A couple of other parameters need to be specified in addition, as they will differen from input to input:

- The location of the downloaded busco lineages (```busco_db_dir```) should be set to the folder containing the busco lineages files downloaded earlier
- The required BUSCO lineage for running the BUSCO analysis needs to set (```busco_lineage``` parameter). Run ```busco --list-datasets``` to get an overview of all available datasets.
- The required oatk lineage for running organelle genome assembly (```oatk_lineage``` parameter). Check https://github.com/c-zhou/OatkDB for an overview of available lineages.
- A boolean value wether the species is plant (for plastid prediction) or not (```oatk_isPlant```; set to either ```True``` or ```False```, without quotation marks)
- The NCBI taxid of your species, required for the decontamination step (```taxid``` parameter)
- The telomere repeat sequence (```curation_telomere_rep```) if performing curation

## Usage and run modes

Before running, make sure to activate the conda environment containing the necessary software: ```conda activate asm_assembly```.
To run the pipeline, run the following command:

```
snakemake --profile config/ --configfile config/asm_params.yaml --snakefile workflow/Snakefile {run_mode}
```

If you invoke the snakemake command in another directory than the one containing the ```workflow``` and ```config``` folders, 
or if the config files (```config.yaml``` and ```asm_params.yaml```) are in another location, you need to specify their correct paths on the command line.

The workflow parameters can be modified in 3 ways:
- Directly modifying the ```config/asm_parameters.yaml``` file
- Overriding the default parameters on the command line: ```--config parameter=new_value```
- Overriding the default parameters using a different yaml file: ```--configfile path_to_parameters.yaml```

The pipeline has different runing modes, and the run mode should always be the last argument on the command line:

- "all" (default): will run the full workflow including pre-assembly (genomescope & smudgeplot), assembly, scaffolding, decontamination, pre-curation, and organelle assembly.
- "all_nocuration": will run the full workflow, but without the final curation step
- "pre_assembly": will run only the pre-assembly steps (genomescope & smudgeplot)
- "assembly": will use hifiasm to produce an assembly based on the input data, and run busco on the assembly. If HiC data is available, the resulting assembly will also be scaffolded.
- "decontamination": will run assembly, scaffolding (if HiC data is provided), and decontamination, and perform BUSCO analysis at different steps.
- "organelles": will run only organnelle genome assembly

## Output

All generated output will be located in the "results" directory, which will be created in the folder from where you invoke the snakemake command.
This results directory contains different subdirectories related to the different steps in the assembly:
- results/pre_assembly: genomescope and smudgeplot output (each in its own subfolder)
- results/assembly: Hifiasm assembly output and corresponding busco results
- results/scaffolding: scaffolding output (if HiC data is provided), separated in two folders:
  - meryl: meryl databases used for filtering HiC reads
  - yahs: scaffolding output, including final scaffolds and their corresponding busco results
- results/decontamination: decontamination output of the final (scaffolded) assembly
- results/organelles: assembled organellar genomes
- results/curation: files necessary for further manual curation of the assembly
- results/stats: basic statistics of the input data

Additionally, a text file containing all software versions will be created in the specified input directory.
The log files of the different steps in the workflow can be found in the ```logs``` directory that will be created.