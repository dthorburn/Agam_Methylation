# Nanopore Methylation Pipeline
Pipelines developed for methylation calls using nanopore reads primarily for the Crisanti Lab using the *Anopheles gambiae* PEST (AgamP4) reference genome, but this pipeline can easily be applied to other organisms. 

## Bismark-Like 
I've developed the pipeline using [nextflow](https://www.nextflow.io/) version 20.10.0 (the most recent version on Imperial HPC; date: 25/01/22). This pipeline only uses ONT-Guppy to call modified bases, and then uses [modbam2bed](https://github.com/epi2me-labs/modbam2bed) to generate an output similar to [Bismakr](https://www.bioinformatics.babraham.ac.uk/projects/bismark/). 

Steps to run this pipeline:

  1. Obtain the ONT-Guppy container from me. If you are not in the Crisanti Lab, build the GPU version of [ONT-Guppy](https://community.nanoporetech.com/protocols/Guppy-protocol/v/gpb_2003_v1_revab_14dec2018/linux-guppy) into a singularity container: `ONT_Guppy_GPU.sif`.
  2. Download a copy of the appropraite reference genome (i.e., *Anopheles gambiae* genome from [vectorbase](https://vectorbase.org/vectorbase/app/record/dataset/DS_2251b21396#description)). 
  3. Clone this repository into your project directory using `git clone`.
  4. Update the project directory path and add required (and optional) arguments to the `Methylation.sh` PBS job submission script. 
  5. Run the pipeline using `qsub Nextflow_Submit.sh`. 

If you are not using a PBS job submission system, you'll need to update all the scripts to reflect the job submission system. See [Nextflow doucmentation](https://www.nextflow.io/docs/latest/executor.html) for help. 

**Important:** for some reason all files have to be in live, and not ephemeral on Imperial HPC. Unsure why, but nextflow fails to load the data correctly and will throw errors if you have the guppy input there.  

The help message from the nextflow script is below:
```
Usage:
  This pipeline was developed to call modified bases in nanopore fast5 data using ONT-Guppy. 

  To use follow these steps:
  1. Update project directory path in Methylation.sh 
  2. Add required (and optional) arguments listed below
  3. Submit pipeline coordinator using qsub Methylation.sh
  
  If you require available HPC jobs for alternative scripts lower job concurrency options. 

  Required arguments:
    --RefGen                                       Path to reference genome. Usage: '--RefGen /path/to/genome.fasta'
    --Container                                    Path to ONT-Guppy singularity container.
    --Input                                        Path to nanopore tarball. 
    --RunID                                        Name of sample. 
  
  Optional arguments:
    --help                                         Show this message

  Concurrency arguments:                           Imperial HPC only permits 50 jobs per user. These options limit the 
                                                   number of concurrent processes running per step. NB. Multiple 
                                                   processes can be running at the same time.
    --GP_Forks                                     Guppy forks. Default: 20
    --BI_Forks                                     Index forks. Default: 20

  Debugging arguments:
    --Skip_Decompress                              Skip decompression of tarball
    --Skip_Guppy                                   Skip calling modified bases with ONT-Guppy
    --Skip_Index                                   Skip indexing output bams
    --Skip_Processing                              Skip processing output bams 
```
