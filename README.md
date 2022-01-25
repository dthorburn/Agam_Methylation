# Agam Methylation
Pipelines developed for methylation calls using nanopore reads primarily for the Crisanti Lab, but this pipeline can easily be applied to other organisms. 

## Bismark-Like 
I've developed the pipeline using [nextflow](https://www.nextflow.io/) version 20.10.0 (the most recent version on Imperial HPC; date: 15/01/22). This pipeline only uses ONT-Guppy to call modified bases, and then uses [modbam2bed](https://github.com/epi2me-labs/modbam2bed) to generate an output similar to [Bismakr](https://www.bioinformatics.babraham.ac.uk/projects/bismark/). 

Steps to run this pipeline:

  1. Obtain the ONT-Guppy container from me. If you are not in the Crisanti Lab, build the GPU version of [ONT-Guppy](https://community.nanoporetech.com/protocols/Guppy-protocol/v/gpb_2003_v1_revab_14dec2018/linux-guppy) into a singularity container: `ONT_Guppy_GPU.sif`.
  2. Download a copy of the *Anopheles gambiae* reference genome from [vectorbase](https://vectorbase.org/vectorbase/app/record/dataset/DS_2251b21396#description). Placing both into the `00_Resources` directory.
  3. Clone this repository into your project directory using `git clone`.
  4. Update the `nextflow.config` file to represent your project requirements. Required paramaters to change are flagged.
  5. Extract the `.fast5` files into the `01_Raw_Input` directory.
  6. Run the pipeline using `qsub Nextflow_Submit`. 

The help message from the nextflow scripts is below:
```
Usage:
  You'll first need to update the paths and config file to reflect your environment and ensure you are in the same directory as the scripts.
  Then, once you place the reference genome and guppy container in 00_Resources/, execute the full workflow with:
  qsub Nextflow_Submit.sh
  
  If the workflow doesn't complete, once errors are resolved, you can resume using optional arguments by updating the Nextflow_Submit.sh script.
  
Directory structure:
  /Project_dir/                                                Project Directory
    | - Nextflow_Submit.sh                                     Pipeline submission script
    | - Methylation.nf                                         Nextflow script
    | - nextflow.config                                        Nextflow config
    | - 00_Resources/                                          Subdirectory of common resources
          | - ONT_Guppy_GPU.sif
          | - VectorBase-54_AgambiaePEST_Genome.fasta
    | - 01_Raw_Input/                                          Raw data input subdirectory
    | - 02_Guppy_Output/                                       
          | - 01_ModBases/
          | - 02_Methylation/
  
Optional arguments:
  --help                                                       Show this message
  --init                                                       Set up directory structure (deprecated)
  --resume                                                     Resume from where error occured
```

## Consensus Approach (*Under development*)
I am using the consensus framework for methylation calling developed into the [METEORE](https://github.com/comprna/METEORE) tool. The nextflow pipeline will be paramaterised to run on the Imperial HPC (psbpro), and will leverage the calls from both ONT-Guppy and [ONT-Megalodon](https://github.com/nanoporetech/megalodon) to create a more robust consensus call when compared to the single caller output above.
