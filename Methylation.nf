#!/usr/bin/env nextflow

// Pipeline developed for methylation calls of Oxford Nanopore reads. 
// Author: Miles Thorburn <d.thorburn@imperial.ac.uk>
// Date last modified: 25/01/2022

def helpMessage() {
  log.info """
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
        """
}

println "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nNanopore Methylation Pipeline v0.1\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
// println "${PWD}, ${HOME}, ${PATH}"


// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Set up directory structure - I wanted to include downloading of useful files, but I can't get wget to work for ref-gen and cant publish the sing' container.
// This is redundant as the publishDir argument will create a directory, but this idea can be used to also download data or set up environments if needed. 
if (params.init) {
    process Init {
      executor = 'local'
      script:
      """
      mkdir $params.ProjectDir/00_Resources
      mkdir $params.ProjectDir/01_Raw_Input
      mkdir $params.ProjectDir/02_Guppy_Output
      mkdir $params.ProjectDir/02_Guppy_Output/01_ModBases
      mkdir $params.ProjectDir/02_Guppy_Output/02_Methylation
      """
    }
}

params.publishDir = './02_Guppy_Output'

// Creating an input channel with the path to the input directory and adding an error message if it's empty. 
Channel
  .fromPath("${params.P1_inDir}/")
  .ifEmpty { error "Cannot find input folder ${params.P1_inDir}" }
  .set { P1_in }

// Channel
//   .fromPath("${params.refGen}", checkIfExists: true)
//   .set { ref_genome }

// The first Guppy process to call modified bases against the reference genome. If statement used to skip steps if jobs crash or don't complete. 
if (params.GupModCall == 1) {
    guppy_config = file( params.P1_gupConf )
    ref_genome = file( params.refGen )

    process GupcallBases { 
      errorStrategy { 'retry' }
      maxRetries 3
      
      // publishDir creates a symlink from temp folders (/ProjectDir/work/) to the output path so places all output files there. Can be modified to only accept specific files using regex
      // i.e., publishDir "params.P1_outDir", mode: 'copy', overwrite: false, pattern: "*.bam"
      publishDir(
        path: "${params.publishDir}/01_ModBases",
        mode: 'copy',
        saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
      )
      
      // Options for submitting to cluster GPU
      // Lower numbers of runners:cpus needed as jobs kept crashing due to too many threads being spawned
      executor = 'pbspro'
      clusterOptions = "-lselect=1:ncpus=${params.P1_threads}:mem=${params.P1_memory}:ngpus=1:gpu_type=${params.P1_GPU} -lwalltime=${params.P1_walltime}:00:00"
     
      // beforeScript 'module load cuda/11.4.2; export SINGULARITY_NV=1'
      // container "${params.Gup_container}"
     
      input:
      path input_dir from P1_in
      path ref_genome
     
      output:
      // Because a wildcard was used, all output files are treated as a single emission from this process and need to be modified later. 
      // Output channels can only be used as input for one other process, so putting all the bams into 2 channels permits use by two processes. 
      path "outdir/pass/*.bam" into all_bams1, all_bams2
              
      script:
      """
      mkdir outdir
      module load cuda/11.4.2
      singularity exec --nv $params.Gup_container \\
      guppy_basecaller \\
        --config "${params.P1_gupConf}" \\
				--device "cuda:0" \\
				--bam_out \\
				--recursive \\
				--compress \\
				--align_ref "${ref_genome}" \\
				-i "${input_dir}" \\
				-s outdir \\
				--gpu_runners_per_device "${params.P1_GPU_runners}" \\
				--num_callers "${params.P1_callers}"
      """
    }
}

if (params.GupIndexBam == 1) {
    // This changes the single emission of all the files into individual files to be processed.
    bams_ch = all_bams2
                .flatMap()  

    process IndexBams {
      // Submitting locally
      executor = 'local'
      
      publishDir "${params.publishDir}/01_ModBases", mode: 'copy'
      
      input: 
      path queryBam from bams_ch
      
      output:
      path "*.bai" into bais_ch
        
      script:
      """
      module load samtools/1.3.1
      samtools index ${queryBam}
      """
    }
}

// This stops ProcBams from starting until all of the processes from IndexBams is complete -- THIS DIDN'T WORK. UNSURE WHY. TRYING .COLLECT() ABOVE INSTEAD. 
//all_bais = bais_ch
//              .flatten()
//              .view()
                
if (params.GupModProc == 1) {
  

  
  ref_genome2 = file( params.refGen )
  
  process ProcBams {
    executor = 'pbspro'
    clusterOptions = "-lselect=1:ncpus=${params.P3_threads}:mem=${params.P3_memory} -lwalltime=${params.P3_walltime}:00:00"
    // Taken from docs, but each failed attempt will cause the sleep to get longer. 
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries 3
    
    publishDir "${params.publishDir}/02_Methylation", mode: 'copy'
    
    input:
    // Needs the declaration block to be all the bams at the same time. 
    path bams from all_bams1
    path bais from bais_ch.collect()
    
    output:
    path "Guppy.cpg.bam"
    
    script:
    """
    module load anaconda3/personal
    source activate modbam2bed
    modbam2bed -e -m ${params.P3_modType} --cpg -t ${params.P3_threads} ${ref_genome2} ${bams} > Guppy.cpg.bam
    """
  }
}
