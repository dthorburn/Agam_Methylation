#!/usr/bin/env nextflow

/*
 *  Pipeline developed for methylation calls of Oxford Nanopore reads. 
 *  Author: Miles Thorburn <d.thorburn@imperial.ac.uk>
 *  Date last modified: 10/03/2022
 */

                                                            // ========================================================
                                                            // Setting help and error messages
                                                            // ========================================================
def helpMessage() {
  log.info """
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
==============================================================================================================================
        """
}
// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

if(!params.RefGen) {
  log.info"""
ERROR: No reference genome path provided! --RefGen /path/to/genome.fasta
==============================================================================================================================
  """
  helpMessage()
  exit 0
}

if(!params.Container) {
  log.info"""
ERROR: No Guppy Container provided! --Container /path/to/Guppy.sif
==============================================================================================================================
  """
  helpMessage()
  exit 0
}

if(!params.Input) {
  log.info"""
ERROR: No input tarball provided! --Input /path/to/ONT-fast5.tar.gz
==============================================================================================================================
  """
  helpMessage()
  exit 0
}

log.info """
==============================================================================================================================
                                              Nanopore Methylation Pipeline v2
==============================================================================================================================

Run ID                : ${params.RunID}
Reference             : ${params.RefGen}
Guppy Container       : ${params.Container}
Input Directory       : ${params.Input}
Results               : ${PWD}/03_Results

==============================================================================================================================
"""
                                                            // ========================================================
                                                            // Setting the value channels (can be read unlimited times)
                                                            // ========================================================
Gup_container = file( params.Container, checkIfExists: true )
run_id = Channel.value( params.RunID )
ref_genome = file( params.RefGen, checkIfExists: true )
ref_name   = ref_genome.getBaseName()
ref_dir    = ref_genome.getParent()
//ref_index  = file( "${ref_dir}/${ref_name}.*" ).flatMap().first()
ref_fai  = file( "${ref_dir}/${ref_name}*.fai", checkIfExists: true ).first()
                                                            // ========================================================
                                                            // Step 1: Decompress tarball
                                                            // ========================================================
if ( params.Skip_Decompress == false ) {
  Channel
    .fromPath("${params.Input}")
    .ifEmpty { error "${params.Input} file not found!" }
    .set { tarball }

  process Decompress {
    errorStrategy { 'retry' }
    maxRetries 3

    executor = 'pbspro'
    clusterOptions = "-lselect=1:ncpus=${params.DC_threads}:mem=${params.DC_memory}gb -lwalltime=${params.DC_walltime}:00:00"

    tag { run_id }

    publishDir(
      path: "${params.TarDir}/",
      mode: 'move',
      saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )

    input:
    path(tb) from tarball
    val run_id

    output:
    path("**.fast5") into fast5s
    //path("**.fast5")
    //path("bad_practises.txt") into fast5s
    //path("${params.TarDir}") into fast5_dir

    script:
    """
    tar -xf ${tb}
    ##echo "" > bad_practises.txt
    """
  }
}
                                                            // ========================================================
                                                            // Step 2: ONT-Guppy Call
                                                            // ========================================================
// The first Guppy process to call modified bases against the reference genome. If statement used to skip steps if jobs crash or don't complete. 
if (params.Skip_Guppy == false) {
  if ( params.Skip_Decompress ){
    Channel
      .fromPath("${params.TarDir}/*.fast5")
      .ifEmpty { error "Cannot find fast5s in ${params.TarDir}" }
      .collect()
      .set { fast5s }
  }
  Channel
    .fromPath("${params.TarDir}")
    .set { fast5_dir }

  process Guppy_Call { 
    errorStrategy { 'retry' }
    maxRetries 3
    maxForks params.GP_Forks
    
    //tag { chunk }
    tag { run_id }
    
    // publishDir creates a symlink from temp folders (/ProjectDir/work/) to the output path so places all output files there. Can be modified to only accept specific files using regex
    // i.e., publishDir "params.Guppy_outDir", mode: 'copy', overwrite: false, pattern: "*.bam"
    publishDir(
      path: "${params.GuppyDir}",
      mode: 'copy',
      //saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) },
    )
    
    // Options for submitting to cluster GPU
    // Lower numbers of runners:cpus needed as jobs kept crashing due to too many threads being spawned
    executor = 'pbspro'
    clusterOptions = "-lselect=1:ncpus=${params.Guppy_threads}:mem=${params.Guppy_memory}gb:ngpus=${params.Guppy_GPUs}:gpu_type=${params.Guppy_GPU} -lwalltime=${params.Guppy_walltime}:00:00"
   
    beforeScript 'module load cuda/11.4.2'
    
    input:
    path ref_genome
    //path ref_index
    path Gup_container
    path fast5_files from fast5s
    path file_in from fast5_dir
   
    output:
    // Because a wildcard was used, all output files are treated as a single emission from this process and need to be modified later. 
    // Output channels can only be used as input for one other process, so putting all the bams into 2 channels permits use by two processes. 
    // path "outdir/pass/*.bam" into all_bams1, all_bams2
    path("*.bam") into all_bams1, all_bams2
            
    script:
    """
    if [ ! -d outdir ]; then mkdir outdir; fi
    n_slots=`expr ${params.Guppy_threads} / 2`
    taskset -c 0-\${n_slots} \\
      singularity exec --nv ${Gup_container} \\
        guppy_basecaller \\
          --config "${params.Guppy_gupConf}" \\
          --device cuda:all \\
          --bam_out \\
          --recursive \\
          --compress \\
          --align_ref ${ref_genome} \\
          -i "${file_in}" \\
          -s outdir \\
          --gpu_runners_per_device "${params.Guppy_GPU_runners}" \\
          --num_callers "${params.Guppy_callers}"
      
    mv outdir/pass/*.bam ./
    """
  }
}
                                                            // ========================================================
                                                            // Step 3: Index Bams
                                                            // ========================================================
if (params.Skip_Index == false) {
  if ( params.Skip_Guppy ){
    Channel
      .fromPath("${params.GuppyDir}/*.bam")
      .ifEmpty { error "No bams found in ${params.GuppyDir}" }
      .set { all_bams1 }
  }
  // This changes the single emission of all the files into individual files to be processed.
  bams_ch = all_bams1
              .map { file -> tuple(file.baseName.replaceAll("^.*_", "").replaceAll(".fast5", ""), file) }
              .flatMap()  

  process IndexBams {
    errorStrategy { 'retry' }
    maxRetries 3
    maxForks params.BI_Forks
  
    executor = 'pbspro'
    clusterOptions = "-lselect=1:ncpus=${params.Index_threads}:mem=${params.Index_memory}gb -lwalltime=${params.Index_walltime}:00:00"
    
    tag { chunk }

    publishDir(
      path: "${params.GuppyDir}",
      mode: 'copy',
    )
    
    beforeScript 'module load samtools/1.2'

    input: 
    set chunk, path(queryBam) from bams_ch
    
    output:
    path("*.bai") into bais_ch
      
    script:
    """
    samtools index ${queryBam}
    """
  }
}
                                                            // ========================================================
                                                            // Step 4: Process modified bases
                                                            // ========================================================
if (params.Skip_Processing == false) {
  if (params.Skip_Guppy) {
    Channel
      .fromPath("${params.GuppyDir}/*.bam")
      .ifEmpty { error "No bams found in ${params.GuppyDir}" }
      .set { all_bams2 }
  }
  if (params.Skip_Index) {
    Channel
      .fromPath("${params.GuppyDir}/*.bai")
      .ifEmpty { error "No bais found in ${params.GuppyDir}" }
      .set { bams_ch }
  }

  process ProcBams {
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries 3
    
    executor = 'pbspro'
    clusterOptions = "-lselect=1:ncpus=${params.BP_threads}:mem=${params.BP_memory} -lwalltime=${params.BP_walltime}:00:00"
  
    tag { run_id }

    publishDir(
      path: "${params.ResDir}",
      mode: 'copy',
    )
    
    beforeScript 'module load anaconda3/personal; source activate modbam2bed'

    input:
    // Needs the declaration block to be all the bams at the same time. 
    path bams from all_bams2
    path bais from bais_ch.collect()
    path ref_genome
    val run_id

    output:
    path "${run_id}.cpg.bam"
    
    script:
    """
    modbam2bed -e -m ${params.P3_modType} --cpg -t ${params.P3_threads} ${ref_genome} ${bams} > ${run_id}.cpg.bam
    """
  }
}