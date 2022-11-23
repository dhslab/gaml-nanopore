version 1.0

workflow ProcessONT {

  input {
    String InDir
    String Name
    String OutDir

    String GuppyConfig
    String HG38
    String CHM13

    String JobGroup

   }

  # recall bases with 5mC (and eventually 5hmC)
  # this probably will be done at runtime eventually,
  # but not for a while I guess
  # note this step should GPU to avoid running for weeks
  call guppy_basecaller {
      input: Dir=InDir,
      Name=Name,
      Reference=Reference,
      jobgroup=JobGroup,
      queue=DragenQueue,
      docker=DragenDocker,
          }

   # we will want FASTQs as output too
   # this task should combine multiple fastqs and then
   # save them to the output directory
   call combine_fastqs {
      input: fastqs=guppy_basecaller.fastqs,
	     outdir=OutDir
	     queue=Queue,
	     jobgroup=JobGroup

   call gather_files as gather_fastqs {
       input: files=combine_fastqs.fastqs
       	      outdir=OutDir + '/fastqs'
	      queue=Queue,
	      jobgroup=JobGroup
   }

   # ideally we will take advantage of CHM13. Do both? or just CHM13?

   scatter(Reference in [HG38,CHM13]){

     # basecalling generates multiple unaligned bams
     # with base mods embedded. Need to align with
     # guppy's version of mm2 to retain the MM and ML tags
     # the aligner can return a sorted and indexed bam. 
     # it does NOT need a GPU
     scatter(ubam in guppy_basecaller.bams){
        call guppy_aligner {
           input: bam=ubam,
	 	reference=Reference,
		queue=Queue,		
      		docker='dhspence/docker-gguppy:latest',
		jobgroup=JobGroup
       }
     }

     # merge individual bams, like normal
     call merge_bams {
    	input: bams=guppy_aligner.bam,
	       name=Samplename,
	       queue=Queue,
	       jobgroup=JobGroup
     }

     # this calls (inherited) variants, phases reads and then
     # recalls and rephases the reads again. The outputs are phased bam
     # (with HP tag and base mods) and VCFs. DV is not supposed to
     # be good for somatic calls so we need to investigate other callers
     # see: https://github.com/kishwarshafin/pepper
     call pepper_margin_deepvariant {
   	 input: bam=merge_bams.bam,
	        reference=Reference,
		name=Samplename,
		queue=GPUQueue,
      		docker='kishwars/pepper_deepvariant:r0.8-gpu',
		jobgroup=JobGroup

     }

     # use ont's modbam2bed to get meth values by haplotype.
     # May need to parallelize this by genome block for speed.
     # could also write a simple python script with the modbampy library.
     call extract_5mc {
     	input: bam=pepper_margin_deepvariant.bam,
	       queue=Queue,
      	       docker='dhspene/docker-longread',
	       jobgroup=JobGroup
     }
     
     # additional tasks that are needed, but methods need to be define/decided:
     #  1) Call somatic small variants. My suggestions are freebayes (I've used on Hifi)
     #      and Octopus
     #  2) Annotate small variants w/ gnomad 3.0 (I have VCFs for this)
     #  3) Call SVs. Haley's suggestions here will be key.
     #  4) Filter potential somatic variants based on phase/haplotype information

     call gather_files as gather_results {
     	  input: files=[<pepper_margin_deepvariant.files,etc>]
	  	 outdir=OutDir + '/' + reference
     }
  }
}

task guppy_basecaller{

  input {  

    String queue
    String docker
    String jobGroup

  }
  command {

# bsub -g /dspencer/adhoc -G compute-dspencer -q general -R 'gpuhost' -gpu "num=2:gmodel=TeslaV100_SXM2_32GB:gmem=16G" -eo testGPU.err -oo testGPU.log -a "docker(dhspence/docker-gguppy)" guppy_basecaller -i ./fast5_pass/ --bam_out -x cuda:all:100% -s ./unaligned_bam -c /opt/ont/guppy/data/dna_r10.4.1_e8.2_260bps_modbases_5mc_cg_sup.cfg --num_callers 2

  }

  runtime {
      docker_image: docker
      queue: queue
      job_group: jobGroup
  }

  output {
      String outdir = "~{OutputDir}"
  }
}

