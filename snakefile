
 # Snakemake Pipeline for Variant calling


rule trim:
   output:
     fq_out1    = "trimmed_fastq/{sample}_1.trimmed.fastq",
     fq_out1un  = "trimmed_fastq/{sample}_1un.trimmed.fastq",
     fq_out2    = "trimmed_fastq/{sample}_2.trimmed.fastq",
     fq_out2un  = "trimmed_fastq/{sample}_2un.trimmed.fastq"
   input:
     fq_in1 = "untrimmed_fastq/{sample}_1.fastq.gz",
     fq_in2 = "untrimmed_fastq/{sample}_2.fastq.gz"
   shell:
     """
       trimmomatic PE -threads 1 \
         {input.fq_in1} {input.fq_in2} \
         {output.fq_out1} {output.fq_out1un} \
         {output.fq_out2} {output.fq_out2un} \
         SLIDINGWINDOW:4:20 \
         MINLEN:25 \
         ILLUMINACLIP:untrimmed_fastq/NexteraPE-PE.fa:2:40:15
     """
