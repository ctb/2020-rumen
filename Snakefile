SAMPLES = ['ERR2027899']

rule all:
    input:
        "outputs/trim/ERR2027899_R1.trim.fq.gz",
        expand("outputs/megahit/tig{tig}.contigs.fa", tig=['32', '52', '85'])

rule adapter_trim:
    input:
        r1 = "inputs/raw/{sample}_1.fastq.gz",
        r2 = 'inputs/raw/{sample}_2.fastq.gz',
        adapters = 'inputs/adapters.fa'
    output:
        r1 = 'outputs/trim/{sample}_R1.trim.fq.gz',
        r2 = 'outputs/trim/{sample}_R2.trim.fq.gz',
        o1 = 'outputs/trim/{sample}_o1.trim.fq.gz',
        o2 = 'outputs/trim/{sample}_o2.trim.fq.gz'
    conda: 'env.yml'
    shell:'''
     trimmomatic PE {input.r1} {input.r2} \
             {output.r1} {output.o1} {output.r2} {output.o2} \
             ILLUMINACLIP:{input.adapters}:2:0:15 MINLEN:25  \
             LEADING:2 TRAILING:2 SLIDINGWINDOW:4:2
    '''

rule kmer_trim_reads:
    input: 
        "outputs/trim/{sample}_R1.trim.fq.gz", 
        "outputs/trim/{sample}_R2.trim.fq.gz"
    output: "outputs/abundtrim/{sample}.abundtrim.fq.gz"
    conda: 'env.yml'
    shell:'''
    interleave-reads.py {input} | 
        trim-low-abund.py -C 3 -Z 18 -M 30e9 -V - -o {output}
    '''

rule index_queries:
     input: "queries/{tig}.fa"
     output:
          expand("queries/{{tig}}.fa.{ext}",
                 ext=['amb', 'ann', 'bwt', 'pac', 'sa'])
     conda: 'env.yml'
     shell: '''
          bwa index {input}
     '''

rule map_nbhd_reads:
     input:
          query = 'queries/{tig}.fa',
          index = 'queries/{tig}.fa.bwt',
          reads = 'rumen_sgc_k31_r1_search_oh0/{tig}.fa.cdbg_ids.reads.fa.gz'
     output:
          'outputs/nbhd_map/{tig}_nbhd_reads.sam'
     conda: 'env.yml'
     threads: 4
     shell: '''
          bwa mem -t {threads} {input.query} {input.reads} > {output}
     '''

rule sam_unmapped_reads:
    input: "outputs/nbhd_map/{sample}_nbhd_reads.sam"
    output: "outputs/nbhd_map/{sample}_unmapped.sam"
    conda: 'env.yml'
    shell:'''
         samtools view -f 4 {input} > {output}
    '''

rule unmapped_reads_to_fasta:
    input: "outputs/nbhd_map/{sample}_unmapped.sam"
    output: "outputs/nbhd_map/{sample}_unmapped.fa"
    conda: 'env.yml'
    shell:'''
         samtools fasta {input} > {output}
    '''

rule megahit_unmapped_reads:
    input: "outputs/nbhd_map/{sample}_unmapped.fa"
    output: "outputs/megahit/{sample}.contigs.fa"
    conda: 'env.yml'
    shell:'''
         megahit -r {input} --min-contig-len 142 \
             --out-dir {wildcards.sample}_megahit \
             --out-prefix {wildcards.sample}
    mv {wildcards.sample}_megahit/{wildcards.sample}.contigs.fa {output}
    rm -rf {wildcards.sample}_megahit
    '''
