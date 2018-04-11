#global settings
shell.prefix("set -o pipefail; ")
shell.prefix("set -e; ")
shell.prefix("set -u; ")
localrules: all
configfile: "config.json"

rule all:
        input:
                expand("00.rawdata/{library}/{library}.foo", library=config["Library"]),
                expand("00.rawdata/{library}/{library}.fastq", library=config["Library"]),
                expand("01.Pre_cutadapt/{library}/{library}.foo2", library=config["Library"]),
                expand("01.Pre_cutadapt/{library}/{library}.trimmed.cutAdapt3.fastq", library=config["Library"]),
                expand("02.hsa_rRNA_bowtie2/{library}/{library}.rRNA.sam",library=config["Library"]),
                expand("02.hsa_rRNA_bowtie2/{library}/{library}.no_rRNA.fq",library=config["Library"]),
                expand("03.ordered_mapping_bowtie2/{library}/00.hg38/{library}.hg38.sam", library=config["Library"]),
                expand("03.ordered_mapping_bowtie2/{library}/00.hg38/{library}.hg38.unAligned.fq", library=config["Library"]),
                expand("03.ordered_mapping_bowtie2/{library}/01.miRNA/{library}.miRNA.unAligned.fq", library=config["Library"]),
                expand("03.ordered_mapping_bowtie2/{library}/01.miRNA/{library}.miRNA.sam", library=config["Library"]),
                expand("03.ordered_mapping_bowtie2/{library}/02.piRNA/{library}.piRNA.unAligned.fq", library=config["Library"]),
                expand("03.ordered_mapping_bowtie2/{library}/02.piRNA/{library}.piRNA.sam", library=config["Library"]),
                expand("03.ordered_mapping_bowtie2/{library}/03.Y_RNA/{library}.Y_RNA.unAligned.fq", library=config["Library"]),
                expand("03.ordered_mapping_bowtie2/{library}/03.Y_RNA/{library}.Y_RNA.sam", library=config["Library"]),
                expand("03.ordered_mapping_bowtie2/{library}/04.snRNA/{library}.snRNA.unAligned.fq", library=config["Library"]),
                expand("03.ordered_mapping_bowtie2/{library}/04.snRNA/{library}.snRNA.sam", library=config["Library"]),
                expand("03.ordered_mapping_bowtie2/{library}/05.srpRNA/{library}.srpRNA.unAligned.fq", library=config["Library"]),
                expand("03.ordered_mapping_bowtie2/{library}/05.srpRNA/{library}.srpRNA.sam", library=config["Library"]),
                expand("03.ordered_mapping_bowtie2/{library}/06.tRNA/{library}.tRNA.unAligned.fq", library=config["Library"]),
                expand("03.ordered_mapping_bowtie2/{library}/06.tRNA/{library}.tRNA.sam", library=config["Library"]),
                expand("03.ordered_mapping_bowtie2/{library}/07.other_lncRNA/{library}.other_lncRNA.unAligned.fq",library=config["Library"]),
                expand("03.ordered_mapping_bowtie2/{library}/07.other_lncRNA/{library}.other_lncRNA.sam",library=config["Library"]),
                expand("03.ordered_mapping_bowtie2/{library}/08.mRNA/{library}.mRNA.unAligned.fq",library=config["Library"]),
                expand("03.ordered_mapping_bowtie2/{library}/08.mRNA/{library}.mRNA.sam",library=config["Library"])

rule fastqc:
        input:
                fastq = "00.rawdata/{library}/{library}.fastq"
        output:
                dir = "00.rawdata/{library}/",
                foo = "00.rawdata/{library}/{library}.foo"
        params:
                jobname = "{library}.fastqc"
        shell:
                """
                fastqc -o {output.dir} {input.fastq}
                """

rule cut_adapter:
        input:
                fastq = "00.rawdata/{library}/{library}.fastq"
        output:
                rmTooShort = "01.Pre_cutadapt/{library}/{library}.trimmed.cutAdapt3.fastq",
                TooShort = "01.Pre_cutadapt/{library}/{library}.trimmed.cutAdapt3.TooShort.fastq"
        params:
                quality_left = config["Cutadapter"]["quality_left"],
                quality_right = config["Cutadapter"]["quality_right"],
                adapter3 = lambda wildcards: config['Sample_adapter3'][wildcards.library],
                miniL = config["Cutadapter"]["miniL"],
                jobname = "{library}.cutadapt"
        log:
                log_cutAdapt3 = "01.Pre_cutadapt/{library}/{library}.cutAdapt3.log"
        shell:
                """
                ~/apps/anaconda2/bin/python ~/local/bin/cutadapt -u -100 -q {params.quality_left},{params.quality_right} -a {params.adapter3} --trim-n  -m {params.miniL} --too-short-output={output.TooShort}  -o {output.rmTooShort} {input.fastq}  >{log.log_cutAdapt3}
                """

rule fastqc2:
        input:
                fastq = "01.Pre_cutadapt/{library}/{library}.trimmed.cutAdapt3.fastq"
        output:
                dir = "01.Pre_cutadapt/{library}/",
                foo2 = "01.Pre_cutadapt/{library}/{library}.foo2"
        params:
                jobname = "{library}.fastqc2"
        shell:
                """
                fastqc -o {output.dir} {input.fastq}
                """

                rule bowtie2_remove_rRNA:
        input:
                fastq = "01.Pre_cutadapt/{library}/{library}.trimmed.cutAdapt3.fastq"
        output:
                no_rRNA = "02.hsa_rRNA_bowtie2/{library}/{library}.no_rRNA.fq",
                sam = "02.hsa_rRNA_bowtie2/{library}/{library}.rRNA.sam"
        params:
                jobname = "{library}.bowtie2_remove_rRNA",
                index = config["References"]["bowtie2_hg38_rRNA_index_dir"],
                cpu = config["Remove_rRNA"]["cpu"]
        shell:
                """
                bowtie2 -p {params.cpu} --sensitive-local --un {output.no_rRNA} -x {params.index} {input.fastq} -S {output.sam}
                """

rule bowtie2_mapping_hg38:
        input:
                fastq = "02.hsa_rRNA_bowtie2/{library}/{library}.no_rRNA.fq"
        output:
                unAligned = "03.ordered_mapping_bowtie2/{library}/00.hg38/{library}.hg38.unAligned.fq",
                sam = "03.ordered_mapping_bowtie2/{library}/00.hg38/{library}.hg38.sam"
        params:
                jobname = "{library}.bowtie2_mapping_hg38",
                index = config["References"]["bowtie2_hg38_genome_index_dir"],
                tmp = "03.ordered_mapping_bowtie2/{library}/00.hg38/",
                cpu = config["Remove_rRNA"]["cpu"]
        shell:
                """
                bowtie2 -p {params.cpu} --sensitive-local --un {output.unAligned} -x {params.index} {input.fastq} -S {output.sam}
                """
                
rule bowtie2_mapping_miRNA:
        input:
                fastq = "02.hsa_rRNA_bowtie2/{library}/{library}.no_rRNA.fq"
        output:
                unAligned = "03.ordered_mapping_bowtie2/{library}/01.miRNA/{library}.miRNA.unAligned.fq",
                sam = "03.ordered_mapping_bowtie2/{library}/01.miRNA/{library}.miRNA.sam"
        params:
                jobname = "{library}.bowtie2_mapping_miRNA",
                index = config["References"]["bowtie2_hg38_miRNA_index_dir"],
                tmp = "03.ordered_mapping_bowtie2/{library}/01.miRNA/",
                cpu = config["Remove_rRNA"]["cpu"]
        shell:
                """
                bowtie2 -p {params.cpu} --sensitive-local --un {output.unAligned} -x {params.index} {input.fastq} -S {output.sam}
                """
                
                rule bowtie2_mapping_piRNA:
        input:
                fastq = "03.ordered_mapping_bowtie2/{library}/01.miRNA/{library}.miRNA.unAligned.fq"
        output:
                unAligned = "03.ordered_mapping_bowtie2/{library}/02.piRNA/{library}.piRNA.unAligned.fq",
                sam = "03.ordered_mapping_bowtie2/{library}/02.piRNA/{library}.piRNA.sam"
        params:
                jobname = "{library}.bowtie2_mapping_piRNA",
                index = config["References"]["bowtie2_hg38_piRNA_index_dir"],
                tmp = "03.ordered_mapping_bowtie2/{library}/02.piRNA/",
                cpu = config["Remove_rRNA"]["cpu"]
        shell:
                """
                bowtie2 -p {params.cpu} --sensitive-local --un {output.unAligned} -x {params.index} {input.fastq} -S {output.sam}
                """

rule bowtie2_mapping_Y_RNA:
        input:
                fastq = "03.ordered_mapping_bowtie2/{library}/02.piRNA/{library}.piRNA.unAligned.fq"
        output:
                unAligned = "03.ordered_mapping_bowtie2/{library}/03.Y_RNA/{library}.Y_RNA.unAligned.fq",
                sam = "03.ordered_mapping_bowtie2/{library}/03.Y_RNA/{library}.Y_RNA.sam"
        params:
                jobname = "{library}.bowtie2_mapping_Y_RNA",
                index = config["References"]["bowtie2_hg38_Y_RNA_index_dir"],
                tmp = "03.ordered_mapping_bowtie2/{library}/03.Y_RNA/",
                cpu = config["Remove_rRNA"]["cpu"]
        shell:
                """
                bowtie2 -p {params.cpu} --sensitive-local --un {output.unAligned} -x {params.index} {input.fastq} -S {output.sam}
                """

rule bowtie2_mapping_snRNA:
        input:
                fastq = "03.ordered_mapping_bowtie2/{library}/03.Y_RNA/{library}.Y_RNA.unAligned.fq"
        output:
                unAligned = "03.ordered_mapping_bowtie2/{library}/04.snRNA/{library}.snRNA.unAligned.fq",
                sam = "03.ordered_mapping_bowtie2/{library}/04.snRNA/{library}.snRNA.sam"
        params:
                jobname = "{library}.bowtie2_mapping_snRNA",
                index = config["References"]["bowtie2_hg38_snRNA_index_dir"],
                tmp = "03.ordered_mapping_bowtie2/{library}/04.snRNA/",
                cpu = config["Remove_rRNA"]["cpu"]
        shell:
                """
                bowtie2 -p {params.cpu} --sensitive-local --un {output.unAligned} -x {params.index} {input.fastq} -S {output.sam}
                """

                rule bowtie2_mapping_srpRNA:
        input:
                fastq = "03.ordered_mapping_bowtie2/{library}/04.snRNA/{library}.snRNA.unAligned.fq"
        output:
                unAligned = "03.ordered_mapping_bowtie2/{library}/05.srpRNA/{library}.srpRNA.unAligned.fq",
                sam = "03.ordered_mapping_bowtie2/{library}/05.srpRNA/{library}.srpRNA.sam"
        params:
                jobname = "{library}.bowtie2_mapping_srpRNA",
                index = config["References"]["bowtie2_hg38_srpRNA_index_dir"],
                tmp = "03.ordered_mapping_bowtie2/{library}/05.srpRNA/",
                cpu = config["Remove_rRNA"]["cpu"]
        shell:
                """
                bowtie2 -p {params.cpu} --sensitive-local --un {output.unAligned} -x {params.index} {input.fastq} -S {output.sam}
                """

rule bowtie2_mapping_tRNA:
        input:
                fastq = "03.ordered_mapping_bowtie2/{library}/05.srpRNA/{library}.srpRNA.unAligned.fq"
        output:
                unAligned = "03.ordered_mapping_bowtie2/{library}/06.tRNA/{library}.tRNA.unAligned.fq",
                sam = "03.ordered_mapping_bowtie2/{library}/06.tRNA/{library}.tRNA.sam"
        params:
                jobname = "{library}.bowtie2_mapping_tRNA",
                index = config["References"]["bowtie2_hg38_tRNA_index_dir"],
                tmp = "03.ordered_mapping_bowtie2/{library}/06.tRNA/",
                cpu = config["Remove_rRNA"]["cpu"]
        shell:
                """
                bowtie2 -p {params.cpu} --sensitive-local --un {output.unAligned} -x {params.index} {input.fastq} -S {output.sam}
                """

rule bowtie2_mapping_other_lncRNA:
        input:
                fastq = "03.ordered_mapping_bowtie2/{library}/06.tRNA/{library}.tRNA.unAligned.fq"
        output:
                unAligned = "03.ordered_mapping_bowtie2/{library}/07.other_lncRNA/{library}.other_lncRNA.unAligned.fq",
                sam = "03.ordered_mapping_bowtie2/{library}/07.other_lncRNA/{library}.other_lncRNA.sam"
        params:
                jobname = "{library}.bowtie2_mapping_other_lncRNA",
                index = config["References"]["bowtie2_hg38_other_lncRNA_All_index_dir"],
                tmp = "03.ordered_mapping_bowtie2/{library}/07.other_lncRNA/",
                cpu = config["Remove_rRNA"]["cpu"]
        shell:
                """
                bowtie2 -p {params.cpu} --sensitive-local --un {output.unAligned} -x {params.index} {input.fastq} -S {output.sam}
                """
                
                rule bowtie2_mapping_mRNA:
        input:
                fastq = "03.ordered_mapping_bowtie2/{library}/07.other_lncRNA/{library}.other_lncRNA.unAligned.fq"
        output:
                unAligned = "03.ordered_mapping_bowtie2/{library}/08.mRNA/{library}.mRNA.unAligned.fq",
                sam = "03.ordered_mapping_bowtie2/{library}/08.mRNA/{library}.mRNA.sam"
        params:
                jobname = "{library}.bowtie2_mapping_mRNA",
                index = config["References"]["bowtie2_hg38_mRNA_index_dir"],
                tmp = "03.ordered_mapping_bowtie2/{library}/08.mRNA/",
                cpu = config["Remove_rRNA"]["cpu"]
        shell:
                """
                bowtie2 -p {params.cpu} --sensitive-local --un {output.unAligned} -x {params.index} {input.fastq} -S {output.sam}
                """
