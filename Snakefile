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
                expand("01.Pre_cutadapt/{library}/{library}.trimmed.cutAdapt3.fastq", library=config["Library"])

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
