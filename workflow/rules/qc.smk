## RSEQC
import pandas as pd
unit = pd.read_csv("config/units.tsv",sep="\t")

rule fastqc:
    input:
        get_fastqs
    output:
        html="results/fastqc/{sample}-{unit}.html",
        zip="results/fastqc/{sample}-{unit}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc/{sample}-{unit}.log"
    threads: 1
    wrapper:
        "0.65.0/bio/fastqc"

rule multiqc:
    input:
        expand(["results/logs/kallisto/quant/{sample}-{unit}.log",
        "results/fastqc/{sample}-{unit}_fastqc.zip",
        "results/trimmed/{sample}-{unit}.fastq.gz"],
        sample= unit["sample"],
		unit=unit["unit"])
    output:
        "results/qc/multiqc_report.html"
    log:
        "logs/multiqc.log"
    wrapper:
        "0.31.1/bio/multiqc"
