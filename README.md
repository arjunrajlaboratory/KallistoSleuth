# KallistoSleuth
Bulk RNA-seq pipeline to perform pseudoalignment and analysis using Lior Pachter's kallisto and sleuth.

## General description:
To process and analyze bulk RNA-seq data, we perform the following steps: (1) quality control on raw sequencing reads, (2) alignment of reads to a reference transcriptome, (3) conversion of aligned reads to a normalized count for each feature (such as a gene or isoform), (4) comparison of datasets to look for enrichment in gene expression.

We tackle these steps using the [Kallisto](https://pachterlab.github.io/kallisto/) + [sleuth](https://pachterlab.github.io/sleuth/) software from Lior Pacter’s group, with a snakemake workflow from Johannes Koester. We based the workflow on his Kallisto+sleuth ([link](https://github.com/snakemake-workflows/rna-seq-kallisto-sleuth)), but found that the sleuth step steps were better done directly in R.

Total estimated time (per sample - assuming ~ 1 GB single-stranded fastq against a mammalian genome ~ 1-2 hours on 1 core, can be parallelized).

## Follow the steps below to set up the pipeline:

1) Download and install [miniconda](https://docs.conda.io/en/latest/miniconda.html) (install the python 3 version) and [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).
Make sure you are using R v3.6, not all dependencies for sleuth have been upgraded for R v4.0. You can download and install older (or newer) versions of R here (link).

2) *** If you are using macOS with a version >10.11 (El Capitan, Sierra, High Sierra, Mojave, or Catalina (as of 2020)). There have been some annoying changes with how R installs some packages, including dependencies of sleuth like the data.table package. You should download and install the 3.6.z version of [r-macos-rtools](https://github.com/rmacoslib/r-macos-rtools/releases/tag/v4.0.0), which will fix that when you load packages in Rstudio.
3) Clone the git repository - https://github.com/arjunrajlaboratory/KallistoSleuth.git  to your local computer, in whatever folder you prefer.
4) Access the top of the folder, you will see empty **data/** folder and a **ref/** folders.
5) Place your zipped fastq files into the data/ folder.
Depending on your reference genome, download the appropriate tar folder from the Pachter lab (https://github.com/pachterlab/kallisto-transcriptome-indices/releases ). If you are not using a model organism you will need to create your own index. *** It is important to record the version of the transcriptome to refer to upon data publication. Your configurations file has a line to record this.
6) Move the tar file above into the ref/ folder and either double click to open the folder or open it via the command line. 
7) In the **config/**, you will see three files - **config.yaml**, **samples.tsv**, and **units.tsv**. Edit those to describe your data.

## Running the pipeline:

We run the pipeline in the terminal. 
1) Make sure you are in the correct conda environment. Often when you download snakemake it will ask you to activate a conda environment named “snakemake”.
> conda activate snakemake

If this works, on the command line you should see “(snakemake)”.

2) Run the following to make sure your pipeline is properly built. It will not run the pipeline, but assure that that your pipeline does not contain hanging files:
> snakemake --use-conda -n

3) Run your pipeline:
> snakemake --use-conda --cores 2

You can change the number of cores as needed.

4) If the pipeline works correctly you should only see green and yellow text, and the final line should read “... (100%) done)”.

5) To view a report on the metrics of each step run:
> open multiqc_report.html

This will open a MultiQC report in a browser window. You can view different metrics of your samples including the duplication rate, presence of adapters, and number of reads aligning to your reference transcriptome, as well as other things.


## Analyzing your data:

The Snakemake team has incorporated sleuth into their pipeline, but we have included it as a separate R file due to inconsistencies and errors attributed with the macOS version.

1) In Rstudio, open the file “workflow/scripts/KallistoSleuth_analysis.R”.

2) Scan the document and adjust various quantities, paths, etc.... as needed. Running Line 6 will hopefully load sleuth and all other necessary packages.

3) Process data and images for comparisons.
