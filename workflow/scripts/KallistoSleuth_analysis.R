## KallistoSleuth_analysis.R

path.home = "~/Documents/Github/KallistoSleuth/"   # path to pipeline

source( paste0(path.home,"workflow/scripts/functions_kallisto_sleuth.R"))
### Paths and files
path.to.results = paste0(path.home,"results/kallisto/")  # path to kallisto results
path.to.metadata = paste0(path.home,"config/")   # path to configuration file

file.samples = paste0(path.to.metadata,"samples.tsv") 
file.units = paste0(path.to.metadata,"units.tsv")

include.batch.effects = TRUE # if you have a column for batch effects in samples.tsv
qval.cutoff = 0.01  # qvalue will be used to determine significant enrichment
number.of.top.genes = 40

save_eps = FALSE # will save images to the figs/ folder if set to TRUE.

### samples and metadata -----
sample_id <- dir(path.to.results) # list sample IDs
sample_id <- sample_id[grep(sample_id, pattern = "idx$",invert = T)]

s2c <- read.table(file.samples, header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::mutate(s2c, path = paste0(path.to.results,sample_id))


### Load gene names assoicated with each transcript ID -------
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "mmusculus_gene_ensembl")
ttg <- biomaRt::getBM(
  attributes = c("ensembl_transcript_id", "transcript_version",
                 "ensembl_gene_id", "external_gene_name", "description",
                 "transcript_biotype"),
  mart = mart)
ttg <- dplyr::rename(ttg, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
ttg <- dplyr::select(ttg, c('target_id', 'ens_gene', 'ext_gene'))

### build sleuth object --------
so <- sleuth_prep(s2c,target_mapping = ttg,
                  aggregation_column = 'ens_gene', 
                  extra_bootstrap_summary = TRUE)

# will counter for batch effects if designated in the header 
if(include.batch.effects){
  so <- sleuth_fit(so, ~condition+batch_effect, 'full')
  so <- sleuth_fit(so, ~batch_effect, 'reduced')
}else{
  so <- sleuth_fit(so, ~condition, 'full')
  so <- sleuth_fit(so, ~1, 'reduced')
}

# perform sleuth tests

models(so) # list models here

so <- sleuth_lrt(so, 'reduced', 'full') # will give the likelihood ratio test
so <- sleuth_wt(so, which_beta = 'conditionlate') # will give the wald test when one of the conditions is specified

tests(so)  # should show the tests you just performed

# write summary tables of tests 
sleuth_table_wt_late_tx <- sleuth_results(so, test_type = "wt", test = 'conditionlate', show_all = TRUE, pval_aggregate = FALSE)
sleuth_table_lrt_tx <- sleuth_results(so, test_type = "lrt", test = 'reduced:full', show_all = TRUE, pval_aggregate = FALSE)


### Plot a PCA of samples based on gene expression ----
sleuth::plot_pca(so, color_by = 'condition')+theme_bw()

sleuth::plot_pc_variance(so,pca_number = 10)

### Plot a comparison of gene expression level vs variability ----
sleuth::plot_mean_var(so, which_model = "full", 
                      point_alpha = 0.4, point_size = 2, 
                      point_colors = c("black", "dodgerblue"),
                      smooth_alpha = 1, 
                      smooth_size = 0.75, 
                      smooth_color = "red")+theme_bw()


### QC: Plot a distribution of reads by different groupings to look for count bias ----
plot_group_density(so, use_filtered = TRUE, units = "est_counts",
                   trans = "log", grouping = "condition", offset = 1)

plot_group_density(so, use_filtered = TRUE, units = "est_counts",
                   trans = "log", grouping = "batch_effect", offset = 1)

### Plot a heatmap of gene expression ----
heatmap = sleuth::plot_transcript_heatmap(so, transcripts = sleuth_table_lrt_tx$target_id[1:number.of.top.genes],
                                          color_high = "orange", color_low = "purple",
                                          main="Top differentially expressed transcripts (TPM)")
dev.off()
plot(heatmap)

# Unfortunately, the plot_trnascript_heatmap does not allow you to add in gene names. 
# We wrote a geom_tile ggplot_object that allows you to do that (see functions_kallisto_sleuth.R)

generate.heatmap(sleuthObject = so, sleuth_table = sleuth_table_lrt_tx, topNgens=number.of.top.genes)

### Plot individual counts as for particular genes for various conditions -----
transcript.of.interest = "ENSMUST00000209735.1"
plot_bootstrap(so, transcript.of.interest, units = "est_counts", color_by = "condition")


### Plot a volcano plot and label top genes.

sleuth_results_oe <- sleuth_results(so, 
                                    test = 'conditionlate', 
                                    show_all = TRUE, pval_aggregate = FALSE)

generate.volcano(sleuthObject = sleuth_results_oe, qvalueCutoff = .01, pvalueCutoff = 10**-5, betaCutoff = 2, topNgenes = number.of.top.genes)
