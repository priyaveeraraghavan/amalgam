# setup script
#### perform DGE
library(tidyr)
library(ggplot2)
library(dplyr)
library(tximport)
library(DESeq2)
library(readxl)
library(pheatmap)
library(RColorBrewer)
library(clusterProfiler)
library(org.Mm.eg.db)
library(oneclust) # colorblind palette
library(ggpubr)
library(ggrepel)
library(knitr)
library(gridExtra)
library(stringr)
library(ggupset)
library(ggrastr)
library(GO.db)
library(transite)
library(rjson)

# args = commandArgs(trailingOnly=TRUE)
# input <- args[1]
input <- "C:/Users/prv892/Documents/Github/amalgam/metadata.json"

file_path_ref <- fromJSON(file = input)

project_dir <- file_path_ref$project_directory


## Default Imports
# get the pre-downloaded transcript to gene map or retrieve it from ensembl
grcm38_101 <- biomaRt::useEnsembl(biomart="ensembl", version=101, 
                                  dataset="mmusculus_gene_ensembl")
                                  
t2g_file <- file.path(project_dir, 
                      "data/external_data/ensembl_GRCm38_101_gene_biotypes.RData")

if (!file.exists(t2g_file)) {


t2g <- biomaRt::getBM(c("ensembl_transcript_id", 
                                  "ensembl_gene_id", 
                                  "external_gene_name", 
                                  "gene_biotype"),
                                mart=grcm38_101)
save(t2g, file=t2g_file)
} else {
  load(t2g_file)
}

# filter the transcripts to get relevant biotypes 

# there are discrepancies in pseudogene annotations, so also cross reference with MGI
# from http://www.informatics.jax.org/downloads/reports/index.html
mgi_biotypes <- read.csv(file.path(project_dir, 
                                   "data/external_data/MGI_ensembl_mapping.csv"), 
                         header=FALSE)
mgi_pseudogenes <- mgi_biotypes %>% subset(grepl("pseudo", V9))

# remove non poly-A RNAs
# eg histones https://www.nature.com/articles/nrg2438
non_polyA <- t2g %>% subset(grepl("^H[1234][abcf]", external_gene_name))

# filter out t2g to remove genes that are not protein coding
# also remove mitochondrially encoded genes
# also remove non poly-A containing genes
t2g_protein <- t2g %>% 
  subset(gene_biotype=="protein_coding" & 
           !grepl("mt-", external_gene_name) & 
           !(ensembl_gene_id %in% mgi_pseudogenes$V6) & 
           !(ensembl_gene_id %in% non_polyA$ensembl_gene_id))
t2g_protein$target_id <- t2g_protein$ensembl_transcript_id

# also create a gene to external gene mapping
ens2ext <- unique(t2g_protein[,c("ensembl_gene_id", "external_gene_name")])


# load metadata
metadata <- read_xlsx(file.path(project_dir, file_path_ref$sequencing$metadata)) %>% 
                        subset(!is_outlier)

salmon_dir <- file.path(project_dir, 
                        file_path_ref$sequencing$salmon_directory) 

# define colors based on CUD colorblind palette
# https://rdrr.io/cran/oneclust/man/cud.html#heading-1
accent_color_light <- cud(1)
accent_color_dark <- cud(3)

background_color <- "gray80"
background_color_dark <- "gray40"
cpn_color <- "#FF6347"
cthpn_color <- "#0A97DD"
negative_color <- "gray30"
white <- "#FFFFFF"
black <- "#000000"
sig_alpha <- 0.05

# function that gets genes associated with a particular GO term
# input: goterm -> character with go ID
# returns: character vector of ensembl gene ID associated with the goterm
get_goterm_genes <- function(goterm) {
  retrieved <- AnnotationDbi::select(org.Mm.eg.db, 
                                     keytype="GOALL", 
                                     keys=goterm, 
                                     columns="ENSEMBL")$ENSEMBL
  
  return(retrieved)
  
}

# helper function that augments a DESeq result with a column ensembl_gene_id
# input: dds_res -> output from DESeq2::results()
# returns: data.frame with contents same columns as input augmented with one
#          additional column ensembl_gene_id
row2col_ensgene <- function(dds_res) {
  
  dds_res <- data.frame(dds_res)
  dds_res$ensembl_gene_id <- rownames(dds_res)
  
  return(dds_res)
  
}

