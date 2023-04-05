# plotting helper scripts
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrastr)
library(ggrepel)
library(ggupset)
library(pheatmap)
library(RColorBrewer)
library(rjson)
library(tidyr)

# do the setup to get color definitions, etc
input <- "C:/Users/prv892/Documents/GitHub/amalgam/metadata.json"

file_path_ref <- fromJSON(file = input)

# import the relevant parameters and references
source(file.path(file_path_ref$project_directory, 
                 file_path_ref$sequencing$setup_script))

# define colors
cpn_color_scale <- colorRampPalette(c("white", cpn_color, "black"))(5)
cpn_color_light <- cpn_color_scale[2]
cpn_color_dark <- cpn_color_scale[4]

cthpn_color_scale <- colorRampPalette(c("white", cthpn_color, "black"))(5)
cthpn_color_light <- cthpn_color_scale[2]
cthpn_color_dark <- cthpn_color_scale[4]


accent_color_scale <- colorRampPalette(c("white", accent_color_dark, "black"))(5)
accent_color_scale_light <- accent_color_scale[2]
accent_color_scale_dark <- accent_color_scale[4]

# provide the PCA based on normalized counts
# provide a data frame of rows as genes and columns as sample types
# will be the top_var most variable genes
get_pcs <- function(norm_counts, top_var=500) {
  
  norm_variances <- sort(apply(norm_counts, 1, var), decreasing=TRUE)[1:top_var]
  
  norm_subset <- data.frame(norm_counts)[names(norm_variances),]
  
  norm_pca <- prcomp(t(norm_subset))
  
  pcs <- data.frame(norm_pca$x)
  
  pcs$sample_id <- rownames(pcs)
  
  pcs <- pcs %>% inner_join(metadata, by="sample_id")
  
  # get variance explained
  variances <- round(norm_pca$sd^2/sum(norm_pca$sd^2)*100)
  
  return(list(pcs=pcs, variances=variances))
}

gene_ratio_compute <- function(string_gene_ratio) {
  gene_counts <- strsplit(string_gene_ratio, "/")[[1]]
  return(as.numeric(gene_counts[1])/as.numeric(gene_counts[2]))
}

plot_goterm_dotplot <- function(goterm_subset, count_thresh=10, ncategory=10) {
  
  goterm_df <- subset(goterm_subset, Count >= count_thresh ) %>%
    rowwise() %>%
    mutate(gene_ratio=gene_ratio_compute(GeneRatio),
           description=str_wrap(Description, width=20)) %>%
    arrange(desc(Cluster), desc(gene_ratio))
  
  descrip_levels <- c()
  for (descrip in goterm_df$description) {
    if (!(descrip %in% descrip_levels)) {
      descrip_levels <- c(descrip_levels, descrip)
    }
    
  }
  goterm_df$description <- factor(goterm_df$description, 
                                  levels=descrip_levels)
  
  plt <- ggplot(goterm_df, aes(str_wrap(gsub("_", " ", Cluster), width=8), description, size=gene_ratio, color=p.adjust)) +
    geom_point() +
    scale_color_gradient(low=accent_color_scale_dark, high=accent_color_scale_light)+
    theme_bw() +
    xlab("") + ylab("") +
    theme(text=element_text(size=8))
  
  return(plt)
  
}