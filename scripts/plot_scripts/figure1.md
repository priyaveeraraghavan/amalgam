Figure1_Analysis
================
Priya Veeraraghavan
2023-03-25

``` r
library(rjson)
input <- "C:/Users/prv892/Documents/GitHub/amalgam/metadata.json"

file_path_ref <- fromJSON(file = input)

# import the relevant parameters and references
source(file.path(file_path_ref$project_directory, 
                 file_path_ref$sequencing$setup_script))
```

    ## Warning: package 'dplyr' was built under R version 4.2.3

``` r
source(file.path(file_path_ref$project_directory,
                 file_path_ref$plots$setup_script))
```

# Main figure

``` r
plot_main_1d <- function(plot_savefile) {
  load(file.path(file_path_ref$project_directory,
                 file_path_ref$sequencing$results$all_dpc23_comparison))
  
  pca_res <- get_pcs(assay(vst_subset))

  plt <- ggplot(pca_res$pcs, 
              aes(PC1, PC2, 
                  color=subtype, shape=location)) + 
         geom_point(size=0.5) +
         xlab(paste0("PC1: ", pca_res$variances[1], "% Variance Explained", collapse="")) +
         ylab(paste0("PC2: ", pca_res$variances[2], "% Variance Explained", collapse="")) +
         theme_classic() + 
    scale_color_manual(values=c(cpn_color, cthpn_color)) +
    theme(legend.position="None",
          text=element_text(size=6))
  
  plot(plt)
  
  
  ggsave(plot_savefile, plt, device="pdf", useDingbats=FALSE, width=46, height=33, units="mm")
  
}

plot_main_1d(file.path(file_path_ref$project_directory,
                       file_path_ref$plots$plot_directory, "main_1d.pdf"))
```

![](figure1_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

The genes chosen for labeling have previously been described as being
localized and/or with localization elements already determined

``` r
plot_main_1e <- function(plot_savefile) {
  load(file.path(file_path_ref$project_directory,
                              file_path_ref$sequencing$results$gc_vs_soma_dpc23_comparison))
  
  res_dge_valid <- subset(res_dge, is_valid_gc & baseMean > 100) 
  
  res_dge_gc_enriched <- subset(res_dge_valid, 
                                log2FoldChange > 0 & padj < 0.05)
  
  res_dge_soma_enriched <- subset(res_dge_valid, 
                                  log2FoldChange < 0 & padj < 0.05)
  
  res_dge_tolabel <- subset(res_dge_valid, 
                            external_gene_name %in% c("Actb", "Nrn1", "Ybx1", "Creb1") |
                              log2FoldChange > 3)
  

  plt <- rasterize(ggplot(res_dge_valid, 
         aes(baseMean, log2FoldChange)) +
    geom_point(color=background_color,
               size=0.5),
    dpi=300) +
    rasterize(geom_point(data=res_dge_gc_enriched,
                         mapping=aes(baseMean, log2FoldChange),
                         color=accent_color_dark,
                         size=0.5),
              dpi=300) +
      rasterize(geom_point(data=res_dge_soma_enriched,
                         mapping=aes(baseMean, log2FoldChange),
                         color=background_color_dark,
                         size=0.5),
              dpi=300) +
    geom_text_repel(data=res_dge_tolabel,
               mapping=aes(baseMean, log2FoldChange, 
                           label=external_gene_name),
               color="black", size=1) +
    scale_x_log10() +
    theme_classic() +
    ylab("Log2 Fold Change GC/Soma") + 
    xlab("Mean Normalized Counts") +
    theme(text=element_text(size=6))
  
  plot(plt)
  
  ggsave(plot_savefile, plot = plt, 
         width = 53, height = 53, units="mm", 
         useDingbats=FALSE)
  
  
  
  
  
}

plot_main_1e(file.path(file_path_ref$project_directory,
                       file_path_ref$plots$plot_directory, "main_1e.pdf"))
```

![](figure1_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
plot_main_1f <- function(plot_savefile) {
  
  load(file.path(file_path_ref$project_directory,
                              file_path_ref$sequencing$results$gc_vs_soma_dpc23_comparison))
  
  go_grouped_simple <- simplify(go_grouped)
  goterm_df <- subset(go_grouped_simple@compareClusterResult, Count >= 10 ) %>%
    rowwise() %>%
    mutate(gene_ratio=gene_ratio_compute(GeneRatio),
           description=str_wrap(Description, width=20)) %>%
    arrange(desc(Cluster), desc(gene_ratio))
  
  goterm_df$description <- factor(goterm_df$description, 
                                  levels=goterm_df$description)
  
  plt <- ggplot(goterm_df, aes(str_wrap(gsub("_", " ", Cluster), width=8), description, size=gene_ratio, color=p.adjust)) +
    geom_point() +
    scale_color_gradient(low=accent_color_scale_dark, high=accent_color_scale_light)+
    theme_bw() +
    xlab("") + ylab("") +
    theme(text=element_text(size=8))
  
  plot(plt)
  
  ggsave(plot_savefile, plot=plt, device="pdf",
         width=29*3, height=61*3, units="mm", useDingbats=FALSE)
  
}

plot_main_1f(file.path(file_path_ref$project_directory,
                       file_path_ref$plots$plot_directory, "main_1f.pdf"))
```

![](figure1_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->