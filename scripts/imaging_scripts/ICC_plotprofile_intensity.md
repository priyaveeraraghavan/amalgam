Proximal/Distal Immunofluoresence in CPN Axons In Vitro
================
Priya Veeraraghavan
2023-04-05

This script is designed to use the mean gray values from along an axon
in an ICC image, where there is at least one channel in addition to
DAPI. Bins the list of mean gray values and returns a proportion
relative to the maximal value.

TODO: Also look at length of neurites at P1 vs P3

``` r
library(dplyr)
```

    ## Warning: package 'dplyr' was built under R version 4.2.3

``` r
library(ggplot2)
```

``` r
project_dir <- "C:/Users/prv892/Documents/GitHub/amalgam"

plotprof_dir <- file.path(project_dir, "data/imaging/icc/plotprofile_quant")

# ideally change this to the illustrator_pdfs folder
plot_dir <-  file.path(project_dir, "data/imaging/icc/plotprofile_quant/plots")
```

``` r
get_meta_from_files <- function(rbp_id) {
  rbp_meta <- data.frame(fname=file.path(plotprof_dir, rbp_id, list.files(file.path(plotprof_dir, rbp_id)))) %>%
    rowwise() %>%
    mutate(stage=strsplit(as.character(basename(fname)), "_")[[1]][2],
           neuron=strsplit(as.character(basename(fname)), "_")[[1]][3],
           channel=gsub(".csv", "", strsplit(as.character(basename(fname)), "_")[[1]][4]))

  return(rbp_meta)
}
```

``` r
get_dapi_profiles <- function(rbp_meta) {
  
  dapi_files <- subset(rbp_meta, channel=="DAPI")$fname
  
  dapi_thresh <- rbp_meta$dapi_thresh[1]

  dapi_profiles <- data.frame(do.call(rbind, 
                     lapply(dapi_files, 
                            function(x) read.csv(x) %>% 
                              rowwise() %>% 
                              mutate(fname=basename(x)))))
  
  plot(ggplot(dapi_profiles, 
         aes(Gray_Value, color=fname)) + 
           geom_density() +
         theme(legend.position="None") +
         scale_x_log10() +
         geom_vline(xintercept = dapi_thresh))
}
```

``` r
get_channel_meangray_bin <- function(file, dapi_file, dapi_thresh, num_bins=20) {

  
  meangray_file_colnames <- c("distance_microns", "gray_value", "is_axon")  

  # read in dapi channel
  res_dapi <- read.csv(dapi_file) 
  colnames(res_dapi) <- meangray_file_colnames

  # read in other channel
  res <- read.csv(file)
  colnames(res) <- meangray_file_colnames
  
  #print(ggplot(res_dapi, aes(distance_microns, gray_value)) + geom_point() + geom_hline(yintercept=400))
  
  # perform the axon computation
  axon_mask_dapi <- res_dapi %>% 
    subset(is_axon == 1 & gray_value < dapi_thresh) %>% 
    dplyr::select(distance_microns)

  
  res_axon <- res %>% inner_join(axon_mask_dapi)
  #print(ggplot(res, aes(distance_microns, gray_value)) + geom_point())



  axon_length_bins <- seq(0, max(res_axon$distance_microns), 
                          by=max(res_axon$distance_microns)/num_bins)
  axon_length_bins_offset <- axon_length_bins + (max(res_axon$distance_microns)/num_bins)/2
  
  
  res_axon$bin <- .bincode(res_axon$distance_microns, 
                           axon_length_bins, include.lowest = TRUE)
  
  res_axon_offset <- res_axon
  res_axon_offset$bin <- .bincode(res_axon$distance_microns,
                                  axon_length_bins_offset, include.lowest = TRUE)+0.5
  
  res_axon <- data.frame(rbind(res_axon, res_axon_offset))
  
  res_axon_binvals <- res_axon %>% 
    group_by(bin) %>% 
    summarize(bin_value_mean=mean(gray_value),
              bin_stdev=sd(gray_value))
  
  res_axon_binvals$bin_name <- sprintf("Bin %s", res_axon_binvals$bin)
  res_axon_binvals$bin_name <- factor(res_axon_binvals$bin_name, 
                                      levels=sprintf("Bin %s", 1:num_bins))
  
  # perform the soma computation
  soma_mask_dapi <- res_dapi %>%
    subset(is_axon == 0) %>% 
    dplyr::select(distance_microns)
  
  res_soma <- res %>% inner_join(soma_mask_dapi)
  res_soma_meangray <- data.frame(bin=0, bin_name="Soma",
                                  bin_value_mean=mean(res_soma$gray_value),
                                  bin_stdev=sd(res_soma$gray_value))


  res_binvals <- data.frame(rbind(res_axon_binvals, res_soma_meangray))
  
  
  res_binvals$fname <- file

  return(res_binvals)
  
}

get_neuron_channel_meangray_bin <- function(neuron_meta) {
  dapi_file <- subset(neuron_meta, channel == "DAPI")$fname
  dapi_thresh <- neuron_meta$dapi_thresh[1]
  
  channel_bins <- data.frame(do.call(rbind, 
                          lapply(neuron_meta$fname, 
                                 function(x) get_channel_meangray_bin(x, dapi_file, dapi_thresh))))
  

  return(channel_bins)

}

calculate_all_neurons <- function(rbp_meta) {
  rbp_bin_result <- data.frame(do.call(rbind,
                            tibble(rbp_meta) %>% 
    group_by(stage, neuron) %>%
    group_map(~get_neuron_channel_meangray_bin(.x)))) %>%
    inner_join(rbp_meta, by="fname")
  
  return(rbp_bin_result)
}

calculate_stats_all_neurons <- function(rbp_neuron_bins) {
  bin_stats <- rbp_neuron_bins %>% 
    subset(!is.na(bin_stdev)) %>%
    group_by(stage, channel, bin, bin_name) %>%
    summarize(mean_gray_value=mean(bin_value_mean),
              stderr_gray_value=(sqrt(sum(bin_stdev^2))/length(bin_stdev))/sqrt(length(bin_stdev)),
              num_values_bin=length(bin_value_mean))
  
  return(bin_stats)
    
}

plot_soma <- function(neuron_bins, plot_savefile, h=80, w=80) {
  
  for (chan in c("FITC", "TRITC", "Cy5")) {
    plt_mean <- ggplot(subset(neuron_bins, bin_name == "Soma" & channel==chan ),
       aes(stage, bin_value_mean, fill=stage)) + geom_boxplot() +
                 ggtitle(sprintf("Channel: %s, Stat: mean", chan)) +
      theme(legend.position="None") +
      theme_classic()

    
    plot(plt_mean)
    
    ggsave(gsub(".pdf", sprintf("_%s.pdf", chan), plot_savefile),
           plot=plt_mean,
           height=h,
           width=w,
           units="mm",
           useDingbats=FALSE)


  }
  
}


plot_axon <- function(neuron_bin_stats, plot_savefile, h=160, w=80) {
  
  for (chan in c("FITC", "TRITC", "Cy5")) {
    
    df_to_plot <- subset(neuron_bin_stats, bin_name != "Soma" & channel == chan & num_values_bin >= 5)
    # need to reorder the bins to get proximal to distal in right order
    maxval_bins <- max(df_to_plot$bin)
    df_to_plot$bin <- maxval_bins+1-df_to_plot$bin
    
    plt <- ggplot(df_to_plot,
                  aes(bin, mean_gray_value, group=stage, color=stage)) +
      geom_line() +
      geom_ribbon(data=df_to_plot,
                  mapping=aes(x=bin, 
                              ymin=mean_gray_value-stderr_gray_value,
                              ymax=mean_gray_value+stderr_gray_value,
                              fill=stage,
                              group=stage),
                  alpha=0.5) +
      ggtitle(sprintf("Channel: %s, Stat: mean", chan)) +
      theme_classic() +
      coord_flip()
    
    plot(plt)
    
    ggsave(gsub(".pdf", sprintf("_%s.pdf", chan), plot_savefile), plot=plt,
           device="pdf",
           height=h, 
           width=w,
           units="mm",
           useDingbats=FALSE)


  }
  
}
```

Calculate the axon stats and plot Weird bump in the middle of P1 Cpeb4
could be real? not a dapi signal. Look at neuron 5 image, for example.

``` r
cpeb4_meta <- get_meta_from_files("Cpeb4")
cpeb4_meta$dapi_thresh <- 500
get_dapi_profiles(cpeb4_meta)
```

![](ICC_plotprofile_intensity_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
cpeb4_neuron_bins <- suppressMessages(calculate_all_neurons(cpeb4_meta))

plot_soma(cpeb4_neuron_bins, file.path(plot_dir, "Cpeb4_soma_intensities.pdf"))
```

![](ICC_plotprofile_intensity_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->![](ICC_plotprofile_intensity_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->![](ICC_plotprofile_intensity_files/figure-gfm/unnamed-chunk-6-4.png)<!-- -->

``` r
cpeb4_neuron_bin_stats <- calculate_stats_all_neurons(cpeb4_neuron_bins)
plot_axon(cpeb4_neuron_bin_stats, file.path(plot_dir, "Cpeb4_axon_profiles.pdf")) 
```

![](ICC_plotprofile_intensity_files/figure-gfm/unnamed-chunk-6-5.png)<!-- -->![](ICC_plotprofile_intensity_files/figure-gfm/unnamed-chunk-6-6.png)<!-- -->![](ICC_plotprofile_intensity_files/figure-gfm/unnamed-chunk-6-7.png)<!-- -->

``` r
rbms1_meta <- get_meta_from_files("Rbms1")
rbms1_meta$dapi_thresh <- 500
get_dapi_profiles(rbms1_meta)
```

    ## Warning: Removed 1 rows containing non-finite values (`stat_density()`).

![](ICC_plotprofile_intensity_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
rbms1_neuron_bins <- suppressMessages(calculate_all_neurons(rbms1_meta))

plot_soma(rbms1_neuron_bins, file.path(plot_dir, "Rbms1_soma_intensities.pdf"))
```

![](ICC_plotprofile_intensity_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->![](ICC_plotprofile_intensity_files/figure-gfm/unnamed-chunk-7-3.png)<!-- -->![](ICC_plotprofile_intensity_files/figure-gfm/unnamed-chunk-7-4.png)<!-- -->

``` r
rbms1_neuron_bin_stats <- calculate_stats_all_neurons(rbms1_neuron_bins)
plot_axon(rbms1_neuron_bin_stats, file.path(plot_dir, "Rbms1_axon_profiles.pdf")) 
```

![](ICC_plotprofile_intensity_files/figure-gfm/unnamed-chunk-7-5.png)<!-- -->![](ICC_plotprofile_intensity_files/figure-gfm/unnamed-chunk-7-6.png)<!-- -->![](ICC_plotprofile_intensity_files/figure-gfm/unnamed-chunk-7-7.png)<!-- -->

``` r
pum1_meta <- get_meta_from_files("Pum1")
pum1_meta$dapi_thresh <- 500
get_dapi_profiles(pum1_meta)
```

![](ICC_plotprofile_intensity_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
pum1_neuron_bins <- suppressMessages(calculate_all_neurons(pum1_meta))

plot_soma(pum1_neuron_bins, file.path(plot_dir, "Pum1_soma_intensities.pdf"))
```

![](ICC_plotprofile_intensity_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->![](ICC_plotprofile_intensity_files/figure-gfm/unnamed-chunk-8-3.png)<!-- -->![](ICC_plotprofile_intensity_files/figure-gfm/unnamed-chunk-8-4.png)<!-- -->

``` r
pum1_neuron_bin_stats <- calculate_stats_all_neurons(pum1_neuron_bins)
plot_axon(pum1_neuron_bin_stats, file.path(plot_dir, "Pum1_axon_profiles.pdf")) 
```

![](ICC_plotprofile_intensity_files/figure-gfm/unnamed-chunk-8-5.png)<!-- -->![](ICC_plotprofile_intensity_files/figure-gfm/unnamed-chunk-8-6.png)<!-- -->![](ICC_plotprofile_intensity_files/figure-gfm/unnamed-chunk-8-7.png)<!-- -->

``` r
pum2_meta <- get_meta_from_files("Pum2")
pum2_meta$dapi_thresh <- 500
get_dapi_profiles(pum2_meta)
```

![](ICC_plotprofile_intensity_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
pum2_neuron_bins <- suppressMessages(calculate_all_neurons(pum2_meta))

plot_soma(pum2_neuron_bins, file.path(plot_dir, "Pum2_soma_intensities.pdf"))
```

![](ICC_plotprofile_intensity_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->![](ICC_plotprofile_intensity_files/figure-gfm/unnamed-chunk-9-3.png)<!-- -->![](ICC_plotprofile_intensity_files/figure-gfm/unnamed-chunk-9-4.png)<!-- -->

``` r
pum2_neuron_bin_stats <- calculate_stats_all_neurons(pum2_neuron_bins)
plot_axon(pum2_neuron_bin_stats, file.path(plot_dir, "Pum2_axon_profiles.pdf")) 
```

![](ICC_plotprofile_intensity_files/figure-gfm/unnamed-chunk-9-5.png)<!-- -->![](ICC_plotprofile_intensity_files/figure-gfm/unnamed-chunk-9-6.png)<!-- -->![](ICC_plotprofile_intensity_files/figure-gfm/unnamed-chunk-9-7.png)<!-- -->

``` r
wave1_meta <- get_meta_from_files("WAVE1")
wave1_meta$dapi_thresh <- 500
get_dapi_profiles(wave1_meta)
```

![](ICC_plotprofile_intensity_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
wave1_neuron_bins <- suppressMessages(calculate_all_neurons(wave1_meta))

wave1_neuron_bins <- wave1_neuron_bins %>%
  rowwise() %>%
  mutate(neuron=sprintf("%s_%s", stage, neuron),
         stage="P3")
wave1_neuron_bin_stats <- calculate_stats_all_neurons(wave1_neuron_bins)
plot_axon(wave1_neuron_bin_stats, file.path(plot_dir, "Wave1_axon_profiles.pdf")) 
```

![](ICC_plotprofile_intensity_files/figure-gfm/unnamed-chunk-10-2.png)<!-- -->![](ICC_plotprofile_intensity_files/figure-gfm/unnamed-chunk-10-3.png)<!-- -->![](ICC_plotprofile_intensity_files/figure-gfm/unnamed-chunk-10-4.png)<!-- -->