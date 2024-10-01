## Quality control analysis of DIA data using [DIA-NN 1.9.1](https://github.com/vdemichev/DiaNN/releases/tag/1.9.1) search results

### General information
In this repository you will find the quarto markdown **`QC_report_DIANN.qmd`** for QC check from DIANN search results and the shiny app **`QC_report_DIANN_app.r`** that will help you to vizualise the results without coding skills. All you need to do is to download the **`QC_report_DIANN_app.r`**, open it in RStudio and click **`Run app`**. The app will start, then you browse the `report.parquet` file from DIANN. Naturally, make sure you have the required libraries installed in order to have the app working properly.

### For the quarto markdown document
Please, note that the column names are slightly different in `report.tsv` and `report.parquet` files. Some importante informations are only available in the `report.parquet`.
You can download the .qmd document and place it in the same folder as the `report.parquet` file.

To have an in depth understanding of what is happening in this script, please check the [DIANN documentation](https://github.com/vdemichev/DiaNN), as well as the original publication of the QuantUMS algorithm by Demichev's group:
>Franziska Kistner, Justus L. Grossmann, Ludwig R. Sinn, Vadim Demichev. QuantUMS: uncertainty minimisation enables confident quantification in proteomics. bioRxiv 2023.06.20.545604; doi: https://doi.org/10.1101/2023.06.20.545604.

Maybe it is worth to check the [diann R package](https://github.com/vdemichev/diann-rpackage) that we are using to extract the matrix of abundance from the DIANN report with MaxLFQ values.

Make sure you have the following packages installed and loaded:
```r
library(diann) # to extract the MaxLFQ matrix from DIANN report
library(arrow)  # to read the report.parquet file
library(here) # to avoid the need for use the path while loading the data
library(tidyverse)  # to do the data wrangling, plots, etc...
library(ggpointdensity) # to reconstruct the m/z density map
library(naniar) # for sparsity analysis
library(limma)  # for statistics (normalization)
library(plotly) # for 3D plot of the QuantUMS scores
```

## Load and filter the data from DIANN 1.9.1 search results
The file recommended to be used is the `report.parquet` generated by the DIANN search.
At the loading, we filter the data by the Q-values of the library, library protein group and protein group. It is important to use the same cut-off for Q-values at the Protein Group and Library Protein Group.

Information from DIANN documentation:
- **Quantity.Quality:** when using QuantUMS is equal to 1.0 / (1.0 + SD), where SD is the standard deviation of the LC-MS-derived error in relative precursor quantification.
- **Empirical.Quality:** when using QuantUMS reflects the agreement of relative precursor quantities obtained using different quantitative features (MS1 / fragment ions).
- **PG.MaxLFQ.Quality:** when using QuantUMS reflects the quality of PG.MaxLFQ.

#### Load the output data from DIANN search using the arrow function `read_parquet()`.

The File.Name column was removed from the output of the `report.parquet`, so we need to recreate it in order to make the file compatible with `dann_matrix()` function.

```r
diann_report <- arrow::read_parquet("report.parquet") %>%
    dplyr::filter(PG.MaxLFQ.Quality >= 0.75 & Lib.PG.Q.Value <= 0.01 & Lib.Q.Value <= 0.01 & PG.Q.Value <= 0.01) %>%
  dplyr::mutate(File.Name = Run)

# extracting the matrix of abundance from DIANN report.parquet file
# observe that if we select the id.header as "Protein.Ids" or "Protein.Group" we will have different results because of the way DIANN makes the protein inference.
unique_genes <- diann::diann_matrix(diann_report,
    id.header = "Protein.Ids",
    quantity.header = "Genes.MaxLFQ.Unique",
    proteotypic.only = TRUE,
    pg.q = .01
)

# count the number of proteins per sample and save it in a new data frame called proteins
proteins <- diann_report %>%
    dplyr::group_by(Run) %>%
    dplyr::summarise(
        n_proteins = n_distinct(Protein.Ids)
    )
```

If it is necessary to use the data in other pipelines, you can write a tsv file filtered for the whole data and for the matrix.
The files will be saved in the working directory and can be inspected in Excel, for instance.

```r
readr::write_tsv(diann_report, "diann_report_QuantUMS.tsv")
readr::write_tsv(unique_genes, "diann_matrix_QuantUMS.tsv")
```

For the reconstruction of the ion chromatograms, the precursor quantity is plotted over the retention time (min) for each sample.

:bulb: You can generate a new column for `condition, group, etc.` and add the `color = condition` to the `aes()`, for instance.

```r
precursor_rt <- diann_report %>%
    ggplot(aes(x = RT, y = Precursor.Quantity)) +
    geom_line(aes(color = Run), show.legend = FALSE) +
    labs(
        x = "Retention time (min)",
        y = "Precursor quantity",
        color = NULL
    ) +
    facet_wrap(~Run)

ggsave("precursor_rt.png",
    path = "plots",
    precursor_rt, width = 15,
    height = 10, units = "in", dpi = 350)
```

Plot the m/z map to show the density of ions collected over the scan range and retention time.
It can be very informative of instability in spray or chromatographic setup in general.

```r
mz_map_density_plot <- diann_report %>%
    ggplot(aes(x = RT, y = Precursor.Mz)) +
    ggpointdensity::geom_pointdensity(size = 0.25) +
    viridis::scale_color_viridis(option = "plasma") +
    scale_x_continuous(limits = c(0, 90)) + # replace 90 with the length of the gradient
    labs(
        x = "Retention time (min)",
        y = " Scan range (m/z)",
        color = NULL
    ) +
    theme(legend.position = "bottom",
        legend.key.width = unit(1.5, "cm")) +
    facet_wrap(~Run)

ggsave("mz_map_density_plot.png",
    path = "plots",
    mz_map_density_plot, width = 15,
    height = 10, units = "in", dpi = 350)
```

The charge state distribution can be informative in experiments using FAIMS, for instance.

```r
precursor_charge_density <- diann_report %>%
    ggplot(aes(x = Precursor.Charge, fill = Run)) +
    geom_density(alpha = 0.7, 
            stat = "density", 
            show.legend = FALSE) +
    labs(
        x = "Precursor charge",
        y = "Density",
        fill = NULL
    ) +
    facet_wrap(~Run)

ggsave("precursor_charge_density.png",
    path = "plots",
    precursor_charge_density, width = 10,
    height = 10, units = "in", dpi = 350)
```

Counting the number of peptides per sample.

```r
peptides_plot <- diann_report %>%
    dplyr::group_by(Run) %>%
    dplyr::summarise(
        n_pepetides = n_distinct(Stripped.Sequence)
    ) %>%
    ggplot(aes(y = Run, x = n_peptides, fill = Run)) +
    geom_bar(stat = "identity", position = "dodge", show.legend = FALSE) +
    geom_vline(xintercept = 40000, linetype = "dashed", color = "black") + # replace the 40000 with a threshold of interest
    geom_text(aes(label = n_peptides),
        color = "white", size = 5,
        hjust = 1, nudge_x = -0.5
        ) +
    labs(y = NULL,
        x = "Number of peptides",
        fill = NULL) +
    theme(axis.text.x = element_text(angle = 90,
                        vjust = 0.5,
                        hjust = 1)
    )

ggsave("peptides_plot.png",
    path = "plots",
    peptides_plot, width = 10,
    height = 10, units = "in", dpi = 350)
```

Counting the number of proteins per sample.

```r
proteins_plot <- proteins %>%
    ggplot(aes(y = Run, x = n_proteins, fill = Run)) +
    geom_bar(stat = "identity", position = "dodge", show.legend = FALSE) +
    geom_vline(xintercept = 4500, linetype = "dashed", color = "black") + # replace the 4500 with a threshold of interest
    geom_text(aes(label = n_proteins),
        color = "white", size = 5,
        hjust = 1, nudge_x = -0.5
        ) +
    labs(y = NULL,
        x = "Number of proteins",
        fill = NULL) +
    theme(axis.text.x = element_text(angle = 90,
                        vjust = 0.5,
                        hjust = 1)
    )

ggsave("proteins_plot.png",
    path = "plots",
    proteins_plot, width = 10,
    height = 10, units = "in", dpi = 350)
```

Evaluate the sparsity profile for each sample. If the sparsity is high in one sample, check the m/z map to understand why. 

```r
sparsity_plot <- unique_genes %>%
    as.data.frame() %>%
    naniar::vis_miss() +
    labs(x = NULL,
        y = "Proteins") +
    theme(text = element_text(color = "black"),
        axis.text.y = element_text(color = "black", vjust = 1),
        axis.text.x = element_text(angle = 90, vjust = 0.5,
                                   hjust = 0, color = "black"),
        line = element_blank()
    )

# save the plot in the working directory
ggsave("sparsity_plot.png",
    path = "plots",
    sparsity_plot, width = 10, bg = "white",
    height = 10, units = "in", dpi = 350)
```

Calculate the proportion of missing values and median abundance per sample and plot the correlation between them.
This step is recommended by Prof. Dr. Clemens Kreutz (Institute of Medical Biometry and Statistics).
You will see that, depending on the normalization method used, the correlation will change.
Try to compare the normalization using scale and quantile, for instance.

```r
sample_abundance_vs_missing <- log2(unique_genes) %>%
    as.data.frame() %>%
    gather(key = "Sample", value = "Intensity") %>%
    dplyr::mutate(missing = is.na(Intensity)) %>%
    dplyr::group_by(Sample) %>%
    dplyr::summarise(
        missing = mean(missing) * 100,
        median_intensity = median(Intensity, na.rm = TRUE)
    )

corr_mean_missing <- sample_abundance_vs_missing %>%
    ggplot(aes(x = missing, y = median_intensity)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_smooth(method = "lm", se = FALSE,
        color = "red") +
    labs(x = "Proportion of missing values (%)",
        y = "Mean abundance",
        color = NULL)

ggsave("corr_mean_missing.png",
    path = "plots",
    corr_mean_missing, width = 3,
    height = 3, units = "in", dpi = 350)
```

Plot the log2 transformed abundance distribution before and after the normalization.

The normalization methods available in `limma` package are: "none", "scale", "quantile", "cyclicloess", "Aquantile", "Gquantile", "Rquantile" or "Tquantile".

From `limma` documentation:
>Scale normalization was proposed by Yang et al (2001, 2002) and is further explained by Smyth and Speed (2003). The idea is simply to scale the log-ratios to have the same median-absolute-deviation (MAD) across arrays. This idea has also been implemented by the maNormScale function in the marray package. The implementation here is slightly different in that the MAD scale estimator is replaced with the median-absolute-value and the A-values are normalized as well as the M-values.

```r
df_long_raw <- unique_genes %>%
    log2() %>%
    as.data.frame() %>%
    gather(key = "Sample", value = "Intensity") %>%
    dplyr::mutate(norm = "Raw matrix")

df_long_norm <- log2(unique_genes) %>%  # if there are Inf values we can add 0.5 to the intensity value before log transformation
    limma::normalizeBetweenArrays(method = "scale") %>% # we can change the method to any other available in limma and to compare the profiles
    as.data.frame() %>%
    gather(key = "Sample", value = "Intensity") %>%
    dplyr::mutate(norm = "MAD normalised")

combined_data <- rbind(df_long_raw, df_long_norm) %>%
    dplyr::mutate(norm = factor(norm, levels = c("Raw matrix", "MAD normalised")))

# plot the intensity distribution on log2 scale
intensity_distribution <- combined_data %>%
    ggplot(aes(
        x = Sample,
        y = Intensity,
        fill = Sample
    )) +
    geom_boxplot() +
    theme(
        axis.text.x = element_text(angle = 90, 
                hjust = 1, vjust = 0.5),
                legend.position = "none"
    ) +
    labs(x = NULL,
        y = "log2(Intensity)",
        fill = NULL) +
  facet_wrap(~norm)

ggsave("intensity_distribution.png",
    path = "plots",
    intensity_distribution, width = 10,
    height = 8, units = "in", dpi = 350)
```

Plot the error of the retention time across the m/z range.
A better understanding of how retention time can be predicted is available in:
>Al Musaimi O, Valenzo OMM, Williams DR. Prediction of peptides retention behavior in reversed-phase liquid chromatography based on their hydrophobicity. J Sep Sci. 2023 Jan;46(2):e2200743. doi: 10.1002/jssc.202200743. Epub 2022 Nov 14. PMID: 36349538; PMCID: PMC10098489.

```r
RT_error <- diann_report %>%
    ggplot(aes(x = Precursor.Mz, 
                y = RT - Predicted.RT)) +
    ggpointdensity::geom_pointdensity(size = 0.25) +
    viridis::scale_color_viridis(option = "plasma") +
    geom_hline(yintercept = c(1,0,-1), linetype = "dashed", color = "black") +
    labs(
        x = "Precursor Mz",
        y = "RT - Predicted RT",
        color = NULL
    ) +
    theme(legend.key.width = unit(2, "cm"),
    legend.position = "bottom") +
    facet_wrap(~Run)

ggsave("RT_error.png",
    path = "plots",
    RT_error, width = 15,
    height = 10, units = "in", dpi = 350)
```

Plot the Posterior Error Probability (PEP) density distribution per Run. This is a estimation of local FDR.

```r
PEP_plot <- diann_report %>%
    ggplot(aes(x = PEP)) +
    geom_density(fill = "tomato") +
    labs(
        x = "PEP",
        y = "Density",
        fill = NULL
    ) +
    theme(legend.position = "none") +
    facet_wrap(~Run)

ggsave("PEP_density.png",
    path = "plots",
    PEP_plot, width = 15,
    height = 10, units = "in", dpi = 350)
```

Plot the distribution of the FWHM per run.

Information from DIANN documentation:
>**FWHM** estimated peak width at half-maximum; note that the accuracy of such estimates sometimes strongly depends on the DIA cycle time and sample injection amount, i.e. they can only be used to evaluate chromatographic performance in direct comparisons with similar settings, including the scan window; another caveat is that FWHM does not reflect any peak tailing.

```r
FWHM_plot <- diann_report %>%
    ggplot(aes(x = FWHM)) +
    geom_density(fill = "tomato") +
    labs(
        x = "FWHM",
        y = "Density",
        fill = NULL
    ) +
    theme(legend.position = "none") +
    facet_wrap(~Run)

ggsave("FWHM_density.png",
    path = "plots",
    FWHM_plot, width = 15,
    height = 10, units = "in", dpi = 350)
```

### evaluate the digestion efficiency by plotting the missed cleavages
We can separate the peptides by the specificity of the enzyme used in the digestion. In this case, we are using Trypsin, so we can evaluate the missed cleavages by the presence of K or R at the C-terminal of the peptide.

```r
missed_cleavages <- diann_report %>%
    dplyr::mutate(specificity = case_when(
        str_detect(Stripped.Sequence, "K$|R$") ~ "Trypsin",
        TRUE ~ "Unespecific")
    ) %>%
    group_by(Run, specificity) %>%
    dplyr::summarise(
        peptides = n()
    ) %>%
    ggplot(aes(x = Run, y = peptides,
                fill = specificity)) +
    geom_col(alpha = 0.7,
            position = "dodge") +
    geom_text(aes(label = peptides),
        position = position_dodge(width = 1),
        vjust = -0.25, size = 3) +
    scale_fill_manual(values = c("Trypsin" = "tomato",
                                "Unespecific" = "steelblue")) +
    labs(
        x = NULL,
        y = "Count",
        fill = "Specificity"
    ) +
    theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90,
                        hjust = 1, vjust = 0.5))

ggsave("missed_cleavages.png",
    path = "plots",
    missed_cleavages, width = 10,
    height = 10, units = "in", dpi = 350)
```

### Evaluate the distribution of the QuantUMS scores in a 3D plot

The QuantUMS scores from DIANN are used to evaluate the quality of the data. DIANN calculates three different scores: quantity quality, empirical quality and PG MaxLFQ quality. Filtering by the empirical quality has the most significant impact on the data, while the other two only impacts in a small proportion.

From the [pre-print](https://www.biorxiv.org/content/10.1101/2023.06.20.545604v1), the QuantUMS is explained as follows:

>"The idea here is that since each feature produces, for each acquisition, an estimate of the precursor quantity, the deviations between these estimates corresponding to different features are indicative of how accurate the quantity estimates are. QuantUMS hence tunes hyperparameters to minimise the empirically measured differences between quantity estimates obtained using different features."

```r
diann_report %>%
    plot_ly(x = ~PG.MaxLFQ.Quality, 
        y = ~Quantity.Quality, 
        z = ~Empirical.Quality,
            color = ~Run, alpha = 0.7,
            type = "scatter3d", mode = "markers") %>%
    layout(scene = list(
        xaxis = list(title = "PG MaxLFQ Quality"),
        yaxis = list(title = "Quantity Quality"),
        zaxis = list(title = "Empirical Quality")
    )) %>%
    subplot()
```

Demichev asking questions about the values of QuantUMS scores said:

>"Basically, low values mean that something is wrong with quantity, high values don't guarantee that it is good though. In QuantUMS, high values mean that LC-MS-related error in the quantity is very likely negligible."

We can check the correlation between MS1 and MS2 chromatograms stratified by a defined empirical quality score. This will show that the correlation is lower for peptides with low empirical scores.

```r
EQScore_cutoff_plot <- diann_report %>%
    dplyr::mutate(EQScore_cutoff = case_when(
        Empirical.Quality >= 0.75 ~ "≥ 0.75",
        TRUE ~ "< 0.75")
    ) %>%
    ggplot() +
    geom_density(aes(x = Ms1.Profile.Corr, fill = EQScore_cutoff), alpha = 0.7) +
    scale_fill_manual(values = c("tomato", "steelblue")) +
    labs(x = "MS1 Profile Correlation",
        y = "Density",
        fill = "Correlation between MS1 and MS2 chromatograms") +
    theme(legend.position = "top") +
    facet_wrap(~Run)

ggsave("EQScore_cutoff_plot.png",
    path = "plots",
    EQScore_cutoff_plot, width = 10,
    height = 10, units = "in", dpi = 350)
```
