## Dashboard QC report for DIANN search results
## Alison Felipe Alencar Chaves
## For detailed explanation, check out the GitHub repository: https://github.com/41ison/QC4DIANN

# Load required libraries
library(shiny)
library(shinydashboard)
library(diann)
library(arrow)
library(tidyverse)
library(ggpointdensity)
library(limma)
library(vegan)
library(lsa)
library(plotly)
library(viridis)
library(ggfortify)

# Increase the maximum file size to 900 MB
options(shiny.maxRequestSize = 900 * 1024^2)

# set the general theme for the plots
theme_set(theme_bw())
theme_update(
  text = element_text(color = "black"),
  axis.text = element_text(color = "black"),
  axis.title = element_text(color = "black", face = "bold"),
  strip.text = element_text(face = "bold"),
  legend.title = element_text(face = "bold", hjust = 0.5),
  legend.title.position = "top"
)

# Define UI for application that reads a parquet file and generates a QC report dashboard
ui <- dashboardPage(

  dashboardHeader(
      title = "QC Reporting dashboard for DIANN search results", titleWidth = "400",
      dropdownMenu(
        type = "messages",
        messageItem(
          from = "Support",
          message = "felipealison@gmail.com",
          icon = icon("envelope")
        )
      )
  ),

  dashboardSidebar(
    sidebarMenu(
      menuItem("QuantUMS filters", tabName = "filters", icon = icon("filter")),
      fileInput(inputId = "report", label = "Choose Parquet File", accept = ".parquet"),
      sliderInput("PG.MaxLFQ.Quality", "PG MaxLFQ Quality score", min = 0, max = 1, value = 0.75, step = 0.05),
      sliderInput("Empirical.Quality", "Empirical Quality score", min = 0, max = 1, value = 0, step = 0.05),
      menuItem("Interactive viewer", tabName = "protein", icon = icon("equalizer", lib = "glyphicon")),
      selectInput("xcol", "X Sample", choices = NULL),
      selectInput("ycol", "Y Sample", choices = NULL),
      downloadButton("download", "Download filtered matrix", class = "butt"),
      tags$head(tags$style(".butt{background:grey;} .butt{color: #337ab7;}"))
    )
  ),
  
  dashboardBody(
    tabItems(
      tabItem(tabName = "filters",
              fluidRow(
                infoBoxOutput("info_box1", width = 12),
                box(title = "Reconstruction of XIC", status = "primary", solidHeader = TRUE, plotOutput("plot1"), collapsible = TRUE),
                box(title = "Density of ions", status = "primary", solidHeader = TRUE, plotOutput("plot2"), collapsible = TRUE),
                box(title = "Retention time error", status = "primary", solidHeader = TRUE, plotOutput("plot3"), collapsible = TRUE),
                box(title = "Charge state distribution", status = "primary", solidHeader = TRUE, plotOutput("plot4"), collapsible = TRUE),
                box(title = "Peptide length", status = "primary", solidHeader = TRUE, plotOutput("plot5"), collapsible = TRUE),
                box(title = "Peptides per sample", status = "primary", solidHeader = TRUE, plotOutput("plot6"), collapsible = TRUE),
                box(title = "Proteins per sample", status = "primary", solidHeader = TRUE, plotOutput("plot7"), collapsible = TRUE),
                box(title = "Sparsity profile", status = "primary", solidHeader = TRUE, plotOutput("plot8"), collapsible = TRUE),
                box(title = "Missing vs median abundance", status = "primary", solidHeader = TRUE, plotOutput("plot9"), collapsible = TRUE),
                box(title = "Abundance before and after MAD normalization", status = "primary", solidHeader = TRUE, plotOutput("plot10"), collapsible = TRUE),
                box(title = "Missed cleavage sites", status = "primary", solidHeader = TRUE, plotOutput("plot11"), collapsible = TRUE),
                box(title = "MS1 Profile Correlation", status = "primary", solidHeader = TRUE, plotOutput("plot12"), collapsible = TRUE),
                box(title = "QuantUMS scores distribution", status = "primary", solidHeader = TRUE, plotOutput("plot13"), collapsible = TRUE),
                box(title = "Gene quantity distribution", status = "primary", solidHeader = TRUE, plotOutput("plot14"), collapsible = TRUE)
              )
      ),
      tabItem(tabName = "protein",
              fluidRow(
                box(title = "Sample correlation - Non-normalized log2(Intensity)", status = "primary", height = 600, solidHeader = TRUE, plotlyOutput("Corr"), collapsible = FALSE),
                tabBox(
                  title = "Similarity metrics", side = "right", height = 600,
                  tabPanel("Cosine similarity", plotOutput("cosine_similarity")),
                  tabPanel("Euclidean distance", plotOutput("euclidean_distance")),
                  tabPanel("Jaccard similarity", plotOutput("jaccard_similarity"))
                ),
                box(title = "QuantUMS score distribution", status = "primary", height = 600, solidHeader = TRUE, plotlyOutput("QuantUMS_dist"), collapsible = FALSE),
                box(title = "Principal Component Analysis", status = "primary", height = 600, solidHeader = TRUE, plotlyOutput("PCA"), collapsible = FALSE)
          )
      )
    )
  )
)

# Define server logic required to read the parquet file and generate the QC report
server <- function(input, output, session) {

  # Reactive expression to read and pre-process the uploaded parquet file
  data <- reactive({
    req(input$report)
    diann_report <- arrow::read_parquet(input$report$datapath) %>%
      dplyr::filter(Lib.PG.Q.Value <= 0.01 & Lib.Q.Value <= 0.01 & PG.Q.Value <= 0.01) %>%
      dplyr::mutate(File.Name = Run,
                    peptide_length = nchar(Stripped.Sequence)) %>%
      dplyr::filter(.$PG.MaxLFQ.Quality >= input$PG.MaxLFQ.Quality & .$Empirical.Quality >= input$Empirical.Quality)
  })

  # Reactive expression to filter number of proteins based on the input filters
    proteins <- reactive({
    req(input$report)
    diann_report <- arrow::read_parquet(input$report$datapath) %>%
      dplyr::filter(Lib.PG.Q.Value <= 0.01 & Lib.Q.Value <= 0.01 & PG.Q.Value <= 0.01) %>%
      dplyr::mutate(File.Name = Run) %>%
      dplyr::filter(.$PG.MaxLFQ.Quality >= input$PG.MaxLFQ.Quality & .$Empirical.Quality >= input$Empirical.Quality) %>%
      dplyr::group_by(Run) %>%
      dplyr::summarise(
        n_proteins = n_distinct(Protein.Ids)
      )
  })

  # Reactive expression to filter number of peptides based on the input filters
    peptides_per_run <- reactive({
    req(input$report)
    diann_report <- arrow::read_parquet(input$report$datapath) %>%
      dplyr::filter(Lib.PG.Q.Value <= 0.01 & Lib.Q.Value <= 0.01 & PG.Q.Value <= 0.01) %>%
      dplyr::mutate(File.Name = Run) %>%
      dplyr::filter(.$PG.MaxLFQ.Quality >= input$PG.MaxLFQ.Quality & .$Empirical.Quality >= input$Empirical.Quality) %>%
      dplyr::group_by(Run) %>%
      dplyr::summarise(
        n_peptides = n_distinct(Stripped.Sequence)
      )
  })

    # Reactive expression to filter matrix of protein abundance based on the input filters
    unique_genes <- reactive({
    req(data())
    diann::diann_matrix(data(),
      id.header = "Protein.Ids",
      quantity.header = "Genes.MaxLFQ.Unique",
      proteotypic.only = FALSE,
      pg.q = .01
    )
  })

  # Observe the uploaded file and update selectInput choices
  observe({
    req(unique_genes())
    colnames <- colnames(unique_genes())
    updateSelectInput(session, "xcol", choices = colnames)
    updateSelectInput(session, "ycol", choices = colnames)
  })

  # Reactive expression to combine the raw and MAD normalized data
combined_data <- reactive({
  req(unique_genes())
  unique_genes() %>%
    log2() %>%
    as.data.frame() %>%
    gather(key = "Sample", value = "Intensity") %>%
    dplyr::mutate(norm = "Raw matrix") %>%
    bind_rows(
      unique_genes() %>%
        log2() %>%
        limma::normalizeBetweenArrays(method = "scale") %>%
        as.data.frame() %>%
        gather(key = "Sample", value = "Intensity") %>%
        dplyr::mutate(norm = "MAD normalised")
    ) %>%
    dplyr::mutate(norm = factor(norm, levels = c("Raw matrix", "MAD normalised")))
})

  QuantUMS_scores <- reactive({
    req(input$report)
    diann_report <- arrow::read_parquet(input$report$datapath) %>%
      dplyr::filter(Lib.PG.Q.Value <= 0.01 & Lib.Q.Value <= 0.01 & PG.Q.Value <= 0.01) %>%
      dplyr::mutate(File.Name = Run) %>%
      dplyr::select(Run, Precursor.Id, PG.MaxLFQ.Quality, Empirical.Quality, Quantity.Quality) %>%
      pivot_longer(-c(Run, Precursor.Id),
                    names_to = "Filter",
                    values_to = "Score")
  })

  gene_quantity <- reactive({
    req(data())
    data() %>%
      dplyr::select(Run, PG.MaxLFQ, Genes.MaxLFQ, Genes.MaxLFQ.Unique) %>%
      pivot_longer(-Run,
                   names_to = "protein_metrics",
                   values_to = "quantity")
  })
  
  missing_vs_mean <- reactive({
    unique_genes() %>%
    log2() %>%
    as.data.frame() %>%
    gather(key = "Sample", value = "Intensity") %>%
    dplyr::mutate(missing = is.na(Intensity)) %>%
    dplyr::group_by(Sample) %>%
    dplyr::summarise(
        missing = mean(missing) * 100,
        median_intensity = median(Intensity, na.rm = TRUE)
    ) %>%
    dplyr::mutate(norm = "Non normalised") %>%
    bind_rows(
      unique_genes() %>%
        log2() %>%
        limma::normalizeBetweenArrays(method = "scale") %>%
        as.data.frame() %>%
        gather(key = "Sample", value = "Intensity") %>%
        dplyr::mutate(missing = is.na(Intensity)) %>%
        dplyr::group_by(Sample) %>%
        dplyr::summarise(
            missing = mean(missing) * 100,
            median_intensity = median(Intensity, na.rm = TRUE)
        ) %>%
        dplyr::mutate(norm = "MAD normalised")
    ) %>%
    dplyr::mutate(norm = factor(norm, levels = c("Non normalised", "MAD normalised")))
  })

output$info_box1 <- renderInfoBox({
    infoBox(title = "The filters Lib.PG.Q.Value ≤ 0.01, Lib.Q.Value ≤ 0.01 and PG.Q.Value ≤ 0.01 are active.",
            paste("PG MaxLFQ Quality score ≥ ", input$PG.MaxLFQ.Quality),
            paste("Empirical Quality score ≥ ", input$Empirical.Quality),
            icon = icon("filter"),
            color = "black"
    )
  })

MS_corr <- reactive({
  req(input$report)
  MS_correlation <- arrow::read_parquet(input$report$datapath) %>%
    dplyr::filter(Lib.PG.Q.Value <= 0.01 & Lib.Q.Value <= 0.01 & PG.Q.Value <= 0.01) %>%
    dplyr::mutate(File.Name = Run)
})

PCA_label <- reactive({
  unique_genes() %>%
  log2() %>%
  na.omit() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample")
})

pca_data <- reactive({
    unique_genes() %>%
    log2() %>%
    na.omit() %>%
    t() %>%
    prcomp(., scale. = TRUE) %>%
    autoplot(data = PCA_label(), colour = "Sample", label = TRUE)
})

# download the filtered matrix
output$download <- downloadHandler(
    filename = function() {
        paste0(input$report, ".tsv")
    },
    content = function(file) {
        readr::write_tsv(as.data.frame(unique_genes()) %>% rownames_to_column("protein_id"), file)
    }
)

# calculate the cosine similarity in the matrix and plot the heatmap
output$cosine_similarity <- renderPlot({
    unique_genes() %>%
    log2() %>%
    na.omit() %>%
    lsa::cosine() %>%
    as.data.frame() %>%
    rownames_to_column(var = "Sample") %>%
    pivot_longer(-Sample, names_to = "Match", values_to = "value") %>%
    dplyr::mutate(Similarity = "Cosine similarity") %>%
    ggplot() +
    geom_tile(aes(x = Sample, y = Match, fill = value)) +
    viridis::scale_fill_viridis(option = "E") +
    theme(text = element_text(size = 15),
        axis.text.x = element_text(angle = 90,
                        hjust = 1, vjust = 0.5),
        axis.text.y = element_text(angle = 0,
                        hjust = 1, vjust = 0.5),
        legend.position = "bottom",
        legend.key.width = unit(2.5, "cm")) +
    labs(x = NULL,
        y = NULL,
        fill = "Cosine similarity")
})

# calculate the euclidean distance in the matrix and plot the heatmap
output$euclidean_distance <- renderPlot({
    unique_genes() %>%
    log2() %>%
    t() %>%
    dist(method = "euclidean") %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column(var = "Sample") %>%
    pivot_longer(-Sample, names_to = "Match", values_to = "value") %>%
    dplyr::mutate(Similarity = "Euclidean distance") %>%
    ggplot() +
    geom_tile(aes(x = Sample, y = Match, fill = value)) +
    viridis::scale_fill_viridis(option = "E") +
    theme(text = element_text(size = 15),
    axis.text.x = element_text(angle = 90,
                        hjust = 1, vjust = 0.5),
        axis.text.y = element_text(angle = 0,
                        hjust = 1, vjust = 0.5),
        legend.position = "bottom",
        legend.key.width = unit(2.5, "cm")) +
    labs(x = NULL,
        y = NULL,
        fill = "Euclidean distance")
        })

# calculate the Jaccard similarity in the matrix and plot the heatmap
output$jaccard_similarity <- renderPlot({
    unique_genes() %>%
    log2() %>%
    t() %>%
    vegan::vegdist(method = "jaccard", na.rm = TRUE) %>%
    as.matrix() %>%
    as.data.frame(as.table(.)) %>% 
    dplyr::mutate(Sample = colnames(.)) %>%
    pivot_longer(-Sample, names_to = "Match", values_to = "value") %>%
    dplyr::mutate(Similarity = "Jaccard similarity") %>%
    ggplot() +
    geom_tile(aes(x = Sample, y = Match, fill = value)) +
    viridis::scale_fill_viridis(option = "E") +
    theme(text = element_text(size = 15),
    axis.text.x = element_text(angle = 90,
                        hjust = 1, vjust = 0.5),
        axis.text.y = element_text(angle = 0,
                        hjust = 1, vjust = 0.5),
        legend.position = "bottom",
        legend.key.width = unit(2.5, "cm")) +
    labs(x = NULL,
        y = NULL,
        fill = "Jaccard similarity")
})

  # Render plots
  output$plot1 <- renderPlot({
    data() %>%
    as.data.frame() %>%
    ggplot(aes(x = RT, y = Precursor.Quantity)) +
    geom_line(color = "darkblue", alpha = 0.7, show.legend = FALSE) +
    labs(
        x = "Retention time (min)",
        y = "Precursor quantity",
        color = NULL
    ) +
    facet_wrap(~Run, scales = "free_y")
  })
  
  output$plot2 <- renderPlot({
    data() %>%
    as.data.frame() %>%
    ggplot(aes(x = RT, y = Precursor.Mz)) +
    ggpointdensity::geom_pointdensity(size = 0.25) +
    viridis::scale_color_viridis(option = "plasma") +
    labs(
        x = "Retention time (min)",
        y = " Scan range (m/z)",
        color = NULL
    ) +
    theme(legend.position = "bottom",
        legend.key.width = unit(1.5, "cm")) +
    facet_wrap(~Run)
  })
  
  output$plot3 <- renderPlot({
    data() %>%
      as.data.frame() %>%
      ggplot(aes(x = Precursor.Mz, 
                 y = RT - Predicted.RT)) +
      ggpointdensity::geom_pointdensity(size = 0.25) +
      viridis::scale_color_viridis(option = "plasma") +
      geom_hline(yintercept = c(1,0,-1), linetype = "dashed", color = "black") +
      labs(
        x = "Precursor m/z",
        y = "RT - Predicted RT (min)",
        color = NULL
      ) +
      theme(legend.key.width = unit(2, "cm"),
            legend.position = "bottom") +
      facet_wrap(~Run)
  })
  
  output$plot4 <- renderPlot({
    data() %>%
    as.data.frame() %>%
    ggplot(aes(x = Precursor.Charge)) +
    geom_density(alpha = 0.7, 
            stat = "density", fill = "darkblue",
            show.legend = FALSE) +
    labs(
        x = "Precursor charge",
        y = "Density",
        fill = NULL
    ) +
    facet_wrap(~Run)
  })

  output$plot5 <- renderPlot({
    data() %>%
      ggplot(aes(x = peptide_length)) +
      geom_histogram(alpha = 0.7,
                     stat = "bin",
                     binwidth = 1,
                     show.legend = FALSE,
                     fill = "tomato") +
      labs(
        x = "Peptide length (a.a.)",
        y = "Count",
        fill = NULL
      ) +
      facet_wrap(~Run, scales = "free")
  })
  
  output$plot6 <- renderPlot({
    peptides_per_run() %>%
    as.data.frame() %>%
    ggplot(aes(y = Run, x = n_peptides)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.7,
            fill = "darkblue", show.legend = FALSE) +
    geom_text(aes(label = n_peptides),
        color = "white", size = 5,
        hjust = 1, nudge_x = -0.5
        ) +
    labs(y = NULL,
        x = "Number of peptides",
        fill = NULL) +
    theme(
      axis.text.x = element_text(angle = 90,
                        vjust = 0.5,
                        hjust = 1)
    )
  })

  output$plot7 <- renderPlot({
    proteins() %>%
    as.data.frame() %>%
    ggplot(aes(y = Run, x = n_proteins)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.7,
            fill = "darkblue", show.legend = FALSE) +
    geom_text(aes(label = n_proteins),
        color = "white", size = 5,
        hjust = 1, nudge_x = -0.5
        ) +
    labs(y = NULL,
        x = "Number of proteins",
        fill = NULL) +
    theme(
      axis.text.x = element_text(angle = 90,
                        vjust = 0.5,
                        hjust = 1)
    )
  })

  output$plot8 <- renderPlot({
    unique_genes() %>%
    as.data.frame() %>%
    gather(key = "Sample", value = "Intensity") %>%
    dplyr::mutate(missing = is.na(Intensity)) %>%
    dplyr::group_by(Sample) %>%
    dplyr::summarise(
        missing = mean(missing) * 100
    ) %>%
    ggplot(aes(x = Sample, y = missing)) +
    geom_col(fill = "darkblue", alpha = 0.7, show.legend = FALSE) +
    geom_text(aes(label = round(missing, 2)),
        vjust = -0.25, size = 3) +
    labs(x = NULL,
        y = "Proportion of missing values per sample (%)",
        fill = NULL) +
    theme(
        axis.text.x = element_text(angle = 90, 
                hjust = 1, vjust = 0.5)
    )
  })

  output$plot9 <- renderPlot({
    missing_vs_mean() %>%
    ggplot(aes(x = missing, y = median_intensity)) +
    geom_point(alpha = 0.7, size = 3) +
    geom_smooth(method = "lm", se = FALSE,
        color = "darkblue") +
    labs(x = "Proportion of missing values (%)",
        y = "Mean log2(abundance)",
        color = NULL) +
    facet_wrap(~norm)
  })

  output$plot10 <- renderPlot({
  combined_data() %>%
    as.data.frame() %>%
    ggplot(aes(
        x = Sample,
        y = Intensity,
        fill = norm
    )) +
    scale_fill_manual(values = c("Raw matrix" = "tomato",
                                "MAD normalised" = "darkblue")) +
    geom_boxplot(alpha = 0.7) +
    theme(
        axis.text.x = element_text(angle = 90, 
                hjust = 1, vjust = 0.5),
                legend.position = "none"
    ) +
    labs(x = NULL,
        y = "log2(Intensity)",
        fill = NULL) +
  facet_wrap(~norm)
  })
  
  output$plot11 <- renderPlot({
    data() %>%
    as.data.frame() %>%
     dplyr::mutate(specificity = case_when(
        str_detect(Stripped.Sequence, "K$|R$") ~ "Trypsin C-termini",
        str_detect(Stripped.Sequence, "E$|D$") ~ "GluC C-termini",
        TRUE ~ "Missed C-termini")
    ) %>%
    group_by(Run, specificity) %>%
    dplyr::summarise(
        peptides = n()
    ) %>%
    ggplot(aes(y = Run, x = peptides,
                fill = specificity)) +
    geom_col(alpha = 0.7,
            position = "stack") +
    scale_fill_manual(values = c("Trypsin C-termini" = "tomato",
                                "GluC C-termini" = "darkblue",
                                "Missed C-termini" = "darkgreen")) +
    labs(
        y = NULL,
        x = "Count",
        fill = "Specificity"
    ) +
    theme(legend.position = "top",
        axis.text.x = element_text(hjust = 0.5, vjust = 1))
  })

output$plot12 <- renderPlot({
    MS_corr() %>% 
    as.data.frame() %>%
    dplyr::mutate(EQScore_cutoff = case_when(
        .$Empirical.Quality >= input$Empirical.Quality ~ "Above threshold",
        TRUE ~ "Below threshold")
    ) %>%
    ggplot(aes(x = Ms1.Profile.Corr, fill = EQScore_cutoff)) +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c("tomato", "darkblue")) +
    labs(x = "MS1 Profile Correlation",
        y = "Density",
        fill = "Correlation between MS1 and MS2 chromatograms") +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 90)) +
    facet_wrap(~Run)
  })

output$plot13 <- renderPlot({
  QuantUMS_scores() %>%
    as.data.frame() %>%
    ggplot() +
    geom_density(aes(x = Score, fill = Filter), alpha = 0.7) +
    scale_fill_manual(values = c("PG.MaxLFQ.Quality" = "tomato",
                                "Empirical.Quality" = "darkblue",
                                "Quantity.Quality" = "darkgreen")) +
    labs(x = "Score",
        y = "Density",
        fill = NULL) +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 90)) +
    facet_wrap(~Run)
  })

output$plot14 <- renderPlot({
  gene_quantity() %>%
    as.data.frame() %>%
    ggplot() +
    geom_density(aes(x = log2(quantity), fill = protein_metrics), alpha = 0.7) +
    scale_fill_manual(values = c("Genes.MaxLFQ" = "tomato",
                                "Genes.MaxLFQ.Unique" = "darkblue",
                                "PG.MaxLFQ" = "darkgreen")) +
    labs(x = "Protein Quantity",
         y = "Density",
         fill = NULL) +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 90)) +
    facet_wrap(~Run)
})

output$QuantUMS_dist <- renderPlotly({
    data() %>%
    plot_ly(x = ~PG.MaxLFQ.Quality, 
            y = ~Quantity.Quality, 
            z = ~Empirical.Quality,
            color = ~Run,
            alpha = 0.6,
            colors = viridis::viridis(256),
            type = "scatter3d",
            mode = "markers") %>%
    layout(scene = list(
        xaxis = list(title = "PG MaxLFQ Quality"),
        yaxis = list(title = "Quantity Quality"),
        zaxis = list(title = "Empirical Quality")
    )) %>%
    subplot()
  })

# plot the sample correlation using the columns from the selectInput
output$Corr <- renderPlotly({
    unique_genes() %>%
    log2() %>%
    as.data.frame() %>%
    ggplot(aes(x = !!sym(input$xcol), y = !!sym(input$ycol))) +
    geom_point(alpha = 0.7, show.legend = FALSE) +
    geom_smooth(method = "lm", se = FALSE,
        color = "darkblue") +
    labs(x = paste0("Log2(", input$xcol, ")"),
        y = paste0("Log2(", input$ycol, ")"))
  })

# plot the PCA
output$PCA <- renderPlotly({
    pca_data()
  })

}

# Run the application
shinyApp(ui = ui, server = server)
