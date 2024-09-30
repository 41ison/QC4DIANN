## Dashboard QC report for DIANN search results

# Load required libraries
library(shiny)
library(shinydashboard)
library(diann)
library(arrow)
library(tidyverse)
library(ggpointdensity)
library(limma)

# Increase the maximum filde size to 200 MB
options(shiny.maxRequestSize = 200 * 1024^2)

# Define UI for application that reads a parquet file and generates a QC report dashboard
ui <- dashboardPage(
  skin = "black",
  dashboardHeader(title = "QC Report for DIANN search results dashboard"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Filters", tabName = "filters", icon = icon("filter")),
      fileInput(inputId = "report", label = "Choose Parquet File", accept = ".parquet"),
      sliderInput("PG.MaxLFQ.Quality", "PG MaxLFQ Quality", min = 0, max = 1, value = 0.75),
      sliderInput("Empirical.Quality", "Empirical Quality", min = 0, max = 1, value = 0)
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "filters",
              fluidRow(
                box(title = "Reconstruction of XIC", status = "primary", solidHeader = TRUE, plotOutput("plot1"), collapsible = TRUE),
                box(title = "Density of ions", status = "primary", solidHeader = TRUE, plotOutput("plot2"), collapsible = TRUE),
                box(title = "Charge state distribution", status = "primary", solidHeader = TRUE, plotOutput("plot3"), collapsible = TRUE),
                box(title = "Proteins per sample", status = "primary", solidHeader = TRUE, plotOutput("plot4"), collapsible = TRUE),
                box(title = "Sparsity", status = "primary", solidHeader = TRUE, plotOutput("plot5"), collapsible = TRUE),
                box(title = "Missing vs median abundance", status = "primary", solidHeader = TRUE, plotOutput("plot6"), collapsible = TRUE),
                box(title = "Abundance", status = "primary", solidHeader = TRUE, plotOutput("plot7"), collapsible = TRUE),
                box(title = "Retention time error", status = "primary", solidHeader = TRUE, plotOutput("plot8"), collapsible = TRUE),
                box(title = "Missed cleavage", status = "primary", solidHeader = TRUE, plotOutput("plot9"), collapsible = TRUE),
                box(title = "Posterior Error Probability", status = "primary", solidHeader = TRUE, plotOutput("plot10"), collapsible = TRUE)
      )
    )
  )
)
)

# Define server logic required to read the parquet file and generate the QC report
server <- function(input, output) {
  theme_set(theme_bw())
  theme_update(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold", hjust = 0.5),
    legend.title.position = "top"
  )

  # Reactive expression to read and pre-process the uploaded parquet file
  data <- reactive({
    req(input$report)
    diann_report <- arrow::read_parquet(input$report$datapath) %>%
      dplyr::filter(Lib.PG.Q.Value <= 0.01 & Lib.Q.Value <= 0.01 & PG.Q.Value <= 0.01) %>%
      dplyr::mutate(File.Name = Run)
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

    # Reactive expression to filter matrix of protein abundance based on the input filters
    unique_genes <- reactive({
    req(data())
    diann::diann_matrix(data() %>%
    dplyr::filter(.$PG.MaxLFQ.Quality >= input$PG.MaxLFQ.Quality & .$Empirical.Quality >= input$Empirical.Quality),
      id.header = "Protein.Ids",
      quantity.header = "Genes.MaxLFQ.Unique",
      proteotypic.only = FALSE,
      pg.q = .01
    )
  })

  # Reactive expression to combine the raw and MAD normalised data
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

  # Render plots
  output$plot1 <- renderPlot({
    data() %>%
    as.data.frame() %>%
    dplyr::filter(.$PG.MaxLFQ.Quality >= input$PG.MaxLFQ.Quality & .$Empirical.Quality >= input$Empirical.Quality) %>%
    ggplot(aes(x = RT, y = Precursor.Quantity)) +
    geom_line(color = "darkblue", alpha = 0.7, show.legend = FALSE) +
    labs(
        x = "Retention time (min)",
        y = "Precursor quantity",
        color = NULL
    ) +
    facet_wrap(~Run)
  })
  
  output$plot2 <- renderPlot({
    data() %>%
    as.data.frame() %>%
    dplyr::filter(.$PG.MaxLFQ.Quality >= input$PG.MaxLFQ.Quality & .$Empirical.Quality >= input$Empirical.Quality) %>%
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
    dplyr::filter(.$PG.MaxLFQ.Quality >= input$PG.MaxLFQ.Quality & .$Empirical.Quality >= input$Empirical.Quality) %>%
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

  output$plot4 <- renderPlot({
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
    theme(axis.text.x = element_text(angle = 90,
                        vjust = 0.5,
                        hjust = 1)
    )
  })

  output$plot5 <- renderPlot({
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

  output$plot6 <- renderPlot({
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
    ggplot(aes(x = missing, y = median_intensity)) +
    geom_point(alpha = 0.7, size = 3) +
    geom_smooth(method = "lm", se = FALSE,
        color = "darkblue") +
    labs(x = "Proportion of missing values (%)",
        y = "Mean abundance",
        color = NULL)
  })

  output$plot7 <- renderPlot({
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

  output$plot8 <- renderPlot({
    data() %>%
    as.data.frame() %>%
    dplyr::filter(.$PG.MaxLFQ.Quality >= input$PG.MaxLFQ.Quality & .$Empirical.Quality >= input$Empirical.Quality) %>%
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
  })

  output$plot9 <- renderPlot({
    data() %>%
    as.data.frame() %>%
    dplyr::filter(.$PG.MaxLFQ.Quality >= input$PG.MaxLFQ.Quality & .$Empirical.Quality >= input$Empirical.Quality) %>%
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
    scale_fill_manual(values = c("Trypsin" = "darkblue",
                                "Unespecific" = "tomato")) +
    labs(
        x = NULL,
        y = "Count",
        fill = "Specificity"
    ) +
    theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90,
                        hjust = 1, vjust = 0.5))
  })

output$plot10 <- renderPlot({
  data() %>%
    as.data.frame() %>%
    dplyr::filter(.$PG.MaxLFQ.Quality >= input$PG.MaxLFQ.Quality & .$Empirical.Quality >= input$Empirical.Quality) %>%
    ggplot(aes(x = PEP)) +
    geom_density(fill = "darkblue", alpha = 0.7) +
    labs(
        x = "PEP",
        y = "Density",
        fill = NULL
    ) +
    theme(legend.position = "none") +
    facet_wrap(~Run)
  })

}

# Run the application
shinyApp(ui = ui, server = server)
