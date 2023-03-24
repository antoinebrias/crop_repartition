#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for crop repartition application
fluidPage(

  titlePanel("Crop repartition"),
  sidebarLayout(
    sidebarPanel(
      numericInput("n_modalites", "Number of crop modalities", value = 6, min = 1, max = 20, step = 1),
      selectInput("str_block_shape", "Block shape (rows x columns)", c("")),
      numericInput("n_rep", "Number of replicates", value = 2, min = 1, max = 20, step = 1),
      selectInput("str_plot_shape", "Plot shape (rows x columns)", c("")),
      selectInput("constraint", "Constraint", c("None", "Orthogonal","Ortho+Diag","Orthogonal 2","Ortho2 + Diag2"),  selected = "Ortho+Diag",),
      actionButton("calculate", "Calculer"),
      downloadButton("downloadCSV", "Exporter en CSV")
    ),
    mainPanel(
      plotOutput("repartition_plot")
    )
  )
)
