require(shiny)
require(ggplot2)

shinyUI(pageWithSidebar(
  
  headerPanel("Stopped Negative Binomial Simulation"),
  
  sidebarPanel(
    sliderInput('pParam', 'p', min=0, max=1,
                value=0.2, step=0.01),
    numericInput("sParam", "s", value=3, min=1, step=1),
    numericInput("tParam", "t", value=20, min=1, step=1),
    actionButton("simButton", "Next Simulation")
  ),
  
  mainPanel(
    h3("Zelterman Plot"),
    plotOutput('zplot'),
    h3("Kane Plot"),
    plotOutput('kplot')
  )
))
