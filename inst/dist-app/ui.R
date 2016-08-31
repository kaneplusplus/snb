require(shiny)

shinyUI(pageWithSidebar(

  headerPanel("Stopped Negative Binomial Distribution"),

  sidebarPanel(
    sliderInput('p_param', 'Responder Probability', min=0, max=1,
                value=0.4, step=0.01),
    numericInput("s_param", "Total Responders", value=4, min=1, step=1),
    numericInput("t_param", "Total Non-responders", value=6, min=1, step=1)
  ),

  mainPanel(
    h3("Density Plot"),
    plotOutput('plot'),
    h3("Endpoint Probabilities"),
    dataTableOutput("endpoint_probs"),
    h3("Density Values"),
    dataTableOutput("tab1")
  )
))

