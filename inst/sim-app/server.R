require(shiny)
require(ggplot2)
require(snm)

flips <- snbFlips(1, 0.2, 3, 20)

shinyServer(function(input, output) {

  flips <- reactive({
    input$simButton
    if (!is.na(input$sParam) && !is.na(input$tParam)) {
      snbFlips(1, input$pParam, input$sParam, input$tParam)
    }
  })

  output$zplot <- renderPlot({
    if (!is.na(input$sParam) && !is.na(input$tParam)) {
      print(zplot(flips(), input$sParam, input$tParam))
    }
  }, height=400)
  
  output$kplot <- renderPlot({
    if (!is.na(input$sParam) && !is.na(input$tParam)) {
      print(kplot(flips(), input$sParam, input$tParam))
    }
  }, height=400)
 
})
