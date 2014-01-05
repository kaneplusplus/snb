require(shiny)
require(ggplot2)
require(snm)

shinyServer(function(input, output) {
  output$plot <- renderPlot({
    if (!is.na(input$sParam) && !is.na(input$tParam)) {
      x <- input$sParam:(input$tParam+input$sParam-1)
      p <- qplot(factor(x), dsnb(x, input$pParam, input$sParam, input$tParam),
        stat="identity", geom="bar", ylab="f(t)", xlab="t")
      print(p)
    }
  }, height=600)
 
})
