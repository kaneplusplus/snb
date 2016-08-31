library(shiny)
library(ggplot2)
library(snb)
library(reshape2)

shinyServer(function(input, output) {
  output$plot = renderPlot({
    s = input$s_param
    t = input$t_param
    p = input$p_param
    if (!is.na(s) && !is.na(t)) {
      x = min(s, t):(s+ t - 1)
      stack = as.data.frame(dsnb_stacked(x, p, s, t))
      stacked_plot(stack$x, stack$s, stack$t)
    }
  })

  output$endpoint_probs = renderDataTable({
    s = input$s_param 
    t = input$t_param 
    density = as.data.frame(dsnb_stacked(min(s,t):(s+t-1), input$p_param, s, t))
    data.frame(list(Outcome=c("Success", "Failure"), 
                    Probability=c(sum(density$s), sum(density$t))))
  }, options=list(searching=FALSE, paging=FALSE))

  output$tab1 = renderDataTable({
    s = input$s_param 
    t = input$t_param 
    density = as.data.frame(dsnb_stacked(min(s,t):(s+t-1), input$p_param, s, t))
    density$TOC = density[,2] + density[,3]
    density[,2] = signif(density[,2], 4)
    density[,3] = signif(density[,3], 4)
    density[,4] = signif(density[,4], 4)
    names(density)=c("Enrollment Number", "Trial Success", "Trial Failure",
                     "Total Outcome Prob.")
    density
    }, 
    options = list(pageLength = 25)
  )

})

