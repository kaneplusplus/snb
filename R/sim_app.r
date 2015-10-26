
#' Run the snb Simulation Web Application
#' 
#' This function launches a shiny web page for visualizing the stopped
#' negative binomial simulation in the browser.
#'
#' @param port The TCP port that the application should listen on. Defaults to port 8100.
#' @param launch.browser If true, the system's default web browser will be launched automatically after the app is started. Defaults to true in interactive sessions only.
#' @param workerId Can generally be ignored. Exists to help some editions of Shiny Server.
#' @importFrom shiny runApp
#' @export
snb_sim_app = function(port = 8100L,
         launch.browser = getOption("shiny.launch.browser", interactive()),
         workerId = "") {
  runApp( system.file("simApp", package='snm'), port=port, 
    launch.browser=launch.browser, workerId=workerId)
}

