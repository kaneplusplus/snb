\name{runSnbDistApp}
\alias{runSnbDistApp}
\title{Run the snb Desity Web Application}
\usage{
  runSnbDistApp(port = 8100L,
    launch.browser = getOption("shiny.launch.browser", interactive()),
    workerId = "")

  runSnbDistApp(port = 8100L,
    launch.browser = getOption("shiny.launch.browser", interactive()),
    workerId = "")
}
\arguments{
  \item{port}{The TCP port that the application should
  listen on. Defaults to port 8100.}

  \item{launch.browser}{If true, the system's default web
  browser will be launched automatically after the app is
  started. Defaults to true in interactive sessions only.}

  \item{workerId:}{Can generally be ignored. Exists to help
  some editions of Shiny Serv}

  \item{port}{The TCP port that the application should
  listen on. Defaults to port 8100.}

  \item{launch.browser}{If true, the system's default web
  browser will be launched automatically after the app is
  started. Defaults to true in interactive sessions only.}

  \item{workerId:}{Can generally be ignored. Exists to help
  some editions of Shiny Serv}
}
\description{
  This function launches a shiny web page for visualizing
  the stopped negative binomial distribution in the
  browser.

  This function launches a shiny web page for visualizing
  the stopped negative binomial distribution in the
  browser.
}

