% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reaction_type.R
\name{dts}
\alias{dts}
\title{Get reaction data from a time, and sort it by frequency of occurrence
NB: it is currently assumed that the data is in timesteps of 20000}
\usage{
dts(time, rd, sr = T)
}
\arguments{
\item{time}{the time}

\item{rd}{the data}

\item{plax}{whether to plot the bottom axis}

\item{mymar}{the figure margins}

\item{cols}{the plot colours (defaults to orange and blue)}

\item{rno}{the run number}
}
\description{
Get reaction data from a time, and sort it by frequency of occurrence
NB: it is currently assumed that the data is in timesteps of 20000
}
