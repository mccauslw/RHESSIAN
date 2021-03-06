% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{spline_draw_cpp}
\alias{spline_draw_cpp}
\title{Call C language spline_draw_c function}
\usage{
spline_draw_cpp(p, m, n_draws)
}
\arguments{
\item{p}{A vector of levels on grid}

\item{m}{A vector of derivatives on grid}

\item{n_draws}{The number of draws to make from the spline density}
}
\value{
The vector of spline density draws
}
\description{
Call C language spline_draw_c function
}
