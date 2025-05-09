% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cormat.R
\name{get.dist}
\alias{get.dist}
\title{Extract value from dist object}
\usage{
get.dist(i, j, dis)
}
\arguments{
\item{i}{row index of dist object}

\item{j}{column index of dist object}

\item{dis}{dist object}
}
\value{
distance value
}
\description{
Convenience function to extract a distance (or correlation) from a dist object
}
\examples{
corrs<-cormat(rbind(
c(2.5,4.2,1.6,6.2,3.1),
c(5.2,3.5,8.2,7.4,5.3),
c(0.3, 0.4, 1.4, 2.7, 3.1),
c(0.4, 0.4, 1.2, 2.7, 3.1),
c(0.5, 0.4, 1.4, 2.4, 3.2)
))
get.dist(3, 4, corrs)

}
