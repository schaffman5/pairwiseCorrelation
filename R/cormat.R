#' Matrix pairwise correlation
#'
#' @description
#' This function will take a matrix from R and do all pairwise
#' correlation calculations for each row in the matrix and
#' return a distance object containing the Pearson correlations.
#'
#' The values returned should be equivalent to running
#' cor(x,y,method="pearson",use="pair") where x and y are the
#' data from two rows being compared.
#'
#' In addition, if the fraction of missing pair data between
#' the two matrix rows is greater than the missingThresh parameter,
#' the correlation is not calculated and an NA value is returned.
#'
#' 2006/02/22 (Mike Schaffer)
#'
#' @param x numeric matrix where rows are variables to compare (e.g. genes) and
#' columns are conditions or samples.
#' @param diag logical to compute the diagonal elements.
#' @param upper logical to show upper part of the matrix.
#' @param missingThresh fraction of paired data that is required for the
#' correlation to be calculated.
#'
#' @returns An R dist object containing the Pearson correlation for each
# 	row comparison.
#'
#' @useDynLib pairwiseCorrelation, .registration=TRUE
#' @export
#'
#' @examples
#'
#' # compute pairwise Pearson correlations across rows
#' cormat(rbind(c(2.5,4.2,1.6,6.2,3.1), c(5.2,3.5,8.2,7.4,5.3), c(0.3, 0.4, 1.4, 2.7, 3.1) ))
#'
#' # compare cormat output to base R cor output:
#' cormat(rbind(c(1,1,3,4,NA),c(1,2,3,4,5) ) )  # 0.9467293
#' cor(c(1,1,3,4,NA),c(1,2,3,4,5),use="pair")   # 0.9467293
#'
cormat <- function (x, diag = FALSE, upper = FALSE, missingThresh=0.50) {
    N <- nrow(x <- as.matrix(x))
    d <- .C("Rcormat", x = as.double(x), nr = N, nc = ncol(x),
        d = double(N * (N - 1)/2), diag = as.integer(FALSE),
        missingThresh = as.double(missingThresh),
        NAOK = TRUE)$d
    attr(d, "Size") <- N
    attr(d, "Labels") <- dimnames(x)[[1]]
    attr(d, "Diag") <- diag
    attr(d, "Upper") <- upper
    attr(d, "method") <- "pearson"
    attr(d, "call") <- match.call()
    class(d) <- "dist"
    return(d)
}


#' Extract value from dist object
#'
#' @description
#' Convenience function to extract a distance (or correlation) from a dist object
#'
#' @param i row index of dist object
#' @param j column index of dist object
#' @param dis dist object
#'
#' @returns distance value
#' @export
#'
#' @examples
#' corrs<-cormat(rbind(c(2.5,4.2,1.6,6.2,3.1), c(5.2,3.5,8.2,7.4,5.3), c(0.3, 0.4, 1.4, 2.7, 3.1) ))
#' get.dist(1, 1, corrs)
#'
get.dist<-function(i, j, dis) {
	if(i==j) {return(1)}
	r<-range(c(i,j))
	i<-r[1]
	j<-r[2]
	dis[attributes(dis)$Size*(i-1) - i*(i-1)/2 + j-i]
}


#' Extract values from dist object for a given element
#'
#' @description
#' Convenience function to extract distances (or correlations) from a dist object
#'
#' @param row index of row in the dist object
#' @param dis dist object
#'
#' @returns vector of distance values
#' @export
#'
#' @examples
#' corrs<-cormat(rbind(c(2.5,4.2,1.6,6.2,3.1), c(5.2,3.5,8.2,7.4,5.3), c(0.3, 0.4, 1.4, 2.7, 3.1) ))
#' get.all.dists(1, corrs)
#'
get.all.dists<-function(row, dis) {
	dat<-unlist(sapply(c(1:attributes(dis)$Size),get.dist,j=row,dis=dis,simplify=F))
	names(dat)<-attributes(dis)$Labels
	dat
}
