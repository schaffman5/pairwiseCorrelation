/* cormat.c: Pairwise Matrix Row Correlation

This function will take a matrix from R and do all pairwise
correlation calculations for each row in the matrix and
return a vector containing the Pearson correlations.

The values returned should be equivalent to running
cor(x,y,method="pearson",use="pair") where x and y are the
data from two rows being compared.

In addition, if the fraction of missing pair data between
the two matrix rows is greater than the missingThresh parameter,
the correlation is not calculated and an NA value is returned.

2006/02/22 (Mike Schaffer)

To compile: from the command line, run:
	R CMD SHLIB cormat.c


To use: from an R script:
	source("cormat.R")
	cormat(testset, diag = FALSE, upper = FALSE, missingThresh=0.50)


Inputs:
	x: a numeric matrix where rows are variables to compare (e.g.
		genes) and columns are conditions.
	diag: a boolean whether the diagonal elements should be
		calculated.
	upper: whether the upper part of the matrix should be
		shown.
	missingThresh: fraction of paired data that is required for
		the correlation to be calculated.


Output:
	An R vector containing the Pearson correlation for each
	row comparison.  It is left up to the user to determine which vector
	elements correspond to a given pair.  A distance object was previously
	returned, but for very large matrices, R was running out of memory.


Example:
	cormat(rbind(c(1,2,3,4),c(1,2,3,4),c(1,2,3,4)) )

	# Test this function
	cormat(rbind(c(1,1,3,4,NA),c(1,2,3,4,5)) )
	# Compare to R's builtin cor function
	cor(c(1,1,3,4,NA),c(1,2,3,4,5),use="pair")
*/

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <math.h>


// Pearson function
// 2006/02/24 MES

double pearsonCor (double *x, int nr, int nc, int i1, int i2, double missingThresh) {
	int j;

	double n=0.0, sY1=0.0, sprodY1Y1=0.0, sY2=0.0;
	double sprodY2Y2=0.0, sprodY1Y2=0.0;
	double sprod=0.0, ssqY1=0.0, ssqY2=0.0;
	double cor;

	for(j=0; j < nc; j++) {
		if(!ISNAN(x[i1 + nr*j]) && !ISNAN(x[i2 + nr*j])) {
			n++;

			sY1 = sY1 + x[i1 + nr*j];
			sprodY1Y1 = sprodY1Y1 + (x[i1 + nr*j]*x[i1 + nr*j]);
			sY2 = sY2 + x[i2 + nr*j];
			sprodY2Y2 = sprodY2Y2 + (x[i2 + nr*j]*x[i2 + nr*j]);
			sprodY1Y2 = sprodY1Y2 + (x[i1 + nr*j]*x[i2 + nr*j]);
		}
	}

	if( (n >= 2) && ( (n/nc) >= missingThresh) ) {
		ssqY1 = sprodY1Y1 - ((sY1*sY1)/n);
		ssqY2 = sprodY2Y2 - ((sY2*sY2)/n);
		sprod = sprodY1Y2 - ((sY1*sY2)/n);
		cor   = sprod/sqrt(ssqY1*ssqY2);
		return cor;

	} else {
		return NA_REAL;
	}
}


// Notes:
// Code to run Pearson on all pairs of genes on an incoming matrix where
// rows are genes and columns are conditions.  The matrix is dealt with as
// a vector in this function and in the Pearson function, so the dimensions
// need to be passed back and forth.  The benefit is speed and not having to
// make copies of the data in memory.  Additionally,the results are returned
// as a vector as R seems to be running out of memory for large matrices when
// Henry was trying to compare 17,000 genes x ~600 conditions.  It is thus
// left up to the user to determine which pair corresponds to which vector
// element (not hard, but requires some thinking as only half of the full
// matrix is returned (it's symmetric, so no need to do the full matrix or
// the diagonal -- unless the diag variable is set to 1 (TRUE)).

void Rcormat(double *x, int *nr, int *nc, double *d, int *diag, double *missingThresh)
{
	int dc, i, j, ij;
	int increment = 10;  // Give user feedback on the progress
					 // in increments of n rows.

	if(*nr > increment) {
		Rprintf("Processing row...\n");
	}

	dc = (*diag) ? 0 : 1; /* diag=1:  we do the diagonal */
	ij = 0;
	for(j = 0 ; j <= *nr ; j++) {

		R_CheckUserInterrupt();  // Check if R user wants to cancel

		if( (*nr > increment) && (( (j+1) % increment) == 0) ) {
			Rprintf("[%d] ",(j+1));
		}

		for(i = j+dc ; i < *nr ; i++) {
			d[ij++] = pearsonCor(x, *nr, *nc, i, j, *missingThresh);
		}
	}

	if(*nr > increment) {
		Rprintf("\nProcessing complete!\n");
	}
}
