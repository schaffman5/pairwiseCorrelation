#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <math.h>
#include <omp.h>


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


void update_progress_bar(int current_step, int total_steps, int bar_width) {
  float progress = (float)current_step / total_steps;
  int pos = (int)(bar_width * progress);

  // start the line with carriage return to overwrite previous output
  Rprintf("\r[");

  // fill in the progress
  for (int i = 0; i < bar_width; ++i) {
    if (i < pos) Rprintf("=");
    else if (i == pos) Rprintf(">");
    else Rprintf(" ");
  }

  // end the progress bar with the percentage
  Rprintf("] %d%%", (int)(progress * 100));

  // flush the output buffer to ensure the progress bar is displayed immediately
  //fflush(stdout);

  R_FlushConsole();
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

  int j;
  int index;

  int dc = (*diag) ? 0 : 1; // when diag==1, include the diagonal by setting dc to zero
  int total = *nr + 1;
	int progress_threshold = 2; // update progress every 2%

	for(j = 0 ; j <= *nr ; j++) {
		R_CheckUserInterrupt();  // Check if R user wants to cancel

		// calculate percent complete
		int percent_complete = (int)(((double)j / total) * 100);

		// update progress display every progress_threshold percent
		if (percent_complete % progress_threshold == 0 && j % progress_threshold == 0) {
		  update_progress_bar(j, total, 50);
		}

    #pragma omp parallel for
		for(int i = j+dc ; i < *nr ; i++) {
		  // determine the index for the upper right triangle of the
		  // correlation matrix (dc==0 indicates if the diagonal is included)
		  index = (int)(j * (*nr - dc) + i - j * (j + 1) / 2) - dc;
		  d[index] = pearsonCor(x, *nr, *nc, i, j, *missingThresh);
		}
	}

	update_progress_bar(total, total, 50);
	Rprintf("\nComplete!\n");
}
