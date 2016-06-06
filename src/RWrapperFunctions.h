

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/Error.h>
#include <cmath>
#include <vector>
#include "DebugConstants.h"
#include <Rembedded.h>
#include <Rdefines.h>

double*** RGETLIST(SEXP sexpobs, int* T, int nsample, int D);

std::vector<std::vector<double> >  RGETMAT(SEXP sexpA, int K, int N);
std::vector<std::vector<int> >  RGETMATINT(SEXP sexpA, int K, int N);
std::vector<std::vector<std::vector<double> > >  RGETMAT(SEXP sexpA, int K, int N, int M);
std::vector<double>  RGETVECT(SEXP sexpA, int K);
std::vector<int>  RGETVECTINT(SEXP sexpA, int K);
SEXP RPREPAREMAT(std::vector<std::vector<double> > transMat);
SEXP RPREPAREMAT(std::vector<std::vector<std::vector<double> > > transMat);
SEXP RPREPAREVECT(std::vector<double> v);
SEXP RPREPAREVECT(std::vector<int> v) ;
SEXP RPREPAREMAT(std::vector<std::vector<int> > transMat);


