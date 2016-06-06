#include <iostream>
#include "RWrapperFunctions.h"
#include <stdio.h>
//

using namespace std;

double*** RGETLIST(SEXP sexpobs, int* T, int nsample, int D) {
    int t,d,n;
    double ***observations = NULL;
    
    if(nsample > 0) {
        int mem = 0;
        observations = (double***)malloc(sizeof(double**)*nsample);
        mem += sizeof(double**)*nsample;
        for(n=0; n<nsample; n++) {
            observations[n] = (double**)malloc(sizeof(double*)*T[n]);
            mem += sizeof(double*)*T[n];
            for(t=0; t<T[n]; t++) {
                observations[n][t] = (double*)malloc(sizeof(double)*D);
                mem += sizeof(double)*D;
                for(d=0; d<D; d++) {
                    observations[n][t][d] = REAL(coerceVector(VECTOR_ELT(sexpobs, n), REALSXP))[t+T[n]*d];
                    //	Rprintf("%f ", observations[n][t][d]);
                }
                //Rprintf("\n");
            }
            //Rprintf("\n\n");
        }
    }
    return observations;
}



std::vector<std::vector<std::vector<double> > >  RGETMAT(SEXP sexpA, int K, int N, int M) {
    int i,j,r;
    std::vector<std::vector<std::vector<double> > > matrix(K, std::vector< std::vector<double> >(N, std::vector<double> (M)));
    //std::vector<std::vector<std::vector<double> > >matrix( K, std::vector<double>(N,std::vector<double>(M)));
    sexpA = coerceVector(sexpA, REALSXP);
    /*
    1,2,3,4,5,6,7,8
    1,1,1,1,2,2,2,2
    0,1,2,3,4,5,6,7
     */
    for(i=0; i<K; i++){
        for(j=0; j<N; j++){
            for(r=0; r<M; r++){
                matrix[i][j][r] = REAL(sexpA)[i*(K*N)+K*j+r];
            }
        }
    }
    return matrix;
}

std::vector<std::vector<double> >  RGETMAT(SEXP sexpA, int K, int N) {
    int i,j;
    std::vector<std::vector<double> > matrix( K, std::vector<double>(N));
    sexpA = coerceVector(sexpA, REALSXP);
    for(i=0; i<K; i++){
        for(j=0; j<N; j++){
            matrix[i][j] = REAL(sexpA)[i+K*j];
        }
    }
    return matrix;
}
std::vector<std::vector<int> >  RGETMATINT(SEXP sexpA, int K, int N) {
    int i,j;
    std::vector<std::vector<int> > matrix( K, std::vector<int>(N));
    sexpA = coerceVector(sexpA, REALSXP);
    for(i=0; i<K; i++){
        for(j=0; j<N; j++){
            matrix[i][j] = REAL(sexpA)[i+K*j];
        }
    }
    
    return matrix;
}
std::vector<double>  RGETVECT(SEXP sexpA, int K) {
    int i;
    std::vector<double> matrix(K);
    for(i=0; i<K; i++){
        matrix[i] = REAL(sexpA)[i];
    }
    
    return matrix;
}
std::vector<int>  RGETVECTINT(SEXP sexpA, int K) {
    int i;
    std::vector<int> matrix(K);
    for(i=0; i<K; i++){
        matrix[i] = INTEGER(sexpA)[i];
    }
    
    return matrix;
}
/////////////////////////////
SEXP RPREPAREMAT(std::vector<std::vector<double> > transMat) {
    SEXP AFit;
    int i,j;
    int K=transMat.size();
    PROTECT(AFit = NEW_NUMERIC(transMat.size()*transMat.at(0).size()));
    int counter=0;
    for(i=0; i<transMat.size(); i++){
        for(j=0; j<transMat.at(0).size(); j++){
            //cout<<"i: "<<i<<"j: "<<j<<" : "<<j+K*i<<" : "<<transMat[i][j]<<endl;
            NUMERIC_POINTER(AFit)[counter] = transMat[i][j];
            counter++;
        }
    }
    UNPROTECT(1);
    return AFit;
}
SEXP RPREPAREMAT(std::vector<std::vector<int> > transMat) {
    SEXP AFit;
    int i,j;
    int K=transMat.size();
    PROTECT(AFit = NEW_NUMERIC(transMat.size()*transMat.at(0).size()));
    int counter=0;
    for(i=0; i<transMat.size(); i++){
        for(j=0; j<transMat.at(0).size(); j++){
            //cout<<"i: "<<i<<"j: "<<j<<" : "<<j+K*i<<" : "<<transMat[i][j]<<endl;
            NUMERIC_POINTER(AFit)[counter] = transMat[i][j];
            counter++;
        }
    }
    UNPROTECT(1);
    return AFit;
}
/////////////////////////////
SEXP RPREPAREMAT(std::vector<std::vector<std::vector<double> > > transMat) {
    SEXP AFit;
    int r,i,j;
    int l=transMat.size();
    int K=transMat.size();
    PROTECT(AFit = NEW_NUMERIC(transMat.size()*transMat.at(0).size()*transMat.at(0).at(0).size()));
    int counter=0;
    for(r=0; r<transMat.size(); r++){
        for(i=0; i<transMat.at(0).size(); i++){
            for(j=0; j<transMat.at(0).size(); j++){
                NUMERIC_POINTER(AFit)[counter] = transMat[r][i][j];
                counter++;
            }
        }
    }
    UNPROTECT(1);
    return AFit;
}
/////////////////////////////
SEXP RPREPAREVECT(std::vector<double> v) {
    SEXP AFit;
    int i,j;
    PROTECT(AFit = NEW_NUMERIC(v.size()));
    for(i=0; i<v.size(); i++){
            NUMERIC_POINTER(AFit)[i] = v[i];
    }
    UNPROTECT(1);
    return AFit;
}
/////////////////////////////
SEXP RPREPAREVECT(std::vector<int> v) {
    SEXP AFit;
    int i,j;
    PROTECT(AFit = NEW_NUMERIC(v.size()));
    for(i=0; i<v.size(); i++){
        NUMERIC_POINTER(AFit)[i] = v[i];
    }
    UNPROTECT(1);
    return AFit;
}
