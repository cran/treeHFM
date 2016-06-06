//
//  MatlabWrapperBM.cpp
//  FGBMA
//
//  Created by Henrik Failmezger on 16.01.15.
//  Copyright (c) 2015 Henrik Failmezger. All rights reserved.
//

#include <stdio.h>
#include "BaumWelch.h"
#include "RWrapperFunctions.h"


extern "C" {
    
SEXP RWrapperBM( SEXP sexpstartProbN, SEXP  sexptransProbSeqN, SEXP sexptransProbDiv, SEXP sexpemProbN, SEXP sexpobservationSequences, SEXP sexpnodeIndices,SEXP sexpparents, SEXP sexpsummationIndices, SEXP sexpindicesXD,SEXP sexptype, SEXP sexpstateIndicesSingle,SEXP sexpnHStates,SEXP sexpnDStates,SEXP sexpnObs,SEXP sexpdimData,SEXP sexpnEStates,SEXP sexpsigma,SEXP sexpmu){
   
    int nHStates=INTEGER(sexpnHStates)[0];
    int nDStates=INTEGER(sexpnDStates)[0];
    int nEStates=INTEGER(sexpnEStates)[0];
    int nObs=INTEGER(sexpnObs)[0];
    int dimData=INTEGER(sexpdimData)[0];
    //Rprintf("Type\n");
    const char *typeC = CHAR(STRING_ELT(sexptype, 0));
    char type=typeC[0];
    //////////////////////// TYPE //////////////////////////////////////////////
    ////////////////////////////////////////////////////////
 
    //transProbSeqM///////////////////////////////////////////////////////
    //Rprintf("startProbN1\n");
    std::vector<std::vector<double> > startProbN=  RGETMAT(sexpstartProbN, nHStates, nHStates);
    //Rprintf("startProbN2\n");
    //printVector(startProbN);
    //Rprintf("startProbN3\n");
    //printVector(startProbN);
    ///////////////////////////////////////////////////////////////////////
    std::vector<std::vector<double> > transProbSeqN= RGETMAT(sexptransProbSeqN, nHStates, nHStates);
    //Rprintf("transProbSeqN\n");
    //////////transProbDiv/////////////////////////////////////////////////
    //cout<<"nDStates: "<<nDStates<<endl;
    std::vector<std::vector<double> > transProbDiv= RGETMAT(sexptransProbDiv, nHStates, nDStates);
    //Rprintf("transProbDiv\n");
    ///////////////////////////////////////////////////////////////////////
    ///////////probO///////////////////////////////////////////////////////
    //cout<<"nObs: "<<nObs<<endl;
    //cout<<"nEStates: "<<nEStates<<endl;
    std::vector<std::vector <double> > emProbN= RGETMAT(sexpemProbN, nHStates, nEStates);
    //printVector(emProbN);
     //Rprintf("emProbN\n");
    ///////////////////////////////////////////////////////////////////////
    
    //Rprintf("observation sequencesN \n");
    
    ///////////////////////////////////////////////////////////////////////
    
    ////////////node indices//////////////////////////////////////
    std::vector <nodeIndices> allNodeIndices(nObs);
    for(int i=0;i<nObs;i++){
        //cout<<"i: "<<i<<endl;
        SEXP g = coerceVector(VECTOR_ELT(sexpnodeIndices, i),INTSXP);
        //cout<<"g: "<<g<<endl;
        std::vector <int>   observationNodeIndices = RGETVECTINT(g,  length(g));
        nodeIndices obsNodeIndices;
        
        //printVector(observationNodeIndices);
        obsNodeIndices.data=observationNodeIndices;
        allNodeIndices.at(i)=obsNodeIndices;
    }
    //printVector(allNodeIndices.at(0).data);
    ////////////parents//////////////////////////////////////
    //Rprintf("parentsN: ");
    std::vector <parentIndices> allParents(nObs);
    for(int i=0;i<nObs;i++){
        SEXP g = coerceVector(VECTOR_ELT(sexpparents, i),INTSXP);
        std::vector <int>   observationParent = RGETVECTINT(g,  length(g));
        parentIndices obsParent;
        obsParent.data=observationParent;
        allParents.at(i)=obsParent;
    }
    //printVector(allParents.at(0).data);
    ///////////////////////////////////////////////////////////////////////
    ////////////observation sequences//////////////////////////////////////
    std::vector <observationSequence> allObservationSequences(nObs);
    std::vector <observationSequenceDiscrete> allObservationSequencesDiscrete(nObs);
    int dimensionData=0;
    if(type=='c'){
        for(int k=0;k<nObs;k++){
            //cout<<"k: "<<k<<endl;
            SEXP g = VECTOR_ELT(sexpobservationSequences, k);
            SEXP Rdim=getAttrib(g,R_DimSymbol);
            int IOS=INTEGER(Rdim)[0];
            int JOS=INTEGER(Rdim)[1];
            observationSequence obsSequence;
            //cout<<"Rdim: IOS: "<<IOS<<"JOS: "<<JOS<<endl;
            std::vector<std::vector <double> >  observationMatrix= RGETMAT(g, IOS, JOS);
            //printVector(observationMatrix);
            //
            //cout<<"observationMatrix:"<<endl;
            //printVector(observationMatrix);
            dimensionData=observationMatrix.size();
            obsSequence.data=observationMatrix;
            allObservationSequences.at(k)=obsSequence;
        }
    }else{
        //Rprintf("observation sequences discrete\n");
        //cout<<"nObs: "<<nObs<<endl;
        for(int k=0;k<nObs;k++){
            //cout<<"k: "<<k<<endl;
            //SEXP g = VECTOR_ELT(sexpobservationSequences, k);
            SEXP g = coerceVector(VECTOR_ELT(sexpobservationSequences, k),INTSXP);
            //cout<<"g: "<<k<<endl;
            //int t = INTEGER(getAttrib(coerceVector(VECTOR_ELT(sexpobservationSequences, k), INTSXP), R_DimSymbol))[0];
            //cout<<"t: "<<t<<endl;
            //SEXP Rdim=getAttrib(g,R_DimSymbol);
            //cout<<"Rdim: "<<k<<endl;
            //int IP=INTEGER(Rdim)[0];
            observationSequenceDiscrete obsSequence;
            //cout<<"observationMatrix: "<<endl;
            std::vector<int>   observationMatrix= RGETVECTINT(g, length(g));
            //
            dimensionData=observationMatrix.size();
            obsSequence.data=observationMatrix;
            allObservationSequencesDiscrete.at(k)=obsSequence;
        }
    }
    //Rprintf("insert: \n");
    ///////////////////////////////////////////////////////////////////////
    ////////////summation indices//////////////////////////////////////
    //Rprintf("summation indices: \n");
    std::vector <std::vector <int> > summationIndices= RGETMATINT(sexpsummationIndices, nHStates, nHStates);
    //Rprintf("subEl: \n");
    summationIndices=subEl(summationIndices,1);
    //printVector(summationIndices);
    ///////////////////////////////////////////////////////////////////////
    ////////////indices XD//////////////////////////////////////
    //Rprintf("indices XD: \n");
    //Rprintf("indices XD2: \n");
    R_len_t IXD = length(sexpindicesXD);
    //INTEGER(RdimXD)[1];
    //cout<<"indices XD: "<<IXD<<endl;
    std::vector <int> indicesXD=RGETVECTINT(sexpindicesXD, IXD);
    indicesXD= subEl(indicesXD,1); //bring them to c++ indices
    //////////// OBS PROB FOR CONTINOUS ///////////////////////////////////
    
    //////////////////////// StateIndicesSingle //////////////////////////////////////////////
    //Rprintf("State indices2: \n");
    SEXP RdimSI=getAttrib(sexpstateIndicesSingle,R_DimSymbol);
    int ISI=INTEGER(RdimSI)[0];
    //cout<<"ISI: "<<ISI<<endl;
    int JSI=INTEGER(RdimSI)[1];
    //cout<<"JSI: "<<JSI<<endl;
    std::vector<std::vector <int> > stateIndicesSingle= RGETMATINT(sexpstateIndicesSingle, ISI, JSI);
    //printVector(stateIndicesSingle);
    //////////////////////// SigmaInit //////////////////////////////////////////////
    std::vector<std::vector<std::vector<double> > > SigmaInit=RGETMAT(sexpsigma,nHStates, dimData, dimData);
    //
    //cout<<"---------------------------------------------------------------"<<endl;
    //printVector(SigmaInit.at(0));
    //cout<<"--------"<<endl;
    //printVector(SigmaInit.at(0));
    //cout<<"---------------------------------------------------------------"<<endl;
    /*
    if(type=='c'){
        for(int i=0;i<nHStates;i++){
            cout<<"i: "<<i<<endl;
            SEXP g = VECTOR_ELT(sexpsigma, i);
            SigmaInit.at(i) = RGETMAT(g, nHStates, nHStates);
        }
    }
     */
    //
    
    ///////////////////////////////////////////////////////////////////////
    //////////////////////// MuInit //////////////////////////////////////////////
    //Rprintf("mu init: \n");
    std::vector<std::vector<double> > muN;
    if(type=='c'){
        muN=RGETMAT(sexpmu, dimData, nHStates);
    }
    //
    ///////////////////////////////////////////////////////////////////////
    //Rprintf("type: \n");
    SEXP HFM, wnames;
    if(type=='c'){
        PROTECT(HFM = NEW_LIST(6));
        PROTECT(wnames = NEW_CHARACTER(6));
        SET_STRING_ELT(wnames, 0, mkChar("transMatSeq"));
        SET_STRING_ELT(wnames, 1, mkChar("transMatDiv"));
        SET_STRING_ELT(wnames, 2, mkChar("initProb"));
        SET_STRING_ELT(wnames, 3, mkChar("loglik"));
        SET_STRING_ELT(wnames, 4, mkChar("mu"));
        SET_STRING_ELT(wnames, 5, mkChar("sigma"));
        SET_NAMES(HFM, wnames);
    }else{
        PROTECT(HFM = NEW_LIST(5));
        PROTECT(wnames = NEW_CHARACTER(5));
        SET_STRING_ELT(wnames, 0, mkChar("transMatSeq"));
        SET_STRING_ELT(wnames, 1, mkChar("transMatDiv"));
        SET_STRING_ELT(wnames, 2, mkChar("initProb"));
        SET_STRING_ELT(wnames, 3, mkChar("loglik"));
        SET_STRING_ELT(wnames, 4, mkChar("emission"));
        SET_NAMES(HFM, wnames);
    }
    
    
   
    //
    SEXP transProSeqbOut;
    SEXP transProbDivOut;
    SEXP startProbOut;
    SEXP llOut;
    if(type=='c'){
        outputC bmOutputC=baumWelch(allObservationSequences,  startProbN, transProbSeqN,  transProbDiv,  emProbN, muN, SigmaInit,  allNodeIndices, summationIndices, stateIndicesSingle, indicesXD,  allParents,  type);
        ///
        
        ///
        transProSeqbOut = RPREPAREMAT(bmOutputC.transProbSeq);
        SET_ELEMENT(HFM, 0, transProSeqbOut);
        ///
        ///
        //Rprintf("transProbDiv: \n");
        transProbDivOut = RPREPAREMAT(bmOutputC.transProbDiv);
        SET_ELEMENT(HFM, 1, transProbDivOut);
        ///
        startProbOut = RPREPAREVECT(bmOutputC.prior.at(0));
        SET_ELEMENT(HFM, 2, startProbOut);
        ///
        llOut = RPREPAREVECT(bmOutputC.allLL);
        SET_ELEMENT(HFM, 3, llOut);
        ///
        ///
        //Rprintf("muN: \n");
        //printVector(bmOutputC.mu);
        SEXP muC = RPREPAREMAT(bmOutputC.mu);
        SET_ELEMENT(HFM, 4, muC);
        ///
        SEXP sigmaC = RPREPAREMAT(bmOutputC.sigma);
        SET_ELEMENT(HFM, 5, sigmaC);
        UNPROTECT(2);
        
    }else{
        //Rprintf("baumWelchDiscrete: ");
        //Rprintf("start: ");

        outputD bmOutputDiscrete=baumWelchDiscrete(allObservationSequencesDiscrete,  startProbN, transProbSeqN,  transProbDiv,  emProbN, allNodeIndices, summationIndices, stateIndicesSingle, indicesXD,  allParents,  type);
        ///
        
        ///
         //Rprintf("transProbSeq: ");
        transProSeqbOut = RPREPAREMAT(bmOutputDiscrete.transProbSeq);
        SET_ELEMENT(HFM, 0, transProSeqbOut);
        ///
         //Rprintf("transProbDiv: ");
        transProbDivOut = RPREPAREMAT(bmOutputDiscrete.transProbDiv);
        SET_ELEMENT(HFM, 1, transProbDivOut);
        ///
        startProbOut = RPREPAREVECT(bmOutputDiscrete.prior.at(0));
        SET_ELEMENT(HFM, 2, startProbOut);
        ///
        llOut = RPREPAREVECT(bmOutputDiscrete.allLL);
        SET_ELEMENT(HFM, 3, llOut);
        //Rprintf("emProbN: ");
        //printVector(bmOutputDiscrete.emProb);
        SEXP emProbOut = RPREPAREMAT(bmOutputDiscrete.emProb);
        SET_ELEMENT(HFM, 4, emProbOut);
        UNPROTECT(2);
        
    }
    ///////////////////////////////////////////////////////////////////////
    
    
    return HFM;
}
}
