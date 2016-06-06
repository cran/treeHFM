HFMviterbi = function(observationSequence,nHStates,dataNodeIndices,priorInit,transmatInit,transInitDiv,type,emissionProb=c(),sigma=c(), mu=c() ){
    
   
    
    
    summationIndices=c()
    endState=nHStates
    startState=1
    for(ns in 1:nHStates){
        summationIndices=rbind(summationIndices,as.integer(startState:endState))
        startState=endState+1
        endState=endState+nHStates
    }
    #create state indices
    stateIndicesSingle=matrix(0,nHStates*nHStates,2)
    countS=1
    for(nS1 in 1:nHStates){
        for (nS2 in 1:nHStates){
            stateIndicesSingle[countS,]=as.integer(c(nS1,nS2));
            countS=countS+1
        }
    }
    ##
    nHStates=as.integer(nHStates)
    nDStates=as.integer(nHStates*nHStates)
    indicesXU=as.integer(1:(nHStates*nHStates))
    indicesXD=c()
    for (iX1 in 1: nHStates){
        iX1N=iX1
        for (iX2 in 1: nHStates){
            indicesXD=as.integer(c(indicesXD,iX1N))
            iX1N=iX1N+nHStates
        }
    }
    #
    parents=as.integer(dataNodeIndices[,1]);
    nodeIndices=as.integer(dataNodeIndices[,2]);
    #
    #

   
    type=as.character(type)
    ###
    if(length(emissionProb)>0){
        nEStates=as.integer(dim(emissionProb)[2])
    }else{
         nEStates=as.integer(0)
    }
    if(type=='d'){
        dimData=as.integer(1);
        ###
        nObs=as.integer(length(observationSequence));
        ###
        dimData=as.integer(dim(observationSequence));
        hfmV_out= .Call("RWrapperMS",priorInit,transmatInit,transInitDiv,emissionProb,observationSequence,nodeIndices,parents,type,nHStates,nDStates,nObs,dimData,c(),c(),summationIndices,indicesXD,nEStates);
    }else{
        ###
        nObs=as.integer(dim(observationSequence)[2]);
        ###
        dimData=as.integer(dim(observationSequence)[1]);
        ####
        hfmV_out= .Call("RWrapperMS",priorInit,transmatInit,transInitDiv,emissionProb,observationSequence,nodeIndices,parents,type,nHStates,nDStates,nObs,dimData,sigma,mu,summationIndices,indicesXD,nEStates);
    }
    hfmViterbi=list()
    hfmViterbi$maxPath = matrix(hfmV_out$maxPath, nrow=nObs, ncol=4, byrow=T)
    hfmViterbi$strainParents = matrix(hfmV_out$strainParents, nrow=hfmV_out$dims[[1]], ncol=hfmV_out$dims[[2]], byrow=T)
    hfmViterbi
}
    