HFMfit = function(data,nHStates,dataNodeIndices,priorInit,transmatInit,transInitDiv,type,emissionProb=c(),SigmaInit=c(), muInit=c() ){
    
   
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
    vSize=6;
    parents=list()
    nodeIndices=list()
    for (eN  in 1:length(dataNodeIndices)){
        parents[[eN]]=as.integer(dataNodeIndices[[eN]][,1])
        nodeIndices[[eN]]=as.integer(dataNodeIndices[[eN]][,2])
    }
    #
    ###
    nObs=as.integer(length(data));
    ###
    dimData=as.integer(dim(data[[1]])[1]);
    type=as.character(type)
    ###
    nHStates=as.integer(nHStates)
    if(length(emissionProb)>0){
        nEStates=as.integer(dim(emissionProb)[2])
    }else{
         nEStates=as.integer(0)
    }
    if(type=='d'){
        dimData=as.integer(1);
        hfm_out= .Call("RWrapperBM",priorInit,transmatInit,transInitDiv,emissionProb, data,nodeIndices,parents,summationIndices,indicesXD,type,stateIndicesSingle,nHStates,nDStates,nObs,dimData,nEStates,c(),c());
    }else{
        hfm_out= .Call("RWrapperBM",priorInit,transmatInit,transInitDiv,emissionProb, data,nodeIndices,parents,summationIndices,indicesXD,type,stateIndicesSingle,nHStates,nDStates,nObs,dimData,nEStates,SigmaInit,muInit);
    }
    hfm=list()
    hfm$transMatSeq = matrix(hfm_out$transMatSeq, nrow=nHStates, ncol=nHStates, byrow=T)
    hfm$transMatDiv = matrix(hfm_out$transMatDiv, nrow=nHStates, ncol=nDStates, byrow=T)
    hfm$ll = hfm_out$loglik
    if(type=='d'){
        #map emissions
        maxEState=nEStates
        emissionM=matrix(hfm_out$emission, nrow=nHStates, ncol=nEStates+1, byrow=T)
        emissionM=emissionM[,1:maxEState]
        #
        hfm$emission = emissionM
    }else{
        hfm$mu = matrix(hfm_out$mu, nrow=nHStates, ncol=nHStates, byrow=T)
        hfm$sigma = array(hfm_out$sigma, c(dimData,dimData,nHStates))
    }
    hfm$initProb = hfm_out$initProb
    hfm
}
    