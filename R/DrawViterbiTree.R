DrawViterbiTree = function(maxPath){
    cellsStrain=matrix(0,dim(maxPath)[1],2);
    cellsStrain[maxPath[,1],1]=maxPath[,2];
    #order the entries in the right direction
    cellsStrain[maxPath[,1],2]=maxPath[,3];
    ######
    parents=c();
    for (st in maxPath[,3]){
        parentStrain=maxPath[maxPath[,3]==st,4];
        stP=unique(parentStrain);
        indexFP=which(stP !=st);
        if(length(indexFP)>0){
            parents=rbind(parents,c(stP[indexFP],st))
        }
    }
    IA = duplicated(parents[,2])
    parents=parents[-which(IA),]
    nodeSubtree=list()
    ####
    if(length(parents)==0){
        submatrix=cellsStrain[,1];
    }else{
        nodes=unique(as.vector(parents));
        visitedNodes=array(0,length(nodes));
        
        startNodes=c(1);
        
        matrix=c();
        stopDrawing=0;
        
        while(!stopDrawing){
            for(subtree in nodes){
                children=parents[parents[,1]==subtree,2];
                if (length(children)==0){
                    matrix=cellsStrain[cellsStrain[,2]==subtree,1];
                    nodeSubtree[[subtree+1]]=matrix;
                    visitedNodes[subtree+1]=1;
                }else{
                    child1=children[1]
                    child2=children[2]
                    if(visitedNodes[child1+1]==1 & visitedNodes[child2+1]==1){
                        if(is.null(dim(nodeSubtree[[child1+1]]))){
                            rowsChild1=1
                            colsChild1=length(nodeSubtree[[child1+1]])
                        }else{
                            rowsChild1=dim(nodeSubtree[[child1+1]])[1]
                            colsChild1=dim(nodeSubtree[[child1+1]])[2]
                        }
                        if(is.null(dim(nodeSubtree[[child2+1]]))){
                            rowsChild2=1
                            colsChild2=length(nodeSubtree[[child2+1]])
                        }else{
                            rowsChild2=dim(nodeSubtree[[child2+1]])[1]
                            colsChild2=dim(nodeSubtree[[child2+1]])[2]
                        }
                        
                       
                        matrixNew=cellsStrain[cellsStrain[,2]==subtree,1]
                        if(is.null(dim(matrixNew))){
                            nrowsM1=1
                            ncolsM1=length(matrixNew)
                        }else{
                            nrowsM1=dim(matrixNew)[1]
                            ncolsM1=dim(matrixNew)[2]
                        }
                        matrixNew=t(as.matrix(matrixNew))
                        matrixJoin=matrix(0,rowsChild1+rowsChild2+nrowsM1,max(colsChild1,colsChild2))
                        matrix2=matrix(0,rowsChild1+rowsChild2+nrowsM1,ncolsM1)
                        matrixJoin[1:rowsChild1,1:colsChild1]=nodeSubtree[[child1+1]];
                        matrixJoin[(rowsChild1+dim(matrixNew)[1]+1):(rowsChild1+1+dim(matrixNew)[1]+(rowsChild2-1)),1:colsChild2]=nodeSubtree[[child2+1]];
                        #add the new strain
                        matrix2[(rowsChild1+1):(rowsChild1+1+dim(matrixNew)[1]-1),1:dim(matrixNew)[2]]=matrixNew;
                        #
                        nodeSubtree[[subtree+1]]=cbind(matrix2,matrixJoin)
                        visitedNodes[subtree+1]=1;
                    }
                }
                if(length(which(visitedNodes==0))==0){
                    stopDrawing=1;
                    break
                }
            }
        }
    }
    submatrix=nodeSubtree[[1]];
    submatrix
}