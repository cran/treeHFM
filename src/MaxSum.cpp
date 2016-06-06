//
//  MaxSum.cpp
//  FGBMA
//
//  Created by Henrik Failmezger on 16.01.15.
//  Copyright (c) 2015 Henrik Failmezger. All rights reserved.
//

#include "MaxSum.h"
using namespace std;
void sendMessageVToFMSHMT( factorGraph &fg, int j, int i){
    ////cout<<"sendMessageVToFMSHMT"<<endl;
    ////cout<<"v: i: "<<i<<" ->"<<"f: j:" <<j<<endl;
    //std::vector<variableNode > v =fg.variableNodes;
    //std::vector<factorNode > f =fg.factorNodes;
    std::vector<double> msg_tmp = zeros(fg.variableNodes.at(i).size);
    //%msg_tmp = ones(1,3);
    //%msg_tmp = 1;
    double logS_tmp = 0.0;
    if (j>i){
        
        if (fg.variableNodes.at(i).observed != 1){
            std::vector<int> diffInd = setdiff(fg.variableNodes.at(i).f_neighbors,j);
            for (int m=0; m<diffInd.size();m++){
                int k=diffInd.at(m);
                //for (k = setdiff(v(i).f_neighbors,j)){
                if (k <= i) {//% down f(k)-->x(i)
                    if(fg.factorNodes.at(k).msg_out_down.size()>0){
                        msg_tmp = addEl(msg_tmp , fg.factorNodes.at(k).msg_out_down.at(findMessage(fg.factorNodes.at(k).msg_out_down,i)).msg);
                        //logS_tmp = logS_tmp + fg.factorNodes.at(k).acc_logS_msg_out_down.at(findMessage(fg.factorNodes.at(k).acc_logS_msg_out_up,i)).msg;
                    }
                }else {//% k > i up f(k)-->x(i)
                    if(fg.factorNodes.at(k).msg_out_up.size() >0){
                        msg_tmp = addEl(msg_tmp,fg.factorNodes.at(k).msg_out_up.at(findMessage(fg.factorNodes.at(k).msg_out_up,i)).msg);
                        //logS_tmp = logS_tmp + fg.factorNodes.at(k).acc_logS_msg_out_up.at(findMessage(fg.factorNodes.at(k).acc_logS_msg_out_up,i)).msg;
                    }
                }
            }
            //
            //msg_tmp_norm=msg_tmp*(1/scale);
            //
            //fg.factorNodes.at(j).msg_in_down_scale.push_back(scale);
            msg msg_in_down;
            msg_in_down.msg=msg_tmp;
            msg_in_down.index=i;
            //
            fg.factorNodes.at(j).msg_in_down.push_back(msg_in_down);
            fg.factorNodes.at(j).N_msg_in_down=fg.variableNodes.at(j).N_msg_in_down+1;
            //
            //cout<<"msg_in_down";
            //f(j).msg_in_down_norm=msg_tmp_norm;
        }else{

            msg_tmp = zeros(fg.variableNodes.at(i).size);
            //msg_tmp = ones(1,3);
            //
            //msg_tmp_norm=msg_tmp*(1/scale);
            //
            //fg.factorNodes.at(j).msg_in_down_scale.push_back(scale);
            msg msg_in_down;
            msg_in_down.msg=msg_tmp;
            msg_in_down.index=i;
            //
            fg.factorNodes.at(j).msg_in_down.push_back(msg_in_down);;
            fg.factorNodes.at(j).N_msg_in_down=fg.variableNodes.at(j).N_msg_in_down+1;
            msgAccLog mA;
            mA.index=i;
            //
            ////cout<<"msg_in_down";
            //f(j).msg_in_down_norm = msg_tmp_norm;
        }
        
    }else{//i>j
        if (fg.variableNodes.at(i).observed != 1) {
            //for (k = setdiff(v(i).f_neighbors,j)){
            std::vector<int> diffInd = setdiff(fg.variableNodes.at(i).f_neighbors,j);
            for (int m=0; m<diffInd.size();m++){
                int k=diffInd.at(m);
                if(k <= i){ // down f(k)-->x(i)
                    //%for m=1:length(f(k).msg_out_down)
                    if(fg.factorNodes.at(k).msg_out_down.size()>0){
                        msg_tmp = addEl(msg_tmp,fg.factorNodes.at(k).msg_out_down.at(findMessage(fg.factorNodes.at(k).msg_out_down,i)).msg);
                        double logValueMsg=0.0;
                        //if(fg.factorNodes.at(k).acc_logS_msg_out_down.size()>0){
                        //}
                        
                    }
                    //%end
	               }else{ // k > i up f(k)-->x(i)
                       //%for m=1:length(f(k).msg_out_up)
                       if(fg.factorNodes.at(k).msg_out_up.size()>0){
                           msg_tmp = addEl(msg_tmp,fg.factorNodes.at(k).msg_out_up.at((findMessage(fg.factorNodes.at(k).msg_out_up,i))).msg);
                       }
                       //%end
                   }
            }
            //%msg_tmp_norm=msg_tmp*(1/scale);
            //
            msg msg_in_up;
            msg_in_up.index=i;
            msg_in_up.msg=msg_tmp;
            //fg.factorNodes.at(j).msg_in_up_scale.at(i) = scale;
            //
            ////cout<<"msg_in_up";
            //
            fg.factorNodes.at(j).msg_in_up.push_back(msg_in_up);
            fg.factorNodes.at(j).N_msg_in_up=fg.variableNodes.at(j).N_msg_in_up+1;
        }else{
            msg_tmp = zeros(fg.variableNodes.at(i).size);
            //%msg_tmp = ones(1,3);
            //%
            double scale=sum(msg_tmp);
            //%msg_tmp_norm=msg_tmp*(1/scale);
            //%
            //fg.factorNodes.at(j).msg_in_up_scale.at(i) = scale;
            msg msg_in_up;
            msg_in_up.index=i;
            msg_in_up.msg=msg_tmp;
            //
            fg.factorNodes.at(j).msg_in_up.push_back(msg_in_up);
            ////cout<<"msg_in_up";
        }
    }
    
    fg.factorNodes.at(j).N_msg_in_sofar = fg.factorNodes.at(j).N_msg_in_sofar + 1;
    fg.variableNodes.at(i).N_msg_out_remain = fg.variableNodes.at(i).N_msg_out_remain - 1;
}

void sendMessageFToVMSHMT( factorGraph &fg, int j, int i,transitionMatrices tM){
     ////cout<<"sendMessageFToVMSHMT"<<endl;
     ////cout<<"f: j: "<<j<<" ->"<<"v: i:" <<i<<endl;
    //UNTITLED3 Summary of this function goes here
    //   Send a message from factor node f(j) to variable node v(i)
    
    //std::vector<variableNode > v =fg.variableNodes;
    //std::vector<factorNode > f =fg.factorNodes;
    std::vector<double > msg_tmp = zeros(fg.variableNodes.at(i).size);
    double logS_tmp = 0.0;
    int strain=0;
    std::vector<std::vector<double > >msg_tmpM;
     std::vector<double > msg_tmpMS;
    if(j>i){
        
        if(fg.factorNodes.at(j).leaf != 1){
            std::vector<int> diffInd = setdiff(fg.factorNodes.at(j).v_neighbors,i);
            for (int m=0; m<diffInd.size();m++){
                int k=diffInd.at(m);
                if(k < j){//  down x(k)-->f(j)
                    msg_tmp =addEl(msg_tmp,fg.factorNodes.at(j).msg_in_down.at(findMessage(fg.factorNodes.at(j).msg_in_down,k)).msg);
                }else{// % k > i up x(k)-->f(j)
                    strain=k;
                    msg_tmp =addEl(msg_tmp,fg.factorNodes.at(j).msg_in_up.at(findMessage(fg.factorNodes.at(j).msg_in_up,k)).msg);
                    
                }
            }
            ////cout<<"---------------------------------------------msg_tmp start-------------------------";
            
            std::vector<double > xM=msg_tmp;
            std::vector<double >xMN=msg_tmp;
           
            //replaces repmat because of speed
            //replaces repmat because of speed
            
            if(fg.factorNodes.at(j).division){
                ////cout<<"division"<<endl;
                strain=0;
                std::vector< std::vector<double > > msg_tmp_strainUp;
                std::vector<int> diffInd = setdiff(fg.factorNodes.at(j).v_neighbors,i);
                for (int m=0; m<diffInd.size();m++){
                    int k=diffInd.at(m);
                    msg_tmp_strainUp.push_back(fg.factorNodes.at(j).msg_in_up.at(findMessage(fg.factorNodes.at(j).msg_in_up,k)).msg) ;
                }
                int variableIndex=0; // i do not need it as everything is summed
                //if(strain>i)
                // variableIndex=2;
                //end
                //UP
                std::vector<int>v_neighbors=fg.factorNodes.at(j).v_neighbors;
                int minIndex=min(fg.factorNodes.at(j).v_neighbors).index;
                v_neighbors.erase(v_neighbors.begin() + minIndex);
                int minN=min(v_neighbors).value;
                std::vector<double >xEU=fg.factorNodes.at(j).probO.at(minN); //'index on + strain is smaller
                
                int maxN=max(v_neighbors).value;
                std::vector<double >xED=fg.factorNodes.at(j).probO.at(maxN); //'index on - strain is larger
                //+ strain
                std::vector<double >oMU;
                for(int c=0;c<fg.variableNodes.at(i).size;c++){
                    for (int l=0;l<fg.variableNodes.at(i).size;l++){
                        oMU.push_back(xEU.at(c));
                    }
                }
                //- strain
                std::vector<double >oMD;
                for(int c=0;c<fg.variableNodes.at(i).size;c++){
                    for (int l=0;l<fg.variableNodes.at(i).size;l++){
                        oMD.push_back(xED.at(c));
                    }
                }
                oMD=arrangeItems(oMD,fg.factorNodes.at(j).indicesXD); //rearrange the entries for X2
                
                //UP
                std::vector<double> xMU;
                for(int c=0;c<fg.variableNodes.at(i).size;c++){
                    for (int l=0;l<fg.variableNodes.at(i).size;l++){
                        xMU.push_back(msg_tmp_strainUp.at(0).at(c));
                    }
                }
                
                //DOWN
                std::vector<double> xMD;
                for (int c=0;c<fg.variableNodes.at(i).size;c++){
                    for (int l=0;l<fg.variableNodes.at(i).size;l++){
                        xMD.push_back(msg_tmp_strainUp.at(1).at(c));
                    }
                }
                xMD=arrangeItems(xMD,fg.factorNodes.at(j).indicesXD); //rearrange the entries for X2
                //
                std::vector<std::vector<double > > xMUrep=repmat(xMU,tM.probDiv.at(variableIndex).size(),1);
                std::vector<std::vector<double > > xMDrep=repmat(xMD,tM.probDiv.at(variableIndex).size(),1);
                std::vector<std::vector<double > > oMUrep= repmat(oMU,tM.probDiv.at(variableIndex).size(),1);
                std::vector<std::vector<double > > oMDrep= repmat(oMD,tM.probDiv.at(variableIndex).size(),1);
                //
                msg_tmpM=addEl(xMUrep,log2M(oMUrep));

                
                msg_tmpM=addEl(msg_tmpM,log2M(tM.probDiv.at(variableIndex)));
                msg_tmpM=addEl(msg_tmpM,log2M(oMDrep));
                msg_tmpM=addEl(msg_tmpM,xMDrep);
                if(j==111){
                    ////cout<<"msg_tmpM"<<endl;
                }
                
            }else{
               
                ////cout<<"straight"<<endl;
                std::vector<double >xE=fg.factorNodes.at(j).probO.at(max(fg.factorNodes.at(j).v_neighbors).value);//';
                std::vector<std::vector<double > >xErep=repmat(xE,tM.probSeq.at(0).size(),1);
                
                std::vector<std::vector<double > >xMrep=repmat(xM,tM.probSeq.at(0).size(),1);
                msg_tmpM=addEl(xMrep,log2M(xErep));
                msg_tmpM =addEl(msg_tmpM,log2M(tM.probSeq.at(0)));
                ////cout<<"msg_tmpM"<<endl;
                if(j==139){
                  //  //cout<<"---------------------------------------------end node-------------------------";
                  //  //cout<<"rep1"<<endl;
                  //  //cout<<"msg_tmpM"<<endl;
                }
            }
            
            
        }
        // transpose message
        //
        //
        //
        if(j==111){
            ////cout<<"msg_tmpM: "<<endl;
        }
        std::vector<double > msg_tmp;
        if(msg_tmpM.size()>1){
            MaxMatrix mMC=maxMatrix(msg_tmpM,1);//the sum is implicit
            std::vector<int > cordsIndicesR=createSequence(msg_tmpM.size(),0,msg_tmpM.size());
            if(j==111){
                ////cout<<"MSG TEMP: "<<endl;
                ////cout<<"cords indices"<<endl;
            }
            std::vector<int > cordsIndicesC=mMC.index;
            if(j==111){
                ////cout<<"cords values"<<endl;
            }
            std::vector<double > msg_tmp_summed1;
            for(int nc=0;nc<msg_tmpM.size();nc++){
                msg_tmp_summed1.push_back(msg_tmpM.at(cordsIndicesR.at(nc)).at(cordsIndicesC.at(nc)));
            }
            if(j==111){
                ////cout<<"msg_tmp_summed1"<<endl;
            }
            msg_tmp=msg_tmp_summed1;
            //
            msgB msgBUP;
            msgBUP.bd=cordsIndicesC;
            
            Max maxIndex=max(setdiff(fg.factorNodes.at(j).v_neighbors,i));
            msgBUP.index=maxIndex.value;
            
            max(setdiff(fg.factorNodes.at(j).v_neighbors,i));
            //
            fg.factorNodes.at(j).backwardUp.push_back(msgBUP); // backward down values
        }else{
            msg_tmp=msg_tmpMS;
            //
            msgB msgBUP;
            msgBUP.bd=createSequence(msg_tmp.size(),0,msg_tmp.size());
            msgBUP.index=i;
            //
            fg.factorNodes.at(j).backwardUp.push_back(msgBUP);
        }
        
        msg msg_out_up;
        msg_out_up.msg=msg_tmp;
        msg_out_up.index=i;
        fg.factorNodes.at(j).msg_out_up.push_back(msg_out_up);
        fg.factorNodes.at(j).N_msg_out_remain = 0;
        ////cout<<"-------------------msg_out_up-----------------"<<endl;
        //fg.factorNodes.at(j).msg_out_up_scale.at(i)=scale;
    }else{
        //j;
        int obs=0;
        std::vector<double >msg_tmpM_summed;
        std::vector<double > msg_out;
        std::vector<double > msg_tmp_strainUp=zeros(fg.variableNodes.at(i).size);
        std::vector<double > msg_tmp_strainDown=zeros(fg.variableNodes.at(i).size);
       // //cout<<"upwards"<<endl;
       // //cout<<fg.factorNodes.at(j).leaf<<endl;
        if (fg.factorNodes.at(j).leaf != 1){
        ////////////////////////////////////
            if(j==20){
                //cout<<"v_neighbors"<<endl;
            }
            std::vector<int> diffInd = setdiff(fg.factorNodes.at(j).v_neighbors,i);
            //cout<<"diffInd"<<endl;
            double logS_tmp = 0.0;
            for (int m=0; m<diffInd.size();m++){
                int k=diffInd.at(m);
                if(k < j){ //% down x(k)-->f(j)
                    if(j==20){
                        //cout<<"k: "<<k<<"j: "<<j<<endl;
                    }
                    msg_tmp = addEl(msg_tmp,fg.factorNodes.at(j).msg_in_down.at(findMessage(fg.factorNodes.at(j).msg_in_down,k)).msg);
                }else{ //% k > i up x(k)-->f(j)
                    if(j==20){
                        //cout<<"k: "<<k<<"j: "<<j<<endl;
                   
                    }
                    msg_tmp = addEl(msg_tmp,fg.factorNodes.at(j).msg_in_up.at(findMessage(fg.factorNodes.at(j).msg_in_up,k)).msg);
                }
            }
            if(j==20){
                //cout<<"msg_tmp"<<endl;
            }
            if(fg.factorNodes.at(j).division){ //check if division =msg_tmp*(1/scale) true and summ states
                
                //cout<<"division msg"<<endl;
                int variableIndex=0;
                if(strain>i){
                    variableIndex=1;
                }
                std::vector<int>v_neighbors=fg.factorNodes.at(j).v_neighbors;
                Min minNN = min(fg.factorNodes.at(j).v_neighbors);
                int minIndex=minNN.index;
                v_neighbors.erase(v_neighbors.begin() + minIndex);
                Min minV =min(v_neighbors);
                int minN=minV.value;
                ////////////////////////////////////////////////////
                std::vector<double >xEU=fg.factorNodes.at(j).probO.at(minN); //'index on + strain is smaller
                Max maxV=max(v_neighbors);
                int maxN=maxV.value;
                std::vector<double >xED=fg.factorNodes.at(j).probO.at(maxN); //'index on - strain is larger
                //+ strain
                std::vector<double >oMU;
                for(int c=0;c<fg.variableNodes.at(i).size;c++){
                    for (int l=0;l<fg.variableNodes.at(i).size;l++){
                        oMU.push_back(xEU.at(c));
                    }
                }
                if(j==20){
                    //cout<<"oMU"<<endl;
                }
                //- strain
                std::vector<double > oMD;
                for(int c=0;c<fg.variableNodes.at(i).size;c++){
                    for (int l=0;l<fg.variableNodes.at(i).size;l++){
                        oMD.push_back(xED.at(c));
                    }
                }
                oMD=arrangeItems(oMD,fg.factorNodes.at(j).indicesXD); //rearrange the entries for X2
                if(j==20){
                    //cout<<"oMD"<<endl;
                }
                ////////////////////////////////////////////////////
                std::vector<double >msg_tmp_strainUp = ones(fg.variableNodes.at(i).size);
                std::vector<double >msg_tmp_strainDown= ones(fg.variableNodes.at(i).size);
                int strain=0;
                std::vector<int> diffInd = setdiff(fg.factorNodes.at(j).v_neighbors,i);
                for (int m=0; m<diffInd.size();m++){
                    int k=diffInd.at(m);
                    if(k < j){ // down x(k)-->f(j)
                        msg_tmp_strainDown = multEl(msg_tmp_strainDown,fg.factorNodes.at(j).msg_in_down.at(findMessage(fg.factorNodes.at(j).msg_in_down,k)).msg);
                    }else{ // k > i up x(k)-->f(j)
                        strain=k;
                        msg_tmp_strainUp = multEl(msg_tmp_strainUp,fg.factorNodes.at(j).msg_in_up.at(findMessage(fg.factorNodes.at(j).msg_in_up,k)).msg);
                    }
                }
                variableIndex=0;
                if(strain<i){
                    variableIndex=1;
                }
                
                ////cout<<"msg_tmp_strainDown"<<endl;
                ////cout<<"msg_tmp_strainUp"<<endl;
                ////////////////////////////////////////////////////
                std::vector<std::vector<double > > xMDown=repmat(msg_tmp_strainDown,tM.probDiv.at(variableIndex).at(0).size(),1);
                std::vector<std::vector<double > > xMDownT=transpose(xMDown);
                //UP
                ////cout<<"UP"<<endl;
                std::vector<double > xMUp;
                for(int c=0;c<tM.probDiv.at(variableIndex).size();c++){
                    for (int l=0;l<msg_tmp_strainUp.size();l++){
                        xMUp.push_back(msg_tmp_strainUp.at(l));
                    }
                }
                ////cout<<"xMUp"<<endl;
                std::vector<std::vector<double > > xMUprep=repmat(xMUp,tM.probDiv.at(variableIndex).size(),1);
                //
                std::vector<std::vector<double > > oMUrep=repmat(oMU,tM.probDiv.at(0).size(),1);
                ////cout<<"oMUrep"<<endl;
                std::vector<std::vector<double > > oMDrep=repmat(oMD,tM.probDiv.at(0).size(),1);
                std::vector<std::vector<double > > transmatObs=addEl(log2M(oMUrep),log2M(tM.probDiv.at(0)));
                ////cout<<"oMDrep"<<endl;
                ////cout<<"transmatObs"<<endl;
                //

                transmatObs=addEl(transmatObs,log2M(oMDrep));
                ////cout<<endl<<variableIndex<<endl;
                if(variableIndex==1){
                    transmatObs=arrangeItems(transmatObs,fg.factorNodes.at(j).indicesXD);
                }
                ////cout<<"transmatObs"<<endl;
                ////cout<<"msg_tmpM"<<endl;
                msg_tmpM=addEl(xMDownT,transmatObs);
               // //cout<<"xMUprep"<<endl;
                msg_tmpM=addEl(msg_tmpM,xMUprep);
                ////cout<<"msg_tmpM"<<endl;
                //
                std::vector<std::vector<int> > indicesToSum=fg.factorNodes.at(j).summationIndices;
                //
                //

                //////////////////////////////////////////////////////////////////
        }else{
            msg_tmp; //'
            
            std::vector<double > xM=msg_tmp;
            std::vector<std::vector<double > > xMrep=repmat(xM,tM.probDiv.at(0).size(),1);
            std::vector<std::vector<double > > xMrepT=transpose(xMrep);
            std::vector<double> xE=fg.factorNodes.at(j).probO.at(i); //' format [p(o|r) p(o|g) ...]
            std::vector<std::vector<double > >xErep=repmat(xE, tM.probDiv.at(0).size(),1);
            
            obs=0;
            msg_tmpM=addEl(xMrepT,log2M(xErep));
            msg_tmpM=addEl(msg_tmpM,log2M(tM.probSeq.at(0)));
            
            /*
            if(j==2){
                //cout<<"xMrepT: "<<endl;
                //cout<<"xErep: "<<endl;
                //cout<<"prob seq"<<endl;
                //cout<<"msg_tmpM"<<endl;
            }
             */
        }
        }else{
                std::vector<double> xE=fg.factorNodes.at(j).probO.at(i);
            
                msg_tmp = addEl(fg.factorNodes.at(j).msg_in_down.at(findMessage(fg.factorNodes.at(j).msg_in_down,1)).msg,log2M(tM.probStart.at(0).at(0)));
           /*
            if(j==1){
                //cout<<"------------------------------------------------"<<endl;
                //cout<<"message"<<endl;
                //cout<<"obs"<<endl;
                //cout<<"xE"<<endl;
                //cout<<endl;
                //cout<<"msg_tmp1"<<endl;
            
            }
            */
                msg_tmpM_summed = addEl(msg_tmp,log2M(xE));
            if(j==1){
             //   //cout<<"msg_tmpM_summed"<<endl;
             //
            }
        }
        std::vector<double> msg_tmp_summed1;
        
        //Summation from downwards variable
        ////cout<<"msg_tmpM"<<endl;
        std::vector<int> maxVj;
        if(msg_tmpM.size()>1){
            if(j==20){
                //cout<<"msg_tmpM"<<endl;
            }
            
            MaxMatrix mMC=maxMatrix(msg_tmpM,2);//the sum is implicit
            ////cout<<"sequence"<<endl;
            std::vector<int > cordsIndicesR=createSequence(msg_tmpM.at(0).size(),0,msg_tmpM.at(0).size());
            if(j==20){
                //cout<<"cordsIndicesR"<<endl;
                //cout<<endl;
            }
            std::vector<int > cordsIndicesC=mMC.index;
            if(j==20){
                //cout<<"cordsIndicesC"<<endl;
            }
          //  //cout<<msg_tmpM.at(0).size()<<endl;
            for(int nc=0;nc<msg_tmpM.at(0).size();nc++){
                if(j==20){
                 //   //cout<<"nc: "<<nc<<endl;
                 //   //cout<<cordsIndicesR.at(nc)<<endl;
                 //   //cout<<cordsIndicesC.at(nc)<<endl;
                 //   //cout<<" m "<<msg_tmpM.at(cordsIndicesC.at(nc)).at(cordsIndicesR.at(nc))<<endl;
                }
                msg_tmp_summed1.push_back(msg_tmpM.at(cordsIndicesC.at(nc)).at(cordsIndicesR.at(nc)));
            }
            std::vector<double>  msg_tmp=msg_tmp_summed1;
            maxVj=cordsIndicesC;
            //
            ////cout<<"assign"<<endl;
            msg_out=msg_tmp_summed1;
            msgB msgBD;
            msgBD.bd=maxVj;
            msgBD.index=max(setdiff(fg.factorNodes.at(j).v_neighbors,i)).value;
            //
            fg.factorNodes.at(j).backwardDown.push_back(msgBD); // backward down values
        }else{
          //  //cout<<"msg_tmpM single:"<<endl;
            msg_tmp=msg_tmpM_summed;
            //
            maxVj=zerosINT(msg_tmpM_summed.size());
            msgB msgBD;
            msgBD.bd=maxVj;
            msgBD.index=i;
            //
            fg.factorNodes.at(j).backwardDown.push_back(msgBD);////?????????
            msg_out=msg_tmp;
        }
        if(j==20){
            //cout<<"msg_tmp_summed1"<<endl;
        }
        //Summation from upwards variable 
        if(fg.factorNodes.at(j).division){ //check if division = true and summ states
           // //cout<<"division: "<<endl;
            std::vector< std::vector<int> > indicesToSum=fg.factorNodes.at(j).summationIndices;
            std::vector<double>  msg_tmp_summed2;
            std::vector<int> cordsMaxValues;
            std::vector<int> maxVjN;
            for(int n=0;n<indicesToSum.size();n++){
                ////cout<<"n :"<<n<<endl;
                Max maxValueA=max( getValuesIndex(msg_tmp_summed1,indicesToSum.at(n)));
                //cout<<"m2"<<endl;
                //cout<<"v: "<<maxValueA.value<<endl;
                msg_tmp_summed2.push_back(maxValueA.value);
                ///
                cordsMaxValues.push_back(indicesToSum.at(n).at(maxValueA.index));
                ///
                maxVjN.push_back(maxValueA.index);
            }
            if(j==20){
                //cout<<"msg_tmp_summed2"<<endl;
            }
            ////cout<<"indices to sum"<<endl;
            //
            Max maxValueBD=max(setdiff(fg.factorNodes.at(j).v_neighbors,i));
            Min minValueBD=min(setdiff(fg.factorNodes.at(j).v_neighbors,i));
            //
            //cout<<"msgBD"<<endl;
            msgB msgBD;
            std::vector<int> arrangedM=arrangeItems(maxVj,cordsMaxValues);
            msgBD.bd=arrangedM;
            msgBD.index=minValueBD.value;
            //
            //
            //cout<<"msgBDMI"<<endl;
            msgB msgBDMI;
            std::vector<int> arrangedMI=arrangeItems(maxVj,cordsMaxValues);
            msgBDMI.bd=arrangedMI;
            msgBDMI.index=maxValueBD.value;
            //
            fg.factorNodes.at(j).backwardDown.push_back(msgBD);
            fg.factorNodes.at(j).backwardDown.push_back(msgBDMI);
            //
            
            msg_out =msg_tmp_summed2;
        }
        msg msg_out_down;
        msg_out_down.index=i;
        msg_out_down.msg=msg_out;
        //cout<<"----------- MSGOUTDOWN --------"<<endl;
        
        fg.factorNodes.at(j).msg_out_down.push_back(msg_out_down);
}
fg.variableNodes.at(i).N_msg_in_sofar = fg.variableNodes.at(i).N_msg_in_sofar + 1;
fg.factorNodes.at(j).N_msg_out_remain = fg.factorNodes.at(j).N_msg_out_remain - 1;
}
/////////
std::vector<std::vector<int> > maxSumAlgorithm(factorGraph fgS, std::vector<int> endNodes, std::vector<int> startNodes, transitionMatrices tM, int N_V){
    factorGraph fg=fgS;
    
    int done = 0;
    //% logS = 0; % log of the scale factor S, where S is the normalization factor of each message
    int counter=0;
    int counterM=0;
    //%%%message passing upwards
    std::string msg;
    //cout<<"MS endNodes"<<endl;
    
    while (endNodes.size()>0){
        int actNode=endNodes.at(0);
        
        //%send messages from the end node until it is not possible any
        //%more
        int stopSending=0;
        if(fg.variableNodes.at(actNode).f_children.size()==0){
            msg="v-->f"; //%ending variable node of a strain
        }else{
            msg="f-->v"; //%division factor node
            if (fg.factorNodes.at(actNode).N_msg_in_sofar >= fg.factorNodes.at(actNode).N_v_neighbors - 1){ //%the node was already visited
                stopSending=1;
                endNodes.erase(endNodes.begin());
            }
        }
        
        while(stopSending==0){
            
            if(msg=="f-->v"){
                std::vector<int> v_parent=fg.factorNodes.at(actNode).v_parents;
                
                if(v_parent.size()==0){
                    if(intersect(endNodes,actNode)==0){
                        endNodes.push_back(actNode);//%%add the actual node to the stack
                    }
                    endNodes.erase(endNodes.begin());//%drop the endNode from the stack
                    stopSending=1;
                    break;
                }else{
                    if (fg.factorNodes.at(actNode).N_msg_in_sofar >= fg.factorNodes.at(actNode).N_v_neighbors - 1){ //if there have not enough messages arrived, stop go to the next end node
                        actNode;
                        v_parent;
                        sendMessageFToVMSHMT(fg,actNode,v_parent.at(0),tM); //only one v_parent possible
                    }else{
                        if(intersect(endNodes,actNode)==0){
                            endNodes.push_back(actNode);//%%add the actual node to the stack
                        }
                        endNodes.erase(endNodes.begin());//%drop the endNode from the stack
                        stopSending=1;
                        break;
                    }
                    actNode=v_parent.at(0); //only one v_parent possible
                    msg="v-->f";
                }
            }else if(msg=="v-->f"){
                if(fg.variableNodes.at(actNode).f_neighbors.size()==0){
                    if(intersect(endNodes,actNode)==0 ){
                        endNodes.push_back(actNode);//%%add the actual node to the stack
                    }
                    endNodes.erase(endNodes.begin());//%drop the endNode from the stack
                    stopSending=1;
                    break;
                }else{
                    int f_parent=min(fg.variableNodes.at(actNode).f_neighbors).value;
                    
                    sendMessageVToFMSHMT(fg,f_parent,actNode);
                    actNode=f_parent;
                    msg="f-->v";
                }
            }
        }
    }
    //%%%message passing downwards
    while(startNodes.size()>0){
        int actNode=startNodes.at(0);
        //%send messages from the end node until it is not possible any
        //%more
        int stopSending=0;
        msg="f-->v"; //%start node is always a factor node
        while (stopSending==0){
            
            if(msg=="f-->v"){
                std::vector<int> v_child=fg.factorNodes.at(actNode).v_children;
                if(v_child.size()>0){ //%%if there have not enough messages arrived, stop go to the next end node
                    int v_childS;
                    if(v_child.size()>1){//%%division
                        if(intersect(startNodes,actNode)==0){
                            startNodes.push_back(actNode); //add the other strain to the startNodes
                            v_childS=max(v_child).value;
                        }else{ //%+ strain was already visited, start with - strain
                            v_childS=min(v_child).value;
                        }
                    }else{
                        v_childS=v_child.at(0);
                    }
                    sendMessageFToVMSHMT(fg,actNode,v_childS,tM);
                    actNode=v_childS;
                    msg="v-->f";
                }else{
                    startNodes.erase(startNodes.begin());//%drop the startNode from the stack
                    stopSending=1;
                    break;
                }
                
            }else if(msg=="v-->f"){
                
                'g';
                if(fg.variableNodes.at(actNode).f_children.size()>0){
                    int f_child=max(fg.variableNodes.at(actNode).f_children).value;
                    sendMessageVToFMSHMT(fg,f_child,actNode);
                    actNode=f_child;
                    msg="f-->v";
                }else{
                    startNodes.erase(startNodes.begin());//drop the startNode from the stack
                    stopSending=1;
                    break;
                }
            }
            
        }
    }
    
    
    //////////////////////////////////////////MAXPATH///////////////////////////////////////
    std::vector<int> visitedNodes=zerosINT(N_V+1);
    std::vector<std::vector<int> > maxPath;
    //cout<<"N_V: "<<N_V<<endl;
    
    //cout<<"calc max path"<<endl;
    //cout<<"end nodes"<<endl;
    int pathWalk=0;
    int maxStateNext;
    //
    std::vector<int> f_parent;
    std::vector<int> v_parent;
    //create queue
    std::queue<int> strains;
    for (int s=0; s<fg.endNodes.size();s++){
        strains.push(fg.endNodes.at(s));
    }
    
    
    int counterMax=0;
    while(strains.size()>0){
        int actNode=strains.front();
        //cout<<"-----------------------------------------------------------------"<<endl;
        //cout<<"actNode: "<<actNode<<endl;
        visitedNodes.at(actNode)=1;
       
        f_parent=fg.variableNodes.at(actNode).f_parents;
        v_parent=fg.variableNodes.at(actNode).v_parents;
        //cout<<"v_parent"<<v_parent.at(0)<<endl;
        //cout<<"msg: "<<findMessage(fg.factorNodes.at(f_parent.at(0)).msg_out_down,actNode)<<endl;
        
        Max maxValue;
        maxValue=max(fg.factorNodes.at(f_parent.at(0)).msg_out_down.at(findMessage(fg.factorNodes.at(f_parent.at(0)).msg_out_down,actNode)).msg);
        counterMax++;
        //
        //cout<<"maxValue: "<<maxValue.value<<endl;
        double pMax=maxValue.value;
        int maxState=maxValue.index;
        //cout<<"maxState: "<<maxState<<endl;
        //
        std::vector<int> maxPathValue;
        maxPathValue.push_back(actNode);
        maxPathValue.push_back(maxState+1);
        maxPathValue.push_back(fg.variableNodes.at(actNode).strain);
        //cout<<"v_parent: "<<v_parent.at(0)<<endl;
        maxPathValue.push_back(fg.variableNodes.at(v_parent.at(0)).strain);
        //
        //cout<<"findMessage(fg.factorNodes.at(f_parent.at(0)).backwardDown,v_parent.at(0)): "<<findMessage(fg.factorNodes.at(f_parent.at(0)).backwardDown,v_parent.at(0))<<endl;
        maxPath.push_back(maxPathValue);
        maxStateNext=fg.factorNodes.at(f_parent.at(0)).backwardDown.at(findMessage(fg.factorNodes.at(f_parent.at(0)).backwardDown,v_parent.at(0))).bd.at(maxState);
        actNode=fg.variableNodes.at(actNode).v_parents.at(0);
        //cout<<"path walk"<<endl;
        //cout<<"actNode: "<<actNode<<endl;
        
        pathWalk=1;
        while(pathWalk==1){
            if(visitedNodes.at(actNode)==0){
                f_parent=fg.variableNodes.at(actNode).f_parents;
                v_parent=fg.variableNodes.at(actNode).v_parents;
                
                maxState=maxStateNext;
                //
                //cout<<"maxPathValue"<<endl;
                std::vector<int> maxPathValue;
                maxPathValue.push_back(actNode);
                maxPathValue.push_back(maxState+1);
                //cout<<fg.variableNodes.at(actNode).strain<<endl;
                maxPathValue.push_back(fg.variableNodes.at(actNode).strain);
                //cout<<"v_parents"<<endl;
                if(fg.variableNodes.at(actNode).v_parents.size()==0){
                    //cout<<"v_parent: 0"<<endl;
                    //cout<<"next: "<<endl;
                    maxPathValue.push_back(fg.variableNodes.at(actNode).strain);
                }else{
                    maxPathValue.push_back(fg.variableNodes.at(v_parent.at(0)).strain);
                }
                maxPath.push_back(maxPathValue);
                //
                visitedNodes.at(actNode)=1;
                
                
                //
                if(fg.variableNodes.at(actNode).v_parents.size()==0){
                    strains.pop();
                    pathWalk=0;
                    break;
                }else{
                    //cout<<"maxStateNext"<<endl;
                    //cout<<"actNode) "<<actNode<<endl;
                    //cout<<"f_parent.at(0) "<<f_parent.at(0)<<endl;
                    //cout<<"fg.factorNodes.at(f_parent.at(0)).division: "<<fg.factorNodes.at(f_parent.at(0)).division<<endl;
                    //cout<<"v_parent.at(0) "<<v_parent.at(0)<<endl;
                    actNode=fg.variableNodes.at(actNode).v_parents.at(0);
                    
                    //cout<<"findMessage(fg.factorNodes.at(f_parent.at(0)).backwardDown,v_parent.at(0))"<<findMessage(fg.factorNodes.at(f_parent.at(0)).backwardDown,v_parent.at(0))<<endl;
                    maxStateNext=fg.factorNodes.at(f_parent.at(0)).backwardDown.at(findMessage(fg.factorNodes.at(f_parent.at(0)).backwardDown,v_parent.at(0))).bd.at(maxState);
                }
            }else{
                strains.pop();
               // strains.erase(strains.begin()+1);
                pathWalk=0;
                break;
            }
            
            }
            pathWalk=0;
    }

    return maxPath;
}
