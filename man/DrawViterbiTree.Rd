\name{DrawViterbiTree}
\alias{DrawViterbiTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Arranges a Viterbi tree}
\description{
The functions arranges the output of HFMviterbi in a way that can be displayed.
}
\usage{
	DrawViterbiTree(maxPath)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{maxPath}{The output of HFMViterbi.}

}
\details{
The function arranges the branches of the tree in a matrix.
}
\value{
A matrix, foreground: Viterbi tree, background:0
}
\author{ Henrik Failmezger, failmezger@googlemail.com}


\examples{
## Fit a continous treeHFM ## 
nHStates = 2;
########create observation sequences########
obs1 <- rbind(runif(10,0,1),runif(10,0,1))
obs2 <- rbind(runif(8,0,1),runif(8,0,1))
data=list()
data[[1]]=obs1
data[[2]]=obs2
######### create guesses for gaussian covariance matrix and means #########
mc2 <- Mclust(t(cbind(obs1,obs2)), G=2)
muInit=mc2$parameters$mean;
SigmaInit=mc2$parameters$variance$sigma
#########create tree topology####################################
nodeIndices1=cbind(c(0,1,2,3,3,4,5,6,7,8),c(1,2,3,4,5,6,7,8,9,10));
nodeIndices2=cbind(c(0,1,2,3,4,4,5,6),c(1,2,3,4,5,6,7,8));
dataNodeIndices=list()
dataNodeIndices[[1]]=nodeIndices1;
dataNodeIndices[[2]]=nodeIndices2;
######### create guesses for prior and transition matrices#########
prior1=array(1,nHStates);
priorInit = array(1,nHStates)/nHStates;
transmatInit = matrix(1,nHStates,nHStates)*(1/nHStates);
transInitDiv= matrix(1,nHStates,nHStates*nHStates)*(1/(nHStates*nHStates))
#
hfm=HFMfit(data,nHStates,dataNodeIndices,priorInit,transmatInit,transInitDiv,'c',
SigmaInit=SigmaInit, muInit=muInit);
# calculate Viterbi path
viterbi=HFMviterbi(obs1,nHStates,nodeIndices1,hfm$initProb,hfm$transMatSeq,hfm$transMatDiv,
'c',sigma=hfm$sigma, mu=hfm$mu);
# draw the path
tree=DrawViterbiTree(viterbi$maxPath)
#plot the tree
#image(t(tree),col=rainbow(10))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
