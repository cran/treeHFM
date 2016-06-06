\name{HFMviterbi}
\alias{HFMviterbi}
%- Also NEED an '\alia
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculates the most probable Hidden State path}
\description{
The function calculates the most probable hidden state path for an observation sequence.
} 
\usage{
HFMviterbi(observationSequence,nHStates,dataNodeIndices,priorInit,
transmatInit,transInitDiv,type,emissionProb=c(),sigma=c(), mu=c())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{observationSequence}{The observation sequence, either a matrix of size KxN, whereas K is the dimension of the data and N the lenght of the sequence or of size 1xK whereas the entries are discrete.}
  \item{nHStates}{Number of hidden states H}
  \item{dataNodeIndices}{The data structure. A list of the toplogy of every data entry of the format [parentIndex childIndex]. }
  \item{priorInit}{The state probabilities }
  \item{transmatInit}{ The transtition probabilities for sequential transitions (H x H) }
  \item{transInitDiv}{ The transtition probabilities for splitting events (H x (H*H)}
\item{type}{Type of the data. ('d' discrete or 'c' continous)}
\item{emissionProb}{Emission probilities for discrete data, a matrix of size E x H, whereas E are the number of discrete observable states.}
\item{sigma}{Covariance matrix for continous data}
\item{mu}{Means for continous data}
}
\details{
  The algorithm calculates the most probable hidden state path for an observation matrix.
}
\value{
A list of two entries:
Path matrix:
		1. Column: Node index 
		2. Column: Hidden state 
		3. Column: Branch of the tree
		4. Column: Parent branch
	
Accossiation of tree branches to parent branches:
	1. Column: Parent branch
	2. Column child branch.
 	
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
