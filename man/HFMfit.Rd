\name{HFMfit}
\alias{HFMfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fit a HFM Model}
\description{
The functions fits an HFM (Hidden Factor Graph Model) to the (tree structured) data.
}
\usage{
	HFMfit(data,nHStates,dataNodeIndices,priorInit,transmatInit,transInitDiv,
	type,emissionProb=c(),SigmaInit=c(), muInit=c() )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{The observations. A list of measurements for every observation. An entry can be a matrix of size KxN, whereas K is the dimension of the data and N the lenght of the sequence or of size 1xK whereas the entries are discrete.}
  \item{nHStates}{Number of hidden states H}
  \item{dataNodeIndices}{The data structure. A list of the toplogy of every data entry of the format [parentIndex childIndex]. }
  \item{priorInit}{The inital state probabilities }
  \item{transmatInit}{ The initial transtition probabilities for sequential transitions (H x H) }
  \item{transInitDiv}{ The initial transtition probabilities for splitting events (H x (H*H)}
\item{type}{Type of the data. ('d' discrete or 'c' continous)}
\item{emissionProb}{Initial emission probilities for discrete data, a matrix of size E x H, whereas E are the number of discrete observable states.}
\item{SigmaInit}{Initial covariance matrix for continous data}
  \item{muInit}{Initial means for continous data}
}
\details{
  The algorithm fits an HFM to sequential data (in this case essentially a Hidden Markov Model is fitted) or to data that is structured as a binary tree. An HFM consists of two transition probability matrices. One for the sequential transitions and one for transitions where a splitting happens.
}
\value{
An HFM:
\item{$initProb}{The initial probabilities (prior).}
\item{$transMatSeq}{The transition matrix for sequential events}
\item{$transMatDiv}{The transition matrix for division events}
\item{$ll}{The log likelihood}
\item{$emission}{The emission probabilities of the hidden states}
}
\author{ Henrik Failmezger, failmezger@googlemail.com}


\examples{
## Fit a discrete treeHFM ## 
nOStates = 3;
nHStates = 2;
#
T = 1;
nex = 10;
obs1=sample(1:nOStates,5,replace=T)
obs2=sample(1:nOStates,5,replace=T)
data=list()
data[[1]]=obs1;
#
nodeIndices1=cbind(c(0,1,2,2,4),c(1,2,3,4,5));
dataNodeIndices=list()
dataNodeIndices[[1]]=nodeIndices1;
prior1=array(1,nHStates);
priorInit = array(1,nHStates)/nHStates;
transmatInit = matrix(1,nHStates,nHStates)*(1/nHStates);
transInitDiv= matrix(1,nHStates,nHStates*nHStates)*(1/(nHStates*nHStates))
obsmatInit = matrix(sample(1:nOStates,nHStates*nOStates,replace=T),nHStates,nOStates);
obsmatInit=obsmatInit/rowSums(obsmatInit,1);
hfm=HFMfit(data,nHStates,dataNodeIndices,priorInit,transmatInit,transInitDiv,
'd',emissionProb=obsmatInit);
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
#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
