\name{cossy}
\alias{cossy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{COSSY function}
\description{COSSY algorithm.}
\usage{
cossy(expression, cls, misset, nmis = 5, pval.ent = F, sig.test = "iqr")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expression}{A data frame containing the gene expression of the samples.}
  \item{cls}{A data frame containg the class labels of the samples.}
  \item{misset}{MIS set to use for training.}
  \item{nmis}{Inteter. Number of MIS to return.}
  \item{pval.ent}{Logical. If True, p-value of entropy, instead of raw entropy score, is used to rank the MISs.}
  \item{sig.test}{Significance test to choose top probes in an mis. It can be either "iqr" or "ttest"}
}


\details{
cossy finds out the important molecular interaction subnetworks from misset. 

}

\value{
Returns a list containing following objects.
\item{topmis}{Top ranked MISs.}
\item{cls}{Class labels.}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~cossy }
