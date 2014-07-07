\name{cossy.v}
\alias{cossy.v}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cross validated cossy.
}
\description{
Find consistent subnetworks using a validation step.
}
\usage{
cossy.v(expression, cls, misset, nmis = seq(1, 15, 2), one.se = T, mis.consistency = T, sig.test = "iqr")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expression}{A data frame containing the gene expression of the samples.}
  \item{cls}{A data frame containg the class labels of the samples.}
  \item{misset}{MIS set to use for training.}
  \item{nmis}{Inteter vector. The possbile number of MISs.}
  \item{one.se}{Logical. Indicates whether one standard error rule should be used.}
  \item{mis.consistency}{Indicates whether the consistent MISs should be used for final model.}
  \item{sig.test}{Significance test to choose top probes in an mis. It can be either "iqr" or "ttest"}
  
}
\details{
  Cossy.v function serves the same purpose as cossy with a validation step. It chooses the number of MISs from the nmis argument which produces the best accuracy in 10-fold cross-validation step. If mis.consistency=TRUE, the final model is built using the most consistent MISs i.e. the most frequent MISs in validation steps. If one.se=TRUE, the final model uses the simples model producing accuracy withing one standard error of the highest accuracy.
}
\value{
Returns a list containing following objects.
\item{topmis}{Top ranked MISs.}
\item{cls}{Class labels.}

}
\keyword{ ~consistent }
