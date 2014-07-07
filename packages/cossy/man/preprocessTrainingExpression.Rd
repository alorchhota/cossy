\name{preprocessTrainingExpression}
\alias{preprocessTrainingExpression}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Preprocess training expression.
}
\description{
Preprocess training expression.
}
\usage{
preprocessTrainingExpression(expression, frank = F, qnorm = T, ztrans = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expression}{
  A data frame returned by readExpression() function.
}
  \item{frank}{
  Boolean: Indicates whether the expression values should be normalized by fuzzy ranking.
}
  \item{qnorm}{
  Boolean: Indicates whether quantile normalization should be applied.
}
  \item{ztrans}{
  Boolean: Indicates whether z-transform should be applied.
}
}
\details{
  Training data often require different kinds of preprocessing (for e.g., quantile normalization, rank normalization, z-transformation etc.). The order of preprocessing is fuzzy rank normalization, quantile normalization, z-transformation. Fuzzy rank normalization is performed as described in Lim et al. (2014).
}
\value{
Returns a list containing the following attribtutes.
\item{expression }{Processed expression.}
\item{qnorm.quantiles }{Quantiles used in quantile normalization.}
\item{ztrans.centers}{Centers used in z-transformation.}
\item{ztrans.scales}{Scales used in z-transformation.}
\item{frank}{Ranks used in fuzzy ranking.}
}
\references{
Lim, K. and Wong. L. (2014) Finding consistent disease subnetworks using PFSNet. Bioinformatics, 20, 189-96.
}
\keyword{ ~normalization }
\keyword{ ~preprocess }
