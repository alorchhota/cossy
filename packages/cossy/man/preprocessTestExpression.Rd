\name{preprocessTestExpression}
\alias{preprocessTestExpression}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Process test data.
}
\description{
The test data must be processed in the same way as the training data were processed.
}
\usage{
preprocessTestExpression(preprocessObj, expression)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{preprocessObj}{
  An R object returned by the preprocessTrainingExpression() function.
}
  \item{expression}{
  Training expression. This is a data frame read using a readExpression() function.
}
}
\details{
  Test data is processed in the same way as the trainign data were processed.
}
\value{
Returns a list containing the following objects:
\item{expression}{Processe expression.}
}
\keyword{ ~preprocess }
\keyword{ ~test data}
