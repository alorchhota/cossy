\name{predict}
\alias{predict}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cossy prediction.
}
\description{
Predicts the class label of a new sample.
}
\usage{
predict(cossyobj, expression)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cossyobj}{
	Cossy model. An R object returned by cossy() or cossy.v() function.
}
  \item{expression}{
	Test expression. A data frame returned by readExpression() function.
}
}
\details{
 Predicts the class label of a new sample using the model built by cossy() or cossy.v() function.
}
\value{
Returns a vector containint two items - 1) the predicted class label, 2) vote for positive class in [0,1] scale.
}
\keyword{ ~test }
\keyword{ ~predict }
