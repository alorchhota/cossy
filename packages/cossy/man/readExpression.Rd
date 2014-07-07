\name{readExpression}
\alias{readExpression}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read gene expression file.
}
\description{
Create R object from a gct file.
}
\usage{
readExpression(gctfile, chipfile = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gctfile}{
  Character string: the path to the gct file. File format is available \href{http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats}{here}.
}
  \item{chipfile}{
  Character string: the path to the microarray chip file. File format is available \href{http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats}{here}. See details.
}
}
\details{
 If chipfile is supplied, then it is used for probe to gene symbols mapping. Otherwise, the DESCRIPTION column in the gct file is used.
}
\value{
A dataframe is returned.
}
\keyword{ ~microarray expression}
\keyword{ ~gct}
