\name{readClass}
\alias{readClass}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read class labels.
}
\description{
Reads the class labels of the samples. 
}
\usage{
readClass(clsfile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clsfile}{
  Character string: the path to the class file. File format is available \href{http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats}{here}.
}
}
\details{
  The number of labels in the class file must be same as that in the gct file.
}
\value{
Returns a data frame.
}
\keyword{ ~cls file }
