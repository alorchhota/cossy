\name{icossy}
\alias{icossy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cossy in Internet.
}
\description{
This package helps web development of cossy.
}
\usage{
icossy(gctfile, chipfile = NA, clsfile, network, nmis, frank = T, qnorm = F, ztrans = F, sig.test = "ttest", mis.consistency=T)
}

\arguments{
  \item{gctfile}{ Character string: path to gct file}
  \item{chipfile}{ Character string or NA. path to chip file.}
  \item{clsfile}{Character string: path to class label file.}
  \item{network}{ Network object.}
  \item{nmis}{Numeric: number of MIS.}
  \item{frank}{Logical: use fuzzy rank normalization.}
  \item{qnorm}{Logical: use quanlite normalization.}
  \item{ztrans}{Logical: use z-transformation.}
  \item{sig.test}{Character string (either 'ttest' or 'iqr')}
  \item{mis.consistency}{Logical: should return consistent mis?}
}
\details{
The format of gctfile, chipfile and clsfile are available \href{http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats}{here}. If chipfile == NA, the mapping from row id (in gct file) to genes is done using the DESCRIPTION column of the gct file, which should contain the corresponding gene symbols.

network is an R object read from a gmt file (\href{http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats}{format}). The package contains sample network objects.

nmis is the number of top MISs to rank. If it is a numeric vector, the final model is built using the #mis that produces highest accuracy in validation step.

}
\value{
Returns a list. If there is no error in running the function, it contains the following items:

\item{status}{OK.}
\item{network}{JSON representation of top molecular interaction subnetowrks.}

If error occurs, the output contains the following items:

\item{status}{ERROR.}
\item{error}{The error message.}

}

\keyword{ ~icossy }
