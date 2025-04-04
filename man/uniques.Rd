% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jamba-string.r
\name{uniques}
\alias{uniques}
\title{apply unique to each element of a list}
\usage{
uniques(
  x,
  keepNames = TRUE,
  incomparables = FALSE,
  useBioc = TRUE,
  useSimpleBioc = FALSE,
  xclass = NULL,
  ...
)
}
\arguments{
\item{x}{input list of vectors}

\item{keepNames}{boolean indicating whether to keep the list element
names in the returned results.}

\item{incomparables}{see \code{\link[=unique]{unique()}} for details, this value is only
sent to \code{S4Vectors::unique()} when the Bioconductor package
\code{S4Vectors} is installed, and is ignored otherwise for efficiency.}

\item{useBioc}{\code{logical}, default TRUE,  indicating whether this
function should try to use \code{S4Vectors::unique()} when the
Bioconductor package \code{S4Vectors} is installed, otherwise it will
use a somewhat less efficient bulk operation.}

\item{useSimpleBioc}{\code{logical}, default FALSE, whether to use a legacy
mechanism with \code{S4Vectors} and is maintained for edge cases where
it might be faster.}

\item{xclass}{\code{character} optional vector of classes, used to invoke
optimized logic when the class is known upfront.}

\item{...}{additional arguments are ignored.}
}
\value{
\code{list} with unique values in each list element.
}
\description{
Apply unique to each element of a list, usually a list of vectors
}
\details{
This function will attempt to use \code{S4Vectors::unique()} which is
substantially faster than any \code{apply} family function, especially
for very long lists. However, when \code{S4Vectors} is not installed,
it applies uniqueness to the \code{unlist}ed vector of values, which is
also substantially faster than the \code{apply} family functions for
long lists, but which may still be less efficient than the
C implementation provided by \code{S4Vectors}.
}
\examples{
L1 <- list(CA=nameVector(LETTERS[c(1:4,2,7,4,6)]),
   B=letters[c(7:11,9,3)],
   C2=NULL,
   D=nameVector(LETTERS[4]));
L1;
uniques(L1);

uniques(L1, useBioc=FALSE);

}
\seealso{
Other jam list functions: 
\code{\link{cPaste}()},
\code{\link{heads}()},
\code{\link{jam_rapply}()},
\code{\link{list2df}()},
\code{\link{mergeAllXY}()},
\code{\link{mixedSorts}()},
\code{\link{rbindList}()},
\code{\link{relist_named}()},
\code{\link{rlengths}()},
\code{\link{sclass}()},
\code{\link{sdim}()},
\code{\link{unnestList}()}
}
\concept{jam list functions}
