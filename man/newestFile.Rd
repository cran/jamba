% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jamba.r
\name{newestFile}
\alias{newestFile}
\title{Return the newest file from a vector of files}
\usage{
newestFile(x, timecol = "mtime", n = 1, ...)
}
\arguments{
\item{x}{\code{character} vector of files, specifying file path where
required.}

\item{timecol}{\code{character} value from the output of \code{base::file.info()}
indicating the time column used to order files. By default \code{"mtime"}
refers to the time the file was last modified.}

\item{n}{\code{integer} number of files to return, in order of the most
recent to the least recent. By default \code{n=1} returns only the one
newest file.}

\item{...}{additional parameters are ignored.}
}
\value{
\code{character} vector \code{length=1} of the most recently modified file
from the input vector \code{x}. Note that any files not found are removed,
using \code{base::file.exists()}, which means invalid symlinks will be ignored.
}
\description{
Return the newest file from a vector of files
}
\details{
This function returns the newest file, defined by the most
recently modified time obtained from \code{base::file.info()}.
}
\examples{
newestFile(list.files());

}
\seealso{
Other jam practical functions: 
\code{\link{breakDensity}()},
\code{\link{call_fn_ellipsis}()},
\code{\link{checkLightMode}()},
\code{\link{check_pkg_installed}()},
\code{\link{colNum2excelName}()},
\code{\link{color_dither}()},
\code{\link{exp2signed}()},
\code{\link{getAxisLabel}()},
\code{\link{isFALSEV}()},
\code{\link{isTRUEV}()},
\code{\link{jargs}()},
\code{\link{kable_coloring}()},
\code{\link{lldf}()},
\code{\link{log2signed}()},
\code{\link{middle}()},
\code{\link{minorLogTicks}()},
\code{\link{printDebug}()},
\code{\link{reload_rmarkdown_cache}()},
\code{\link{renameColumn}()},
\code{\link{rmInfinite}()},
\code{\link{rmNA}()},
\code{\link{rmNAs}()},
\code{\link{rmNULL}()},
\code{\link{setPrompt}()}
}
\concept{jam practical functions}
