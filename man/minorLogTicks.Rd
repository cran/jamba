% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jamba-axis.R
\name{minorLogTicks}
\alias{minorLogTicks}
\title{Calculate major and minor tick marks for log-scale axis}
\usage{
minorLogTicks(
  side = NULL,
  lims = NULL,
  logBase = 2,
  displayBase = 10,
  logStep = 1,
  minorWhich = c(2, 5),
  asValues = TRUE,
  offset = 0,
  symmetricZero = (offset > 0),
  col = "black",
  col.ticks = col,
  combine = FALSE,
  logAxisType = c("normal", "flip", "pvalue"),
  verbose = FALSE,
  ...
)
}
\arguments{
\item{side}{\code{integer} value indicating which axis to produce tick
marks, 1=bottom, 2=left, 3=top, 4=right.}

\item{lims}{\code{numeric} vector length=2, indicating specific numeric
range to use for tick marks.}

\item{logBase}{\code{numeric} value indicating the logarithmic base, assumed
to be applied to the numeric \code{lims} limits, or the axis range,
previously.}

\item{displayBase}{\code{numeric} value indicating the base used to position
axis labels, typically \code{displayBase=10} is used to draw labels
at typical positions.}

\item{logStep}{\code{integer} value indicating the number of log steps
between major axis label positions. Typically \code{logStep=1} will
draw a label every log position based upon \code{displayBase}, for
example \code{displayBase=10} and \code{logStep=1} will use \code{c(1,10,100,1000)};
and \code{displayBase=10} and \code{logStep=2} would use \code{c(1,100,10000)}.}

\item{minorWhich}{\code{integer} vector of values to label, where those
integer values are between 1 and \code{displayBase}, for example
\code{displayBase=10} may label only \code{c(2,5)}, which implies minor
tick labels at \code{c(2, 5, 20, 50, 200, 500)}. Any minor labels
which would otherwise equal a major tick position are removed.
By default, when \code{displayBase=2}, \code{minorWhich=c(1.5)} which has the
effect of drawing one minor label between each two-fold
major tick label.}

\item{asValues}{\code{logical} indicating whether to create exponentiated
numeric labels. When \code{asValues=FALSE}, it creates \code{expression} objects
which include the exponential value. Use \code{asValues=FALSE} and
\code{logAxisType="pvalue"} to draw P-value labels.}

\item{offset}{\code{numeric} value added during log transformation, typically
of the form \code{log(1 + x)} where \code{offset=1}. The offset is used to
determine the accurate numeric label such that values of \code{0} are
properly labeled by the original numeric value.}

\item{symmetricZero}{\code{logical} indicating whether numeric values
are symmetric around zero. For example, log fold changes should
use \code{symmetricZero=TRUE} which ensures a log2 value of \code{-2} is
labeled \code{-4} to indicate a negative four fold change. If
\code{symmetricZero=FALSE} a log2 value of \code{-2} would be labeled
\code{0.0625}.}

\item{col, col.ticks}{\code{character} color used for the axis label, and
axis tick marks, respectively, default "black".}

\item{combine}{\code{logical}, default FALSE, whether to combine major and
minor ticks into one continuous set of major tick marks.}

\item{logAxisType}{\code{character} string indicating the type of log axis:
\itemize{
\item normal: typical axis style and orientation
\item flipped: used for reverse orientation
\item pvalue: used for \code{-log10(pvalue)} orientation.
}}

\item{verbose}{logical indicating whether to print verbose output.}

\item{...}{additional parameters are ignored.}
}
\value{
\code{list} of axis tick positions, and corresponding labels,
for major and minor ticks. Note that labels may be \code{numeric},
\code{character}, or \code{expression}. Specifically when \code{expression}
the \code{graphics::axis()} must be called once per label.
\itemize{
\item majorTicks: \code{numeric} position of each major tick mark
\item minorTicks: \code{numeric} position of each minor tick mark
\item allTicks: \code{numeric} position of each major tick mark
\item majorLabels: label to show for each tick mark
\item minorLabels: label to show for each tick mark
\item minorSet: the \code{numeric} steps requested for minor ticks
\item minorWhich: the \code{numeric} steps requested for minor labels
\item allLabelsDF: \code{data.frame} with all tick marks and labels, with
colname \code{"use"} indicating whether the label is displayed beside
each tick mark.
}
}
\description{
Calculate major and minor tick marks for log-scale axis
}
\details{
This function is called by \code{minorLogTicksAxis()}, and
it may be better to use that function, or \code{logFoldAxis()}
or \code{pvalueAxis()} which has better preset options.

This function calculates log units for the axis of an
existing base R plot.
It calculates appropriate tick and label positions for:
\itemize{
\item major steps, which are typically in log steps; and
\item minor steps, which are typically a subset of steps at one
lower log order.
}

For example, log 10 steps would be: \code{c(1, 10, 100, 1000)},
and minor steps would be \code{c(2, 5, 20, 50, 200, 500, 2000, 5000)}.
\subsection{Motivation}{

This function is motivated to fill a few difficult cases:
\enumerate{
\item Label axis ticks properly
when used together with \code{offset}. For example \code{log2(1 + x)}
uses \code{offset=1}. Other offsets can be used as relevant.
\item Create axis labels which indicate negative fold change
values, for example \code{-2} in log2 fold change units would
be labeled with fold change \code{-4}, and not \code{0.0625}.
\item Use symmetric tick marks around x=0 when applied to log fold changes.
\item Display actual P-values when plotting \code{log10(Pvalue)}, which
is common for volcano plots.
}
}
}
\examples{
## This example shows how to draw axis labels manually,
## but the function minorLogTicksAxis() is easier to use.
xlim <- c(0,4);
nullPlot(xlim=xlim, doMargins=FALSE);
mlt <- minorLogTicks(1,
   logBase=10,
   offset=1,
   minTick=0);
maj <- subset(mlt$allLabelsDF, type \%in\% "major");
graphics::axis(1, las=2,
   at=maj$tick, label=maj$text);
min <- subset(mlt$allLabelsDF, type \%in\% "minor");
graphics::axis(1, las=2, cex.axis=0.7,
   at=min$tick, label=min$text,
   col="blue");
graphics::text(x=log10(1+c(0,5,50,1000)), y=rep(1.7, 4),
   label=c(0,5,50,1000), srt=90);

nullPlot(xlim=c(-4,10), doMargins=FALSE);
abline(v=0, lty=2)
graphics::axis(3, las=2);
minorLogTicksAxis(1, logBase=2, displayBase=10, symmetricZero=TRUE);

nullPlot(xlim=c(-4,10), doMargins=FALSE);
graphics::axis(3, las=2);
minorLogTicksAxis(1, logBase=2, displayBase=10, offset=1);
x2 <- stats::rnorm(1000) * 40;
d2 <- stats::density(log2(1+abs(x2)) * ifelse(x2<0, -1, 1));
lines(x=d2$x, y=normScale(d2$y)+1, col="green4");

nullPlot(xlim=c(0,10), doMargins=FALSE);
graphics::axis(3, las=2);
minorLogTicksAxis(1, logBase=2, displayBase=10, offset=1);
x1 <- c(0, 5, 15, 200);
graphics::text(y=rep(1.0, 4), x=log2(1+x1), label=x1, srt=90, adj=c(0,0.5));
graphics::points(y=rep(0.95, 4), x=log2(1+x1), pch=20, cex=2, col="blue");

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
\code{\link{newestFile}()},
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
