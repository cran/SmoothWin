\name{plot.SmoothWin}
\alias{plot.SmoothWin}

\title{
 Plot function for the SmoothWin object
}
\description{
  This function plots a SmoothWin object
}
\usage{
  \method{plot}{SmoothWin}(x,
                           ylab   = 'Response'              ,
                           xlab   = 'Time (continuous)'     ,
                           sub    = NULL                    ,
                           col    = NULL                    ,
                           digits = 2                       ,
                           ...
                           )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  SmoothWin object
}
  \item{ylab}{
  Label on the y axis. Default `Response`
}
  \item{xlab}{
  Label on the x axis. Default `Time (continuous)`
}
\item{sub}{
  See the `sub` parameter in `plot()` function. If left NULL then some information about the final window will be shown. Default NULL
}
  \item{col}{
   Colour parameter for the points. Set to NULL to use the default colouring (spectrum colouring). Default NULL
  }
  \item{digits}{
   The number of visible digits for l, k and SWS. Default 2
  }
  \item{\dots}{
  Optional parameters that can be passed to the `plot`/`qqPlot` function. See `car` package for the qqPlot function
}
}
\value{
  Called mainly for its side effect of drawing a scatterplot of the response against time, overlaying the scaled window weights and vertical lines at the modes.

  If a final model is available, invisibly returns a list with components \code{weight} (numeric vector of full window weights along the time axis), \code{k} and \code{l} (the sharpness and bandwidth selection results from the fit, as stored in \code{x$final.k} and \code{x$final.l}), and \code{object} (the original \code{SmoothWin} object).

  If windowing failed and no model is stored, prints a message and does not plot; the return value is \code{NULL} in that case.
}
\author{
	Hamed Haselimashhadi <hamedhaseli@gmail.com>
}

\seealso{
	\code{\link{SmoothWin}}
}

\examples{
	example(SmoothWin)
}
