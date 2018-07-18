\name{expWeight}
\alias{expWeight}
\title{
  This function computes the smooth windowing weights
}
\description{
  This function computes the exponential weights (smooth windowing function) for different shapes (k)  and bandwidth (l) and provides some plots
}
\usage{
  expWeight(t, k, l, m = 0, plot = FALSE, ...)
}
\arguments{
  \item{t}{
    Vector. a vector of positive continuous values 
  }
  \item{k}{
    A single value for sharpness
  }
  \item{l}{
    A single value for bandwidth
  }
  \item{m}{
    Vector. The location of the modes on t (modes are the peak of the windows)
  }
  \item{plot}{
    Logical flag. Set to true shows a plot of the weights.
  }
  \item{\dots}{
    The parameters that can be passed to the plot() function
  }
}
\value{
  A vector of weights
}

\author{
  Hamed Haselimashhadi <hamedhm@ebi.ac.uk>
}

\seealso{
  \code{\link{SmoothWin}}
}
\examples{
  par(mfrow = c(3, 1))
  # Example 1 - no merging in windows
  weight = expWeight(
    t = 1:100                                       ,
    k = 5                                           ,
    l = 10                                          ,
    m = c(25, 50, 75)                               ,
    plot = TRUE                                     ,
    ### Passed parameters to the plot function
    type = 'l'                                      ,
    lty = 2                                         ,
    lwd = 3                                         ,
    main = 'If windows do not intersect, then wont merge!'
  )
  # Example 2 - merging in windows
  weight = expWeight(
    t = 1:100                                       ,
    k = 5                                           ,
    l = 15                                          ,
    m = c(25, 50, 75)                               ,
    plot = TRUE                                     ,
    ### Passed parameters to the plot function
    type = 'l'                                      ,
    lty = 2                                         ,
    lwd = 3                                         ,
    main = 'If windows intersect, then merge!'
  )
  # Example 3 - partial merging in windows
  weight = expWeight(
    t = 1:100                                       ,
    k = 1                                           ,
    l = 12                                          ,
    m = c(25, 50, 75)                               ,
    plot = TRUE                                     ,
    ### Passed parameters to the plot function
    type = 'l'                                      ,
    lty = 2                                         ,
    lwd = 3                                         ,
    main = 'If windows intersect with small k, then partially merge!'
  )
}

