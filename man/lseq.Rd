\name{lseq}
\alias{lseq}
\title{Logarithmically spaced sequence}
\description{
  Generates \code{length.out} positive values whose logarithms are equally
  spaced from \code{log(from)/adj} to \code{log(to)} (equivalently, values
  follow a geometric progression on that log interval). Used inside
  \code{\link{SmoothWin}} to build default grids for bandwidth (\code{l})
  and sharpness (\code{k}).
}
\usage{
lseq(from = 1, to = 5, length.out = 6, adj = 1)
}
\arguments{
  \item{from}{
    Positive numeric; start of the sequence on the original scale when
    \code{adj = 1} (first generated value equals \code{from}). Must be
    strictly positive so that \code{log(from)} is defined.
  }
  \item{to}{
    Positive numeric; end of the sequence on the original scale when
    \code{adj = 1} (last generated value equals \code{to}). Must be
    strictly positive.
  }
  \item{length.out}{
    Non-negative integer; desired length of the result. Passed to
    \code{\link[base]{seq}} as \code{length.out}.
  }
  \item{adj}{
    Positive numeric; divisor applied to \code{log(from)} in the lower
    endpoint of the underlying arithmetic sequence in log space:
    \code{seq(log(from)/adj, log(to), length.out = length.out)}. The
    default \code{1} gives a geometric grid from \code{from} to \code{to}.
    Values other than \code{1} change the effective lower endpoint before
    exponentiation (so the first value is \code{exp(log(from)/adj)} rather
    than \code{from} when \code{from != 1}).
  }
}
\details{
  The implementation is
  \code{exp(seq(log(from)/adj, log(to), length.out = length.out))}.
  With \code{adj = 1}, values are equally spaced on a logarithmic scale
  between \code{from} and \code{to}, which is natural for multiplicative
  search grids.
}
\value{
  A numeric vector of length \code{length.out} with strictly positive
  entries (assuming \code{from}, \code{to}, and \code{adj} are positive and
  \code{length.out} is at least 1).
}
\author{
  Hamed Haselimashhadi <hamedhaseli@gmail.com>
}
\seealso{
  \code{\link{SmoothWin}} for use in default \code{l} and \code{k} grids.
}
\examples{
## Geometric sequence from 1 to 5 (compare to linear seq)
lseq(from = 1, to = 5, length.out = 6)

## Shorter grid similar to default sharpness values in SmoothWin
lseq(from = 0.5, to = 10, length.out = 8)
}
\keyword{math}
