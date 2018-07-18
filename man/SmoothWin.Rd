\name{SmoothWin}
\alias{SmoothWin}
\title{
Implementation of the soft windowing on the linear models
}
\description{
  Implementation of the soft windowing on a range of methods by imposing weights on the model. The parameters k and l control the shape and bandwidth of the windowing function. 
}
\usage{
SmoothWin(object                                          ,
         data                                             ,
         t                                                ,
         m                                                ,
         l = function(ignore.me.in.default) {
           r = SmoothWin:::lseq(
             from = 1                                          ,
             to = max(abs(t[m] - min(t)), abs(t[m] - max(t)),1),
             length.out = min(200, max(1,diff(range(t))))
           )
           r = unique(round(r))
           return(r)
         }                                                ,
         k = SmoothWin:::lseq(from = .5                   ,
                              to = 10                     ,
                              length.out = 30)            ,
         min.obs   = function(ignore.me.in.default) {
           lutm = length(unique(t[m]))
           r = ifelse(lutm > 1, 35, max(pi * sqrt(length(t)), 35))
           r = max(r * lutm, length(m), na.rm = TRUE)
           r = min(r, length(t), na.rm = TRUE)
           return(r)
         }                                                ,
         weightFUN = function(x) {
           x
         }                                                ,
         check = 2                                        ,
         threshold = 10 ^ -12                             ,
         messages  = TRUE,
         seed      = 123456                               ,
         simple.output = FALSE                            ,
         ...)
    }
\arguments{
  \item{object}{
  The fitted model. The object must support "update(weights =)". See examples
}
  \item{data}{
  data.frame. Input data that is used to fit the initial model
}
  \item{t}{
  Vector of (continues) time values.
}
  \item{m}{
  Vector of integers. Mode indices on the time component. For example 10, 11, 12. Note that it is different from t[10], t[11], t[12]
}
  \item{l}{
  Vector of continues values. Bandwidth parameter. The default uses the maximum distance of the modes (t[m]) from the time boundaries, max(max(t)-t[m],t[m]-min(t)) split on 200 points on the logarithmic scale.
}
  \item{k}{
  Vector of continues values. Shape parameter. The default uses 30 split of the values from 0.5 to 10 on the logarithmic scale.
}
  \item{min.obs}{
  Single value. The minimum observations that must be in the window. The default uses the following steps.\cr
    1. If there are more than one modes in the data, then: 35*(the number of modes)\cr
    2. If there is a single mode in the data, then: max(pi * sqrt(length(t)), 35)\cr
  * min.obs must not be less than the number of observations on the mode times\cr
  ** min.obs must be less than the total observations
}
  \item{weightFUN}{
  Weight function. By default, a vector of weights call "ModelWeight" is passed to this function. See the examples. 
}
  \item{check}{
  Single integer in \{0,1,2\}. \cr
  check=1, the function selects the times (t) with more than one observations. Further, the function only selects the values with weights greater than the threshold parameter. Mostly useful in fitting linear mixed models\cr
  check=2 (default), the function only selects the values with weights greater than the threshold parameter\cr 
  check=0, disable all checks
}
  \item{threshold}{
  Single continues value. The minimum value for weights before pushing them into zero. Default 10^-12
}
  \item{messages}{
  Logical value. Set to TRUE (default) to see the messages and errors 
}
  \item{seed}{
  seed. Default 123456
}
  \item{simple.output}{
  Logical flag. Setting to TRUE leads to not exporting the list of models for l and k. Useful to prevent memory leak.
}
  \item{\dots}{
  Other parameters that can be passed to the weightFUN()
}
}


\value{
  \item{final.k, final.l }{Final values for k and l}
  \item{finalModel, model.l, model.k }{List of models for l, k and the final model. }
  \item{others}{The input parameters such as x, y, t and so on}
}

\author{
Hamed Haselimashhadi <hamedhm@ebi.ac.uk>
}


\seealso{
  \code{\link{expWeight}}
}

\examples{
# All examples use the Orthodont dataset in the nlme package
library(nlme)
# Sort the data on the time component (age)
Orthodont =  Orthodont[order(Orthodont$age),]
# Modes
mode = which(Orthodont$age \%in\%  c(12))
# Time component
time = Orthodont$age

#################### Examples ####################
### Example 1. Linear model
# Method 1 (recommanded)
f  = formula(distance ~ Sex)
lm = do.call('lm', list(formula = f, data = Orthodont))
rm(f)

# Method 2 (can cause error if you pass the formula to the lm function)
# lm = lm(distance ~ Sex, data = Orthodont)

lm.result = SmoothWin(
  object = lm,
  data = Orthodont,
  t = time,
  m = mode,
  check = 0,
  weightFUN = function(x) {
    x
  }
)
plot(
  lm.result,
  col = Orthodont$Sex,
  pch = as.integer(Orthodont$Sex),
  main = 'Simple liner model'
)


#### Example 2. Linear Model Using Generalized Least Squares
# Method 1 (recommanded)
f   = formula(distance ~ Sex)
gls = do.call('gls', list(model = f, data = Orthodont))
rm(f)

# Method 2 (can cause error if you pass the formula to the gls function)
# gls = gls(distance ~ Sex, data = Orthodont)

gls.result = SmoothWin(
  object = gls,
  data = Orthodont,
  t = time,
  m = mode,
  check = 0,
  weightFUN = function(ignore.me) {
    varFixed( ~ 1 / ModelWeight) #nlme package uses the inverse weights
  }
)
plot(
  gls.result,
  col = Orthodont$Sex,
  pch = as.integer(Orthodont$Sex),
  main = 'Linear model using GLS'
)

#### Example 3. Linear mixed model
# Method 1 (recommanded)
fixed  = formula(distance ~ Sex)
random = formula( ~ 1 | Subject)
lme = do.call('lme', list(
  fixed = fixed,
  random = random,
  data = Orthodont
))
rm(fixed, random)

# Method 2 (can cause error if you pass the formula to the lme function)
# lme = lme(fixed = distance ~ Sex, random=~1|Subject , data = Orthodont)

lme.result = SmoothWin(
  object = lme,
  data = Orthodont,
  t = time,
  m = mode,
  # Remove zero weights as well as single observation dates
  check = 1,
  weightFUN = function(ignore.me) {
    varFixed( ~ 1 / ModelWeight)
  }
)
plot(
  lme.result,
  col = Orthodont$Sex,
  pch = as.integer(Orthodont$Sex),
  main = 'Linear mixed model'
)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{"Windowing", "SmoothWin" }% use one of  RShowDoc("KEYWORDS")
