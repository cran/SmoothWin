gridSearchModel = function(object            ,
                           data              ,
                           t                 ,
                           m = mean(t)       ,
                           l = 1             ,
                           k = 1             ,
                           weightFUN = function(x, ...) {
                             x
                           }                   ,
                           threshold = 10 ^ -18,
                           check = 1           ,
                           messages = TRUE     ,
                           onlyOne  = FALSE    ,
                           ...) {
  lk = length(k)
  ll = length(l)
  m  = unique(m)
  ## outputs
  lmodel = lweight = list()
  rmat   = matrix(0, ncol = 4       , nrow = ll * lk)
  colnames(rmat)          = c('Ind',
                              'ObsInInterval',
                              'k',
                              'l')
  counter = countProgress = 1
  pb = txtProgressBar(
    min = 0,
    max = lk * ll,
    style = 3,
    width = 50
  )
  for (lp in l) {
    for (kp in k) {
      weight = expWeight(
        t = t            ,
        k = kp           ,
        l = lp           ,
        m = m            ,
        plot = 0
      )
      
      wi      = checkWeightsN(
        w         = weight    ,
        threshold = threshold ,
        t         = t         ,
        check     =  check
      )
      inn                 = wi$NoZeWe
      newdata             = data[wi$wInd,]
      newdata$ModelWeight = wi$w
      slm = tryCatch(
        expr = {
          do.call('update',
                  list(
                    object                        ,
                    weights  = weightFUN (newdata$ModelWeight)                ,
                    data     = droplevels(newdata),
                    ...
                  ))
        },
        error = function(err, messsages = messages) {
          if (messsages)
            message('\n Error: ', err)
          slm = NULL
        } ,
        warning = function(war, messsages = messages) {
          if (messsages)
            message('\n Warning: ', war)
          slm = NULL
        }
      )
      # Progress bar
      setTxtProgressBar(pb = pb, value = countProgress)
      countProgress = countProgress + 1
      # pass or skip
      if (is.null(slm)) {
        rmat = rmat[-nrow(rmat), ]
        next
      }
      
      lmodel[[counter]]   = slm
      lweight[[counter]]  = wi$w
      rmat  [counter, ]    = c(counter ,
                               sum(inn),
                               kp      ,
                               lp)
      counter           = counter  + 1
    }
  }
  close(pb)
  cat('\n')
  
  
  return(list(
    output = as.data.frame(rmat),
    models = if (onlyOne && length(lmodel) < 2) {
      lmodel[[1]]
    } else{
      lmodel
    },
    weights = if (onlyOne && length(lweight) < 2) {
      lweight[[1]]
    } else{
      lweight
    }
  ))
}
