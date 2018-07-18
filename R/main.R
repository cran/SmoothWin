SmoothWin = function(object                                           ,
                     data                                             ,
                     t                                                ,
                     m                                                ,
                     l = function(ignore.me.in.default) {
                       r = SmoothWin:::lseq(
                         from = 1                                          ,
                         to = max(abs(t[m] - min(t)), abs(t[m] - max(t)), 1),
                         length.out = min(200, max(1, diff(range(
                           t
                         ))))
                       )
                       r = unique(round(r))
                       return(r)
                     },
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
                     ...) {
  set.seed(seed)
  min.obs = ceiling(ifelse(is.function(min.obs), min.obs(), min.obs))
  l = if (is.function(l)) {
    l()
  } else{
    l
  }
  k = if (is.function(k)) {
    k()
  } else{
    k
  }
  k       = sort(k[!is.na(k)], decreasing = TRUE)
  l       = sort(l[!is.na(l)])
  argg    = c(as.list(environment()), list())
  
  if (length(unique(t[m])) > 10)
    message('More than 10 modes detected. The entire procedure can take a long time!')
  
  if (length(m) > min.obs) {
    stop('min.obs is less than the number of treatments!')
  } else {
    msg(argg)
  }
  ### 1. Determining l
  message('\n 1|3 Searching for l ...\n')
  rl = gridSearchModel(
    object = object           ,
    data = data               ,
    weightFUN = weightFUN     ,
    check = check             ,
    t = t                     ,
    m = t[m]                  ,
    l = l                     ,
    k = max(k)                ,
    threshold = threshold     ,
    messages = messages       ,
    onlyOne  = FALSE          ,
    ...
  )
  finall = tv.test(rl, argg, 'l')
  if (is.null(finall))
    finall = max(l)
  ### 2. Determining k
  message('\n 2|3 Searching for k ... \n')
  rk = gridSearchModel(
    object = object           ,
    data = data               ,
    weightFUN = weightFUN     ,
    check = check             ,
    t = t                     ,
    m = t[m]                  ,
    l = finall                ,
    k = k                     ,
    threshold = threshold     ,
    messages = messages       ,
    onlyOne  = FALSE          ,
    ...
  )
  finalk = tv.test(rk, argg, 'k')
  if (is.null(finalk))
    finalk = max(k)
  
  ##### final model
  message('\n 3|3 Forming the final model \n')
  finalr = gridSearchModel(
    object = object           ,
    data = data               ,
    weightFUN = weightFUN     ,
    check = check             ,
    t = t                     ,
    m = t[m]                  ,
    l = finall                ,
    k = finalk                ,
    threshold = threshold     ,
    messages = messages       ,
    onlyOne  = TRUE           ,
    ...
  )
  if (simple.output) {
    rk = rl = NULL
  }
  out = list(
    object  = object,
    data    = data  ,
    final.k = finalk,
    final.l = finall,
    finalModel = finalr,
    model.l = rl       ,
    model.k = rk       ,
    min.obs = min.obs  ,
    input   = argg
  )
  class(out) = 'SmoothWin'
  return(out)
}



plot.SmoothWin = function(x, ylab = 'Response', ...) {
  t = x$input$t
  y = x$data[, all.vars(formula(x$finalModel$models))[1]]
  ly = length(y)
  m = x$input$m
  plot(
    t,
    y,
    xlab = 'Time',
    ylab = ylab,
    sub = paste(
      'l=',
      round(x$final.l, 2),
      ', k=',
      round(x$final.k, 2),
      ', #=',
      x$finalModel$output$ObsInInterval,
      ' [',
      round(x$finalModel$output$ObsInInterval / ly * 100),
      '%]',
      ', MaxBW=',
      max(x$input$l, na.rm = TRUE),
      ', MinObs=',
      x$min.obs,
      sep = ''
    ),
    ...
  )
  
  abline(v = unique(t[m]),
         lty = 2,
         col = 'gray')
  wp = SmoothWin::expWeight(
    t = t,
    k = x$final.k,
    l = x$final.l,
    m = unique(t[m]),
    plot = 0
  )
  lines(
    t,
    min(y) + wp * (max(y) - min(y)),
    col = 'gray',
    lty = 4,
    lwd = 4
  )
  #####
}