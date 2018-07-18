.onAttach <- function(lib, pkg) {
  packageStartupMessage(
    paste0(
      '\n ========================================================================',
      '\n If you have any question about this package contact me on               ',
      '\n      hamedhm@ebi.ac.uk                                                  ',
      '\n *** This project is supported by European Bioinformatic Institute (EBI) ',
      '\n *** https://www.mousephenotype.org/                                     ',
      '\n ========================================================================'
    ),
    domain = NULL,
    appendLF = TRUE
  )
}
###########
expF = function(x, k, l, m) {
  m  = unique(m)
  r = plogis(x,
             location = m - l,
             scale = 1 / k) *
    (1 - (plogis(
      x,
      location = m + l,
      scale = 1 / k
    )))
  return(r)
}
###########
expWeight = function(t, k, l, m = 0, plot = FALSE , ...) {
  m  = unique(m)
  lm = length(m)
  if (lm > 10)
    message('more than 10 modes detected. It can take a long time to finish!')
  r  = matrix(0, ncol = length(m), nrow = length(t))
  for (i in 1:lm) {
    r[, i] = expF (x = t,
                   k = k,
                   l = l,
                   m = m[i])
  }
  s = rowSums(r)
  if (lm > 1) {
    for (j in 2:lm) {
      cm  = combn(lm, j)
      for (i in 1:ncol(cm)) {
        s = s + (-1) ^ (j + 1) * apply(r, 1, function(x) {
          prod(x[cm[, i]])
        })
      }
    }
  }
  s = s / ifelse(max(s, na.rm = TRUE) != 0, max(s, na.rm = TRUE), 1)
  if (plot) {
    plot(
      t,
      s,
      xlab = 'Time',
      ylab = 'Weight',
      xaxt = 'n',
      xlim = c(min(t, na.rm = TRUE) * .85, max(t, na.rm = TRUE) * 1.5),
      ...
    )
    abline(
      v = c(m),
      lty = 1,
      col = 1:lm,
      lwd = 2
    )
    abline(v = c(m + l, m - l),
           lty = 2,
           col = 1:lm)
    axis(
      side = 1,
      at = c(m, m - l, m + l),
      labels = c(m, m - l, m + l),
      lwd.ticks = 2
    )
    legend(
      'topright',
      title = 'm/m+l/m-l',
      legend = paste(round(m, 3), round(m + l, 3), round(m - l, 3), sep = '/'),
      col = 1:lm,
      fill = 1:lm
    )
  }
  return(s)
}
###########
checkWeightsN <- function(w             ,
                          t             ,
                          check = 1     ,
                          #1 all+nosingle, #2 all without nosingle #0 off
                          threshold = 10 ^ -18) {
  n  = length(t)
  if (!is.null(w) && check > 0) {
    #&& var(w) > threshold) {
    r      = 1:n
    zw     = which(abs(w) > threshold)
    #### need at least 2 observations in a group
    if (check == 1) {
      t      = as.character(t)
      tz     = table(t[zw])
      MorTh1 = names(tz)[which(tz > 1)]
      zw     = zw[t[zw] %in% MorTh1] # no singleDay-singleData
    }
    if (length(zw) > 0) {
      r  = r[zw]
      w  = w[zw] / sum(w[zw])
    } else{
      message(
        '\n * Model weights are ignored! ** weight may be all close to zero *** there may be all dates with single weight [can cause error in the mixed model] **** Setting check = 1 or check = 2 may solve the problem.\n'
      )
      w = (w * 0 + 1) / n
      r = 1:n
    }
  } else{
    r = 1:n
  }
  return(list(
    w = w,
    wInd = r,
    NoZeWe = sum(w > threshold)
  ))
}
###########
lseq = function (from = 1,
                 to = 5,
                 length.out = 6,
                 adj = 1)
{
  r = exp(seq(log(from) / adj, log(to), length.out = length.out))
  return(r)
}
###########
msg = function(args, ...) {
  message(
    '\nMinimum observations in the model:',
    '\n Min observations around each experiment date: ',
    args$min.obs / length(unique(args$t[args$m])),
    '\n Min observations required: ' ,
    args$min.obs + length(args$m),
    '\n  The number of modes: ',
    length(unique(args$t[args$m])),
    '\n  The time      range: ',
    paste(round(range(args$t), 3), collapse = ' to '),
    '\n  The bandwith  range: ',
    paste(round(range(args$l), 3), collapse = ' to '),
    ' [',
    length(args$l),
    ' split]',
    '\n  The shape     range: ',
    paste(round(range(args$k), 3), collapse = ' to '),
    ' [',
    length(args$k),
    ' split]',
    '\n'
  )
}
###########
tv.test = function(obj, args, name = 'parameter', ...) {
  if (length(obj$models) > 1) {
    df = lapply(obj$models, resid)
    tt = c()
    for (i in 2:length(obj$models)) {
      vtl = var.test(df[[i]][obj$weights[[i]] > args$threshold], df[[i - 1]][obj$weights[[i - 1]] > args$threshold])$p.value
      ttl = t.test  (df[[i]][obj$weights[[i]] > args$threshold], df[[i - 1]][obj$weights[[i - 1]] > args$threshold])$p.value
      tt[i - 1] = vtl * ttl + vtl  + ttl
    }
    #####
    
    al       = data.frame(
      t.pval = c(tt, Inf),
      ol     = obj$output$ObsInInterval,
      l      = obj$output[, name]
    )
    al       = al[al$ol > args$min.obs + length(args$m),]
    final   = al$l[which.min(al$t.pval)[1]]
    # hist(al$t.pval)
    # abline(v=final)
    if (length(final) < 1 || is.na(final)) {
      message(paste('\n An optimal', name, 'is not found. Max value will be used.'))
      final = NULL
    }
  } else{
    final = NULL
  }
  return(final)
}
