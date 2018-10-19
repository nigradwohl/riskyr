## Interfacing table calculation with the riskyr function (development):



riskyr <- function(scen.lbl = "",  ## WAS: txt$scen.lbl,
                   scen.lng = txt$scen.lng,
                   scen.txt = txt$scen.txt, popu.lbl = txt$popu.lbl,
                   cond.lbl = txt$cond.lbl,
                   cond.true.lbl = txt$cond.true.lbl, cond.false.lbl = txt$cond.false.lbl,
                   dec.lbl = txt$dec.lbl,
                   dec.pos.lbl = txt$dec.pos.lbl, dec.neg.lbl = txt$dec.neg.lbl,
                   hi.lbl = txt$hi.lbl, mi.lbl = txt$mi.lbl,
                   fa.lbl = txt$fa.lbl, cr.lbl = txt$cr.lbl,
                   prev = NA, ## e.g., nprev = 1 - prev (specify in function call?)
                   ppod = NA,
                   sens = NA,
                   spec = NA, fart = NA,
                   N = NA,  ## WAS: freq$N,
                   hi = NA, mi = NA,
                   fa = NA, cr = NA,
                   scen.src = txt$scen.src,
                   scen.apa = txt$scen.apa
                   ## TODO: Deprecate but keep parameters that are not used anymore:
                   ## e.g., convert fart to FNR internally.
                   ) {

  ## Note:  this interface is specific to the diagnostic case (alias: riskyr.diagnostic):


  ## TODO: Decide on flexibility!
  ## Matrix for condition-decision case:
  ## This one should probably not be flexible to have clear variable names for internal use:
  cnd_dec_var <- matrix(c("hi", "mi", "cond.true", "sens", "FPR",
                          "fa", "cr", "cond.false", "FNR", "spec",
                          "dec.pos", "dec.neg", "N", "ppod", "pned",  # pned: proportion negative decisions.
                          "PPV", "FOR", "prev", "acc", NA,
                          "FDR", "NPV", "nprev", NA, NA),
                        ## TODO: Find good variable names for the missing quantities.
                        nrow = 5, ncol = 5)

  ## for user-defined display an additional flexible matrix may be used:
  cnd_dec_lab <- matrix(c(hi.lbl, mi.lbl, cond.true.lbl, "sens", "FPR",
                      fa.lbl, cr.lbl, cond.false.lbl, "FNR", "spec",
                      dec.pos.lbl, dec.neg.lbl, "N", "ppod", "1-ppod",
                      "PPV", "FOR", "prev", "acc", NA,
                      "FDR", "NPV", "1-prev", NA, NA),
                    nrow = 5, ncol = 5)

  ## TODO: Here we need reasonable defaults, but should allow users to customize!


  ## (0): Initialize some stuff: ------
  num_mat <- matrix(c(hi, mi, hi+mi, sens, NA,
                      fa, cr, fa+cr, fart, spec,
                      hi+fa, mi+cr, N, ppod, NA,
                      NA, NA, prev, NA, NA,
                      NA, NA, NA, NA, NA,
                      NA, NA, NA, NA, NA),
                    nrow = 5, ncol = 5
                    )

  ## TODO: Provide inputs for composite variables!

  numeric <- calc_tab(num_mat)  # calculate the table object (maybe use try catch?)

  ## TODO (in calc tab) name table types.

  object <- list(numeric = numeric,  # slot for all numeric content (how to split for different types?)

            vars = cnd_dec_var,  # slot for variable names to access elements
            global_labels = list(scen.lbl = scen.lbl, scen.lng = scen.lng, scen.txt = scen.txt,
                                 popu.lbl = popu.lbl, cond.lbl = cond.lbl),  # labels not pertaining to cells.
            source = list(scen.src = scen.src, scen.apa = scen.apa)  # source info.
            )

  # object <- list(scen.lbl = scen.lbl, scen.lng = scen.lng, scen.txt = scen.txt,
  #                popu.lbl = popu.lbl, cond.lbl = cond.lbl,
  #                cond.true.lbl = cond.true.lbl, cond.false.lbl = cond.false.lbl,
  #                dec.lbl = dec.lbl, dec.pos.lbl = dec.pos.lbl, dec.neg.lbl = dec.neg.lbl,
  #                hi.lbl = hi.lbl, mi.lbl = mi.lbl, fa.lbl = fa.lbl, cr.lbl = cr.lbl,
  #                prev = probs[1],
  #                sens = probs[2],
  #                spec = probs[4], fart = probs[5],
  #                N = N,
  #                hi = freqs$hi, mi = freqs$mi,
  #                fa = freqs$fa, cr = freqs$cr,
  #                scen.src = scen.src, scen.apa = scen.apa)
  #
  ## Add class riskyr:
  class(object) <- c("riskyr.diagnostic")

  return(object)

}


riskyr(hi = 1, mi = 3, fa = 4, cr = 5)


## testing method:
## Setting a general riskyr class:
  setClass("riskyr",
           slots = c(numeric = "numeric",
                     vars = "character",
                     labels = "character",
                     source = "character")
           )
## Setting a specific class for the diagnostic case:
  setClass("riskyr.diagnostic", contains = "riskyr",  # let riskyr.diagnostic inherit riskyr.
           slots = c(numeric = "numeric",
                     vars = "character",
                     labels = "character",
                     source = "character")
  )

## Set method for the generic class:
setMethod("summary", signature = "riskyr",
          definition = function(object){return(round(object$numeric[[1]], 2))})



## Overwrite the cite method to obtain source information:
setMethod("cite", signature = "riskyr",
          definition = function(keys){return(keys$source)})


## Overwrite the summary method to obtain a summary of the table:
setMethod("summary", signature = c("riskyr.diagnostic"),
          definition = function(object, relf = FALSE){

            ix <- ifelse(relf, 2, 1)  # decide whether to report relative frequencies.

            nums <- round(object$numeric[[ix]], 2)
            vars <- object$vars
            out <- paste0(vars, " = ", nums)

            head <- ifelse(relf, "relf", "freq")
            cat(head, "\n")
            return(out)
            })


## Test the class and method:
tst <- riskyr(hi = 1, mi = 3, fa = 4, cr = 5)

summary(tst, relf = TRUE)  # summary now returns the source info (obviously this is just a test).
cite(tst)  # cite now returns source information (also for object of class riskyr.diagnostic, inheriting riskyr)!
