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
                          "PPV", "FOR", "prev", "acc", "",
                          "FDR", "NPV", "nprev", "", ""),
                        ## TODO: Find good variable names for the missing quantities.
                        nrow = 5, ncol = 5)

  ## for user-defined display an additional flexible matrix may be used:
  cnd_dec_lab <- matrix(c(hi.lbl, mi.lbl, cond.true.lbl, "sens", "FPR",
                      fa.lbl, cr.lbl, cond.false.lbl, "FNR", "spec",
                      dec.pos.lbl, dec.neg.lbl, "N", "ppod", "1-ppod",
                      "PPV", "FOR", "prev", "acc", "",
                      "FDR", "NPV", "1-prev", "", ""),
                    nrow = 5, ncol = 5)

  ## TODO: Here we need reasonable defaults, but should allow users to customize!


  ## (0): Initialize some stuff: ------
  ## Basic: dim 1 (cnd) and dim 2 (dec)
    num_mat1 <- matrix(c(hi, mi, hi+mi, sens, NA,
                        fa, cr, fa+cr, fart, spec,
                        hi+fa, mi+cr, N, ppod, NA,
                        NA, NA, prev, NA, NA,
                        NA, NA, NA, NA, NA,
                        NA, NA, NA, NA, NA),
                      nrow = 5, ncol = 5
                      )

    ## TODO: Provide inputs for composite variables!


    numeric1 <- calc_tab(num_mat1)  # calculate the table object (maybe use try catch?)

    ## TODO: In which order to calculate the tables?

  ## Calculate dim1 (cnd) and dim 3 (acc) --------------------------
    num_mat2 <- matrix(NA,
                      nrow = 5, ncol = 5
    )
    #num_mat2[4, 3] <- numeric1[[1]][4, 4]  # enter accuracy instead of ppod.
    num_mat2[3, 4] <- numeric1[[1]][3, 4]  # enter prevalence (needed, so that the table is specified!).
    num_mat2[4, 1:2] <- diag(numeric1[[1]][c(4, 5), c(1, 2)])  # enter sens and spec (note, that this reverses the columns!).
    num_mat2[3, 3] <- numeric1[[1]][3, 3]  # enter N; CAUTION: may not always work.

    numeric2 <- calc_tab(num_mat2)  # calculate dim 3 (acc).

    ## Important: Columnns need to be reveactually the reversal is sufficient to do the trick!

    cnd_acc_var <- matrix(c("hi", "mi", "cond.true", "sens", "FPR",
                            "cr", "fa", "cond.false", "spec", "FNR",
                            "cor", "inc", "N", "acc", "err",  # pned: proportion negative decisions.
                            "hi.cor", "mi.err", "prev", "ppod", "",
                            "cr.cor", "fa.err", "nprev", "", ""),
                          ## TODO: Find good variable names for the missing quantities.
                          nrow = 5, ncol = 5)

  ## TODO (in calc tab) name table types.

  ## Calculate dim 2 (dec) and dim 3 (acc):
    num_mat3 <- matrix(NA,
                       nrow = 5, ncol = 5
    )

    tnum1 <- t(numeric1[[1]])  # transpose whole table!
    t(cnd_dec_var)  # for testing: transposed variables.

    num_mat3[3, 4] <- tnum1[3, 4]  # enter ppod (needed, so that the table is specified!).
    num_mat3[4, 1:2] <- diag(tnum1[c(4, 5), c(1, 2)])  # enter PPV and NPV (note, that this reverses the columns!).
    num_mat3[3, 3] <- tnum1[3, 3]  # enter N; CAUTION: may not always work.

    numeric3 <- calc_tab(num_mat3)  # calculate dim 3 (acc).

    dec_acc_var <- matrix(c("hi", "fa", "dec.yes", "PPV", "FOR",
                            "cr", "mi", "dec.no", "NPV", "FDR",
                            "cor", "inc", "N", "acc", "err",  # pned: proportion negative decisions.
                            "hi.cor", "fa.err", "ppod", "prev", "",
                            "cr.cor", "mi.err", "pned", "", ""),
                          ## TODO: Find good variable names for the missing quantities.
                          nrow = 5, ncol = 5)


  ## TODO: Decide where to put numeric and related label info.  Current idea:
        ## - one list for each dimension (numeric, vars, labels)
        ## - global dictionary, which dimension is which? (alternatively: by order)
        ## - plus: global attributes; including colors for specific variables (color scheme); potentially labels

  ## Bind into one object:
  object <- list(numeric = list(
                            dim1_2 = numeric1,
                            dim1_3 = numeric2,
                            dim2_3 = numeric3
                          ),  # slot for all numeric content (how to split for different types?)

                 ## TODO: Use specific or generic names in the object?
                 ## I currently seem to have solved a rather generic case...
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



### Skeleton for a general function: --------------------------------

riskyr_gen <- function(x,  # some input: table, list of frequencies / probs, generically named parameters.
                       info,  # scenario information (probably several parameters.)
                       type = "diagnostic",  # overall type ("generic", "treatment", ...). Controls methods.
                       ## Dimension names of the table; defaults to diagnostic case (necessary?).
                       dim1 = "cnd",
                       dim2 = "dec",
                       dim3 = "acc"
                        ## Potentially also include labels.
){

}

