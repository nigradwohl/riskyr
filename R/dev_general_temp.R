## TEMPORARY FILE FOR DEVELOPING GENERALIZATION FUNCTIONS:


## A. Calculating frequencies: --------------------------

library(riskyr)

## TODO: Function does not calculate sums currently! Necessary?

## TODO: Is a separation of rows and cols useful or is one sufficient with transposing?

## 1. Function to calculate rows: --------------------------
  crow <- function(ftab) {
    ## Get indices of NA rows:
    rcalc_ix <- which(is.na(ftab), arr.ind = TRUE)
    #rcalc_ix <- rcalc_ix[order(rcalc_ix[, 1]), ]
    dupl_ix <- !duplicated(rcalc_ix[, 1]) & !duplicated(rcalc_ix[, 1], fromLast = TRUE)
    # index for duplicates.
    rcalc_ix <- rbind(rcalc_ix[dupl_ix, ])  # remove duplicated rows.

    ## Get all missings not in the last (sum) column:
    if (any(rcalc_ix[,2] < 3)) {

      rmix1 <- rbind(rcalc_ix[rcalc_ix[, 2] < 3, c(1, 2)])  # cbind ensures matrix.
      rmix1

      ## Insert missings into appropriate places; incalculable columns remain NA:
      ftab[rmix1] <- ftab[rmix1[, 1], 3] - rowSums(rbind(ftab[rmix1[, 1], -3]), na.rm = TRUE)
    }

    ## Get all missings in the last (sum) column:
    if (any(rcalc_ix[,2] == 3)) {

      rmix2 <- rbind(rcalc_ix[rcalc_ix[, 2] == 3, c(1, 2)])  # cbind ensures matrix.
      rmix2

      ftab[rmix2] <- rowSums(rbind(ftab[rmix2[, 1], -3]))
    }

    return(ftab)
  }

## 2. Function to calculate columns: ---------------------
  ccol <- function(ftab) {
    ## Get indices of NA rows:
    ccalc_ix <- which(is.na(ftab), arr.ind = TRUE)
    dupl_ix <- !duplicated(ccalc_ix[, 2]) & !duplicated(ccalc_ix[, 2], fromLast = TRUE)
    # index for duplicates.
    ccalc_ix <- rbind(ccalc_ix[dupl_ix, ])  # remove duplicated rows.
    # ccalc_ix <- rcalc_ix[order(rcalc_ix[, 1]), ]

    ## Get all missings not in the last (sum) column:
    if (any(ccalc_ix[,1] < 3)) {

      cmix1 <- rbind(ccalc_ix[ccalc_ix[, 1] < 3, c(1, 2)])  # cbind ensures matrix.
      # mix1

      ## Insert missings into appropriate places; incalculable columns remain NA:
      ftab[cmix1] <- ftab[3, cmix1[, 2]] - colSums(cbind(ftab[-3, cmix1[, 2]]), na.rm = TRUE)
    }

    ## Get all missings in the last (sum) column:
    if (any(ccalc_ix[,2] == 3)) {

      cmix2 <- rbind(ccalc_ix[ccalc_ix[, 1] == 3, c(1, 2)])  # cbind ensures matrix.
      # mix2

      ftab[cmix2] <- colSums(cbind(ftab[-3, cmix2[, 2]]))
    }

    return(ftab)
  }

## 3. Function to complete a frequency table: --------------------
  comp_ftab <- function(ftab) {

    ## Get the dimensions:
      n_rowf <- nrow(ftab)
      n_colf <- ncol(ftab)

      ## (1) Check for rows and columns that can be calculated: -----
      ## Logical index matrix for NA values:
      mf_na_log <- is.na(ftab)  # get matrix of NAs for the frequency case.

      ## Stopping point 1:
      if (!any(mf_na_log)) {
        return(ftab)  # return, when table is complete.
      }


    ## (a) Rows: -----
      if(any(rowSums(mf_na_log) == 1)) {

        ## If any row can be calculated:
        ftab <- crow(ftab)  # calculate the rows.
      }


    ## (b) Columns: -----

      ## Renew log_na:
      mf_na_log <- is.na(ftab)  # get matrix of NAs

      ## Stopping point 2:
      if (!any(mf_na_log)) {
        return(ftab)  # return, when table is complete.
      }

      if(any(colSums(mf_na_log) == 1)) {
        ftab <- ccol(ftab)  # calculate the columns.
      }



    ## XX: Return statement and stopping rules: --------
    mf_na_log <- is.na(ftab)  # get matrix of NAs after one run.

    ## Recursive STOPPING rule 1:
    if (!any(mf_na_log)) {
      return(ftab)  # return, when table is complete.

      ## STOPPING RULE 2:
    } else if(!any(rowSums(mf_na_log) == 1 | colSums(mf_na_log) == 1)) {
      ## if no single value in neither dimension.

      # stop("The information you specified is not sufficient to calculate the full table.\nI need at least one more value."      )
      ## Do not throw the error in this case!
      return(ftab)  # return the incomplete table!

      ## CONTINUE if not complete:
    } else {

      return(comp_ftab(ftab))
    }

  }

## B. Probability table: ------------------------------

  ## 1. If the frequency table is complete, calculate probabilities:

  ## (a) Calculate all complements (for 2 x n matrix): ------
    compr_pcomp <- function(p_tab, transpose = FALSE, tol = 0.001) {

      ## Checking for completeness: -------

      if (transpose) p_tab <- t(p_tab)  # shouuld the table be transposed?

      ## One probability table as input:
        mp_na_log <- is.na(p_tab)  # get NA values.
        p_inc <- any(mp_na_log)  # check for any missings.

      ## Calculate complemets for rows: -----
        if(p_inc) {  # if row elements are missing.

          row_na <- which(mp_na_log, arr.ind = TRUE)  # get rows with NA-values.
          row_na <- rbind(
            row_na[!(duplicated(row_na[,1]) | duplicated(row_na[,1], fromLast = TRUE)), ]
            )  # remove duplicated rows (which can't be calculated).

          ## Complete row matrix:
          p_tab[row_na] <- 1 - rowSums(rbind(p_tab[row_na[, 1], ]), na.rm = TRUE)

          ## Check summing up to one:
          rows_eq1 <- isTRUE(all.equal(
            rowSums(rbind(p_tab[row_na[,1], ])),
            rep(1, nrow(row_na)),
            tolerance = tol))

          ## Warn (or provide error, if they do not):
          if(any(!rows_eq1)) warning(paste0("Probabilities in row ", which(!rows_eq1), " do not sum up to 1."))
        }

      # retranspose table if it has been transposed:
        if(transpose) p_tab <- t(p_tab)

      return(p_tab)
    }


  # compr_pcomp(test_p2)  # use an example.



  ## (b) Calculate probs from other probs if necessary: ------------------

      ## Testing stuff:
        # p_tabs$p_row <- compr_pcomp(p_tabs$p_row)  # calculate complements first.
        # p_tabs$p_col <- compr_pcomp(p_tabs$p_col, transpose = TRUE)
        #
        # pr <- p_tabs$p_row
        # pc <- p_tabs$p_col
        #
        # pr; pc
        #
        # pr[1,] <- NA
        #
        # ## Redefine as table 1 and 2 which then can be exchanged:
        # ptb1 <- pr  # t(a$p_col)
        # ptb2 <- t(pc)  # a$p_row
        # ptb1; ptb2

    ## Function:
      comp_ptab <- function(ptb1, ptb2, transp_relf = FALSE,  # should relf be transposed?
                            warn = FALSE) {

       ## A. Checking:

         ## Bring matrices into same form:
           transp <- FALSE  # set marker for whether the table has been transposed.

           if(!all(dim(ptb1) == dim(ptb2))) {
             ptb2 <- t(ptb2)  # transpose, if necessary.
             transp <- TRUE

             if(!all(dim(ptb1 == dim(ptb2)))) stop("Input tables need to have the same dimensions. ")
           }

      ## B. Calculate table:

          ## 1. Try to calculate relf_tab:
           relf_tab1 <- ptb1[1:2, ] * ptb2[3, ]

           rlf <- FALSE

           if(!any(is.na(relf_tab1))){  # stopping condition 1?

             rlf <- TRUE

             # ## Complete the table:
             #   relf_tab <- comp_rftab(relf_tab)
             #   pcol <- t(t(relf_tab[1:2, ]) / relf_tab[3, ])
             #   prow <- relf_tab[, 1:2] / relf_tab[, 3]
             #
             # return(list(rftab = relf_tab, prow = prow, pcol = pcol))
           }

         ## 2. Test prerequisite for possible calculation:

           c1 <- !all(is.na(rowSums(ptb1)))  # not all rows in ptb1 may be NA.

           c2 <- sum(is.na(rowSums(ptb2))) < 2  & !any(is.na(ptb2[3, ]))
           ## less than 2 (1) may be NA, one needs to be an unconditional probability (in row 3).

           ## Note: for optimally specified table these are false!

         if(!(c1 & c2) & !rlf) {

           ## Warn, if not turned off:
           if(warn) warning("These probability tables do not allow for any calculation")

           return(list(rftab = FALSE))   # return FALSE (to be passed on to other functions).


         } else {

         ## 3. Try to calculate any other probability:
           # tb1r1 <- tb2[1:2, 1] * tb1[3, ] / tb2[3, 1]  # tb1 row 1.
           # tb1r2 <- tb2[1:2, 2] * tb1[3, ] / tb2[3, 2]  # tb1 row 2.

           ptb1r12 <- t(ptb2[1:2, 1:2] * ptb1[3, ]) / ptb2[3, 1:2]  # alternative for both rows?
           ptb1r3 <- diag(ptb1[1:2, 1:2] * ptb2[3, 1:2] / ptb2[1:2, 1:2])  # tb1 row 3.

           ## Bind the CPs and UCPs:
           ptb1_temp <- rbind(ptb1r12, ptb1r3, deparse.level = 0)   # set deparse.level to avoid naming.

           ptb1_temp <- comb_tabs(ptb1, ptb1_temp)  # combine new information with old information.

           ## Calculate missing complements:
             if(any(rowSums(is.na(ptb1_temp)) == 1)) {  # if any row is incomplete:
               ptb1_temp <- compr_pcomp(ptb1_temp)  # complete the table with complement (if necessary).
             }

        ## C. Testing: ----------------

          ## 1. Test whether inputs are consistent:
           ## Calculate inverse table (but: what about cases with recursion?):
             relf_tab2 <- t(ptb2[1:2, ] * ptb1_temp[3, ])

           ## TODO: Also calculate other values within the function?

           ## Test, whether tables (without NA are equal):
             onetab_na <- is.na(relf_tab1) | is.na(relf_tab2)
             rlf_eq <- isTRUE(all.equal(relf_tab1[!onetab_na], relf_tab2[!onetab_na]))

           if(!rlf_eq) {  # in case the two tables are unequal:

             ## Test, where the problem is:
             # which(relf_tab1 != relf_tab2, arr.ind = TRUE)

             stop("Specified probabilities are inconsistent!")

             ## Use a handler in the final function to localize the problematic inputs!

           } else {

             ## If tables are equal (outer) and one table was complete:
               if (rlf) {

                 relf_tab <- if (transp_relf) t(comp_rftab(relf_tab1)) else comp_rftab(relf_tab1)

                 pcol <- t(t(relf_tab[1:2, ]) / relf_tab[3, ])
                 prow <- relf_tab[, 1:2] / relf_tab[, 3]

                 return(list(rftab = relf_tab, prow = prow, pcol = pcol))

               } else {  ## If tables are equal (outer) but none was complete:

                 relf_tab <- comb_tabs(relf_tab1, relf_tab2)  # combine tables.

                 ## If only one element is missing:
                 if(sum(is.na(relf_tab)) == 1) {
                   relf_tab[is.na(relf_tab)] <- 1 - sum(relf_tab, na.rm = TRUE)
                 }

                 ## HERE!

                 relf_tab <- comp_rftab(relf_tab)  # complete the table.
                 ptb1 <- relf_tab[, 1:2] / relf_tab[, 3]
                 ptb2 <- t(t(relf_tab[1:2, ]) / relf_tab[3, ])

               }


           ## Go into new recursive loop:
            return(comp_ptab(ptb1, ptb2, transp_relf = transp_relf, warn = warn))  # when to do this?

             ## TODO: Can this recursive case happen without overflow?
            }

         }
      }

#
#     ## Test the function:
#      comp_ptab(p_row, p_col)  # get a relative frequency table.
#         ## TODO: Is that sufficient?
#      comp_ptab(t(pc), pr)  # this does (correctly) not work, as pc cannot be calculated.



  ## (c) Testing complements: ------

     test_pcomp <- function(ptab) {

       test <- all.equal(ptab[, 1], 1 - ptab[, 2])  # pretty economic!

       return(isTRUE(test))
     }




## C. Mixed functions: ------------------------------

  ## 1. Probs from frequencies:------------------------

    p_from_f <- function(ftab) {
      p_row <- ftab[1:3, c(1,2)] / ftab[1:3, 3]
      p_col <- t(t(ftab[c(1,2), 1:3]) / ftab[3, 1:3])
      p_dia <- sum(diag(ftab)[1:2]) / diag(ftab)[3]  # add the accuracy metric.

      return(list(p_row = p_row, p_col = p_col, p_dia = p_dia))
    }


  ## 2. Frequencies from probs:-------

      ## a) Function to complete the relative frequency table:
        comp_rftab <- function(rftab) {
          rftab <- rbind(rftab, colSums(rftab))
          rftab <- cbind(rftab, rowSums(rftab))
          return(rftab)
        }

      ## b) Calculate N from feletive frequencies and any given frequencies:
        calc_N <- function(rftab, ftab) {

          if(all(dim(rftab) < dim(ftab))) rftab <- comp_rftab(rftab)  # check whether matrices are equal.

          N <- as.integer(round(ftab / rftab))
            ## round, to obtain (hopefully) the correct numbers.

            ## TODO: Test whether rounding works in all cases!

          N <- unique(N[!is.na(N)])  # get the (hopefully) unique value.

          if (length(N) > 1) warning("There is something off with the calculation of N!")

          return(N)
        }


      ## c) Calculate frequencies from relative frequencies:
        f_from_rf <- function(rftab, ftab = NULL, N = NULL) {

          N_null <- is.null(N)
          ftab_null <- is.null(ftab)

          if(N_null & ftab_null) stop("Either N or a frequency table must be specified. ")

          if(N_null){N <- calc_N(rftab, ftab)}

          ftab <- round(rftab * N, 0)  # calculate frequency table.
          ## TODO: Monitor!

          if(!is_freq(ftab)) warning("Frequency table does not consist of frequencies!")

          ## Complete the table:
          ## Given no margin sums exist, caculate them:
          if(any(dim(ftab) < c(3, 3))) {
            ftab <- rbind(ftab, colSums(ftab))
            ftab <- cbind(ftab, rowSums(ftab))
          }


          return(ftab)

        }

  ## 3. Function to calculate frequencies from prob and frequency (condprob * freq = ucprob): -----------

        ## TODO: Probably an early function!
        ## Test case (definition below):
        tst <- test_p23  # [1, 3] can be calculated from [1, 2] / [1, 5] = 25/.595!


        tst <- complete_tab

        ## Note: Acts on the whole table!  needs all inputs.

        ## Function:
        f_from_pf <- function(tab) {
          ## Note: For rows only (analogous for columns by transposing).
          ## Sum frequency:
          rsfreq <- tab[1:3, 1:2] / tab[1:3, 4:5]
          ## calculates sum frequency from probability.
          sfreq <- comb_tabs(rsfreq[, 1], rsfreq[, 2])

          ## Cell frequency:
          cfreq <- tab[1:3, 3] * tab[1:3, 4:5]
          ## Calculates cell frequency from probability.
          ## Bind and combine:
          rfreq <- cbind(cfreq, sfreq)

          ## TODO:
          ## ISSUE: Round frequencies / throw an error/warning if they are not integer / just leave it.

          # out <- comb_tabs(tst[1:3, 1:3], rfreq)  # returns frequency part only!  This may be changed.

          ## Error if combining tabs does not work: ------------------
            tab[1:3, 1:3] <- tryCatch({
              ## Execute the function:
              comb_tabs(tab[1:3, 1:3], rfreq)
            },
            ## If a warning occurs, throw an informative error message:
            warning = function(w) {
              warning(w)

              comp <- tab[1:3, 1:3] != rfreq  # get comparison.
              conf <- which(comp & !is.na(comp), arr.ind = TRUE)

              ## TODO (Potential): if semantics are included, return also the cell name.

              ## get conflicting cells:
              cls <- apply(conf, 1, paste, collapse = ",")

              ## Collapse cells of interest:
              msg <- paste0("for cell(s) ",
                            paste0("[", cls, "]", collapse = ", "))

              ## Provide error message:
              stop(paste0("Your input is overspecified.\nProbabilities and frequencies imply different values ",
                          msg))
            })

          ## If everything is fine, return: -----------------------
          return(tab)
        }

        ## Test:
        f_from_pf(tst)



## D. Helper functions: --------------------------------

    ## 1. Function to combine two tables: ----------
        comb_tabs <- function(tab1, tab2) {

          ## Check, whether both tables have same type and dimensions:
          if(!isTRUE(all.equal(dim(tab1), dim(tab2)))) {
            stop("Tables do not have the same dimensions.  Execution halted. ")
          }


          ## Way 1:
            # ifelse(is.na(p_tabs$p_row),  # if an entry the first matrix is NA.
            #        ifelse(is.na(p_row), NA, # and if the entry in the second matrix is NA there as well, code NA
            #               p_row),  # else use the value from the second matrix.
            #        ifelse(is.na(p_row), p_tabs$p_row,  # if the entry in the first matrix is not NA,
            #               # but the second is, use the value of the frst matrix.
            #               # otherwise compare the two values.
            #               p_tabs$p_row
            #        )

          ## Way 2: Alternative strategy:
            ## Get indices for NA-values:
            ixna1 <- is.na(tab1)
            ixna2 <- is.na(tab2)
            ix12 <- (!ixna1 & ixna2) | (ixna1 & !ixna2)  # statement for "either NA in 1 or NA in 2".

            ## Test tables for equality where none is NA:
            tab_eq <- isTRUE(all.equal(tab1[!(ixna1 | ixna2)], tab2[!(ixna1 | ixna2)]))
              ## get whether tables are equal everywhere they are not NA.
              ## If any value is false, they do not match (sufficiently)

          if(!all(tab_eq)){

            ## To allow for general usage, throw a warning in case of inconsistency:
            warning("Provided inputs do not match within tolerance. In conflicting cases the first input is used. ")

            }

          ref_tab <- tab1  # set final matrix to reference matrix.

          tab1[ix12 & ixna1] <- 0  # set NAs where only tab1 is NA to 0.
          tab2[ix12 & ixna2] <- 0  # set NAs where only tab2 is NA to 0.
          ref_tab[ix12] <- tab1[ix12] + tab2[ix12]

          return(ref_tab)
        }


## XX. Test cases: -------------------------------------

    ## 0. Full table: ----------

        complete_freq <- rbind(c(17, 25),
                               c(28, 27)
        )
        ## Sums (frequency table):
        complete_freq <- cbind(complete_freq, rowSums(complete_freq))
        complete_freq <- rbind(complete_freq, colSums(complete_freq))

        rel_freq <- complete_freq / complete_freq[dim(complete_freq)[1], dim(complete_freq)[2]]

        ## Probs:
        p_row <- t(apply(complete_freq, 1, function(X) X/X[3]))[,1:2]
        p_col <- apply(complete_freq, 2, function(X) X/X[3])[1:2, ]

        ## Complete table (freqs and probs):
        complete_tab <- cbind(complete_freq, p_row)  # add the row probabilities.
        complete_tab <- rbind(complete_tab, cbind(p_col, NA, NA))
        complete_tab[4, 4] <- (complete_tab[1, 1] + complete_tab[2, 2]) / complete_tab[3, 3]
        complete_tab

    ## 1. Completing frequencies is sufficient:  ----------------
        test_freq <- complete_freq

        test_freq[, 2] <- NA
        test_freq[2, ] <- NA

    ## 2. Frequencies cannot be completed: ----------------
        test_fp <- complete_tab
        ix <- cbind(c(1, 1, 2, 3, 3), c(1, 3, 2, 1, 3))
        test_fp[ix] <- NA
        test_fp[1:3, 1:3]  # frequency portion.
        test_fp[4:5, 1:3]  # col probs.
        test_fp[1:3, 4:5]  # row probs.

    ## 3. Probabilities need complements: ------------------
        test_p1 <- test_fp
        test_p1[4, ] <- NA
        test_p1[, 4] <- NA

    ## 4. Frequencies cannot be calculated, probabilities need complements, and some are missing: ------
        test_p2 <- test_p1
        ix <- cbind(c(5, 5, 3), c(1, 2, 5))
        test_p2[ix] <- NA

tab <- test_p2

    ## 5. Create probability table (also use for translation function):
        prev <- 0.75; sens <- 0.8; spec = 0.78; N = 120

        ## A table for mapping input:
        map_tab <- matrix(c("hi", "fa", "nyes", "PPV", "FOR",
                            "mi", "cr", "nno", "FDR", "NPV",
                          "cond", "nocond", "N", "prev", "inv_prev",
                          "sens", "fart", "ppod", "acc", NA,
                          "FPR", "spec", "inv_ppod", NA, NA),
                          ncol = 5, nrow = 5
        )

        test_p3 <- matrix(NA, ncol = 5, nrow = 5)

        test_p3[map_tab %in% c("sens", "spec", "N", "prev")] <- c(N, prev, sens, spec)


        tab <- test_p3
        # calc_tab(test_p3)



## D. Full function: --------------------------------

calc_tab <- function(tab) {

  otab <- tab  # save original input.

  ## (A) Preparations: -----
    ## (1) Get and check the dimensions of the input: ---------
      n_row <- nrow(tab)
      n_col <- ncol(tab)

    ## Currently all calculations are implemented for a 2x2 frequency table + sums and probabilities (5x5) only!
      ## Throw an exception, if the table does not match:
        if(!all(c(n_row, n_col) == c(5, 5))) {

          ## Create matrix for display:
            mt1 <- matrix(c(rep("freq", 9), rep("prob", 6)), ncol = 5)
            mt2 <- matrix(c(rep("prob", 7), rep(" NA ", 3)), ncol = 5)
            mt <- rbind(mt1, mt2)

            mt_rnd <- paste0(apply(mt, 1, paste0, collapse = " | "), collapse = "\n")  # rendered table.

            ## Some additional text:
            errtxt <- "Table dimensions do not fit.  I need a 5x5 table of the following form:"

            stop(paste0(errtxt, "\n\n", mt_rnd, "\n\n(Some values may be NA)"))
        }


    ## (2) Check, whether there are any or only NAs: -----

        ## Note, that some cells are NA by default!

          mat_hlp <- matrix(TRUE, ncol = n_col, nrow = n_row)  # helper matrix where to check.
          mat_hlp[5, c(4, 5)] <- FALSE
          mat_hlp[4, 5] <- FALSE
          na_log <- is.na(tab) & mat_hlp  # get matrix of NAs overall.

        ## If there are no NAs in the relevant portion, return the table:
          if (!any(na_log)) {
            return(tab)  # return, if table is already complete.
          }

        ## If no values have been provided:
          if(all(na_log)) {

            ## Create matrix for display:
              mt1 <- matrix(c(rep("freq", 9), rep("prob", 6)), ncol = 5)
              mt2 <- matrix(c(rep("prob", 7), rep(" NA ", 3)), ncol = 5)
              mt <- rbind(mt1, mt2)

              mt_rnd <- paste0(apply(mt, 1, paste0, collapse = " | "), collapse = "\n")  # rendered table.

            ## Some additional text:
              errtxt <- "No values have been provided.  I need a 5x5 table of the following form:"

            ## thow exception:
              stop(paste0(errtxt, "\n\n", mt_rnd, "\n\n(Some values may be NA)"))

          }



  ## (B) Calculate frequencies from frequencies and probabilities: ----------------

          ## Note: This is done first as:
              ## - it may already indicate overspecification if results don't match and
              ## - it then allows to potentially calculate all other frequencies.

    ## Calculate frequencies from provided probabilities * frequencies :
      rtab <- f_from_pf(tab = tab)  # for rows.
      ctab <- t(f_from_pf(tab = t(tab)))  # for columns (tranpose and retranspose).

      ## TODO: New problem: calculation of freqs from probs may require rounding!
      ## How to proceed? Loss of precision resulting in errors possible; relative frequencies as solution?
      ## Intuition: Non-integer frequencies are a kind of misspecification and should cause an error
        ## But: interest in cases with known probabilities...

      ## Test case (frombelow):
        ## tst_smp[1,1] <- NA

      tab <- tryCatch({   #define output table.
        comb_tabs(rtab, ctab)  # combine to table.
      },
      warning = function(w) {
        message(w)

        ## Get conflicting cell(s):
        neql <- rtab != ctab  # test inequality.
        conf <- which(neql & !is.na(neql), arr.ind = TRUE)

        ## TODO (Potential): if semantics are included, return also the cell name.

        ## get conflicting cells:
        cls <- apply(conf, 1, paste, collapse = ",")

        ## Collapse cells of interest:
        msg <- paste0("for cell(s) ",
                      paste0("[", cls, "]", collapse = ", "))

        ## Provide error message:
        stop(paste0("Your input is overspecified.\nRow and column probabilities imply different values ",
                    msg))
      }
      )

      ## Here also errors due to inconsistency can occur!
      ## Test case: tst_smp

  ## (C) Calculate frequencies from frequencies: --------

      ftab_us <- tab[1:3, 1:3]  # get the frequency proportion of the table.

      ftab <- comp_ftab(ftab_us)  # calculate the frequency table.

      ## Note: Even without any frequencies a decent table may be provided by calculating an N!

      ## Check output:
      ## CONTINUE HERE!
      ## TODO: Here the function stumbles across probabilities and one frequency only!
      ## This is due to NA values not being equal to the sums.  Use comb_tabs and tryCatch?
        # comb_tabs(rowSums(ftab[, 1:2]), ftab[, 3])
        # comb_tabs(colSums(ftab[1:2, ]), ftab[3, ])

      rsum <- rowSums(ftab[, 1:2])
      csum <- colSums(ftab[1:2, ])  # columns.

      ## Check only non-NA cells (where any sums are):
      rtst <- isTRUE(all.equal(rsum[!is.na(rsum)], ftab[, 3][!is.na(rsum)]))  # rows.
      ctst <- isTRUE(all.equal(csum[!is.na(csum)], ftab[3, ][!is.na(csum)]))  # columns

      ## Catch any inconsistent sums in the table:
      if(any(!c(rsum, csum) & !is.na(c(rsum, csum)))) {

        ## Your inputs in row / col imply a different frequency for cell...
        ## But: already calculated frequencies from probs... test both?

        ## Test user input frequencies:
          conftab <- comp_ftab(otab)

          rsum <- which(rowSums(conftab[, 1:2]) != conftab[, 3])
          csum <- which(colSums(conftab[1:2, ]) != conftab[3, ])

          ## TODO: Test for NA!

          conf_row <- ifelse(length(rsum) > 0,
                             paste0("for row ", rsum, collapse = ", "),
                             "")
          conf_col <- ifelse(length(csum) > 0,
                             paste0("for column ", csum, collapse = ", "),
                             "")

          if(nchar(conf_row) > 0 | nchar(conf_col) > 0) {  # if either applies:

            sep <- ifelse(nchar(conf_row) > 0 && nchar(conf_col) > 0, " and ", "")  # set separator.

            stop("Your frequency input is inconsistent.  It implies different frequencies ",
                 paste(conf_row, conf_col, sep = sep), ".")
          }

        ## Test frequencies already derived from probabilities otherwise:
          ## Test case:
            ## tst_smp[1,1] <- NA; tst_smp[2,1] <- 10

          ## I am not sure, whether this can even occur!  I din't manage to produce this error so far!
          stop("TEST! The probabilities and frequencies you specified imply different frequencies.  Please check your inputs.")

      }



  ## (D) Calculate probabilities from frequencies: -------

      p_tabs <- p_from_f(ftab)  # get probability information from the frequency table.
        ## Will output all NAs, if no frequncies have been provided.

    ## Assemble with probability table provided by the user:
      p_row_us <- tab[1:3, 4:5]  # get row probability table.
      p_col_us <- tab[4:5, 1:3]  # get column probabiilty table.

    ## Combine the matrices and capture inconsistencies:
      ## For row probabilities:
      p_row <- tryCatch({
          ## Execute the function:
          comb_tabs(p_tabs$p_row, p_row_us)
        },
        ## If a warning occurs, throw an informative error message:
        warning = function(w) {
          message(w)
          stop("Provided frequencies imply probabilities different from provided probabilities. ")
        })

      ## For col probabilities:
      p_col <- tryCatch({
          comb_tabs(p_tabs$p_col,  p_col_us)
        },
        warning = function(w) {
          message(w)
          stop("Provided frequencies imply probabilities different from provided probabilities. ")
        })

  ## (D) Calculate probabilities from probabilities: -----

    ## (1) Calculate all possible probability complements: ------
        p_row <- compr_pcomp(p_row)
        p_col <- compr_pcomp(p_col, transpose = TRUE)
          ## transpose to use row function for columns.

        ## Check for consistency of complements; throw error otherwise:
          if(!all(test_pcomp(p_row), test_pcomp(t(p_col)))) {

            ## Gather information to output an informative error message:
              row_num <- which(rowSums(p_row) != 1)  # get row numbers with errors.

              ## Paste rows with complements not addig up to 1:
              if(length(row_num) > 0) {
                row_err <- paste0("Row probabilities in row ", paste0(row_num, collapse = ", "))
              } else {
                row_err <- ""  # if no rows exist, leave empty.
              }


              col_num <- which(colSums(p_col) != 1)  # get col numbers with errors.
              col_err <- ""

              ## Paste columns with complements not addig up to 1:
              if(length(col_num) > 0) {
                if(length(row_num) > 0) {
                  col_err <- paste0(" and column probabilities in column ", paste0(col_num, collapse = ", "))
                } else {  # if no rows are erroneous:
                  col_err <- paste0("Column probabilities in column ", paste0(col_num, collapse = ", "))
                }

            }

            ## Compose error message:
              comp_const <- paste0(row_err, col_err, " do not add up to 1. ")

            ## Throw an error message:
            stop(comp_const)
          }

      ## TODO: HERE!!!###
    ## (2) Calculate missing probabilities from Bayes' theorem: ------
        ## For both directions test for inconsistent probabilities
        ## if the function throws an error, search for the reason to output an informative error message.

      ## (a) Calculation from the direction of the row-table: -----
          relfr <- tryCatch({

            ## Try to calculate the probability tables:
              comp_ptab(p_row, p_col)  # calculate relf in one way ...

            },
            ## If an error occurs:
            error = function (e) {

              ## Start some testing:
                ## Is the problem in user probability input?
                  ## Get user inputs:
                  us_r <- compr_pcomp(p_row_us)
                  us_c <- t(compr_pcomp(t(p_col_us)))
                  us <- NULL

                  ## Test whether user inputs are okay:
                  us <- tryCatch({
                    comp_ptab(us_r, us_c)  # try function on user inputs.
                  },

                  ## If an error occurs, this is evidence that user inputs already were flawed:
                  error = function(e) {

                    # message(e)  # output the original error.
                    stop("Probably: Your specified probabilities are inconsistent and the table is overspecified.
                         Please enter only two probabilities in each dimension,
                         containing at least one unconditional probability (prev, ppod, ...;
                         in case of redundant probabilities you may enter more)")

                  ## TODO: Make this handler scenario dependent by specifying some identifying return value!

                  }
                  )  # end inner tryCatch().

                  if(length(us) > 0) {  # as soon as there is any output to us:
                    ## If the try catch on user input does not produce an error due to inconsistency,
                    ## it has to be dependent on the frequencies specified.

                    # message(e)  # output the original error.
                    stop("Probabilities derived from your frequencies do not match the specified probabilities.
                         Please specify whether you want the table to be based on frequencies or probabilities
                         by deleting one of them. ")
                  }


            })  # end testing relfr.

      ## (b) Calculation from the direction of the column-table: -----
          relfc <- tryCatch({

            ## Try to calculate the probability tables:
            comp_ptab(t(p_col), p_row, transp_relf = TRUE)  # calculate relf in one way ...

            },
            ## If an error occurs:
            error = function (e) {

              ## Start some testing:
              ## Is the problem in user probability input?
              ## Get user inputs:
              us_r <- compr_pcomp(p_row_us)
              us_c <- t(compr_pcomp(t(p_col_us)))
              us <- NULL  # set flag to NULL.

              ## Test whether user inputs are okay:
              us <- tryCatch({
                comp_ptab(t(us_c), us_r, transp_relf = TRUE)  # try function on user inputs.
              },

              ## If an error occurs, this is evidence that user inputs already were flawed:
              error = function(e) {

                # message(e)  # output the original error.
                stop("Probably: Your specified probabilities are inconsistent and the table is overspecified.
                     Please enter only two probabilities in each dimension,
                     containing at least one unconditional probability (prev, ppod, ...;
                     in case of redundant probabilities you may enter more)")

                ## TODO: Make this handler scenario dependent by specifying some identifying return value!

              }
                )  # end inner tryCatch() for relfc.


              ## If the try catch on user input does not produce an error due to inconsistency,
              ## it has to be dependent on the frequencies specified.

              # message(e)  # output the original error.
              stop("Probabilities derived from your frequencies do not match the specified probabilities.
                   Please specify whether you want the table to be based on frequencies or probabilities
                   by deleting one of them. ")

          })  # end testing relfc.

        ## TODO: Do relfr and relfc always produce the same result, 1 no result, 2 no result, or an error?
        ## Or are inconsistent outputs still possible?:
        ## Apparently, you can get no error in one but an error in the other.  So always both should be checked.

        ## TODO:
        ## Apparently, some value combinations are not possible as well!
            ## E.g., for ppod = .3; spec = .52; the NPV cannot be .5 (min ~.54)
            ## TODO: Do these problems always coincide with negative values? Probably yes!


      ## (c) Testing the results: ------------------------------------
        ## TODO!
        ## Test results for consistency:
            ## Potential problems:
              ## No table can be calculated
              ## negative values
              ##

        ## Note: is.logical() does not work for testing, as it encounters a list both times!

        ## (i) Test whether the table can be calculated from at least one direction: ------
            ## For testing use output length!
            lenr <- length(relfr) == 1
            lenc <- length(relfc) == 1

            if (all(lenr, lenc)) {

              ## TODO: Is it possible to calculate some (relative) frequencies and then continue?
                ## Related to below: may allow to calculate N!

              ## Therefore, comp_ptab should return everything it can calculate, in every case!

              ## TODO: Here one may also intersect probabilities and frequencies!
                ## See test_p23: calculate frequency from corresponding probability.
                ## Traces back to obervation that 1 frequenccy and 1 probability per row is sufficient.

              ## For sum frequencies:
                ftab[, 1:2] / p_row

                ftab[1:2, ] / p_col

              ## For cell frequencies:
                ftab[, 3] * p_row
                t(ftab[3, ] * t(p_col))

                ## Currently HERE !!! ##

              # stop("The information you specified is not sufficient to produce the table.
              #      I need at leas one more input. ")

              ## Alternatively:
                warning("Cannot calculate the complete table from given frequencies and probabilities.\nI return everything possible."
                )

                tab[1:3, 1:3] <- ftab
                tab[1:3, 4:5] <- p_row
                tab[4:5, 1:3] <- p_col
                tab[4, 4] <- p_tabs$p_dia

                return(tab)
            }


              ## If at least one table can be calculated, continue:

        ## (ii) Given that both tables could be calculated test and combine them: --------
            if (!any(lenr, lenc)) {

              ## Compare them:
              if(isTRUE(all.equal(relfr, relfc))) {

                relf <- relfr  # combine them if equal.

              } else {  # if both tables were calculated but are inconsistent:

                ## Issue an error message!
                stop("Inputs are inconsistent! ")  ## TODO: Make error message more informative!

              }

            } else {  # if one could not be calculated:

              relf <- list(relfr, relfc)[!c(lenr, lenc)][[1]]  # get calculated table.

            }

        ## (iii) Test for negative values: ---------
        neg_log <- any(unlist(relf) < 0 | unlist(relf) > 1)  # logical index.

          ## Negative values occur if values are specified which are not possible according to the functions
            ## e.g, prev (uc) of 0.3, spec (cc) of 0.49, and an NPV of 0.52 (cc; min for all sens 0.54).

        if (neg_log) {  # if negative values occur:

          ## Identify where negative values stem from:


          stop("Value error: Your inputs produced negative values.
               This indicates that you entered an impossible case. Please revise your inputs. ")

          ## TODO: make the messge more informative!

        }

    ## (E) Mixed calculations -------------------
        ## Calculate frequencies from relf:

        ## HERE!

        # relf <- comp_rftab(relf)  ## complete the table.

        N <- tab[3, 3]  # get provided N.

        if(is.na(N)) {  # if no N was provided

          if (all(is.na(otab[1:3, 1:3]))) {  # test, whether any frequency was provided, if not:

            ## Note, that this case seems very unlikely...

            ## Testing table:
            # [,1]      [,2]      [,3]
            # [1,] 0.1752577 0.2577320 0.4329897
            # [2,] 0.2886598 0.2783505 0.5670103
            # [3,] 0.4639175 0.5360825 1.0000000

            ## 97 is optimal (all integer frequencies)!
            ## One might test what is an optimal factor...

            # riskyr:::factors_min_diff

            prec <- 5
            N <- 10^0  # set N.
            rftb <- round(relf$rftab, prec) * N  # multiply by 1.
            while(!all(rftb > 1)) {
              N <- N * 10
              rftb <- rftb * N
            }

            # rftb

            round(rftb * N, 1)

            ## TODO: Find a nice method to calculate N (with controlled rounding, as well).

          }

          ## If other frequencies are available:
            ftab2 <- f_from_rf(relf$rftab, ftab)
            ## TODO: Also check against N?

            N <- calc_N(relf$rftab, ftab)

        } else {  # if an N is avaliable use it:

            ftab2 <- f_from_rf(relf$rftab, N = N)
        }

          ## TODO: Allow for some testing, whether N is identified (via calc_N?)
          ## TODO: Test somewhere that summing up works (it does not necessarily! weird rounding?)
              ## Case: tab <- test_p2; N <- 100
              ## This may indicate a general scaling problem, especially for smaller numbers (1000 is fine)...


        ## Check provided frequencies for equality:
          ftab_eq <- ftab == ftab2

          if(!all(ftab_eq[!is.na(ftab_eq)])) {
           stop("Frequencies calculated from probabilities do not match provided frequencies. ")

            ## TODO: Problem occurs here!
              ## If frequency inputs do not match probability inputs this results in infinite recursion on ftab!
              ## In case tst_smp one may adjust hits, N, or the probabilities.
              ## Maybe one may simply throw an error...


          } else {
            ftab <- ftab2
          }


    ## (F) General consistency checks:

        ## TODO: Finish!
        ## (1) Test whether probabilities are complements:
          test_pcomp(p_row)
          test_pcomp(t(p_col))


    ## (G) Finishing the table:----

        ## Reassemble tables:
          tab[1:3, 1:3] <- ftab
          tab[1:3, 4:5] <- relf$prow
          tab[4:5, 1:3] <- relf$pcol
          tab[4, 4] <- sum(diag(relf$rftab)[1:2])  # calculate accuracy from relative frequencies!

          tabrf <- tab
          tabrf[1:3, 1:3] <- relf$rftab  # also output table of relative frequencies.

        ## Update NA values:
          na_log <- is.na(tab) & mat_hlp

        ## Test, whether table is complete:
          if(any(na_log)) {  # If there are NAs,

            ## TODO: Condition for stopping:

            ## else: go into a recursive next trial.
            return(calc_tab(tab))

          } else {

          ## Return result:
            return(list(tab, tabrf))
          }

}  # eof calc_tab.

## Testing:
        calc_tab(test_p2)
        all.equal(complete_tab, calc_tab(test_p2)[[1]])

        test_p22 <- test_p2
        test_p22[1, ] <- NA  # remove first line.

        system.time({calc_tab(test_p22)})  # still works.

        test_p22[5, 3] <- 0.7  # changing 1- prev.
        calc_tab(test_p22)  # does not work properly anymore for 0.7 (smaller values are okay).
          ## Problem: values depend on each other, but only one is changed (update function?)
          ## Some givens still seem to change...

        ## Simplified world:
        tst_smp <- test_p3  # use a minimal set.
        calc_tab(tst_smp)  # not overspecified!
        # here everything works fine, as no dependent inputs were specified.
        tst_smp[4, 3] <- 0.01

        calc_tab(tst_smp)

        tst_smp[4, 1]  <- 0.7  # add PPV; should lead to overspecification.
        calc_tab(tst_smp)  # error message due to overspecifiaction.

        tst_smp[4, 1]  <- NA
        tst_smp[1, 1]  <- 30
        calc_tab(tst_smp)  # including another frequency: error for inconsistency.

        test_p23 <- test_p2
        test_p23[, 3] <- NA
        calc_tab(test_p23)  # this doesn't work; more informative error message?
        ## TODO: This would work! In row 1 one could calculate col 3 from prob and (rel)freq!
        ## [1,3] * 0.595 = 25 <=> [1, 3] = 25 / 0.595 = 42

        test_p23[3, 3] <- 32
        calc_tab(test_p23)  # note: N of 32 is stupid and should lead to a warning
         ## (as not in line with other frequencies)

        test_p23[3, 3] <- 97
        calc_tab(test_p23)  # problem: if I change the N, the rest does not follow suit!
        ## The inputs become inconsistent (maybe beyond repair?)
        ## Note, tables cannot be altered in this fashion!
        ## TODO: Manage the warning!  Output looks decent!


        ## The calculations are pretty fast!
        system.time(
          replicate(100, calc_tab(test_p1))
        )

## TODO:
        ## Check the function thoroughly! Where may be problems?
        ## Include stopping condition; is recursiveness in the wrapper even warranted?
        ## What about inconsistent input (non-matching probabilities)?
        ## What to hold constant, if something is changed (frequencies, NPV, sens ...)?


## E. Default tables and translation

      ## TODO: How to interface named methods with the table calculation?
        ## Directly or via general function?


## F. Conceptual stuff: --------------------
    ## 1. Thoughts on the calculation of probabilities: =====


      ## a. Patterns: -----
        ## Rules for being able to calculate a given probability: -----
        ## 1. 3 probabilities: 2 conditional and 1 converse unconditional probability are provided
        ## (then one can also calculate the relative frequency table):
        rw <- matrix(c("a", "b", NA, "1 - a", "1 - b", NA), ncol = 2)
        cl <- matrix(c(NA, NA, "c", NA, NA, "1 - c"), nrow = 2, byrow = TRUE)
        rw; cl
        rw <- !is.na(rw); cl <- !is.na(cl)
        rw; t(cl)  # this pattern works.
        ## --> transposing one and combining them provides a complete 2 x 3 table.

        ## 2. 4 probabilities: Transposing one table either:
        ## a) Same pattern in both matrices; exactly one NA:
        rw <- cbind(c(T, F, T), c(F, F, F))
        cl <- rbind(c(T, F, T), c(F, F, F))
        rw; cl
        ## OR:
        rw <- cbind(c(F, F, F), c(F, T, T))
        cl <- rbind(c(F, F, F), c(F, T, T))
        rw; cl
        rw == t(cl)

        ## b)
        rw <- cbind(c(F, F, F), c(T, F, T))
        cl <- rbind(c(F, T, T), c(F, F, F))
        rw; cl
        ## OR:
        t(rw); t(cl)
        ## Given that one is transposed:
        rw; t(cl)
        ## OR:
        t(rw); cl

    ## b. Other stuff: ----------------------
        ## 1. General procedure for calculating unconditional probabilities:

        ## Test, which are missing:
        ucr_na <- all(is.na(p_tabs$p_row[3, ]))  # are all entries in the last row (unconditional) NA?
        ucc_na <- all(is.na(p_tabs$p_col[, 3]))  # are all entries in the last col (unconditional) NA?

        if (ucc_na | ucr_na & !(ucc_na & ucr_na)) {  # only if some unconditional probabilities are missing.

          ## Calculate sub-matrices:
          a <- pc[1:2,1:2]
          b <- pr[1:2,1:2]
          ## Note, there may not be more than two rows/cols NA!

          ## In case the row probs (typically prev) are NA
          if (ucr_na) {
            ucr <- b * pc[,3] / a  # unconditional row probabilities.
            pr[3,] <- unique(ucr)[rowSums(!is.na(ucr)) > 0, ]
          }

          ## In case the col probs are NA (typically ppod):
          if (ucc_na){
            ucc <- a * pr[3, ] / b  # unconditional column probabilities;
            pc[, 3] <- unique(ucc)[colSums(!is.na(ucr)) > 0, ]  # TODO!  Check colsums!
            # this is the opposite operation (not working here).
          }
        }

        pr <- compr_pcomp(pr)

        ## 2. Approach via relative frequency table:
        relf_tab <- pr[1:2, ] * t(pc)[3, ]  # relative frequencies from c_row & uc_col.

        rlf_tb <- t(t(pc)[1:2, ] * t(pr)[, 3])

        ## NOTE: Condition is that pr[1:2,] or pc[, 1:2] does not contain NAs.

        comb_tabs(relf_tab, rlf_tb)


        ## The example with pr[1,] <- NA amounts to a case b) (see below).
        pr; pc

        ## Redefine as table 1 and 2 which then can be exchanged:
        tb1 <- t(pc)  # t(a$p_col)
        tb2 <- pr  # a$p_row
        tb1; tb2

        which(is.na(tb1), arr.ind = TRUE)  # test only relevant indices?

        ## Prerequisites:
        ## 1 UCP and 2 CPs (either both in other table or 1 in same other in other)
        ## To be able to calculate anything from a given (transposed) table: max 1 row NA.

        ## Calculate conditional probabilities: -----
        ## CP1 = CP2 * UCP2 / UCP1:
        tb1[1, 1] = tb2[1, 1] * tb1[3, 1] / tb2[3, 1]  ## OR:
        tb1[1, 2] = tb2[2, 1] * tb1[3, 2] / tb2[3, 1]

        ## tb1 and tb2 can simply be exchanged to calculate the respective other probability.

        tb1[2, 1] = tb2[1, 2] * tb1[3, 1] / tb2[3, 2]
        tb1[2, 2] = tb2[2, 2] * tb1[3, 2] / tb2[3, 2]

        ## Calculate UCPs: -----
        ## UCP1 = CP1 * UCP2 / CP2:
        tb1[3, 1] = tb1[1, 1] * tb2[3, 1] / tb2[1, 1]
        tb1[3, 1] = tb1[2, 1] * tb2[3, 2] / tb2[1, 2]

        tb1[3, 2] = tb1[1, 2] * tb2[3, 1] / tb2[2, 1]
        tb1[3, 2] = tb1[2, 2] * tb2[3, 2] / tb2[2, 2]

        ## In matrix calculations this now amounts to:
        tb2[1:2, 1] * tb1[3, ] / tb2[3, 1]  # tb1 row 1.
        tb2[1:2, 2] * tb1[3, ] / tb2[3, 2]  # tb1 row 2.

        tb1[1:2, 1] * tb2[3, 1:2] / tb2[1, 1:2]  # tb1 row 3, col 1.
        tb1[1:2, 2] * tb2[3, 1:2] / tb2[2, 1:2]  # tb1 row 3, col 2.

        ## These are 4 instead of 8 equations.

        diag(tb1[1:2, 1:2] * tb2[3, 1:2] / tb2[1:2, 1:2])  # tb1 row 3.

        ## One more equation may be saved.

        ## Final result:
        tb2[1:2, 1] * tb1[3, ] / tb2[3, 1]  # tb1 row 1.
        tb2[1:2, 2] * tb1[3, ] / tb2[3, 2]  # tb1 row 2.
        diag(tb1[1:2, 1:2] * tb2[3, 1:2] / tb2[1:2, 1:2])  # tb1 row 3.


        tpc <- t(pc)  # create a transposed version.
        pr; tpc

        rw <- cbind(c(F, F, F), c(T, F, T))
        cl <- rbind(c(F, T, T), c(F, F, F))
        tcl <- t(cl)
        rw; tcl  # either these may be used or they may be swapped.

        ## 1 - PPV  = ((1 - spec) / DP:) * CN:
        pr[1, ] <- tpc[tcl][1] / tpc[tcl][2] * pr[rw]  # get upper row in pr!

        ## 1 - sens = ((1 - NPV) / CP:) * DN:
        pr[tcl][1] / pr[tcl][2] * tpc[rw]  # get first column in pc!

        ## TODO: Does this hold for all cases?
        ## TODO: Do all possible calculations?

        ## OR:
        t(rw); t(cl)
        ## 1 - sens = ((1 - NPV) / CP:) * DN
        pr[t(cl)][1] / pr[t(cl)][2] * pc[t(rw)]

        sum(is.na(pr[rw]), is.na(pc[cl]))  == 1  # this is the condition to calculate the NA.

        ## Alternatively:
        pc[cl][1] / pc[cl][2] * pr[t(cl)]  # get upper row!


        ## Get the complement(s):
        pr <- compr_pcomp(pr)


        ## From frequencies?:
        # comb_tabs(ftab / 97, relf_tab)

        ## Note, that in the relative frequency table the margin sums
        ## are the corresponding unconditional probabilities!
