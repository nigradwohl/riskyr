## TEMPORARY FILE FOR DEVELOPING GENERALIZATION FUNCTIONS:


## A. Calculating frequencies: --------------------------

## TODO: Function does not calculate sums currently!

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


  ## Set tolereance:
  tol <- 0.001

  ## OLD STUFF:
  ## (a) Mute the frequency proportion:
  ptab <- exp_tab

  ## Set the dimensions of the frequency part:
  frow <- 3
  fcol <- 3
  ## This may be done dynamically in the future.
  ## TODO!
  ptab[1:frow, 1:fcol] <- NA  # mute the frequency part.

    ## Bayes rule: P(A|B) = (P(B|A) * P(A)) / P(B)

  ## (a) Calculate all complements: ------
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

      p_tabs$p_row <- compr_pcomp(p_tabs$p_row)  # calculate complements first.
      p_tabs$p_col <- compr_pcomp(p_tabs$p_col, transpose = TRUE)

      ## Calculate PPV as example:
        ## assume that sens and spec are in the rows.
        ## if unconditional probabilites are missing, it doesn't work.

        ## Abbreviations for testing:
        pr <- p_tabs$p_row
        pc <- p_tabs$p_col

        pr; pc

        pr[1,] <- NA

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


      ## TODO!!! Important note: if one can calculate the relative frequency table, this is sufficient !!! ##
        ## TODO: But is it necessary?

        ## TODO: Finish this!
        ## Plan:
          ## Calculate relftab first
          ## If not possible, calculate unconditional probability
          ## Calculate complements
          ## retry relftab

  ## Function to calculate probabilities: ---------

    ## Testing stuff:
            pr <- p_tabs$p_row
            pc <- p_tabs$p_col

            pr; pc

            pr[1,] <- NA

            ## Redefine as table 1 and 2 which then can be exchanged:
            ptb1 <- pr  # t(a$p_col)
            ptb2 <- t(pc)  # a$p_row
            ptb1; ptb2

     comp_ptab <- function(ptb1, ptb2) {

       transp <- FALSE  # set marker for whether the table has been transposed.

       if(!all(dim(ptb1 == dim(ptb2)))) {
         ptb2 <- t(ptb2)  # transpose, if necessary.
         transp <- TRUE

         if(!all(dim(ptb1 == dim(ptb2)))) stop("Input tables need to have the same dimensions. ")
       }

       ## 1. Try to calculate relf_tab:
          relf_tab <- ptb1[1:2, ] * ptb2[3, ]

       ## 2. Try to calculate any other probability:
          # tb1r1 <- tb2[1:2, 1] * tb1[3, ] / tb2[3, 1]  # tb1 row 1.
          # tb1r2 <- tb2[1:2, 2] * tb1[3, ] / tb2[3, 2]  # tb1 row 2.

          ptb1r12 <- t(ptb2[1:2, 1:2] * ptb1[3, ]) / ptb2[3, 1:2]  # alternative for both rows?
          ptb1r3 <- diag(ptb1[1:2, 1:2] * ptb2[3, 1:2] / ptb2[1:2, 1:2])  # tb1 row 3.

          ## Bind the CPs and UCPs:
          ptb1_temp <- rbind(ptb1r12, ptb1r3, deparse.level = 0)   # set deparse.level to avoid naming.

          ptb1 <- comb_tabs(ptb1, ptb1_temp)  # combine new information with old information.

          ## TODO: Condition:
          ptb1 <- compr_pcomp(ptb1)  # complete the table with complement (if necessary).

     }

        ## Rules for being able to calculate a given probability:
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


    ## Note: ----
        ## Testing complements:
        all.equal(p_tabs$p_row[, 1], 1 - p_tabs$p_row[, 2])  # pretty economic!



## C. Mixed functions: ------------------------------

  ## 1. Probs from frequencies:------------------------

  p_from_f <- function(ftab) {
    p_row <- ftab[1:3, c(1,2)] / ftab[1:3, 3]
    p_col <- t(t(ftab[c(1,2), 1:3]) / ftab[3, 1:3])
    p_dia <- sum(diag(ftab)[1:2]) / diag(ftab)[3]  # add the accuracy metric.

    return(list(p_row = p_row, p_col = p_col, p_dia = p_dia))
  }




        ## TODO!
  ## 2. Frequencies from probs:-------


        ## Is one frequency sufficient to calculate all others?
        complete_tab[1,1]/(complete_tab[1,1]/97)  # yes, by this logic!
        ## If the relative frequency table can be completely calculated, this is a piece of cake.

        ## Knowing then N, one can multiply the relative frequencies to obtain the final table.

        ## TODO: Function taking relative frequencies!

        p_tabs

        pr <- p_tabs$p_row
        pc <- p_tabs$p_col

        relf_tab <- pr[1:2, ] * pc[, 3]
        relf_tab <- rbind(relf_tab, colSums(relf_tab))
        relf_tab <- cbind(relf_tab, rowSums(relf_tab))
        relf_tab

        N <- as.integer(round(ftab / relf_tab))  # round, to obtain (hopefully) the correct numbers.
        N <- unique(N[!is.na(N)])

        if (length(N) > 1) warning("There is something off with the calculation of N!")

        relf_tab * N  # this is correct, but not integer...

        as.integer(relf_tab * N)  # changing to integers comes with some problems!

        a <- round(relf_tab * N, 0)  # rounding seems to do the job (in this case)...
        ## TODO: Monitor!

        p_from_f(a)$p_row == p_tabs$p_row  # test calculating the probabilities from there.

        ftab <- relf_tab * N  # calculate the frequency table.

        is_freq(round(relf_tab * N, 0))

        ## Then do some checking (rowsums, colsums...)


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
          tab_eq <- (tab1 == tab2)[!(ixna1 | ixna2)]
          # get whether tables are equal everywhere they are not NA.
          ## If any value is false, they do not match (sufficiently)

          if(!all(tab_eq)){warning("Provided inputs do not match within tolerance.")}

          ## Set NA values to -1; adding two NAs will produce -2 (to be set NA);
          ## all non-NA values will remain above and +1 can be added.
          ref_tab <- tab1  # set final matrix to reference matrix.

          tab1[ix12 & ixna1] <- 0  # set NAs where b is fine to 0.
          tab2[ix12 & ixna2] <- 0
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

## D. Full function: --------------------------------

calc_tab <- function(tab) {

  ## OLD STUFF: ----------
    ## Get the dimensions:
    n_row <- nrow(tab)
    n_col <- ncol(tab)

    ## (0) Check for any NAs:
    mat_hlp <- matrix(TRUE, ncol = n_col, nrow = n_row)
    mat_hlp[5, c(4, 5)] <- FALSE
    mat_hlp[4, 5] <- FALSE
    na_log <- is.na(tab) & mat_hlp  # get matrix of NAs overall.

    if (!any(na_log)) {
      return(tab)  # return, if table is already complete.
    }

    ## (A) Calculate frequencies from frequencies: -----

      ftab <- tab[1:3, 1:3]  # get the frequency proportion.
      ## Get the dimensions:
      n_rowf <- nrow(ftab)
      n_colf <- ncol(ftab)

      ftab <- comp_ftab(ftab)  # calculate the frequency table.

    ## (B) Calculate probabilities from frequencies: -------

      ## Calculate probability table from frequencies.
      ## TODO: Check first, whether necessary?
        p_tabs <- p_from_f(ftab)

      ## Assemble with probability table provided by the user:
        p_row <- tab[1:3, 4:5]  # get row probability table.
        p_col <- tab[4:5, 1:3]  # get column probabiilty table.

      ## Combine the matrices:
        p_tabs$p_row <- comb_tabs(p_tabs$p_row, p_row)  # reference table needs to be in the first place.
        p_tabs$p_col <- comb_tabs(p_tabs$p_col, p_col)


    ## (C) Calculate probabilities from probabilities: -----
      ## (1) Calculate all possible probability complements (see above?): ------
        p_tabs$p_row <- compr_pcomp(p_tabs$p_row)
        p_tabs$p_col <- compr_pcomp(p_tabs$p_col, transpose = TRUE)
          ## transpose to use row function for columns.

      ## (2) Calculate missing probabilities from Bayes' theorem: ------
        ## TODO: HERE!

    ## (D) Mixed calculations? -------------------
        ## TODO: Even perform them (if necessary)?
    ## Note: As long, as there are two conditional probabilities
    ## and the corresponding unconditional probability, all can be calculated in a 2x2 table.
    ## However, this is redundant with the pathway over the frequencies:
    ## - as soon, as any frequency is given--may be an arbitrary N--a row or column can be calculated;
    ## This allows to calculate the col or row sums
    ## From these one can calculate the cell frequencies via the corresponding conditional ps!

    ## ANY frequency + appropriate probabilities will do!
    ## More generally: one needs one frequency and all (i.e., n_p - 1) probabilities.

    ## This accounts for using the frequency + probability algorithm and the frequency algorithm.


    ##+++HERE+++##

    ## (2) Calculating frequencies from probabilities: ------
      ## Note: as long as there is one probability and one frequency in any row or column it can be calculated.

      ## Mapping:
      ## Within each row or column the other dimension maps 1 -> 4, 2 -> 5 in a 2x2 table.
      ## More generally it maps according to cell_ix -> cell_ix + dim + 1  (dim - 1 + 2).

      ## TODO: When to set N = 1 (or other arbitrary number)?

    ## (D) General consistency checks:
      ## TODO? Calculate the table in different ways and compare the results?

      ## Potential ways of checkign:
      ## - is_valid_prob_pair {riskyr} (passing tol!)
      ## - is_prob

      ## TODO: Compare non-NA values of input and corresponding output values.
        ## Give precedence to table from frequencies!

    ## (E) Finishing the table:

        ## Reassemble tables:
        tab[!is.na(ptab)] <- ptab[!is.na(ptab)]

}
