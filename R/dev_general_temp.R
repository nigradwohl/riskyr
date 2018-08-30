## TEMPORARY FILE FOR DEVELOPING GENERALIZATION FUNCTIONS:


## A. Calculating frequencies: --------------------------

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
  comp_pcomp <- function(p_tabs, tol = 0.001) {

    ## Checking for completeness: -------
    ## Version with table input:
      mp_r <- ptab[1:frow, -(1:fcol)]  # matrix of probabilities calculated from rows.
      mp_r_na_log <- is.na(mp_r)
      mp_c <- ptab[-(1:frow), 1:fcol]  # matrix of probabiilities calculated from columns.
      mp_c_na_log <- is.na(mp_c)



    # if (!any(c(pr_inc, pc_inc))) # TODO: this is a stopping condition for the probabilitiy case.

    ## Version with list input:
      mp_r <- p_tabs$p_row
      mp_r_na_log <- is.na(mp_r)
      mp_c <- p_tabs$p_col
      mp_c_na_log <- is.na(mp_c)

      pr_inc <- any(mp_r_na_log)  # check rows.
      pc_inc <- any(mp_c_na_log)  # check cols.

    ## i. Calculate complemets for rows: -----
      if(pr_inc) {  # if row elements are missing.

        row_na <- which(mp_r_na_log, arr.ind = TRUE)  # get rows with NA-values.
        row_na <- rbind(
          row_na[!(duplicated(row_na[,1]) | duplicated(row_na[,1], fromLast = TRUE)), ]
          )  # remove duplicates.
        ## Complete row matrix:
        mp_r[row_na] <- 1 - rowSums(rbind(mp_r[row_na[, 1], ]), na.rm = TRUE)

        ## Check summing up to one:
        rows_eq1 <- isTRUE(all.equal(
          rowSums(rbind(mp_r[row_na[,1], ])),
          rep(1, nrow(row_na)),
          tolerance = tol))

        ## Warn (or provide error, if they do not):
        if(any(!rows_eq1)) warning(paste0("Probabilities in row ", which(!rows_eq1), " do not sum up to 1."))

        ## Include into original table:
        p_tabs$p_row <- mp_r
      }


    ## ii. Calculate complemets for columns: ------
      if (pc_inc) {  # if column elements are missing.

        col_na <- which(mp_c_na_log, arr.ind = TRUE)  # get cols with NA-values.
        col_na <- rbind(
          col_na[!(duplicated(col_na[,2]) | duplicated(col_na[,2], fromLast = TRUE)), ]
          ) # remove duplicates.
        ## Complete column matrix:
        mp_c[col_na] <- 1 - colSums(cbind(mp_c[, col_na[, 2]]), na.rm = TRUE)

        ## Check summing up to one:
        cols_eq1 <- isTRUE(all.equal(
          colSums(cbind(mp_c[, col_na[,2]])),
          rep(1, nrow(col_na)),
          tolerance = tol))

        ## Warn (or provide error, if they do not):
        if(any(!cols_eq1)) {
          warning(paste0("Probabilities in column ", which(!cols_eq1), " do not sum up to 1.\n"))
        }

        ## Include into original table:
        p_tabs$p_col <- mp_c
      }

    return(p_tabs)
  }

  comp_pcomp(test_p2)  # use an example.



  ## (b) Calculate probs from other probs if necessary: ------------------
    ## Understand the calculations in matrix style:

      ## Calculate PPV as example:
        b <- t(t(p_tabs$p_col[, 1:2]) * p_tabs$p_row[3, ])
        b

        b[1, 1] / (b[1, 1] + b[1, 2])  ## PPV.
        b[1, 2] / (b[1, 1] + b[1, 2])  ## 1 - PPV.

        b[2, 1] / (b[2, 1] + b[2, 2])  ## 1- NPV.
        b[2, 2] / (b[2, 1] + b[2, 2])  ## NPV.

      ## Repeat for sens:
        d <- t(p_tabs$p_row[1:2, ] * p_tabs$p_col[, 3])  # just use a slightly different table!
        d
        d == b

        d[1, 1] / (d[1, 1] + d[1, 2])  ## sens.
        d[1, 2] / (d[1, 1] + d[1, 2])  ## 1 - sens.

        d[2, 1] / (d[2, 1] + d[2, 2])  ## 1- spec.
        d[2, 2] / (d[2, 1] + d[2, 2])  ## spec.

        ## Note that you can typically only calculate one of both!

      ## Calculate prevalence:
        p_tabs
        p_tabs$p_row[1,1] * p_tabs$p_col[1, 3] / (p_tabs$p_col[1,1])  # prevalence!
        d[1, 1] / p_tabs$p_col[1,1]  # TP/N / sens.

## C. Mixed functions: ------------------------------

  ## 1. Probs from frequencies:------------------------

  p_from_f <- function(ftab) {
    p_row <- ftab[1:3, c(1,2)] / ftab[1:3, 3]
    p_col <- t(t(ftab[c(1,2), 1:3]) / ftab[3, 1:3])
    p_dia <- sum(diag(ftab)[1:2]) / diag(ftab)[3]  # add the accuracy metric.

    return(list(p_row = p_row, p_col = p_col, p_dia = p_dia))
  }


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
                               c(28, 30)
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
      ## TODO: Check first, whwether necessary?
        p_tabs <- p_from_f(ftab)

      ## Assemble with probability table provided by the user:
        p_row <- tab[1:3, 4:5]  # get row probability table.
        p_col <- tab[4:5, 1:3]  # get column probabiilty table.

      ## Combine the matrices:
        p_tabs$p_row <- comb_tabs(p_tabs$p_row, p_row)  # reference table needs to be in the first place.
        p_tabs$p_col <- comb_tabs(p_tabs$p_col, p_col)





    ## (C) Calculate probabilities from probabilities: -----
      ## (1) Calculate all possible probabilities (see above?): ------
        ## (a) Mute the frequency proportion:
        ptab <- tab

        ## Set the dimensions of the frequency part:
        frow <- 3
        fcol <- 3
        ## This may be done dynamically in the future.
        ## TODO!
        ptab[1:frow, 1:fcol] <- NA  # mute the frequency part.

    ## TODO: Proper naming!!!

    ## (C) Mixed calculations? -------------------
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
