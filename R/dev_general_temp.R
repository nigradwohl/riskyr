## TEMPORARY FILE FOR DEVELOPING GENERALIZATION FUNCTIONS:

## A. Calculating frequencies: --------------------------

## Function to calculate rows: --------------------------
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

      stop("The information you specified is not sufficient to calculate the full table.\nI need at least one more value."
      )

      ## CONTINUE if not complete:
    } else {

      return(comp_ftab(ftab))
    }

  }

## B. Probability table: ------------------------------

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

    ## Checking for completeness:
    mp_r <- ptab[1:frow, -(1:fcol)]  # matrix of probabilities calculated from rows.
    mp_r_na_log <- is.na(mp_r)
    mp_c <- ptab[-(1:frow), 1:fcol]  # matrix of probabiilities calculated from columns.
    mp_c_na_log <- is.na(mp_c)

    pr_inc <- any(mp_r_na_log)  # check rows.
    pc_inc <- any(mp_c_na_log)  # check cols.

    if (!any(c(pr_inc, pc_inc))) # TODO: this is a stopping condition for the probabilitiy case.

    ## Calculate compelemets for rows: -----
      if(pr_inc) {  # if row elements are missing.

        row_na <- which(mp_r_na_log, arr.ind = TRUE)  # get rows with NA-values.
        ## Complete row matrix:
        mp_r[row_na] <- 1 - rowSums(rbind(mp_r[row_na[, 1], ]), na.rm = TRUE)

        ## Check summing up to one:
        rows_eq1 <- isTRUE(all.equal(rowSums(mp_r), rep(1, nrow(mp_r)), tolerance = tol))

        ## Warn (or provide error, if they do not):
        if(any(!rows_eq1)) warning(paste0("Probabilities in row ", which(!rows_eq1), " do not sum up to 1."))

        ## Include into original table:
        ptab[1:frow, -(1:fcol)] <- mp_r
      }


    ## Calculate compelemets for columns: ------
    if (pc_inc) {  # if column elements are missing.

      col_na <- which(mp_c_na_log, arr.ind = TRUE)  # get cols with NA-values.
      ## Complete column matrix:
      mp_c[col_na] <- 1 - colSums(cbind(mp_c[, col_na[, 2]]), na.rm = TRUE)

      ## Check summing up to one:
      cols_eq1 <- colSums(mp_c) == 1

      ## Warn (or provide error, if they do not):
      if(any(!comp_col)) {
        warning(paste0("Probabilities in column ", which(!cols_eq1), " do not sum up to 1.\n"))
      }

      ## Include into original table:
      ptab[-(1:frow), 1:fcol] <- mp_c
    }



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

    ## (B) Calculate probabilities from probabilities: -----
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

    ## (E) Finishing the table:

        ## Reassemble tables:
        tab[!is.na(ptab)] <- ptab[!is.na(ptab)]

}
