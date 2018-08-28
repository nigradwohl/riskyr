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

## D. Full function: --------------------------------

calc_tab <- function(tab) {

}
