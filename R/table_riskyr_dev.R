## 3. table.riskyr function + helpers: ----------------

## 3a. riskyr_table function: -------------------

## Function takes a riskyr object and transforms it into a matrix: ---------

tst <- riskyr(N = 10, prev = 0.3, sens = 0.8, spec = 0.6)  # riskyr object for testing.

riskyr_table <- function(object, type = "cd") {


  ## Note:  this interface is specific to the diagnostic case (alias: riskyr.diagnostic):

  ## TODO: Decide on flexibility!
  ## Matrix for condition-decision case:
  ## This one should probably not be flexible to have clear variable names for internal use:
  # cnd_dec_var <- matrix(c("hi", "mi", "cond.true", "sens", "FPR",
  #                         "fa", "cr", "cond.false", "FNR", "spec",
  #                         "dec.pos", "dec.neg", "N", "ppod", "pned",  # pned: proportion negative decisions.
  #                         "PPV", "FOR", "prev", "acc", "",
  #                         "FDR", "NPV", "nprev", "", ""),
  #                       ## TODO: Find good variable names for the missing quantities.
  #                       nrow = 5, ncol = 5)
  ## Variable names are commented out for current purposes.

  ## for user-defined display an additional flexible matrix may be used:
  cnd_dec_lab <- with(object,
                      matrix(c("hi", "mi", "cond.true", "Sens", "FNR",
                               "fa", "cr", "cond.false", "FPR", "Spec",
                               "dec.pos", "dec.neg", "N", "ppod", "ppnd",
                               "PPV", "FOR", "Prev", "Acc", "",
                               "FDR", "NPV", "1-prev", "", ""),
                             nrow = 5, ncol = 5)
  )


  # Treatment vs. condition:
  trt_cnd_lab <- with(object,
                      matrix(c("hi", "mi", "treated", "Sens", "FPR",
                               "fa", "cr", "untreated", "FNR", "Spec",
                               "condition", "no condition", "N", "prev", "1 - prev",
                               "PPV", "FOR", "ptreated", "Acc", "",
                               "FDR", "NPV", "puntreated", "", ""),
                             nrow = 5, ncol = 5)
  )

  lab_mat <- switch (type,
                     cd = cnd_dec_lab,
                     tc = trt_cnd_lab)



  ## Basic: dim 1 (cnd) and dim 2 (dec)
  probs <- comp_prob_prob(prev = object$prev, sens = object$sens, spec = object$spec)
  freqs <- comp_freq_prob(N = object$N, prev = object$prev, sens = object$sens, spec = object$spec)


  ## Create a numeric matrix:
  num_mat <- with(object,
                  matrix(c(freqs$hi, freqs$mi, freqs$cond_true, probs$sens, probs$mirt,
                           freqs$fa, freqs$cr, freqs$cond_false, probs$fart, probs$spec,
                           freqs$dec_pos, freqs$dec_neg, freqs$N, probs$ppod, 1-probs$ppod,
                           probs$PPV, probs$FOR, probs$prev, probs$acc, NA,
                           probs$FDR, probs$NPV, 1-probs$prev, NA, NA),
                         nrow = 5, ncol = 5
                  )
  )


  ## Descriptive information:
  desc <- list(
    scen_lbl = object$scen_lbl,
    popu_lbl = object$popu_lbl,
    sdt_lbl = object$sdt_lbl,

    cond_lbl = object$cond_lbl,
    dec_lbl = object$dec_lbl,

    scen_txt = object$scen_txt,
    scen_src = object$scen_src,
    scen_apa = object$scen_apa,
    scen_lng = object$scen_lng
    )

  obj <- list(dsc = desc, lbl = lab_mat, num = num_mat)

  ## TODO: For now limit to 1st and 2nd dimension!

  class(obj) <- c("riskyr.tabular")


  return(obj)
}


# -------------------------------------------

## Plotting of tabular objects: ----------

plot.riskyr.tabular <- function(x = NULL,        # require riskyr scenario
                        type = "prism",  # default type
                        # by = "cddc",   # default perspective
                        ...) {


  type <- tolower(type)  # ensure lowercase

  # Test type argument:
  if (!type %in% c(#
    # plot_prism:
    "prism", "fprism", "tree", "ftree", "net", "fnet", "network",
    # plot_area:
    "area", "farea", "mosaic",
    # plot_tab:
    "tab", "table", "ftab", "ctab",
    # plot_icons:
    "icon", "icons", "iconarray",
    # plot_bar:
    "bar", "bars", "barplot", "fbar",
    # plot_curve:
    "curve", "curves",
    # plot_plane:
    "plane", "planes", "cube")) {

    message("Unknown plot type (in plot.riskyr): Using type = 'prism'.")
    type <- "prism"

  }


## Get central inputs:
  nums <- x$num  # numerical information.
  lbls <- x$lbl  # adjusted label information.

  # TODO: Ultimately, functions should rely on the tabular representation to be more flexible.

  x_txt <- init_txt(scen_lbl = x$dsc$scen_lbl,

                    popu_lbl = x$dsc$popu_lbl,
                    N_lbl = lbls[3, 3],

                    cond_lbl = x$dsc$cond_lbl,
                    cond_true_lbl = lbls[3, 1],
                    cond_false_lbl = lbls[3, 2],

                    dec_lbl  = x$dsc$dec_lbl,
                    dec_pos_lbl = lbls[1, 3],
                    dec_neg_lbl = lbls[2, 3],

                    acc_lbl = lbls[4, 4],
                    dec_cor_lbl = "",  # TODO!
                    dec_err_lbl = "x$dec_err_lbl",

                    sdt_lbl = x$dsc$sdt_lbl,
                    hi_lbl = lbls[1, 1],
                    mi_lbl = lbls[2, 1],
                    fa_lbl = lbls[1, 2],
                    cr_lbl = lbls[2, 2],

                    scen_txt = x$dsc$scen_txt,
                    scen_src = x$dsc$scen_src,
                    scen_apa = x$dsc$scen_apa,
                    scen_lng = x$dsc$scen_lng
  )


  if ((substr(type, 1, 3) == "tab") || (type == "ftab") || (type == "ctab")) {

    plot_tab(prev = nums[3, 4],
             sens = nums[4, 1], mirt = nums[5, 1],
             spec = nums[5, 2], fart = nums[4, 2],
             N = nums[3, 3],
             # Options:
             lbl_txt = x_txt,
             title_lbl = x_txt$scen_lbl,
             ...
    )
  } # if (type == "tab")

  ## 2. Area / mosaic plot:
  if ((substr(type, 1, 4) == "area") || (type == "farea") ||
      (substr(type, 1, 6) == "mosaic")) {  # "mosaic"

    plot_area(prev = nums[3, 4],
              sens = nums[4, 1], mirt = nums[5, 1],
              spec = nums[5, 2], fart = nums[4, 2],
              N = nums[3, 3],
              # Options:
              lbl_txt = x_txt,
              title_lbl = x_txt$scen_lbl,
              ...
    )

  } # if (type == "area")

  ## 3. Icon array:
  if (substr(type, 1, 4) == "icon") {

    plot_icons(prev = nums[3, 4],
               sens = nums[4, 1], mirt = nums[5, 1],
               spec = nums[5, 2], fart = nums[4, 2],
               N = nums[3, 3],
               # Options:
               lbl_txt = x_txt,
               title_lbl = x_txt$scen_lbl,
               ...
    )

  } #  if (type == "icon")

  ## 4. Prism plot:
  if ((substr(type, 1, 5) == "prism") || (substr(type, 1, 6) == "fprism") ||
      (substr(type, 1, 3) == "net")   || (substr(type, 1, 4) == "fnet")   ||
      (substr(type, 1, 4) == "tree")  || (substr(type, 1, 5) == "ftree")) {

    plot_prism(prev = nums[3, 4],
               sens = nums[4, 1], mirt = nums[5, 1],
               spec = nums[5, 2], fart = nums[4, 2],
               N = nums[3, 3],
               # Options:
               lbl_txt = x_txt,
               title_lbl = x_txt$scen_lbl,
               ...
    )

  } # if (type == "prism")

  ## 5. Bar plot / frequency bars:
  if ((substr(type, 1, 3) == "bar") || (substr(type, 1, 4) == "fbar")) {

    plot_bar(prev = nums[3, 4],
             sens = nums[4, 1], mirt = nums[5, 1],
             spec = nums[5, 2], fart = nums[4, 2],
             N = nums[3, 3],
             # Options:
             lbl_txt = x_txt,
             title_lbl = x_txt$scen_lbl,
             ...
    )

  } # if (type == "bar")

  ## 6. Curve of probabilities:
  if (substr(type, 1, 5) == "curve") {

    plot_curve(prev = nums[3, 4],
               sens = nums[4, 1], mirt = nums[5, 1],
               spec = nums[5, 2], fart = nums[4, 2],
               N = nums[3, 3],
               # Options:
               lbl_txt = x_txt,
               title_lbl = x_txt$scen_lbl,
               ...
    )
  } # if (type == "curve")

  ## 7. Plane/cube of probabilities:
  if ((substr(type, 1, 5) == "plane") || (substr(type, 1, 4) == "cube")) {

    plot_plane(prev = nums[3, 4],
               sens = nums[4, 1], mirt = nums[5, 1],
               spec = nums[5, 2], fart = nums[4, 2],
               N = nums[3, 3],
               # Options:
               lbl_txt = x_txt,
               title_lbl = x_txt$scen_lbl,
               ...
    )
  } # if (type == "plane")

}



## Test examples: --------
a <- riskyr_table(tst)
plot(a)

b <- riskyr_table(tst, type = "tc")
plot(b, f_lbl = "nam")
plot(b, type = "tree", by = "cdac", f_lbl = "nam", p_lbl = "nam")
plot(b, type = "icons", lbl_txt = b$lbl)  # TODO: does not work!

# TODO: Here switching perspectives becomes interesting!

## Summarising riskyr objects: -----------------

## 3b. nice_table function: --------------------

nice_tab <- function(smr, space = 2) {  # TODO: Give object class summary.

  max_nchar <- max(nchar(smr))
  padding <- paste0(rep(" ", space), collapse = "")  # repeat for padding.
  cellwd <- max_nchar + space * 2  # cellwidth.

  ## Prolong each cell to max (repeat for each cell and collapse):
  naddspace <- c(max_nchar - nchar(smr))  # get the number of spaces to add.
  addspace <- sapply(naddspace, function(x) {paste0(rep(" ", times = x), collapse = "")})

  pad_vec <- rep(padding, length(smr))
  #pad_vec[1] <- paste0(padding, " ", collapse = "")

  cells <- paste0(pad_vec, smr, addspace, padding)
  mcells <- matrix(cells, nrow = 5, ncol = 5)  # bring back to form.

  ctab <- paste0(apply(mcells, 1, paste0, collapse = "|"), "\n")  # only with colseps.

  csep <- paste0(rep("-", cellwd), collapse  = "")  # horizontal separator for each cell.
  hsep <- paste0(paste0(rep(csep, 5), collapse = "|"), "\n")  # horizontal separator.

  tab <- c("\n", hsep, rbind(ctab, rep(hsep, length(ctab))))

  cat(tab)


  ## TODO: Add headers!
  ## TODO: First row is off; maybe due to the \n separator.

}


## 3c. as.table.riskyr method: --------------------
## Takes a riskyr object and prints out a nice table:

#' Create table of risk information.
#'
#' \code{as.table.riskyr} provides a \code{as.table} method for objects of class "riskyr" to allow for readable output.
#'
#' @format A "riskyr" object.
#'
#' @param x An object of class "riskyr", usually a result of a call to \code{riskyr}.
#'
#' @param sep By which character should the entries be separated from the labels? (default " = ").
#'
#' @examples
#' table(scenarios$n4)
#'
#' @family summary functions
#'
#' @export

as.table.riskyr <- function(x, sep = " = ") {

  tab <- riskyr_table(x)  # convert riskyr to table.

  nums <- round(tab$num, 2)
  lbl <- tab$lbl
  nums[is.na(nums)] <- ""
  out <- matrix(paste0(lbl, sep, nums), nrow = 5, ncol = 5)
  out[out == sep] <- ""  # replace empty cells.

  class(out) <- c("table.riskyr")

  return(out)


  #
  #             ## TODO: Allow the transposed forms as well.
  #
  #             ix <- ifelse(relf, 2, 1)  # decide whether to report relative frequencies.
  #             head <- ifelse(relf, "relf", "freq")
  #             cat(head, "\n")
  #             return(out)

}

## Test the class and method:
# tst <- riskyr(hi = 1, mi = 3, fa = 4, cr = 5)
#

## 3d. print.table.riskyr method: ----------------

#' Printing tabular risk information.
#'
#' \code{print.table.riskyr} provides a \code{print} method for objects of class "table.riskyr".
#'
#' @format Printed output of a "table.riskyr" object.
#'
#' @param x An object of class "table.riskyr", usually a result of a call to \code{table.riskyr}.
#'
#' @param space How much space should be included left and right of the table entries?
#'
#' @examples
#' table(scenarios$n4)
#'
#' @family print functions
#'
#' @export

print.table.riskyr <- function(x, space = 2) {
  nice_tab(x, space = space)
}


