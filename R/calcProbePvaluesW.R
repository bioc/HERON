calcProbePValuesW<-function(c_mat, colData_in, sd_shift,
                            abs_shift, paired) {

    if (paired) {
        stop("paired with wilcox not supported")
    } else {
        if (!is.na(sd_shift)) {
            stop("sd_shift with wilcox not supported")
        }
        calcProbePValuesWUnpaired(c_mat, colData_in, abs_shift = abs_shift)
    }
}

#' Calculate Probe p-values using a two-sample wilcoxon test
#'
#' @param probe_mat numeric matrix or data.frame of values
#' @param colData_in design data.frame
#' @param exact a logical indicating whether an exact p-value should be
#' computed (see wilcox.test for details)
#' @param abs_shift absolute shift to use when calculating p-values
#'
#' @return matrix of p-values on the post columns defined in the colData matrix
#' @export
#'
#' @importFrom stats wilcox.test
#'
#' @examples
#' data(heffron2021_wuhan)
#' colData_wu <- colData(heffron2021_wuhan)
#' pval_res <- calcProbePValuesWUnpaired(assay(heffron2021_wuhan), colData_wu)
calcProbePValuesWUnpaired<-function(
    probe_mat,
    colData_in,
    exact = NULL,
    abs_shift = 0
) {
    message("exact:", exact)
    if (is.na(abs_shift)) {abs_shift <- 0}
    pre_cols <- rownames(colData_in)[tolower(colData_in$visit) == "pre"]
    post_cols <- rownames(colData_in)[tolower(colData_in$visit) == "post"]
    ans <- matrix(1, nrow=nrow(probe_mat), ncol=length(post_cols))
    colnames(ans) <- post_cols
    rownames(ans) <- rownames(probe_mat)
    for (probe in rownames(ans)) {
        x_val <- as.numeric(probe_mat[probe,pre_cols])
        cexact <- exact # Prevent a warning if there are ties
        if (is.null(cexact) && any(duplicated(x_val))) {
            cexact <- FALSE
        }
        for (post_col in post_cols) {
            y_val <- as.numeric(probe_mat[probe, post_col])

            test_res <- wilcox.test(
                x = x_val,
                mu = y_val - abs_shift,
                paired = FALSE,
                exact = cexact,
                alternative = "less"
            )
            ans[probe, post_col] <- test_res$p.value
        }
    }
    ans <- as.data.frame(ans, stringsAsFactors=FALSE)
    return(ans)
}
