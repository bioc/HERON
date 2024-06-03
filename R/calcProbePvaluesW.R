#' Title
#'
#' @param probe_mat
#' @param colData_in
#' @param exact
#' @param abs_shift
#'
#' @return
#' @export
#'
#' @importFrom stats wilcox.test
#'
#' @noRd
calcProbePValuesWUnpaired<-function(
    probe_mat,
    colData_in,
    exact = FALSE,
    abs_shift = 0
) {
    pre_cols <- rownames(colData_in)[tolower(colData_in$visit) == "pre"]
    post_cols <- rownames(colData_in)[tolower(colData_in$visit) == "post"]
    ans <- matrix(1, nrow=nrow(probe_mat), ncol=length(post_cols))
    colnames(ans) <- post_cols
    rownames(ans) <- rownames(probe_mat)
    for (probe in rownames(ans)) {
        x_val <- as.numeric(probe_mat[probe,pre_cols])
        for (post_col in post_cols) {
            y_val <- as.numeric(probe_mat[probe, post_col])

            test_res <- suppressWarnings(
                wilcox.test(
                    x = x_val,
                    mu = y_val - abs_shift,
                    paired = FALSE,
                    exact = exact,
                    alternative = "less"
                )
            )
            ans[probe, post_col] <- test_res$p.value
        }
    }
    ans <- as.data.frame(ans, stringsAsFactors=FALSE)
    return(ans)
}
