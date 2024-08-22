calcProbePValuesT<-function(c_mat, colData_in, sd_shift,
                            abs_shift, paired) {

    if (paired) {
        calcProbePValuesTPaired(c_mat, colData_in, sd_shift, abs_shift)
    } else {
        calcProbePValuesTUnpaired(c_mat, colData_in, sd_shift, abs_shift)
    }
}

#' Calculate Probe p-values using a differential unpaired t-test
#'
#' @param probe_mat numeric matrix or data.frame of values
#' @param colData_in design data.frame
#' @param sd_shift standard deviation shift to use when calculating p-values
#' Either sd_shift or abs_shift should be set
#' @param abs_shift absolute shift to use when calculating p-values
#'
#' @return matrix of p-values on the post columns defined in the colData matrix
#' @export
#'
#' @examples
#' data(heffron2021_wuhan)
#' colData_wu <- colData(heffron2021_wuhan)
#' pval_res <- calcProbePValuesTUnpaired(assay(heffron2021_wuhan), colData_wu)
calcProbePValuesTUnpaired<-function(
        probe_mat,
        colData_in,
        sd_shift=NA,
        abs_shift=NA
) {
    if (!is.na(sd_shift) && !is.na(abs_shift)) {
        stop("Either sd or abs can be set, not both.")
    }
    pre_cols <- rownames(colData_in)[tolower(colData_in$visit) =="pre"]
    post_cols <- rownames(colData_in)[tolower(colData_in$visit) == "post"]

    pre_means <- rowMeans(probe_mat[,pre_cols])
    pre_sds <- matrixStats::rowSds(as.matrix(probe_mat[,pre_cols]))
    post_means <- rowMeans(probe_mat[,post_cols])
    post_sds <- matrixStats::rowSds(as.matrix(probe_mat[,post_cols]))
    pre_var <- matrixStats::rowVars(as.matrix(probe_mat[,pre_cols]))
    n <- length(pre_cols)
    pre_stderr <- sqrt(pre_var / n)
    dfree <- n-1

    pars <- data.frame(
        pre_mean = pre_means, pre_sds = pre_sds,
        post_mean = post_means, post_sds = post_sds,
        pre_stderr = pre_stderr, diff_mean = post_means - pre_means,
        dfree = rep(dfree, nrow(probe_mat)),
        stringsAsFactors=FALSE
    )
    rownames(pars) <- rownames(probe_mat)
    post_mat <- probe_mat[,post_cols]
    post_tv <- getPostTVal(post_mat, pre_means, pre_stderr,
        pre_sds, sd_shift, abs_shift)
    colnames(post_tv) <- post_cols

    ans <- matrix(NA, nrow = nrow(probe_mat),ncol=length(post_cols))
    rownames(ans) <- rownames(probe_mat)
    colnames(ans) <- post_cols


    for (post_col in post_cols) {
        ans[,post_col] <-
            stats::pt(q=post_tv[,post_col], df=dfree, lower.tail=FALSE)
    }
    ans <- as.data.frame(ans,stringsAsFactors=FALSE)

    pars <- cbind(pars, post_tv)
    attr(ans, "pars") <- pars
    return(ans)
}



#' Calculate Probe p-values using a differential paired t-test
#'
#' @param probe_mat numeric matrix or data.frame of values
#' @param colData_in design data.frame
#' @param sd_shift standard deviation shift to use when calculating p-values.
#' Either sd_shift or abs_shift should be set
#' @param abs_shift absolute shift to use when calculating p-values.
#' @param debug print debugging information
#'
#' @return matrix of p-values on the post columns defined in the colData matrix.
#' Attributes of the matrix are:
#'
#' pars - data.frame parameters used in the paired t-test for each row
#' (e.g. df, sd)
#'
#' mapping - data.frame of mapping used for pre-post column calculation
#' diff_mat - data.frame
#' containing the post-pre differences for each sample (column) and probe (row)
#'
#' @export
#'
#' @examples
#' data(heffron2021_wuhan)
#' colData_wu <- colData(heffron2021_wuhan)
#' pre_idx = which(colData_wu$visit == "pre")
#' ## Make some samples paired
#' colData_post = colData_wu[colData_wu$visit == "post",]
#' new_ids = rownames(colData_post)[seq_len(5)]
#' colData_wu$ptid[pre_idx[seq_len(5)]] = new_ids
#' exprs <- assay(heffron2021_wuhan, "exprs")
#' pval_res <- calcProbePValuesTPaired(exprs, colData_wu)
calcProbePValuesTPaired <- function(
        probe_mat,
        colData_in,
        sd_shift = NA,
        abs_shift = NA,
        debug = FALSE
) {
    if (!is.na(sd_shift) && !is.na(abs_shift)) {
        stop("Either sd or abs can be set. Not both.")
    }
    mapping <- getPairedMapping(colData_in)
    ans <- matrix(NA, nrow = nrow(probe_mat), ncol=nrow(mapping))
    diff_mat <- probe_mat[,mapping$post] - probe_mat[,mapping$pre]
    colnames(diff_mat) <- mapping$ptid
    rownames(diff_mat) <- rownames(probe_mat)
    pars <- matrix(data = NA, nrow=nrow(probe_mat), ncol = 5)
    colnames(pars) <- c("diff_mean", "diff_sd", "diff_var",
                        "diff_stderr", "dfree")
    rownames(pars) <- rownames(probe_mat)
    for (r_idx in seq_len(nrow(probe_mat))) {
        x <- c(t(probe_mat[r_idx,mapping$post] - probe_mat[r_idx,mapping$pre]))
        nx <- length(x)
        mx <- mean(x, na.rm=TRUE)
        sx <- stats::sd(x, na.rm=TRUE)
        vx <- stats::var(x, na.rm=TRUE)
        stderr <- sqrt(vx/nx)
        dfree <- sum(!is.na(x))
        pars[r_idx, "diff_mean"] <- mx
        pars[r_idx, "diff_var"] <- vx
        pars[r_idx, "diff_sd"] <- sx
        pars[r_idx, "diff_stderr"] <- stderr
        pars[r_idx, "dfree"] <- dfree
        for (c_idx in seq_len((nrow(mapping)))) {
            tp <- getPTP(x[c_idx], stderr, sd_shift, sx, abs_shift, dfree)
            ans[r_idx, c_idx] <- tp
        }
    }
    ans <- as.data.frame(ans, stringsAsFactors=FALSE)
    colnames(ans) <- mapping$ptid
    rownames(ans) <- rownames(probe_mat)
    attr(ans, "pars") <- pars
    attr(ans, "mapping") <- mapping
    attr(ans, "diff_mat") <- diff_mat
    return(ans)
}

getPostTVal <- function(
        post_mat, pre_means,
        pre_stderr, pre_sds,
        sd_shift, abs_shift) {
    no_shift <- is.na(sd_shift) && is.na(abs_shift)
    if (no_shift) {
        post_tvalues <- (post_mat - pre_means)/pre_stderr
    } else if (!is.na(sd_shift)) {
        post_tvalues <- (post_mat - pre_means-sd_shift*pre_sds)/pre_stderr
    } else {
        post_tvalues <- (post_mat - pre_means-abs_shift)/pre_stderr
    }
    return(post_tvalues)
}

