#' Calculate Probe-level p-values
#'
#' calculates p-values on the matrix (can be sequence as well), using either a
#' z-score global, t-stat differential, or a combination of both.
#' input a matrix that has samples on the columns and sequence probes on the
#' rows.
#'
#' @param probe_mat numeric matrix, rows as seqs and columns as samples
#' @param colData_in experiment design data.frame
#' @param d_sd_shift sd shift multiplier for diff test
#' @param d_abs_shift abs shift to use for the diff test
#' @param d_paired do paired t-test
#' @param g_sd_shift sd shift multiplier for global test
#' @param use which p-value method(s) to use
#'
#' @return matrix of calculated p-values
#' @noRd
calcProbePValuesProbeMat<-function(
        probe_mat,
        colData_in,
        d_sd_shift = NA,
        d_abs_shift = NA,
        d_paired = FALSE,
        g_sd_shift=0,
        use = "both") {

    diff_fxn <- switch(
        use,
        "z" = NULL,
        "t" = calcProbePValuesT,
        "tz" = calcProbePValuesT,
        "w" = calcProbePValuesW,
        NULL
    )
    g_fxn <- switch(
        use,
        "z" =  calcProbePValuesZ,
        "tz" = calcProbePValuesZ,
        NULL
    )
    if (is.null(diff_fxn) && is.null(g_fxn)) {
        stop("Unknown use parameter:", use)
    }
    if (missing(colData_in) || is.null(colData_in)) {
        c_mat <- probe_mat
    } else {
        c_mat <- probe_mat[,rownames(colData_in)]
    }
    praw <- NULL
    if (!is.null(diff_fxn)) {
        praw <- diff_fxn(c_mat, colData_in, d_sd_shift, d_abs_shift, d_paired)
    }
    if (!is.null(g_fxn)) {
        pg <- g_fxn(c_mat, colData_in, g_sd_shift)
        if (is.null(praw)) {
            praw <- pg
        } else {
            praw <- combinePValueMatrix(praw, pg)
        }
    }
    praw[is.na(praw)] <- 1.0 # Conservatively set NAs to p-value = 1
    ans <- praw
    return(ans)
}

combinePValueMatrix<-function(pmat1, pmat2) {
    use_cols <- intersect(colnames(pmat1), colnames(pmat2))
    if (length(use_cols) != ncol(pmat1)) {
        warning("Combining p-values, some columns are mismatched")
    }
    ans <- pmat1[,use_cols]
    for (col in use_cols) {
        ans[,col] <- pmax(ans[,col], pmat2[,col])
        ##(Wilkinson's max p-value)
        ans[,col] <- stats::pbeta(ans[,col], 2, 1)
    }
    return(ans)
}

#' Calculate p-values using the "exprs" assay
#'
#' @param obj HERONSequenceDataSet or HERONProbeDataSet
#' @param colData_in optional column DataFrame (default: NULL => colData(obj)))
#' @param d_sd_shift standard deviation shift for differential test
#' @param d_abs_shift absolute shift for differential test
#' @param d_paired run paired analysis
#' @param g_sd_shift standard deviation shift for global test
#' @param use use global-test ("z"), differential-test using t.test ("t"),
#' differential-test using wilcox ("w"), or both global and differential ("tz")
#' @param p_adjust_method method for adjusting p-values
#'
#' @return HERONSequenceDataSet/HERONProbeDataSet with the pvalue assay added
#' @export
#'
#' @examples
#' data(heffron2021_wuhan)
#' seq_pval_res <- calcCombPValues(heffron2021_wuhan)
calcCombPValues<-function(
    obj,
    colData_in = NULL,
    d_sd_shift = NA,
    d_abs_shift = NA,
    d_paired = FALSE,
    g_sd_shift = 0,
    use = "tz",
    p_adjust_method = "BH"
) {
    stopifnot(is(obj, "HERONSequenceDataSet") || is(obj, "HERONProbeDataSet"))
    if (is.null(colData_in)) {colData_in <- colData(obj)}
    pval <- calcProbePValuesProbeMat(
        probe_mat = assay(obj, "exprs"),
        colData_in = colData_in,
        d_sd_shift = d_sd_shift,
        d_abs_shift = d_abs_shift,
        d_paired = d_paired,
        g_sd_shift = g_sd_shift,
        use = use
    )

    res <- addPValues(obj, pval)
    res <- p_adjust_ds(res)
    return(res)
}

#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment assay<-
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors metadata
#' @importFrom S4Vectors metadata<-
#' @noRd
addPValues<-function(obj, pval) {
    res <- obj
    if (all(dim(assay(obj, "exprs")) == dim(pval)) &&
        all(colnames(assay(obj,"exprs")) == colnames(pval)) &&
        all(rownames(assay(obj,"exprs")) == rownames(pval))
    ) {
        assay(res, "pvalue") <- pval
    } else {
        # Post rows/columns
        exprs_old <- assay(obj, "exprs")
        colData_old <- colData(obj)
        colData_new <- colData_old[(colnames(pval)),]
        rownames(colData_new) <- colData_new$ptid
        exprs_new <- exprs_old[rownames(pval), colnames(pval)]
        colnames(exprs_new) <- colData_new$ptid
        if (is(obj, "HERONSequenceDataSet")) {
            res <- HERONSequenceDataSet(
                exprs = exprs_new,
                colData = colData_new
            )
            metadata(res) <- metadata(obj)
        } else if (is(obj, "HERONProbeDataSet")) {
            rowRanges_new <- rowRanges(obj)
            res <- HERONProbeDataSet(
                assays = list(exprs = exprs_new),
                colData = colData_new,
                rowRanges = rowRanges_new
            )
        }
        colnames(pval) <- colData_new$ptid
        assay(res, "pvalue") <- pval
    }
    return(res)
}

#' Calculate Global p-values Using Normal (z) Distribution
#'
#' @param probe_mat numeric matrix or data.frame of values
#' @param colData_in design data.frame
#' @param sd_shift standard deviation shift to use when calculating p-values
#' (default 0)
#' @param all calculate p-values for all samples, otherwise report the p-values
#' for just the post samples
#'
#' @return matrix of "p-values"
#'
#' @examples
#' data(heffron2021_wuhan)
#' colData_wu <- colData(heffron2021_wuhan)
#' pval_res <- calcProbePValuesZ(assay(heffron2021_wuhan), colData_wu)
#' @noRd
calcProbePValuesZ<-function(
        probe_mat,
        colData_in,
        sd_shift = 0,
        all = FALSE
) {

    if (all || missing(colData_in) || is.null(colData_in)) {
        message("No colData or all asked for, calculating on all columns")
        all_cols <- colnames(probe_mat)
        post_cols <- all_cols
    } else {
        pre_cols <- rownames(colData_in)[tolower(colData_in$visit)=="pre"]
        post_cols <- rownames(colData_in)[tolower(colData_in$visit) == "post"]
        all_cols <- c(pre_cols, post_cols)
    }

    vals <- unlist(probe_mat[,all_cols])
    g_mean <- mean(vals, na.rm=TRUE)
    g_sd <- stats::sd(vals, na.rm=TRUE)

    zvalues <- matrix(NA, nrow = nrow(probe_mat), ncol=length(post_cols))
    rownames(zvalues) <- rownames(probe_mat)
    colnames(zvalues) <- post_cols
    zvalues[,post_cols] <- as.matrix(probe_mat[,post_cols] - g_mean) / g_sd

    pvalues <- apply(zvalues - sd_shift, 2, stats::pnorm, lower.tail=FALSE)

    pars <- c("mean" = g_mean, "sd" = g_sd)
    attr(pvalues,"pars") <- pars
    attr(pvalues,"zscore") <- zvalues

    return(pvalues)
}


getPairedMapping<-function(colData) {
    pre_df <- colData[tolower(colData$visit) =="pre",, drop=FALSE]
    post_df <- colData[tolower(colData$visit) == "post",, drop=FALSE]

    mapping <- data.frame(
        ptid = pre_df$ptid,
        pre = rownames(pre_df),
        post = rep(NA, nrow(pre_df)),
        stringsAsFactors=FALSE
    )
    rownames(mapping) <- mapping$ptid
    for (idx in seq_len(nrow(post_df))) {
        post_ptid <- post_df$ptid[idx]
        if (post_ptid %in% rownames(mapping)) {
            mapping[post_ptid,"post"] <- rownames(post_df)[idx]
        }
    }

    mapping <- stats::na.omit(mapping)
    return(mapping)
}

getPTP<-function(x, stderr, sd_shift, sx, abs_shift, dfree) {
    no_shift <- is.na(sd_shift) & is.na(abs_shift)
    if (no_shift) {
        tstat <- (x)/stderr
    } else if (!is.na(sd_shift)) {
        tstat <- (x-sd_shift*sx)/stderr
    } else {
        tstat <- (x-abs_shift)/stderr
    }
    ans <- stats::pt(tstat, dfree, lower.tail=FALSE)
    return(ans)
}




#' Calculate protein-level p-values
#'
#' @param epitope_ds HERONEpitopeDataSet with the "pvalue" assay
#' @param metap_method meta p-value method to use
#' @param p_adjust_method p.adjust method to use
#'
#' @details
#' see calcEpitopePValues for a list of meta p-value methods supported
#' by HERON. the protein should be one that requires at least one of the
#' epitope p-values to be small (e.g. wmax1).
#'
#' @return HERONProteinDataSet with the "pvalue" and "padj" assays
#' @export
#'
#' @seealso [stats::p.adjust()] for p_adjust_parameter.
#' @seealso [calcEpitopePValues()] for meta p-value methods
#'
#' @examples
#' data(heffron2021_wuhan)
#' pval_seq_res <- calcCombPValues(heffron2021_wuhan)
#' pval_pr_res <- convertSequenceDSToProbeDS(pval_seq_res)
#' calls_res <- makeProbeCalls(pval_pr_res)
#' segments_res <- findEpitopeSegments(calls_res, "unique")
#' epval_res <- calcEpitopePValues(calls_res, segments_res)
#' ppval_res <- calcProteinPValues(epval_res)
calcProteinPValues<-function(
        epitope_ds,
        metap_method = "wmin1",
        p_adjust_method = "BH"

) {
    stopifnot(is(epitope_ds, "HERONEpitopeDataSet"))
    stopifnot("pvalue" %in% assayNames(epitope_ds))
    epi_pvalues <- assay(epitope_ds, "pvalue")
    protein_pvalues <- calcMetaPValuesMat(
        pvalues_mat = epi_pvalues,
        by_vec = getEpitopeProtein(rownames(epi_pvalues)),
        method = metap_method
    )
    res <- HERONProteinDataSet(pvalue = protein_pvalues)
    colData(res) <- colData(epitope_ds)
    metadata(res) <- metadata(epitope_ds)
    res <- p_adjust_ds(res, p_adjust_method)
    return(res)
}


#' Calculate epitope-level p-values
#'
#' @param probe_pvalues matrix of probe p-values
#' @param epitope_ids vector of epitope ids
#' @param metap_method meta p-value method to use
#'
#' @return matrix of epitope-level p-values
#' @noRd
calcEpitopePValuesMat<-function(
    probe_pvalues,
    epitope_ids,
    metap_method = "wmax1"
) {

    epi_probe_df <- getEpitopeIDsToProbeIDs(epitope_ids)
    #Remove any probes that are not in the probe_pvalues matrix
    epi_probe_df <- epi_probe_df[epi_probe_df$PROBE_ID %in%
                                    rownames(probe_pvalues),]

    #Order the p-values to match the epi_probe_df
    probe_pvalues_mat <- probe_pvalues[epi_probe_df$PROBE_ID,]

    #by vec will be the epitope_id associated with the probe.
    by_vec <- epi_probe_df$Epitope_ID

    #Do the meta p-value calculation.
    epitope_pvalues <- calcMetaPValuesMat(
        pvalues_mat = probe_pvalues_mat,
        by_vec = by_vec,
        method = metap_method
    )
    return(epitope_pvalues)
}

#' Calculate epitope-level p-values
#'
#' @param probe_pds HERONProbeDataSet with the "pvalue" assay
#' @param epitope_ids vector of epitope ids
#' @param metap_method meta p-value method to use (see below)
#' @param p_adjust_method what p.adjust method to use.
#'
#' @details
#' The meta p-value methods supported by \code{calcEpitopePValues} are:
#' min_bonf*,
#' min*,
#' max*,
#' fischer/sumlog,
#' hmp/harmonicmeanp,
#' wilkinsons_min1/tippets,
#' wilkinsons_min2/wmin2,
#' wilkinsons_min3,
#' wilkinsons_min4,
#' wilkinsons_min5,
#' wilkinsons_max1/wmax1,
#' wilkinsons_max2/wmax2,
#' and cct.
#'
#' When choosing a p-value method, keep in mind that the epitope p-value should
#' be one that requires most of the probe p-values to be small (e.g. *wmax1*)
#' Other p-value methods such as the*cct* and the *hmp* have been shown to be
#' more accurate with p-value that have dependencies.
#'
#' @seealso [stats::p.adjust()] for p_adjust_parameter.
#'
#' @return HERONEpitopeDataSet with "pvalue" and "padj" assays
#' @export
#'
#' @examples
#' data(heffron2021_wuhan)
#' pval_seq_res <- calcCombPValues(heffron2021_wuhan)
#' pval_pr_res <- convertSequenceDSToProbeDS(pval_seq_res)
#' calls_res <- makeProbeCalls(pval_pr_res)
#' segments_res <- findEpitopeSegments(calls_res, "unique")
#' epval_res <- calcEpitopePValues(calls_res, segments_res)
#' @importFrom SummarizedExperiment assayNames
#' @importFrom SummarizedExperiment colData<-
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment assays<-
#' @importFrom SummarizedExperiment mcols
calcEpitopePValues<-function(
        probe_pds,
        epitope_ids,
        metap_method = "wmax1",
        p_adjust_method = "BH"
) {
    stopifnot(is(probe_pds, "HERONProbeDataSet"))
    stopifnot("pvalue" %in% assayNames(probe_pds))
    pvalues_mat <- calcEpitopePValuesMat(
        probe_pvalues = assay(probe_pds, "pvalue"),
        epitope_ids = epitope_ids,
        metap_method = metap_method
    )

    eds <- HERONEpitopeDataSet(pvalue = pvalues_mat)
    colData(eds) <- colData(probe_pds)
    metadata(eds) <- metadata(probe_pds)
    metadata(eds)$probe_meta <- mcols(rowRanges(probe_pds))
    metadata(eds)$epitope_ids <- epitope_ids
    eds <- p_adjust_ds(eds, p_adjust_method)
    return(eds)
}

#' Adjust a matrix of p-values column-by-column
#'
#' @param pvalues_mat matrix of p-values
#' @param method what adjustment algorithm to use (see p.adjust)
#'
#' @return matrix of column-by-column adjusted p-values
#'
#' @examples
#' mat = matrix(runif(25), nrow=5)
#' rownames(mat) = paste0("A;",seq_len(nrow(mat)))
#' p_adjust_mat(mat)
#' @noRd
p_adjust_mat<-function(pvalues_mat, method = "BH") {
    ans <- pvalues_mat
    for (col_idx in seq_len(ncol(ans))) {
        ans[,col_idx] <- stats::p.adjust(ans[,col_idx], method=method)
    }
    return(ans)
}

#' Adjust an assay of p-values column-by-column
#'
#' @param obj SummarizedExperiment with a "pvalue" assay
#' @param method what adjustment algorithm to use (see p.adjust)
#'
#' @return SummarizedExperiment with the "padj" assay added
#' @noRd
p_adjust_ds <- function(obj, method = "BH") {
    stopifnot(is(obj, "SummarizedExperiment"))
    stopifnot("pvalue" %in% assayNames(obj))
    assay(obj, "padj") <- p_adjust_mat(assay(obj, "pvalue"), method = method)
    return(obj)
}

