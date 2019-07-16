#' Estimate cross-condition mediation probability at the gene-level.
#'
#' For each cis-trans gene pair, estimate the probability of mediation.
#'
#' @param cis_res list object returned by running Primo on gene-level cis-association statistics.
#' @param condcorr_res list object returned by running Primo on conditional correlation statistics.
#' @param cisID vector of gene names corresponding that match the identifiers of of \code{cis_res} (e.g. rows in
#' list element \code{post_prob}).
#' @param condcorrID data.frame with two columns ("cis", "trans") that match the identifiers
#' of \code{condcorr_res} (i.e. rows in list element \code{post_prob}).
#' @param ncond integer of the number of conditions required for mediation.
#'
#' @return A data.frame with the following columns:
#' \tabular{ll}{
#' \code{cis} \tab cis-gene identifier\cr
#' \code{trans} \tab trans-gene identifier\cr
#' \code{ncond} \tab number of conditions required for mediation\cr
#' \code{prob_cis_all} \tab posterior probability of cis-association in all conditions\cr
#' \code{prob_condcorr_ncond} \tab posterior probability of conditional correlation
#' in at least \code{ncond} conditions\cr
#' \code{prob_med} \tab probability of mediation in at least \code{ncond} conditions
#' }
#'
#' The main element of interest for inference is the column \code{prob_med}.
#'
#' @export
#'
CCmed_gene <- function(cis_res,condcorr_res,cisID,condcorrID,ncond){

  ## extract posterior probability matrices
  PP_cis <- cis_res$post_prob
  PP_condcorr <- condcorr_res$post_prob

  ## obtain probabilities used in P_med
  prob_cis_all <- PP_cis[,ncol(PP_cis)]

  prob_condcorr_byN <- Primo::collapse_pp_num(PP_condcorr,prefix="pp_ge_")
  prob_condcorr_ncond <- prob_condcorr_byN[,paste0("pp_ge_",ncond)]

  ## merge by IDs
  prob_cis_all <- data.frame(cis=cisID, prob_cis_all)
  prob_condcorr_ncond <- data.frame(condcorrID, prob_condcorr_ncond)
  PP_merged <- merge(prob_cis_all,prob_condcorr_ncond,by="cis")

  ## calculate probability of mediation
  PP_merged$prob_med <- PP_merged$prob_cis_all * PP_merged$prob_condcorr_ncond
  PP_merged$ncond <- ncond

  ## consider editing so that more than just the final probability is returned
  return(PP_merged[,c("cis","trans","ncond","prob_cis_all","prob_condcorr_ncond","prob_med")])
}

#' Estimate cross-condition mediation probability for a GWAS SNP.
#'
#' For each SNP-cis-trans gene trio, estimate the probability of mediation.
#' Used to identify trans-associations of GWAS SNPs.
#' .
#'
#' @param cis_res list object returned by running Primo on GWAS SNP cis-association statistics.
#' @param condcorr_res list object returned by running Primo on conditional correlation statistics.
#' @param cisID data.frame with two columns ("variant_id", "cis") that match the identifiers of
#' \code{cis_res} (i.e. rows in list element \code{post_prob}).
#' @param condcorrID data.frame with two columns ("cis", "trans") that match the identifiers
#' of \code{condcorr_res} (i.e. rows in list element \code{post_prob}).
#' @param ncond integer of the number of conditions required for mediation.
#'
#' @return data.frame with the following columns:
#' \tabular{ll}{
#' \code{variant_id} \tab GWAS SNP identifier\cr
#' \code{cis} \tab cis-gene identifier\cr
#' \code{trans} \tab trans-gene identifier\cr
#' \code{ncond} \tab number of conditions required for mediation\cr
#' \code{prob_med} \tab probability of mediation in at least \code{ncond} conditions
#' }
#'
#' The main element of interest for inference is the column \code{prob_med}.
#'
#' @export
#'
CCmed_gwas <- function(cis_res,condcorr_res,cisID,condcorrID,ncond){

  ## extract posterior probability matrices
  PP_cis <- cis_res$post_prob
  PP_condcorr <- condcorr_res$post_prob

  ## obtain probabilities used in P_med
  PP_cis_byGroup <- Primo::collapse_pp_combin(PP_cis, combos=ncond, prefix="pp_")
  PP_condcorr_byGroup <- Primo::collapse_pp_combin(PP_condcorr, combos=ncond, prefix="pp_")

  ## merge by IDs
  PP_cis_byGroup <- data.frame(cisID, PP_cis_byGroup)
  PP_condcorr_byGroup <- data.frame(condcorrID, PP_condcorr_byGroup)
  PP_merged <- merge(PP_cis_byGroup,PP_condcorr_byGroup,by="cis")

  ## calculate probability of mediation
  ## consider recoding (second +2 is because of trans gene column location in data.frame)
  PP_med_mat <- as.matrix(PP_merged[,3:ncol(PP_cis_byGroup)]) * as.matrix(PP_merged[,(ncol(PP_cis_byGroup)+2):ncol(PP_merged)])
  PP_merged$prob_med <- matrixStats::rowMaxs(PP_med_mat)
  PP_merged$ncond <- ncond

  ## consider editing so that more than just the final probability is returned
  return(PP_merged[,c("variant_id","cis","trans","ncond","prob_med")])
}
