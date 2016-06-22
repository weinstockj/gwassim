#' Applies the sidaks pvalue correction
#' @export
sidaks = function(gwas){
  result = min(gwas$pvalue)
  result = 1 - (1 - result) ^ nrow(gwas)
  return(result)
}

fisherHelper = function(x){
  k = nrow(x)
  chi = -2 * sum(log(x$pvalue))
  pval = pchisq(chi, df = 2 * k, lower.tail = F)
  return(pval)
}

#' Applies Fisher's pvalue correction
#'
#' @export
fisher = function(gwas){
  result = fisherHelper(gwas)
  return(result)
}

#' Applies Manuel Ferreira's canoncical correlation test
#' @export
ccaTest = function(gene, pheno){
  result = as.numeric(geneCCA(gene, pheno$phenotype))
  return(result)
}

#' @export
vegas = function(gene, gwas){
  N = 1000
  pval = vegasHelper(gene, gwas, N)
  if(pval < 0.1){
    N = 10 ^ 4
    pval = vegasHelper(gene, gwas, N)
#     if(pval < 0.001){
#       N = 10 ^ 6
#       pval = vegasHelper(gene, gwas, N)
#     }
  }
  return(pval)
}

vegasHelper = function(gene, gwas, N) {
  chi = sum(qchisq(gwas$pvalue, 1, lower.tail = FALSE))
  # R = mvtnorm::rmvnorm(N, mean = rep(0, ncol(gene)), sigma = cor(gene), method ='chol')
  R = rmvnorm_(N, mu = rep(0, ncol(gene)), sigma = cor(gene))
  R = R ^ 2
  R = rowSums(R)
  p = sum(R > chi) / N
  return(p)
}

#' @export
gates = function(gene, gwas, config, retKeySnp = FALSE){
  gwas = gwas[order(gwas$pvalue), ] # order by ascending pvalue
  M = gatesHelper(gene, gwas$snp, config)
  N_LIMITS = 35
  if(nrow(gwas) > N_LIMITS){
    gwas = gwas[1:N_LIMITS, ]
  }
  indices = 1:nrow(gwas)
  gwas$m = sapply(indices,
    function(x) gatesHelper(gene, gwas$snp[1:x], config))
  gwas$adjpval = M * gwas$pvalue / gwas$m
  if(!retKeySnp) {
    return(min(gwas$adjpval))
  } else {
    vec = c(min(gwas$adjpval), which.min(gwas$adjpval))
    names(vec) = c("pvalue", "index")
    return(vec)
  }
}

gatesHelper = function(gene, snp_names, config){
  # stopifnot(inherits(gene, "genosim"))
  markers = gene[, snp_names, drop = F]
  mat = addCorNoise(cor(markers), config)
  M = gatesHelper_(mat)
  return(M)
}


#' @export
hyst = function(gene, gwas, config){
  blocks = blockFromSnp(colnames(gene), config)
  res = data.frame(pvalue = numeric(length = length(unique(blocks))),
    key = numeric(length = length(unique(blocks))),
    block = numeric(length = length(unique(blocks))))
  index = 1
  for(b in unique(blocks)) {
    snps = which(blocks == b)
    geneSub = gene[, snps]
    gwasSub = gwas[snps, ]
    res[index, 1:2] = gates(geneSub, gwasSub, config, TRUE)
    res[index, 3] = b
    index = index + 1
  }
  if(nrow(res) > 1) {
    pval = scaleTest(res, gene, config)
  } else {
    pval = res$pvalue
  }
  return(pval)
}

scaleTest = function(res, gene, config) {
  N_BLOCKS = config[["N_BLOCKS"]]
  N_MARKERS = config[["N_MARKERS"]]
  chi = -2 * sum(log(res$pvalue))
  markers = N_MARKERS * (res$block - 1) + res$key
  subLD = cor(gene[, markers])
  scaleVal = 0
  N_COLS = ncol(subLD)
  for(i in 1:N_COLS) {
    for(j in 1:N_COLS) {
      if(i < j){
        rprime = subLD[i, j]
        scaleVal = scaleVal + rprime * (3.25 + 0.75 * rprime)
      }
    }
  }
  scaleVal = 1 + scaleVal / (2 * N_COLS)
  df = (2 * N_COLS) / scaleVal
  pval = pchisq(chi, df = df, ncp = scaleVal, lower.tail = F)
  return(pval)
}

#' @export
skat = function(gene, pheno, config) {
  # Z is matrix of SNPs, in dosage format
  # SKAT fails if class(obj) != "matrix"
  PHENO_DIST = config[["PHENO_DIST"]]
  SKAT_DIST_LABEL = switch(PHENO_DIST,
                           "gaussian" = "C",
                           "binomial" = "D")
  Z = gene
  Z_classes = class(Z)
  NECESSARY_CLASS = "matrix"
  # class(Z) = c(NECESSARY_CLASS, setdiff(Z_classes, NECESSARY_CLASS))
  class(Z) = NECESSARY_CLASS
  obj = SKAT::SKAT_Null_Model(pheno$phenotype ~ 1, out_type = SKAT_DIST_LABEL)
  mod = SKAT::SKAT_CommonRare(Z, obj)
  return(mod$p.value)
}


#' @export
MAGMA = function(phenotype, SNP_matrix, prune = .001){
  pc = princomp(SNP_matrix, cor = FALSE, scores = TRUE)
  pc2 = cumsum(pc$sdev ^ 2 / sum(pc$sdev ^ 2))
  k = length(pc2[pc2 <= 1 - prune])
  fit = lm(phenotype ~ pc$scores[, 1:k])
  f = summary(fit)$fstatistic
  p = pf(f[1], f[2], f[3], lower.tail = F)
  return(p)
}
