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

fisher = function(gwas){
  result = fisherHelper(gwas)
  return(result)
}

ccaTest = function(gene, pheno){
  result = as.numeric(geneCCA(gene, pheno$phenotype))
  return(result)
}

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
  R = mvtnorm::rmvnorm(N, mean = rep(0, ncol(gene)), sigma = cor(gene), method ='chol')
  R = R ^ 2
  R = apply(R, 1, sum)
  p = length(R[R > chi]) / N
  return(p)
}

gates = function(gene, gwas, retKeySnp = FALSE){
  gwas = gwas[order(gwas$pvalue), ] # order by ascending pvalue
  M = gatesHelper(gene, gwas$snp)
  N_LIMITS = 35
  if(nrow(gwas) > N_LIMITS){
    gwas = gwas[1:N_LIMITS, ]
  }
  indices = 1:nrow(gwas)
  gwas$m = sapply(indices, function(x) gatesHelper(gene, gwas$snp[1:x]))
  gwas$adjpval = M * gwas$pvalue / gwas$m
  if(!retKeySnp) {
    return(min(gwas$adjpval))
  } else {
    vec = c(min(gwas$adjpval), which.min(gwas$adjpval))
    names(vec) = c("pvalue", "index")
    return(vec)
  }
}

gatesHelper = function(gene, snp_names){
  # stopifnot(inherits(gene, "genosim"))
  markers = gene[, snp_names, drop = F]
  mat = cor(markers)
  changeR = function(x) {
    0.2982 * x ^ 6 - 0.0127 * x ^ 5 + 0.0588 * x ^ 4 +
      0.0099 * x ^ 3 + 0.628 * x ^ 2 - 0.0009 * x
  }
  mat = apply(mat, 1:2, function(x) changeR(x))
  decomp = eigen(mat)
  eigenvals = decomp$values
  eigenvals = eigenvals[eigenvals > 1]
  val = sum(eigenvals - 1)
  M = length(snp_names) - val
  return(M)
}

hyst = function(gene, gwas, config, est_blocks = FALSE){
  if(est_blocks) {
    est_n_blocks = sum(eigen(cor(gene)) > 1)
    pca_decomp = princomp(gene)
    hap_blocks = apply(pca_decomp, 1, which.max) # chooses block to load on to
    gwas$block = hap_blocks
  } else {
    blocks = blockFromSnp(colnames(gene), config)
    res = data.frame(pvalue = numeric(length = length(unique(blocks))),
      key = numeric(length = length(unique(blocks))),
      block = numeric(length = length(unique(blocks))))
    index = 1
    for(b in unique(blocks)) {
      snps = which(blocks == b)
      geneSub = gene[, snps]
      gwasSub = gwas[snps, ]
      res[index, 1:2] = gates(geneSub, gwasSub, TRUE)
      res[index, 3] = b
      index = index + 1
    }
    pval = scaleTest(res, gene, config)
    return(pval)
  }
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
