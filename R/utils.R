setConfiguration_ = function(N_MARKERS, LD, simple_interface,
  sigma, N_COV, ALLELEFRQ, N_STRANDS, N_ALLELES = N_MARKERS * 2,
  N_BLOCKS, BLOCK_COR, EFFECT_SIZE = .01, N_CAUSAL = 5, N_CAUSAL_PER_BLOCK = 1,
  PHENO_SD = 3, PHENO_DIST = "gaussian", COR_NOISE_VAR = 0){
  list(N_MARKERS = N_MARKERS, N_ALLELES = N_ALLELES, simple_interface = simple_interface,
    LD = LD, sigma = sigma, N_COV = N_COV,
    ALLELEFRQ = ALLELEFRQ, N_STRANDS = N_STRANDS, N_BLOCKS = N_BLOCKS,
    BLOCK_COR = BLOCK_COR, EFFECT_SIZE = EFFECT_SIZE,
    N_CAUSAL = N_CAUSAL, N_CAUSAL_PER_BLOCK = N_CAUSAL_PER_BLOCK,
    PHENO_SD = PHENO_SD, PHENO_DIST = PHENO_DIST, COR_NOISE_VAR = COR_NOISE_VAR)
}

setConfiguration = function(LD, N_BLOCKS, N_MARKERS) {
  stopifnot(LD %in% c("low", "medium", "high"))
  LD = switch(LD,
              "low" = .3,
              "medium" = .5,
              "high" = .7)
  N_COV = 1
  sigma = matrix(LD, nrow = N_MARKERS, ncol = N_MARKERS)
  diag(sigma) = 1
  N_BLOCKS = 1
  ALLELEFRQ = runif(N_MARKERS, .05, .95)
  N_STRANDS = 2000
  BLOCK_COR = .15
  N_CAUSAL = 1
  EFFECT_SIZE = .01
  PHENO_DIST = "gaussian"
  COR_NOISE_VAR = .02
  config = setConfiguration(N_MARKERS = N_MARKERS,
                            LD = LD, sigma = sigma, N_COV = N_COV, ALLELEFRQ = ALLELEFRQ,
                            N_STRANDS = N_STRANDS,
                            N_BLOCKS = N_BLOCKS,
                            BLOCK_COR = BLOCK_COR, N_CAUSAL = N_CAUSAL,
                            PHENO_DIST = PHENO_DIST,
                            EFFECT_SIZE = EFFECT_SIZE, COR_NOISE_VAR = COR_NOISE_VAR)

}

simCov = function(config, limitZ = 1){
  N_MARKERS = config[["N_MARKERS"]]
  N_COV = config[["N_COV"]]
  sigma = config[["sigma"]]
  simple_interface = config[["simple_interface"]]
  if(simple_interface) {
    return(sigma)
  } else {
      z = Inf
      while(z > limitZ){
        sim = rWishart(N_COV, N_MARKERS, sigma)[, , 1]
        z = psych::cortest(cor(sim), sigma, n1 = N_MARKERS ^ 2)$z
      }
      sim = as.matrix(matrix::nearPD(sim)$mat)
      return(sim)
    }
}


simBlockR = function(block1, block2){
  res = cancor(block1, block2)$cor
  return(max(res))
}

simBlock = function(config){
  N_STRANDS = config[["N_STRANDS"]]
  sigma = config[["sigma"]]
  ALLELEFRQ = config[["ALLELEFRQ"]]
  N_MARKERS = config[["N_MARKERS"]]
  chrom1 = bindata::rmvbin(N_STRANDS, margprob=ALLELEFRQ, sigma=sigma)
  chrom2 = bindata::rmvbin(N_STRANDS, margprob=ALLELEFRQ, sigma=sigma)
  block = chrom1 + chrom2
  return(block)
}

simBlockSet = function(config){
  N_BLOCKS = config[["N_BLOCKS"]]
  BLOCK_COR = config[["BLOCK_COR"]]
  blockList = list()
  for(i in 1:N_BLOCKS){
    if(!length(blockList)){
      blockCur = simBlock(config)
    } else {
      blockNew = simBlock(config)
      while(simBlockR(blockCur, blockNew) > BLOCK_COR){
        blockNew = simBlock(config)
      }
      blockCur = blockNew
    }
    blockList[[length(blockList) + 1]] = blockCur
  }
  blockSet = Reduce(cbind, blockList)
  NCOL = ncol(blockSet)
  colnames(blockSet) = sprintf("SNP_%d", 1:NCOL)
  class(blockSet) = c("genosim", class(blockSet))
  return(blockSet)
}

getBlock = function(gene, config, block){
  N_MARKERS = config[["N_MARKERS"]]
  indices = (((block - 1) * N_MARKERS) + 1):(block * N_MARKERS)
  if(!all(indices %in% seq_len(ncol(gene)))){
    stop("Indices are invalid; block number may be misspecified")
  }
  return(gene[, indices])
}

simGWAS = function(gene, pheno, config, parallel = F, cl){
  PHENO_DIST = config[["PHENO_DIST"]]
  df = data.frame(gene, pheno = pheno$phenotype)
  if(!parallel) {
    res = lapply(names(df)[-ncol(df)], function(x) snpTest(x, df, PHENO_DIST))
  } else {
    stopifnot(!missing(cl) | !inherits(cl, "cluster"))
    parallel::clusterExport(cl, "snpTest")
    res = parallel::parLapply(cl, names(df)[-ncol(df)], function(x) snpTest(x, df, PHENO_DIST))
  }
  res = Reduce(rbind, res)
  rownames(res) = 1:nrow(res)
  res$snp = as.character(res$snp)
  res$block = blockFromSnp(res$snp, config)
  class(res) = c("gwasresult", class(res))
  return(res)
}

snpTest = function(snpName, dat, dist){
  SNP = dat[, snpName]
  pheno = dat[, "pheno"]
  if(dist == "gaussian"){
    lmod = summary(lm(pheno ~ SNP))
    beta = lmod$coefficients[, "Estimate"][2]
    se = lmod$coefficients[, "Std. Error"][2]
    pvalue = lmod$coefficients[, "Pr(>|t|)"][2]
  } else if(dist == "binomial"){
    lmod = summary(glm(pheno ~ SNP, family = binomial(link = "logit")))
    # choosing to not exponentiate
    beta = lmod$coefficients[, "Estimate"][2]
    se = lmod$coefficients[, "Std. Error"][2]
    pvalue = lmod$coefficients[, "Pr(>|z|)"][2]
  } else {
    stop("Phenotype distribution needs to be gaussian or binomial")
  }
  return(data.frame(snp = snpName, beta = beta, se = se, pvalue = pvalue))
}

plot.genosim = function(x) {
  if(requireNamespace("LDheatmap", quietly = TRUE)) {
    rgb.palette = colorRampPalette(rev(c("blue", "red")), space = "rgb")
    LDheatmap::LDheatmap(cor(x),  flip = T, color=rgb.palette(18))
  }
}

plot.gwasresult = function(x, type = "manhattan"){
  if(type == "manhattan"){
    ggplot2::ggplot(data = data.frame(x, number = 1:nrow(x)),
      ggplot2::aes(x = number, y = -log(pvalue), color = as.factor(block))) +
      ggplot2::geom_point(size = 4) + ggthemes::theme_tufte() +
      ggplot2::geom_hline(mapping = ggplot2::aes(yintercept = 3), linetype = "dashed")+
      ggplot2::xlab("SNP number in gene") +
      ggplot2::ylab(expression(-log[10](italic(p)))) +
      ggplot2::labs(colour = "block") +
      ggplot2::ggtitle("Manhattan plot of simulated data")
  } else if(type == "qqplot") {
    print(qq(x$pvalue))
  } else {
    stop("type needs to be manhattan or qqplot")
  }
}

qq = function(pvector, title="Quantile-quantile plot of p-values") {
  # http://gettinggeneticsdone.blogspot.com/2009/11/qq-plots-of-p-values-in-r-using-ggplot2.html
  # Thanks to Daniel Shriner at NHGRI for providing this code for creating expected and observed values
  o = -log10(sort(pvector,decreasing=F))
  e = -log10( 1:length(o)/length(o) )
  plot=ggplot2::qplot(e,o, xlim=c(0,max(e)), ylim=c(0,max(o))) + ggplot2::stat_abline(intercept=0,slope=1, col="red")
  plot=plot+ggplot2::ggtitle(title) + ggthemes::theme_tufte()
  plot=plot+ggplot2::xlab(expression(Expected~~-log[10](italic(p))))
  plot=plot+ggplot2::ylab(expression(Observed~~-log[10](italic(p))))
}

blockFromSnp = function(snpname, config){
  N_BLOCKS = config[["N_BLOCKS"]]
  N_MARKERS = config[["N_MARKERS"]]
  if(is.numeric(snpname)) {
    snpnum = snpname
  } else if(is.character(snpname)) {
    snpnum = as.numeric(substr(snpname, 5, 7))
  } else {
    stop("snpname is not numeric or character")
  }
  block = ceiling(snpnum / N_MARKERS)
  if(block > N_BLOCKS){
    stop("invalid snpname; this block does not exist")
  } else if(block < 1){
    stop("invalid snpname; this block does not exist")
  }
  return(block)
}

blockFromSnp = Vectorize(blockFromSnp, vectorize.args = "snpname")

addCorNoise = function(mat, config) {
  N_ROWS = nrow(mat)
  N_COLS = ncol(mat)
  if(N_ROWS != N_COLS) {
    stop("matrix is not square")
  }
  COR_NOISE_VAR = config[["COR_NOISE_VAR"]]
  mu = rep(0, N_COLS)
  sigma = diag(N_COLS) * COR_NOISE_VAR
  noise = MASS::mvrnorm(N_COLS, mu, sigma)
  noise[lower.tri(noise)] = t(noise)[lower.tri(noise)]
  return(mat + noise)
}
