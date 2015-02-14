#' Add together two numbers.
#'
#' @param config A configuration object.
#' @param gene A genosim object.
#' @return An object of class phenosim. Includes simulated phenotype,
#' causal blocks, causal markers, and structure


simPhenotype = function(config, gene, verbose = F){
  N_CAUSAL_PER_BLOCKS = config[["N_CAUSAL_PER_BLOCK"]]
  N_CAUSAL = config[["N_CAUSAL"]]
  N_BLOCKS = config[["N_BLOCKS"]]
  N_MARKERS = config[["N_MARKERS"]]
  EFFECT_SIZE = config[["EFFECT_SIZE"]]
  PHENO_SD = config[["PHENO_SD"]]
  PHENO_DIST = config[["PHENO_DIST"]]
  N_CAUSAL_BLOCKS = floor(N_CAUSAL / N_CAUSAL_PER_BLOCKS)
  # ASSUMING EFFECTS ARE INDEPENDENT
  EFFECT_SIZE_PER_CAUSAL = EFFECT_SIZE / N_CAUSAL
  COR = sqrt(EFFECT_SIZE_PER_CAUSAL)
  chosenMarkers = vector("integer", length = N_CAUSAL)
  slopes = vector("numeric", length = N_CAUSAL)
  chosenBlocks = vector("integer", length = N_CAUSAL_BLOCKS)
  countBlock = 1
  countMarker = 1
  while(countBlock <= N_CAUSAL_BLOCKS){
    block = sample(1:N_BLOCKS, size = 1)
    if(block %in% chosenBlocks) {
      next
    } else {
      chosenBlocks[countBlock] = block
      countBlock = countBlock + 1
      indices = (((block - 1) * N_MARKERS) + 1):(block * N_MARKERS)
      marker = sample(indices, size = N_CAUSAL_PER_BLOCKS, replace = F)
      for(m in marker){
        # choose sign of effect of marker
        sign = sample(c(-1, 1), size = 1)
        COR = COR * sign
        chosenMarkers[countMarker] = m
        markerSd = sd(gene[, m])
        # SLOPE = COR * Sy / Sx
        slope = COR * PHENO_SD / markerSd
        slopes[countMarker] = slope
        countMarker = countMarker + 1
      }
    }
  }
  causalString = paste(sprintf("%.2f * SNP_%d", slopes, chosenMarkers),
    collapse = " + ")
  pheno = eval(parse(text = causalString), envir = as.data.frame(gene))
  phenonew = pheno + rnorm(n = length(pheno), mean = 0,
    sqrt(PHENO_SD ^ 2 - var(pheno)))
  if(EFFECT_SIZE == 0){
    chosenBlocks = NULL
    chosenMarkers = NULL
    causalString = "SNPs do not explain phenotype; no causal string"
  }
  if(verbose){
    cat(sprintf("\nChosen SNPs: %s", paste(chosenMarkers, collapse = ", ")))
    cat(sprintf("\nChosen Blocks: %s", paste(chosenBlocks, collapse = ", ")))
    cat(sprintf("\nCausal structure is: %s", causalString))
  }
 if(PHENO_DIST == "gaussian"){

  } else if(PHENO_DIST == "binomial") {
    # pass through sigmoid
    phenonew = round(1 / (1 + exp(-phenonew)))
  } else {
    stop("Phenotype distribution needs to be gaussian or binomial")
  }
  phenoList = list(phenotype = phenonew, blocks = chosenBlocks,
    markers = chosenMarkers, causalString = causalString)
  class(phenoList) = "phenosim"
  return(phenoList)
}

print.phenosim = function(x){
  cat("\nHere are the properties of the simulated phenotype:\n")
  cat(sprintf("\nChosen SNPs: %s", paste(x$markers, collapse = ", ")))
  cat(sprintf("\nChosen Blocks: %s", paste(x$blocks, collapse = ", ")))
  cat(sprintf("\nCausal structure is: %s", x$causalString))
  cat(sprintf("\nDistribution is: %s", x$dist))
  cat(sprintf("\nmean is: %.4f", mean(x$phenotype)))
  cat(sprintf("\nstandard deviation is: %.4f", sd(x$phenotype)))

}

summary.phenosim = function(x){
  print(x)
}

