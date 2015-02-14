library("gwassim")
library("parallel")

N_MARKERS = 8
LD = sqrt(.9)
N_COV = 1
sigma = matrix(LD, nrow = N_MARKERS, ncol = N_MARKERS)
diag(sigma) = 1
N_BLOCKS = 15
ALLELEFRQ = runif(N_MARKERS, .05, .95)
N_STRANDS = 2000
# Number of individuals * 2 (DNA strands).
BLOCK_COR = .15
N_CAUSAL = 4
EFFECT_SIZE = .0
PHENO_DIST = "gaussian"
config = setConfiguration(N_MARKERS = N_MARKERS,
  LD = LD, sigma = sigma, N_COV = N_COV, ALLELEFRQ = ALLELEFRQ,
  N_STRANDS = N_STRANDS,
  N_BLOCKS = N_BLOCKS,
  BLOCK_COR = BLOCK_COR, N_CAUSAL = N_CAUSAL,
  PHENO_DIST = PHENO_DIST,
  EFFECT_SIZE = EFFECT_SIZE)

totalSim = function(config){
  gene1 = simBlockSet(config)
  pheno = simPhenotype(config, gene1)
  gwas = simGWAS(gene1, pheno, config)
  EFFECT_SIZE = config[["EFFECT_SIZE"]]
  N_MARKERS = config[["N_MARKERS"]]
  N_BLOCKS = config[["N_BLOCKS"]]
  LD = config[["LD"]]
  res = data.frame(EFFECT_SIZE, N_MARKERS, N_BLOCKS, LD)
  res$sidaks = sidaks(gwas)
  res$ccaTest = ccaTest(gene1, pheno)
  res$fisher = fisher(gwas)
  res$vegas = vegas(gene1, gwas)
  res$gates = gates(gene1, gwas)
  return(res)
}

ptm <- proc.time()
N_SIM = 500
cl = makeCluster(3)
clusterExport(cl, varlist = c("totalSim", "config"))
clusterEvalQ(cl, {library(gwassim)})
res = parLapply(cl, 1:N_SIM, function(x) totalSim(config))
stopCluster(cl)
res = Reduce("rbind", res)
proc.time() - ptm
colMeans(res)
lapply(res[, 5:ncol(res)], function(x) prop.table(table(x < .05)))
# 2.26 seconds..per sim

