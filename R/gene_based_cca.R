
geneCCA = function(gene, phen, verbose = FALSE) {

  dna=na.omit(gene)

  nind=nrow(dna)
  nsnps=ncol(dna)

  # This version is for univariate analysis only (ie. single quantitative trait)
  ntraits = 1

  ### ========================================
  ### Remove SNPs that cause multicollinearity
  keep.dna=check_mc(dna)
  dna.pruned=dna[, keep.dna]

  ### ===================
  ### Gene-based analysis
  p.gene = gene.assoc(dna.pruned, phen)

  if(verbose){
    cat(paste("\nN individuals:",nind))
    cat(paste("\nN SNPs before pruning:", ncol(dna)))
    cat(paste("\nN SNPs after pruning:", ncol(dna.pruned)))
    cat(paste("\nCCA gene-based P-value =", p.gene))
  }
  invisible(data.frame(pvalue = p.gene))
}

### ================
### Helper functions
###

### Gene-based test
gene.assoc = function(dna,phen)
{
  nsnps=max(ncol(dna),1,na.rm=T); ntraits=max(ncol(phen),1,na.rm=T);
  if (nsnps == 1)
    nind = length(dna)
  else
    nind=nrow(dna)

  all.can = as.numeric(cancor(dna,phen)$cor)
  all.can[all.can>0.99]=0.99

  # F-test
  fapp = function(nsnps,ntraits,nind,all.can)
  {
    wilks = prod(1-all.can^2)
    w = nind - 1 - 0.5*(nsnps+ntraits+1)
    t = sqrt( (nsnps^2 * ntraits^2 - 4) / (nsnps^2 + ntraits^2 - 5) )
    if (nsnps*ntraits == 2) t=1
    df1 = nsnps * ntraits
    df2 = w*t - 0.5*nsnps*ntraits + 1
    f = ( (1-wilks^(1/t)) / (wilks^(1/t)) ) * df2/df1
    p = pf(f,df1,df2,lower.tail=F)
    p
  }
  p = fapp(nsnps,ntraits,nind,all.can)
  p
}



### Function to CHECK FOR MULTI-COLLINEARITY
### Thanks to Shaun Purcell

check_mc <- function(traits)
{

  maxvif=2	# VIF threshold
  maxcor=0.8	# r2 threshold

  ### Remove SNPs based on high pairwise correlation
#   tnames=1:ncol(traits)
#   rem = numeric(0)
#   for (i1 in 1:ncol(traits))
#     for (i2 in 1:ncol(traits))
#       if (i1<i2)
#         if( cor( traits[,i1],traits[,i2],use="complete.obs" )^2 > maxcor | !is.finite(cor(traits[,i1],traits[,i2],use="complete.obs")) )
#           rem = c(rem,i1)
#   rem = as.integer(which(apply(cor(gene1) ^ 2 > .38 & cor(gene1) < 1, 2, function(x) any(x))))
  tnames = 1:ncol(traits)
#   removeIdx = function(x){
#     vec = cor(traits) ^ 2
#     vec = vec[1:x, x]
#     if(any(vec > maxcor & vec < 1)){
#       return(x)
#     } else return()
#   }
#   removeIdx = Vectorize(removeIdx, vectorize.args = "x")
#   rem = unlist(removeIdx(tnames))
  rem = unique(traverse_cor(cor(traits) ^ 2, maxcor))

  rem = unique(rem)
  if(length(rem)>0)
  {
    cat(paste("\nBased on pairwise correlation: dropping ",length(rem),"of ",ncol(traits)," SNP(s) "))
    for(rm in rem)
      cat(rm," ")
    tnames=tnames[-rem]
    cat("\n")
  }


  ### Test based on variance inflation factor (VIF)
  mc=1
  da=0
  w=1
  if (length(tnames)<2)
    mc=0

  while(mc)
  {
    w=w+1
    k = length(tnames)
    data = data.frame(traits[,tnames])
    data = data[!apply(is.na(traits),1,any) ,]
    dv = rnorm(length(data[,1]))
    varinf = vif(lm( dv ~ . , data=data))
    if(any(!is.finite(varinf)))
    {
      cat("\n\n** WARNING ** \n** WARNING **  NaN produced in VIF routine:")
      cat(" reduce marker set.\n** WARNING **\n")
    }
    if ( max(abs(varinf),na.rm=T)>maxvif )
    {
      tmp = tnames
      ntmp = which.max(varinf)

#       if(da==0)
#         cat("\nBased on VIF: dropping SNP(s) ")

      da = 1
#       cat(tnames[ntmp]," ")
      tnames = tnames[-1*ntmp]
      if (length(tnames==1))
        mc <- 0
    }

    else
      mc <- 0
  }

#   if (da==1)
#     cat("\nKeeping",length(tnames),"SNPs\n\n")
  tnames

}




#######################################################
# VIF function tken from the CAR library, by J. Fox
# Generalized Variance-Inflation Factors
#######################################################

vif<-function(mod){
  #last modified 13 Dec 2000 by J. Fox
  UseMethod("vif")
}

vif.lm<-function(mod) {
  #last modified 2 Dec 2003 by J. Fox
  if (!is.null(weights(mod))) stop("requires unweighted lm")
  if(!has.intercept(mod)) stop("requires model with intercept.")
  terms<-term.names(mod)[-1]
  n.terms<-length(terms)
  if (n.terms < 2) stop("model contains fewer than 2 terms")
  R<-cor(model.matrix(mod)[,-1])
  detR<-det(as.matrix(R))
  result<-matrix(0,n.terms,3)
  rownames(result)<-terms
  colnames(result)<-c("GVIF","Df","GVIF^(1/2Df)")
  assign<-mod$assign
  for (term in 1:n.terms){
    subs<-which(assign==term)-1
    result[term,1]<-det(as.matrix(R[subs,subs]))*
      det(as.matrix(R[-subs,-subs]))/detR
    result[term,2]<-length(subs)
  }
  if (all(result[,2]==1)) result<-result[,1]
  else result[,3]<-result[,1]^(1/(2*result[,2]))
  result
}

vif.default<-function(mod){
  #last modified 13 Dec 2000 by J. Fox
  stop("requires lm object")
}


# Utility functions (J. Fox)
# last modified 19 Nov 04 by J. Fox

inv<-function(x) solve(x)

has.intercept<-function (model, ...) {
  UseMethod("has.intercept")
}

has.intercept.default<-function(model, ...) {
  any(names(coefficients(model))=="(Intercept)")
}

term.names<-function (model, ...) {
  UseMethod("term.names")
}

term.names.default<-function (model, ...) {
  term.names<-labels(terms(model))
  if (has.intercept(model)) c("(Intercept)", term.names)
  else term.names
}


