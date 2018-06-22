#' Breakpoints location
#'
#'This function is an auxiliary tool designed to identifies breakpoints along a dissimilarity profile table.
#' @param smw,pool an object either from class \code{smw} or from class \code{pool}, resulted from \code{\link{SMW.test}}   and \code{\link{PoolSMW}}, respectively.
#' @param seq.sig specifies the length of a sequence of significant dissimilarity values that will be considered in defining the community breakpoints. Defaults to \code{seq.sig=3};
#' @param peaks.choice defines if the breakpoints should be chosen as those samle positions corresponding to the maximum dissimilarity in the sequence  (\code{peak.choice="max"}) or as those sample positions corresponding to the median position of the sequence. Defaults to \code{peak.choice="max"};
#' @return  The function returns a subset of the original \code{DP} table containing only the samples suggested as breakpoints
#' @examples
#' data(sim1)
#' sim1.o<-OrdData(sim1)
#' sim1.smw<-SMW.test(y.ord=sim1.o$y,w=20)
#' BPlocation(sim1.smw)
#' @author Danilo Candido Vieira
#' @export
BPlocation <- function(smw,  seq.sig=3, peaks.choice=c("max","median")) UseMethod("BPlocation")

#' @export
BP<-function (smw,  seq.sig=3, peaks.choice=c("max","median"))
{

  peaks.choice<-match.arg(peaks.choice,c("max","median"))
  sig.name<-paste("sig:",smw$params$sig.test,sep="")
  SMWs<-smw$DP
  #which(unlist(lapply(SMWs,is.logical)))

  #grouping
  re<-split(SMWs, cumsum(c(1, diff(SMWs$sig) !=0)))
  groups<-re[unlist(lapply(re,function(x) sum(x$sig)!=0))]
  groups<-groups[ unlist(lapply(groups,nrow))>=seq.sig]
  point<-NULL
  if(peaks.choice=="max")
  { point<-unlist(lapply(groups,function(x) x[which(x[,3]==max(x[,3])),1]))}
  if(peaks.choice=="median"){point<-round(unlist(lapply(lapply(lapply(groups, "[",1),unlist),median)))}
  names(point)<-NULL
  BPs<-SMWs[SMWs[,1]%in%point,]
  return(BPs)

}

#' @describeIn BPlocation object from \code{\link{SMW.test}}
#' @export
BPlocation.smw<-function (smw,  ...){BP(smw,...)}

#' @export
#' @describeIn BPlocation object from \code{\link{PoolSMW}}
BPlocation.pool<-function (pool,  ...) {BP(smw,...)}
