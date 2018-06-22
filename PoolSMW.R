#' Pooled split-moving-window analysis (SMW) with permutation tests
#Identifies community breakpoints using a pooled dissimilarity profile. The function performs SMW analyses using several windows sizes and averages the normalized dissimilarities (z-scores).
#' @param y.ord  the ordered community matrix;
#' @param windows a vector for the window sizes (all sizes have to be even);
#' @param n.rand the number of randomizations for significance computation;
#' @param z the critical value for significance test;
#' @return  An object of class \code{pool}, which is a list consisting of:
#' \itemize{
#' \item \code{DP}: the pooled dissimilarity profile. Dataframe giving the positions, labels, averaged z-score and the significance test for each sample.
#' \item \code{result.pooled}: a list of of \code{smw} objects (returned from \code{\link{SMW.test}}) corresponding to each window size analysed;
#' \item \code{params}: the input arguments.
#' }
#' @author Danilo Candido Vieira
#' @references
#' \enumerate{
#'   \item Erdos, L., Z. Bátori, C. S. Tölgyesi, and L. Körmöczi. 2014. The moving split window (MSW) analysis in vegetation science - An overview. Applied Ecology and Environmental Research 12:787–805.
#'   \item Cornelius, J. M., and J. F. Reynolds. 1991. On Determining the Statistical Significance of Discontinuities with Ordered Ecological Data. Ecology 72:2057–2070.
#' }
#' @examples
#' data(sim1)
#' sim1.o<-OrdData(sim1)
#' sim1.pool<-PoolSMW(y.ord=sim1.o$y,windows=c(20,30,40))
#' @export
PoolSMW<-function(y.ord,windows,n.rand=99,z=1.85)
{
  sig.test<-"AverageZ"
  argg <- c(as.list(environment()), list())
  if(is.null(rownames(y.ord))){rownames(y.ord)<-1:nrow(y.ord)}
  if(any(windows%%2==1)){stop("all Window sizes must be enven")}
  smw.tab<-data.frame(positions=(windows[1]/2):(nrow(y.ord)-(windows[1]/2)))
  smw.tab$SiteID<-rownames(y.ord)[smw.tab[,1]]
  Z.table<-z.scores<- diss<-smw.tab
  result.tests<-list()
  result.tests[[1]]<-SMW.test(y.ord,windows[1],progress=F, n.rand=n.rand)
   z.scores[which(smw.tab[,1]%in%result.tests[[1]]$DP[,1]),1+2]<- result.tests[[1]]$DP$zscore
  diss[which(smw.tab[,1]%in%result.tests[[1]]$DP[,1]),1+2]<- result.tests[[1]]$DP$diss
  message(paste("Progress", (round(1/length(windows)*100)), "%"))
  for( i in 2:length(windows))
  {
    w<-windows[i]
    result.tests[[i]]<-SMW.test(y.ord,w,progress=F, n.rand=n.rand)
    z.scores[which(smw.tab[,1]%in%result.tests[[i]]$DP[,1]),i+2]<- result.tests[[i]]$DP$zscore
    diss[which(diss[,1]%in%result.tests[[i]]$DP[,1]),i+2]<-result.tests[[i]]$DP$diss
    message(paste("Progress", (round(i/length(windows)*100)), "%"))
  }

   Z.table$AverageZ<-apply(z.scores[,-(1:2)],1,mean,na.rm=T)
   Z.table$sig<-Z.table$AverageZ>z
   colnames(Z.table)[ colnames(Z.table)=="sig"]<-paste("sig:",sig.test,sep="")
  smw.p<-list(DP=Z.table,result.pooled=result.tests,params=argg)
  class(smw.p) <- "pool"
  return(smw.p)
}

