#' Plotting together the dissimilarity profiles from different windows sizes.
#'
#' Plot the dissimilarity profile using different windows sizes. This function applies only to
#' objects of class \code{"pool"}.  Through this function,  non-transformed profiles values
#'   from the different windows sizes are drawn on the same plot. The function returns an invisible frequency table of breakpoints, which can be displayed in the plot.
#'   See \code{\link{PoolSMW}}

#' @param pool An  object of class\code{"pool"} from \code{\link{PoolSMW}}
#' \itemize{
#'   \item \code{'Z'}: consider normalized dissimilarity (z-scores) discontinuities  that exceed a "z" critical value.
#'   \item \code{'SD'}: consider dissimilarity discontinuities that exceed mean plus one standard deviation
#'   \item \code{'SD2'}: consider dissimilarity discontinuities that exceed mean plus two standard deviation
#'   \item \code{'Tail1'}: consider dissimilarity discontinuities that exceed 95 percent confidence limits
#' }
#' @param title An overall title for the plot
#' @param legend Logical. Should the legend be displayed?
#' @param count Logical, Should the counting table appear in the plot?
#' @param xlab The label for the x axis
#' @param ylab The label for the y axis
#' @param lwd Line width
#' @param colors Vector of colors for the lines. If \code{NULL}, the function uses the colour palette \code{\link[grDevices]{rainbow}} from R stats
#' @param L.adj numeric (0-1). Adjustment of the position of the legend along the axis "x"
#' @param ... Further arguments passed from BPlocation function
#' @return Returns an invisible table with the counts of breakpoints
#' @author Danilo Candido Vieira
#' @examples
#' data(sim1)
#' sim1.o<-OrdData(sim1)
#' sim1.pool<-PoolSMW(y.ord=sim1.o$y,windows=c(20,30,40))
#' count.breaks<-Dprofile2(sim1.pool)
#' count.breaks
#' @export
Dprofile2<-function(pool, which=c('diss','zscore'),title="Pooled SMW",count=TRUE,
                    legend=T,colors=NULL,xlab=NULL,ylab=NULL,lwd=2,L.adj=-0.4, ...)
{
  smw.p<-pool
  which<-match.arg(which,c('diss','zscore'))
  if(class(smw.p)!="pool"){stop("smw.p should be of class 'pooled'")}
  y.ord<-smw.p$params$y.ord
  diss.list<-lapply(smw.p$result.pooled, function(x)(x$DP[,c("positions",which)]))
  cps.list<-lapply(smw.p$result.pooled,BPlocation)
  cps.list<-cps.list[unlist(lapply(cps.list, nrow))>0]
  freq.cps<-unlist(lapply(cps.list,function(x) x[,1]))
  freq.cps<-table(freq.cps)
  windows<-smw.p$params$windows
  y.limits<-range(na.omit(unlist(lapply(diss.list, function(x) x[,2]))))
  if(legend==T){
  par(mar=c(5.1,4.1,4.1,8.1), xpd=TRUE)}
  plot(1:nrow(y.ord),seq(min(y.limits),max(y.limits),length.out = nrow(y.ord)), type="n",las=1, xlab=" ",
       ylab="",main=title)
  if(is.null(ylab)){ylab=colnames(diss.list[[1]])[2]}
  if(is.null(xlab)){xlab="sample positions"}
  mtext(xlab,1,3)
  mtext(ylab,2,3)
  if(is.null(colors)==T){cols<-rainbow(length(windows))} else{ cols=colors}

  cps.list<-lapply(smw.p$result.pooled,BPlocation)
  for( i in 1:length(windows))
  { lines(diss.list[[i]], col=cols[i],lwd=lwd)
  if(count==T){
    if(sum(unlist(lapply(cps.list, nrow)))==0){stop("no breakpoints identified in any of the windows")}
    points(cps.list[[i]][,c("positions",which)], pch=16,col=cols[i])}}
  if(legend==T)
  { legend("topright" , inset=c(L.adj,0), col=cols,lty=1,legend=paste("w",windows,sep=""),y.intersp=.75, bty="n")
    if(count==T){
      legend("bottomright", inset=c(L.adj,0),legend=paste("sample", names(freq.cps),"=", freq.cps),
             cex=.8, bty="o",title="Counts of breakpoints", xjust=.5, yjust=0)}}
return(invisible(freq.cps))
}

