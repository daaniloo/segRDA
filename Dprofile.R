#' Plot the dissimilarity profiles
#'
#' Generates a plot with dissimilarity values vs. the location of the window midpoint. If present,
#' significant dissimilarities and their peaks can be indicated.
#' @param smw,pool It can be either an object of class \code{'smw'} or an object of class \code{'pool'}
#'  (resulted from the functions \code{\link{SMW.test}} and \code{\link{PoolSMW}}, respectively).
#' @param BPs optional. A vector with the breakpoints location;
#' @param which the column from the dissimilarity profile table to be used in the plot. For objects of class \code{smw}: \code{"diss"} for dissimilarity values and \code{"zscore"} for Z-scores; for objects of class \code{pool}: the function automatically define \code{which="AverageZ"}
#' @param pch vector of length 3 specifying the symbols in plotting dissimilarity values, significant dissimilarities and breakpoints. Defaults to \code{pch = c(16,16,17)};
#' @param colors vector of length 3 specifying the colors of the plot in the same way as the pch argument. Defaults to \code{colors=c("black","red","blue")}
#' @param bg.colors vector of colors for the background according to the breakpoints. If \code{NULL} , the function uses the colour palette \code{\link[grDevices]{rainbow}} from R stats;
#' @param xlab The label for the x axis
#' @param yalb The label for the y axis
#' @param ... Further graphical parameters
#' @author Danilo Candido Vieira
#' @examples
#' data(sim1)
#' sim1.o<-OrdData(sim1)
#' sim1.smw<-SMW.test(y.ord=sim1.o$y,w=20)
#' Dprofile(sim1.smw)
#'
#' sim1.pool<-PoolSMW(y.ord=sim1.o$y,windows=c(20,30,40))
#' Dprofile(sim1.pool)
#' @export
Dprofile <- function(smw, BPs=NA,which=c("zscore","diss","AverageZ"),pch=c(16,16,17),colors=c("black","red","blue"),bg=NULL,bg.colors=NULL,
                     legend=T,...) UseMethod("Dprofile")
#' @export
DP<-function(smw, BPs=NA,which=c("zscore","diss"),pch=c(16,16,17),colors=c("black","red","blue"),bg=NULL,bg.colors=NULL,
                       legend=T,...)
{
  if(class(smw)=="smw")
  { which<-match.arg(which,c("zscore", "diss"))
  } else{which<-"AverageZ" }
  cps<-BPs
  positions=1:nrow(smw$params$y.ord)
  SMWs<-smw$DP
  values<-seq(min(SMWs[,which]),max(SMWs[,which]),length.out = length(positions))

  plot(positions,values, type="n",pch=pch[1],  col=colors[1],...)

  if(is.null(bg)==F)
  { if(anyNA(cps)){stop("BPs should be specified when bg is not NULL")}
    if(is.null(bg.colors)) {bg.colors<-rainbow(length(BPs)+1)}
    lim <- par("usr")
    bc<-c(lim[1],cps,lim[2])
    bg.adj<-adjustcolor(bg.colors, alpha.f = bg)
    for(i in 1:(length(bc)-1)){rect(bc[i], bc[1]-100, bc[i+1], bc[2], border = "white", col = bg.adj[i])}
    rect(lim[1],lim[3],lim[2],lim[4])}
  lines(SMWs[,c("positions",which)], type="b", pch=pch[1], col=colors[1])


  sig.col<-which((unlist(lapply(SMWs,is.logical)))==TRUE)
  points(SMWs[SMWs[,(sig.col)],c("positions",which)], pch=pch[2],col=colors[2])
  if(any(is.na(cps)!=T))
  { points(SMWs[SMWs$positions%in%cps,c("positions",which)], pch=pch[3],col=colors[3],cex=1.3)
    text(SMWs[SMWs$positions%in%cps,c("positions",which)], labels=paste(paste("B",1:length(cps),sep=""),"=", cps),cex=.7, adj=0,pos=2)  }
  sig.test<-smw$params$sig.test
  if(sig.test=="Z")
  {  sig.title<-paste("Significance test: \n",sig.test,"-score","(>=",smw$params$z,")", sep="")} else {
    sig.title<-paste("Significance test:\n",sig.test)
  }
  if(legend==T)
  {legend("topleft",legend=sig.title,bty="n", cex=.8)}



}


#' @export
#' @describeIn Dprofile object from \code{\link{SMW.test}}
Dprofile.smw<-function(smw,which=c("zscore","diss"),...){DP(smw,which=c("zscore","diss"),...)}


#' @export
#' @describeIn Dprofile object from \code{\link{PoolSMW}}
Dprofile.pool<-function(pool,which="'AverageZ'",...){DP(pool,which=c('AverageZ'),...)}




