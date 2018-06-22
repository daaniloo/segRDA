#' Split moving window analysis (SMW) with randomization tests
#'
#' This function generates a dissimilarity profile with randomizations tests using a single
#' window size \code{w}.
#' @param y.ord  the ordered community matrix;
#' @param w a numeric value for window size (has to be even);
#' @param dist dissimilarity index used in \code{\link[vegan]{vegdist}}. Defaults to \code{'bray'}.
#' @param n.rand the number of randomizations for significance computation;
#' @param rand The type of randomization for significance computation (Erdös et.al, 2014):
#' \itemize{
#'   \item \code{'shift'} : restricted randomization in which data belonging to the same species are randomly shifted along the data series ("Random shift");
#'   \item \code{'plot'} : unrestricted randomization: each sample is randomly repositioned along the data series ("Random plot").
#' }
#' @param sig.test significance test used to test whether a detected discontinuity differs significantly from those appearing in a random pattern.
#' The following tests are considered with default to \code{sig.test="Z"}:
#' \itemize{
#'   \item \code{'Z'}: consider normalized dissimilarity (z-scores) discontinuities  that exceed a "z" critical value.
#'   \item \code{'SD'}: consider dissimilarity discontinuities that exceed mean plus one standard deviation;
#'   \item \code{'SD2'}: consider dissimilarity discontinuities that exceed mean plus two standard deviation;
#'   \item \code{'Tail1'}: Consider dissimilarity discontinuities that exceed 95 percent confidence limits.
#' }
#' @param z the critical value for the significance computation Defaults to \code{'z=1.85'} (Erdös et.al, 2014);
#' @param progress logical. If \code{TRUE} displays the progress of the analysis.
#' @return  An object of class \code{"smw"}, which is a list consisting of:
#' \enumerate{
#' \item \code{DP}: the dissimilarity profile.  Dataframe giving the positions, labels, values of dissimilarity, z-scores and,  significance test for each sample;
#' \item \code{DP.rand}: dataframe containing the randomized dissimilarity profile;
#' \item \code{SD} standard deviation for each sample position;
#' \item \code{D.overall}: overall expected mean dissimilarity;
#' \item \code{SD.overall}: average standard deviation for the dissimilarities;
#' \item \code{params}: list with input arguments.
#' }
#' @author Danilo Candido Vieira. Fabio Cop, Gustavo Fonseca
#' @references
#' \itemize{
#'   \item Erdos, L., Z. Bátori, C. S. Tölgyesi, and L. Körmöczi. 2014. The moving split window (MSW) analysis in vegetation science - An overview. Applied Ecology and Environmental Research 12:787–805.
#'   \item Cornelius, J. M., and J. F. Reynolds. 1991. On Determining the Statistical Significance of Discontinuities with Ordered Ecological Data. Ecology 72:2057–2070.
#' }
#' @examples
#' data(sim1)
#' sim1.o<-OrdData(sim1)
#' sim1.smw<-SMW.test(y.ord=sim1.o$y,w=20)
#' @export
SMW.test<-function(y.ord,w,dist="bray",n.rand=99,rand=c("shift","plot"),
                   sig.test=c("Z","SD","SD2",'Tail1'), z=1.85,
                  progress=F)
{
  if(n.rand<2){stop ("number of randomizations not alowed")}
  rand<-match.arg(rand,c("shift","plot"))
  sig.test<-match.arg(sig.test,c("Z","SD","SD2",'Tail1'))
  argg <- c(as.list(environment()), list())
  DPtable<-smw(y.ord,w,dist)
  OB<-DPtable[,3]
  sig.name<-paste("sig:",sig.test,sep="")
  DP.rand<-matrix(NA,length(OB),n.rand)
  for( b in 1:n.rand)
  {
    if(rand=='shift')
    { comm.rand<-apply(y.ord,2,function(sp)sp[sample(1:length(sp))])
      DP.rand[,b]<-smw(data.frame(comm.rand),w,dist)[,3]}
    if(rand=='plot')
    { comm.rand<-t(apply(y.ord,1,function(sp)sp[sample(1:length(sp))]))
      DP.rand[,b]<-smw(comm.rand,w,dist)[,3]}
    if(progress==T)
    {message(paste("Progress", (round(b/n.rand*100)), "%"))}}
  rownames(DP.rand)<-DPtable[,1]
  Dmean<-apply(DP.rand,1,mean)
  SD<-apply(DP.rand,1,sd)
  D.overall<-sum(Dmean)/(nrow(y.ord)-w)
  SD.overall<-sum(SD)/(nrow(y.ord)-w)
  Dz<-(OB-D.overall)/SD.overall
  DPtable$zscore<-Dz
    #significance
  Sigif.Tests<-DPtable[,1:2]
  Sigif.Tests$SD<- OB>Dmean+(SD)
  Sigif.Tests$SD2<-  OB>Dmean+(2*SD)
  Sigif.Tests$Z<-Dz>z
  p.value<-NULL
  for(i in 1:length(OB))
  {p.value[i]<-pnorm(OB[i],mean=Dmean[i], sd=SD[i], lower.tail = F)}
  Sigif.Tests$Tail1<-p.value<0.05
  Sigif.Tests$p.value<-round(p.value,10)
  DPtable$sig<-  Sigif.Tests[,which(colnames(Sigif.Tests)==sig.test)]
  colnames(DPtable)[ colnames(DPtable)=="sig"]<-sig.name
  if(sig.test=="Tail1"){DPtable$p.value<-round(p.value,10)}
  smw.o<-list(DP= DPtable,DP.rand=DP.rand,Dmean=Dmean,
                   SD=SD,D.overall=D.overall,SD.overall=SD.overall, params=argg)
  class(smw.o) <- "smw"
  return(smw.o)
}
#' @export
smw<-function(y.ord,w, dist="bray")
{
  if(w%%2==1){stop('window size should be even')}
  diss<-NULL
  y.ord<-data.frame(y.ord)
  for(i in 1:(nrow(y.ord)-w+1)){
    wy.ord<-y.ord[i:(i+(w-1)),]

    half.a<-apply(wy.ord[1:(nrow(wy.ord)/2),],2,sum)
    half.b<-apply(wy.ord[-c(1:(nrow(wy.ord)/2)),],2,sum)
    d<-vegdist(rbind(half.a,half.b),dist)
    diss[i]<-d
    k<-(w/2)
    for(i in 1:((nrow(y.ord)-w)))
    {k[i+1]<-(w/2)+i}
  }
  result<-data.frame(positions=k,sampleID= rownames(y.ord)[k],diss=diss)
  return(invisible(result))
}

