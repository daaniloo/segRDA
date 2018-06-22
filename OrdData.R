#' Data ordering
#'
#' Ordinates both community and explanatory matrices based on on the first RDA score.
#' @param x explanatory matrix;
#' @param y community matrix;
#' @param ... parameters passed to vegan::\code{\link[vegan]{rda}};
#' @return a list consisting of:
#' \enumerate{
#' \item{envi}: the ordered explanatory matrix
#' \item{com}: the ordered community matrix
#' }
#' @author Danilo Candido Vieira
#' @examples
#' data(sim1)
#' sim1.o<-OrdData(x=sim1[[1]], y=sim1[[2]])
#' @export
OrdData<-function(x,y,...)
{

  if(is.null(rownames(y))){ rownames(y)<-1:nrow(y)}
  if(is.null(colnames(y))){ colnames(y)<-1:ncol(y)}
  if(is.null(rownames(x))){ rownames(x)<-1:nrow(x)}
  if(is.null(colnames(x))){ colnames(x)<-1:ncol(x)}
  model<-rda(y~.,data.frame(x),...)
  site.order<-order(scores(model,1,"sites"))
  sp.order<-order(scores(model,1,"species"))
    y.ord<-y[site.order,sp.order]
  x.ord<-x[site.order,]
  return(list(x=x.ord,y=y.ord))
}
