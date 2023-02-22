#' Basic description
#'
#' @description Given the sampless, this function returns the mean and weight matrix.
#' @usage Weight(G)
#' @param G a list of all groups
#' @return the mean and weight matrix
#' @examples
#' G=list()
#' set.seed(1)
#' n1=20;n2=40
#' N=n1+n2;
#' G[[1]]=c(1:n1);G[[2]]=c((n1+1):(n1+n2));
#' Weight(G)
#' @export
Weight=function(G){
  N=sum(lengths(G))
  m=function(x,y) return(2*x*y/N)
  m24=function(x,y) return(2*x*y/N+2*x*y*(x+y-2)/(N*(N-1))+4*x*(x-1)*y*(y-1)/(N*(N-1)))
  W=mu=matrix(NA,length(G),length(G))
  for(i in 1:(length(G)-1))
   for(j in (i+1):length(G)){
    mu[i,j]=m(length(G[[i]]),length(G[[j]]))
    W[i,j]=m24(length(G[[i]]),length(G[[j]]))-(mu[i,j])^2
                            }

  W=1/sqrt(W)
return(list(weight=W,mean=mu))
}


