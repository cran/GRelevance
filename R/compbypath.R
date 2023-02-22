#' Basic description
#'
#' @description Given the groups and the shortest Hamiltonian path, this function returns the number of edges that connect nodes between samples.
#' @usage compbypath(G,re.path)
#' @param G a list of all groups
#' @param re.path the shortest Hamiltonian path returned from the function Hpath
#' @return the number of edges that connect nodes between samples
#' @examples
#' d=100;n1=20;n2=30;n3=40;
#' N=n1+n2+n3
#' mu1=rep(0,d)
#' mu2=mu1
#' mu3=mu2+0.1
#' cov1=0.2^(abs(outer(1:d,1:d,"-")))
#' cov2=0.2^(abs(outer(1:d,1:d,"-")))
#' cov3=0.4^(abs(outer(1:d,1:d,"-")))
#' sam1=MASS::mvrnorm(n=n1,mu=mu1,Sigma=cov1)
#' sam2=MASS::mvrnorm(n=n2,mu=mu2,Sigma=cov2)
#' sam3=MASS::mvrnorm(n=n3,mu=mu3,Sigma=cov3)
#' Data=rbind(sam1,sam2,sam3)
#' Dist=philentropy::distance(Data, method = "euclidean")
#' Dist[lower.tri(Dist)] <- NA
#' Dist[diag(Dist)] <- NA
#' G=list()
#' G[[1]]=c(1:n1);G[[2]]=c((n1+1):(n1+n2));G[[3]]=c((n1+n2+1):(n1+n2+n3));
#' compbypath(G,Hpath(1,N,Dist))
#' @import MASS
#' @import philentropy
#' @seealso Hpath
#' @export
compbypath=function(G,re.path){
  Edge.c=matrix(NA,length(G),length(G))
  n0=dim(re.path)[1]+1
  for(g1 in 1:(length(G)-1))
   for(g2 in (g1+1):length(G)){
    Edge.c[g1,g2]=0;
    for(i in 1:dim(re.path)[1]){
      if((re.path[i,1]%in% G[[g1]]) && (re.path[i,2]%in% G[[g2]])) Edge.c[g1,g2]=Edge.c[g1,g2]+1;
      if((re.path[i,1]%in% G[[g2]]) && (re.path[i,2]%in% G[[g1]])) Edge.c[g1,g2]=Edge.c[g1,g2]+1;
                               }
                              }
return(list(EC=Edge.c))
}



