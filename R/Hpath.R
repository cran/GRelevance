#' Basic description
#'
#' @description Applies the path.kruskal function based on the nodes and edge.cost (sorts the weights from minimum to maximum). Given the starting node, ending node, and the distance matrix, this function returns the list of nodes of each edge from the shortest Hamiltonian path. We have the Hamiltonian path from path.kruskal
#' @usage Hpath(n1,n2,mat)
#' @param n1 starting node
#' @param n2 ending node
#' @param mat distance matrix (distance type is determined by the reader)
#' @return list of nodes of each edge from the shortest Hamiltonian path
#' @examples
#' G=list()
#' set.seed(1)
#' n1=20;n2=40
#' N=n1+n2;
#' G[[1]]=c(1:n1);G[[2]]=c((n1+1):(n1+n2));
#' d=10
#' mu1=rep(0,d)
#' mu2=mu1+0.1
#' true.cov1=0.4^(abs(outer(1:d,1:d,"-")))
#' true.cov2=0.4^(abs(outer(1:d,1:d,"-")))
#' sam1=MASS::mvrnorm(n=n1,mu=mu1,Sigma=true.cov1)
#' sam2=MASS::mvrnorm(n=n2,mu=mu2,Sigma=true.cov2)
#' Data=rbind(sam1,sam2)
#' Dist=philentropy::distance(Data, method = "euclidean")
#' Dist[lower.tri(Dist)] <- NA
#' Dist[diag(Dist)] <- NA
#' Hpath(1,N,Dist)
#' @seealso path.kruskal
#' @import MASS
#' @import philentropy
#' @export
Hpath=function(n1,n2,mat){
  n0=n2-n1+1
  edge.cost=matrix(NA,nrow=n0*(n0-1)/2,ncol=3)

  temp=1;
  for(i in n1:(n2-1))
    for(j in (i+1):(n2))
    {
      edge.cost[temp,3]=mat[i,j];edge.cost[temp,1]=i-n1+1;edge.cost[temp,2]=j-n1+1;temp=temp+1;}

  edge.cost=edge.cost[sort.list(edge.cost[,3]), ]
  #based on path

  opt.path=path.kruskal(c(1:n0),edge.cost)
  return(opt.path$edge.cost[,1:2])
}

