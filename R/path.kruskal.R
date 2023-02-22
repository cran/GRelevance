#' Basic description
#'
#' @description Calculates the shortest Hamiltonian path based on the sorted edge weights and the nodes
#' @usage path.kruskal(nodes,edge_cost)
#' @param nodes sequence of nodes 1,...,n from the graph which is based on the high-dimensional data that is provided by the reader
#' @param edge_cost sorted edge weights
#' @return the shortest Hamiltonian path
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
#' mat=Dist
#' n1=1; n2=N; n0=n2-n1+1
#' edge.cost=matrix(NA,nrow=n0*(n0-1)/2,ncol=3)
#' temp=1;
#' for(i in n1:(n2-1))
#'  for(j in (i+1):(n2))
#'    {
#'     edge.cost[temp,3]=mat[i,j];edge.cost[temp,1]=i-n1+1;edge.cost[temp,2]=j-n1+1;temp=temp+1;}
#' edge.cost=edge.cost[sort.list(edge.cost[,3]), ]
#' path.kruskal(c(1:n0),edge.cost)
#' @seealso Hpath
#' @import MASS
#' @import philentropy
#' @export
path.kruskal=function (nodes, edge_cost)
{   n0=length(nodes)
cost=0
components <- matrix(c(nodes, nodes), ncol = 2)
shp <- matrix(ncol = 3)[-1, ]
edge_cost=rbind(matrix(ncol = 3)[-1, ],edge_cost)
i <- 1
degrees=rep(0,length(nodes))
while (nrow(shp) < n0 - 1) {
  min.shp <- edge_cost[i, ]
  iComp <- components[components[, 1] == min.shp[1], 2]
  jComp <- components[components[, 1] == min.shp[2], 2]
  if (iComp != jComp && max(c(degrees[min.shp[1]],degrees[min.shp[2]]))<2) {
    shp <- rbind(shp, min.shp)
    cost=cost+edge_cost[i, 3]
    components[components[, 2] == jComp, 2] <- iComp
    degrees[min.shp[1]]=degrees[min.shp[1]]+1;
    degrees[min.shp[2]]=degrees[min.shp[2]]+1;
  }
  i <- i + 1
}

return(list(edge.cost=shp))
}
