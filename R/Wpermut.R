#' Basic description
#'
#' @description Given the groups, the weight matrix and the observed statistic, this function returns the pvalue.
#' @usage Wpermut(G,W,obs)
#' @param G a list of all groups
#' @param W the weight matrix
#' @param obs the observed statistic
#' @return the pvalue
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
#' counts=compbypath(G,Hpath(1,N,Dist))
#' W=Weight(G)
#' #W[i,j]=0 #if we donot consider this relevance between sample i and sample j
#' C=counts$EC
#' WC=W$weight*C
#' WS=sum(WC[!is.na(WC)])
#' Wpermut(G,W$weight,WS)
#' @import MASS
#' @import philentropy
#' @export
Wpermut=function(G,W,obs){
  N=sum(lengths(G))
  m=function(x,y) return(2*x*y/N)
  m24=function(x,y) return(2*x*y/N+2*x*y*(x+y-2)/(N*(N-1))+4*x*(x-1)*y*(y-1)/(N*(N-1)))
  m14=function(x,y,z) return(2*x*y*z*(2*z-1)/(N*(N-1)))#z is the common group
  m04=function(x,y,z,u) return(4*x*y*z*u/(N*(N-1)))
  mu=sigma2=0
  for(g1 in 1:(length(G)-1))
   for(g2 in (g1+1):length(G))
    mu=mu+W[g1,g2]*m(length(G[[g1]]),length(G[[g2]]))

  for(i1 in 1:(length(G)-1))
   for(j1 in (i1+1):length(G))
     for(i2 in 1:(length(G)-1))
       for(j2 in (i2+1):length(G)){
         m2=m(length(G[[i1]]),length(G[[j1]]))*m(length(G[[i2]]),length(G[[j2]]))
         v1=c(i1,j1);v2=c(i2,j2);v12=intersect(v1,v2)
         if(length(v12)==2)
          sigma2=sigma2+W[i1,j1]*W[i2,j2]*(m24(length(G[[i1]]),length(G[[j1]]))-m2)
         if(length(v12)==1)
          sigma2=sigma2+W[i1,j1]*W[i2,j2]*(m14(length(G[[v1[v1!=v12]]]),length(G[[v2[v2!=v12]]]),length(G[[v12]]))-m2)
         if(length(v12)==0)
          sigma2=sigma2+W[i1,j1]*W[i2,j2]*(m04(length(G[[i1]]),length(G[[j1]]),length(G[[i2]]),length(G[[j2]]))-m2)
                                  }
return(list(pvalue=stats::pnorm((obs-mu)/sqrt(sigma2))))
}


