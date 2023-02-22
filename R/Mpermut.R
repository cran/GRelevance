#' Basic description
#'
#' @description Given the groups and the observed statistic, this function returns the pvalue.
#' @usage Mpermut(G,W,obs)
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
#' Z=(C-W$mean)*W$weight
#' obs=min(Z[!is.na(Z)])
#' Mpermut(G,W$weight,obs)
#' @import mvtnorm
#' @import MASS
#' @import philentropy
#' @export
Mpermut=function(G,W,obs){
  obs=rep(obs, length(G)*(length(G)-1)/2)
  L=function(i,j) return(j-i+(2*length(G)-i)*(i-1)/2)
  inv.L=function(ell){
         for(i in 1:(length(G)-1))
          for(j in (i+1):length(G))
            if(L(i,j)==ell) return(c(i,j));
                      }
  N=sum(lengths(G))
  m=function(x,y) return(2*x*y/N)
  m24=function(x,y) return(2*x*y/N+2*x*y*(x+y-2)/(N*(N-1))+4*x*(x-1)*y*(y-1)/(N*(N-1)))
  m14=function(x,y,z) return(2*x*y*z*(2*z-1)/(N*(N-1)))#z is the common group
  m04=function(x,y,z,u) return(4*x*y*z*u/(N*(N-1)))
  mu=rep(0,length(G)*(length(G)-1)/2)
  Sigma=matrix(NA,length(G)*(length(G)-1)/2,length(G)*(length(G)-1)/2)
  for(i in 1:(length(G)*(length(G)-1)/2)){
   v1=inv.L(i)
   if(W[v1[1],v1[2]]==0) obs[i]=-Inf
   if(W[v1[1],v1[2]]!=0) obs[i]=obs[i]/W[v1[1],v1[2]]
   for(j in 1:(length(G)*(length(G)-1)/2)){
    v2=inv.L(j);v12=intersect(v1,v2);
    m2=m(length(G[[v1[1]]]),length(G[[v1[2]]]))*m(length(G[[v2[1]]]),length(G[[v2[2]]]))
    if(length(v12)==2)
      Sigma[i,j]=m24(length(G[[v1[1]]]),length(G[[v1[2]]]))-m2
    if(length(v12)==1)
      Sigma[i,j]=m14(length(G[[v1[v1!=v12]]]),length(G[[v2[v2!=v12]]]),length(G[[v12]]))-m2
    if(length(v12)==0)
      Sigma[i,j]=m04(length(G[[v1[1]]]),length(G[[v1[2]]]),length(G[[v2[1]]]),length(G[[v2[2]]]))-m2
                                          }
                                         }
  prob=1-mvtnorm::pmvnorm(lower=obs,upper=rep(Inf, length(G)*(length(G)-1)/2), mean=mu, sigma=Sigma)
return(list(pvalue=prob[1],mean=as.vector(mu), sigma=Sigma))
}


