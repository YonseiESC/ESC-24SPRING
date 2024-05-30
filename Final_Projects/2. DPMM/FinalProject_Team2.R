###<Data Generating>###
library("MASS")
set.seed(11)
n<-60 #각각 60개씩 sampling

#upper right
m1<-c(1.5, 1.5)
S1<-matrix(c(0.3, 0.05, 0.05, 0.3), ncol=2)
clus1<-mvrnorm(n, m1, S1)

#lower right
m2<-c(1.5, -1.5)
S2<-matrix(c(0.5, -0.08, -0.08, 0.2), ncol=2)
clus2<-mvrnorm(n, m2, S2)

#upper left
m3<-c(-1.5, 1.5)
S3<-matrix(c(0.1, 0.03, 0.03, 0.1), ncol=2)
clus3<-mvrnorm(n, m3, S3)

#lower left
m4<-c(-1.5, -1.5)
S4<-matrix(c(0.8, 0.50, 0.50, 0.8), ncol=2)
clus4<-mvrnorm(n, m4, S4)

datc<-rbind(clus1, clus2, clus3, clus4) # 60X4=>총 240개

#visualization
plot(clus1, col="grey",xlim=c(-3, 3), ylim=c(-3, 3), xlab="X", ylab="Y")
points(clus2, col="red")
points(clus3, col="blue")
points(clus4, col='green')



#Appendix B
###<Dirichlet Process Mixture Model>###
###<DPMM>###
library(mvtnorm)
crp_gibbs<-function(data, alpha, mu0, sigma0, sigma_y, c_init, maxIters=1000){
  data_dim<-ncol(data) #dimension of data points
  N<-nrow(data) #number of data points
  #data= NXD matrix of data points
  #sigma_y: measurement error of data y, assumed known, same for all clusters
  #mu0, sigma0: prior mean and variance around unknown mean mu
  #c_init: initial assignments of data points to clusters
  
  #Prior
  #y~N(mu, sigma_y)
  tau0<-solve(sigma0) #prior precision on mu, inverse of prior covariance
  tau_y<-solve(sigma_y)
  #CRP Gibs Sampler 시작
  z<-c_init
  n_k<-as.vector(table(z)) #각 클러스터에 속한 데이터 포인트의 수
  Nclust<-length(n_k) #현재 클러스터의 수
  #CRP Gibs Sampler
  res<-matrix(NA, nrow=N, ncol=maxIters) #cluster membership storage
  pb<-txtProgressBar(min=0, max=maxIters, style=3)
  for (iter in 1:maxIters) {
    for (n in 1:N) {
      c_i<-z[n] #n번째 사람의 table assignment가 어디인지
      n_k[c_i]<-n_k[c_i]-1 #remove nth person from table
      if (n_k[c_i]==0) { #If 테이블이 비면, replace new empty cluster with last cluster
        n_k[c_i]<-n_k[Nclust] #Nclust는 항상 already occupied tables
        loc_z<-(z==Nclust) #Nclust+1은 next new table to be created
        z[loc_z]<-c_i
        n_k<-n_k[-Nclust]
        Nclust<-Nclust-1
      }
      z[n]<- -1
      logp<-rep(NA, Nclust+1) #log probabilities for clusters, add unoccupied table
      for (c_i in 1:Nclust) { #loop over already occupied tables
        tau_p<-tau0 + n_k[c_i] * tau_y #calculate tau_p, cluster precision
        sig_p<-solve(tau_p) #cluster variance, inverse of tau_p
        loc_z<-which(z==c_i)
        if (length(loc_z)>1){
          sum_data<-colSums(data[z==c_i, ])
        }else{
          sum_data<-data[z==c_i, ]
        }
        #calculate posterior mean, mean_p
        mean_p<-sig_p %*% (tau_y %*% sum_data + tau0 %*% t(mu0))
        logp[c_i]<-log(n_k[c_i]) + dmvnorm(data[n,], mean=mean_p, 
                                           sigma=sig_p + sigma_y, log=TRUE)}
      #calculate log probability of unoccupied table
      #calculate probability for newly created cluste
      logp[Nclust+1]<-log(alpha) + dmvnorm(data[n,], mean=mu0, 
                                           sigma=sigma0 + sigma_y, log=TRUE)
      #transform unnormalized log probability to probability
      max_logp<-max(logp)
      logp<-logp - max_logp
      loc_probs<-exp(logp)
      loc_probs<-loc_probs/sum(loc_probs)
      #draw a sample of which cluster this new customer should belong to
      newz<-sample(1:(Nclust + 1), 1, replace=TRUE, prob=loc_probs)
      if (newz==Nclust + 1){
        n_k<-c(n_k, 0)      
        Nclust<-Nclust + 1
      }
      z[n]<-newz
      n_k[newz]<-n_k[newz] + 1
    }
    setTxtProgressBar(pb, iter)
    res[, iter]<-z
  }
  close(pb)
  invisible(res)
}



###<Find the Maximum A Posteriori (posterior mode) of the class membership estimates>###
alpha<-0.01
mu0<-matrix(rep(0,2), ncol=2, byrow=TRUE)
sigma0<-diag(2)*3^2
sigma_y<-diag(2)*1
c_init<-rep(1, nrow(datc))
results<-crp_gibbs(data=datc, alpha=alpha, mu0=mu0, sigma0=sigma0, sigma_y=sigma_y,
                   c_init=rep(1, nrow(datc)))
tab<-apply( results, 1, FUN=function(x){
  tab<-table(x)
  ans<-names(tab[which.max(tab)])
  return(ans)})
table(tab)



###<Effect of input values of parameters on CRP results>###
set.seed(11)
alpha<-c(0.01, 1.00, 3.00, 5.00)
sigma_y<-c(0.50, 1.00, 1.50, 3.00)
for (i in 1:length(alpha)){
  for (j in 1:length(sigma_y)){
    results<-crp_gibbs(data=datc, alpha=alpha[i], 
                       mu0=matrix(rep(0, 2), ncol=2, byrow=TRUE), 
                       sigma0 = diag(2)*3^2, sigma_y=diag(2)*sigma_y[j],
                       c_init=rep(1, nrow(datc)))
    tab<-apply(results, 1, FUN=function(x){
      tab<-table(x)
      ans<-names(tab[which.max(tab)])
      return(ans)})
    cat("alpha= ", alpha[i], "sigma_y= ", sigma_y[j], "\n")
    print(table(tab))
  }
}



###<Information on misclassifications>###
set.seed(11)
results<-crp_gibbs(data=datc, alpha=1.0, mu0=matrix(rep(0,2), ncol=2, byrow=TRUE),
                   sigma0=diag(2)*3^2, sigma_y=diag(2)*1.0, c_init=rep(1, nrow(datc)))
c_modal<-apply(results, 1, FUN=function(x){
  tab<-table(x)
  ans<-names(tab[which.max(tab)])
  return(ans)
})
c_true<-rep(1:4, each=60)
table(c_true, c_modal)



###<Estimating cluster means and covariance>###
mu0<-matrix(rep(0, 2), ncol=2, byrow=T)
sigma0<-matrix(c(3, 1, 1, 3), ncol=2)
sigma_y<-matrix(c(1, 0.2, 0.2, 1), ncol=2)
tau0<-solve(sigma0)
tau_y<-solve(sigma_y)
n_k<-table(c_modal)
for (c_i in 1:length(n_k)){
  tau_p<-tau0 + n_k[c_i]*tau_y
  sig_p<-solve(tau_p)
  y<-datc[which(c_modal==names(n_k)[c_i]), ]
  sum_data<-colSums(y)
  mean_p<-sig_p %*% (tau_y %*% sum_data + tau0 %*% t(mu0))
  y_ss<-sweep(y, MARGIN=2, STATS=mean_p, FUN="-")
  covariance<-(sigma0 + crossprod(y_ss))/(n_k[c_i] + 2-2-1)
  cat("\n cluster ", names(n_k)[c_i], "\n")
  print(round(mean_p, 3))
  print(round(covariance, 3))
}