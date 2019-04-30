set.seed(12345)
N <-500
J <- 20
a.true <- runif(J,0.5,2)
b.true <- rnorm(J,0,1)
theta.true <- rnorm(N,0,1)
tau.true <- rnorm(N,0,1)
beta.true <- runif(1,0.5,2)
gamma.true <- rnorm(J,0,1)
alpha.true <- runif(J,0.5,2)
sigma.true <- 1/alpha.true
delta.true <- runif(1,0.5,2)

set.seed(3)
#runif(1,0.5,1.5)
mean.lognorm <- outer(-delta.true*tau.true-beta.true*theta.true,beta.true*b.true+gamma.true,"+")
log.response.time <- rnorm(N*J,0,1)*matrix(rep(sigma.true,N),N,byrow=TRUE)+mean.lognorm
response.time <-exp(log.response.time)
time.sum <- rowSums(response.time)
C <-as.vector( quantile(time.sum)[4] )
D <- matrix(0,N,J)
cum.response.time <- t(apply(response.time ,1,cumsum))
D[cum.response.time <= C] <- 1
response.time.observed <- matrix(NA,N,J)
response.time.observed[cum.response.time <=C] <- response.time[cum.response.time <=C]

response.probability <- 1/(1+exp(-outer(theta.true,a.true,"*")+matrix(rep(a.true*b.true,N),N,byrow=TRUE)))
response.data.com <-matrix(rbinom(N*J, 1, response.probability) ,N,J)
response.data <- response.data.com 
response.data[D<1] <- 9
XIJ <- response.data
XIJ[response.data>1] <- 1
YIJ <- 1-response.data
YIJ[response.data>1] <- 1
par<- list(c(),c(),c(),c(),c(),c())
names(par) <- c("a","b","alpha","beta","gamma","delta")
M <- 500
A <- B <- ALPHA <- GAMMA <-  matrix(NA,M+1,J)
BETA <- DELTA <- rep(NA,M+1)
par$a <-rep(1,J) 
#a.true
A[1,] <- par$a

par$b <- sort(rnorm(J,0,1))[rank(colMeans(YIJ))]

B[1,] <- par$b
par$alpha <-  rep(1,J)
# alpha.true 
ALPHA[1,] <- par$alpha

par$beta <-  runif(1,0.5,1.5)
# beta.true
BETA[1] <- par$beta
par$delta <- runif(1,0.5,1.5)
#delta.true

DELTA[1] <- par$delta

par$gamma <-  sort(rnorm(J,0,1))[rank(colMeans(response.time.observed,na.rm=TRUE) )]
GAMMA[1,] <- par$gamma

K <- 20
theta.max <- 4
delta.theta <- 2*theta.max/K
theta.seq <- (seq(-theta.max,theta.max,delta.theta)+delta.theta/2)[1:K]
tau.seq <- theta.seq
pi.theta <- dnorm(tau.seq,0,1)
omega.tau <- pi.theta


TIJ <- response.time.observed
TIJ[is.na(TIJ)] <- 1
response.time.observed.cum <-  t(apply(response.time.observed ,1,cumsum))

V <- rowSums(D)+1
S <- rowSums(D)
response.time.leave <- rep(0,N)
for(i in 1:N){
  if (S[i]>0){
    response.time.leave[i] <- C-response.time.observed.cum[i,S[i]]}
  else {response.time.leave[i] <- C}
}

response.time.leave.J <-response.time.leave [V<J+1]
log.response.time.leave <- rep(0,N)
log.response.time.leave[V<J+1] <- log(response.time.leave.J)
VV <- V
VV[VV>J] <- 0
VJ <- VV[VV>0]
VVV <- VV%x%rep(1,K)

V.MATRIX <-t(apply(1-D,1,cumsum)) 
V.MATRIX[V.MATRIX>1] <- 0

V.MATRIX <-t(apply(1-D,1,cumsum)) 
V.MATRIX[V.MATRIX>1] <- 0
step.size.a <- step.size.b <-  step.size.alpha <- step.size.gamma <- rep(1,J)
round.a <-round.b <- round.alpha <-round.beta <- round.gamma <-rep(0,J)
step.size.delta <- step.size.beta <- 1

responseitem <- colSums(D)
weightitem <- 1/responseitem
weightall <- 1/J
H1 <- function(x1,x2,j){
  pkj<- 1/(1+exp(-x1*theta.seq+x1*x2))
  pkj[pkj<0.0001] <- 0.0001
  pkj[pkj>0.9999] <- 0.9999
  L1 <-sum(((outer(XIJ[,j],log(pkj),"*")+outer(1-XIJ[,j],log(1-pkj),"*"))*post.theta*matrix(rep(D[,j] ,K),N,byrow=FALSE)) )*weightitem[j]
  return(L1)
}
H2 <- function(y1,y2,y3,y4,y5,y6,j){
  pkj2<- 1/(1+exp(-y1*theta.seq+y1*y2))
  pkj2[pkj2<0.0001] <- 0.0001
  pkj2[pkj2>0.9999] <- 0.9999
  pnorm.t <- pnorm(y3*(log.response.time.leave-y4*y2-y5)%x%matrix(1,K,K)+ matrix(rep(t(outer(y4*theta.seq,y6*tau.seq,"+")),N),ncol=K,byrow=T),0,1)
  pnorm.t[pnorm.t>0.9999] <- 0.9999
  L2 <-(sum(((outer(XIJ[,j],log(pkj2),"*")+outer(1-XIJ[,j],log(1-pkj2),"*"))*post.theta*matrix(rep(D[,j] ,K),N,byrow=FALSE)) )-0.5*y3*y3*sum(STAND.POST*(D[,j]%x%matrix(1,K,K))*( (log(TIJ[,j])-y4*y2-y5)%x%matrix(1,K,K)+ matrix(rep(t(outer(y4*theta.seq,y6*tau.seq,"+")),N),ncol=K,byrow=T))*( (log(TIJ[,j])-y4*y2-y5)%x%matrix(1,K,K)+ matrix(rep(t(outer(y4*theta.seq,y6*tau.seq,"+")),N),ncol=K,byrow=T)))+sum(STAND.POST*(V.MATRIX[,j]%x%matrix(1,K,K))*log(1-pnorm.t)))*weightitem[j]
  return(L2)
}

H3 <- function(z2,z3,z4,z5,z6,j){
  pnorm.t3 <- pnorm(z3*(log.response.time.leave-z4*z2-z5)%x%matrix(1,K,K)+ matrix(rep(t(outer(z4*theta.seq,z6*tau.seq,"+")),N),ncol=K,byrow=T),0,1)
  pnorm.t3[pnorm.t3>0.9999] <- 0.9999
  L3 <- (-0.5*z3*z3*sum(STAND.POST*(D[,j]%x%matrix(1,K,K))*( (log(TIJ[,j])-z4*z2-z5)%x%matrix(1,K,K)+ matrix(rep(t(outer(z4*theta.seq,z6*tau.seq,"+")),N),ncol=K,byrow=T))*( (log(TIJ[,j])-z4*z2-z5)%x%matrix(1,K,K)+ matrix(rep(t(outer(z4*theta.seq,z6*tau.seq,"+")),N),ncol=K,byrow=T)))+sum(STAND.POST*(V.MATRIX[,j]%x%matrix(1,K,K))*log(1-pnorm.t3)))*weightitem[j]+log(z3)
  return(L3)
}

H4 <- function(t2,t3,t4,t5,t6){
  L4 <- rep(NA,J)
  for (j in 1:J){
    L4[j] <- H3(t2[j],t3[j],t4,t5[j],t6,j)
  }
  LL4 <- weightall* sum(L4)
  return(LL4)
}


rho <- 0.5


for (m in 1:M){
  LEA <- matrix(1,N*K,K)
  
  post.theta <- matrix(0,N,K)
  post.tau <- matrix(0,N,K)
  
  P.K1J <-  1/(1+exp(-outer(theta.seq,par$a,"*")+matrix(rep(par$a*par$b,K),K,byrow=TRUE)))
  
  P.NK1J <- (matrix(rep(t(P.K1J),N),ncol=J,byrow=TRUE))*(XIJ%x%rep(1,K))+(matrix(rep(t(1-P.K1J),N),ncol=J,byrow=TRUE))*((YIJ)%x%rep(1,K))
  PROB.NK1 <- matrix(rep(apply(P.NK1J, 1, prod),K),ncol=K,byrow=F) #
  alphasquare <- par$alpha*par$alpha
  
  
  P.N.K1K2 <- exp(-0.5*rowSums(D*(matrix(rep(alphasquare,N),N,byrow=TRUE)))%x%matrix(rep(tau.seq*tau.seq,K),K,byrow = T)*par$delta*par$delta-(rowSums(D*(matrix(rep(alphasquare *par$beta*par$delta,N),N,byrow=TRUE))))%x%(outer(theta.seq,tau.seq))-0.5*rowSums(D*(matrix(rep(alphasquare *par$beta*par$beta,N),N,byrow=TRUE)))%x%matrix(rep(theta.seq*theta.seq,K),K,byrow =F)-rowSums(D*(matrix(rep(alphasquare ,N),N,byrow=TRUE)*log(TIJ)- matrix(rep((par$beta*par$b+par$gamma)*alphasquare ,N),N,byrow=TRUE)))%x%(outer(par$beta*theta.seq,par$delta*tau.seq,"+")))
  
  LEA[VVV>0,] <- 1-pnorm(matrix(rep(c(t(matrix(rep(par$alpha[VJ]*(log(response.time.leave.J)-par$beta*par$b[VJ]-par$gamma[VJ]),K),ncol=K,byrow = F)+outer(par$alpha[VJ],par$beta*theta.seq)
  )),K),ncol=K,byrow =F)+ outer(par$alpha[VJ],par$delta*tau.seq)%x%rep(1,K),0,1)
  
  
  JOINT.POST <- PROB.NK1*P.N.K1K2*LEA*matrix(rep(t(outer(pi.theta,omega.tau)),N),ncol=K,byrow=TRUE)
  STAND.POST <- matrix(0,N*K,K)
  
  for (i in 1:N){
    AA <- JOINT.POST[(i*K-K+1):(i*K),]
    AA[AA>1e+300]<-1e+300
    total <- sum(AA)
    STAND.POST[(i*K-K+1):(i*K),] <- AA/total
    post.theta[i,]<- rowSums(AA)/total
    post.tau[i,] <- colSums(AA)/total
  }
  HH3 <- rep(NA,1,J)
  roud.beta <- rep(NA,1,J)
  round.delta <- rep(NA,1,J)
  for (j in 1:J){
    
    
    log.mean.tj <- (log(TIJ[,j])-par$beta*par$b[j]-par$gamma[j])%x%matrix(1,K,K)+ matrix(rep(t(outer(par$beta*theta.seq,par$delta*tau.seq,"+")),N),ncol=K,byrow=T)
    leave.tj <- (log.response.time.leave-par$beta*par$b[j]-par$gamma[j])%x%matrix(1,K,K)+ matrix(rep(t(outer(par$beta*theta.seq,par$delta*tau.seq,"+")),N),ncol=K,byrow=T)
    fnor <- 1-pnorm(par$alpha[j]*leave.tj,0,1)
    fnor[fnor<0.0001] <- 0.0001
    nor <- dnorm(par$alpha[j]*leave.tj,0,1)/fnor
    
    round.a[j] <-(sum(outer(XIJ[,j],P.K1J[,j],"-")*(matrix(rep(theta.seq,N),N,byrow=TRUE)-par$b[j])*post.theta*matrix(rep(D[,j] ,K),N,byrow=FALSE)))*weightitem[j] 
    HH1 <- H1(par$a[j],par$b[j],j)
    while ( H1((par$a[j]+step.size.a[j]*round.a[j]),par$b[j],j)-HH1<0.5*step.size.a[j]*round.a[j]*round.a[j]) {
      step.size.a[j] <- rho*  step.size.a[j]
    }
    
    round.b[j] <- (sum(-par$a[j]*outer(XIJ[,j],P.K1J[,j],"-")*post.theta*matrix(rep(D[,j] ,K),N,byrow=FALSE))+alphasquare[j]*par$beta*sum(STAND.POST*(D[,j]%x%matrix(1,K,K))*log.mean.tj) +par$alpha[j]*par$beta*sum(STAND.POST*(V.MATRIX[,j]%x%matrix(1,K,K))*nor))*weightitem[j]
    HH2 <- H2(par$a[j],par$b[j],par$alpha[j],par$beta,par$gamma[j],par$delta,j)
    while ( H2(par$a[j],(par$b[j]+step.size.b[j]*round.b[j]),par$alpha[j],par$beta,par$gamma[j],par$delta,j)-HH2<0.5*step.size.b[j]*round.b[j]*round.b[j]) {
      step.size.b[j] <- rho*  step.size.b[j]
    }
    
    HH3[j] <- H3(par$b[j],par$alpha[j],par$beta,par$gamma[j],par$delta,j)
    
    round.alpha[j] <-1/par$alpha[j]- (par$alpha[j]*sum(STAND.POST*(D[,j]%x%matrix(1,K,K))*log.mean.tj*log.mean.tj)-sum(STAND.POST*(V.MATRIX[,j]%x%matrix(1,K,K))*nor*leave.tj))*weightitem[j] 
    
    while(par$alpha[j]+step.size.alpha[j]*round.alpha[j]<0){
      step.size.alpha[j] <- rho*  step.size.alpha[j]
    }
    while ( H3(par$b[j],(par$alpha[j]+step.size.alpha[j]*round.alpha[j]),par$beta,par$gamma[j],par$delta,j)-HH3[j]<0.5*step.size.alpha[j]*round.alpha[j]*round.alpha[j]) {
      step.size.alpha[j] <- rho*  step.size.alpha[j]
    }  
    
    round.gamma[j] <- (alphasquare[j]*sum(STAND.POST*(D[,j]%x%matrix(1,K,K))*log.mean.tj) +par$alpha[j]*sum(STAND.POST*(V.MATRIX[,j]%x%matrix(1,K,K))*nor))*weightitem[j]
    while ( H3(par$b[j],par$alpha[j],par$beta,(par$gamma[j]+step.size.gamma[j]*round.gamma[j]),par$delta,j)-HH3[j]<0.5*step.size.gamma[j]*round.gamma[j]*round.gamma[j]) {
      step.size.gamma[j] <- rho*  step.size.gamma[j]
    }
    
    
    round.beta[j] <- (alphasquare[j]*sum(STAND.POST*(par$b[j]-matrix(rep(theta.seq,N*K),ncol=K,byrow = F))*(D[,j]%x%matrix(1,K,K))*log.mean.tj) +par$alpha[j]*sum(STAND.POST*(par$b[j]-matrix(rep(theta.seq,N*K),ncol=K,byrow = F))*(V.MATRIX[,j]%x%matrix(1,K,K))*nor))*weightitem[j]
    
    round.delta[j] <- -(alphasquare[j]*sum(STAND.POST*(matrix(rep(tau.seq,N*K),ncol=K,byrow = T))*(D[,j]%x%matrix(1,K,K))*log.mean.tj)-par$alpha[j]*sum(STAND.POST*(matrix(rep(tau.seq,N*K),ncol=K,byrow = T))*(V.MATRIX[,j]%x%matrix(1,K,K))*nor))*weightitem[j]
    
  }
  
  
  HH3.all <- weightall*sum(HH3)
  round.beta.all <- weightall*sum(round.beta)
  
  while ( H4(par$b,par$alpha,(par$beta+step.size.beta* round.beta.all),par$gamma,par$delta)-HH3.all<0.5*step.size.beta*round.beta.all*round.beta.all) {
    step.size.beta<- rho*  step.size.beta
  }
  
  
  
  round.delta.all <- weightall*sum(round.delta)
  while ( H4(par$b,par$alpha,par$beta,par$gamma,(par$delta+step.size.delta*round.delta.all))-HH3.all<0.5*step.size.delta*round.delta.all*round.delta.all) {
    step.size.delta<- rho*  step.size.delta
  }
  
  par$a <- A[m,]+step.size.a * round.a 
  A[m+1,] <- par$a 
  
  par$b <- B[m,]+step.size.b * round.b
  B[m+1,] <- par$b
  
  par$alpha <- ALPHA[m,]+step.size.alpha * round.alpha 
  ALPHA[m+1,] <- par$alpha
  
  
  par$gamma <- GAMMA[m,]+step.size.gamma * round.gamma
  GAMMA[m+1,] <- par$gamma
  
  par$beta <- BETA[m]+step.size.beta * round.beta.all
  BETA[m+1] <- par$beta
  par$delta <- DELTA[m]+step.size.delta * round.delta.all
  DELTA[m+1] <- par$delta
}

t2 <- Sys.time()
t2-t1

mean(abs(a.true-par$a))
mean(abs(b.true-par$b))
mean(abs(alpha.true-par$alpha))
mean(abs(gamma.true-par$gamma))
par$beta-beta.true
par$delta - delta.true

i.a <-i.b <- i.alpha <-i.beta <- i.gamma <-i.beta <- i.delta <- rep(NA,J)
for (j in 1:J){
  log.mean.tj1 <- (log(TIJ[,j])-par$beta*par$b[j]-par$gamma[j])%x%matrix(1,K,K)+ matrix(rep(t(outer(par$beta*theta.seq,par$delta*tau.seq,"+")),N),ncol=K,byrow=T)
  leave.tj1 <- (log.response.time.leave-par$beta*par$b[j]-par$gamma[j])%x%matrix(1,K,K)+ matrix(rep(t(outer(par$beta*theta.seq,par$delta*tau.seq,"+")),N),ncol=K,byrow=T)
  fnor1 <- 1-pnorm(par$alpha[j]*leave.tj1,0,1)
  fnor1[fnor<0.0001] <- 0.0001
  nor1 <-  dnorm(par$alpha[j]*leave.tj1,0,1)/fnor1
  
  MA <-nor1*nor1- par$alpha[j]*leave.tj1*nor1
  
  i.a[j] <-(sum(matrix(rep(P.K1J[,j]*(1-P.K1J[,j]),N),N,byrow=TRUE)*(matrix(rep(theta.seq,N),N,byrow=TRUE)-par$b[j])*(matrix(rep(theta.seq,N),N,byrow=TRUE)-par$b[j])*post.theta*matrix(rep(D[,j] ,K),N,byrow=FALSE)))
  i.b[j] <- (sum(par$a[j]*matrix(rep(P.K1J[,j]*(1-P.K1J[,j]),N),N,byrow=TRUE)*post.theta*matrix(rep(D[,j] ,K),N,byrow=FALSE)))+responseitem[j]*par$alpha[j]*par$alpha[j]*par$beta*par$beta+par$alpha[j]*par$alpha[j]*par$beta*par$beta*sum(STAND.POST*(V.MATRIX[,j]%x%matrix(1,K,K))*MA)
  i.gamma[j] <- responseitem[j]*par$alpha[j]*par$alpha[j]+par$alpha[j]*par$alpha[j]*sum(STAND.POST*(V.MATRIX[,j]%x%matrix(1,K,K))*MA)
  i.alpha[j] <- responseitem[j]*log(par$alpha[j])+sum(STAND.POST*(D[,j]%x%matrix(1,K,K))*log.mean.tj1*log.mean.tj1)+sum(STAND.POST*(V.MATRIX[,j]%x%matrix(1,K,K))*MA*leave.tj1*leave.tj1)
  i.delta[j] <- par$alpha[j]*par$alpha[j]*sum(STAND.POST*(matrix(rep(tau.seq*tau.seq,N*K),ncol=K,byrow = T))*(D[,j]%x%matrix(1,K,K)))+par$alpha[j]*par$alpha[j]*sum(STAND.POST*(matrix(rep(tau.seq*tau.seq,N*K),ncol=K,byrow = T))*(V.MATRIX[,j]%x%matrix(1,K,K))*MA)
  i.beta[j] <- par$alpha[j]*par$alpha[j]*sum(STAND.POST*(par$b[j]-matrix(rep(theta.seq,N*K),ncol=K,byrow = F))*(par$b[j]-matrix(rep(theta.seq,N*K),ncol=K,byrow = F))*(D[,j]%x%matrix(1,K,K)))+par$alpha[j]*par$alpha[j]*sum(STAND.POST*(par$b[j]-matrix(rep(theta.seq,N*K),ncol=K,byrow = F))*(par$b[j]-matrix(rep(theta.seq,N*K),ncol=K,byrow = F))*(V.MATRIX[,j]%x%matrix(1,K,K))*MA)
  
  
}
i.delta.all <- sum(i.delta)
i.beta.all <- sum(i.beta)


SE <- matrix(NA,J,6)
colnames(SE) <- c("a","b","alpha","gamma","beta","delta")
SE[,1] <- i.a
SE[,2] <- i.b
SE[,3] <- i.alpha
SE[,4] <- i.gamma
SE[1,5] <- i.beta.all
SE[1,6] <- i.delta.all
SE.RE <- 1/sqrt(SE)
ESTIMATE <- matrix(NA,J,18)
colnames(ESTIMATE) <- c("a.true","b.true","alpha.true","gamma.true","beta.true","delta.true","a.est","b.est","alpha.est","gamma.est","beta.est","delta.est","a.se","b.se","alpha.se","gamma.se","beta.se","delta.se")
ESTIMATE[,13:18] <- SE.RE
ESTIMATE[,1] <- a.true
ESTIMATE[,2] <- b.true
ESTIMATE[,3] <- alpha.true
ESTIMATE[,4] <- gamma.true
ESTIMATE[1,5] <- beta.true
ESTIMATE[1,6] <- delta.true
ESTIMATE[,7] <-  par$a
ESTIMATE[,8] <-  par$b
ESTIMATE[,9] <-  par$alpha
ESTIMATE[,10] <-  par$gamma
ESTIMATE[1,11] <-  par$beta
ESTIMATE[1,12] <- par$delta

write.csv(ESTIMATE, file = "500ESTIMATE3.csv")

