## Author: Genevieve Lefebvre
## 2018-2019

if (getRversion() >=  '3.6.0') {
  RNGkind(sample.kind="Rounding")
}

# Load some R libraries
library(copula)
library(extraDistr)
library(truncnorm)
library(truncdist)

# Generate exposure
expit<-function(lp){
  p<-exp(lp)/(exp(lp)+1)
  return(p)
}

DS<-NULL

Y1<-NULL
Y0<-NULL

# Nb monte Carlo replications (that is, nb of datasets)
rep<-100

ate<-rep(0,rep)

# Sample size
n<-1000

# Nb of covariates generated
ncov<-200

# Random covariate indices (required later)
set.seed(26769)
Ind<-sample(1:ncov,ncov)

# load coefficient regressions for treatment and outcome models (resp)
load("pcoef_scenario2.Rdata")
load("ycoef_scenario2.Rdata")

pcoef<-pcoef_scenario2
ycoef<-ycoef_scenario2

# Start looping to generate dataset and results
set.seed(46767)

for (k in 1:rep){

print(k)

DSInd<-rep(k,n)

# Gaussian copula
tau<-c(0.1,0.2,0.4,0.6,0.8)
a<- sin(pi*tau*0.5)

norm.cop1 <- normalCopula(a[1],dim = 10,dispstr = "ex")
N1a<-rCopula(n, norm.cop1)

norm.cop2 <- normalCopula(a[2],dim = 8,dispstr = "ex")
N2a<-rCopula(n, norm.cop2)

norm.cop3 <- normalCopula(a[3],dim = 6,dispstr = "ex")
N3a<-rCopula(n, norm.cop3)

norm.cop4 <- normalCopula(a[4],dim = 4,dispstr = "ex")
N4a<-rCopula(n, norm.cop4)

norm.cop5 <- normalCopula(a[5],dim = 2,dispstr = "ex")
N5a<-rCopula(n, norm.cop5)

N1_1<-qnorm(N1a[,1],1,0.5)
N1_2<-qt(N1a[,2],6,5)
N1_3<-qgamma(N1a[,3],8,1)
N1_4<-qlnorm(N1a[,4],0.3,1.2)
N1_5<-qbeta(N1a[,5],2,1)
N1_6<-qbinom(N1a[,6],1,0.05)
N1_7<-qbinom(N1a[,7],1,0.3)
N1_8<-qbinom(N1a[,8],1,0.5)
N1_9<-qbinom(1-N1a[,9],1,0.07)
N1_10<-qbinom(1-N1a[,10],1,0.6)

N1<-cbind(N1_1,N1_2,N1_3,N1_4,N1_5,N1_6,N1_7,N1_8,N1_9,N1_10)

N2_1<-qnorm(1-N2a[,1],0,1)
N2_2<-qt(N2a[,2],5,8)
N2_3<-qgamma(1-N2a[,3],2,25)
N2_4<-qgamma(N2a[,4],10,25)
N2_5<-qbinom(N2a[,5],1,0.55)
N2_6<-qbinom(N2a[,6],1,0.1)
N2_7<-qbinom(1-N2a[,7],1,0.4)
N2_8<-qbinom(1-N2a[,8],1,0.87)

N2<-cbind(N2_1,N2_2,N2_3,N2_4,N2_5,N2_6,N2_7,N2_8)

N3_1<-qnorm(N3a[,1],2,3)
N3_2<-qt(N3a[,2],3)
N3_3<-qbeta(1-N3a[,3],2,25)
N3_4<-qlnorm (N3a[,4],0.2,1)
N3_5<-qbinom(N3a[,5],1,0.2)
N3_6<-qbinom(1-N3a[,6],1,0.35)

N3<-cbind(N3_1,N3_2,N3_3,N3_4,N3_5,N3_6)

N4_1<-qnorm(N4a[,1],1,2)
N4_2<-qt(N4a[,2],10)
N4_3<-qt(N4a[,3],4,25)
N4_4<-qbinom(1-N4a[,4],1,0.35)
 
N4<-cbind(N4_1,N4_2,N4_3,N4_4)

N5_1<-qbinom(N5a[,1], 1,0.26)
N5_2<-qbinom(N5a[,2], 1,0.06)
N5<-cbind(N5_1,N5_2)

N<-cbind(N1,N2,N3,N4,N5)

# Student copula
tau2<-c(0.05,0.15,0.25,0.35,0.45)
a2<- sin(pi*tau2*0.5)

t.cop1 <- tCopula(a2[1],dim = 10,dispstr = "ex",df = 10,df.fixed = TRUE)
T1a<-rCopula(n, t.cop1)

t.cop2 <- tCopula(a2[2],dim = 8,dispstr = "ex",df = 3,df.fixed = TRUE)
T2a<-rCopula(n, t.cop2)

t.cop3 <- tCopula(a2[3],dim = 6,dispstr = "ex",df = 5,df.fixed = TRUE)
T3a<-rCopula(n, t.cop3)

t.cop4 <- tCopula(a2[4],dim = 4,dispstr = "ex",df = 8,df.fixed = TRUE)
T4a<-rCopula(n, t.cop4)

t.cop5 <- tCopula(a2[5],dim = 2,dispstr = "ex",df = 1,df.fixed = TRUE)
T5a<-rCopula(n, t.cop5)

T1_1<-qnorm(T1a[,1],2,3)
T1_2<-qt(T1a[,2],11)
T1_3<-qgamma(T1a[,3],8,2)
T1_4<-qlnorm(T1a[,4],0,1.4)
T1_5<-qbeta(1-T1a[,5],2,2)
T1_6<-qbinom(T1a[,6],1,0.06)
T1_7<-qbinom(T1a[,7],1,0.35)
T1_8<-qbinom(T1a[,8],1,0.09)
T1_9<-qbinom(1-T1a[,9],1,0.6)
T1_10<-qpois(T1a[,10],2)

T1<-cbind(T1_1,T1_2,T1_3,T1_4,T1_5,T1_6,T1_7,T1_8,T1_9,T1_10)

T2_1<-qnorm(T2a[,1],-4,1)
T2_2<-qt(1-T2a[,2],7,20)
T2_3<-qgamma(T2a[,3],3,29)
T2_4<-qgamma(1-T2a[,4],11,25)
T2_5<-qbinom(T2a[,5],1,0.55)
T2_6<-qbinom(1-T2a[,6],1,0.11)
T2_7<-qbinom(1-T2a[,7],1,0.4)
T2_8<-qbinom(T2a[,8],1,0.8)

T2<-cbind(T2_1,T2_2,T2_3,T2_4,T2_5,T2_6,T2_7,T2_8)

T3_1<-qnorm(T3a[,1],1,0.6)
T3_2<-qt(T3a[,2],3,1)
T3_3<-qbeta(T3a[,3],2,8)
T3_4<-qlnorm (1-T3a[,4],0.1,0.7)
T3_5<-qbinom(1-T3a[,5],1,0.22)
T3_6<-qbinom(T3a[,6],1,0.38)

T3<-cbind(T3_1,T3_2,T3_3,T3_4,T3_5,T3_6)

T4_1<-qunif(1-T4a[,1],0,5)
T4_2<-qinvchisq(T4a[,2],6)
T4_3<-qt(T4a[,3],9,13)
T4_4<-qzip(1-T4a[,4],1,0.7)
 
T4<-cbind(T4_1,T4_2,T4_3,T4_4)

T5_1<-qbinom(T5a[,1],1,0.36)
T5_2<-qpois(1-T5a[,2],5)

T5<-cbind(T5_1,T5_2)

T<-cbind(T1,T2,T3,T4,T5)

# Gumbel copula 
tau<-c(0.1,0.2,0.4,0.6,0.8)
thetaGumbel <- copGumbel@iTau(tau)
family <- "Gumbel" 

Gumbel.cop1<-onacopulaL(family,list(theta=thetaGumbel[1],1:10))
G1a<-rCopula(n,Gumbel.cop1)

Gumbel.cop2<-onacopulaL(family,list(theta=thetaGumbel[2],1:8))
G2a<-rCopula(n,Gumbel.cop2)

Gumbel.cop3<-onacopulaL(family,list(theta=thetaGumbel[3],1:6))
G3a<-rCopula(n,Gumbel.cop3)

Gumbel.cop4<-onacopulaL(family,list(theta=thetaGumbel[4],1:4))
G4a<-rCopula(n,Gumbel.cop4)

Gumbel.cop5<-onacopulaL(family,list(theta=thetaGumbel[5],1:2))
G5a<-rCopula(n,Gumbel.cop5)

G1_1<-qnorm(G1a[,1],-1,0.1)
G1_2<-qt(1-G1a[,2],18,-0.5)
G1_3<-qgamma(G1a[,3],1,5)
G1_4<-qlnorm(G1a[,4], 0.03,0.2)
G1_5<-qbeta(G1a[,5], 3,1.5)
G1_6<-qbinom(1-G1a[,6],1,0.90)
G1_7<-qbinom(G1a[,7],1,0.33)
G1_8<-qbinom(1-G1a[,8],1,0.66)
G1_9<-qbinom(1-G1a[,9],1,0.17)
G1_10<-qpois(G1a[,10],3)

G1<-cbind(G1_1,G1_2,G1_3,G1_4,G1_5,G1_6,G1_7,G1_8,G1_9,G1_10)

G2_1<-qtruncnorm(1-G2a[,1],11,a=9,b=Inf)
G2_2<-qt(1-G2a[,2],2.9,10)
G2_3<-qunif(G2a[,3],2.5,250)
G2_4<-qinvchisq(G2a[,4],20,25)
G2_5<-qzip(G2a[,5],5,0.55)
G2_6<-qnbinom(1-G2a[,6],7,0.3) 
G2_7<-qbinom(G2a[,7],1,0.14)
G2_8<-qbinom(1-G2a[,8],1,0.67)

G2<-cbind(G2_1,G2_2,G2_3,G2_4,G2_5,G2_6,G2_7,G2_8)

G3_1<-qnorm(1-G3a[,1],2,0.3)
G3_2<-qt(G3a[,2],7.3)
G3_3<-qgamma(1-G3a[,3],5,9)
G3_4<-qlnorm(G3a[,4],1,0.2)
G3_5<-qnbinom(G3a[,5],15,0.9)
G3_6<-qbinom(G3a[,6],1,0.02)

G3<-cbind(G3_1,G3_2,G3_3,G3_4,G3_5,G3_6)

G4_1<-qtruncnorm(G4a[,1],-10,10,a=-Inf,b=5)
G4_2<-qt(1-G4a[,2],14)
G4_3<-qinvchisq(G4a[,3],24,3)
G4_4<-qzip(1-G4a[,4],0.2,0.8)
 
G4<-cbind(G4_1,G4_2,G4_3,G4_4)

G5_1<-qbinom(G5a[,1],1,0.62)
G5_2<-qbinom(1-G5a[,2],1,0.34)

G5<-cbind(G5_1,G5_2)

G<-cbind(G1,G2,G3,G4,G5)

# Frank copula 
tau2<-c(0.05,0.15,0.25,0.35,0.45)
thetaFrank <- copFrank@iTau(tau2)
family <- "Frank" 

Frank.cop1<-onacopulaL(family,list(theta=thetaFrank[1],1:10))
F1a<-rCopula(n, Frank.cop1)

Frank.cop2<-onacopulaL(family,list(theta=thetaFrank[2],1:8))
F2a<-rCopula(n, Frank.cop2)

Frank.cop3<-onacopulaL(family,list(theta=thetaFrank[3],1:6))
F3a<-rCopula(n, Frank.cop3)

Frank.cop4<-onacopulaL(family,list(theta=thetaFrank[4],1:4))
F4a<-rCopula(n, Frank.cop4)

Frank.cop5<-onacopulaL(family,list(theta=thetaFrank[5],1:2))
F5a<-rCopula(n, Frank.cop5)

F1_1<-qnorm(F1a[,1],-10,50)
F1_2<-qt(1-F1a[,2],22,50)
F1_3<-qgamma(F1a[,3],1000,5)
F1_4<-qlnorm(F1a[,4],2,2)
F1_5<-qbeta(F1a[,5],30,0.5)
F1_6<-qbinom(F1a[,6],1,0.95)
F1_7<-qnbinom(1-F1a[,7],15,0.33)
F1_8<-qbinom(F1a[,8],1,0.46)
F1_9<-qdunif(1-F1a[,9],2,39)
F1_10<-qpois(1-F1a[,10],6)

F1<-cbind(F1_1,F1_2,F1_3,F1_4,F1_5,F1_6,F1_7,F1_8,F1_9,F1_10)

F2_1<-qtruncnorm(F2a[,1],11,30,a=9,b=Inf)
F2_2<-qtrunc(F2a[,2],"t",4,a=-20,b=29)
F2_3<-qunif(F2a[,3],-30,25)
F2_4<-qgamma(1-F2a[,4],20,250)
F2_5<-qzip(1-F2a[,5],10,0.08)
F2_6<-qnbinom(F2a[,6],7,0.5) 
F2_7<-qbinom(1-F2a[,7],1,0.44)
F2_8<-qbinom(F2a[,8],1,0.47)

F2<-cbind(F2_1,F2_2,F2_3,F2_4,F2_5,F2_6,F2_7,F2_8)

F3_1<-qinvchisq(1-F3a[,1],6.1,6.1)
F3_2<-qtrunc(1-F3a[,2],"t",1,a=-100,b=29)
F3_3<-qgamma(F3a[,3],15,9)
F3_4<-qlnorm(F3a[,4],1,0.004)
F3_5<-qnbinom(1-F3a[,5],2,0.07)
F3_6<-qbinom(F3a[,6],1,0.18)

F3<-cbind(F3_1,F3_2,F3_3,F3_4,F3_5,F3_6)

F4_1<-qtruncnorm(1-F4a[,1],0,1,a=-1,b=1)
F4_2<-qt(F4a[,2],11)
F4_3<-qinvchisq(F4a[,3],13,23)
F4_4<-qzip(1-F4a[,4],12,0.05)
 
F4<-cbind(F4_1,F4_2,F4_3,F4_4)

F5_1<-qbinom(F5a[,1],1,0.62)
F5_2<-qzip(F5a[,2],1,0.34)

F5<-cbind(F5_1,F5_2)

F<-cbind(F1,F2,F3,F4,F5)

# Joe copula 
tau3<-c(0.15,0.25,0.45,0.65,0.85)
thetaJoe <- copJoe@iTau(tau3)
family <- "Joe" 

Joe.cop1<-onacopulaL(family,list(theta=thetaJoe[1],1:10))
J1a<-rCopula(n,Joe.cop1)

Joe.cop2<-onacopulaL(family,list(theta=thetaJoe[2],1:8))
J2a<-rCopula(n,Joe.cop2)

Joe.cop3<-onacopulaL(family,list(theta=thetaJoe[3],1:6))
J3a<-rCopula(n,Joe.cop3)

Joe.cop4<-onacopulaL(family,list(theta=thetaJoe[4],1:4))
J4a<-rCopula(n,Joe.cop4)

Joe.cop5<-onacopulaL(family,list(theta=thetaJoe[5],1:2))
J5a<-rCopula(n,Joe.cop5)

J1_1<-qnorm(J1a[,1],40,9)
J1_2<-qtrunc(1-J1a[,2],"t",7,a=-30,b=4)
J1_3<-qgamma(1-J1a[,3],3,80)
J1_4<-qlnorm(J1a[,4],0.15,0.32)
J1_5<-qbeta(J1a[,5],0.5,100)
J1_6<-qbinom(1-J1a[,6],1,0.75)
J1_7<-qnbinom(1-J1a[,7],15,0.83)
J1_8<-qbinom(1-J1a[,8],1,0.12)
J1_9<-qdunif(J1a[,9],-12,59)
J1_10<-qpois(J1a[,10],16)

J1<-cbind(J1_1,J1_2,J1_3,J1_4,J1_5,J1_6,J1_7,J1_8,J1_9,J1_10)

J2_1<-qtruncnorm(J2a[,1],7,50,a=1,b=70)
J2_2<-qtrunc(J2a[,2],"t",3,a=-5,b=45)
J2_3<-qunif(1-J2a[,3],-38,-15)
J2_4<-qgamma(J2a[,4],200,250)
J2_5<-qzip(1-J2a[,5],15,0.28)
J2_6<-qnbinom(J2a[,6],17,0.75) 
J2_7<-qbinom(1-J2a[,7],1,0.09)
J2_8<-qbinom(1-J2a[,8],1,0.91)

J2<-cbind(J2_1,J2_2,J2_3,J2_4,J2_5,J2_6,J2_7,J2_8)

J3_1<-qinvchisq(1-J3a[,1],40.1,6.1)
J3_2<-qtrunc(J3a[,2],"t",2,a=-10,b=6)
J3_3<-qgamma(J3a[,3],1500,0.9)
J3_4<-qlnorm(1-J3a[,4],0.1,0.006)
J3_5<-qnbinom(J3a[,5],2,0.5)
J3_6<-qbinom(J3a[,6],1,0.52)

J3<-cbind(J3_1,J3_2,J3_3,J3_4,J3_5,J3_6)

J4_1<-qtruncnorm(J4a[,1],0,1,a=-0.78,b=0.99)
J4_2<-qt(1-J4a[,2],21)
J4_3<-qinvchisq(J4a[,3],60,23)
J4_4<-qzip(1-J4a[,4],5,0.15)
 
J4<-cbind(J4_1,J4_2,J4_3,J4_4)

J5_1<-qbinom(J5a[,1],1,0.5)
J5_2<-qnbinom(J5a[,2],2,0.5)

J5<-cbind(J5_1,J5_2)

J<-cbind(J1,J2,J3,J4,J5)

# Clayton copula 
thetaClayton <- copClayton@iTau(tau2)
family <- "Clayton" 

Clayton.cop1<-onacopulaL(family,list(theta=thetaClayton[1],1:10))
C1a<-rCopula(n,Clayton.cop1)

Clayton.cop2<-onacopulaL(family,list(theta=thetaClayton[2],1:8))
C2a<-rCopula(n,Clayton.cop2)

Clayton.cop3<-onacopulaL(family,list(theta=thetaClayton[3],1:6))
C3a<-rCopula(n,Clayton.cop3)

Clayton.cop4<-onacopulaL(family,list(theta=thetaClayton[4],1:4))
C4a<-rCopula(n,Clayton.cop4)

Clayton.cop5<-onacopulaL(family,list(theta=thetaClayton[5],1:2))
C5a<-rCopula(n,Clayton.cop5)

C1_1<-qnorm(C1a[,1],-40,19)
C1_2<-qtrunc(C1a[,2],"t",1,20,a=-30,b=4)
C1_3<-qgamma(1-C1a[,3],183,80)
C1_4<-qlnorm(1-C1a[,4], 0.25,0.2)
C1_5<-qbeta(C1a[,5], 500,10)
C1_6<-qbinom(C1a[,6],1,0.71)
C1_7<-qnbinom(1-C1a[,7],5,0.8)
C1_8<-qbinom(C1a[,8],1,0.21)
C1_9<-qdunif(C1a[,9],-112,-59)
C1_10<-qpois(C1a[,10],9)

C1<-cbind(C1_1,C1_2,C1_3,C1_4,C1_5,C1_6,C1_7,C1_8,C1_9,C1_10)

C2_1<-qtruncnorm(C2a[,1],17,5,a=7,b=70)
C2_2<-qtrunc(1-C2a[,2],"t",18,a=-5,b=2)
C2_3<-qunif(C2a[,3], 38,85)
C2_4<-qgamma(1-C2a[,4],2,1250)
C2_5<-qzip(C2a[,5],7,0.05)
C2_6<-qnbinom(C2a[,6],170,0.15) 
C2_7<-qbinom(C2a[,7],1,0.29)
C2_8<-qbinom(1-C2a[,8],1,0.43)

C2<-cbind(C2_1,C2_2,C2_3,C2_4,C2_5,C2_6,C2_7,C2_8)

C3_1<-qinvchisq(1-C3a[,1],4.1,10.1)
C3_2<-qtrunc(1-C3a[,2],"t",2,a=-20,b=16)
C3_3<-qgamma(C3a[,3],0.4,0.9)
C3_4<-qlnorm(1-C3a[,4],0.18,0.06)
C3_5<-qnbinom(C3a[,5],12,0.45)
C3_6<-qbinom(C3a[,6],1,0.62)

C3<-cbind(C3_1,C3_2,C3_3,C3_4,C3_5,C3_6)

C4_1<-qtruncnorm(C4a[,1],0,0.41,a=-0.5,b=0.5)
C4_2<-qt(C4a[,2],20)
C4_3<-qinvchisq(1-C4a[,3],40,24)
C4_4<-qzip(1-C4a[,4],6,0.19)
 
C4<-cbind(C4_1,C4_2,C4_3,C4_4)

C5_1<-qbinom(C5a[,1],1,0.15)
C5_2<-qnbinom(1-C5a[,2],2,0.54)

C5<-cbind(C5_1,C5_2)

C<-cbind(C1,C2,C3,C4,C5)

# Independent copula (20 covariables)

Ia<- matrix(runif(n*20),ncol=20)
I_1<-qnorm(Ia[,1],20,4)
I_2<-qt(Ia[,2],13,ncp=-4)
I_3<-qgamma(Ia[,3],1,8)
I_4<-qlnorm(Ia[,4], 1.5,1.4)
I_5<-qbeta(Ia[,5],62,1)
I_6<-qtrunc(Ia[,6],"t",4,a=-4,b=45)
I_7<-qgamma(Ia[,7],80,2000)
I_8<-qlnorm(Ia[,8], 0.15,0.05)
I_9<-qbeta(Ia[,9], 4,10)
I_10<- qtruncnorm(Ia[,10],8,30,a=-49,b=170)

I_11<-qbinom(Ia[,11],1,0.07)
I_12<-qbinom(Ia[,12],1,0.14)
I_13<-qnbinom(Ia[,13],3,0.29)
I_14<-qbinom(Ia[,14],1,0.21)
I_15<-qpois(Ia[,15],4)
I_16<-qbinom(Ia[,16],1,0.28)
I_17<-qbinom(Ia[,17],1,0.35)
I_18<-qbinom(Ia[,18],1,0.42)
I_19<-qzip(Ia[,19],2,0.3)
I_20<-qpois(Ia[,20],27)

I<-cbind(I_1,I_2,I_3,I_4,I_5,I_6,I_7,I_8,I_9,I_10,I_11,I_12,I_13,I_14,I_15,I_16,I_17,I_18,I_19,I_20)

# Combine data
Va<-cbind(N,T,G,F,J,C,I)

# Shuffle columns according to covariates' indices generated previously (before starting loop)
V<-data.frame(Va[,Ind])

for (i in 1:ncol(V)) {
   colnames(V)[i] <- paste("V", i, sep="")
}

### Generate treatment

V1sqrt<-sqrt(V[,1])
sdV1sqrt<-4*sd(V1sqrt)

sdV5<-sd(V[,5])
sdV32<-sd(V[,32])
intV5V32<-V[,5]*V[,32]
sdintV5V32<-4*sd(intV5V32)

sdV70<-sd(V[,70])
V70square<-V[,70]^2
sdV70square<-4*sd(V70square)

IndV101<-V[,101]>2.5

IndcontV179<-(V[,179]+1.5)*(V[,179]< -0.5)
sdIndcontV179<-4*sd(IndcontV179)

Vtreatment<-cbind(1,V1sqrt,V[,5],V[,32],intV5V32,V[,70],V70square,IndV101,V[,150],IndcontV179)

probexpint<- apply(t(t(Vtreatment)*pcoef),1,sum)
expitp<-expit(probexpint)
A<-rbinom(n,1,expitp)

### Generate Y

VYa<-cbind(V1sqrt,V[,5],V[,23],V[,32],intV5V32,V[,70],V70square,IndV101,V[,150],IndcontV179,V[,199])
VY<-data.frame(VYa)

meany<- apply(t(t(VY)*ycoef),1,sum)+0.8*A

err<-rnorm(n,sd=1)

Y<-meany+err 

meany11<- apply(t(t(VY)*ycoef),1,sum)+0.8*1
meany10<- apply(t(t(VY)*ycoef),1,sum)+0.8*0

Y11<-meany11+err
Y10<-meany10+err 

atea<-Y11-Y10
ate[k]<-mean(atea)

datatest<-cbind(Y,A,V)

CS<-data.frame(datatest)

# Create CSV files
csvFileName <- paste("CHDScenario2DS",k,".csv",sep="")
write.csv(CS, file=csvFileName,row.names=FALSE) 

DSa<-data.frame(DSInd,CS)

DS<-rbind(DS,DSa)

Y1<-c(Y1,Y11)

Y0<-c(Y0,Y10)

} # close loop

DS_scenario2<-DS
save(DS_scenario2,file="DS_scenario2.Rdata")

ate_scenario2<-ate
save(ate_scenario2,file="ate_scenario2.Rdata")

# > summary(ate)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.8     0.8     0.8     0.8     0.8     0.8 