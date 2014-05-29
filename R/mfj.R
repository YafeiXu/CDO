mfj <-
function(theta1Input,theta2Input,weightInput,M,dateInput){
#vMixedTrancheGsT=function(theta1Input,theta2Input,weightInput,M,dateInput)
#theta1Input=.15
#theta2Input=.08
#weightInput=.03
#MInput=100
#dateInput=c("2007/10/23")
#---define variables for using later
#M=MInput # Monte Carlo runsssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
d=125 # # of entities
#rho=rhoInput # correlation parameter in copulaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
U=matrix(NA,M,d)
norm.cop=list()
tao=matrix(NA,d,M)
# Settlement Date  
CompDate=list()
CompDate=as.Date(dateInput) # 设定计算日Settlement Dateeeeeeeeeeeeeeeeeeeeeeeeeeeee
#lambda=lambdaInput
lambda=numeric()

lambdaVector=read.csv("C:/defIntensity.csv")
lambda=lambdaVector[which(as.Date(lambdaVector[,1])==as.Date(CompDate)),3]

diffTime=numeric()
timeTable=data.frame()
Time=numeric()
TT=numeric()
BETA=numeric()
payDay=read.csv("C:/payday.csv")[,1]
BETA=numeric()

# Flat Recovery Rate		40% 设定损失回收率R
R=numeric()
R=0.4#recoveryInput
# Risk-free interest rate		5.0% 设定无风险回报率
IntRate=numeric()
IntRate=0.03#intRateInput









################################################################################
#---check needed packages 检查软件包是否安装
is.installed=list()
is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1])

if(is.installed("copula")){library(copula)
}else{
install.packages("copula")
library(copula)}

if(is.installed("matrixcalc")){library(matrixcalc)
}else{
install.packages("matrixcalc")
library(matrixcalc)}


#---1) 2)Compute tao 计算违约时间  copulaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa


frank.cop1 <- frankCopula(iTau(frankCopula(), theta1Input), dim = d)


gumbel.cop <- gumbelCopula(iTau(gumbelCopula(), theta2Input), dim = d)

rN=ifelse((M*weightInput)<1,1,floor(M*weightInput))

U=t(rbind(rCopula(rN, frank.cop1 ),rCopula(M-rN,gumbel.cop))[sample(M,M,replace=FALSE),])


tao=(-1/lambda)*log(U)

#---3) 4)Compute Time 求Time(term structure)
diffTime=as.numeric(as.Date(payDay)-as.Date(CompDate))
timeTable=data.frame(as.Date(payDay),diffTime)
if(timeTable[1,2]<0){timeTable=timeTable[-which(timeTable[,2]<0),]}else{timeTable=timeTable}

Time=timeTable[,1] # term structure 未来实际的payDay
 
TT=length(timeTable[,1]) #  
tPercent=numeric()
tPercent=cumsum(as.numeric(as.Date(c(CompDate,Time))[-1]-as.Date(c(CompDate,Time))[-length(c(CompDate,Time))])/365)
ONE=matrix(NA,125*M,TT)



#---6) Compute ONE matrix  求ONE矩阵
taoVector=matrix(NA,125*M,1)
taoVector=vec(tao)
ONE=apply(taoVector,1,function(x)as.numeric(x<=tPercent))
ONE=t(ONE)

#---5) Compute BETA 求beta 折现因子
BETA=exp(-IntRate*(tPercent))

#---6) L matrix 求Lt矩阵
L=matrix(NA,M,TT)
for(i in 1:M){
   for(j in 1:TT){
      L[i,j]=sum(ONE[(i*125-124):(i*125),j])
   }
}
L=((1-R)/d)*L

#---7) Lj matrix 求Lj矩阵
l1=numeric()
u1=numeric()
l2=numeric()
u2=numeric()
l3=numeric()
u3=numeric()
l4=numeric()
u4=numeric()
l5=numeric()
u5=numeric()
l1=0
u1=0.03
l2=0.03
u2=0.06
l3=0.06
u3=0.09
l4=0.09
u4=0.12
l5=0.12
u5=0.22
fct=function(x){

L1=numeric()
L2=numeric()
L3=numeric()
L4=numeric()
L5=numeric()

L1=0.03*as.numeric(x>u1)+(x-l1)*as.numeric(x>l1 & x<=u1)+0*as.numeric(x<=l1)
L2=0.03*as.numeric(x>u2)+(x-l2)*as.numeric(x>l2 & x<=u2)+0*as.numeric(x<=l2)
L3=0.03*as.numeric(x>u3)+(x-l3)*as.numeric(x>l3 & x<=u3)+0*as.numeric(x<=l3)
L4=0.03*as.numeric(x>u4)+(x-l4)*as.numeric(x>l4 & x<=u4)+0*as.numeric(x<=l4)
L5=0.1*as.numeric(x>u5)+(x-l5)*as.numeric(x>l5 & x<=u5)+0*as.numeric(x<=l5)
return(c(L1,L2,L3,L4,L5))
}

LjResult=matrix(NA,5,M*TT)
LjResult=sapply(t(L),fct)
LjResult=t(LjResult)

#---8) LjResultMinus1
minusIndex=numeric()
LjResultTemp=matrix(NA,M*TT,5)
minusIndex=c(1:M)*TT
LjResultTemp=LjResult
LjResultTemp[minusIndex,]=0
LjResultMinus1=LjResultTemp
LjResultMinus1=LjResultMinus1[-length(LjResultMinus1[,1]),]
LjResultMinus1=rbind(seq(0,0,length.out=5),LjResultMinus1)

#---9) LjDiff
LjDiff=matrix(NA,M*TT,5)
LjDiff=LjResult-LjResultMinus1

#---10) BETAMatrix
betaMatrix=matrix(NA,M*TT,5)
betaMatrix=matrix(BETA,M*TT,5)

#---11) 12)upStar, upSumStar, eUpSumStar
upStar=matrix(NA,M*TT,5)
upSumStar=numeric()
eUpSumStar=numeric()
upStar=betaMatrix*LjDiff
upSumStar=colSums(upStar)
eUpSumStar=upSumStar*(1/M)

#---13) FjResult
FjResult=matrix(NA,M*TT,5)
FjResult=cbind((u1-l1)-LjResult[,1],(u2-l2)-LjResult[,2],(u3-l3)-LjResult[,3],(u4-l4)-LjResult[,4],(u5-l5)-LjResult[,5])

#---14) FjResultMinus1
FjResultMinus1=matrix(NA,M*TT,5)
FjResultTemp=matrix(NA,M*TT,5)
FjResultTemp=FjResult
FjResultTemp[minusIndex,]=0
FjResultMinus1=FjResultTemp[-length(FjResult[,1]),]
FjResultMinus1=rbind(seq(0,0,length.out=5),FjResultMinus1)

#---15) FjAdd
FjAdd=matrix(NA,M*TT,5)
FjAddHalf=matrix(NA,M*TT,5)
FjAdd=FjResult+FjResultMinus1
FjAddHalf=0.5*FjAdd

#---16) deltaMatrix
deltaT=numeric()
deltaT=as.numeric(as.Date(c(CompDate,Time))[-1]-as.Date(c(CompDate,Time))[-length(c(CompDate,Time))])/365
deltaMatrix=matrix(NA,M*TT,5)
deltaMatrix=matrix(deltaT,M*TT,5)

#---17) lowStar, lowSumStar, eLowSumStar
lowStar=matrix(NA,M*TT,5)
lowSumStar=numeric()
eLowSumStar=numeric()
lowStar=FjAddHalf*betaMatrix*deltaMatrix
lowSumStar= colSums(lowStar)
eLowSumStar=(1/M)*lowSumStar

#---18) tranches' spreads
trancheSpread1=(eUpSumStar[1]-0.05*eLowSumStar[1])/0.03
trancheSpread2=eUpSumStar[2]/eLowSumStar[2]
trancheSpread3=eUpSumStar[3]/eLowSumStar[3]
trancheSpread4=eUpSumStar[4]/eLowSumStar[4]
trancheSpread5=eUpSumStar[5]/eLowSumStar[5]

#---19) return value
return(c(trancheSpread1,trancheSpread2,trancheSpread3,trancheSpread4,trancheSpread5))

}
