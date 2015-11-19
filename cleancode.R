if(!("moments" %in% installed.packages()[,1])){install.packages("moments")}
if(!("CCA" %in% installed.packages()[,1])){install.packages("CCA")}
if(!("ggplot2" %in% installed.packages()[,1])){install.packages("ggplot2")}
library(moments)
library(CCA)
library(ggplot2)

setwd("C:/Users/MACIEK/Desktop/analizaWW/AWWPD1")
zdjecia<-read.csv("zdjecia.csv",sep=";",header = TRUE,colClasses=c(rep("numeric",11)),row.names=4)
zdjecia0<-zdjecia
skalarad<-c("b.slabo","slabo","przecietnie", "dobrze", "b.dobrze","rewelacyjnie")
#loop across 3-factor var. with the same levels
factors<-sapply(c(5,8,9),function(x) zdjecia[,x]<<-factor(zdjecia[,x],labels=skalarad))
# factor - kolor
zdjecia$kolor<-factor((zdjecia$kolor),labels=c("fioletowy","zielony", "niebieski","inny"))
#factor - miejsce
zdjecia$miejsce<-factor((zdjecia$miejsce),labels=c("Polska","inny kraj"))
#factor - odbicia
zdjecia$odbicia<-factor((zdjecia$odbicia),labels=c("brak","wys. odbicia"))



zdjecia02<-zdjecia0
sap1<-sapply(c("miejsce","odbicia","kolor"),function(x) zdjecia02[,x]<<-as.factor(zdjecia02[,x]))
for(i in c("miejsce","odbicia","kolor")){
  sap2<-sapply(2:length(levels(zdjecia02[,i])),function(x) zdjecia02[paste0(i,x)]<<-as.numeric(zdjecia02[,i]==levels(zdjecia02[,i])[x]))
}
zdjecia02<-zdjecia02[,!names(zdjecia02) %in% c("miejsce","odbicia","kolor")]


head(zdjecia0,3)
head(zdjecia,3)
head(zdjecia02,3)


komPolska<-zdjecia[which(zdjecia$miejsce=="Polska"),4]
kominne<-zdjecia[which(zdjecia$miejsce=="inny kraj"),4]
t.test(komPolska,kominne)



zdjecianie<-zdjecia[which(zdjecia[,"kolor"]=="niebieski"),]
zdjecianieosoby<-zdjecianie$osoby
zdjecianiekoty<-zdjecianie$koty
t.test(zdjecianieosoby,zdjecianiekoty)



tableao<-as.matrix(table((zdjecia$artyzm),(zdjecia$ostrosc),deparse.level=2))
tableao
sumr<-apply(tableao,1,sum)
sumc<-apply(tableao,2,sum)
tableao<-rbind(tableao,sumc)
tableao<-cbind(tableao,sumr=c(sumr,sum(sumr)))
sapply(1:6,function(x) tableao[,x]/tableao[7,x])
#correalation and significance of it
cor.test(as.numeric(zdjecia0$artyzm),as.numeric(zdjecia0$ostrosc),method = c("spearman"))


set.seed(308914)
randomrows<-sort(sample(1:nrow(zdjecia0), 781, replace = FALSE, prob = NULL))
zdjecia.random<-zdjecia0[randomrows,]
summary(zdjecia.random)


kurtosis<-sapply(1:ncol(zdjecia.random),function(x) kurtosis(zdjecia.random[,x]))
skewness<-sapply(1:ncol(zdjecia.random),function(x) skewness(zdjecia.random[,x]))
sd<-sapply(1:ncol(zdjecia.random),function(x) sd(zdjecia.random[,x]))
mean<-sapply(1:ncol(zdjecia.random),function(x) mean(zdjecia.random[,x]))
bound<-mean +3*sd
jarque.stat<-sapply(1:ncol(zdjecia.random),function(x) jarque.test(zdjecia.random[,x])$statistic)
stats<-matrix(1:(ncol(zdjecia0)*6),ncol=ncol(zdjecia0),nrow=6)
stats<-rbind(kurtosis,skewness,jarque.stat,sd,mean,bound)
colnames(stats)<- colnames(zdjecia0)
stats


par(mfrow=c(3,4), mar=c(4,4,2,1), oma=rep(2,4))
for(i in 1:ncol(zdjecia0)){
  hist(zdjecia0[,i],xlab = colnames(zdjecia0)[i],main = "")
}



round(cor(zdjecia0, method = c("spearman")),2)



round(cor(zdjecia0, method = c("pearson")),2)



cor.test(as.numeric(zdjecia0$artyzm),as.numeric(zdjecia0$komentarze),method = c("pearson"))
cor.test(as.numeric(zdjecia0$ocena),as.numeric(zdjecia0$dzien),method = c("spearman"))



zbior1<-c("artyzm", "ostrosc", "ocena")
zbior2<-c("kolor2","kolor3","kolor4","koty", "miejsce2", "odbicia2", "dzien", "osoby")
matY<-zdjecia02[,zbior1]
matX<-zdjecia02[,zbior2]
corrall<-matcor(matX,matY)
corrall$XYcor



cc1 <- cc(matX,matY)
# sklad oraz rodzaje wynikow dla funkcji cc
str(cc1)



cc1$cor
cc1$ycoef
cc1$xcoef



wyniki<-list()
wyniki[[1]] <- diag(sqrt(diag(cov(matY)))) %*% cc1$ycoef
rownames(wyniki[[1]])<-rownames(cc1$ycoef)
wyniki[[2]] <- diag(sqrt(diag(cov(matX)))) %*% cc1$xcoef
rownames(wyniki[[2]])<-rownames(cc1$xcoef)
wyniki



WILKSL<-function(matX,matY,cc1){
  ev <- (1 - cc1$cor^2)
  n <- dim(matX)[1]
  p <- length(matX)
  q <- length(matY)
  k <- min(p, q)
  m <- n - 3/2 - (p + q)/2
  w <- rev(cumprod(rev(ev)))
  # initialize
  d1 <- d2 <- f <- vector("numeric", k)
  
  for (i in 1:k) {
    s <- sqrt((p^2 * q^2 - 4)/(p^2 + q^2 - 5))
    si <- 1/s
    d1[i] <- p * q
    d2[i] <- m * s - p * q/2 + 1
    r <- (1 - w[i]^si)/w[i]^si
    f[i] <- r * d2[i]/d1[i]
    p <- p - 1
    q <- q - 1
  }
  pv <- pf(f, d1, d2, lower.tail = FALSE)
  dmat <- cbind(WilksL = w, F = f, df1 = d1, df2 = d2, p = pv)
  return(dmat)
  ## source: http://www.ats.ucla.edu/stat/r/dae/canonical.htm
}


WILKS2<-function(cc1){
  dmat2<-matrix(0,nrow=ncol(matY),ncol=2)
  sapply(1:ncol(matY),function(i) dmat2[i,1]<<-prod((1-cc1$cor^2)[i:(ncol(matY))]))
  return(dmat2)
}

WILKSL(matX,matY,cc1)
WILKS2(cc1)



REDUNT<-function(matX,matY,cc1){
  eigenmatY<-cc1$cor
  vector1<-vector(,length(matY))
  sapply(1:ncol(matY),function(i) vector1[i]<<-sqrt(eigenmatY[i]))
  names1<-c("opposite variance","own variance")
  names2<-c("prop stdvar v","prop stdvar u")
  matim<<-list()
  for(i in (1:ncol(matY))){
    a<-round(sum((cc1$scores$corr.Y.xscores[,i])^2)/ncol(matY),3)
    b<-round(a/vector1[i],3)
    c<-round(sum((cc1$scores$corr.X.yscores[,i])^2)/ncol(matX),3)
    d<-round(c/vector1[i],3)
    assign(paste0("amat",i),matrix(c(c,d,a,b),byrow=TRUE,nrow=2,ncol=2,dimnames=list(names2,names1)),inherit=TRUE)
    matim[[i]]<<-get(paste0("amat",i))
  }
  return(matim)
}

REDUNT(matX,matY,cc1)

par(mfrow=c(3,4), mar=c(4,4,2,1), oma=rep(2,4))
for(i in 1:ncol(zdjecia0)){
  boxplot(zdjecia0[,i],xlab = colnames(zdjecia0)[i],main = "")
}

outliers<-list()
for(i in 1:ncol(zdjecia0)){
  tt<-zdjecia0[which(zdjecia0[,i]>stats[6,i]),]
  outliers[[i]]<-tt
}
names(outliers)<-c(colnames(zdjecia0))

outliers[1:4]

