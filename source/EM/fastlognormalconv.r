library(pracma)
# Fast lognormal convolution by Mark van de Wiel
fastlognormalconv <- function(mus,intlength=0.5,maxim=100){ 
  #mus <- c(log(50),log(50),0.5,0.5)
  a<- seq(intlength/2,maxim+intlength/2,by=intlength)
  mu1 <- mus[1];mu2<-mus[2];s1<-mus[3];s2 <- mus[4]
  a1 <- dlnorm(a,mu1,s1); 
  b1 <- plnorm(a,mu2,s2) - plnorm(a-intlength,mu2,s2)
  pm <- polymul(a1,b1)
  return(intlength*cumsum(pm[1:length(a)])[-length(a)]) #only return cdf upto maxim
}

