##### weighted mean function 
# v, mean 
# s, error of the mean

myweighted_means_err<- function(v, s) {
  vv<- sum(v*(s^(-2)))
  ss<-  sum(s^(-2)) 
  return( c(mean=vv/ss, err=sum((v- vv/ss)^2)/(length(v)- 1) ))
  }

