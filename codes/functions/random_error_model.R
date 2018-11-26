#### change the i to different iteration numbers , default here is 50 

fisher_mean<- function (rep1, rep2, rep3, se1,se2,se3) {
  unweigthed_mean= (rep1+rep2+rep3)/3 
  se0= 0.5* ((rep1-unweigthed_mean)^2 + (rep2-unweigthed_mean)^2 +(rep3-unweigthed_mean)^2)
  weighted_mean= (rep1*(se0+se1^2)^(-1)+ rep2*(se0+se2^2)^(-1) +rep3*(se0+se3^2)^(-1) )/((se0+se1^2)^(-1)+ (se0+se2^2)^(-1) +(se0+se3^2)^(-1))
  se<- c()
  se[1]<- se0
  weighted_means<- c()
  weighted_means[1]<- weighted_mean
  for (i in 2:50) {
    se[i]= se[i-1]*((se[i-1]+se1^2)^(-2)*(rep1-weighted_means[i-1])^2 + (se[i-1]+se2^2)^(-2)*(rep2-weighted_means[i-1])^2 +  (se[i-1]+se3^2)^(-2)*(rep3-weighted_means[i-1])^2)/((se[i-1]+se1^2)^(-1) + (se[i-1]+se2^2)^(-1) +  (se[i-1]+se3^2)^(-1) - ((se[i-1]+se1^2)^(-2) + (se[i-1]+se2^2)^(-2) +  (se[i-1]+se3^2)^(-2))/((se[i-1]+se1^2)^(-1) + (se[i-1]+se2^2)^(-1) +  (se[i-1]+se3^2)^(-1)))
    weighted_means[i]= (rep1*(se[i]+se1^2)^(-1)+ rep2*(se[i]+se2^2)^(-1) +rep3*(se[i]+se3^2)^(-1))/((se[i]+se1^2)^(-1) + (se[i]+se2^2)^(-1) +  (se[i]+se3^2)^(-1))
  }
  return(c(weighted_mean= weighted_means[50], 
              se_final= 1/ ((se[50]+ se1^2)^(-1) + (se[50]+ se2^2)^(-1) + (se[50]+ se3^2)^(-1))^0.5))
}

#### test the cycle of iteration 
