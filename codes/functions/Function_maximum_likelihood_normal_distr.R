#### maximum likelyhood in the biological relevant range for each rep 

## mu = mean values to rescale 
## sigma, standard deviation of the mean value to rescale 
## mean_to_estimate: kth bin value within the range. 

maximum_likelihood<- function(mu, sigma, mean_to_estimate){
  probability<- (sigma*(2*pi)^0.5)^-1* exp(-0.5*((mu-mean_to_estimate)/sigma)^2)
  likelihood<- sum(probability* mean_to_estimate)/sum(probability)
  return(likelihood)
} 
