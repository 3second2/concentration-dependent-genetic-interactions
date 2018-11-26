### please change the path to the source file 
### please note that the input and output is linear values (not log transformed)
### load the loess model to predict protein dimer to the total protein concentration and vice versa
### for all the functions, paramters maximum GFP and autorluorecence needs to be provided. 
### for the paper, max GFP =3470.67, auto= 23.24.
### for the paper, total protein amount, at low expression= 5.5e-7, at high expression =8.4e-7


load("path/loess_protein_total_predic_free_dimer.RData") # please change the path
 ######## 1.     from delta_deltaG_folding to GFP signal,
      
      deltaG_deltaG_folding_to_gfp_signal<- function(deltaG_deltaG_folding, total_protein, max_gfp,auto) { 
      # defining parameters
      
      R= 1.98*10^(-3) # kcal/mol
      Temp= 310.15
      Ka= 5 *10^7 # ka for dimmer formation 
      DNA= 10^(-9)
      wt_deltaG= -2.908485
      
      # Ackers model
      relations=list(configarations=seq(1:8), deltaGs=c(0, -11.7, -10.1, -10.1, -23.8, -21.8, -22.2, -33.9) ,dimer_number=c(0,1,1,1,2,2,2,3))  
      # folding model
      fraction_folded<- exp(-(deltaG_deltaG_folding+ wt_deltaG)/(R*Temp))/(1+exp(-(deltaG_deltaG_folding+ wt_deltaG)/(R*Temp)))
      function_protein<- total_protein*fraction_folded
      
      free_dimer= predict(loess_linear, function_protein) 
      configaration<-c() 
      sum_config<- c()
      r_formd<- c() 
      for (i in 1:8) {
      configaration[i]<- exp(-relations[[2]][i]/(R*Temp))*free_dimer^relations[[3]][i]
      sum_config<- sum(configaration)
      config_prob<- configaration/sum_config 
      r_formd<- sum(config_prob[i]*relations[[3]][i])*2*DNA
      }
      Propor_sup= config_prob[2]+ config_prob[3]+ config_prob[5]+ config_prob[6]+ config_prob[7]+ config_prob[8]
      return(gfp=max_gfp-(max_gfp-auto)*Propor_sup)
      }
      
####### 2. from protein amount to GFP 
      
      protein_to_gfp_signal<- function(function_protein, max_gfp,auto) { 
      
      R= 1.98*10^(-3) # kcal/mol
      Temp= 310.15
      Ka= 5 *10^7 # ka for dimmer formation 
      DNA= 10^(-9)
      wt_deltaG= -2.908485
      
      relations=list(configarations=seq(1:8), deltaGs=c(0, -11.7, -10.1, -10.1, -23.8, -21.8, -22.2, -33.9) ,dimer_number=c(0,1,1,1,2,2,2,3))
      
      free_dimer= predict(loess_linear, function_protein) 
      configaration<-c() 
      sum_config<- c()
      r_formd<- c() 
      for (i in 1:8) {
      configaration[i]<- exp(-relations[[2]][i]/(R*Temp))*free_dimer^relations[[3]][i]
      sum_config<- sum(configaration)
      config_prob<- configaration/sum_config 
      r_formd<- sum(config_prob[i]*relations[[3]][i])*2*DNA
      }
      Propor_sup= config_prob[2]+ config_prob[3]+ config_prob[5]+ config_prob[6]+ config_prob[7]+ config_prob[8]
      return(gfp=max_gfp-(max_gfp-auto)*Propor_sup)
      }
      
##### 3. from deltaG folding to GFP signal. 
      ### this is the same as first function, in case that absolute deltaG folding is used for calculating GFP
      
      deltaG_folding_to_gfp_signal<- function(deltaG_folding, total_protein, max_gfp,auto) { 
      
      R= 1.98*10^(-3) # kcal/mol
      Temp= 310.15
      Ka= 5 *10^7 # ka for dimmer formation 
      DNA= 10^(-9)
      
      relations=list(configarations=seq(1:8), deltaGs=c(0, -11.7, -10.1, -10.1, -23.8, -21.8, -22.2, -33.9) ,dimer_number=c(0,1,1,1,2,2,2,3)) 
      fraction_folded<- exp(-deltaG_folding/(R*Temp))/(1+exp(-deltaG_folding/(R*Temp)))
      function_protein<- total_protein*fraction_folded
      
      free_dimer= predict(loess_linear, function_protein) 
      configaration<-c() 
      sum_config<- c()
      r_formd<- c() 
      for (i in 1:8) {
      configaration[i]<- exp(-relations[[2]][i]/(R*Temp))*free_dimer^relations[[3]][i]
      sum_config<- sum(configaration)
      config_prob<- configaration/sum_config 
      r_formd<- sum(config_prob[i]*relations[[3]][i])*2*DNA
      }
      Propor_sup= config_prob[2]+ config_prob[3]+ config_prob[5]+ config_prob[6]+ config_prob[7]+ config_prob[8]
      return(gfp=max_gfp-(max_gfp-auto)*Propor_sup)
      }
      
#####  Reverse processes ########
##### 4. GFP to the protein amount 
      
      
      output_to_protein<- function(auto, max_gfp, output) { # output is variant GFP, in linear scale
      suppression=1-(output-auto)/(max_gfp-auto) 
      
      R= 1.98*10^(-3) # kcal/mol
      Temp= 310.15
      Ka= 5 *10^7 # ka for dimmer formation 
      DNA= 10^(-9)
      
      relations=list(
      configarations=seq(1:8), 
      deltaGs=c(0, -11.7, -10.1, -10.1, -23.8, -21.8, -22.2, -33.9) ,
      dimer_number=c(0,1,1,1,2,2,2,3))
      
      
      Sup2<- function(free_dimer){  # free_dimer in log scale 
      k_for_each<- c() 
      for (i in 1:8) {
      k_for_each[i]<- exp(-relations[[2]][i]/(R*Temp)) 
      }
      return(1-(1+k_for_each[3]*exp(free_dimer))/(1+(k_for_each[2]+k_for_each[3]+k_for_each[4])*exp(free_dimer) + (k_for_each[5]+k_for_each[6]+k_for_each[7])*(exp(free_dimer))^2 + k_for_each[8]*(exp(free_dimer))^3)) 
      } 
      
      inverse = function (f, lower = 10^(-40), upper = 10^(-3)) { # minimum interval and muaximum interval
      function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper, extendInt = "yes")[1]
      }
      free_dimer_from_prop= inverse(Sup2, 10^(-40), 10^(-3) )
      free_dimer= exp(unlist(free_dimer_from_prop(suppression)))
      
      configaration= c()
      
      for (i in 1:8) {
      configaration[i]<- exp(-relations[[2]][i]/(R*Temp))*free_dimer^relations[[3]][i]
      sum_config<- sum(configaration)
      config_prob<- configaration/sum_config 
      r_formd<- sum(config_prob[i]*relations[[3]][i])*2*DNA
      }
      Propor_sup= config_prob[2]+ config_prob[3]+ config_prob[5]+ config_prob[6]+ config_prob[7]+ config_prob[8]
      Total_protein= (free_dimer/Ka)^0.5+ 2*free_dimer+ r_formd
      return(Total_protein)
      }
      
##### 5. total protein to fraction folded 
      fr_folded<- function(v, total_protein) { # v, is the total protein calculated from GFP for each variant
      return(v/total_protein)
      }
      
##### 6. from fraction folded to the deltaG of folding 
      fraction_folded_to_deltaG_folding<- function(fr) {
      R= 1.98*10^(-3) # kcal/mol
      Temp= 310.15 
      deltaG= -log(fr/(1-fr))*R*Temp
      return(deltaG)
      }
      
      
#### alternative submodel #######
###################################
      
########
####### 7. linear assumption of the regulatory interaction 
      ### for the paper, total protein amount, at low expression= 5.5e-7, at high expression =8.4e-7
      ### myslop = -0.52, myinter = -17
      
      
      lin_assumption_deltaG_from_gfp<- function(x, myslop, myinter,wt_gfp){ # Here, for this function, all inputs are in the log2 scale 
     
      R= 1.98*10^(-3) # kcal/mol
      Temp= 310.15
      log2_pr= (x- myinter)/myslop
      log2_wt_pr= (wt_gfp- myinter)/myslop
      fr= (2^log2_pr/2^log2_wt_pr)/114*115
      if (fr>0 & fr<1) {
      deltaG= -log(fr/(1-fr))*R*Temp
      }
      return(deltaG)
      }
      
###### reverse process of function 7, to calculate GFP from deltaG
##### 8.  from delta G to the output GFP with the regulatory interaction linear assumption 
      ### for the paper, total protein amount, at low expression= 5.5e-7, at high expression =8.4e-7
      ### myslop = -0.52, myinter = -17
      
      lin_assumption_gfp_from_deltaG<- function(deltaG, myslop, myinter,wt_protein) { # slope, intercept in log2 scale, wt_protein in linear scale
      R= 1.98*10^(-3) # kcal/mol
      Temp= 310.15
      fraction_folded<- exp(-deltaG/(R*Temp))/(1+exp(-deltaG/(R*Temp)))
      function_protein<- fraction_folded*wt_protein/114*115
      log2_gfp= myslop* log2(function_protein) + myinter
      return(log2_gfp)
      }
      