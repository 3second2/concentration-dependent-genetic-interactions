#### model based on three nonlinear models of protein expression - ftiness effects 
### for fig5, extended data fig5
### the code includes additional plots as reference 
### any modification to the plots were performed on the Adobe illustrator 
### 'save' option is not included. 
## please save or export plots or data files by modifiying the code


####
#### make simple functions linking protein expression to the fitness

myfunction_lists<- list( # 
  inc<- function(x){ # the maximum fitness set as 1, 
    return(x/(0.1+ x)) # 0.1~ fitness=0.5;  1~ fitness=0.91 ; 10~ finess=0.99 
  }, 
  dec<- function(x){ # the maximum fitness set as 1, 
    return(1/(1+(x/10))) # 0.1~ fitness=0.99;  1~ fitness=0.91 ; 10~ finess=0.5
  }, 
  bel<-  function(x){
    return(1.2*x/(0.1+ x)*(1/(1+x/10))) # 1.2* inc*dec # multiply to 1.2 so that max is 1.  
  }# 1 ~ fitness=0.99; 0.1 and 10 fitness equals to 0.6
)
names(myfunction_lists)<- c("incre", "decre", "peaked")

####
#### define protein folding model  

fraction_folded<- function(deltaG){
  R= 1.98*10^(-3) # kcal/mol
  Temp= 310.15
  return(exp(-(deltaG)/(R*Temp))/(1+exp(-(deltaG)/(R*Temp))))
}

#####
#####
#####   generate random data, evenly spaced protein folding energies 
#####
#####
mydeltaGs<- seq(-5, 5, by=0.1) # 51, even distribution 
fraction_folded<- fraction_folded(mydeltaGs)
toy_stability<- data.frame(mydeltaGs, fraction_folded)

####
####
####  different wild  type protein stability 
####  in combination of different protein concentration
####
#### 
deltaG_wts<- c(-3, -1.6, -1, 0) #i
wts<- c(0.1, 1, 10) # k 
### a range of deltaGs - (-3, -1.6, -1, 0) 

fig_list_diff_typs<- list(
  deltG_min3= list(
    inc=list(), # l and h
    dec=list(), # l and h
    bel=list()
  ), 
  deltG_min1.6= list(
    inc=list(), # l and h
    dec=list(), # l and h
    bel=list()
  ), 
  deltG_min1= list(
    inc=list(), # l and h
    dec=list(), # l and h
    bel=list()
  ), 
  deltG_min0= list(
    inc=list(), # l and h
    dec=list(), # l and h
    bel=list()
  )
)


for (i in 1:4) { # wt deltaG 
  for (j in 1:3) { # protein expression - fitness functions
    for (k in 1:3){ # protein expression level 
      
      fig_list_diff_typs[[i]][[j]][[k]] =  expand.grid(amut=toy_stability$mydeltaGs, bmut=toy_stability$mydeltaGs)
      fig_list_diff_typs[[i]][[j]][[k]]$ab_mut<-fig_list_diff_typs[[i]][[j]][[k]]$amut + fig_list_diff_typs[[i]][[j]][[k]]$bmut - deltaG_wts[i]
      wt_fit<- myfunction_lists[[j]](fraction_folded(deltaG_wts[i])* wts[k]) 
      a_fit<- myfunction_lists[[j]](fraction_folded(fig_list_diff_typs[[i]][[j]][[k]]$amut)* wts[k]) 
      b_fit<- myfunction_lists[[j]](fraction_folded(fig_list_diff_typs[[i]][[j]][[k]]$bmut)* wts[k])
      fig_list_diff_typs[[i]][[j]][[k]]$ab_fit<- log10(myfunction_lists[[j]](fraction_folded(fig_list_diff_typs[[i]][[j]][[k]]$ab_mut)* wts[k])) 
      fig_list_diff_typs[[i]][[j]][[k]]$ab_exp_fit<- log10(a_fit* b_fit/ wt_fit)
      fig_list_diff_typs[[i]][[j]][[k]]$epis<- fig_list_diff_typs[[i]][[j]][[k]]$ab_fit - fig_list_diff_typs[[i]][[j]][[k]]$ab_exp_fit 
    }
  }
}

fig_lists_l_h<- list( 
  deltaG_min3= list(), 
  deltaG_min1.6= list() , 
  deltaG_min1= list() , 
  deltaG_min0= list() 
)

#### prepare for the plot 

for (i in 1:4) {
  for (j in 1:3) {
    fig_lists_l_h[[i]][[j]] = rbind(rbind (fig_list_diff_typs[[i]][[j]][[1]], fig_list_diff_typs[[i]][[j]][[2]]),fig_list_diff_typs[[i]][[j]][[3]] )
    fig_lists_l_h[[i]][[j]]$wt_amounts= c(rep("0.1", 2601),  rep("1", 2601),rep("10", 2601))
  }
}

### change the deltaG to the delta deltaG
for (i in 1:4) {
  for (j in 1:3) {
    fig_lists_l_h[[i]][[j]]$amut= fig_lists_l_h[[i]][[j]]$amut - deltaG_wts[i]
    fig_lists_l_h[[i]][[j]]$bmut= fig_lists_l_h[[i]][[j]]$bmut - deltaG_wts[i]
  }
}

##### subset the data to plot, because delta deltaG >5 or < -1 do not have much information and depending on the wild type protein folding energy
### the delta deltaG reanges are different. 
# only take the delta deltaG values existing in all the wildtype cases were collected to plot. 
# only two out of the three examined protein expression level cases were presented in the paper. 

library(ggplot2)
library(ggpubr)


for (i in 1:4) {
  for (j in 1:3){
    fig_lists_l_h[[i]][[j]]= fig_lists_l_h[[i]][[j]][
      fig_lists_l_h[[i]][[j]]$amut>= -1 &fig_lists_l_h[[i]][[j]]$amut<= 5 & 
        fig_lists_l_h[[i]][[j]]$bmut>= -1 &fig_lists_l_h[[i]][[j]]$bmut<= 5, ]
    names(fig_lists_l_h[[i]][[j]])[5]<- "ab_exp"
  }
}

########## plot the fitness as a function in relation to the folding energy changes 
########

increas<- ggplot(fig_lists_l_h[[2]][[1]])+ geom_line(aes(x= ab_mut, y=ab_fit, col=wt_amounts)) + 
  scale_colour_manual(values=c("blue","black", "red")) + theme_classic() + facet_grid(.~as.factor(wt_amounts))+ xlim(-5,10)
decreas<- ggplot(fig_lists_l_h[[2]][[2]])+ geom_line(aes(x= ab_mut, y=ab_fit, col=wt_amounts)) + 
  scale_colour_manual(values=c("blue","black", "red")) + theme_classic()+ facet_grid(.~as.factor(wt_amounts))+ xlim(-5,10)
bell<- ggplot(fig_lists_l_h[[2]][[3]])+ geom_line(aes(x= ab_mut, y=ab_fit, col=wt_amounts)) + 
  scale_colour_manual(values=c("blue","black", "red")) + theme_classic()+ facet_grid(.~as.factor(wt_amounts))+ xlim(-5,10)

ggarrange(increas, decreas, bell, nrow=3, ncol=1)


####
##### tile plot 
####

myepis<- list(deltaG_min3= list(), 
              deltaG_min1.6= list() , 
              deltaG_min1= list() , 
              deltaG_min0= list() )  # epistasis 

for (i in 1:4){
  for (j in 1:3) {
    myepis[[i]][[j]]<- ggplot(fig_lists_l_h[[i]][[j]]) +
      geom_tile(aes(amut, y = bmut, fill = epis)) +
      scale_fill_gradient2(low = "blue", mid = "gray88", high = "magenta") + 
      facet_grid(wt_amounts~.)+ 
      labs(x="delta_deltaG", y="delta_deltaG") + 
      geom_vline(xintercept=0, linetype=2, col="gray") + 
      geom_hline(yintercept=0, linetype=2, col="gray") + theme_bw() 
  }
}

yobs<- list(deltaG_min3= list(), 
            deltaG_min1.6= list() , 
            deltaG_min1= list() , 
            deltaG_min0= list() )  # observed fitness 

for (i in 1:4){
  for (j in 1:3) {
    myobs[[i]][[j]]<- ggplot(fig_lists_l_h[[i]][[j]]) +
      geom_tile(aes(amut, y = bmut, fill = ab_fit)) +
      scale_fill_gradient2(low="dark green", mid = "gray88", high = "magenta") + 
      facet_grid(wt_amounts~.)+ 
      labs(x="delta_deltaG", y="delta_deltaG") + 
      geom_vline(xintercept=0, linetype=2, col="gray") + 
      geom_hline(yintercept=0, linetype=2, col="gray")+ theme_bw()
  }
}

myexp<- list(deltaG_min3= list(), 
             deltaG_min1.6= list() , 
             deltaG_min1= list() , 
             deltaG_min0= list() )  # expected fitness 

for (i in 1:4){
  for (j in 1:3) {
    myexp[[i]][[j]]<- ggplot(fig_lists_l_h[[i]][[j]]) +
      geom_tile(aes(amut, y = bmut, fill = ab_exp)) +
      scale_fill_gradient2(low = "dark green", mid = "gray88", high = "magenta") + 
      facet_grid(wt_amounts~.)+ 
      labs(x="delta_deltaG", y="delta_deltaG") + 
      geom_vline(xintercept=0, linetype=2, col="gray") + 
      geom_hline(yintercept=0, linetype=2, col="gray") + theme_bw()
  }
}

####### arrange in one figure 
############################# 
ggarrange(ncol = 3, nrow = 3, 
          myepis[[1]][[1]], 
          myepis[[1]][[2]],
          myepis[[1]][[3]], 
          myobs[[1]][[1]], 
          myobs[[1]][[2]],
          myobs[[1]][[3]], 
          myexp[[1]][[1]], 
          myexp[[1]][[2]],
          myexp[[1]][[3]]) 

ggarrange(ncol = 3, nrow = 3, myepis[[2]][[1]], 
          myepis[[2]][[2]],
          myepis[[2]][[3]], 
          myobs[[2]][[1]], 
          myobs[[2]][[2]],
          myobs[[2]][[3]], 
          myexp[[2]][[1]], 
          myexp[[2]][[2]],
          myexp[[2]][[3]]) 

ggarrange(ncol = 3, nrow = 3, myepis[[3]][[1]], 
          myepis[[3]][[2]],
          myepis[[3]][[3]], 
          myobs[[3]][[1]], 
          myobs[[3]][[2]],
          myobs[[3]][[3]], 
          myexp[[3]][[1]], 
          myexp[[3]][[2]],
          myexp[[3]][[3]]) 

ggarrange(ncol = 3, nrow = 3, myepis[[3]][[1]], 
          myepis[[4]][[2]],
          myepis[[4]][[3]], 
          myobs[[4]][[1]], 
          myobs[[4]][[2]],
          myobs[[4]][[3]], 
          myexp[[4]][[1]], 
          myexp[[4]][[2]],
          myexp[[4]][[3]]) 

########
#########



