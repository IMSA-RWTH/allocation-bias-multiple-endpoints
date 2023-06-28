#' Function to calculate the actual mean FWER and T1E and the probabilities of obtaining an actual FWER and T1E that not exceeds the  significance level [P(FWER<=0.05)and P(T1E<=0.05)] regarding different randomization procedures in clinical trials that are evaluated with the Sidak and all-or-none procedure, respectivel.

library(randomizeR)


#' @param N total sample size
#' @param m number of endpoints
#' @param randomization one of the randomization procedures "RAR", "PBR", "CHEN", "CR", "BSD", "MP", "EBC", "UD" (default parameter for this randomization procedures: PBR(4), BSD(3), EBC(0.67), CHEN(2,0.67), MP(3))
#' @param r_number number of randomization lists that be used to calculate mean FWER and T1E and the probabilities of obtaining an actual FWER and T1E that exceeds the 5% significance level
#' @return object of the class 'randSeq' of the package randomizeR
random_list_samp<-function(N,m,randomization,r_number,CHEN_par1=2,CHEN_par2=0.67,BSD_par=3,EBC_par=0.67,MP_par=3,PBR_par=4){

  if( randomization=='CR'){
    Rand<-crPar(N ,K = 2, ratio = rep(1, 2), groups = LETTERS[1:2])
    myPar<- genSeq(Rand,r_number)
    return(myPar)
  }else if(randomization=='BSD'){
    Rand<- bsdPar(N, BSD_par, groups = LETTERS[1:2])
    myPar<- genSeq(Rand,r_number)
    return(myPar)
  }else if(randomization=='CHEN'){  
    Rand<- chenPar(N,CHEN_par1,CHEN_par2,groups = LETTERS[1:2])
    myPar<- genSeq(Rand,r_number)
    return(myPar)
  }else if(randomization=='EBC'){
    Rand<-ebcPar(N, EBC_par, groups = LETTERS[1:2])
    myPar<- genSeq(Rand,r_number)
    return(myPar)
  }else if(randomization=='MP'){
    Rand<-mpPar(N ,mti=MP_par, ratio = rep(1, 2), groups = LETTERS[1:2])
    myPar<- genSeq(Rand,r_number)
    return(myPar)
  }else if(randomization=='PBR'){
    if(N %% PBR_par==0){
      Rand<-pbrPar(rep(PBR_par, times= N/PBR_par), K = 2, groups = LETTERS[1:2])
      myPar<- genSeq(Rand,r_number)
      return(myPar)
    }else{stop("N is not divisible by PBR_par")}
  }else if(randomization=='RAR'){
    Rand<-rarPar(N,K=2,groups = LETTERS[1:2])
    myPar<- genSeq(Rand,r_number)
    return(myPar)
  }else{stop("randomization procedure is not defined")}
  
}

#'@param Rand_list matrix that rows contains randomization lists with entries "A" (treatment group) and "B" (control group)
#'@return vector, with the numbers of allocations to the control group

get_numb_contr<-function(Rand_list){
  NC<-rep(NULL,dim(Rand_list)[1])
  
  for(j in 1:dim(Rand_list)[1]){
    NC[j]=length(which(Rand_list[j,]=="B"))
  }
  return(NC)
}

#'@param Rand_list matrix that rows contains randomization lists with entries "A" (treatment group) and "B" (control group)
#'@return vector, with the numbers of allocations to the treatment group

get_numb_treat<-function(Rand_list){
  NE<-rep(NULL,dim(Rand_list)[1])
  
  for(j in 1:dim(Rand_list)[1]){
    NE[j]=length(which(Rand_list[j,]=="A"))
  }
  return(NE)
}

#' @param m number of endpoints
#' @param sigma standard deviation of the outcome of the different endpoints, vector of size m
#' @param Rand_list matrix that rows contains randomization lists with entries "A" (treatment group) and "B" (control group)
#' @return a named list, with transformed variances named val and the transformation matrix of the patients named vec.


get_transf<-function(m,sigma,Rand_list, p=NULL){
  
  T<-diag(sigma)
  R<-matrix(p,m,m)+diag(rep(1-p,m))
  S<- T%*%R%*%T
  ev<-eigen(x=S)
  
  return(ev)
}

#'@param RandSeq list with the information of the randomization of the patients [entries "A" (treatment group) and "B" (control group)]
#'@param sigma standard deviation of the outcome of the different endpoints, vector of size m
#'@param NE number of allocations to the treatment group
#'@param NC number of allocation to the control group
#'@param tar vector of size m that contains the allocation bias effects for each patient 

get_noncentralparam<-function(RandSeq,sigma,NE,NC,tau){
  N=NE+NC
  if(length(which(RandSeq=='A'))==(NE+NC)){
    
    delta=0
    lambda= 1/sigma^2*(sum(tau^2)-N*mean(tau[RandSeq=="A"])^2)
    
  }else if(length(which(RandSeq=='B'))==(NE+NC)){
    delta=0
    lambda= 1/sigma^2*(sum(tau^2) - N*mean(tau[RandSeq=="B"])^2 )
    
  }else{
    
    delta=1/sigma*sqrt((NE*NC)/N)*(mean(tau[RandSeq=="A"]) - mean(tau[RandSeq=="B"]))
    lambda= 1/sigma^2*(sum(tau^2)-NE*mean(tau[RandSeq=="A"])^2 - NC*mean(tau[RandSeq=="B"])^2 )
  }
  
  return(data.frame(delta=delta,lambda=lambda))
}



#' Explanation of the inputs of @function multEndp:
#' @param N total sample size
#' @param m number of endpoints
#' @param sigma standard deviation of the outcome of the different endpoints, vector of size m
#' @param randomization one of the randomization procedures "RAR", "PBR", "CHEN", "CR", "BSD", "MP", "EBC", "UD" (default parameter for this randomization procedures: PBR(4), BSD(3), EBC(0.67), CHEN(2,0.67), MP(3))
#' @param eta endpoint-specific allocation bias effects, vector of the size m
#' @param r_number number of randomization lists that be used to calculate mean FWER and T1E and the probabilities of obtaining an actual FWER and T1E that exceeds the 5% significance level
#' @param procedure one of the procedures "Sidak", "AllorNone" that evaluates multiple endpoints in clinical trials
#' @param alpha significance level (default: 0.05)
#' @param p we assume for the correlation matrix of the endpoints compound symmetry, then the correlation between the endpoints is p (default case: uncorrelated endpoints)

#' @return for @param proceedure = "Sidak" we get a dataframe with variables @param N, @param m, @param eta, @param p, mean actual FWER and P(FWER<=0.05)
#'         for @param procedure ="AllorNone" we get a dataframe with variables @param N, @param m, @param eta, mean actual T1E and P(T1E<=0.05)


multEndp<- function(N,m,sigma, randomization, eta, r_number,procedure,alpha=0.05,CHEN_par1=2,CHEN_par2=0.67,BSD_par=3,EBC_par=0.67,MP_par=3,PBR_par=4, p=NULL){
  set.seed(1)
  
  #'Checking that the inputs of @function multEndp are in the correct form.
  if(N<=0 | m<=0 | r_number<=0){
    stop('Incorrect input for N, m or r_number')
  }
  if( length(eta)!=m && length(which(eta<0))!=0){
    stop("Check the size of eta and if eta>=0")
  }
  if(!is.null(p)){
    if( p<0 | p>1){
    stop('Incorrect input for p')
  }}
  
  #' Simulation of the actual mean FWER and the probabilities of obtaining an actual FWER that maintain the  significance level of the Sidak procedure
  if(procedure=='Sidak'){
    #'According to the randomization procedure defined in @param randomization, a Monte Carlo sample of randomization lists in the size of @param r_number is generated. 
    myPar<-random_list_samp(N,m,randomization,r_number,CHEN_par1=CHEN_par1,CHEN_par2=CHEN_par2,BSD_par=BSD_par,EBC_par=EBS_par,MP_par=MP_par,PBR_par=PBR_par)
    Rand_list<- getRandList(myPar)
    
    #' Determining the sample size of the treatment and control group regarding the randomization lists of the Monte Carlo sample
    NE<-get_numb_treat(Rand_list)
    NC<-get_numb_contr(Rand_list)
    
    #' Define the adjusted Sidak significance level
    alpha_adj<-1-(1-alpha)^(1/m)
    
    FWER<-rep(NULL, r_number)
    
    #' Distinguishing between the cases of uncorrelated and correlated endpoints: We start with the uncorrelated case
    if(!is.null(p) && p!=0){ 
      #' principal component analysis of the covariance matrix to transform the model of the patient responses 
      transf=get_transf(m,sigma,Rand_list, p=p)
      sigma_new=sqrt(transf$val)
      A=transf$vec
      
      biasSeq<-getExpectation(myPar, selBias("CS", 1, "exact")) 
      
      for(j in 1:r_number){
        
        er=1
        
        #' Defining the allocation bias effect according to the biasing policy and transform it with the principal component transformation matrix
        tau<-matrix(,nrow=m, ncol=N)
        for(k in 1:N){
          tau[,k]=biasSeq[j,k]*eta
        }
        tau_new=t(A)%*%tau
        
        #' Calculating the non-centality parameter of the doubly non-central t-distribution of the t-statistic with the transformed patient responses under the null hypotheses
        delta=matrix(, nrow=r_number,ncol=m)
        lambda=matrix(,nrow=r_number, ncol=m)
        
        for (i in 1:m){
          
          param=get_noncentralparam(Rand_list[j,],sigma_new[i],NE[j],NC[j],tau_new[i,])
          delta[j,i]=param$delta
          lambda[j,i]=param$lambda
          
          #Calculating the FWER regarding each randomization list of the Monte Carlo sample 
          ub_ac<-ceiling(lambda[j,i]/2 + qpois(.995, lambda[j,i]/2))
          lb_ac<-max(floor(lambda[j,i]/2 - qpois(.995, lambda[j,i]/2)), 0)
          er<-er*(1-(doublyT(qt(alpha_adj/2,N-2),N-2,delta[j,i],lambda[j,i],lb=lb_ac,ub=ub_ac)+doublyT(qt(alpha_adj/2,N-2),N-2,-delta[j,i], lambda[j,i],lb=lb_ac,ub=ub_ac)))
        }
         
        FWER[j]<- 1-er 
      }
      
      #'Determining the actual mean FWER and the probabilities of obtaining an actual FWER that exceeds the significance level regarding the Monte Carlo sample of randomization lists
      mean_FWER<- mean(FWER)
      go_nogo<-length(which(FWER<=0.05))/(r_number)
      
      eta_mat<-matrix(eta,nrow=1, byrow=TRUE)
      out<-data.frame(N=N,m=m,eta=matrix(eta,nrow=1, byrow=TRUE), randomization=randomization, mean_FWER=mean_FWER , go_nogo=go_nogo,p=p)
      
    }else{
    #'Case of uncorrelated endpoints
    #'Defining the allocation bias effect according to the biasing policy
      biasSeq<-getExpectation(myPar, selBias("CS", 1, "exact")) 
      
      for(j in 1:r_number){
        
        er=1
        
        tau<-matrix(,nrow=m, ncol=N)
        for(k in 1:N){
          tau[,k]=biasSeq[j,k]*eta
        }
        
       
        
        delta=matrix(, nrow=r_number,ncol=m)
        lambda=matrix(,nrow=r_number, ncol=m)
        
        for (i in 1:m){
          
          #' Calculating the non-centality parameter of the doubly non-central t-distribution of the t-statistic under the null hypotheses
          param=get_noncentralparam(Rand_list[j,],sigma[i],NE[j],NC[j],tau[i,])
          delta[j,i]=param$delta
          lambda[j,i]=param$lambda
          
          #Calculating the FWER regarding each randomization list of the Monte Carlo sample
          ub_ac<-ceiling(lambda[j,i]/2 + qpois(.995, lambda[j,i]/2))
          lb_ac<-max(floor(lambda[j,i]/2 - qpois(.995, lambda[j,i]/2)), 0)
          er<-er*(1-(doublyT(qt(alpha_adj/2,N-2),N-2,delta[j,i],lambda[j,i],lb=lb_ac,ub=ub_ac)+doublyT(qt(alpha_adj/2,N-2),N-2,-delta[j,i], lambda[j,i],lb=lb_ac,ub=ub_ac)))
        }
        
        FWER[j]<- 1-er
      }
      
      #'Determining the actual mean FWER and the probabilities of obtaining an actual FWER that exceeds the significance level regarding the Monte Carlo sample of randomization lists
      mean_FWER<- mean(FWER)
      go_nogo<-length(which(FWER<=0.05))/(r_number)

      out<-data.frame(N=N,m=m,eta=matrix(eta,nrow=1, byrow=TRUE), randomization=randomization, mean_FWER=mean_FWER , go_nogo=go_nogo,p=0)
      
    }
    
#' Simulation of the actual mean T1E and the probabilities of obtaining an actual T1E that maintain the  significance level of the all-or-none procedure
  }else if(procedure=='AllorNone'){ 
    
    myPar<-random_list_samp(N,m,randomization,r_number,CHEN_par1=CHEN_par1,CHEN_par2=CHEN_par2,BSD_par=BSD_par,EBC_par=EBS_par,MP_par=MP_par,PBR_par=PBR_par)
    Rand_list<- getRandList(myPar)
    
    NE<-get_numb_treat(Rand_list)
    NC<-get_numb_contr(Rand_list)
    
    T1E<-rep(NULL, r_number)
    
    delta=matrix(, nrow=r_number,ncol=m)
    lambda=matrix(,nrow=r_number, ncol=m)
    
    biasSeq<-getExpectation(myPar, selBias("CS", 1, "exact")) 
    for(j in 1:r_number){
      
      er<-numeric(0)
      
      tau<-matrix(,nrow=m, ncol=N)
      for(k in 1:N){
        tau[,k]=biasSeq[j,k]*eta
      }
      
      for (i in 1:m){
        
        param=get_noncentralparam(Rand_list[j,],sigma[i],NE[j],NC[j],tau[i,])
        delta[j,i]=param$delta
        lambda[j,i]=param$lambda
        
        ub_ac<-ceiling(lambda[j,i]/2 + qpois(.995, lambda[j,i]/2))
        lb_ac<-max(floor(lambda[j,i]/2 - qpois(.995, lambda[j,i]/2)), 0)
        er<-c(er,doublyT(qt(alpha,N-2),N-2,-delta[j,i], lambda[j,i],lb=lb_ac,ub=ub_ac))
      }
      
      T1E[j]<- max(er)
    }
    
    mean_T1E<- mean(T1E)
    go_nogo<-length(which(T1E<=0.05))/r_number
    
    
    eta_mat<-matrix(eta,nrow=1, byrow=TRUE)
    out<-data.frame(N=N,m=m,eta=eta_mat, randomization=randomization, mean_T1E=mean_T1E , go_nogo=go_nogo)
    
  }else if(procedure== 'Lauter'){
    stop('under construction')
  } else if(procedure=='GLs'){
    stop('under construction')
  }else if( procedure=='OLS'){
    stop('under construction')
  } else{
    stop( 'Incorrect input for procedure')
  }
  return(out)
}

