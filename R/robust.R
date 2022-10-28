#' Covert p-star threshold values to the probabilities of responding in each category
#'
#' This function turns matrix of p-star threshold values into probabilities of responding in each category
#' @param pstar A matrix of p-star threshold values
#' @return The probabilities of responding in each category
pstar_to_p<-function(pstar){
  #input is an array with nrows = nitems, ncol = nthresh, dim[3] = nsubjects
  if(!is.na(dim(pstar)[3])){ #if there is only one subject
    P<-array(dim = c(dim(pstar)[1], dim(pstar)[2]+1, dim(pstar)[3]))
    for(f in 1:dim(pstar)[3]){
      P[,1,f]<-1-pstar[,1,f]
      for(g in 2:dim(pstar)[2]){
        P[,g,f]<-pstar[,g-1,f]-pstar[,g,f]
      }
      P[,(dim(pstar)[2]+1),f]<-pstar[,dim(pstar)[2],f]
    }
  }else{ #if there is more than one subject
    P<-matrix(nrow = dim(pstar)[1], ncol =dim(pstar)[2]+1)
    P[,1]<-1-pstar[,1]
    for(g in 2:dim(pstar)[2]){
      P[,g]<-pstar[,g-1]-pstar[,g]
    }
    P[,(dim(pstar)[2]+1)]<-pstar[,dim(pstar)[2]]
  }
  return(P)
}

#' Response Probability Calculation (2PL)
#'
#' This function calculates the probabilities of getting items correct given person's abilities and item parameters
#' @param thetas A vector of length p containing students’ latent traits, where p is the number of test dimensions
#' @param a A p x n matrix containing fixed item slope parameters
#' @param d A vector of length n containing fixed intercept parameters, for n test items
#' @return A vector of probabilities of endorsing each of n test items
probs.gen <- function(thetas, a, d){
  thetas <- matrix(thetas, nrow = dim(a)[2], byrow = T)
  r <- a%*%thetas+d
  P<- 1/(1+exp((-1.7)*r))
  return(list(residual=r, P=matrix(P, ncol=1)))
}

#' Response Probability Calculation (2PL) and response category probabilities (P)
#'
#' This function generates threshold probabilities (p-star) and response category probabilities (P)
# using parameters needed for the GRM
#' @param thetas An array of person abilities
#' @param a A matrix of slope parameters
#' @param b Array of intercept parameters
#' @return list of threshold probabilities (p-star) and the corresponding response category probabilities (P)
probs.gen.grm<-function(thetas, a, b){
  l<-length(thetas) #number of subjects
  n<-length(a) # test length
  nthresh<-length(b)/n
  #generate threshold probabilities based on theta
  exponent<-apply(matrix(thetas), 1, function(theta) a*(theta-b))
  exponent<-array(exponent, dim = c(n,nthresh,l))
  pstar<-1/(1+exp(-1.7*exponent))
  P<-pstar_to_p(pstar)
  return(list(pstar = pstar, P = P))
}

#' Bisquare Weighting Function
#'
#' This function returns the value of the bisquare weight given a residual and bisquare tuning paramter
#' @param r A residual that measures the inconsistency of a response from the subject's assumed response model, on one item
#' @param B Bisquare tuning parameter
#' @return Bisquare weight value
bisquare<-function(r, B){
  w <- r
  for(i in 1:length(r)){
    if(abs(r[i]) <= B){
      w[i] <- (1-(r[i]/B)^2)^2
    }else{
      w[i] <- 0.0000001
    }
  }
  return(w)
}

#' Huber Weighting Function
#'
#' This function returns the value of the Huber weight given a residual and Huber tuning paramter
#' @param r A residual that measures the inconsistency of a response from the subject's assumed response model, on one item
#' @param H Huber tuning parameter
#' @return Huber weight value
huber<-function(r, H){
  w <- r
  for(i in 1:length(r)){
    if(abs(r[i]) <= H){
      w[i] <- 1
    }else{
      w[i] <- H/abs(r[i])
    }
  }
  return(w)
}

#' Data generation function from probabilities of responding in each category
#'
#' This function generates data from probabilities of responding in each category
#' @param H Huber tuning parameters
#' @return The probabilities of responding in each category
data.gen<-function(P){
  if(is.array(P)){ #if more than one subject
    dat<-matrix(nrow = nrow(P), ncol = dim(P)[3])
    for(i in 1:dim(P)[3]){
      dat[,i]<-apply(P[,,i], 1, function(x) sample(c(1:dim(P)[2]), 1, replace = T, prob=x))
    }
  }else{#if only one subject
    dat<-apply(P, 1, function(x) sample(c(1:dim(P)[2]), 1, replace = T, prob=x))
  }
  return(dat)
}

#' Ability Estimation Function Using Robust Estimation
#'
#' This function return the list of ability estimations based on the given weighting function
#' @param dat A \eqn{n \times m} matrix of dichotomously-coded data where 1 represents an endorsed response and 0 represents a non-endorsed response. \emph{n} is the test length and \emph{m} is the number of subjects.
#' @param a A \eqn{p\timesn} matrix containing fixed item slope parameters for \emph{n} items and \emph{p} dimensions
#' @param d A vector of length \emph{n} containing fixed intercept parameters for n items
#' @param iter Maximum number of iterations for the Newton-Rapshon method. Default is 100.
#' @param cutoff Threshold value to terminate the iteration when the likelihood changes below this value, which means that the estimation is converged. Default is 0.01.
#' @param theta.initial A vector of length \emph{p} containing initial latent trait values for the maximum likelihood estimation. The default is 0 for each of the \emph{p} dimensions.
#' @param weight.category The weighting strategy to use: "equal", "bisquare", or "Huber". Default is "equal", which is equally weighted as in standard maximum likelihood estimation.
#' @param tuning.par The tuning parameter for "bisquare" or "Huber" weighting functions. Greater tuning parameters result in less downweighting.
#' @details The goal of robust estimation is to downweigh potentially aberrant responses to lessen their impact on the estimation of \eqn{\theta}. Robust estimates resist the harmful effects of response disturbances and tend to be less biased estimates of true ability than maximum likelihood estimates.
#' The contribution of item \emph{i} to the overall log-likelihood for one subject is weighted with a weight \eqn{\omega(r_i)} as a function of a residual \eqn{r_i} for the item \emph{i}:
#' \deqn{\sum_i^n \omega(r_i) \frac{\partial}{\partial\boldsymbol\theta} \ln L(\boldsymbol\theta;x_i) = 0 }
#' The residual, which measures the inconsistency of a response from the subject's assumed response model, is \deqn{r_i = \textbf a_i\boldsymbol\theta + d_i } in the multidimensional case.
#' Two types of weight functions are used: Tukey's bisquare weighting function (Mosteller & Tukey, 1977)
#' \deqn{\omega(r_i)=\begin{cases}[1-(r_i/B)^2]^2, & \text{if} |r_i|\leq B.\\0, & \text{if} |r_i|>B.\end{cases}}
#' and the Huber weighting function (Huber, 1981)
#' \deqn{\omega(r_i)=\begin{cases}1, & \text{if} |r_i|\leq H.\\H/|r_i|, & \text{if} |r_i|>H.\end{cases}}
#' Both functions are effective in estimating more accurate scores with aberrant data, although the bisquare weight function may lead to nonconvergence when using data containing a high proportion of incorrect responses (Schuster & Yuan, 2011).
#' Huber, P. (1981) \emph{Robust Statistics}. Wiley, New York. https://doi.org/10.1002/0471725250
#' Mosteller, F., & Tukey, J. W. (1977). \emph{Data Analysis and Regression: A Second Course in Statistics}. Reading, MA: Addison-Wesley Pub Co.
#' Schuster, C., & Yuan, K.-H. (2011). Robust Estimation of Latent Ability in Item Response Models. \emph{Journal of Educational and Behavioral Statistics}, 36(6), 720–735. https://doi.org/10.3102/1076998610396890
#' @return An array of estimated person abilities
#' @export

theta.est<-function(dat, a, d, iter=30, cutoff=.01, theta0=rep(0,ncol(a)), weight.type="equal", tuning.par=NULL){
  # first to check if the turning parameter is given when the weight.type is not "normal"
  if (weight.type != "equal") {
    if (is.null(tuning.par)) {
      stop(paste("The turning parameter cannot be null when the weight.type is ", weight.type, sep = ""))
    }
  }
  dim<-ncol(a) #number of dimensions
  n<-length(d) #number of items
  N<-ncol(dat) #number of subjects
  D<-1.7

  theta.out<-matrix(nrow = N, ncol = dim)
  convergence<-matrix(0, nrow=N, ncol = dim)
  theta.progression<-array(NA, dim = c(dim, iter, N))
  residual<-matrix(data=NA, nrow = n, ncol = N)

  for (i in 1:N){ #looping over subjects if more than one
    theta<-matrix(theta0, nrow = ncol(a), byrow = T)
    P0<-0 #starting probability for calculating likelihood for convergence
    for(j in 1:iter){ #prespecified iterations of Newton-Raphson Algorithm
      ex<-a%*%theta+d #residual
      P<-1/(1+exp(-D*ex))
      Hess<-matrix(ncol = dim, nrow = dim)
      D1<-matrix(nrow=dim)
      for(k in 1:dim){
        weighting.term <- NULL
        if (weight.type == "bisqure") {
          weighting.term <- bisquare(ex, tuning.par)
        } else if (weight.type == "Huber") {
          weighting.term <- huber(ex, tuning.par)
        } else {
          weighting.term <- 1
        }
        if (is.null(weighting.term)) {
          stop("Cannot determine the weighting function.")
        }

        D1[k] <- D*sum(weighting.term*a[,k]*(dat-P)) #first derivative
        for(l in 1:dim){
          Hess[k,l] <- (-(D)^2)*sum(weighting.term*a[,k]*a[,l]*(1-P)*(P)) #Hessian matrix of 2nd derivatives
        }
      }
      if(is.singular.matrix(Hess)){ #check if Hessian matrix is singular
        convergence[i,]<-rep("Singular", dim)
        for(h in 1:dim){ #bound theta estimates between -3 and 3
          if(theta[h]>3){theta[h]<-NA}else if(theta[h]<(-3)){theta[h]<-NA}
        }
        break
      }
      if(any(is.na(theta-solve(Hess)%*%D1))){
        theta.out[i,]<-theta<-theta-solve(Hess)%*%D1
        convergence[i,]<-as.numeric(is.na(theta.out[i,]))
        break # break if noncomverging
      }
      theta<-theta-solve(Hess)%*%D1 #update theta
      residual[,i]<-ex
      log_like<-sum(log(P))-sum(log(P0)) #compare log-likelihood with log-likelihood from previous iteration for convergence criterion
      if(abs(log_like)<cutoff){
        break
      }
      P0<-P
      theta.progression[,j,i]<-theta
      #print(theta) #view convergence of theta for debugging purposes
    }
    if(j==iter){#if theta never converged to tolerance
      theta.out[i,]<-theta<-rep(NA, dim)
      convergence[i,]<-rep(1, dim)
    }
    for(h in 1:dim){ #bound theta estimates between -3 and 3
      if(!is.na(theta[h]) & theta[h]>3){theta[h]<-3}
      if(!is.na(theta[h]) & theta[h]<(-3)){theta[h]<-(-3)}
    }
    theta.out[i,]<-theta
  }
  return(list(Theta = theta.out, Convergence = convergence, theta.progression = theta.progression, residual=residual))
}

#' Ability Estimation Function Using Robust Estimation
#'
#' This function return the list of ability estimations based on the given weighting function
#' @param dat Response data to load
#' @param a Matrix of slope parameters
#' @param b Array of intercept parameters
#' @param iter Max number of iterations. Default is 100
#' @param cutoff Threshold value to terminate the iteration when the likelihood changes below this value, which means that the estimation is converged.
#' @param theta.initial Array of initial thetas. Default is 0s
#' @param weight.category The weighting strategy to use, "equal", "bisquare" and "Huber". Default is "equal", which is equally weighted.
#' @param tuning.par The tuning parameter for "bisquare" or "Huber" weighting functions
#' @details Bisquqre weighting function
#' \deqn{\omega(r_i)=\begin{cases}[1-(r_i/B)^2]^2, & \text{if} |r_i|\leq B.\\0, & \text{if} |r_i|>B.\end{cases}}
#' Huber weighting function
#' \deqn{\omega(r_i)=\begin{cases}1, & \text{if} |r_i|\leq H.\\H/|r_i|, & \text{if} |r_i|>H.\end{cases}}
#' where \eqn{r_i} measures the inconsistency of a response of item i from the subject's assumed response model. Mislevy and Bock proposed a residual for the unidimentional case, \deqn{r_i=a_i(\theta-b_i)}
#' @return Array of estimated person abilities
#' @export
theta.est.grm<-function(dat, a, b, iter=30, cutoff=0.01, theta.initial=rep(0,ncol(a)), weight.type="equal", tuning.par=NULL){
  # first to check if the turning parameter is given when the weight.type is not "normal"
  if (weight.type != "equal") {
    if (is.null(tuning.par)) {
      stop(paste("The turning parameter cannot be null when the weight.type is ", weight.type, sep = ""))
    }
  }
  l<-dim(dat)[2] #number of subjects
  n<-dim(dat)[1] # test length
  nthresh<-length(b)/n #number of threshold parameters

  dat.b<-array(dim = list(n, l, 2))
  theta.est2<- matrix(data=NA, nrow=l)
  convergence<-matrix(0, nrow=l)
  theta.progression<-matrix(NA, nrow = l, ncol = iter)
  residual<-matrix(data=NA, nrow = n, ncol = l)
  for (i in 1:l){
    for(j in 1:n){
      #extract threshold values of interest according to each response & stored in dat.b
      # e.g. if the response is in category "3", we will use the threshold parameters from the previous threshold and the threshold corresponding to that category
      if(dat[j,i]==1){
        dat.b[j,i,1]<- -1000
        dat.b[j,i,2]<-b[j,1]
      }else if(dat[j,i]>nthresh){
        dat.b[j,i,1]<-b[j,dat[j,i]-1]
        dat.b[j,i,2]<-1000
      }else{
        dat.b[j,i,1]<-b[j,dat[j,i]-1]
        dat.b[j,i,2]<-b[j,dat[j,i]]
      }
    }
  }
  for(i in 1:l){
    theta<-init.val
    P0<-0
    for (k in 1:iter){
      exponent_0<-a*(theta-dat.b[,i,1])
      exponent_1<-a*(theta-dat.b[,i,2])
      #P*k-1 and P*k
      ps0<-1/(1+exp(-1.7*exponent_0))
      ps1<-1/(1+exp(-1.7*exponent_1))
      qs0<-1-ps0
      qs1<-1-ps1
      #P_k
      P<-ps0-ps1

      pstar<-1/(1+exp(-1.7*(a*(theta-b))))
      probs<-pstar_to_p(pstar)
      expected.value<-probs%*%matrix(c(1:(nthresh+1))) #expected response
      residual[,i]<-(dat[,i]-expected.value)/nthresh #residual based on the difference between the bserved and expected response
      weighting.term <- NULL
      if (weight.type == "bisqure") {
        weighting.term <- bisquare(residual[,i], tuning.par)
      } else if (weight.type == "Huber") {
        weighting.term <- huber(residual[,i], tuning.par)
      } else {
        weighting.term <- 1
      }
      if (is.null(weighting.term)) {
        stop("Cannot determine the weighting function.")
      }

      D1<-sum(1.7*a*weighting.term*(ps0*qs0-ps1*qs1)/P) #first derivative
      D2<-sum(1.7^2*a^2*weighting.term*( (ps0*qs0*(qs0-ps0)-ps1*qs1*(qs1-ps1))/P - (ps0*qs0-ps1*qs1)^2/P^2 )) #second derivative
      if(is.na(theta-D1/D2)){ #check if we received NAs
        theta.est2[i]<-theta<-NA
        convergence[i,1]<-1
        break
      }
      theta<-theta.progression[i,k]<-theta-D1/D2
      log_like<-sum(log(P))-sum(log(P0)) #loglikelihood for convergence criterion
      if(abs(log_like)<cutoff){
        break
      }
      P0<-P
    }
    if(k==iter){#if theta never converged to within tolerance
      theta.est2[i]<-theta<-NA
      convergence[i,1]<-1
    }else if(!is.na(theta) & theta< -3){ #replace thetas that converged outside [-3, 3]
      theta<--3
    }else if(!is.na(theta) & theta>3){
      theta<-3
    }
    theta.est2[i]<-theta
  }
  return(list(theta = theta.est2, convergence=convergence, theta.progression = theta.progression, residual=residual))
}
