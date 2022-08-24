#' Response Probability Calculation (2PL)
#'
#' This function calculates the probabilities of get items correct given person's abilities and item parameters
#' @param thetas An array of person abilities
#' @param a A matrix of slope parameters
#' @param d Array of intercept parameters
#' @return Array of probabilities of having correct response
probs.gen <- function(thetas, a, d){
  thetas <- matrix(thetas, nrow = dim(a)[2], byrow = T)
  r <- a%*%thetas+d
  P<- 1/(1+exp((-1.7)*r))
  return(list(residual=r, P=matrix(P, ncol=1)))
}

#' Bisquare Weighting Function
#'
#' This function return the value of bisquare weighting value given Residual and bisquare tuning paramter
#' @param r Residual that measures the inconsistency of a response from the subject's assumed response model
#' @param B Bisquare tuning parameters
#' @return Bisquare weighting function
bisquare<-function(r, B){
  w <- r
  for(i in 1:length(r)){
    if(abs(r[i]) <= B){
      w[i] <- (1-(r[i]/B)^2)^2
    }else{
      w[i] <- 0
    }
  }
  return(w)
}

#' Huber Weighting Function
#'
#' This function return the value of Huber weighting value given Residual and Huber tuning paramter
#' @param r Residual that measures the inconsistency of a response from the subject's assumed response model
#' @param H Huber tuning parameters
#' @return Huber weighting function
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


#' Ability Estimation Function using MLE
#'
#' This function return the list of ability estimations based on the given weighting function
#' @param dat Response data to load
#' @param a Matrix of slope parameters
#' @param d Array of intercept parameters
#' @param iter Max number of iterations. Default is 100
#' @param theta.initial Array of initial thetas. Default is 0s
#' @param weight.category The weighting strategy to use, "normal", "bisquare" and "Huber". Default is "normal", which is equally weighted.
#' @param tuning.par The tuning parameter for "bisquare" or "Huber" weighting functions
#' @details Huber weighting function
#' @return Array of estimated person abilities
#' @export
theta.recov <- function(dat, a, d, iter=100, theta.initial=rep(0,ncol(a)), weight.category="normal", tuning.par=NULL){
  # first to check if the turning parameter is given when the weight.category is not "normal"
  if (weight.category != "normal") {
    if (is.null(tuning.par)) {
      stop(paste("The turning parameter cannot be null when the weighting strategy is ", weight.category, sep = ""))
    }
  }
  dim <- ncol(a) #number of dimensions
  n <- length(d) #number of items
  N <- length(dat)/n #number of subjects
  D <- 1.7
  theta.out <- matrix(nrow = N, ncol = dim)
  for (i in 1:N){ #looping over subjects if more than one
    theta <- matrix(theta.initial, nrow = ncol(a), byrow = T)
    for(j in 1:iter){ #prespecified iterations of Newton-Raphson Algorithm
      ex <- a%*%theta+d #residual
      P <- 1/(1+exp(-D*ex))
      Hess <- matrix(ncol = dim, nrow = dim)
      D1 <- matrix(nrow=dim)
      for(k in 1:dim){
        weighting.term <- NULL
        if (weight.category == "bisqure") {
          weighting.term <- bisquare(ex, tuning.par)
        } else if (weight.category == "Huber") {
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
      if(any(is.na(theta-solve(Hess)%*%D1))){
        break # break if noncomverging
      }
      theta <- theta-solve(Hess)%*%D1 #update theta
      # print(theta) #view convergence of theta for debugging purposes
      for(h in 1:dim){ #bound theta estimates between -3 and 3
        if(theta[h] > 3){theta[h] <- 3}
        if(theta[h] < -3){theta[h] <- -3}
      }
    }
    theta.out[i,] <- theta
  }
  return(theta.out)
}
