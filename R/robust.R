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
      w[i] <- 0
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


#' Ability Estimation Function Using Robust Estimation
#'
#' This function return the list of ability estimations based on the given weighting function
#' @param dat A \eqn{n \times m} matrix of dichotomously-coded data where 1 represents an endorsed response and 0 represents a non-endorsed response. \emph{n} is the test length and \emph{m} is the number of subjects.
#' @param a A \eqn(p\timesn} matrix containing fixed item slope parameters for \emph{n} items and \emph{p} dimensions
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
theta.est <- function(dat, a, d, iter=100, cutoff=0.01, theta.initial=rep(0,ncol(a)), weight.type="equal", tuning.par=NULL){
  # first to check if the turning parameter is given when the weight.type is not "normal"
  if (weight.category != "equal") {
    if (is.null(tuning.par)) {
      stop(paste("The turning parameter cannot be null when the weight.type is ", weight.category, sep = ""))
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
        break # break if not converging
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
