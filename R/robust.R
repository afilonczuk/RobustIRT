#' Convert \eqn{P^*} threshold values to the probabilities of responding in each category
#'
#' This function calculates \eqn{P}, the probabilities of responding in each category, from the \eqn{P^*} threshold values using the graded response model.
#' The probability that a subject responds in or above a category \eqn{k} for item \eqn{i} is \eqn{P^*_{ik}(\theta) = \frac{1}{1+ e^{-a_i (\theta-b_{ik})}}},
#  for \eqn{K} categories and \eqn{K-1} threshold parameters  ($b_{i1}, ... b_{i,K-1}$), where $b_{ik}$ separates response category \eqn{k} and \eqn{k+1} (\eqn{k=1,..K-1}) (Embretson & Reise 2000).
#  \eqn{a_i} is the item discrimination parameter.  The probability of endorsing exactly category \eqn{k} is \eqn{P_{ik}(\theta) = P^*_{i,k}(\theta) - P^*_{i,k+1}(\theta),} where \eqn{P^*_{i1}(\theta) \equiv 1.0} and \eqn{P^*_{iK}(\theta) \equiv 0.0.}
#' @param pstar A \eqn{n \times K-1 \times m} array of \eqn{P^*} threshold probability values, for \eqn{K} categories over \eqn{n} items and \eqn{m} subjects
#' @return The probabilities \eqn{P} of responding in each category
pstar_to_p<-function(pstar){
  if(!is.na(dim(pstar)[3])){ 
    # If there is only one subject
    # Initialize array for probabilities
    P<-array(dim = c(dim(pstar)[1], dim(pstar)[2]+1, dim(pstar)[3]))
    for(f in 1:dim(pstar)[3]){
      # Calculate P_{j1} as 1 - P*_{j1}
      P[,1,f]<-1-pstar[,1,f]
      for(g in 2:dim(pstar)[2]){
        # Calculate P_{jk} as P*_{j,k-1} - P*_{jk}
        P[,g,f]<-pstar[,g-1,f]-pstar[,g,f]
      }
      # Calculate P_{jK} as P*_{j,K-1}
      P[,(dim(pstar)[2]+1),f]<-pstar[,dim(pstar)[2],f]
    }
  }else{ 
    # If there is more than one subject
    # Initialize matrix for probabilities
    P<-matrix(nrow = dim(pstar)[1], ncol =dim(pstar)[2]+1)
    # Calculate P_{j1} as 1 - P*_{j1}
    P[,1]<-1-pstar[,1]
    for(g in 2:dim(pstar)[2]){
      # Calculate P_{jk} as P*_{j,k-1} - P*_{jk}
      P[,g]<-pstar[,g-1]-pstar[,g]
    }
    # Calculate P_{jK} as P*_{j,K-1}
    P[,(dim(pstar)[2]+1)]<-pstar[,dim(pstar)[2]]
  }
   # Return the matrix/array of category response probabilities
  return(P)
}

#' Response Probability Calculation (2PL)
#'
#' This function calculates the probabilities of getting items correct given the item parameters and a person's abilities
#  \eqn{\boldsymbol\theta = ({\theta}_1, {\theta}_2, ... {\theta}_p)'} on \emph{p} dimensions, according to the multidimensional 2PL IRT model:
#  \eqn{P(X=1 | \boldsymbol{\theta}) =  \frac{1}{1+ e^{-1.7(\boldsymbol{a_i}\boldsymbol\theta+d_i)}}}.
#  An item \emph{i} has slope parameters \eqn{\boldsymbol{a_i}=(a_{1i}, a_{2i}, ..., a_{pi})'} and intercept parameter \eqn{d_i}.
#' @param thetas A vector of length \emph{p} containing students’ latent traits, where \emph{p} is the number of test dimensions
#' @param a A \eqn{n \times p} matrix containing fixed item slope parameters for \emph{p} dimensions and \emph{n} test items
#' @param d A vector of length \emph{n} containing fixed intercept parameters, for \emph{n} test items
#' @return P A vector of probabilities for endorsing each of \emph{n} test items
#' @return residual A vector of residuals detecting misfit between the subject and \emph{n} test items, according to the formula \eqn{r_i = \textbf a_i\boldsymbol\theta + d_i} (see theta.est)
probs.gen <- function(thetas, a, d){
  thetas <- matrix(thetas, nrow = dim(a)[2], byrow = T)
  r <- a%*%thetas+d
  P<- 1/(1+exp((-1.7)*r))
  return(list(residual=r, P=matrix(P, ncol=1)))
}

#' Response Probability Calculation (GRM) and response category probabilities (P)
#'
#' This function generates \eqn{P^*}, the threshold probability, and \eqn{P}, the probability of responding in each category for a vector of latent traits
#  using parameters needed for the GRM
#  The probability that a subject responds in or above a category \eqn{k} for item \eqn{i} is \eqn{P^*_{ik}(\theta) = \frac{1}{1+ e^{-a_i (\theta-b_{ik})}}},
#  for \eqn{K} categories and \eqn{K-1} threshold parameters  ($b_{i1}, ... b_{i,K-1}$), where $b_{ik}$ separates response category \eqn{k} and \eqn{k+1} (\eqn{k=1,..K-1}) (Embretson & Reise 2000).
#  \eqn{a_i} is the item discrimination parameter.  The probability of endorsing exactly category \eqn{k} is \eqn{P_{ik}(\theta) = P^*_{i,k}(\theta) - P^*_{i,k+1}(\theta),} where \eqn{P^*_{i1}(\theta) \equiv 1.0} and \eqn{P^*_{iK}(\theta) \equiv 0.0.}
#' @param thetas A vector of latent traits for \eqn{m} subjects
#' @param a A vector of discrimination parameters for \eqn{n} test items
#' @param b A matrix of category threshold parameters, with the number of rows corresponding to \eqn{n} test items and the number of columns corresponding to \eqn{K-1} thresholds between \eqn{K} categories
#' @return pstar A \eqn{n \times (K-1) \times m} array of category threshold probabilities for \eqn{n} items, \eqn{K} categories, and \eqn{m} subjects
#' @return P a \eqn{n \times K \times m} array of probabilities of responding in each of \eqn{K} categories for \eqn{n} items, across \eqn{m} subjects
#' @examples
#' thetas <- rnorm(15)  # Example latent traits for 15 subjects
#' a <- runif(10, 0.5, 1.5)  # Example discrimination parameters for 10 items
#' b <- t(apply(matrix(runif(10*4, -2.5,2.5), nrow = 10, ncol = 4), 1, sort))  # Example threshold parameters for 10 items and 4 thresholds (5 categories)
#' result <- probs.gen.grm(thetas, a, b)  # Calculate probabilities
probs.gen.grm<-function(thetas, a, b){
  # Number of subjects
  l<-length(thetas)
  # Number of items
  n<-length(a) 
  # Number of thresholds between categories
  nthresh<-length(b)/n
  
  # Generate threshold probabilities based on theta
  exponent <- apply(matrix(thetas), 1, function(theta) a*(theta-b))
  exponent <- array(exponent, dim = c(n,nthresh,l))
  # Calculate P* using the GRM formula
  pstar <- 1/(1+exp(-1.7*exponent))
  # Convert P* to category probabilities using the pstar_to_p function
  P<-pstar_to_p(pstar)

  # Return both P* and P
  return(list(pstar = pstar, P = P))
}

#' Bisquare Weighting Function
#'
#' This function returns the value of the bisquare weight given a residual and bisquare tuning paramter
#' @param r A residual that measures the inconsistency of a response from the subject's assumed response model, on one item. Residuals that are NA are given a weight of 0.
#' @param B Bisquare tuning parameter
#' @return Bisquare weight value
bisquare<-function(r, B){
   w<-ifelse(is.nan(r), 0, 
            ifelse(abs(r) <= B, (1-(r/B)^2)^2, 0.0000001))
  
  return(w)
}

#' Huber Weighting Function
#'
#' This function returns the value of the Huber weight given a residual and Huber tuning paramter
#' @param r A residual that measures the inconsistency of a response from the subject's assumed response model, on one item. Residuals that are NA are given a weight of 0.
#' @param H Huber tuning parameter
#' @return Huber weight value
huber<-function(r, H){
  w<-ifelse(is.nan(r), 0, 
            ifelse(abs(r) <= H, 1, H/abs(r)))
  return(w)
}

#' Likert-type data generation function from probabilities of responding in each category
#'
#' This function generates Likert-type data from probabilities of responding in each category
#' @param P A \eqn{n \times K \times m} array of probabilities of responding in each of \emph{K} categories over \emph{n} items for \emph{m} subjects
#' @return dat A \eqn{n \times m} matrix of randomly generated Likert-type data using the sample() function for \emph{m} subjects over \emph{n} items
#' @examples
#' thetas <- rnorm(15)  # Example latent traits for 15 subjects
#' a <- runif(10, 0.5, 1.5)  # Example discrimination parameters for 10 items
#' b <- t(apply(matrix(runif(10*4, -2.5,2.5), nrow = 10, ncol = 4), 1, sort))  # Example threshold parameters for 10 items and 4 thresholds (5 categories)
#' probs <- probs.gen.grm(thetas, a, b)  # Calculate probabilities
#' data <- data.gen(probs$P)  # Generate Likert-type data

data.gen<-function(P){
  if(is.array(P)){ 
    # If more than one subject
    # Initialize matrix for generated data
    dat<-matrix(nrow = nrow(P), ncol = dim(P)[3])
    for(i in 1:dim(P)[3]){
      # Loop over each subject
      # For each item, sample a response category based on the probabilities
      dat[,i]<-apply(P[,,i], 1, function(x) sample(c(1:dim(P)[2]), 1, replace = T, prob=x))
    }
  }else{
    # If only one subject
    # For each item, sample a response category based on the probabilities
    dat<-apply(P, 1, function(x) sample(c(1:dim(P)[2]), 1, replace = T, prob=x))
  } 
  # Return the generated Likert-type data
  return(dat)
}

#' Ability Estimation Function Using Robust Estimation - Dichotomous Data
#'
#' This function return the list of ability estimations based on the given weighting function
#' @param dat A \eqn{n \times m} matrix of dichotomously-coded data where 1 represents an endorsed response and 0 represents a non-endorsed response. \emph{n} is the test length and \emph{m} is the number of subjects.
#' @param a A \eqn{p \times n} matrix containing fixed item slope parameters for \emph{n} items and \emph{p} dimensions
#' @param d A vector of length \emph{n} containing fixed intercept parameters for n items. If the model is unidimensional, this vector should contain difficulty parameters corresponding to \emph{b} in the unidimensional 2PL model.
#' @param iter Maximum number of iterations for the Newton-Raphson method. Default is 100.
#' @param cutoff Threshold value to terminate the iteration when the likelihood changes below this value, which means that the estimation is converged. Default is 0.01.
#' @param init.val A vector of length \emph{p} containing initial latent trait values for the maximum likelihood estimation. The default is 0 for each of the \emph{p} dimensions.
#' @param weight.category The weighting strategy to use: "equal", "bisquare", or "Huber". Default is "equal", which is equally weighted as in standard maximum likelihood estimation.
#' @param tuning.par The tuning parameter for "bisquare" or "Huber" weighting functions. Greater tuning parameters result in less downweighting.
#' @details The goal of robust estimation is to downweigh potentially aberrant responses to lessen their impact on the estimation of \eqn{\theta}. Robust estimates resist the harmful effects of response disturbances and tend to be less biased estimates of true ability than maximum likelihood estimates.
#' Using the multidimensional IRT model for dichotomous data, the probability of a correct response (e.g., "1") on item \emph{i} is given by the formula \eqn{P(X=1 | \boldsymbol{\theta}) =  \frac{1}{1+ e^{-1.7(\boldsymbol{a_i}\boldsymbol\theta+d_i)}},} where \eqn{\boldsymbol{a}} is a \eqn{p \times n} matrix of item slope parameters for \emph{n} items and \emph{p} dimensions. Note that \emph{d} here is considered the intercept parameter and itself is not the item difficulty as is \emph{b} in the 2PL model for unidimensional dichotomous data.
#' Unidimensional robust estimations can also be calculated from this model when input parameters have \emph{p}=1 dimensions. The 2PL probability model for unidimensional data is \eqn{P(X=1 | \theta) = \frac{1}{1+ e^{-1.7a_{i}(\theta-b_{i})}}} where \emph{a} is the same as the multidimensional case where \emph{p=1}.
#' The contribution of item \emph{i} to the overall log-likelihood for one subject is weighted with a weight \eqn{\omega(r_i)} as a function of a residual \eqn{r_i} for the item \emph{i}:
#' \deqn{\sum_i^n \omega(r_i) \frac{\partial}{\partial\boldsymbol\theta} \ln L(\boldsymbol\theta;x_i) = 0 }
#' The residual, which measures the inconsistency of a response from the subject's assumed response model, is \deqn{r_i = \textbf a_i\boldsymbol\theta + d_i } in the multidimensional case.
#' Two types of weight functions are used: Tukey's bisquare weighting function (Mosteller & Tukey, 1977)
#' \deqn{\omega(r_i)=\begin{cases}[1-(r_i/B)^2]^2, & \text{if} |r_i|\leq B.\\0, & \text{if} |r_i|>B.\end{cases}}
#' and the Huber weighting function (Huber, 1981)
#' \deqn{\omega(r_i)=\begin{cases}1, & \text{if} |r_i|\leq H.\\H/|r_i|, & \text{if} |r_i|>H.\end{cases}}
#' Both functions are effective in estimating more accurate scores with aberrant data, although the bisquare weight function may lead to nonconvergence when using data containing a high proportion of incorrect responses (Schuster & Yuan, 2011).
#' @references Huber, P. (1981) \emph{Robust Statistics}. Wiley, New York. https://doi.org/10.1002/0471725250
#' @references Mosteller, F., & Tukey, J. W. (1977). \emph{Data Analysis and Regression: A Second Course in Statistics}. Reading, MA: Addison-Wesley Pub Co.
#' @references Schuster, C., & Yuan, K.-H. (2011). Robust Estimation of Latent Ability in Item Response Models. \emph{Journal of Educational and Behavioral Statistics}, 36(6), 720–735. https://doi.org/10.3102/1076998610396890
#' @return Theta A \eqn{m \times p} matrix of ability estimates for \emph{m} subjects and \emph{p} dimensions
#' @return Convergence \eqn{m \times p} matrix containing indicators of convergence for \emph{m} subjects: a “0” indicates the value converged; a “1” indicates the maximum likelihood estimation did not converge to any value; and “Singular” indicates the Hessian matrix was singular and could not be used to continue the maximum likelihood estimation.
#' @return theta.progression A \eqn{p \times iter \times m} array for \emph{p} dimensions, \emph{iter} number of iterations supplied to the input, and \emph{m} subjects. Each column provides the updated theta estimate at each iteration of the Newton-Raphson algorithm until the change in log-likelihood for that subject reaches the cutoff value, infinite values (nonconverged), or encounters a singular matrix error.
#' @return residual A \eqn{n \times m} matrix with residuals corresponding to the ability estimate for \emph{m} subjects respective to the \emph{n} test items
#' @export
#'
#' @examples
#' # Test length
#' n <- 30
#'
#' ## Number of iterations of Newton's method
#' iter <- 15
#'
#' ## Number of dimensions
#' dim <- 3
#'
#' ## Correlation between one person's thetas
#' cor <- .6
#'
#' ## Covariance matrix for generating thetas
#' sigma <- matrix(cor, ncol = dim, nrow = dim)
#' diag(sigma) <- 1
#' s <- rep(1, dim)
#' sigma <- s*sigma*s
#'
#' ## Generate real thetas
#' thetas <- matrix(mvrnorm(n = 1, rep(0,dim), Sigma = sigma), ncol = dim)
#'
#' ## Generate slope parameters
#' a <- matrix(runif(n*dim, .5, 1.5), nrow = n, ncol = dim)
#' ## Generate intercept parameters
#' d <- matrix(rnorm(n), ncol = 1)
#' ## Calculate probabilities
#' probs <- probs.gen(thetas, a, d)
#' ## Generate data from probabilities
#' dat <- apply(probs$P, c(1,2), function(x) rbinom(1,1,x))
#'
#' ## Estimate thetas
#' theta.est(dat, a, d, iter = 15, theta0=matrix(rep(0,dim)), weight.type="Huber", tuning.par=1)


theta.est<-function(dat, a, d, iter=30, cutoff=.01, init.val=rep(0,ncol(a)), weight.type="equal", tuning.par=NULL){
  # first to check if the turning parameter is given when the weight.type is not "equal"
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
  residual<-array(data=NA, dim = c(n, N, iter))

  for (i in 1:N){ #looping over subjects if more than one
    if(length(init.val)>dim){
      theta<-matrix(init.val[i,], nrow = ncol(a), byrow = T)
    }else{
      theta<-matrix(init.val, nrow = ncol(a), byrow = T)
    }
    P0<-0 #starting probability for calculating likelihood for convergence
    for(j in 1:iter){ #prespecified iterations of Newton-Raphson Algorithm
      if(dim==1){
        theta<-as.vector(theta)
        ex<-a*(theta-d)
        }else{ex<-a%*%theta+d} #residual
      P<-1/(1+exp(-D*ex))
      Hess<-matrix(ncol = dim, nrow = dim)
      D1<-matrix(nrow=dim)
      for(k in 1:dim){
        weighting.term <- NULL
        if (weight.type == "bisquare") {
          weighting.term <- bisquare(ex, tuning.par)
        } else if (weight.type == "Huber") {
          weighting.term <- huber(ex, tuning.par)
        } else {
          weighting.term <- 1
        }
        if (is.null(weighting.term)) {
          stop("Cannot determine the weighting function.")
        }

        D1[k] <- D*sum(weighting.term*a[,k]*(dat[,i]-P)) #first derivative
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
      residual[,i,j]<-ex
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

#' Ability Estimation Function Using Robust Estimation (GRM)
#'
#' This function returns ability estimations based response data, item parameters for the graded respose model, and the given weighting function
#' @param dat A \eqn{n \times m} matrix of polytomously-scored data (e.g., Likert-type) for \emph{n} items and \emph{m} subjects.
#' @param a Vector of slope parameters for \emph{n} items
#' @param b A \eqn{n \times (K-1)} matrix of category threshold parameters for \emph{K} categories
#' @param iter Max number of iterations. Default is 100
#' @param cutoff Threshold value to terminate the iteration when the likelihood changes below this value, which means that the estimation is converged.
#' @param init.val Vector of initial latent trait for the maximum likelihood estimation for \emph{m} subjects. If a single value is provided, that initial value will be used for all subjects. Default is 0.
#' @param weight.category The weighting strategy to use, "equal", "bisquare" and "Huber". Default is "equal", which is equally weighted as in standard maximum likelihood estimation.
#' @param tuning.par The tuning parameter for "bisquare" or "Huber" weighting functions. Greater tuning parameters result in less downweighting in robust estimation.
#' @details The goal of robust estimation is to downweigh potentially aberrant responses to lessen their impact on the estimation of \eqn{\theta}. Robust estimates resist the harmful effects of response disturbances and tend to be less biased estimates of true ability than maximum likelihood estimates.
#' Under the graded response model (GRM), the probability that a subject responds in or above a category $k$ for item $i$ is \eqn{P^*_{ik}(\theta) = \frac{1}{1+ e^{-a_i (\theta-b_{ik})}}}  (Embretson, 2000). \eqn{a_i} is the item discrimination parameter. There are \emph{K} categories and \eqn{K-1} threshold parameters (\eqn{b_{i1}, ... b_{i,K-1}}), where \eqn{b_{ik}} separates response category \eqn{k} and \eqn{k+1} (\eqn{k=1,..K-1}).
#' The probability of endorsing exactly category \eqn{k} is therefore: \eqn{P_{ik}(\theta) = P^*_{i,k}(\theta) - P^*_{i,k+1}(\theta),} where \eqn{P^*_{i1}(\theta) \equiv 1.0} and \eqn{P^*_{iK}(\theta) \equiv 0.0.}
#' The contribution of item \emph{i} to the overall log-likelihood for one subject is weighted with a weight \eqn{\omega(r_i)} as a function of a residual \eqn{r_i} for the item:
#' \deqn{\sum^n_{i=1} \omega(r_i) \sum^K_{k=1} u_{ik}\text{log}P_{ik} = 0 }
#’ \eqn{u_{ik}} is an indicator function: \deqn{u_{ik} = \begin{cases}
#      1 & \text{if } X_i = k; \\
#      0 & \text{otherwise}.
#   \end{cases} }
#' The residual, which measures the inconsistency of a response from the subject's assumed response model, is \deqn{r_i = \frac{1}{K-1}\left[X_i - E(X_i|\hat{\theta})\right]} for the GRM.
#  The difference in fit is determined between the observed response \eqn{X_i} and expected score \eqn{E(X_i|\theta) = \sum_{k=1}^KkP_{ik}(\theta)}
#' Two types of weight functions are used: Tukey's bisquare weighting function (Mosteller & Tukey, 1977)
#' \deqn{\omega(r_i)=\begin{cases}[1-(r_i/B)^2]^2, & \text{if} |r_i|\leq B.\\0, & \text{if} |r_i|>B.\end{cases}}
#' and the Huber weighting function (Huber, 1981)
#' \deqn{\omega(r_i)=\begin{cases}1, & \text{if} |r_i|\leq H.\\H/|r_i|, & \text{if} |r_i|>H.\end{cases}}
#' Both functions are effective in estimating more accurate scores with aberrant data, although the bisquare weight function may lead to nonconvergence when using data containing a high proportion of incorrect responses (Schuster & Yuan, 2011).
#' @references Huber, P. (1981) \emph{Robust Statistics}. Wiley, New York. https://doi.org/10.1002/0471725250
#' @references Mosteller, F., & Tukey, J. W. (1977). \emph{Data Analysis and Regression: A Second Course in Statistics}. Reading, MA: Addison-Wesley Pub Co.
#' @references Schuster, C., & Yuan, K.-H. (2011). Robust Estimation of Latent Ability in Item Response Models. \emph{Journal of Educational and Behavioral Statistics}, 36(6), 720–735. https://doi.org/10.3102/1076998610396890
#' @return theta Ability estimates for \emph{m} subjects. NAs replace values that did not converge to any value. Estimates that converged to values less than -3.0 were replaced with -3.0, while estimates that converged to values greater than 3.0 were replaced with 3.0.
#' @return convergence Indicators of convergence for \emph{m} subjects: a “0” indicates the value converged, while a “1” indicates the maximum likelihood estimation did not converge to any value.
#' @return standard.error Standard errors of the theta estimates for \emph{m} subjects, given by the square root of the reciprocal of the Fisher information. NAs replace nonconverging values. 
#' @return theta.progression A matrix with rows corresponding to each subject and columns corresponding to the number of iterations supplied to the input. Each column provides the updated theta estimate at each iteration of the Newton-Raphson algorithm until the change in log-likelihood for that subject reaches the cutoff value or the value is nonconverged (reaches infinite values).
#' @return residual A \eqn{n \times m \times iter} array containing residuals corresponding to the ability estimate for \emph{m} subjects respective to the \emph{n} test items at each iteration until convergence, nonconvergence, or singular matrix is reached.
#' @export
#' @examples
#' # Test Length
#' n<-30
#' 
#' # Number of thresholds (5-point Likert scale)
#' nthresh<-4
#' 
#' # Number of iterations of Newton's method
#' iter <- 15
#' 
#' # Set critical value for convergence criteria
#' crit.val<-0.01
#' 
#' # Set real thetas - 5 subjects
#' thetas<-c(-2,-1,0,1,2)
#' 
#' # Set item slope
#' a<-runif(n, .90, 2.15)
#' 
#' # Set threshold parameters
#' b<-t(apply(matrix(runif(n*4, -2.5,2.5), nrow = n, ncol =4), 1, sort))
#' 
#' # Calculate Probabilities
#' probs<-probs.gen.grm(thetas, a, b)
#' 
#' # Generate Likert data
#' dat<-data.gen(probs$P)
#' 
#' # Make the data aberrant by reverse coding 20% items
#' ab.prop<-0.2
#' index<-sample(c(1:n), ab.prop*n)
#' ab.dat<-dat
#' ab.dat[index, ]<-apply(matrix(dat[index,]), c(1,2), function(x) return(nthresh+2-x))
#'
#' 
#' # Calculate MLE (Non-robust)
#' mle<-theta.est.grm(ab.dat, a, b, iter, crit.val, init.val=0, weight.type="equal")
#' 
#' # Use MLE as starting value, or 0 if NA
#' start.val<-apply(mle$theta, c(1,2), function(x) ifelse(is.na(x), 0, x))
#' 
#' # Calculate bisquare- and Huber-weighted robust estimates
#' b.est<- theta.est.grm(ab.dat, a, b, iter, crit.val, init.val=start.val, weight.type="bisquare", tuning.par=4)
#' h.est<-theta.est.grm(ab.dat, a, b, iter, crit.val, init.val=start.val, weight.type="Huber", tuning.par=1)
#' 
#' # Compare robust ability estimates with MLE
#' b.est$theta
#' h.est$theta
#' mle$theta
#' 
theta.est.grm<-function(dat, a, b, iter=30, cutoff=0.01, init.val=0, weight.type="equal", tuning.par=NULL){
  # Check if the turning parameter is given when the weight.type is not "equal"
  if (weight.type != "equal") {
    if (is.null(tuning.par)) {
      stop(paste("The turning parameter cannot be null when the weight.type is ", weight.type, sep = ""))
    }
  }
  # Number of subjects
  l<-ncol(dat) 
  # Test length
  n<-nrow(dat) 
  # Number of threshold parameters
  nthresh<-ncol(b) 

  # Array to store threshold values of interest for each response
  dat.b<-array(dim = list(n, l, 2))
  
  # Arrays/matrices to store final values of theta estimates, standard errors, nonconvergence rates, and the residuals and thetas at each iteration
  theta.est2<- standard.error<- matrix(data=NA, nrow=l)
  convergence<-matrix(0, nrow=l)
  theta.progression<-matrix(NA, nrow = l, ncol = iter)
  residual<-matrix(data=NA, nrow = n, ncol = l)

  # Select threshold values (b) based on response categories
  # b_{jk} and b_{j,k+1} are used to calculate the probability of the response k on item jfor (i in 1:l){
    for(j in 1:n){
      if(dat[j,i]==1){
        # For responses in lowest category
        dat.b[j,i,1]<- -100000 # very small b_{jk} will give probability 1
        dat.b[j,i,2]<-b[j,1]
      }else if(dat[j,i]>nthresh){
        # For responses in highest category
        dat.b[j,i,1]<-b[j,dat[j,i]-1]
        dat.b[j,i,2]<-100000 # very large b_{j,k+1} will give probability 0
      }else{
        # For responses in middle categories
        dat.b[j,i,1]<-b[j,dat[j,i]-1]
        dat.b[j,i,2]<-b[j,dat[j,i]]
      }
    }
  }

  # Loop to estimate theta for each subject 
  for(i in 1:l){
    # Initialize theta value
    theta<-ifelse(length(init.val)>1, init.val[i], init.val)
    P0<-0

    # Iterative loop for maximum likelihood estimation of theta
    for (k in 1:iter){

      # Compute item response probability for the response category and next category
      exponent_0<-a*(theta-dat.b[,i,1])
      exponent_1<-a*(theta-dat.b[,i,2])
      
      # Calculate P*_k and P*_{k+1}
      ps0<-1/(1+exp(-1.7*exponent_0))
      ps1<-1/(1+exp(-1.7*exponent_1))
      qs0<-1-ps0
      qs1<-1-ps1
      
      # Compute P_k, the probability of response k
      P<-ps0-ps1

      # Calculate expected response
      pstar<-1/(1+exp(-1.7*(a*(theta-b))))
      probs<-pstar_to_p(pstar)
      expected.value<-probs%*%matrix(c(1:(nthresh+1))) 

      # Calculate standardized residual
      residual[,i]<-(dat[,i]-expected.value)/sqrt(rowSums(apply(matrix(1:(nthresh+1)), 1, function(x) x-expected.value)^2*probs))  
      
      # Compute item response weights based on specified weight function (bisquare, Huber, equal)
      weighting.term <- NULL
      if (weight.type == "bisquare") {
        weighting.term <- bisquare(residual[,i], tuning.par)
      } else if (weight.type == "Huber") {
        weighting.term <- huber(residual[,i], tuning.par)
      } else {
        weighting.term <- 1
      }
      # Check if weighting term is determined
      if (is.null(weighting.term)) {
        stop("Cannot determine the weighting function.")
      }
      
      # First and second derivatives of the log-likelihood
      D1<-sum(1.7*a*weighting.term*(ps0*qs0-ps1*qs1)/P) 
      D2<-sum(1.7^2*a^2*weighting.term*( (ps0*qs0*(qs0-ps0)-ps1*qs1*(qs1-ps1))/P - (ps0*qs0-ps1*qs1)^2/P^2 )) 
      
      # Check for NAs & record nonconvergence
      if(is.na(theta-D1/D2)){
        theta.est2[i]<-theta<-NA
        convergence[i,1]<-1
        break
      }

      # Update and store theta for this iteration
      theta<-theta.progression[i,k]<-theta-D1/D2
      
      # Stop Newton-Raphson method if log-likelihood difference is converged / less than cutoff
      log_like<-sum(log(P))-sum(log(P0))
      if(abs(log_like)<cutoff){
        break
      }

      # Update probability for loglikelihood comparison
      P0<-P
    }

    # Store final theta estimate for subject
    theta.est2[i]<-theta
    
    # Compute standard error via test information
    pstar<-1/(1+exp(-1.7*(a*(theta-b))))
    probs<-pstar_to_p(pstar)
    test.info<-sum(a^2*probs*(1-probs))
    standard.error[i]<-ifelse(is.na(1/sqrt(test.info)), NA, 1/sqrt(test.info))

    # Handle cases where theta didn't converge withing number of iterations
    if(k==iter){
      theta.est2[i]<-standard.error[i]<-NA
      convergence[i,1]<-1
    }else if(!is.na(theta) & theta< -3){ 
      # Replace thetas that converged outside [-3, 3]
      theta<--3
      pstar<-1/(1+exp(-1.7*(a*(theta-b))))
      probs<-pstar_to_p(pstar)
      test.info<-sum(a^2*probs*(1-probs))
      standard.error[i]<-ifelse(is.na(1/sqrt(test.info)), NA, 1/sqrt(test.info))
      theta.est2[i]<-theta
    }else if(!is.na(theta) & theta>3){
      theta<-3
      pstar<-1/(1+exp(-1.7*(a*(theta-b))))
      probs<-pstar_to_p(pstar)
      test.info<-sum(a^2*probs*(1-probs))
      standard.error[i]<-ifelse(is.na(1/sqrt(test.info)), NA, 1/sqrt(test.info))
      theta.est2[i]<-theta
    }
  }
  return(list(theta = theta.est2, convergence=convergence, standard.error = standard.error, theta.progression = theta.progression, residual=residual))
}

#' Plot histogram of residuals along plot of weight (dependent on TuCo) vs residuals
#'
#' This function plots a histogram of residuals along the graph of the weighting function (dependent on the tuning parameter) as a function of the residual
#' @param r A vector of residuals
#' @param H Huber tuning parameter
#' @param B Bisquare tuning parameter
#' @details The goal of this plot is to visualize the proportion of residuals that are downweighted based on the tuning parameter and allow the researcher to choose a tuning parameter that suits their data well.
#               For a set of residuals with larger variance, a larger tuning parameter should be used.
#               Generally, the tail end of the weighting function should approach the tail end of the distribution of residuals.
#               To increase the downweighting applied in estimation, use a smaller tuning parameter. To decrease the amount of downweighting, use a greater tuning parameter.
#               The function will plot the histogram of residuals below (1) the Huber weight curve if \emph{H} is supplied to the function, (2) the bisquare weight curve if \emph{B} is supplied, or (3) both the Huber and bisquare weight curves if both tuning parameters are supplied.
#               If no tuning parameter is supplied, just the histogram of residuals is generated.
#' @return Histogram plot of residuals beneath a graph of the weight functions vs. the residuals
#' @export
#' @examples
#' ########## choose.tuco example - Unidimensional IRT ######
#' n=40
#' # Generate real thetas
#' thetas<-matrix(seq(0,2, by=.05), ncol=1)
#' # Set item slope and intercept
#' a<-matrix(runif(n, .5, 1.5), ncol=1) 
#' d<-rnorm(n)
#'
#' # Introduce response disturbances: working at a suboptimal level (theta minus 1 standard deviation), for last 40% of items
#' theta.drop<-1
#' chng.pt<-0.6
#' probs<-rbind(apply(thetas, 1, function(x) probs.gen(x, matrix(a[1:(chng.pt*n)], ncol=1), d[1:(chng.pt*n)])$P), apply(thetas-theta.drop, 1, function(x) probs.gen(x, matrix(a[(chng.pt*n+1):n,], ncol=1), d[(chng.pt*n+1):n])$P))
#' dat<-apply(probs, c(1, 2), function(x) rbinom(1, 1, x))
#' 
#' Estimate thetas
#' example<-theta.est(dat, a, d, iter=30, cutoff=.01, init.val=rep(0,ncol(a)), weight.type="equal", tuning.par=NULL)
#' choose.tuco(r=matrix(na.omit(example$residual), ncol=1), B=4)
#'
#' ########### choose.tuco example - GRM ######
#' n=40
#' nthresh<-4
#' # Generate real thetas
#' thetas<-seq(-2,2.1, by=.1)
#'
#' # Set item slope
#' a<-runif(n, .90, 2.15) 
#' # Set category threshold parameters
#' b<- matrix(runif(n*nthresh, -2.5,2.5), nrow = n, ncol =nthresh)
#' b<-t(apply(b, 1, sort)) 
#' 
#' # Calculate response probabilities and generate data
#' probs<-probs.gen.grm(thetas, a, b)
#' dat<-data.gen(probs$P)
#' 
#' Introduce response disturbance: random guessing for latter 40% of the exam
#' abdat<-dat
#' chng.pt<-.6 
#' abdat[(chng.pt*n+1):n, ]<-sample(c(1:(nthresh+1)), length(thetas)*(n-chng.pt*n), replace = T)
#' Calculate ability estimates and residuals
#' mle<-theta.est.grm(dat, a, b, iter=30, cutoff=0.01, init.val=0, weight.type="equal")
#' choose.tuco(matrix(mle$residual), H=.1, B=.8)
#'
#' ######## choose.tuco example - MIRT ##########
#' data(SAT12)
#' SAT12[SAT12 == 8] <- NA #set 8 as a missing value
#'
#' # Correct answer key
#' key <- c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5)
#' scoredSAT12 <- key2binary(SAT12, key)
#' specific <- c(2, 3, 2, 3, 3, 2, 1, 2, 1, 1, 1, 3, 1, 3, 1, 2, 1, 1, 3, 3, 1, 1, 3, 1, 3, 3, 1, 3, 2, 3, 1,2) #which factor each item loads on
#' b_mod1 <- mirt(scoredSAT12, specific)
#' ipars<-matrix(unlist(coef(b_mod1))[1:(32*6)], nrow = length(key), byrow=T) #item parameters
#'
#' ## Set Parameters
#' a <- ipars[,1:3]
#' d<- ipars[,4]
#' # Remove vectors with missing data
#' dat<-scoredSAT12[!is.na(rowSums(scoredSAT12)),] 
#' colnames(dat)<-NULL
#'
#' # Calculate theta estimates and residuals
#' out<-theta.est(t(dat), a, d, iter=30, cutoff=.01, weight.type="equal")
#' choose.tuco(matrix(out$residual[,,2]), H=1, B=4)

choose.tuco<-function(r, H=NULL, B=NULL){
  # r is a vector of residuals
  residuals<-data.frame(Residual =r)
  hist.out<-ggplot(residuals, aes(x=Residual))+geom_histogram(aes(y = ..density..), bins=50)+ ylab("Density")
  if(!is.null(H) & !is.null(B)){
    weight.out<-ggplot()+stat_function(fun=function(x) huber(x, H), aes(colour = "Huber"))+
      stat_function(fun=function(x) bisquare(x, B),  aes(colour = "Bisquare"))+
      xlim(min(r, na.rm =T), max(r, na.rm =T)) + ylab("Weight")+
      scale_color_manual(name = "Function", breaks=c('Bisquare', 'Huber'),
                         values=c('Bisquare'="darkcyan", 'Huber'='firebrick')) +theme(legend.position = c(.9, .74))+
      ggtitle("Weights Applied in Estimation")
    return(do.call(ggarrange, c(list(weight.out, hist.out+ggtitle("Histogram of Residuals")), ncol = 1, nrow = 2)))
  }else if(is.null(H) & !is.null(B)){
    weight.out<-ggplot()+stat_function(fun=function(x) bisquare(x, B), aes(colour = "Bisquare"))+
      xlim(min(r, na.rm =T), max(r, na.rm =T)) + ylab("Weight")+
      scale_color_manual(name = "Function", breaks=c('Bisquare'),
                         values=c('Bisquare'="darkcyan")) +theme(legend.position = c(.9, .74))+
      ggtitle("Weights Applied in Estimation")
    return(do.call(ggarrange, c(list(weight.out, hist.out+ggtitle("Histogram of Residuals")), ncol = 1, nrow = 2)))
  }else if(is.null(B) & !is.null(H)){
    weight.out<-ggplot()+stat_function(fun=function(x) huber(x, H), aes(colour = "Huber"))+
      xlim(min(r, na.rm =T), max(r, na.rm =T)) + ylab("Weight")+
      scale_color_manual(name = "Function", breaks=c('Huber'),
                         values=c('Huber'='firebrick')) +theme(legend.position = c(.9, .74))+
      ggtitle("Weights Applied in Estimation")
    return(do.call(ggarrange, c(list(weight.out, hist.out+ggtitle("Histogram of Residuals")), ncol = 1, nrow = 2)))
  }else {
    return(hist.out+ggtitle("Histogram of Residuals"))
  }
}

#' Plot to compare robust estimates with MLE
#'
#' This function plots robust estimates vs the MLE
#' @param dat \eqn{n \times m} matrix of response data for \emph{n} items and \emph{m} subjects
#' @param a \eqn{n \times p} matrix of slope parameters for \emph{n} items and \emph{p} dimensions (\emph{p=1} if using the GRM or unidimensional 2PL model)
#' @param b If type = “GRM”, an \eqn{n \times (K-1)} matrix of intercept parameters
#' @param d If type = “MIRT”, a vector of discrimination parameters for \emph{n} items
#' @param iter Maximum number of iterations. Default is 100
#' @param cutoff Threshold value to terminate the iteration when the likelihood changes below this value, which means that the estimation is converged.
#' @param H Huber tuning parameter
#' @param B Bisquare tuning parameter
#' @param same.plot.dim If TRUE and type = “MIRT”, estimates across all \emph{p} dimensions will be plotted on the same graph. If FALSE (default) and type = “MIRT”, one plot per dimension will be generated.
#' @param same.plot If TRUE (default) and both \emph{H} and \emph{B} are supplied, generates both the Huber and bisquare plots in the same image frame. If FALSE, the Huber and bisquare plots are generated on separate images.
#' @param type MIRT or GRM
#' @details When the data is not disturbed, robust estimates should not differ greatly from the maximum likelihood estimate (MLE).
#                                       By plotting the robust estimates against the MLE, the user can identify possible aberrant trends, if the robust estimates are far from the MLE, as indicated by the distance from the \eqn{y=x} identity line.
#                                       Larger discrepancies between the point plotted for a subject and the identity line suggest there may be some disturbance in this subject’s data that the robust estimation may be correcting.
#                                       At least one tuning parameter \emph{H} or \emph{B} must be supplied to the function; if both are supplied, the function will return separate plots for both weighting systems.

#' @return ‘Summary Statistics (Huber)’ If \emph{H} is supplied, a dataframe where each row provides a subject’s ID, MLE, the Huber-weighted robust estimate, the minimum distance between the point on the plot and the identity line, and their response vector. The subjects are organized by greatest to least distance.
#' @return ‘Summary Statistics (Bisquare)’ If \emph{B} is supplied, a dataframe where each row provides a subject’s ID, MLE, the bisquare-weighted robust estimate, the minimum distance between the point on the plot and the identity line, and their response vector. The subjects are organized by greatest to least distance.
#' @return Plots If same.plot = TRUE and both \emph{H} and \emph{B} are supplied, each robust ability estimate is plotted against the MLE; graphs for each of the Huber- and bisquare-weighted estimates are generated separately but on the same image frame. The identity line \eqn{y=x} is plotted as a reference line.
#' @return `Huber Plot` If same.plot = FALSE or \emph{B} is not supplied, each Huber-weighted robust ability estimate is plotted against the MLE with the identity line \eqn{y=x} as reference.
#' @return `Bisquare Plot` If same.plot = FALSE or \emph{H} is not supplied, each bisquare-weighted robust ability estimate is plotted against the MLE with the identity line \eqn{y=x} as reference.
#' @export
#' @example
#' ## Test length
#' n <- 30
#'
#' ## Number of iterations of newton's method
#' iter <- 15
#'
#' ## Number of thresholds
#' nthresh <- 4
#'
#' ## Generate real thetas
#' thetas <- seq(-2, 2, by=.1)
#'
#' ## Generate item slope
#' a <- runif(n, .90, 2.15)
#'
#' ## Generate category threshold parameters
#' b <- matrix(runif(n*nthresh, -2.5,2.5), nrow = n, ncol =nthresh)
#' b <- t(apply(b, 1, sort))
#'
#' ## Calculate probabilities
#' probs <- probs.gen.grm(thetas, a, b)
#'
#' ## Generate input data from probabilities
#' abdat <- data.gen(probs$P)
#'
#' ## Introduce aberrant responses: random guessing for latter 40% of the exam
#' chng.pt <- .6
#' abdat[(chng.pt*n+1):n, ] <- sample(c(1:(nthresh+1)), length(thetas)*(n-chng.pt*n), replace = T)
#'
#' ## Plot the GRM
#' out<-theta_plots(abdat, a, b=b, iter=30, cutoff=0.01, H=.1, B=1, same.plot = F, type="GRM")
#' ## Check bisquare plot
#' out$`Bisquare Plot`
#' ## Check Huber summary
#' out$`Summary Statistics (Huber)`
theta_plots<-function(dat, a, d=NULL, b=NULL, iter=30, cutoff=0.01, H=NULL, B=NULL, same.plot.dim = F, same.plot = T, type){
  if(type != "Dichotomous" & type != "GRM"){
    return(print("Please enter a valid type of model (e.g., 'Dichotomous' or 'GRM')."))
  }else if(type == "Dichotomous"){
    if(is.null(d)){ return(print("Please enter a vector of intercept values (d)."))}
    h.plots<-b.plots<-list()
    dim<-ncol(a)
    if(!same.plot.dim){
      theta_estimate=theta.est(dat, a,d, iter, cutoff, init.val=rep(0,ncol(a)), weight.type = "equal")$Theta
      dim<-ncol(theta_estimate) #number of dimensions
      n<-nrow(theta_estimate) #number of subjects

      if(!is.null(H) & !is.null(B)){
        huber_theta_estimate=theta.est(dat, a,d, iter, cutoff, init.val=rep(0,ncol(a)), weight.type = "Huber", tuning.par = H)$Theta
        bisquare_theta_estimate=theta.est(dat, a,d, iter, cutoff, init.val=rep(0,ncol(a)), weight.type = "bisquare", tuning.par = B)$Theta
        pnt.h<-matrix(apply(cbind(matrix(theta_estimate), matrix(huber_theta_estimate)), 1, function(x) sum(x)/2), ncol = dim)
        pnt.b<-matrix(apply(cbind(matrix(theta_estimate), matrix(bisquare_theta_estimate)), 1, function(x) sum(x)/2), ncol = dim)
        Distance.h = sqrt((theta_estimate-pnt.h)^2+(huber_theta_estimate-pnt.h)^2)
        Distance.b = sqrt((theta_estimate-pnt.b)^2+(bisquare_theta_estimate-pnt.b)^2)
        stats.h<-data.frame(Dis=apply(Distance.h, 1, function(x) mean(x, na.rm=T)),
                                     ID = 1:nrow(theta_estimate))
        stats.b<-data.frame(Dis=apply(Distance.b, 1, function(x) mean(x, na.rm=T)),
                            ID = 1:nrow(theta_estimate))
        h.names<-b.names<-c("Dis", "ID")
        for(i in 1:dim){
          stats.h<-cbind(stats.h, MLE = theta_estimate[,i],
                Huber = huber_theta_estimate[,i],
                Distance =Distance.h[,i])
          stats.b<-cbind(stats.b, MLE = theta_estimate[,i],
                Bisquare = bisquare_theta_estimate[,i],
                Distance =Distance.b[,i])
          h.names<-c(h.names, paste0("MLE", i), paste0("Huber", i), paste0("Distance", i))
          b.names<-c(b.names, paste0("MLE", i), paste0("Bisquare", i), paste0("Distance", i))
        }
        colnames(stats.h)<-h.names
        colnames(stats.b)<-b.names
        sum.stats.h<-cbind(stats.h, t(dat))%>%arrange(desc(Dis))
        sum.stats.b<-cbind(stats.b, t(dat))%>%arrange(desc(Dis))

        for(i in 1:dim){
          #message(i)
          h.plots[[i]] <- local({
            i <- i
            huberplot<- ggplot(mapping = aes (x = theta_estimate[,i], y = huber_theta_estimate[,i]))+ geom_abline(color = "red", slope = 1) +
              geom_point() + xlab(bquote(hat(theta)[MLE])) + ylab(bquote(hat(theta)[Huber] ~ " " ~ (H== .(H) ))) +
              ggtitle(bquote("Huber-Weighted Robust Estimates of " ~ theta ~ " vs. the MLE, Dimension" ~ .(i) ))
          })
          b.plots[[i]] <- local({
            i <- i
            bisquareplot<- ggplot(mapping = aes (x = theta_estimate[,i], y = bisquare_theta_estimate[,i]))+ geom_abline(color = "red", slope = 1) +
              geom_point()+ xlab(bquote(hat(theta)[MLE])) + ylab(bquote(hat(theta)[Bisquare] ~ " " ~ (B== .(B) ))) +
              ggtitle(bquote("Bisquare-Weighted Robust Estimates of " ~ theta ~ " vs. the MLE, Dimension" ~ .(i) ))
          })
        }
        return(list("Summary Statistics (Huber)" = sum.stats.h[,-1], "Summary Statistics (Bisquare)" = sum.stats.b[,-1], "Huber Plot" = do.call(ggarrange, c(h.plots, ncol = 1, nrow = dim, common.legend = T)), "Bisquare Plot" = do.call(ggarrange, c(b.plots, ncol = 1, nrow = dim, common.legend = T))))
      }else if(is.null(H) & !is.null(B)){
        bisquare_theta_estimate=theta.est(dat, a,d, iter, cutoff, init.val=rep(0,ncol(a)), weight.type = "bisquare", tuning.par = B)$Theta

        pnt.b<-matrix(apply(cbind(matrix(theta_estimate), matrix(bisquare_theta_estimate)), 1, function(x) sum(x)/2), ncol = dim)
        Distance.b = sqrt((theta_estimate-pnt.b)^2+(bisquare_theta_estimate-pnt.b)^2)
        stats.b<-data.frame(Dis=apply(Distance.b, 1, function(x) mean(x, na.rm=T)),
                            ID = 1:nrow(theta_estimate))
        b.names<-c("Dis", "ID")
        for(i in 1:dim){
          stats.b<-cbind(stats.b, MLE = theta_estimate[,i],
                         Bisquare = bisquare_theta_estimate[,i],
                         Distance =Distance.b[,i])
          b.names<-c(b.names, paste0("MLE", i), paste0("Bisquare", i), paste0("Distance", i))
        }
        colnames(stats.b)<-b.names
        sum.stats.b<-cbind(stats.b, t(dat))%>%arrange(desc(Dis))

        for(i in 1:dim){
          b.plots[[i]] <- local({
            i <- i
            bisquareplot<- ggplot(mapping = aes (x = theta_estimate[,i], y = bisquare_theta_estimate[,i]))+ geom_abline(color = "red", slope = 1) +
              geom_point()+ xlab(bquote(hat(theta)[MLE])) + ylab(bquote(hat(theta)[Bisquare] ~ " " ~ (B== .(B) ))) +
              ggtitle(bquote("Bisquare-Weighted Robust Estimates of " ~ theta ~ " vs. the MLE, Dimension" ~ .(i) ))
          })
        }
        return(list("Summary Statistics" = sum.stats.b[,-1], "Bisquare Plots" =do.call(ggarrange, c(b.plots, ncol = 1, nrow = dim, common.legend = T))))
      }else if(!is.null(H) & is.null(B)){
        huber_theta_estimate=theta.est(dat, a,d, iter, cutoff, init.val=rep(0,ncol(a)), weight.type = "Huber", tuning.par = H)$Theta

        pnt.h<-matrix(apply(cbind(matrix(theta_estimate), matrix(huber_theta_estimate)), 1, function(x) sum(x)/2), ncol = dim)
        Distance.h = sqrt((theta_estimate-pnt.h)^2+(huber_theta_estimate-pnt.h)^2)
        stats.h<-data.frame(Dis=apply(Distance.h, 1, function(x) mean(x, na.rm=T)),
                            ID = 1:nrow(theta_estimate))
        h.names<-c("Dis", "ID")
        for(i in 1:dim){
          stats.h<-cbind(stats.h, MLE = theta_estimate[,i],
                         Huber = huber_theta_estimate[,i],
                         Distance =Distance.h[,i])
          h.names<-c(h.names, paste0("MLE", i), paste0("Huber", i), paste0("Distance", i))
        }
        colnames(stats.h)<-h.names
        sum.stats.h<-cbind(stats.h, t(dat))%>%arrange(desc(Dis))

        for(i in 1:dim){
          h.plots[[i]] <- local({
            i <- i
            huberplot<- ggplot(mapping = aes (x = theta_estimate[,i], y = huber_theta_estimate[,i]))+ geom_abline(color = "red", slope = 1) +
              geom_point() + xlab(bquote(hat(theta)[MLE])) + ylab(bquote(hat(theta)[Huber] ~ " " ~ (H== .(H) ))) +
              ggtitle(bquote("Huber-Weighted Robust Estimates of " ~ theta ~ " vs. the MLE, Dimension" ~ .(i) ))
          })

        }
        return(list("Summary Statistics (Huber)" = sum.stats.h[,-1], "Huber Plots" = do.call(ggarrange, c(h.plots, ncol = 1, nrow = dim, common.legend = T))))
      }else{ return(print("A valid tuning parameter is needed."))}
    }else{ # if not same plot
      theta_estimate=theta.est(dat, a,d, iter, cutoff, init.val=rep(0,ncol(a)), weight.type = "equal")
      dim<-ncol(theta_estimate$Theta) #number of dimensions
      n<-nrow(theta_estimate$Theta) #number of subjects

      if(!is.null(H) & !is.null(B)){
        huber_theta_estimate=theta.est(dat, a,d, iter, cutoff, init.val=rep(0,ncol(a)), weight.type = "Huber", tuning.par = H)
        bisquare_theta_estimate=theta.est(dat, a,d, iter, cutoff, init.val=rep(0,ncol(a)), weight.type = "bisquare", tuning.par = B)

        pnt.h<-matrix(apply(cbind(matrix(theta_estimate), matrix(huber_theta_estimate)), 1, function(x) sum(x)/2), ncol = dim)
        pnt.b<-matrix(apply(cbind(matrix(theta_estimate), matrix(bisquare_theta_estimate)), 1, function(x) sum(x)/2), ncol = dim)
        Distance.h = sqrt((theta_estimate-pnt.h)^2+(huber_theta_estimate-pnt.h)^2)
        Distance.b = sqrt((theta_estimate-pnt.b)^2+(bisquare_theta_estimate-pnt.b)^2)
        stats.h<-data.frame(Dis=apply(Distance.h, 1, function(x) mean(x, na.rm=T)),
                            ID = 1:nrow(theta_estimate))
        stats.b<-data.frame(Dis=apply(Distance.b, 1, function(x) mean(x, na.rm=T)),
                            ID = 1:nrow(theta_estimate))
        h.names<-b.names<-c("Dis", "ID")
        for(i in 1:dim){
          stats.h<-cbind(stats.h, MLE = theta_estimate[,i],
                         Huber = huber_theta_estimate[,i],
                         Distance =Distance.h[,i])
          stats.b<-cbind(stats.b, MLE = theta_estimate[,i],
                         Bisquare = bisquare_theta_estimate[,i],
                         Distance =Distance.b[,i])
          h.names<-c(h.names, paste0("MLE", i), paste0("Huber", i), paste0("Distance", i))
          b.names<-c(b.names, paste0("MLE", i), paste0("Bisquare", i), paste0("Distance", i))
        }
        colnames(stats.h)<-h.names
        colnames(stats.b)<-b.names
        sum.stats.h<-cbind(stats.h, t(dat))%>%arrange(desc(Dis))
        sum.stats.b<-cbind(stats.b, t(dat))%>%arrange(desc(Dis))

        dat<-data.frame(MLE=matrix(theta_estimate$Theta),
                        Huber=matrix(huber_theta_estimate$Theta),
                        Bisquare=matrix(bisquare_theta_estimate$Theta),
                        Dimension=as.factor(rep(1:dim, each=n)))

        huberplot<- ggplot(data = dat, aes(x=MLE, y=Huber)) + geom_abline(color = "red", slope = 1) +
          geom_point(aes(color=Dimension)) + xlab(bquote(hat(theta)[MLE])) +
          ylab(bquote(hat(theta)[Huber] ~ " " ~ (H== .(H) ))) + ggtitle(bquote("Huber-Weighted Robust Estimates of " ~ theta ~ " vs. the MLE"))

        bisquareplot<- ggplot(data = dat, aes(x=MLE, y=Bisquare)) + geom_abline(color = "red", slope = 1) +
          geom_point(aes(color=Dimension)) + xlab(bquote(hat(theta)[MLE])) +
          ylab(bquote(hat(theta)[Bisquare] ~ " " ~ (B== .(B) ))) + ggtitle(bquote("Bisquare-Weighted Robust Estimates of " ~ theta ~ " vs. the MLE"))

        return(list("Summary Statistics (Huber)" = sum.stats.h[,-1], "Summary Statistics (Bisquare)" = sum.stats.b, "Huber Plots" = huberplot, "Bisquare Plots" = bisquareplot))
      }else if(is.null(H) & !is.null(B)){
        bisquare_theta_estimate=theta.est(dat, a,d, iter, cutoff, init.val=rep(0,ncol(a)), weight.type = "bisquare", tuning.par = B)

        pnt.b<-matrix(apply(cbind(matrix(theta_estimate), matrix(bisquare_theta_estimate)), 1, function(x) sum(x)/2), ncol = dim)
        Distance.b = sqrt((theta_estimate-pnt.b)^2+(bisquare_theta_estimate-pnt.b)^2)
        stats.b<-data.frame(Dis=apply(Distance.b, 1, function(x) mean(x, na.rm=T)),
                            ID = 1:nrow(theta_estimate))
        b.names<-c("Dis", "ID")
        for(i in 1:dim){
          stats.b<-cbind(stats.b, MLE = theta_estimate[,i],
                         Bisquare = bisquare_theta_estimate[,i],
                         Distance =Distance.b[,i])
          b.names<-c(b.names, paste0("MLE", i), paste0("Bisquare", i), paste0("Distance", i))
        }
        colnames(stats.b)<-b.names
        sum.stats.b<-cbind(stats.b, t(dat))%>%arrange(desc(Dis))


        dat<-data.frame(MLE=matrix(theta_estimate$Theta),
                        Bisquare=matrix(bisquare_theta_estimate$Theta),
                        Dimension=as.factor(rep(1:dim, each=n)))
        bisquareplot<- ggplot(data = dat, aes(x=MLE, y=Bisquare)) + geom_abline(color = "red", slope = 1) +
          geom_point(aes(color=Dimension)) + xlab(bquote(hat(theta)[MLE])) +
          ylab(bquote(hat(theta)[Bisquare] ~ " " ~ (B== .(B) ))) + ggtitle(bquote("Bisquare-Weighted Robust Estimates of " ~ theta ~ " vs. the MLE"))
        return(list("Sumary Statistics (Bisquare)" = sum.stats.b[,-1], "Bisquare Plot" = bisquareplot))
      }else if(!is.null(H) & is.null(B)){
        huber_theta_estimate=theta.est(dat, a,d, iter, cutoff, init.val=rep(0,ncol(a)), weight.type = "Huber", tuning.par = H)

        pnt.h<-matrix(apply(cbind(matrix(theta_estimate), matrix(huber_theta_estimate)), 1, function(x) sum(x)/2), ncol = dim)
        Distance.h = sqrt((theta_estimate-pnt.h)^2+(huber_theta_estimate-pnt.h)^2)
        stats.h<-data.frame(Dis=apply(Distance.h, 1, function(x) mean(x, na.rm=T)),
                            ID = 1:nrow(theta_estimate))
        h.names<-c("Dis", "ID")
        for(i in 1:dim){
          stats.h<-cbind(stats.h, MLE = theta_estimate[,i],
                         Huber = huber_theta_estimate[,i],
                         Distance =Distance.h[,i])
          h.names<-c(h.names, paste0("MLE", i), paste0("Huber", i), paste0("Distance", i))
        }
        colnames(stats.h)<-h.names
        sum.stats.h<-cbind(stats.h, t(dat))%>%arrange(desc(Dis))

        dat<-data.frame(MLE=matrix(theta_estimate$Theta),
                        Huber=matrix(huber_theta_estimate$Theta),
                        Dimension=as.factor(rep(1:dim, each=n)))

        huberplot<- ggplot(data = dat, aes(x=MLE, y=Huber)) + geom_abline(color = "red", slope = 1) +
          geom_point(aes(color=Dimension)) + xlab(bquote(hat(theta)[MLE])) +
          ylab(bquote(hat(theta)[Huber] ~ " " ~ (H== .(H) ))) + ggtitle(bquote("Huber-Weighted Robust Estimates of " ~ theta ~ " vs. the MLE"))

        return(list("Summary Statistics (Huber)" = sum.stats.h[,-1], "Huber Plots" = huberplot))
      }else{return(print("A valid tuning parameter is needed."))}
    }
  }else if(type == "GRM"){

    theta_estimate=theta.est.grm(dat, a,b, iter, cutoff, 0, weight.type="equal")$theta
    if(!is.null(H)& !is.null(B)){
      huber_theta_estimate=theta.est.grm(dat, a,b, iter, cutoff, 0, weight.type="Huber", tuning.par=H)$theta
      bisquare_theta_estimate=theta.est.grm(dat, a,b, iter, cutoff, 0, weight.type="bisquare", tuning.par=B)$theta

      h.plot<- ggplot(mapping = aes (x = theta_estimate, y = huber_theta_estimate))+ geom_abline(color = "red", slope = 1) +
        geom_point() + xlab(bquote(hat(theta)[MLE])) + ylab(bquote(hat(theta)[Huber] ~ " " ~ (H== .(H) ))) +
        ggtitle(bquote("Huber-Weighted Robust Estimates of " ~ theta ~ " vs. the MLE" ))

      b.plot<- ggplot(mapping = aes (x = theta_estimate, y = bisquare_theta_estimate))+ geom_abline(color = "red", slope = 1) +
        geom_point()+ xlab(bquote(hat(theta)[MLE])) + ylab(bquote(hat(theta)[Bisquare] ~ " " ~ (B== .(B) ))) +
        ggtitle(bquote("Bisquare-Weighted Robust Estimates of " ~ theta ~ " vs. the MLE" ))

      pnt.h<-apply(cbind(theta_estimate, huber_theta_estimate), 1, function(x) sum(x)/2)
      pnt.b<-apply(cbind(theta_estimate, bisquare_theta_estimate), 1, function(x) sum(x)/2)
      sum.stats.h<-cbind(data.frame(ID = 1:nrow(theta_estimate),
                                    MLE = theta_estimate,
                                    Huber = huber_theta_estimate,
                                    Distance = sqrt((theta_estimate-pnt.h)^2+(huber_theta_estimate-pnt.h)^2)), t(dat))%>%arrange(desc(Distance))
      sum.stats.b<-cbind(data.frame(ID = 1:nrow(theta_estimate),
                                    MLE = theta_estimate,
                                    Bisquare = bisquare_theta_estimate,
                                    Distance = sqrt((theta_estimate-pnt.b)^2+(bisquare_theta_estimate-pnt.b)^2)), t(dat)) %>%arrange(desc(Distance))
      if(same.plot){ #allows user to have sperate plots or see both Huber and Bisquare at same time
        return(list("Summary Statistics (Huber)" = sum.stats.h, "Summary Statistics (Bisquare)" = sum.stats.b, "Plots" = do.call(ggarrange, c(list(h.plot, b.plot), ncol = 1, nrow = 2))))
      }else{
        return(list("Summary Statistics (Huber)" = sum.stats.h, "Summary Statistics (Bisquare)" = sum.stats.b, "Huber Plot" = h.plot, "Bisquare Plot" =b.plot))
      }
    }else if(!is.null(H)& is.null(B)){
      huber_theta_estimate=theta.est.grm(dat, a,b, iter, cutoff, 0, weight.type="Huber", tuning.par=H)$theta

      h.plot<- ggplot(mapping = aes (x = theta_estimate, y = huber_theta_estimate))+ geom_abline(color = "red", slope = 1) +
        geom_point() + xlab(bquote(hat(theta)[MLE])) + ylab(bquote(hat(theta)[Huber] ~ " " ~ (H== .(H) ))) +
        ggtitle(bquote("Huber-Weighted Robust Estimates of " ~ theta ~ " vs. the MLE" ))

      pnt.h<-apply(cbind(theta_estimate, huber_theta_estimate), 1, function(x) sum(x)/2)
      sum.stats.h<-cbind(data.frame(ID = 1:nrow(theta_estimate),
                                    MLE = theta_estimate,
                                    Huber = huber_theta_estimate,
                                    Distance = sqrt((theta_estimate-pnt.h)^2+(huber_theta_estimate-pnt.h)^2)), t(dat))%>%arrange(desc(Distance))

      return(list("Summary Statistics (Huber)" = sum.stats.h,  "Huber Plot" = h.plot))
    }else if(is.null(H)& !is.null(B)){
      bisquare_theta_estimate=theta.est.grm(dat, a,b, iter, cutoff, 0, weight.type="Huber", tuning.par=H)$theta

      b.plot<- ggplot(mapping = aes (x = theta_estimate, y = bisquare_theta_estimate))+ geom_abline(color = "red", slope = 1) +
        geom_point()+ xlab(bquote(hat(theta)[MLE])) + ylab(bquote(hat(theta)[Bisquare] ~ " " ~ (B== .(B) ))) +
        ggtitle(bquote("Bisquare-Weighted Robust Estimates of " ~ theta ~ " vs. the MLE" ))

      pnt.b<-apply(cbind(theta_estimate, bisquare_theta_estimate), 1, function(x) sum(x)/2)
      sum.stats.b<-cbind(data.frame(ID = 1:nrow(theta_estimate),
                                    MLE = theta_estimate,
                                    Bisquare = bisquare_theta_estimate,
                                    Distance = sqrt((theta_estimate-pnt.b)^2+(bisquare_theta_estimate-pnt.b)^2)), t(dat)) %>%arrange(desc(Distance))
      return(list("Summary Statistics (Bisquare)" = sum.stats.b, "Bisquare Plot" = b.plot))

    }else{return(print("A valid tuning parameter is needed."))}

  }
}



