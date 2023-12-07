## set of test functions
#
# The location of the global optimum, glob_xstar, can be retrieved 
# by calling the function with the argument get_glob_xstar = TRUE and any 
# xx of the right dimension. Returns NA's when glob_xstar is not known.
#


# Ackley function. Global minimum at glob_xstar
ackley <- function(xx, get_glob_xstar=FALSE, a=20, b=0.2, c=2*pi){
  glob_xstar <- rep(0,length(xx))
  if (get_glob_xstar) return(glob_xstar) else {
    xx <- xx - glob_xstar
    aa <- 6.4*xx
    d <- length(aa)
    sum <- -a*exp(-b*sqrt((1/d)*sum(aa*aa))) - exp((1/d)*sum(cos(c*aa)))+ a + exp(1)
    y <- sum 
    return(y)    
  }
}

## Rastrigin function. Global minimum at glob_xstar
rastrigin <- function(xx, get_glob_xstar=FALSE){
  glob_xstar <- rep(0,length(xx))
  if (get_glob_xstar) return(glob_xstar) else {
    scaling_fact <- 1.024
    xx <- xx - glob_xstar
    aa <- scaling_fact*xx
    d <- length(aa)
    sum <- sum(aa*aa - 10*cos(2*pi*aa))
    y <- sum + 10*d
    return(y)
  }
}

# glob_xstar at rep(3.209687,d) where fstar=0
## Schwefel function.
schwefel <- function(xx, get_glob_xstar=FALSE){
  glob_xstar <- rep(3.209687,length(xx))
  if (get_glob_xstar) return(glob_xstar) else {
    xx <- xx + 1
    aa <- 100*xx
    d <- length(aa)
    sum <- 418.9829*d - sum(aa*sin(sqrt(abs(aa))))
    y <- sum
    return(y)
  }
}

## Sphere function. 
sphere <- function(xx, get_glob_xstar=FALSE){
  d <- length(xx)
  glob_xstar<-1:d
  if (get_glob_xstar) return(glob_xstar) else {
    aa <- xx - glob_xstar
    y <- sum(aa*aa)
    return(y)
  }
}

##### Michalewicz function #########
#glob min in 2D: -1.83, 5D: -4.71 , 10D: -9.64
# not sure what glob_xstar is so returning NA's
michalewicz <- function(xx, get_glob_xstar=FALSE, m=10){
  if (get_glob_xstar) return(rep(NA,length(xx))) else {
    aa <- 0.1*pi*(xx + 5)
    #   aa <- xx
    d <- length(aa)
    i <- 1:d
    sum <- sum(sin(aa)*sin((i*aa*aa)/pi)^(2*m))
    y <- -sum
    return(y + 2)
  }
}

##### quadratic function #########"
quadratic <- function(xx,cond.no=20,get_glob_xstar=FALSE){
  glob_xstar <- rep(0,length(xx))
  if (get_glob_xstar) return(glob_xstar) 
  else {
    aa <- 1.024*xx
    d <- length(aa)
    xstar <- glob_xstar
    lambdas <- diag(seq(1, cond.no, ,d))
    # matrix with arbitrary orientation. The seed number decides the orientation.
    if (!exists("glob_umat")) {
      set.seed(1) # change this seed to change the orientation of the quadratic function
      glob_umat <<- qr.Q(qr(matrix(runif(d*d),nrow=d,ncol=d)))
      glob_umat <<- glob_umat[,sample(seq(1,d))]
      #  glob_umat <<- diag(1, d) # to generate a quadratic function whose principal axes are aligned with coordinates
    }
    H <- glob_umat%*%lambdas%*%t(glob_umat)  
    y <- 0.5*t(aa - xstar)%*%H%*%(aa - xstar)
    return(y)    
  }
}

quadratic_ill <- function(xx,get_glob_xstar=FALSE){
  return(quadratic(xx,cond.no=10000,get_glob_xstar = get_glob_xstar))
}

#### Tunnel function #######"
tunnel <- function(xx, get_glob_xstar=FALSE, b=0.5){
  glob_xstar <- rep(0,length(xx))
  if (get_glob_xstar) return(glob_xstar) 
  else {
    d <- length(xx)
    xx <- xx - glob_xstar
    in_tunnel <- TRUE
    if (d > 1) {
      for(i in 2:d){
        if ((xx[i] < -b) | (xx[i] > b)){
          in_tunnel <- FALSE
        }
      }
      if (in_tunnel) {
        y <- -exp(-norm(as.matrix(xx), type="F")/10)
      } else{
        y <- t(xx)%*%xx
      }
    } else{
      y <- xx^2
    }
    return(y)
  }
}

##### function with hierarchical sensitivities and local optima #########
quad_wave <- function(xx, get_glob_xstar=FALSE){
  if (get_glob_xstar) return(rep(NA,length(xx))) 
  else {
    cond.no <- 1.e1
    d <- length(xx)
    lambdas <- diag(seq(1, cond.no, ,d))
    glob_umat <<- diag(1, d) # to generate a quadratic function whose 
    # principal axes are aligned with coordinates. Could be made into a rotating matrix.
    H <- glob_umat%*%lambdas%*%t(glob_umat)  
    y <- 1+0.5*t(xx)%*%H%*%(xx)-0.2*sum(cos(4*pi*xx[c(1,2)]))
    return(y)
  }
}


##### non differentiable and convex L1 norm #########
L1norm <- function(xx, get_glob_xstar=FALSE){
  if (get_glob_xstar) return(rep(0,length(xx))) 
  else {
    y <- sum(abs(xx))
    return(y)
  }
}

rosen <- function(xx, get_glob_xstar=FALSE)
{
  ##########################################################################
  #
  # ROSENBROCK FUNCTION
  # global optimum (glob_xstar) at (1,...,1)
  # copied and slightly changed from
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  ##########################################################################
  if (get_glob_xstar) return(rep(1,length(xx))) 
  else {
    d <- length(xx)
    xi <- xx[1:(d-1)]
    xnext <- xx[2:d]
    y <- sum(100*(xnext-xi^2)^2 + (xi-1)^2)
    return(y)
  }
}

