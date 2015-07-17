\name{PDM-package}
\alias{PDM-package}
\alias{PDM}
\docType{package}
\title{
Principal Direction of Mediation
}
\description{
This package provides functions that calculate the estimates of the Principal Direction of Mediations (PDMs) and corresponding path coefficients, provided treatment (e.g. thermal pain), response (e.g. reported pain), and mediation data (e.g. measurements of fMRI data).
}
\details{
\tabular{ll}{
Package: \tab PDM\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-06-24\cr
Imports:\tab MASS\cr
License: \tab GPL (>= 2)\cr
}

PDM_1 calculates the first PDM and corresponding pathway coefficients, see help(PDM_1) for details; PDM_2 calculates the first two PDMs and corresponding pathway coefficients, see help(PDM_2) for details; and PDM_k calculates the first k PDMs and corresponding pathway coefficients, see help(PDM_k) for details.
}
\author{
Oliver Chén

Maintainer: Oliver Chén <oliver@jhmi.edu>
}
\references{
Chén O, Ogburn E, Crainiceanu C, Caffo B, Lindquist M. Principal Direction of Mediation. To submit.
}

\keyword{ Principal Direction of Mediation; Indirect Effect; Pathway Coefficients }
\seealso{ }

\examples{ 
library(MASS)

k = 10
mu_1 = 10
s_1 = 3
mu_2 = 8
s_2 = 2

Posdef <- function (n, ev = runif(n, 0, 10)) 
{
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  #O <- Q %*% diag(ph)
  
  O <- tcrossprod(Q, t(diag(ph)))
  #Z <- t(O) %*% diag(ev) %*% O
  Z <- crossprod(O,   tcrossprod( diag(ev), t(O)))
   
  return(Z)
}

set.seed(2)

Sigma<- round(Posdef(n=(k+1), ev=1:(k+1)),3)

n=1000 

v=mvrnorm(n, rnorm((k+1), mu_1, s_1), Sigma) 

m1=v[,1:k]
x1=v[, (k+1) ]

l=rnorm(k, mu_2, s_2)

w11=l/sqrt(t(tcrossprod(l, t(l))))

b0=0.4
b1=0.2
r1=0.5

u_m = c(2, 3, 4, 5, 4, 6, 2, 5, 8, 1)

s_mx = Sigma[, (k+1) ][1:k]

u_x=3

s_x=Sigma[(k+1),(k+1)]

#a0=w11%*% (u_m - s_mx %*%solve(s_x)%*%u_x)

a0 = tcrossprod(w11, t(u_m - tcrossprod(s_mx,  t ( tcrossprod(solve(s_x), t(u_x)) )) ) )

#a1=w11%*% (s_mx %*%solve(s_x))

a1 = tcrossprod(w11, t(tcrossprod(s_mx, t(solve(s_x)))))


e1=rnorm(n,0,1)
#y1=b0+r1*x1+b1*m1%*%w11+e1

y1 = b0+r1*x1+b1*tcrossprod(m1, t(w11))


################################################################################################
###  Estimate of the first two PDMs and the corresponding pathway coefficients
################################################################################################

output <- PDM_2(x=x1, y=y1, m=m1, theta_1=rep(1,5), theta_2=rep(1,6), w1=rep(1,length(w11
)), w2=rep(1,length(w11)), 
               imax=100,  tol_1=10^-{7}, tol_2=10^-{3}, tol_3 =20, tol_4=100,
              interval_1 = 10^6, step_1 = 10^4, interval_2 = 10^6,step_2 = 10^4 )

# The first PDM
output$w1

# The second PDM
output$w2

# Pathway coefficient relating to the first PDM
output$theta_1

# Pathway coefficient relating to the second PDM
output$theta_2

               
}
