#' Generate hyperparameters for gibbs sampling.
#'
#' @param vec A named list containing custom hyperparameter values. Any missing components will be filled in with the default value. The components are:
#' \describe{
#' \item{alpha_w}{Shape parameter for p(w_{pk}). Default is 1.}
#' \item{beta_w}{Rate parameter for p(w_{pk}). Default is 5.}
#' \item{c_h}{Shape parameter for non-null scores. Default is 10.}
#' \item{epsilon}{Mean for null scores. Default is 0.25.}
#' \item{a_mu}{Shape parameter for p(mu). Default is 2.}
#' \item{b_mu}{Scale parameter for p(mu). Default is 20.}
#' \item{a_a}{Alpha parameter for p(pi_ps). Default is 1.}
#' \item{b_a}{Beta parameter for p(pi_pis). Default is 1.}
#' \item{beta_0}{Prior mean for intercept of probit stick breaking parameter. Default is 1.5.}
#' \item{tau_0}{Prior precision for probit stick breaking parameters. DEfault is 2.}
#' }
#'
#' @return A list containing hyperparameters to use with the gibbs sampler.
#' @export
#'
#' @examples
#'
#' #Use Default hyperparameters
#' hyperparams <- initialize_zinmf()
#'
#' # Change alpha_w parameter:
#' hyperparams_2 <- initialize_zinmf(list('alpha_w'=2))
#'
initialize_zinmf <- function(vec=list()){
  default_priors = c(
    "alpha_w" = 1,
    'beta_w' = 5,
    "c_h" = 10,
    "epsilon" = 0.25,
    "a_mu" = 2,
    "b_mu" = 20,
    "a_a" = 1,
    "b_a" = 1,
    "beta_0" = 1.5,
    "tau_0" = 2
  )

  out <- vector("numeric")
  for(value in names(default_priors)){
    if(is.null(vec[[value]])){
      out[[value]] <- default_priors[value]
    } else {
      out[[value]] <- vec[[value]]
    }
  }
  out <- c(out)
  return(out)
}

#' Helper function to calculate cosine distance between columns of two matrices.
#'
#' @param x A matrix of dimensions \eqn{P \times N_1}.
#' @param y A matrix of dimensions \eqn{P \times N_2}.
#'
#' @return A matrix containing the pairwise cosine distances between x and y.
#' @export
#'
cosineDist <- function(x,y) {
  (x%*%t(y))/(sqrt(rowSums(x^2) %*% t(rowSums(y^2))))
}

#' Gibbs Sampler for Zero-Inflated Bayesian Infinite Multi-Study NMF
#'
#' The function implements the gibbs sampler described in Hansen et al. (2025+)
#' The code will return samples from the posterior distribution.
#'
#'
#' @param M_s A list of length \eqn{S}, containing counts matrices of dimensions \eqn{P \times N_s}.
#' @param x_s A list of length \eqn{S}, containing covariate matrices of dimensions \eqn{N_s \times D_s}, where \eqn{D_s} is the number of covariates in study s. If no covariates are provided, an intercept will be used for each study.
#' @param K The maximum number of patterns. Default is 10.
#' @param L The maximum number of clusters within pattern scores. Default is 5.
#' @param niter The number of posterior samples. Default is 2000.
#' @param burnin The number of burnin samples. Default is 10000.
#' @param thin Interval between saved samples. Default is 10.
#' @param priors A list containing the hyperparameter values generated from \code{initialize_zinmf}.
#'
#' @return A list containing posterior samples.  The components of the list are:
#' \item{W}{A tensor of dimension \eqn{P \times K \times niter/thin.}}
#' \item{H}{A list of length S, each item containing a tensor of dimension \eqn{K \times N_s \times niter/thin}.}
#' \item{A}{A list of length S, each item containing a tensor of dimension \eqn{P \times N_s \times niter/thin}.}
#' \item{PI}{A list of length S, each item containing a matrix of dimension \eqn{P\times niter/thin}.}
#' \item{D}{A list of length S, each item containing a tensor of dimension \eqn{K \times N_s \times L \times niter/thin}.}
#' \item{Theta}{A tensor of dimension \eqn{K \times L \times niter/thin}.}
#' \item{Beta}{A list of length S, each item containing a tensor of dimension \eqn{K \times L \times Q_s \times niter/thin}.}
#' @importFrom abind abind
#' @importFrom truncnorm rtruncnorm
#' @import LaplacesDemon
#' @import progress
#' @export
#' @examples
#' \dontrun{
#' # Generate Data
#' S <- 3
#' P <- 20
#' K <- 5
#' N_s <- c(50, 60, 70)
#'
#' W <- matrix(rgamma(P*K, shape=1, scale=2), nrow=P, ncol=K)
#' H_s <- lapply(1:S, function(s){matrix(rgamma(K*N_s[s], shape=10, rate=10/20), nrow=K, ncol=N_s[s])})
#' M_s <- lapply(1:S, function(s) matrix(rpois(N_s[s]*P, W %*% H_s[[s]], nrow=K, ncol=N_s[s]))
#'
#' # Run CAVI
#' gibbs_zinmf(M_s=M_s, K=5)
#' }
#'

gibbs_zinmf <- function(M_s,
                        x_s = NULL,
                        K = 10,
                        L = 5,
                        niter = 2000,
                        burnin = 10000,
                        thin = 10,
                        priors = initialize_zinmf()
) {

  # Check dimensions of data and ensure studies contain the same P
  S = length(M_s)
  if(length(unique(sapply(M_s, nrow)))!=1){
    errorCondition("The data must contain the same number of rows in each study!")
  } else {
    P <- unique(sapply(M_s, nrow))
  }
  N_s <- sapply(M_s, ncol)

  if(is.null(x_s)){
    x_s <- lapply(1:S, function(s) matrix(1, ncol=1, nrow=N_s[s]) )
  }
  Q_s <- sapply(x_s, ncol)

  tot_iter <- niter + burnin

  ################################################################################
  # INITIALIZE PARAMETERS
  ################################################################################
  W <- sapply(1:K, function(k){
    rgamma(P,
           priors['alpha_w'],
           priors['beta_w']
    )
  })

  # index beta_{skil} with beta[[s]][k, i, l, ]
  beta <- lapply(1:S, function(s){
    out <- array(dim=c(Q_s[s], L-1, K))
    out[,,] <- sapply(1:K, function(k){
      sapply(1:(L-1), function(l){
        array(rmvn(1, mu=c(priors['beta_0'],rep(0, Q_s[s]-1)), Sigma = priors['tau_0']^(-1)*diag(Q_s[s])), dim=c(Q_s[s],1))
      })
    }, simplify="array")
    aperm(out, c(3, 2, 1))
  })

  omega <- lapply(1:S, function(s){
    sapply(1:L, function(l){
      sapply(1:N_s[s], function(i){
        sapply(1:K, function(k){
          ifelse(l==L, 1, pnorm(x_s[[s]][i,, drop=FALSE] %*% beta[[s]][k,l,], mean=0, sd=1)) * ifelse(l == 1, 1, prod(sapply(1:(l-1), function(r) 1- pnorm(x_s[[s]][i,, drop=FALSE] %*% beta[[s]][k,r,], mean=0, sd=1) )))
        })
      }, simplify='array')
    }, simplify="array")
  })

  d <- lapply(1:S, function(s){
    out <- sapply(1:N_s[s], function(i){
      sapply(1:K, function(k){
        rmultinom(1, 1, prob=omega[[s]][k,i,])
      })
    }, simplify='array')
    aperm(out, c(2,3,1))
  })

  c_l <- c(1, rep(priors['c_h'], L-1))
  theta <- t(
    sapply(1:K, function(k){
      c(priors['epsilon'],  rinvgamma(L-1, priors['a_mu'], priors['b_mu']))
    })
  )

  H <- lapply(1:S, function(s) {
    sapply(1:N_s[s], function(i) sapply(1:K, function(k){
      rgamma(1, shape=2, rate=2/5)
    }))
  })

  pi_sa <- lapply(1:S, function(s) {
    array(rbeta(P, priors['a_a'], priors['b_a']), dim = c(P))
  })
  a_s <- lapply(1:S, function(s)
    matrix(rbern(P * N_s[s], pi_sa[[s]]), nrow = P, ncol = N_s[s]))

  Z_s <- lapply(1:S, function(s){
    out <- sapply(1:N_s[[s]], function(i){
      sapply(1:P, function(p){
        if (a_s[[s]][p, i] == 1) {
          array(rep(0, K), dim=c(K))
        } else {
          array(rmultinom(1, M_s[[s]][p, i], (W[p, ] * H[[s]][, i] +
                                                1e-10)), dim=c(K))
        }
      }, simplify="array")
    }, simplify="array")
    aperm(out, c(2,1,3))
  })

  ################################################################################
  # PERFORM SAMPLING
  ################################################################################
  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                         total = tot_iter,
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = FALSE,    # If TRUE, clears the bar when finish
                         width = 100)      # Width of the progress bar

  for (iter in 1:tot_iter) {

    if(iter == burnin){
      # Initialize chains
      out.length <- niter / thin
      W_chain <- array(dim = c(P, K, out.length))
      H_chain <- lapply(1:S, function(s) array(dim = c(K, N_s[s], out.length)))
      A_chain <- lapply(1:S, function(s) array(dim = c(P, N_s[s], out.length)))
      PI_chain <- lapply(1:S, function(s) array(dim = c(P, out.length)))
      D_chain = lapply(1:S, function(s) array(dim=c(K, N_s[s], L, out.length)))
      Theta_chain = array(dim=c(K, L, out.length))
      Beta_chain = lapply(1:S, function(s) array(dim=c(K, L-1, Q_s[s], out.length)))
    }


    pb$tick()

    # Drop duplicated patterns if they have a high cosine similarity
    if(iter < 0.9*burnin & iter > 0.25*burnin){
      k <- 1
      while(k<K) {
        dupes <- which(cosineDist(t(W[,k, drop=FALSE]), t(W[,-(1:k), drop=FALSE])) > 0.9)

        if(length(dupes)>0){
          dupes <- dupes + k

          W <- W[,-dupes,drop=FALSE]
          theta <- theta[-dupes,,drop=FALSE]
          for(s in 1:S){
            H[[s]][k,] <- H[[s]][k,] +  apply(H[[s]][dupes,,drop=FALSE], 2, sum)
            #Z_s[[s]][,k,] <-  Z_s[[s]][,k,] +  apply(Z_s[[s]][,dupes,, drop=FALSE], c(1,3), sum)

            H[[s]] <- H[[s]][-dupes,,drop=FALSE]
            Z_s[[s]] <- Z_s[[s]][,-dupes,,drop=FALSE]
            omega[[s]] <- omega[[s]][-dupes,,,drop=FALSE]
            d[[s]] <- d[[s]][-dupes,,,drop=FALSE]
            y[[s]] <- y[[s]][-dupes,,,drop=FALSE]
            beta[[s]] <- beta[[s]][-dupes,,,drop=FALSE]
          }
          K <- K - length(dupes)
        }
        k <- k + 1
      }

    }

    if(iter > 0.25*burnin & iter < 0.9*burnin){
      # Drop Patterns observed in less than 1% of total data
      drops <- which(apply(sapply(1:S, function(s) apply(d[[s]][,,2:L], c(1), sum) ), c(1), sum) < 0.01*sum(N_s))
      if(length(drops)>0){
        W <- W[,-drops,drop=FALSE]
        theta <- theta[-drops,,drop=FALSE]
        for(s in 1:S){
          H[[s]] <- H[[s]][-drops,,drop=FALSE]
          Z_s[[s]] <- Z_s[[s]][,-drops,,drop=FALSE]
          omega[[s]] <- omega[[s]][-drops,,,drop=FALSE]
          d[[s]] <- d[[s]][-drops,,,drop=FALSE]
          y[[s]] <- y[[s]][-drops,,,drop=FALSE]
          beta[[s]] <- beta[[s]][-drops,,,drop=FALSE]
        }
        K <- K - length(drops)
      }
    }


    # Sample Latent counts
    Z_s <- lapply(1:S, function(s){
      out <- array(0, dim=c(P, K, N_s[s]))
      nonzeros <- which(M_s[[s]]!=0, arr.ind = TRUE)

      for(idx in 1:nrow(nonzeros)){
        p <- nonzeros[idx, 1]
        i <- nonzeros[idx, 2]
        prob_vec <- (W[p, ] * H[[s]][, i]) + 1e-20
        out[p,,i] <- rmultinom(1, M_s[[s]][p,i], prob_vec)
      }

      out

    })

    # Sample pattern matrix
    W <- sapply(1:K, function(k){
      rgamma(P,
             priors['alpha_w'] + apply(sapply(1:S, function(s) apply(Z_s[[s]][, k, ], 1, sum)), 1, sum),
             priors['beta_w'] + (1- do.call(cbind, a_s)) %*% do.call(cbind, H)[k,]
             )
    })

    # Sample Exposures
    H <- lapply(1:S, function(s){
      L_s <- apply(d[[s]], c(1,2), function(x) which(x==1))


      matrix(
        rgamma(n = K*N_s[s],
               shape = c_l[L_s] +  apply(Z_s[[s]], c(2,3), sum),
               rate = c_l[L_s] / t(sapply(1:K, function(k) theta[k,L_s[k,]] ))  + t(W) %*% (1-a_s[[s]])
        ),
        nrow=K, ncol=N_s[s]
      )
    })


    # Sample a, dirac 0's
    a_s <- lapply(1:S, function(s) {
      matrix(rbern(P * N_s[s], (1 + exp(
        matrix(
          log(1 - pi_sa[[s]]) - log(pi_sa[[s]]),
          nrow = P,
          ncol = N_s[s]
        ) + dpois(M_s[[s]], lambda = W %*% H[[s]] + 1e-10, log = TRUE) - ifelse(M_s[[s]] ==
                                                                                  0, 0, -Inf)
      ))^(-1)), nrow = P, ncol = N_s[s])
    })

    # Sample probability of A
    pi_sa <- lapply(1:S, function(s) {
      sapply(1:P, function(p) {
        rbeta(1, priors['a_a'] + sum(a_s[[s]][p, ]), priors['b_a'] + N_s[s] - sum(a_s[[s]][p, ]))
      })
    })

    # Calculate cluster weights
    omega <-  lapply(1:S, function(s){
      mat_1 <- abind( apply(beta[[s]], c(1,2), function(beta, x) pnorm(x %*% beta), x=x_s[[s]] ), array(1, dim=c(N_s[s], K, 1)) )
      mat_2 <- abind(array(1, dim=c(N_s[s], K, 1)), 1-mat_1[,,-L])
      mat_2 <- aperm(apply(mat_2, c(1,2), cumprod), c(2,3,1))
      out_mat <- mat_1 * mat_2
      aperm(out_mat, c(2,1,3))
    })

    # Sample cluster membership
    d <- lapply(1:S, function(s){
      out <- sapply(1:N_s[s], function(i){
        sapply(1:K, function(k){
          rmultinom(n = 1, size = 1, prob= (omega[[s]][k,i,]) * dgamma(H[[s]][k,i], shape=c_l, rate=c_l/theta[k,]) + 1e-20)
        })
      }, simplify='array')
      aperm(out, c(2,3,1))
    })

    # Sample cluster atoms
    theta <- t(
      sapply(1:K, function(k){
        c(priors['epsilon'],
          sapply(2:L, function(l) {
            rinvgamma(1,
                      shape = priors['a_mu'] + c_l[l]*sum(sapply(1:S, function(s){sum(d[[s]][k,,l]==1)})),
                      scale = priors['b_mu'] +  c_l[l]*sum(sapply(1:S, function(s) sum(H[[s]][k, which(d[[s]][k,,l] == 1)]) ))
            )
          })
        )
      })
    )

    # Sample probits
    y <- lapply(1:S, function(s){
      l_mat <- 1*(d[[s]][,,-L] == 1)
      r_mat <- 1*aperm(apply(l_mat, c(1,2), function(x) (cumsum(x)==0)), c(2,3,1))
      mean_y <- aperm(apply(beta[[s]], c(1,2), function(beta, x) x %*% beta, x=x_s[[s]] ), c(2,1,3))

      lower_bounds <- c(-Inf, 0)
      upper_bounds <- c(Inf, 0)

      out <- array(NA_real_, dim=dim(mean_y))

      ids_to_sample <- l_mat | r_mat

      out[ids_to_sample] <- rtruncnorm(n=mean_y[ids_to_sample], a=lower_bounds[l_mat[ids_to_sample]+1], b=upper_bounds[r_mat[ids_to_sample]+1], mean=mean_y[ids_to_sample], sd=1)
      out
    })

    # Sample covariate effects
    beta <- lapply(1:S, function(s){
      out <- array(dim=c(Q_s[s], L-1, K))
      out[,,] <- sapply(1:K, function(k){
        sapply(1:(L-1), function(l){
          s_l <- apply(d[[s]][k,,], MARGIN=1, function(x) which(x==1)) >=l

          if(sum(s_l)>0){
            x_s_tilde = x_s[[s]][s_l,,drop=FALSE]
            y_s_tilde = matrix(y[[s]][k,s_l,l], ncol=1)
          } else {
            x_s_tilde = matrix(0, nrow=1, ncol=Q_s[s])
            y_s_tilde = matrix(0, nrow=1, ncol=1)
          }

          beta_covar = chol2inv(chol(priors['tau_0']*diag(Q_s[s]) + crossprod(x_s_tilde)))
          beta_covar = round(beta_covar, 5)
          beta_mu = c(beta_covar %*% (priors['tau_0']*diag(Q_s[s]) %*% matrix(c(priors['beta_0'],rep(0, Q_s[s]-1)),ncol=1) + t(x_s_tilde) %*% y_s_tilde))
          rmvn(1, mu=beta_mu, Sigma = beta_covar)
        })
      }, simplify="array")
      aperm(out, c(3, 2, 1))
    })

    if ((iter > burnin) & (iter %% thin == 0)) {
      idx <- (iter - burnin) / thin
      W_chain[, , idx] <- W
      Theta_chain[,,idx] <- theta
      for (s in 1:S) {
        H_chain[[s]][,,idx] <- H[[s]]
        A_chain[[s]][,,idx] <- a_s[[s]]
        PI_chain[[s]][, idx] <- pi_sa[[s]]
        D_chain[[s]][,,,idx] <- d[[s]]
        Beta_chain[[s]][,,,idx] <- beta[[s]]
      }
    }
  }

  return(list(
    "W" = W_chain,
    "H" = H_chain,
    "A" = A_chain,
    "PI" = PI_chain,
    "D" = D_chain,
    "Theta" = Theta_chain,
    "Beta" = Beta_chain
  ))
}
