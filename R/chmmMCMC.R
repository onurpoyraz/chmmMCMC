#' A Reference Class to represent a single chain in HMM or CHMM
#'
#' This class corresponds to chain of the HMM.
#' It performs all the MCMC sampling procedures for the sigle chain.
#' @field observations N x L matrix where rows corresponds to
#'     number of observations while columns corresponds to sequence lengths
#'     MUST be provided by user
#' @field params Parameter list that defines the model parameters.
#'     MUST be supplied by user.
#'     Details of that field are in the package description.
#' @field latent_sequences sampled_latent sequence. It is necessary for model to continue.
#'     It is not necessary for any further action.
#' @field log_likelihoods Log likelihoods of each sequence at each MCMC iteration.
#'     Rows corresponds to different sequences and columns corresponds to iterations
#' @field pi0s Sampled pi0s. Rows corresponds to number of states and columns corresponds iteration number.
#' @field betas Sampled beta parameters. Rows corresponds to beta values and columns corresponds iteration number.
#' @field transitions A list that model needs. Those transitions are converted from corresponding beta iteration.
#'     It is not necessary for any further action.
#' @field emissions Sampled emissions. It is 3 dimensional array and dim[1] corresponds to n_states,
#'     dim[2] corresponds to n_observations and dim[3] corresponds to number of MCMC samples.
#' @field covariates The latent sequences of the other chains. MUST be provided in the updating the chain.
#' @field forward_alpha Forward table. It will be not necessary for any further action.
#' @field diagnosis List that includes MCMC diagnosis. It is the result of the diagnostics() function.
#' @field posterior List that includes posterior probabilities of pi0, transition and emission.


Chain <- setRefClass(
  "Chain",
  fields = list(
    observations = "ANY",
    latent_sequences = "ANY",
    log_likelihoods = "ANY",
    pi0s = "ANY",
    betas = "ANY",
    transitions = "ANY",
    emissions = "ANY",
    params = "ANY",
    covariates = "ANY",
    forward_alpha = "ANY",
    diagnosis = "ANY",
    posterior = "ANY"
  ),
  methods = list(
    initialize = function(observations, params)
    {
      #' Load observations
      .self$observations <<- observations
      # Load model parameters
      .self$params <<- params
      # Calculate step sizes and dimensions according to params and coupling_mode
      initialize_model()
      # Initialize the pi0 probabilities
      pi0s <<- array(data = NA, dim=c(.self$params$n_states,
                                      .self$params$n_iter))
      pi0s[,1] <<- .self$params$pi0_prior / sum(.self$params$pi0_prior)
      # Initialize the beta parameters
      betas <<- array(data = NA, dim=c(.self$params$dim_beta,
                                       .self$params$n_iter))
      betas[,1] <<- 0
      # Calculate the transition probabilities based on beta values
      transitions <<- calculate_transition(betas[, 1]) ## some functions
      # Initialize the emisssion probabilities
      emissions <<- array(data = NA, dim=c(.self$params$n_states,
                                           .self$params$n_observation,
                                           .self$params$n_iter))
      emissions[,,1] <<- .self$params$emission_prior / rowSums(.self$params$emission_prior)
      # Initialize the covariates
      covariates <<- array(1, dim = c(.self$params$N,
                                      .self$params$L,
                                      .self$params$n_covariates))
      # Initialize the log-likelihood
      log_likelihoods <<- array(data = NA, dim=c(.self$params$N,
                                                 .self$params$n_iter))
      # Initialize forward_alpha
      forward <- forward_filtering(transitions, 1)
      forward_alpha <<- forward$alpha
      log_likelihoods[,1] <<- forward$log_prob
      # Initialize latent sequences
      latent_sequences <<- array(NA, dim = c(.self$params$N,
                                             .self$params$L))
      backward_sampling()
      # Set diagnostics
      diagnosis <<- list()
      diagnosis$acceptance_count_warmup <<- 0
      diagnosis$acceptance_count <<- 0
      # Set posterior list
      posterior <<- list()
    },
    initialize_model = function()
    {
      # Set dimension of the beta and step sizes
      intercept_params <- rep(params$step_size_intercept, params$n_states * (params$n_states - 1))
      if (params$model == "Single"){
        step_size <- intercept_params
        params$n_covariates <<- 0
      } else if (params$model == "Or"){
        covariate_params <- rep(params$step_size_coefficient, params$n_states * (params$n_states - 1))
        covariates_params <- rep(covariate_params, params$n_states - 1)
        step_size <- c(intercept_params, covariates_params)
        calculate_grid()
      } else if (params$model == "Additive"){
        covariate_params <- rep(params$step_size_coefficient, params$n_states * (params$n_states - 1)**2)
        covariates_params <- rep(covariate_params, params$n_covariate)
        step_size <- c(intercept_params, covariates_params)
        calculate_grid()
      } else if (params$model == "Full"){
        covariate_params <- rep(params$step_size_coefficient, params$n_states * (params$n_states - 1))
        covariates_params <- rep(covariate_params, params$n_states ** params$n_covariates - 1)
        step_size <- c(intercept_params, covariates_params)
        calculate_grid()
      }
      params$step_size <<- step_size
      params$dim_beta <<- length(step_size)
      params$N <<- dim(observations)[1] ## number of sequence
      params$L <<- dim(observations)[2] ## length of sequence
    },
    calculate_grid = function()
    {
      l <- list()
      for (i in 1:params$n_covariates) {
        l[[i]] <- 1:params$n_states
      }
      ll <- expand.grid(l)
      params$grid <<- array(as.numeric(unlist(ll)), dim=dim(ll))
    },
    update_chain = function(new_covariates, iter)
    {
      covariates <<- new_covariates
      sample_pi0(iter)
      sample_emission(iter)
      sample_beta(iter)
      backward_sampling()
    },
    sample_pi0 = function(iter)
    {
      # Count the number of times each state appears as the first.
      # Draw a sample from Dirichlet distribution with given counts and given prior
      # Returns p0_sample (n_states)
      counts <- tabulate(latent_sequences[,1], nbins=params$n_states)
      dirichlet_params <- counts + params$pi0_prior
      pi0s[, iter] <<- c(LaplacesDemon::rdirichlet(1, dirichlet_params))
    },
    sample_emission = function(iter)
    {
      # Count the number of emissions from all hidden states to all data states.
      # Draw a sample from Dirichlet distribution with given counts and given prior
      # Returns emission_sample (n_states, n_observation)
      counts <- array(0, dim = c(params$n_states, params$n_observation))
      for(n in 1:params$N) {
        for(l in 1:params$L) {
          i <- latent_sequences[n, l]
          j <- observations[n, l]
          counts[i, j] <- counts[i, j] + 1
        }
      }
      dirichlet_params <- counts + params$emission_prior
      emissions[, , iter] <<- LaplacesDemon::rdirichlet(params$n_states, dirichlet_params)
    },
    forward_filtering = function(transition_list, iter)
    {
      emission_tensor <- tensorA::to.tensor(emissions[, , iter])
      alpha <- tensorA::to.tensor(NA, c(S=params$n_states, L=params$L, N=params$N))
      log_prob <- tensorA::to.tensor(0, dim=c(N=params$N))

      if (params$model == "Single") {
        transition_tensor <- tensorA::to.tensor(array(transition_list[[1]], dim = c(params$n_states,
                                                                                    params$n_states,
                                                                                    params$N)))
        names(transition_tensor) <- c("S", "C", "N")
      }

      for(t in 1:params$L) {
        if (t == 1) {
          alpha[[L=1]] <- pi0s[, iter]
        }
        else {
          if (params$model != "Single") {
            covariates_t <- covariates[,t-1,,drop=FALSE]
            dim(covariates_t) <- c(dim(covariates_t)[1], dim(covariates_t)[3])
            selector <- apply(covariates_t, 1, .self$enumerator)
            transition <- transition_list[selector]
            transition_tensor <- tensorA::to.tensor(array(as.numeric(unlist(transition)),
                                                          dim = c(params$n_states,
                                                                  params$n_states,
                                                                  params$N)))
            names(transition_tensor) <- c("S", "C", "N")
          }
          alpha[[L=t]] <- tensorA::einstein.tensor(alpha[[L=t-1]], transition_tensor, by="N")
        }
        not_na <- !is.na(observations[,t])
        if(!all(not_na==FALSE)) {
          alpha[,t,not_na] <- emission_tensor[,observations[,t][not_na]] * alpha[,t,not_na]
        }
        prob <- tensorA::mean.tensor(alpha[[L=t]], along="S") * params$n_states
        alpha[[L=t]] <- alpha[[L=t]] / prob
        log_prob <- log_prob + log(prob)
      }
      return(list(alpha=alpha, log_prob=log_prob))
    },
    sample_seqs = function(weights, y)
    {
      x <- runif(nrow(weights))
      cumul.w <- weights %*% upper.tri(diag(ncol(weights)), diag = TRUE) /
        rowSums(weights)
      i <- rowSums(x > cumul.w) + 1L
      return (y[i])
    },
    backward_sampling = function()
    {
      for(t in (params$L):1) {
        if (t == params$L) {
          sampling_distribution <- forward_alpha[[L=t]]
        } else {
          if (params$model == "Single") {
            transition_tensor <- tensorA::to.tensor(transitions[[1]][,latent_sequences[,t+1]])
          } else {
            covariates_t <- covariates[,t,,drop=FALSE]
            dim(covariates_t) <- c(dim(covariates_t)[1], dim(covariates_t)[3])
            selector <- apply(covariates_t, 1, .self$enumerator)
            transition <- transitions[selector]
            transition_array <- mapply(function(x, y) x[,y], transition, latent_sequences[,t+1])
            transition_tensor <- tensorA::to.tensor(transition_array)
          }
          names(transition_tensor) <- c("S", "N")
          density <- transition_tensor * forward_alpha[[L=t]]
          sampling_distribution <- density / (tensorA::mean.tensor(density, along="S") * params$n_states)
        }
        latent_sequences[,t] <<- sample_seqs(t(sampling_distribution), 1:params$n_states)
      }
    },
    softmax = function (x)
    {
      return(exp(x - matrixStats::rowLogSumExps(x)))
    },
    calculate_transition = function(beta_values)
    {
      transition_list <- list()
      dim(beta_values) <- c(params$n_states,
                            params$n_states - 1,
                            length(beta_values) / (params$n_states * (params$n_states - 1)))
      beta <- abind::abind(beta_values, -apply(beta_values, c(1,3), FUN=sum), along=2)
      intercept <- beta[,,1]
      transition_list[[1]] <- softmax(intercept)
      if (params$model == "Or"){
        for (n in 2:dim(params$grid)[1]) {
          states <- params$grid[n,]
          combine <- intercept
          for (i in 2:params$n_states) {
            if (any(states == i)) {
              combine <- combine + beta[,,i]
            }
          }
          transition_list[[n]] <- softmax(combine)
        }
      } else if (params$model == "Additive"){
        for (n in 2:dim(params$grid)[1]) {
          states <- params$grid[n,]
          combine <- intercept
          j <- 0
          for (i in states) {
            k <- 1 + j * (params$n_states - 1)
            if (i > 1) {
              combine <- combine + beta[,,i + k - 1]
            }
            j <- j + 1
          }
          transition_list[[n]] <- softmax(combine)
        }
      } else if (params$model == "Full"){
        for (n in 2:dim(params$grid)[1]) {
          transition_list[[n]] <- softmax(intercept + beta[,,n])
        }
      }
      return(transition_list)
    },
    enumerator = function(covariate){
      y <- 1
      for (i in 1:params$n_covariates) {
        y <- y + (covariate[i] - 1) * params$n_states**(i-1)
      }
      return(y)
    },
    sample_beta = function(iter)
    {
      num_params <- length(betas[,iter - 1])
      if (iter > params$warmup) {
        c_square <- 5.76 / num_params
        covariance <- cov(t(betas[,round(params$warmup/2):(iter - 1)])) * c_square
      } else {
        c_square <- params$step_size/num_params
        covariance <- diag(length(betas[,iter - 1])) * c_square
      }
      beta_values <- c(mvtnorm::rmvnorm(1, betas[,iter - 1], covariance, method = "chol"))
      transition_list <- calculate_transition(beta_values)

      metropolis(transition_list, beta_values, iter)
    },
    metropolis = function(transition_list, beta_values, iter)
    {
      proposal <- calculate_log_posterior(transition_list, beta_values, iter)
      previous <- calculate_log_posterior(transitions, betas[, iter-1], iter)
      ratio <- exp(proposal$log_posterior - previous$log_posterior)
      acceptance_rate <- runif(1)
      if (ratio > acceptance_rate){
        transitions <<- transition_list
        betas[, iter] <<- beta_values
        forward_alpha <<- proposal$alpha
        log_likelihoods[,iter] <<- proposal$log_prob
        if (iter > params$warmup) {
          diagnosis$acceptance_count <<- diagnosis$acceptance_count + 1
        } else {
          diagnosis$acceptance_count_warmup <<- diagnosis$acceptance_count_warmup + 1
        }
      } else {
        betas[, iter] <<- betas[, iter-1]
        forward_alpha <<- previous$alpha
        log_likelihoods[,iter] <<- previous$log_prob
      }
    },
    calculate_log_posterior = function(transition_list, beta_values, iter)
    {
      forward <- forward_filtering(transition_list, iter)
      log_likelihood <- sum(forward$log_prob)
      log_prior <- calculate_log_prior(beta_values)
      forward$log_posterior <- log_likelihood + log_prior
      return(forward)
    },
    calculate_log_prior = function(beta_values){
      intercept_params <- 1:(params$n_states * (params$n_states - 1))
      log_prior_intercept <- sum(dnorm(beta_values[intercept_params], 0, params$sd_intercept, log = TRUE))
      log_prior_coefficient <- prior_density(beta_values[-intercept_params])
      log_prior <- log_prior_intercept + log_prior_coefficient
      return(log_prior)
    },
    prior_density = function(values)
    {
      if (params$n_covariates == 0) {
        y <- 0
      } else if (params$prior_distribution == "Gaussian") {
        y <- sum(dnorm(values, 0, params$sd_coefficient, log = TRUE))
      } else if (params$prior_distribution == "Laplace"){
        y <- sum(LaplacesDemon::dlaplace(values, 0, params$sd_coefficient, log = TRUE))
      } else if (params$prior_distribution == "Horseshoe"){
        y <- sum(dhorseshoe(values, 0, params$sd_coefficient, log = TRUE))
      }
      return(y)
    },
    dhorseshoe = function(x, mu, scale, log=TRUE) {
      xx <- ((x - mu) / scale)**2 / 2
      g <- 0.5614594835668851  # exp(-0.5772156649015328606)
      b <- 1.0420764938351215   # sqrt(2 * (1-g) / (g * (2-g)))
      h_inf <- 1.0801359952503342  #  (1-g)*(g*g-6*g+12) / (3*g * (2-g)**2 * b)
      q <- 20. / 47. * xx**1.0919284281983377
      h <- 1. / (1 + xx**(1.5)) + h_inf * q / (1 + q)
      c <- 0.5 * log(2 * pi**3) + log(g * scale)
      log_density <- log(log(1 + g / (xx + 1e-12) - (1 - g) / (h + b * xx)**2)) - log(1 + (1 - g) / g * exp(-xx / (1 - g))) - c
      if (log) {
        output <- log_density
      } else {
        output <- exp(log_density)
      }
      return(output)
    },
    diagnostics = function()
    {
      diagnosis$waic_score <<- loo::waic(t(log_likelihoods[,(params$warmup + 1):params$n_iter]))
      diagnosis$loo_score <<- loo::loo(t(log_likelihoods[,(params$warmup + 1):params$n_iter]))
    },
    full_transition = function(select=(params$warmup+1):params$n_iter)
    {
      return (apply(betas[,select], 2, function(x) abind::abind(calculate_transition(x), rev.along=0)))
    },
    calculate_posteriors = function()
    {
      if (params$n_iter > params$warmup) {
        select <- (params$warmup+1):params$n_iter
      } else {
        select <- 1:params$n_iter
      }
      posterior$pi0 <<- rowMeans(pi0s[,select])
      posterior$emission <<- apply(emissions[,,select], c(1,2), FUN=mean)
      t_array <- rowMeans(full_transition(select))
      dim(t_array) <- c(params$n_states, params$n_states, length(transitions))
      posterior$transitions <<- lapply(1:length(transitions), function(i) t_array[,,i])
      posterior$beta <<- rowMeans(betas[,select])
    },
    clear_garbage = function()
    {
      latent_sequences <<- NULL
      transitions <<- NULL
      covariates <<- NULL
      forward_alpha <<- NULL
    }
  )
)
#' A Reference Class to CHMM
#'
#' This class corresponds to CHMM training.
#' It basically manage the Chain class.
#' @field params Parameter list that defines the model parameters.
#'     MUST be supplied by user.
#'     Details of that field are in the package description.
#' @field chain_names Optional parameter to give the chains a name.
#' @field chains A list that store Chain class objects. It has the size of  number of chains.
#' @field log_likelihoods Log likelihoods of each sequence at each MCMC iteration.
#'     It is the sum of log_likelihoods of all chains.
#'     Rows corresponds to different sequences and columns corresponds to iterations
#' @field n_chains Number of chains. Class calculate it from length of observation list.
#' @field runtime Total time spend on training
#' @field diagnosis Model diagnostics and model scores.


Chmm_Mcmc <- setRefClass(
  "Chmm_Mcmc",
  fields = list(
    params = "ANY",
    chain_names = "ANY",
    chains = "ANY",
    log_likelihoods = "ANY",
    n_chains = "ANY",
    runtime = "ANY",
    diagnosis = "ANY"
  ),
  contains = "Chain",
  methods = list(
    initialize = function(observations, params, chain_names=NULL)
    {
      .self$params <<- params
      .self$n_chains <<- length(observations)
      if (is.null(chain_names)) {
        .self$chain_names <<- sapply(1:n_chains, function(x) paste("C", toString(x), sep=""))
      } else {
        .self$chain_names <<- chain_names
      }
      .self$chains <<- list()
      for (chain in 1:n_chains) {
        .self$chains[[chain]] <<- Chain$new(observations=observations[[chain]], params=params)
      }
      .self$diagnosis <<- list()
    },
    train = function()
    {
      ptm <- proc.time()
      for (iter in 2:params$n_iter) {
        for (chain in 1:n_chains) {
          new_covariates <- concat_covariates(chains[(-chain)])
          chains[[chain]]$update_chain(new_covariates, iter)
        }
      }
      log_likelihood_list <- lapply(chains, function(x) x$log_likelihoods)
      log_likelihoods <<- Reduce('+', log_likelihood_list)
      runtime <<- proc.time() - ptm
      for (chain in 1:n_chains) {
        chains[[chain]]$calculate_posteriors()
        chains[[chain]]$clear_garbage()
      }
      names(chains) <<- chain_names
    },
    concat_covariates = function(covariate_chains)
    {
      covariate_list <- lapply(covariate_chains, function(x) x$latent_sequences)
      new_covariates <- abind::abind(covariate_list, rev.along=0)
      return(new_covariates)
    },
    diagnostics = function()
    {
      diagnosis$waic_score <<- loo::waic(t(log_likelihoods[,(params$warmup + 1):params$n_iter]))
      diagnosis$loo_score <<- loo::loo(t(log_likelihoods[,(params$warmup + 1):params$n_iter]))
    }
  )
)
