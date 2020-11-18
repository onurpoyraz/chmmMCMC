\name{chmmMCMC}
\alias{chmmMCMC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Coupled Hidden Markov Models with MCMC Sampling
}
\description{
Calculate the CHMM and HMM model parameters using backward sampling. Takes list of observations where each observation corresponds to a chain and model parameters list as an argument.
}
\usage{
model <- Chmm_Mcmc(observations, params)
model$train()
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{observations}{
  List of observations. Observations should be in matrix form where columns corresponds to observation length and rows corresponds to different sequence.
  }
  \item{params}{
  List of model parameters. Params MUST have following items:
  }
  \item{n_states}{
  Number of latent states in the model
  }
  \item{n_observation}{
  Number of observations states in the model.
  }
  \item{model}{
  Defines the model. It could be one of the following: "Single", "Or", "Additive" or "Full"
  }
  \item{n_covariates}{
  Number of covariate chains. It should be 0 if the model is "Single".
  }
  \item{prior_distribution}{
  Choses the prior distribution. Should be one of the following: "Gaussian", "Laplace", or "Horseshoe". Horseshoe provides sparse parameters.
  }
  \item{step_size_intercept}{
  Step size of intercept term. Should be relatively small number. 0.1 is recommended.
  }
  \item{step_size_coefficients}{
  Step size of coefficient term. Should be relatively small number. 0.1 is recommended.
  }
  \item{sd_intercept}{
  Prior standard deviation of the intercept term. Should be high for weak prior ad small for strong prior.
  }
  \item{sd_coefficient}{
  Prior standard deviation of the coefficient term. Should be high for weak prior ad small for strong prior.
  }
  \item{n_iter}{
  Number of iteration. At least 2000 iteration is recommended.
  }
  \item{warmup}{
  Warmup length. At least 2000 iteration is recommended.
  }
  \item{pi0_prior}{
  Initial state prior. Must have the length of n_states. It is the prior for dirichlet.
  }
  \item{emission_prior}{
  Emission prior. Rows represents 'n_states' while columns represents 'n_observations'. Must be the matrix of N_states x N_observations. Each row of that matrix is prior for dirichlet distribution.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Return an Chmm-Mcmc object. It includes the overall model evaluation of the whole chains in diagnosis. Additionallyy Chmm-Mcmc object has "chains" field which includes all the necessary information about each chain and their parameters.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Onur Poyraz
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Dummy observations from different chains
observation1=matrix(c(1,2,3,3,NA,3,2,1,
                      1,2,3,3,2,3,2,1,
                      1,2,3,3,3,3,2,1), byrow=TRUE, nrow=3)

observation2=matrix(c(1,2,3,3,NA,3,2,1,
                      3,2,3,1,3,1,2,1,
                      3,2,3,2,NA,2,2,1), byrow=TRUE, nrow=3)

## Create observations list
observations <- list()
observations[[1]] <- observation1
observations[[2]] <- observation2

## Define the model parameters
n_states <- 3
n_observation <- 3
certainty <- 0.9
pi0_prior <- rep(1, n_states)
emission_prior <- array((1-certainty) / (n_observation-1), dim = c(n_states, n_observation))
diag(emission_prior) <- certainty
params <- list(n_states = n_states,
               n_observation = n_observation,
               n_covariates = 1, ## number of the remaining chain
               model="Additive",
               prior_distribution = "Laplace",
               step_size_intercept = 0.01,
               step_size_coefficient = 0.01,
               sd_intercept = 5,
               sd_coefficient = 5,
               n_iter = 10,
               warmup = 5,
               pi0_prior = pi0_prior,
               emission_prior = emission_prior
)

# Initialize the model
model <- Chmm_Mcmc(observations, params)

# Train the model
model$train()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~coupled hidden markov model }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~markov chain monte carlo }% __ONLY ONE__ keyword per line
