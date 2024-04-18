# genAB generates two vectors, alpha and beta, based on the given parameters.
# seed: A seed value for reproducibility of random generation.
# TrueM: The number of true values in alpha and beta generated from a normal distribution.
# p1, p2, p3: Parameters to control the length and properties of segments in alpha and beta.
genAB <- function(seed, TrueM = 50, p1 = 150, p2 = 100, p3 = 300) {
  # Set the seed for reproducibility
  set.seed(seed) 
  
  # Generate the 'alpha' vector:
  # - The first TrueM elements are drawn from a normal distribution with mean 0 and standard deviation 1.5.
  # - The next p1 elements are zeros.
  # - The following p2 elements are drawn from a normal distribution (mean 0, sd 1.5).
  # - The last p3 elements are zeros.
  alpha <- c(rnorm(TrueM, mean = 0, sd = 1.5),  # True values
             rep(0, p1),                       # Zero padding
             rnorm(p2, mean = 0, sd = 1.5),    # More true values
             rep(0, p3))                       # Zero padding
  
  # Generate the 'beta' vector:
  # - The first TrueM elements are drawn from a normal distribution with mean 0 and standard deviation 1.5.
  # - The next p1 elements are also drawn from a normal distribution (same parameters).
  # - The remaining p2 + p3 elements are zeros.
  beta <- c(rnorm(TrueM, mean = 0, sd = 1.5),  # True values
            rnorm(p1, mean = 0, sd = 1.5),     # More true values
            rep(0, p2 + p3))                   # Zero padding
  
  # Return a list containing both alpha and beta vectors
  return(list(alpha = alpha, beta = beta))
}