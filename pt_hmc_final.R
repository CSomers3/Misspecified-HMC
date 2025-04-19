# Bayesian Misspecification Analysis via (Parallel Tempered) HMC

library(coda)
library(mvtnorm)
library(parallel)
library(pryr)
library(gridExtra)
library(ggplot2)
library(corrplot)


# ------------------ DATA GENERATION ------------------
set.seed(3499)
n <- 100
p <- 10

# Create design matrix with controlled correlation structure
Z <- matrix(nrow = n, ncol = p)
Z[,1] <- rep(1, n) # intercept column
noise <- matrix(rnorm(n*(p-1)), nrow = n, ncol = p-1)
Z[,2] <- noise[,1]
Z[,3] <- noise[,1] + 0.2*noise[,2]
Z[,4] <- 0.5*noise[,3]
Z[,5] <- noise[,4]
Z[,6] <- 2*noise[,5] + 20*noise[,6]
Z[,7] <- noise[,6]
Z[,8] <- 0.5 * (noise[,7] + noise[,4] + noise[,8] + noise[,1])
Z[,9] <- noise[,8] + 10*noise[,4]
Z[,10] <- noise[,5] + 0.5*noise[,9]

# Create true parameter values
beta0 <- seq(-0.2, 0.2, length = p)
gamma0 <- 1/4

# Create predictor with bimodal distribution
x <- rep(NA, n)
x[1:(n/2)] <- runif(n/2, min = 0.5, max = 1.5)
x[(n/2+1):n] <- runif(n/2, min = 4, max = 8)


# ------------------------------------------------------
# Skew-Normal Distribution Functions
# ------------------------------------------------------
# Generate samples from Skew-Normal distribution SN(m, sigma, alpha)
rskew_normal <- function(n, m = 0, sigma = 1, alpha = 0) {
  Z <- rnorm(n, mean = m, sd = sigma)
  probs <- pnorm(alpha * (Z - m) / sigma)
  b <- 2 * (runif(n) <= probs) - 1
  return(b * Z)
}

# Density function for Skew-Normal
dskew_normal <- function(x, m = 0, sigma = 1, alpha = 3) {
  2 * dnorm(x, mean = m, sd = sigma) * pnorm(alpha * (x - m) / sigma)
}

# Log-density function for Skew-Normal
ldskew_normal <- function(x, m = 0, sigma = 1, alpha = 3) {
  log(2) + dnorm(x, mean = m, sd = sigma, log = TRUE) + 
    pnorm(alpha * (x - m) / sigma, log.p = TRUE)
}

# Generate error terms from SN(0, 1, 3)
epsilon <- rskew_normal(n, m = 0, sigma = 1, alpha = 3)

# Generate response variable from quadratic model
y <- gamma0 * x^2 + Z %*% beta0 + epsilon


# Predictor correlation mat
corrplot(cor(Z[, -1]), method = 'number', type = "lower", tl.cex = 0.8, tl.col = "black")

# Plot the skew-normal distribution
plot_skew_normal <- function() {
  ggplot(data.frame(epsilon), aes(x = epsilon)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50,
                   fill = "lightblue", color = "white", linewidth = 0.5) +
    stat_function(fun = function(x) 2 * dnorm(x) * pnorm(3 * x),
                  color = "red", linewidth = 1) +
    labs(title = "Skew-Normal Error Distribution",
         x = "Value", y = "Density") +
    theme_classic(base_size = 22)
}

# Plot the relationship between x and y
plot_xy_relationship <- function() {
  ggplot(data.frame(x, y), aes(x, y)) +
    geom_point(color = "blue", alpha = 0.6, size = 1.5) +
    stat_function(fun = function(x) gamma0 * x^2 + beta0[1], color = "red", linewidth = 1) +
    labs(title = "Data with True Quadratic Relationship", x = "x", y = "y") +
    theme_classic(base_size = 22)
}




# ------------------ POSTERIOR FUNCTIONS ------------------
log_prior <- function(theta) {
  mvtnorm::dmvnorm(theta, log = TRUE)
}

log_posterior <- function(theta, model = "linear") {
  gamma <- theta[1]
  beta <- theta[-1]
  
  mu <- if(model == "linear") gamma * x + Z %*% beta else gamma * x^2 + Z %*% beta
  
  residuals <- y - mu
  log_lik <- sum(ldskew_normal(residuals))
  
  log_lik + log_prior(theta)
}

grad_log_posterior <- function(theta, model = "linear") {
  gamma <- theta[1]
  beta <- theta[-1]
  
  if(model == "linear") {
    mu <- gamma * x + Z %*% beta
    x_term <- x
  } else {
    mu <- gamma * x^2 + Z %*% beta
    x_term <- x^2
  }
  
  residuals <- y - mu
  
  alpha <- 3
  alpha_res <- alpha * residuals
  ratio_term <- dnorm(alpha_res) / pmax(pnorm(alpha_res), 1e-10)
  
  score <- residuals - alpha * ratio_term
  
  grad_gamma <- sum(score * x_term)
  grad_beta <- as.vector(t(Z) %*% score)
  
  c(grad_gamma, grad_beta) - theta  # Prior gradient (assuming standard normal)
}



# ------------------ PT-HMC IMPLEMENTATION------------------
# Implements parallel tempering HMC + automatic step size adaptation
# - Multi-chain temperature swapping for rugged posteriors
# - Chunked warmup with step size adaption
# - Skew-normal likelihood integration


run_hmc <- function(model = "linear", n_samples = 50000, warmup = 5000, 
                    L = 15, epsilon = 0.0015, n_temps = 5, temp_ladder = NULL, 
                    swap_every = 10, target_accept = 0.8, n_cores = detectCores() - 1,
                    x_input = NULL, Z_input = NULL, y_input = NULL) {
  start_time <- Sys.time()
  mem_before <- mem_used()
  
  p_dim <- p + 1
  x_local <- if(!is.null(x_input)) x_input else x
  Z_local <- if(!is.null(Z_input)) Z_input else Z
  y_local <- if(!is.null(y_input)) y_input else y
  
  
  # Regular HMC if n_temps=1, PT-HMC otherwise
  is_regular_hmc <- (n_temps == 1)
  
  if(is_regular_hmc) {
    temp_ladder <- 1
    message("Running HMC...")
  } else {
    temp_ladder <- if(is.null(temp_ladder)) exp(seq(0, log(5), length.out = n_temps)) else temp_ladder
    n_temps <- length(temp_ladder)
    message("Running parallel tempering HMC with ", n_temps, " temperatures...")
  }
  
  # Set up parallel processing (if needed)
  if(n_temps > 1 && n_cores > 1) {
    cl <- makeCluster(min(n_temps, n_cores))
    on.exit(stopCluster(cl))
    
    clusterExport(cl, c("x", "Z", "y", "p", "log_posterior", "grad_log_posterior", 
                        "log_prior", "ldskew_normal", "mem_used", "effectiveSize"), 
                  envir = environment())
  } else {
    cl <- NULL
  }
  
  tempered_log_posterior <- function(theta, model, temperature) {
    if(temperature == 1) return(log_posterior(theta, model))
    
    gamma <- theta[1]
    beta <- theta[-1]
    mu <- if(model == "linear") gamma * x + Z %*% beta else gamma * x^2 + Z %*% beta
    
    residuals <- y - mu
    log_lik <- sum(ldskew_normal(residuals))
    
    log_lik / temperature + log_prior(theta)
  }
  
  tempered_grad_log_posterior <- function(theta, model, temperature) {
    grad_log_posterior(theta, model) / temperature
  }
  
  chain_states <- lapply(1:n_temps, function(i) rnorm(p_dim, sd = 0.1 * i))
  samples <- if(is_regular_hmc) matrix(0, nrow = n_samples, ncol = p_dim) else matrix(0, nrow = n_samples, ncol = p_dim)
  
  accepts <- rep(0, n_temps)
  swap_attempts <- swap_accepts <- matrix(0, nrow = n_temps, ncol = n_temps)
  
  # ---- Hamiltonian Dynamics Core ----
  run_hmc_steps <- function(i_temp, n_steps) {
    #  Leapfrog integrator with Metropolis acceptance
    temp <- temp_ladder[i_temp]
    theta <- chain_states[[i_temp]]
    local_accepts <- 0
    trajectory <- matrix(NA, nrow = n_steps, ncol = p_dim)
    
    log_post <- function(t) tempered_log_posterior(t, model, temp)
    grad_post <- function(t) tempered_grad_log_posterior(t, model, temp)
    
    for(i in 1:n_steps) {
      r <- rnorm(p_dim)                                 # Momentum: r ~ N(0, I)
      current_H <- -log_post(theta) + 0.5 * sum(r^2)    # Hamiltonian: -logpost + K(p)
      theta_prop <- theta
      r_prop <- r
      
      # Leapfrog integration
      for(j in 1:L) {
        r_prop <- r_prop + 0.5 * epsilon * grad_post(theta_prop)  # 0.5 step momentum
        theta_prop <- theta_prop + epsilon * r_prop               # full step position
        r_prop <- r_prop + 0.5 * epsilon * grad_post(theta_prop)  # 0.5 step momentum
      }
      
      proposed_H <- -log_post(theta_prop) + 0.5 * sum(r_prop^2)
      # MH acceptance probability
      log_alpha <- min(0, current_H - proposed_H) 
      
      if(is.finite(log_alpha) && log(runif(1)) < log_alpha) {
        theta <- theta_prop
        local_accepts <- local_accepts + 1
      }
      
      trajectory[i, ] <- theta
    }
    
    list(theta = theta, accepts = local_accepts, trajectory = trajectory)
  }

  # Warmup run
  warmup_chunks <- 10
  chunk_size <- ceiling(warmup / warmup_chunks)

  message("Starting (chunked) warmup phase...")
  for(chunk in 1:warmup_chunks) {
    chunk_accepts <- rep(0, n_temps)  # Reset chunk acceptances
    
    for(step in 1:ceiling(chunk_size / swap_every)) {
      # Run HMC steps
      if(is_regular_hmc || is.null(cl)) {
        results <- lapply(1:n_temps, function(i_temp) run_hmc_steps(i_temp, swap_every))
      } else {
        results <- parLapply(cl, 1:n_temps, function(i_temp) run_hmc_steps(i_temp, swap_every))
      }
      
      for(i_temp in 1:n_temps) {
        chain_states[[i_temp]] <- results[[i_temp]]$theta
        accepts[i_temp] <- accepts[i_temp] + results[[i_temp]]$accepts
        chunk_accepts[i_temp] <- chunk_accepts[i_temp] + results[[i_temp]]$accepts
      }
      
      # Temperature swaps (only for PT-HMC)
      if(!is_regular_hmc) {
        # Adjacent temperature swaps using MH ratio
        # only swap neighbours to maintain reasonable acceptance 
        for(i_temp in 1:(n_temps-1)) {
          j_temp <- i_temp + 1
          theta_i <- chain_states[[i_temp]]
          theta_j <- chain_states[[j_temp]]
          
          # cross temperature posterior checks
          log_post_i_at_i <- tempered_log_posterior(theta_i, model, temp_ladder[i_temp])
          log_post_i_at_j <- tempered_log_posterior(theta_i, model, temp_ladder[j_temp])
          log_post_j_at_i <- tempered_log_posterior(theta_j, model, temp_ladder[i_temp])
          log_post_j_at_j <- tempered_log_posterior(theta_j, model, temp_ladder[j_temp])
          
          # mh acceptance probability
          log_alpha <- (log_post_i_at_j + log_post_j_at_i) - (log_post_i_at_i + log_post_j_at_j)
          swap_attempts[i_temp, j_temp] <- swap_attempts[i_temp, j_temp] + 1
          
          if(is.finite(log_alpha) && log(runif(1)) < min(0, log_alpha)) {
            tmp <- chain_states[[i_temp]]
            chain_states[[i_temp]] <- chain_states[[j_temp]]
            chain_states[[j_temp]] <- tmp
            swap_accepts[i_temp, j_temp] <- swap_accepts[i_temp, j_temp] + 1
          }
        }
      }
    }
    
    # ---- Stepsize adaptation after each chunk ----
    chunk_accept_rate <- sum(chunk_accepts) / (chunk_size * n_temps)
    factor <- exp(2 * (chunk_accept_rate - target_accept))
    epsilon <- epsilon * factor
    
    message(sprintf("Chunk %d | acceptance rate: %.3f | adapted epsilon: %g",
                    chunk, chunk_accept_rate, epsilon))
  }
  
  
  # ---- Final warmup summary ----
  warmup_accept_rate <- sum(accepts) / (warmup * n_temps)
  message(sprintf("Final epsilon after warmup: %g (overall acceptance rate: %.3f)",
                  epsilon, warmup_accept_rate))
  
  # Reset metrics
  accepts <- rep(0, n_temps)
  swap_attempts <- swap_accepts <- matrix(0, nrow = n_temps, ncol = n_temps)
  
  
  # Main sampling phase
  message("Starting main sampling phase...")
  for(i in 1:ceiling(n_samples/swap_every)) {
    # Run HMC steps (parallel or sequential)
    if(is_regular_hmc || is.null(cl)) {
      results <- lapply(1:n_temps, function(i_temp) run_hmc_steps(i_temp, swap_every))
    } else {
      results <- parLapply(cl, 1:n_temps, function(i_temp) run_hmc_steps(i_temp, swap_every))
    }
    
    for(i_temp in 1:n_temps) {
      chain_states[[i_temp]] <- results[[i_temp]]$theta
      accepts[i_temp] <- accepts[i_temp] + results[[i_temp]]$accepts
      
      # For regular HMC, store all samples
      # For PT-HMC, store only coldest chain (T=1)
      if(is_regular_hmc || i_temp == 1) {
        start_idx <- (i-1)*swap_every + 1
        end_idx <- min(i*swap_every, n_samples)
        if(end_idx >= start_idx) {
          samples_to_store <- min(end_idx - start_idx + 1, swap_every)
          samples[start_idx:end_idx, ] <- results[[i_temp]]$trajectory[1:samples_to_store, ]
        }
      }
    }
    
    # Temperature swaps (only for PT-HMC)
    if(!is_regular_hmc) {
      for(i_temp in 1:(n_temps-1)) {
        j_temp <- i_temp + 1
        theta_i <- chain_states[[i_temp]]
        theta_j <- chain_states[[j_temp]]
        
        log_post_i_at_i <- tempered_log_posterior(theta_i, model, temp_ladder[i_temp])
        log_post_i_at_j <- tempered_log_posterior(theta_i, model, temp_ladder[j_temp])
        log_post_j_at_i <- tempered_log_posterior(theta_j, model, temp_ladder[i_temp])
        log_post_j_at_j <- tempered_log_posterior(theta_j, model, temp_ladder[j_temp])
        
        log_alpha <- (log_post_i_at_j + log_post_j_at_i) - (log_post_i_at_i + log_post_j_at_j)
        swap_attempts[i_temp, j_temp] <- swap_attempts[i_temp, j_temp] + 1
        
        if(is.finite(log_alpha) && log(runif(1)) < min(0, log_alpha)) {
          tmp <- chain_states[[i_temp]]
          chain_states[[i_temp]] <- chain_states[[j_temp]]
          chain_states[[j_temp]] <- tmp
          swap_accepts[i_temp, j_temp] <- swap_accepts[i_temp, j_temp] + 1
        }
      }
    }
  }
  
  # Performance metrics
  runtime <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  memory_mb <- (mem_used() - mem_before) / (1024*1024)
  
  # Calculate swap rates (only relevant for PT-HMC)
  swap_rates <- NULL
  if(!is_regular_hmc) {
    swap_rates <- matrix(0, nrow = n_temps, ncol = n_temps)
    for(i in 1:n_temps) {
      for(j in 1:n_temps) {
        if(swap_attempts[i,j] > 0) swap_rates[i,j] <- swap_accepts[i,j] / swap_attempts[i,j]
      }
    }
  }
  
  # results
  result <- list(
    samples = samples,
    acceptance = accepts/n_samples,
    epsilon = epsilon,
    runtime = runtime,
    memory_mb = memory_mb,
    algorithm = if(is_regular_hmc) "HMC" else "PT-HMC"
  )
  
  # + PT-specific information
  if(!is_regular_hmc) {
    result$swap_rates <- swap_rates
    result$temp_ladder <- temp_ladder
  }
  
  return(result)
}



# ------------------ BAYESBAG IMPLEMENTATION ------------------
bayesbag_hmc <- function(B = 50, n_cores = detectCores() - 1,
                         n_samples = 50000, warmup = 5000, 
                         L = 10, epsilon = 0.002) {
  start_time <- Sys.time()
  mem_before <- mem_used()
  
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl))
  
  # Export required functions and data
  clusterExport(cl, c("x", "Z", "y", "p", "log_posterior", 
                      "grad_log_posterior", "log_prior", "as.mcmc",
                      "ldskew_normal", "run_hmc", "mem_used", "effectiveSize",
                      "quantile", "parLapplyLB", "detectCores"), 
                envir = environment())
  
  clusterEvalQ(cl, { library(coda); library(mvtnorm); library(parallel) })
  
  # Worker function for a single bootstrap
  process_bootstrap <- function(b) {
    set.seed(3499 + b)
    n <- length(y)
    indices <- sample(n, replace = TRUE)
    x_b <- x[indices]
    Z_b <- Z[indices, , drop = FALSE]
    y_b <- y[indices]
    
    result <- run_hmc(
      model = "linear",
      n_samples = n_samples,
      warmup = warmup,
      L = L,
      epsilon = epsilon,
      x_input = x_b,
      Z_input = Z_b,
      y_input = y_b
    )
    
    list(
      index = b,
      samples = result$samples,
      acceptance_rate = result$acceptance,
      epsilon = result$epsilon,
      ess = effectiveSize(as.mcmc(result$samples))
    )
  }
  
  # Run all bootstraps in parallel
  results <- parLapplyLB(cl, 1:B, process_bootstrap)
  
  # Progress reporting (optional)
  cat(sprintf("Completed %d / %d bootstrap ensembles\n", B, B))
  flush.console()
  
  # Extract posterior means from each bootstrap sample
  posterior_means <- do.call(rbind, lapply(results, function(x) colMeans(x$samples)))
  
  # Calculate ESS values
  ess_values <- lapply(results, function(x) x$ess)
  min_ess <- sapply(ess_values, min)
  med_ess <- sapply(ess_values, median)
  
  credible_intervals <- apply(posterior_means, 2, 
                              function(param_samples) {
                                quantile(param_samples, probs = c(0.025, 0.975))
                              })
  
  memory_mb <- (mem_used() - mem_before) / (1024*1024)
  
  list(
    posterior_means = posterior_means,
    bagged_means = colMeans(posterior_means),
    bagged_vars = apply(posterior_means, 2, var),
    credible_intervals = credible_intervals,
    min_ess = min_ess,
    med_ess = med_ess,
    acceptance_rates = sapply(results, function(x) x$acceptance_rate),
    runtime = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
    memory_mb = memory_mb
  )
}




# ------------------ RWM IMPLEMENTATION ------------------
rwm <- function(model = "linear", n_samples = 50000, warmup = 5000, sigma_proposal = 0.01) {
  start_time <- Sys.time()
  mem_before <- mem_used()
  
  # get posterior
  log_post <- function(theta) log_posterior(theta, model)
  
  # initialise
  p_dim <- p + 1
  chain <- matrix(0, nrow = n_samples, ncol = p_dim)
  theta <- rep(0, p_dim)
  accepts <- 0
  
  # Warmup + proposal adaptation
  for (i in 1:warmup) {
    # Propose new value
    theta_prop <- theta + rnorm(p_dim, mean = 0, sd = sigma_proposal)
    
    # Metropolis acceptance step
    log_alpha <- log_post(theta_prop) - log_post(theta)
    
    if (is.finite(log_alpha) && log(runif(1)) < min(0, log_alpha)) {
      theta <- theta_prop
      accepts <- accepts + 1
    }
    
    # Adapt proposal variance
    if (i %% 100 == 0) {
      accept_rate <- accepts / i
      sigma_proposal <- sigma_proposal * ifelse(accept_rate > 0.3, 1.1, 
                                                ifelse(accept_rate < 0.2, 0.9, 1))
    }
  }
  
  # Main sampling
  accepts <- 0
  for (i in 1:n_samples) {
    # Propose new value
    theta_prop <- theta + rnorm(p_dim, mean = 0, sd = sigma_proposal)
    
    # Metropolis acceptance step
    log_alpha <- log_post(theta_prop) - log_post(theta)
    
    if (is.finite(log_alpha) && log(runif(1)) < min(0, log_alpha)) {
      theta <- theta_prop
      accepts <- accepts + 1
    }
    
    chain[i, ] <- theta
    
  }
  
  # performance metrics
  runtime <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  memory_mb <- (mem_used() - mem_before) / (1024*1024)
  
  list(
    samples = chain, 
    acceptance = accepts/n_samples, 
    sigma_proposal = sigma_proposal,
    runtime = runtime,
    memory_mb = memory_mb
  )
}


# ---------------- DATA GENERATION ---------------

plot_skew_normal()
plot_xy_relationship()


# ------------------ MODEL RUNS ------------------

nits <- 20000
warm <- 8000


# ---- MAIN MODELS ----

cat("Running Linear PT-HMC model...\n")
linear_fit <- run_hmc("linear", n_samples = nits, L = 15,
                      warmup = warm, n_temps = 5, epsilon = 0.0015)

cat("Running Quadratic PT-HMC model...\n")
quadratic_fit <- run_hmc("quadratic", n_samples = nits, L=15,
                         warmup = warm, n_temps = 5, epsilon = 0.0015)


cat("Running BayesBag HMC model...\n")
bayesbag_fit <- bayesbag_hmc(B = 50, n_samples = 20000, warmup = 8000)



# ---- STRAWMAN/BASELINE MISSPECIFIED MODELS ----

# RWM 
cat("Running Linear RWM model...\n")
linear_rwm_fit <- rwm("linear", n_samples = nits)

# HMC (NO Tempering)
cat("Running Linear PT-HMC model...\n")
hmc_line_fit <- run_hmc("linear", n_samples = nits,
                        warmup = warm, n_temps = 1)




# ------------------ DIAGNOSTICS ------------------

# main
linear_mcmc <- as.mcmc(linear_fit$samples)
quadratic_mcmc <- as.mcmc(quadratic_fit$samples)

linear_ess <- effectiveSize(linear_mcmc)
quadratic_ess <- effectiveSize(quadratic_mcmc)

# ESS Table
ess_table <- data.frame(
  Parameter = c("γ", paste0("β", 0:9)),
  Linear_HMC_ESS = round(linear_ess, 1),
  Quadratic_HMC_ESS = round(quadratic_ess, 1)
)

knitr::kable(ess_table, caption = "Effective Sample Size by Parameter", 
             col.names = c("Parameter", "Linear HMC ESS", "Quadratic HMC ESS"))

# comparisons
hmc_line_mcmc <- as.mcmc(hmc_line_fit$samples)
rwm_line_mcmc <- as.mcmc(linear_rwm_fit$samples)

rwm_ess <- effectiveSize(rwm_line_mcmc)
hmc_ess <- effectiveSize(hmc_line_mcmc)
geweke_linear <- geweke.diag(linear_mcmc)
geweke_quadratic <- geweke.diag(quadratic_mcmc)

# Extract posterior means and calculate fitted values
post_means_linear <- colMeans(linear_fit$samples)
post_means_quad <- colMeans(quadratic_fit$samples)
post_means_bayesbag <- bayesbag_fit$bagged_means

fitted_linear <- post_means_linear[1] * x + Z %*% post_means_linear[-1]
fitted_quadratic <- post_means_quad[1] * x^2 + Z %*% post_means_quad[-1]
fitted_bayesbag <- post_means_bayesbag[1] * x + Z %*% post_means_bayesbag[-1]

# Calculate MSE
mse_linear <- mean((y - fitted_linear)^2)
mse_quadratic <- mean((y - fitted_quadratic)^2)
mse_bayesbag <- mean((y - fitted_bayesbag)^2)

# 95% credible intervals
credible_interval_linear <- apply(linear_fit$samples, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
credible_interval_quadratic <- apply(quadratic_fit$samples, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
credible_interval_bayesbag <- bayesbag_fit$credible_intervals

# Helper function to format credible intervals
format_credible_intervals <- function(ci_matrix) {
  apply(ci_matrix, 2, function(ci) paste0("[",round(ci[1], 5), ", ", round(ci[2], 5), "]"))
}

parameter_table <- data.frame(
  Parameter = c("γ", paste0("β", 0:(length(post_means_linear) - 2))),
  True = c(gamma0, beta0),
  Linear_HMC = round(post_means_linear, 5),
  Linear_HMC_CI = format_credible_intervals(credible_interval_linear),
  Quadratic_HMC = round(post_means_quad, 5),
  Quadratic_HMC_CI = format_credible_intervals(credible_interval_quadratic),
  BayesBag = round(post_means_bayesbag, 5),
  BayesBag_CI = format_credible_intervals(credible_interval_bayesbag)
)

# MSE 
mse_row <- data.frame(
  Parameter = "MSE",
  True = NA,
  Linear_HMC = round(mse_linear, 5),
  Linear_HMC_CI = NA,
  Quadratic_HMC = round(mse_quadratic, 5),
  Quadratic_HMC_CI = NA,
  BayesBag = round(mse_bayesbag, 5),
  BayesBag_CI = NA
)

# Print with kable
parameter_table_with_mse <- rbind(parameter_table, mse_row)
knitr::kable(parameter_table_with_mse, 
             caption = "Posterior Means and Credible Intervals by Parameter (plus MSE)")


# performance metrics table
diagnostics_table <- data.frame(
  Metric = c("Acceptance Rate", "Min ESS", "Med ESS", "Runtime (sec)", "Memory (MB)"),
  Linear_PT_HMC = c(
    round(mean(linear_fit$acceptance), 3),
    round(min(linear_ess), 1),
    round(median(linear_ess), 1),
    round(linear_fit$runtime, 2),
    round(linear_fit$memory_mb, 2)
  ),
  Quadratic_PT_HMC = c(
    round(mean(quadratic_fit$acceptance), 3),
    round(min(quadratic_ess), 1),
    round(median(quadratic_ess), 1),
    round(quadratic_fit$runtime, 2),
    round(quadratic_fit$memory_mb, 2)
  ),
  Linear_RWM = c(
    round(linear_rwm_fit$acceptance, 3),
    round(min(rwm_ess), 1),
    round(median(rwm_ess), 1),
    round(linear_rwm_fit$runtime, 2),
    round(linear_rwm_fit$memory_mb, 2)
  ),
  Linear_HMC = c(
    round(hmc_line_fit$acceptance, 3),
    round(min(hmc_ess), 1),
    round(median(hmc_ess), 1),
    round(hmc_line_fit$runtime, 2),
    round(hmc_line_fit$memory_mb, 2)
  ),
  BayesBag = c(
    round(mean(bayesbag_fit$acceptance_rates), 3),
    round(mean(bayesbag_fit$min_ess), 1),
    round(mean(bayesbag_fit$med_ess), 1),
    round(bayesbag_fit$runtime, 2),
    round(bayesbag_fit$memory_mb, 2)
  )
)

# print
knitr::kable(diagnostics_table, caption = "Model Performance", align = c('l', 'r', 'r', 'r', 'r', 'r'))


# density plots
create_density_plots <- function(linear_samples, quadratic_samples, bayesbag_means) {
  param_names <- c("γ", paste0("β", 0:(ncol(linear_samples)-2)))
  true_values <- c(gamma0, beta0)
  
  plot_data <- rbind(
    data.frame(
      Model = "Linear",
      Parameter = rep(param_names, each = nrow(linear_samples)),
      Value = as.vector(linear_samples)
    ),
    data.frame(
      Model = "Quadratic",
      Parameter = rep(param_names, each = nrow(quadratic_samples)),
      Value = as.vector(quadratic_samples)
    )
  )
  
  # Create density plot
  density_plot <- ggplot(plot_data, aes(x = Value, fill = Model)) +
    geom_density(alpha = 0.5) +
    geom_vline(data = data.frame(Parameter = param_names, true_value = true_values),
               aes(xintercept = true_value), linetype = "dashed") +
    geom_vline(data = data.frame(Parameter = param_names, bayes_value = bayesbag_means),
               aes(xintercept = bayes_value), linetype = "dotted", color = "purple") +
    facet_wrap(~ Parameter, scales = "free", ncol = 4) +
    theme_minimal() +
    labs(title = "Posterior Distributions",
         subtitle = "Dashed: True value, Dotted: BayesBag estimate") +
    scale_fill_brewer(palette = "Set1") +
    theme(legend.position = "bottom")
  
  return(density_plot)
}

density_plot <- create_density_plots(linear_fit$samples, quadratic_fit$samples, bayesbag_fit$bagged_means)
density_plot
#ggsave("posterior_densities.png", density_plot, width = 10, height = 6, dpi = 300, bg = "white")




analyse_fit <- function(fit, true_values, model_name) {
  n_params <- ncol(fit$samples)
  
  # Create parameter names
  param_names <- vector("list", n_params)
  param_names[[1]] <- expression(gamma)
  for(i in 2:n_params) {
    param_names[[i]] <- bquote(beta[.(i-2)])
  }
  
  # Create a wider plot
  #png(paste0(model_name, "_diagnostics.png"), width = 1800, height = 1200)
  
  # Set up 6 columns, 4 rows (2 rows for trace, 2 rows for ACF)
  layout_matrix <- matrix(1:(6*4), nrow=4, ncol=6, byrow=TRUE)
  layout(layout_matrix)
  
  par(mar = c(4, 4, 3, 1), cex.main = 3, cex.lab = 2, cex.axis = 1.5)
  
  # Plot traceplots for all parameters
  for(i in 1:n_params) {
    plot(1:nrow(fit$samples), fit$samples[, i], 
         type = "l", col = "darkblue",
         main = bquote(paste("Traceplot: ", .(param_names[[i]]))),
         xlab = "Iteration", ylab = "Value")
    if(i <= length(true_values)) {
      abline(h = true_values[i], col = "red", lty = 2)
    }
  }
  
  # Fill remaining trace plots with empty plots
  for(i in (n_params+1):(6*2)) {
    plot(1, type="n", axes=FALSE, xlab="", ylab="")
  }
  
  # Plot ACF for all parameters
  for(i in 1:n_params) { 
    acf(fit$samples[, i], 
        main = bquote(paste("ACF: ", .(param_names[[i]]))),
        lag.max = 50, col = "blue")
  }
  
  # Fill remaining ACF plots with empty plots
  for(i in (n_params+1):(6*2)) {
    plot(1, type="n", axes=FALSE, xlab="", ylab="")
  }
  
  #dev.off()
}

true_values <- c(gamma0, beta0)
analyse_fit(linear_fit, true_values, "Linear")
analyse_fit(quadratic_fit, true_values, "Quadratic")
