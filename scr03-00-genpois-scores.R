#!/usr/bin/env Rscript

if (!require('yaml')) install.packages('yaml'); library('yaml')
if (!require('doMC')) install.packages('doMC'); library('doMC')
if (!require('stats')) install.packages('stats'); library('stats')


args <- commandArgs(trailingOnly=TRUE)
if(length(args) == 0){
        cat("Please provide a YAML file. \n")
        cat("Usage: Rscript --vanilla scr04-00-enrichment-scores.R sample_library.yaml \n")
        stop()
}

params <- read_yaml(file.path(args[1]))
#params <- read_yaml(file.path("meta", "sample_library_mystat.yaml")) # For testing keep it commented in production


# Functions
lt1 <- 1.0 - .Machine$double.eps

gp_mle <- function(y) {

  fun <- function(par, y, y1, n, sy) {
    theta <- exp( par[1] ) ;  lambda <- 1 / ( 1 + exp(- par[2]) )
    f <- n * log( theta ) + sum( y1 * log( theta + lambda * y ) ) -
         n * theta - sy * lambda
    -f
  }

  sy <- sum(y)
  n <- length(y)
  my <- sy / n
  #s2 <- Rfast::Var(y)
  s2 <- var(y)
  y1 <- y - 1

  theta <- log( sqrt(my^3 / s2 ) )
  f <- 1 - sqrt(my / s2)

  lamb_check <- f / (1 - f)
  if(lamb_check > 0){
          lambda <- log(f / (1 - f))
          ini <- c( theta, lambda )

          mod <- optim(ini, fun, y = y, y1 = y1, n = n, sy = sy, control = list(maxit = 10000))
          lik1 <- mod$value
          mod <- optim(mod$par, fun, y = y, y1 = y1, n = n, sy = sy, control = list(maxit = 10000) )
          lik2 <- mod$value
          while ( lik1 - lik2 > 1e-6 ) {
            lik1 <- lik2
            mod <- optim(mod$par, fun, y = y, y1 = y1, n = n, sy = sy, control = list(maxit = 10000) )
            lik2 <- mod$value
          }
          lik <-  - mod$value - sum( lgamma(y + 1) )
          par <- mod$par
          par <- c( exp( par[1] ), 1 / ( 1 + exp(- par[2]) ), lik )
          names(par) <- c("theta", "lambda", "loglik")
          return(par)
  } else {
          par <- c(NA, NA, NA)
          names(par) <- c("theta", "lambda", "loglik")
          return(par)
          
  }
}

estimate_gp_distributions <- function(input_counts, output_counts, uniq_input_values) {
        cat("Estimating GP distributions\n")
        df_gp_dist <- data.frame()
        uniq_input_values <- sort(uniq_input_values)
        for (input_value in uniq_input_values) {
                curr_counts <- output_counts[input_counts == input_value, 1]
                if (length(curr_counts) < 50) next

                gp_mle_res <- gp_mle(curr_counts)
                cat("The result is:", input_value, gp_mle_res, "\n")
                lambda <- gp_mle_res[["lambda"]]
                theta <- gp_mle_res[["theta"]]
                temp_df <- data.frame("idxs" = input_value, "theta" = theta, "lambda" = lambda)
                df_gp_dist <- rbind(df_gp_dist, temp_df)
        }
        return(df_gp_dist)
}


lambda_theta_regression <- function(df_gp_dist) {
        cat("Conducting lambda theta regression\n")
        lambda_ <- mean(df_gp_dist$lambda)
        lambda_fit <- function(x) lambda_
        coeffs <- lm(df_gp_dist$theta ~ df_gp_dist$idxs)$coefficients
        theta_fit <- function(x) coeffs[2] * x + coeffs[1]
  
        return(list(lambda_fit=lambda_fit, theta_fit=theta_fit))
}


log_GP_pmf <- function(x, theta, lambd) {
        # Ensure arguments to log are positive
        term1 <- ifelse(theta > 0, log(theta), -Inf)
        term2 <- ifelse(theta + x * lambd > 0, (x - 1) * log(theta + x * lambd), -Inf)
        term3 <- theta + x * lambd
        term4 <- lgamma(x + 1)
        result <- term1 + term2 - term3 - term4
        return(result)
}

log_GP_sf <- function(x, theta, lambd) {
        count <- x + 1
        accum <- log_GP_pmf(count, theta, lambd)
        repeat {
                count <- count + 1
                new_term <- log_GP_pmf(count, theta, lambd)
                new <- -Inf  # Default if both terms are -Inf
                # Implement log-sum-exp more safely
                if (is.finite(accum) & is.finite(new_term)) {
                max_term <- max(accum, new_term)
                new <- log(exp(accum - max_term) + exp(new_term - max_term)) + max_term
                }
                # Check for non-finite values before comparison
                if(!is.finite(new) | !is.finite(accum)){
                    break
                }else if(new - accum < 1e-6) {
                break
                }
                accum <- new
                }
        return(accum)
}


precompute_pvals <- function(regression_results, uniq_combos) {
        cat("Computing p-values\n")
        df_log10pval_hash <- data.frame(matrix(nrow =0, ncol = 3))
        colnames(df_log10pval_hash) <- c(colnames(uniq_combos), "log10pval_hash")
        temp_df <- df_log10pval_hash

        for (i in 1:nrow(uniq_combos)) {
                combo <- uniq_combos[i,]
                ic <- combo[[1]]
                oc <- combo[[2]]
                log_pval <- log_GP_sf(oc, regression_results$theta_fit(ic)[[1]], regression_results$lambda_fit(ic))
                log10pval_hash <- log_pval * log10(exp(1)) * -1
                temp_df[1,1] <- ic
                temp_df[1,2] <- oc
                temp_df[1,3] <- log10pval_hash
                df_log10pval_hash <- rbind(df_log10pval_hash, temp_df)
        }
        return(df_log10pval_hash)
}



user_input <- function(prompt) {
  if (interactive()) {
    return(invisible(readline(prompt)))
  } else {
    cat(prompt)
    return(invisible(readLines("stdin", n=1)))
  }
}

# Define input folders
output_normalized_counts_dir <- params$output_normalized_counts_dir


# Create folders
output_generalized_poisson_p_vals_dir <- params$output_generalized_poisson_p_vals_dir

if(dir.exists(output_generalized_poisson_p_vals_dir)){
        cat("Output directory already exists, output maybe over written.: ", output_generalized_poisson_p_vals_dir, "\n")
        user_input("Press [enter] to continue or [ctrl+z] to quit.")
}else{
        cat("Creating output dir: ", output_generalized_poisson_p_vals_dir, "\n")
        dir.create(output_generalized_poisson_p_vals_dir, recursive = TRUE, showWarnings = FALSE)
}

# Load the input_dat from the input file
files <- list.files(path = output_normalized_counts_dir, full.names = TRUE, pattern = "*.tsv")
registerDoMC(detectCores(all.tests = FALSE, logical = FALSE) - 1)
null_out <- foreach(file_i = 1:length(files)) %dopar% {
        input <- files[file_i]
        cat("Processing:", basename(input), "\n")
        input_dat <- read.table(input, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names = FALSE)

        # Extract clones, input counts, and output counts
        input_counts <- as.integer(input_dat[,2])
        output_counts <- as.matrix(input_dat[, -c(1, 2)]) + 1 # Adding pseudocounts

        # Estimate generalized Poisson distributions for every input count
        uniq_input_values <- unique(input_counts)
        distribution_results <- estimate_gp_distributions(input_counts, output_counts, uniq_input_values)
        distribution_results <- na.omit(distribution_results)

        # Regression on all of the theta and lambda values computed
        regression_results <- lambda_theta_regression(distribution_results)

        # Precompute CDF for possible input-output combinations
        uniq_combos <- unique(input_dat[,2:3])
        log10pval_results <- precompute_pvals(regression_results, uniq_combos)

        merge_log10_id <- merge(input_dat, log10pval_results, all.x = TRUE)
        output_dat <- merge_log10_id[order(merge_log10_id$id),]

        samp_name <- colnames(output_dat)[2]

        output_dat <- output_dat[, c("id", "log10pval_hash")]
        colnames(output_dat) <- c("id", samp_name)

        output_file_name <- file.path(output_generalized_poisson_p_vals_dir, 
                                      paste0(tools::file_path_sans_ext(basename(input)), "_mlxp.tsv"))

        write.table(output_dat, file = output_file_name, 
                    sep = "\t", row.names = FALSE, quote = FALSE)
}

