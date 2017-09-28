rankbs <- function(Y, n, H, d, iter, alpha, CPU, seed, resampling){

  a <- length(Y)
  N <- sum(n)

  w <- function(Y, l, i){
    Y_new <- rbind(Y[[l]]$response, Y[[i]]$response)
    Ranks <- apply(Y_new, 2, rank)

    Ra <- list()
    Ra[[l]] <- Ranks[1:n[l], ]
    Ra[[i]] <- Ranks[(n[l]+1):(n[l]+n[i]), ]
    w <- 1/n[l]*(colMeans(Ra[[i]]) - (n[i]+1)/2)
    return(w)
  }

  p <- matrix(NA, ncol = d, nrow = a)

  for (i in 1:a){
    w_tmp <- matrix(NA, nrow = a, ncol = d)
    for (l in 1:a){
      w_tmp[l, ] <- w(Y, l, i)
    }
    p[i, ] <- colMeans(w_tmp)
  }

  p_vec <- c(t(p))

  # Test statistic
  TS <- N* t(p_vec)%*%H%*%p_vec
#-------------------------------------- Bootstrap ------------------
  boot <- function(i, ...){
    index <- list()
    for(i in 1:a){
      index[[i]] <- sample(1:n[i], size = n[i], replace = TRUE)
    }

    Y_tmp <- do.call(rbind, Y)

    Y_b <- list()
    for(i in 1:a){
      Y_b[[i]] <- Y[[i]][index[[i]], ]
    }

    p_boot <- matrix(NA, ncol = d, nrow = a)

    for (i in 1:a){
      w_tmp <- matrix(NA, nrow = a, ncol = d)
      for (l in 1:a){
        w_tmp[l, ] <- w(Y_b, l, i)
      }
      p_boot[i, ] <- colMeans(w_tmp)
    }

    p_vec_boot <- c(t(p_boot))

    TSs <- N* t(p_vec_boot - p_vec)%*%H%*%(p_vec_boot - p_vec)
    return(TSs)
  }

  #-------------------------------------- Wild Bootstrap ------------------
  Z <- function(s, r){
    Y_new <- rbind(Y[[s]]$response, Y[[r]]$response)
    R_sr <- apply(Y_new, 2, rank)
    R_r <- apply(Y[[r]]$response, 2, rank)

    n_s <- n[s]
    n_r <- n[r]

    Ra_r <- R_sr[(n_s+1):(n_s+n_r), ]    # passende n_r Ranks
    Z_sr <- 1/n_s*((Ra_r- R_r) - (colMeans(Ra_r)-(n_r+1)/2))
    return(Z_sr)
  }

  # Wild BS
  wildBS <- function(i, ...){
    epsilon <- function(l){
      epsilon <- matrix(2*rbinom(n[l], 1, 1/2)-1, n[l], d)
      return(epsilon)
    }

    Z_star <- function(l, i){
      Z_star <- colMeans(epsilon(i)*Z(l, i))- colMeans(epsilon(l)*Z(i, l))
      return(Z_star)
    }

    p_boot <- matrix(NA, ncol = d, nrow = a)
    for (i in 1:a){
      Z_tmp <- matrix(NA, nrow = a, ncol = d)
      for (l in 1:a){
        Z_tmp[l, ] <- Z_star(l, i)
      }
      p_boot[i, ] <- colMeans(Z_tmp)
    }

    p_vec_boot <- c(t(p_boot))

    TSs <- N* t(p_vec_boot)%*%H%*%p_vec_boot
    return(TSs)
  }


  #-----------------------------------------------------------------------#

  cl <- parallel::makeCluster(CPU)
  if(seed != 0){
    parallel::clusterSetRNGStream(cl, iseed = seed)
  }

  if(resampling == "bootstrap"){
    bs <- parallel::parSapply(cl, 1:iter, FUN = boot)
  } else if (resampling == "WildBS"){
    bs <- parallel::parSapply(cl, 1:iter, FUN = wildBS)
  }
  quant <- quantile(unlist(bs), prob = 1-alpha)
  ecdf <- ecdf(unlist(bs))
  p_value <- 1-ecdf(TS)

  parallel::stopCluster(cl)

  result <- list()
  result$statistic <- c(TS, p_value)
  result$p <- p

  return(result)
}
