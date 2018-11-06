rankbs <- function(Y, n, H, d, iter, alpha, CPU, seed, resampling){

  a <- length(Y)
  N <- sum(n)

  w <- function(Y, l, i){
    if (d==1){
      Y_new <- c(Y[[l]]$response, Y[[i]]$response)
    Ranks <- rank(Y_new)
    Ra <- list()
    Ra[[l]] <- Ranks[1:n[l]]
    Ra[[i]] <- Ranks[(n[l]+1):(n[l]+n[i])]
    w <- as.matrix(1/n[l]*(mean(Ra[[i]]) - (n[i]+1)/2))

      } else {
    Y_new <- rbind(Y[[l]]$response, Y[[i]]$response)
    Ranks <- apply(Y_new, 2, rank)
    Ra <- list()
    Ra[[l]] <- Ranks[1:n[l], ]
    Ra[[i]] <- Ranks[(n[l]+1):(n[l]+n[i]), ]
    w <- 1/n[l]*(colMeans(Ra[[i]]) - (n[i]+1)/2)
    }
    return(w)
  }

  w_t <- sapply(1:a, function(x) sapply(1:a, function(y) w(Y, y, x), simplify = "array"), simplify = "array")
  # array: w(l, i) = w_t[, l, i]
  if(d==1){
    p_vec <- colMeans(w_t)
  } else {
  p_vec <- c(apply(w_t, 3, rowMeans))
}

  # Test statistic
  TS <- N* t(p_vec)%*%H%*%p_vec
#-------------------------------------- Bootstrap ------------------
  boot <- function(i, ...){
    index <- list()
    for(i in 1:a){
      index[[i]] <- sample(1:n[i], size = n[i], replace = TRUE)
    }

    Y_b <- list()
    for(i in 1:a){
      Y_b[[i]] <- Y[[i]][index[[i]], ]
    }


    w_t_boot <- sapply(1:a, function(x) sapply(1:a, function(y) w(Y_b, y, x), simplify = "array"), simplify = "array")

    if(d==1){
      p_vec_boot <- colMeans(w_t_boot)
    } else {
    p_vec_boot <- c(apply(w_t_boot, 3, rowMeans))
}
    TSs <- N* t(p_vec_boot - p_vec)%*%H%*%(p_vec_boot - p_vec)
    return(TSs)
  }

  #-------------------------------------- Wild Bootstrap ------------------
  Z <- function(s, r){

    if(d == 1){
      Y_new <- c(Y[[s]]$response, Y[[r]]$response)
      R_sr <- rank(Y_new)
      R_r <- rank(Y[[r]]$response)

      n_s <- n[s]
      n_r <- n[r]

      Ra_r <- R_sr[(n_s+1):(n_s+n_r)]    # passende n_r Ranks
      Z_sr <- as.matrix(1/n_s*(Ra_r- R_r) - w_t[r, s])
    } else {

    Y_new <- rbind(Y[[s]]$response, Y[[r]]$response)
    R_sr <- apply(Y_new, 2, rank)
    R_r <- apply(Y[[r]]$response, 2, rank)

    n_s <- n[s]
    n_r <- n[r]

    Ra_r <- R_sr[(n_s+1):(n_s+n_r), ]    # passende n_r Ranks
    Z_sr <- 1/n_s*(Ra_r- R_r) - w_t[, r, s]
    }
    return(Z_sr)
  }


  epsilon <- function(l){
    eps <- 2*rbinom(n[l], 1, 1/2)-1
    return(eps)
  }

  Z_star <- function(l, i, epsi){
    Z_star <- colMeans(unlist(epsi[l])*Z(i, l))- colMeans(unlist(epsi[i])*Z(l, i))
    return(Z_star)
  }

  # Wild BS
  wildBS <- function(...){

    epsi <- lapply(1:a, epsilon)
    Z_t <- sapply(1:a, function(x) sapply(1:a, function(y) Z_star(y, x, epsi)), simplify = "array")

    if(d==1){
      p_vec_boot <- colMeans(Z_t)
    } else{
    p_vec_boot <- c(apply(Z_t, 3, rowMeans))
    }
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
  result$p <- matrix(p_vec, ncol = d, byrow = TRUE)

  return(result)
}
