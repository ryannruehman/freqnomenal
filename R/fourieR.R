

frequencies <- function(signal, time){
  if (anyNA(signal)){
    stop("signal input contains NA values")
  }
  if (anyNA(time)){
    stop("time input contains NA values")
  }
  if (length(signal) != length(time)){
    stop("signal and time must have the same number of data entries")
  }
  checkunif <- diff(time)
  if (max(abs(checkunif - mean(checkunif))) > 0.001 * mean(checkunif)){
    stop("signal must be uniformly sampled")
  }
  n <- length(signal)
  tstar <- time[2] - time[1]
  Fund_freq <- 1 / (n * tstar)
  N <- n / 2
  N_list <- seq(1, N)
  f_0 <- 0
  f_ny <- (n/2) * Fund_freq
  f_klist <- seq(1, N - 1)
  freqs <- numeric(N - 1)
  for(i in seq_along(f_klist)){
    freqs[i] <- f_klist[i] * Fund_freq
  }
  frequencies <- data.frame(Frequencies = c(f_0, freqs, f_ny))
  return(frequencies)
}


dft <- function(signal, time){
  if (anyNA(signal)){
    stop("signal input contains NA values")
  }
  if (anyNA(time)){
    stop("time input contains NA values")
  }
  if (length(signal) != length(time)){
    stop("signal and time must have the same number of data entries")
  }
  checkunif <- diff(time)
  if (max(abs(checkunif - mean(checkunif))) > 0.001 * mean(checkunif)){
    stop("signal must be uniformly sampled")
  }
  n <- length(signal)
  tstar <- time[2] - time[1]
  Fund_freq <- 1 / (n * tstar)
  N <- n / 2
  N_list <- seq(1, N)
  f_0 <- 0
  f_ny <- (n/2) * Fund_freq
  f_klist <- seq(1, N - 1)
  freqs <- numeric(N - 1)
  for(i in seq_along(f_klist)){
    freqs[i] <- f_klist[i] * Fund_freq
  }
  chi_cos <- matrix(NA, nrow = length(time), ncol = length(freqs))
  for (i in 1:nrow(chi_cos)){
    for (j in 1:ncol(chi_cos)){
      chi_cos[i, j] <- cos(2 * pi * freqs[j] * time[i])
    }
  }

  chi_sin <- matrix(NA, nrow = length(time), ncol = length(freqs))
  for (i in 1:nrow(chi_sin)){
    for (j in 1:ncol(chi_sin)){
      chi_sin[i, j] <- sin(2 * pi * freqs[j] * time[i])
    }
  }

  chi_ny <- seq(1:n)
  for (i in seq_along(time)){
    if (i %% 2 == 0) chi_ny[i] <- -1
    else chi_ny[i] <- 1
  }

c_0 <- sum(signal) / length(signal)
a_k <- colSums(signal * chi_cos) / colSums(chi_cos * chi_cos)
b_k <- colSums(signal * chi_sin) / colSums(chi_sin * chi_sin)
c_ny <- sum(signal * chi_ny) / sum(chi_ny * chi_ny)

dftresults <- data.frame(c_0 = c_0,
                         a_k = a_k,
                         b_k = b_k,
                         c_ny = c_ny)
return(dftresults)

}

psd <- function(signal, time){
  if (anyNA(signal)){
    stop("signal input contains NA values")
  }
  if (anyNA(time)){
    stop("time input contains NA values")
  }
  if (length(signal) != length(time)){
    stop("signal and time must have the same number of data entries")
  }
  checkunif <- diff(time)
  if (max(abs(checkunif - mean(checkunif))) > 0.001 * mean(checkunif)){
    stop("signal must be uniformly sampled")
  }
  n <- length(signal)
  tstar <- time[2] - time[1]
  Fund_freq <- 1 / (n * tstar)
  N <- n / 2
  N_list <- seq(1, N)
  f_0 <- 0
  f_ny <- (n/2) * Fund_freq
  f_klist <- seq(1, N - 1)
  freqs <- numeric(N - 1)
  for(i in seq_along(f_klist)){
    freqs[i] <- f_klist[i] * Fund_freq
  }


  chi_cos <- matrix(NA, nrow = length(time), ncol = length(freqs))
  for (i in 1:nrow(chi_cos)){
    for (j in 1:ncol(chi_cos)){
      chi_cos[i, j] <- cos(2 * pi * freqs[j] * time[i])
    }
  }

  chi_sin <- matrix(NA, nrow = length(time), ncol = length(freqs))
  for (i in 1:nrow(chi_sin)){
    for (j in 1:ncol(chi_sin)){
      chi_sin[i, j] <- sin(2 * pi * freqs[j] * time[i])
    }
  }

  chi_ny <- seq(1:n)
  for (i in seq_along(time)){
    if (i %% 2 == 0) chi_ny[i] <- -1
    else chi_ny[i] <- 1
  }

  c_0 <- sum(signal) / length(signal)
  a_k <- colSums(signal * chi_cos) / colSums(chi_cos * chi_cos)
  b_k <- colSums(signal * chi_sin) / colSums(chi_sin * chi_sin)
  c_ny <- sum(signal * chi_ny) / sum(chi_ny * chi_ny)


  psdf_0 <- (c_0)^2 / n
  psdf_k <- ((a_k)^2 + (b_k)^2) / n
  psdf_ny <- (c_ny)^2 / n


psdresult <- data.frame(PSD = c(psdf_0, psdf_k, psdf_ny),
                        Frequency = c(f_0, freqs, f_ny))


return(psdresult)
}






















