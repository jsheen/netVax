# From Lee's code
t.test.calcs <- function(sd1, sd2, alpha=0.05, N=NULL, delta=NULL, pow=NULL) {
  if (is.null(N)+is.null(delta)+is.null(pow)+is.null(alpha) != 1) {
    print("Error: Need exactly one of N, delta, power, or alpha to be null")
    return(NA)
  } else if (is.null(N)) {
    approx <- ((sd1^2+sd2^2)*(qnorm(1-alpha/2)+qnorm(pow))^2)/(delta^2)
    if (approx < 1) {
      return(1)
    } else {
      return(uniroot(f=function(x) ((sd1^2+sd2^2)*(qt(p=1-alpha/2, df=2*x-2)+qt(p=pow, df=2*x-2))^2)/(delta^2)-x,
                     lower=approx, upper=approx*5)$root)
    }
  } else if(is.null(pow)) {
    return(pow <- pt(q=sqrt((N*delta^2)/(sd1^2+sd2^2))-qt(p=1-alpha/2, df=2*N-2), df=2*N-2))
  } else if(is.null(delta)) {
    return(sqrt(((sd1^2+sd2^2)*(qt(p=1-alpha/2, df=2*N-2)+qt(p=pow, df=2*N-2))^2)/N))
  } else {
    return(2*(1-pt(q=sqrt((N*delta^2)/(sd1^2+sd2^2))-qt(p=pow, df=2*N-2), df=2*N-2)))
  }
}

# Function to get final size
getFinalSize <- function(q, i0, r0, gamma, phi) {
  # Make relation
  s_init <- 1 - i0 - q
  s_infs <- seq(0, s_init, 0.00001)
  r0_func <- function(q, i0, gamma, phi, s_inf) {
    num <- (s_init) * log(((1 - gamma) * (s_init)) / s_inf)
    denom <- ((1 - gamma) * (s_init) - s_inf) + gamma * (s_init) * (1 - (s_inf / ((1 - gamma) * (s_init)))^phi) + i0
    return(num / denom)
  }
  r0s <- sapply(s_infs, FUN=r0_func, q=q, i0=i0, gamma=gamma, phi=phi)
  # Get final size proportions
  if (r0 == 0) { # If the infection has an r0 of 0, then it will not spread
    s_final <- 1 - i0 - q # prop. suscep = 1 - prop. init. infect - prop. init. recov.
  } else {
    s_final <- s_infs[which(r0s < r0)[1]] # prop. of the initially susceptible who get infected
  }
  v_final <- gamma * (s_init) * (s_final / ((1 - gamma) * (s_init)))^phi
  # Prop. recovereds of:
  # (1) Pre-existing recovereds
  # (2) naive susceptibles who got infecteds
  # (3) vaccinateds who got infecteds
  # (4) initial infections
  r_final <- 1 - s_final - v_final
  return(c(s_final, v_final, r_final))
}

# Function to get difference in final sizes
getDiffFinalSize <- function(q, i0, r0t, r0w, alpha, r0v, phi_trad, phi_trans) {
  if(q + alpha + i0 > 1) {
    return('q + alpha + i0 cannot be greater than 1.')
  } else {
    # 0) Set up global variables
    s_init = 1 - i0 - q
    # 1) Get final size of traditional vaccine
    trad_vax_res <- getFinalSize(q=q, i0=i0, r0=r0w, gamma=alpha / s_init, phi=phi_trad)
    # 2a) Get final size of "vaccine transmission." When phi=0, 
    #     The recovered prop. is actually the prop. vaccinated. 
    trans_vax_pre <- getFinalSize(q=q, i0=alpha, r0=r0v, gamma=0, phi=0)
    prop_pop_vax <- trans_vax_pre[3] - q # recovereds (made up of vaccine and init recovered) - init. recov = prop. of entire pop. vaccinated.
    # 2b) Get final size after wildtype transmission
    trans_vax_res <- getFinalSize(q=q, i0=i0, r0=r0w, gamma=prop_pop_vax / s_init, phi=phi_trans)
    s_inf_trad <- trad_vax_res[1]
    s_T_inf_trad <- trad_vax_res[2]
    s_inf_trans <- trans_vax_res[1]
    s_T_inf_trans <- trans_vax_res[2]
    # 3) Define f_trad, f_trans, f_v, f_trad_vax, f_trans_vax, f_trad_unvax, f_trans_unvax
    f_trad <- 1 - s_T_inf_trad - s_inf_trad
    f_trans <- 1 - s_T_inf_trans - s_inf_trans
    f_v <- prop_pop_vax
    f_trad_vax <- alpha - s_T_inf_trad
    f_trans_vax <- f_v - s_T_inf_trans
    f_trad_unvax <- f_trad - f_trad_vax
    f_trans_unvax <- f_trans - f_trans_vax
    # Define estimands
    ## Delta 1
    k_delta1 <- 1
    delta1 <- (f_trad / k_delta1) - (f_trans / k_delta1)
    ## Delta 4
    delta4 <- (f_trad_vax / alpha) - (f_trans_vax / f_v)
    if (delta4 < 1e-10) {
      delta4 <- 0
    }
    ## Delta 5
    delta5 <- (f_trad_unvax / (1 - alpha)) - (f_trans_unvax / (1 - alpha))
    ## Delta 2
    delta2 <- (f_trad_unvax / (1 - alpha)) - (f_trans_vax / f_v)
    ## Delta 3
    delta3 <- (f_trad_unvax / (1 - alpha)) - ((f_trans - f_trans_vax) / (1 - f_v))
    return(c(delta1, delta2, delta3, delta4, delta5,
             f_trad, f_trans, f_v, f_trad_vax, f_trans_vax, f_trad_unvax, f_trans_unvax, #12
             f_v - alpha, 1 - f_v))
  }
}

getSampSize <- function(n, q, i0, r0t, r0w, alpha, r0v, phi_trad, phi_trans, sampN, sigma_trad, diff_init_vax_unvax) {
  sigma_trans = sigma_trad
  s <- sampN
  # Get final sizes
  res <- getDiffFinalSize(q, i0, r0t, r0w, alpha, r0v, phi_trad, phi_trans)
  delt1 <- res[1]
  delt2 <- res[2]
  delt3 <- res[3]
  delt4 <- res[4]
  delt5 <- res[5]
  f_trad <- res[6]
  f_trans <- res[7]
  f_v <- res[8]
  f_trad_vax <- res[9]
  f_trans_vax <- res[10]
  f_trad_unvax <- res[11]
  f_trans_unvax <- res[12]
  # Delta 1. Overall.
  var_trad_delt1 <- (n - s) / (s*(n-1)) * (f_trad * (1 - f_trad)) + 
    (n / (n - 1)) * ((s - 1) / s) * as.numeric(sigma_trad)
  var_trans_delt1 <- (n - s) / (s*(n-1)) * (f_trans * (1 - f_trans)) + 
    (n / (n - 1)) * ((s - 1) / s) * as.numeric(sigma_trans)
  req_samp_delt1 <- 2*ceiling(t.test.calcs(sd1=sqrt(var_trad_delt1), sd2=sqrt(var_trans_delt1), alpha=0.05, N=NULL, delta=delt1, pow=0.8))
  # Delta 4. Init vax.
  if (diff_init_vax_unvax) { # if can differentiate, then can sample 's' from 'k'
    s_delt4 <- s
    if (s_delt4 > n * alpha) { # if larger than 'k', sample all 'k'
      s_delt4 <- n * alpha
    }
  } else { # if cannot, can only sample 's*alpha' from 'k'
    s_delt4 <- s * alpha
  }
  f_trad4 <- f_trad_vax / alpha
  f_trans4 <- f_trans_vax / f_v
  var_trad_delt4 <- ((n * alpha) - s_delt4) / (s_delt4 * (n * alpha - 1)) * (f_trad4 * (1 - f_trad4)) +
    (n * alpha) / ((n * alpha) - 1) * ((s_delt4 - 1) / s_delt4) * var_trad_delt1 * (f_trad4 / f_trad)^2
  var_trans_delt4 <- ((n * alpha) - s_delt4) / (s_delt4 * (n * alpha - 1)) * (f_trans4 * (1 - f_trans4)) +
    (n * alpha) / ((n * alpha) - 1) * ((s_delt4 - 1) / s_delt4) * var_trans_delt1 * (f_trans4 / f_trans)^2
  if (delt4 == 0) {
    req_samp_delt4 <- 1e9
  } else {
    req_samp_delt4 <- 2*ceiling(t.test.calcs(sd1=sqrt(var_trad_delt4), sd2=sqrt(var_trans_delt4), alpha=0.05, N=NULL, delta=delt4, pow=0.8))
  }
  # Delta 5. Init unvax.
  if (diff_init_vax_unvax) {
    s_delt5 <- s
  } else {
    s_delt5 <- s * (1 - alpha)
  }
  f_trad5 <- f_trad_unvax / (1 - alpha)
  f_trans5 <- f_trans_unvax / (1 - alpha)
  var_trad_delt5 <- (n * (1 - alpha) - s_delt5) / (s_delt5 * (n * (1 - alpha) - 1)) * (f_trad5 * (1 - f_trad5)) +
    (n * (1 - alpha)) / (n * (1 - alpha) - 1) * ((s_delt5 - 1) / s_delt5) * var_trad_delt1 * (f_trad5 / f_trad)^2
  var_trans_delt5 <- (n * (1 - alpha) - s_delt5) / (s_delt5 * (n * (1 - alpha) - 1)) * (f_trans5 * (1 - f_trans5)) +
    (n * (1 - alpha)) / (n * (1 - alpha) - 1) * ((s_delt5 - 1) / s_delt5) * var_trans_delt1 * (f_trans5 / f_trans)^2
  req_samp_delt5 <- 2*ceiling(t.test.calcs(sd1=sqrt(var_trad_delt5), sd2=sqrt(var_trans_delt5), alpha=0.05, N=NULL, delta=delt5, pow=0.8))
  # Delta 2. Indir. vax.
  if (diff_init_vax_unvax) { # this s_delt2 is specific to the trans. arm.
    s_delt2 <- s * ((f_v - alpha) / (1 - alpha))
  } else {
    s_delt2 <- s * (f_v - alpha)
  }
  f_trans2 <- f_trans_vax / f_v
  var_trans_delt2 <- (n * (1 - alpha) - s_delt2) / (s_delt2 * (n * (1 - alpha) - 1)) * (f_trans2 * (1 - f_trans2)) +
    (n * (1 - alpha)) / (n * (1 - alpha) - 1) * ((s_delt2 - 1) / s_delt2) * var_trans_delt1 * (f_trans2 / f_trans)^2
  req_samp_delt2 <- 2*ceiling(t.test.calcs(sd1=sqrt(var_trad_delt5), sd2=sqrt(var_trans_delt2), alpha=0.05, N=NULL, delta=delt2, pow=0.8))
  # Delta 3. Never vax.
  if (diff_init_vax_unvax) { # this s_delt3 is specific to the trans. arm.
    s_delt3 <- s * ((1 - f_v) / (1 - alpha))
  } else {
    s_delt3 <- s * (1 - f_v)
  }
  f_trans3 <- f_trans_vax / f_v
  var_trans_delt3 <- (n * (1 - alpha) - s_delt3) / (s_delt3 * (n * (1 - alpha) - 1)) * (f_trans3 * (1 - f_trans3)) +
    (n * (1 - alpha)) / (n * (1 - alpha) - 1) * ((s_delt3 - 1) / s_delt3) * var_trans_delt1 * (f_trans3 / f_trans)^2
  req_samp_delt3 <- 2*ceiling(t.test.calcs(sd1=sqrt(var_trad_delt5), sd2=sqrt(var_trans_delt3), alpha=0.05, N=NULL, delta=delt3, pow=0.8))
  
  return(req_samp_delt1)
}

n = 1000
q = 0.1
i0 = 1/1000
r0t = 0
r0w = 2
alpha = 0.05
r0v = 0.9
phi_trad = 0.2
phi_trans = 0.2
sampN = 100
sigma_trad = 0.01
diff_init_vax_unvax = FALSE

# Panel 1
pan1_xs <- seq(0.5, 1.15, by=0.05)
pan1res1 <- sapply(pan1_xs, getSampSize, n = n, q = q, i0 = i0, r0t = r0t,
                   r0w = r0w, alpha = alpha, phi_trad = phi_trad,
                   phi_trans = phi_trans, sampN = sampN, sigma_trad = sigma_trad, diff_init_vax_unvax = FALSE)

# Panel 2
pan2_xs <- seq(10, 100, by=1)
pan2res1 <- sapply(pan2_xs, getSampSize, n = n, q = q, i0 = i0, r0t = r0t,
                   r0w = r0w, alpha = alpha, r0v = 0.9, phi_trad = phi_trad,
                   phi_trans = phi_trans, sigma_trad = sigma_trad, diff_init_vax_unvax = FALSE)
pan2res2 <- sapply(pan2_xs, getSampSize, n = n, q = q, i0 = i0, r0t = r0t,
                   r0w = r0w, alpha = alpha, r0v = 1.1, phi_trad = phi_trad,
                   phi_trans = phi_trans, sigma_trad = sigma_trad, diff_init_vax_unvax = FALSE)

# Panel 3
pan3_xs <- c(0.0001, 0.001, 0.01, 0.1)
pan3res1 <- sapply(pan3_xs, getSampSize, n = n, q = q, i0 = i0, r0t = r0t,
                   r0w = r0w, alpha = alpha, r0v = 0.9, phi_trad = phi_trad,
                   phi_trans = phi_trans, sampN = sampN, diff_init_vax_unvax = FALSE)
pan3res2 <- sapply(pan3_xs, getSampSize, n = n, q = q, i0 = i0, r0t = r0t,
                   r0w = r0w, alpha = alpha, r0v = 1.1, phi_trad = phi_trad,
                   phi_trans = phi_trans, sampN = sampN, diff_init_vax_unvax = FALSE)

# Plot fig
png(file="~/netVax/code_output/figs/fig3.png",
    width=8, height=3, units="in", res=400)
par(mfrow = c(1, 3))
plot(pan1_xs, pan1res1, type='l', 
     ylab='Required number of clusters',
     xlab=expression('R'['0,v']),
     col='black', lwd=3, cex.lab=1.5)
plot(pan2_xs, pan2res1, type='l', 
     ylab='',
     xlab='Number sampled from each cluster',
     col='blue', lwd=3, ylim=c(2, 30), cex.lab=1.25)
lines(pan2_xs, pan2res2, col='red', lwd=3)
legend(x='topright', legend=c(expression('R'['0,v']*'=0.9'), expression('R'['0,v']*'=1.1')),
       col=c("blue", "red"), lty=1, lwd=3, cex=1.2, bty='n')
plot(pan3_xs, pan3res1, type='l', 
     ylab='',
     xlab='Between cluster variance',
     col='blue', lwd=3, cex.lab=1.5)
lines(pan3_xs, pan3res2, col='red', lwd=3)
legend(x='topleft', legend=c(expression('R'['0,v']*'=0.9'), expression('R'['0,v']*'=1.1')),
       col=c("blue", "red"), lty=1, lwd=3, cex=1.2, bty='n')
dev.off()

