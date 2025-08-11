library(Rcpp)
library(RcppArmadillo)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridisLite)
setwd("~/Library/CloudStorage/Box-Box/Xi/abbvie/AKSA")

Rcpp::cppFunction('
double ks2_fast(NumericVector x, NumericVector y) {
  int nx = x.size(), ny = y.size();
  std::sort(x.begin(), x.end());
  std::sort(y.begin(), y.end());
  int i = 0, j = 0;
  double cdf_x = 0.0, cdf_y = 0.0, dmax = 0.0;

  while (i < nx && j < ny) {
    if (x[i] <= y[j]) {
      ++i;   cdf_x = double(i) / nx;
    } else {
      ++j;   cdf_y = double(j) / ny;
    }
    dmax = std::max(dmax, std::fabs(cdf_x - cdf_y));
  }
  return dmax;
}
')

grid <- expand.grid(pi = c(0,seq(0.1,0.8,0.1)), N=c(60,90,120), deltaA=C(3,4,5), deltaB=c(5),
                    pattern = c("monotone"),
                    KEEP.OUT.ATTRS = FALSE,
                    stringsAsFactors = FALSE)
nperm  <- 1000         # permutations
alpha  <- 0.05         # test level
Rreps <- 1000
result.p <- list()

system.time({
  set.seed(2026)
  for (i in seq_len(nrow(grid))) {
    par <- grid[i, ]
    print(par)
    pi <- par$pi
    N <- par$N
    deltaA <- par$deltaA
    deltaB <- par$deltaB
    pattern <- par$pattern
    pAs <- c()
    pBs <- c()
    pAKSAs <- c()
    
    for(r in 1:Rreps){
      # simulate data
      nZero <- round(N * pi); nZero
      B     <- c(rep(0, nZero), runif(N - nZero))
      B     <- sample(B)
      Tt    <- sample(rep(0:1, length.out = N))
      Y     <- rnorm(N, sd = 1)
      
      # simulate treatment effect pattern
      if (pattern == "spike"){
        Y[B == 0 & Tt == 1] <- Y[B == 0 & Tt == 1] + deltaA
      }
      if (pattern == "monotone") {
        rk <- rank(B, ties.method = 'min') / N
        Y[Tt == 1 & B > 0] <- Y[Tt == 1 & B > 0] + deltaB * rk[Tt == 1 & B > 0]
      }
      if (pattern == "mix"){
        Y[B == 0 & Tt == 1] <- Y[B == 0 & Tt == 1] + deltaA
        rk <- rank(B, ties.method = 'min') / N
        Y[Tt == 1 & B > 0] <- Y[Tt == 1 & B > 0] + deltaB * rk[Tt == 1 & B > 0]
      }
      
      ## ----- Part A -----
      safe_diff <- function(a, b) {
        if (length(a) == 0 || length(b) == 0) return(0)
        abs(mean(a) - mean(b))
      }
      if(all(B>0)){
        pA <- 1
      } else{
        G     <- Tt + 2L * (B > 0)
        obsA  <- safe_diff(Y[G == 1L], Y[G == 0L])
        permA <- vapply(1:nperm, function(b) {
          Y_perm <- Y_hat + perm_res[, b]          
          safe_diff(Y_perm[G == 1L], Y_perm[G == 0L])
        }, numeric(1))
        pA    <- (1 + sum(permA >= obsA)) / (nperm + 1)
      }
      pAs[r] <- pA
      
      ## ----- Part B -----
      pos   <- which(B > 0)
      if (length(pos) < 5) {
        pB <- 1
      } else {
        ord  <- order(B[pos])
        Ypos <- Y[pos][ord];  Tpos <- Tt[pos][ord];  npos <- length(Ypos)
        ks_vec <- sapply(1:(npos - 1), function(k) {
          yk <- Ypos[1:k]; tk <- Tpos[1:k]
          if (all(tk == 0) || all(tk == 1)) return(0)
          ks2_fast(yk[tk == 1], yk[tk == 0])
        })
        AKSAobs <- mean(ks_vec)

        permB <- vapply(1:nperm, function(b) {
          Y_perm <- Y_hat + perm_res[, b]        
          Yp     <- Y_perm[pos][ord]
          d <- sapply(1:(npos - 1), function(k) {
            yk <- Yp[1:k]; tk <- Tpos[1:k]         
            if (all(tk == 0) || all(tk == 1)) return(0)
            ks2_fast(yk[tk == 1], yk[tk == 0])
          })
          mean(d)
        }, numeric(1))
        pB <- (1 + sum(permB >= AKSAobs)) / (nperm + 1)
      }
      pBs[r] <- pB
      
      #------ regular AKSA -------------
      ord  <- order(B)
      Ypos <- Y[ord];  Tpos <- Tt[ord];  npos <- length(Ypos)
      ks_vec <- sapply(1:(npos - 1), function(k) {
        yk <- Ypos[1:k]; tk <- Tpos[1:k]
        if (all(tk == 0) || all(tk == 1)) return(0)
        ks2_fast(yk[tk == 1], yk[tk == 0])
      })
      AKSAobs <- mean(ks_vec)

      permB <- vapply(1:nperm, function(b) {
        Y_perm <- Y_hat + perm_res[, b]        
        Yp     <- Y_perm[ord]
        d <- sapply(1:(npos - 1), function(k) {
          yk <- Yp[1:k]; tk <- Tpos[1:k]         
          if (all(tk == 0) || all(tk == 1)) return(0)
          ks2_fast(yk[tk == 1], yk[tk == 0])
        })
        mean(d)
      }, numeric(1))
      pAKSA <- (1 + sum(permB >= AKSAobs)) / (nperm + 1)
      pAKSAs[r] <- pAKSA
    }
    
    # raw fisher
    pFs <- 1 - pchisq(-2 * log(pAs*pBs), df = 4)
    
    # correction fisher
    rho   <- suppressWarnings(cor(pAs, pBs, method = "spearman", use = "complete.obs"))
    if (is.na(rho)){
      rho <- 0  # constant p-values ⇒ treat as independent
    }
    cfac  <- 1 + rho
    df_b  <- 4 / cfac       
    pCFs  <- 1 - pchisq(-2 * log(pAs*pBs)/cfac, df = df_b)
    
    result <- data.frame(pA=pAs, pB=pBs, pAKSA=pAKSAs, pF=pFs, pCF=pCFs)
    key <- paste(as.character(par), collapse = "_")
    result.p[[key]] <- result
  }
})

grid$pAKSA <- sapply(result.p, function(df) mean(df$pAKSA <= 0.05, na.rm = T))
grid$pA <- sapply(result.p, function(df) mean(df$pA <= 0.05, na.rm = T))
grid$pB <- sapply(result.p, function(df) mean(df$pB <= 0.05, na.rm = T))
grid$pF <- sapply(result.p, function(df) mean(df$pF <= 0.05, na.rm = T))
grid$pCF <- sapply(result.p, function(df) mean(df$pCF <= 0.05, na.rm = T))
grid$pF[grid$pi==0] <- grid$pAKSA[grid$pi==0]
grid$pCF[grid$pi==0] <- grid$pAKSA[grid$pi==0]

#write.csv(grid, 'result/mono_all.csv', row.names = F)
#saveRDS(result.p, 'result/mono_all.rds')

grid <- read.csv('result/mono_all.csv')
# r.p <- readRDS('result/mono_all.rds')
# result.p <- c(result.p, r.p)


heat_df <- grid %>% 
  select(pi, N, deltaB, AKSA = pAKSA, Fisher = pF, Brown = pCF) %>% 
  pivot_longer(AKSA:Brown,
               names_to  = "Method",
               values_to = "Power") %>% 
  mutate(
    pi     = as.numeric(pi),          # treat as discrete for tidy tiles
    N      = as.numeric(N),
    deltaB = factor(deltaB, 
                    levels = sort(unique(deltaB))),  # ensures row order
    Method = factor(Method, levels = c("AKSA", "Fisher", "Brown"))
  )

# plot fisher
plot_df <- heat_df %>% 
  filter(Method != "Brown") %>%              # drop Brown if not plotting
  mutate(se   = sqrt(Power * (1 - Power) / Rreps),        # binomial SE
         lo95 = pmax(0, Power - qnorm(0.975) * se),      # lower limit
         hi95 = pmin(1, Power + qnorm(0.975) * se)) %>%  # upper limit
  mutate(
    N      = factor(N,      levels = c(30, 60, 90, 120)),     # row order
    deltaB = factor(deltaB, levels = sort(unique(deltaB))),
    Method = factor(Method, levels = c("AKSA", "Fisher"))
  )

ggplot(plot_df[plot_df$N!=30,], aes(pi, Power, colour = Method, linetype = Method)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95, fill = Method),
             alpha = 0.3, colour = NA, show.legend = FALSE) +
  geom_line(aes(linewidth = Method)) +
  geom_point(aes(size = Method)) +
  facet_grid(rows = vars(N),
             cols  = vars(deltaB),
             labeller = labeller(deltaB = \(x) paste0("Δ = ", x),
                                 N      = \(x) paste("N =", x))) +
  scale_x_continuous(breaks = seq(0, 0.8, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  
  scale_colour_manual(values = c(Fisher = "#CC79A7",
                                 AKSA   = "grey50")) +
  scale_linetype_manual(values = c(Fisher = "solid",
                                   AKSA   = "dashed")) +
  scale_linewidth_manual(values = c(AKSA   = 0.8,
                                    Fisher = 1), guide = "none") +
  scale_size_manual(values = c(AKSA   = 1.8,
                               Fisher = 2), guide = "none") +
  scale_fill_manual(values = c(Fisher = "#CC79A7",   # ribbon colour
                                 AKSA   = "grey50")) +
  labs(x = 'Zero inflation rate',
       y = "Power",
       colour   = "Method",
       linetype = "Method") +
  theme_bw(base_size = 16) +
  guides(linetype = guide_legend(override.aes = list(linewidth=1.5, size = 2))) +
  geom_hline(yintercept = 0.8, linetype = "dashed") +
  theme(panel.spacing   = unit(1, "lines"),
        legend.position = "bottom",
        legend.title    = element_blank())

# plot brown
plot_df <- heat_df %>% 
  filter(Method != "Fisher") %>%              # drop Brown if not plotting
  mutate(se   = sqrt(Power * (1 - Power) / Rreps),        # binomial SE
         lo95 = pmax(0, Power - qnorm(0.975) * se),      # lower limit
         hi95 = pmin(1, Power + qnorm(0.975) * se)) %>%  # upper limit
  mutate(
    N      = factor(N,      levels = c(30, 60, 90, 120)),     # row order
    deltaB = factor(deltaB, levels = sort(unique(deltaB))),
    Method = factor(Method, levels = c("AKSA", "Brown"))
  )

ggplot(plot_df[plot_df$N!=30,], aes(pi, Power, colour = Method, linetype = Method)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95, fill = Method),
              alpha = 0.3, colour = NA, show.legend = FALSE) +
  geom_line(aes(linewidth = Method)) +
  geom_point(aes(size = Method)) +
  facet_grid(rows = vars(N),
             cols  = vars(deltaB),
             labeller = labeller(deltaB = \(x) paste0("Δ = ", x),
                                 N      = \(x) paste("N =", x))) +
  scale_x_continuous(breaks = seq(0, 0.8, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  
  scale_colour_manual(values = c(Brown = "#56B4E9",
                                 AKSA   = "grey50")) +
  scale_linetype_manual(values = c(Brown = "solid",
                                   AKSA   = "dashed")) +
  scale_linewidth_manual(values = c(AKSA   = 0.8,
                                    Brown = 1), guide = "none") +
  scale_size_manual(values = c(AKSA   = 1.8,
                               Brown = 2), guide = "none") +
  scale_fill_manual(values = c(Brown = "#56B4E9",   # ribbon colour
                               AKSA   = "grey50")) +
  labs(x = 'Zero inflation rate',
       y = "Power",
       colour   = "Method",
       linetype = "Method") +
  theme_bw(base_size = 16) +
  guides(linetype = guide_legend(override.aes = list(linewidth=1.5, size = 2))) +
  geom_hline(yintercept = 0.8, linetype = "dashed") +
  theme(panel.spacing   = unit(1, "lines"),
        legend.position = "bottom",
        legend.title    = element_blank())


