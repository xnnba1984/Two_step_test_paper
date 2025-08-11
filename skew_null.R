library(Rcpp)
library(RcppArmadillo)
library(ggplot2)
library(dplyr)
library(tidyr)
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

grid <- expand.grid(pi = seq(0,0.8,0.1), N=c(60,90,120), deltaA=1, deltaB=3,
                    pattern = c("null"), skew=c('mild','moderate','strong'),
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
    skew <- par$skew
    pAs <- c()
    pBs <- c()
    pAKSAs <- c()
    
    for(r in 1:Rreps){
      # simulate data
      nZero <- round(N * pi); nZero
      if(skew=='mild'){
        B     <- c(rep(0, nZero), rbeta(N - nZero, shape1 = 2, shape2 = 5))
      }
      if(skew=='moderate'){
        B     <- c(rep(0, nZero), rbeta(N - nZero, shape1 = 1, shape2 = 4))
      }
      if(skew=='strong'){
        B     <- c(rep(0, nZero), rbeta(N - nZero, shape1 = 0.5, shape2 = 3))
      }
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
        permA <- replicate(nperm, {
          gperm <- sample(G)
          safe_diff(Y[gperm == 1L], Y[gperm == 0L])
        })
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
        permB   <- replicate(nperm, {
          tperm <- sample(Tpos)
          d <- sapply(1:(npos - 1), function(k) {
            yk <- Ypos[1:k]; tk <- tperm[1:k]
            if (all(tk == 0) || all(tk == 1)) return(0)
            ks2_fast(yk[tk == 1], yk[tk == 0])
          })
          mean(d)
        })
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
      permB   <- replicate(nperm, {
        tperm <- sample(Tpos)
        d <- sapply(1:(npos - 1), function(k) {
          yk <- Ypos[1:k]; tk <- tperm[1:k]
          if (all(tk == 0) || all(tk == 1)) return(0)
          ks2_fast(yk[tk == 1], yk[tk == 0])
        })
        mean(d)
      })
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

write.csv(grid, 'result/type_I_skew.csv', row.names = F)
saveRDS(result.p, 'result/result_p_type_I_skew.rds')

# visualize type I control
df_long <- grid[grid$N!=30,] %>%     
  select(pi, N, skew,AKSA = pAKSA, Fisher = pF, Brown = pCF) %>% 
  pivot_longer(AKSA:Brown, names_to = "Method", values_to = "err")

df_long <- df_long %>% 
  mutate(se  = sqrt(err * (1 - err) / Rreps),
         lo  = pmax(0, err - 1.96 * se),
         hi  = pmin(1, err + 1.96 * se))

df_long$Method <- factor(df_long$Method,
                         levels = c("AKSA", "Fisher", "Brown"))
df_plot <- df_long[df_long$skew=='strong',]

ggplot(df_plot, aes(x = pi, y = err)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25, fill = "#7baac7") +
  geom_line(size = 1.2, colour = "#1f78b4") +
  geom_point(size = 2.5, colour = "#0d4f7d") +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey30") +
  facet_grid(rows = vars(N),          # three rows: N = 60, 90, 120
             cols = vars(Method),     # three columns: AKSA, Fisher, Brown
             labeller = labeller(
               N = function(x) paste("N =", x))) +
  labs(x = '\nZero inflation rate (π)',
       y = "Type-I error\n") +
  #scale_x_continuous(breaks = seq(0,0.08,0.02)) +
  scale_y_continuous(limits = c(0.01,0.09), breaks = seq(0.01, 0.09, by = 0.02)) + 
  theme_minimal(base_size = 16) +
  theme(strip.text = element_text(face = "bold", size = 16),
        panel.spacing.x = unit(2, "lines"),
        panel.spacing.y = unit(2.5, "lines"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"))

