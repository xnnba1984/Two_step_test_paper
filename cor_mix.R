library(Rcpp)
library(RcppArmadillo)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridisLite)
library(scales)
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

grid <- expand.grid(pi = seq(0.1,0.8,0.2), N=c(60, 90, 120), deltaA=0.8, deltaB=2,
                    pattern = c("mix"), k=c(0.5, 1, 2, 3),
                    KEEP.OUT.ATTRS = FALSE,
                    stringsAsFactors = FALSE)
nperm  <- 1000         # permutations
alpha  <- 0.05         # test level
Rreps <- 1000
result.p <- list()
#tune_sigma <- function(delta, rho)  delta * rho / 2

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
    rho_target <- par$rho
    k <- par$k
    
    for(r in 1:Rreps){
      # simulate data
      Z <- abs(rnorm(1))
      nZero <- round(N * pi); nZero
      B     <- c(rep(0, nZero), runif(N - nZero))
      B     <- sample(B)
      Tt    <- sample(rep(0:1, length.out = N))
      Y     <- rnorm(N, sd = 1)
      
      # simulate treatment effect pattern
      if (pattern == "spike"){
        Y[B == 0 & Tt == 1] <- Y[B == 0 & Tt == 1] + deltaA + k * Z
        Y[B > 0 & Tt == 1] <- Y[B > 0 & Tt == 1] + k * Z
      }
      if (pattern == "monotone") {
        rk <- rank(B, ties.method = 'min') / N
        Y[Tt == 1 & B > 0] <- Y[Tt == 1 & B > 0] + deltaB * rk[Tt == 1 & B > 0] + k*Z
        Y[Tt == 1 & B==0] <- Y[Tt == 1 & B==0] + k*Z 
      }
      if (pattern == "mix"){
        Y[B == 0 & Tt == 1] <- Y[B == 0 & Tt == 1] + deltaA + k*Z
        rk <- rank(B, ties.method = 'min') / N
        Y[Tt == 1 & B > 0] <- Y[Tt == 1 & B > 0] + deltaB * rk[Tt == 1 & B > 0] + k*Z
      }
      
      ## ----- Part A -----
      safe_diff <- function(a, b) {
        if (length(a) == 0 || length(b) == 0) return(0)
        abs(mean(a) - mean(b))
        #mean(a) - mean(b)
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

# compute cor
sapply(result.p, function(df) cor(df$pA, df$pB, method = 'spearman'))
grid$cor <- sapply(result.p, function(df) cor(df$pA, df$pB, method = 'spearman'))
grid$cor_round <- round(grid$cor, 1)

# compute fisher power based on adjusted pvalue cutoff
p.c <- read.csv('result/cor_null.csv')
p.c <- p.c[p.c$rho!=0, c('rho','critF')]
pF.adj <- c()
for(i in 1:length(result.p)){
  rho <- grid$cor_round[i]
  cutoff <- p.c$critF[p.c$rho==rho]
  pF.adj[i] <- mean(result.p[[i]]$pF<cutoff)
}
grid$pF_adj <- pF.adj
grid$diff <- grid$pCF - grid$pF_adj
table(grid$cor_round[grid$N==60])
table(grid$cor_round[grid$N==90])
table(grid$cor_round[grid$N==120])

g1 <- grid
r.p1 <- result.p

write.csv(grid, 'result/cor_mix.csv', row.names = F)
saveRDS(result.p, 'result/cor_mix_p.rds')

summary_df <- g1 %>%                                    
  group_by(N, pattern, cor_round) %>%                             # 1) group by N and rounded ρ
  summarise(                                             # 2) average the power-gap
    diff = mean(diff, na.rm = TRUE),                     #    (use na.rm if needed)
    .groups = "drop"                                     # 3) drop grouping metadata
  )
summary_df <- summary_df[summary_df$cor_round%in%c(.1,.2,.3,.4,.5),]

ggplot(summary_df, aes(x = cor_round, y = diff*100, group = 1)) +                               
  geom_hline(yintercept = 0, linetype = "dashed") +    # diff = 0 reference
  geom_line(colour = "#1f78b4", size = 1) +
  geom_point(colour = "#1f78b4", size = 2.5) +
  facet_grid(rows = vars(N), cols = vars(pattern)) +   # 3 rows × 1 col
  labs(x = "\nSpearman correlation",
       y = "Power difference\n(Brown – Fisher)\n") +
  scale_x_continuous(breaks = sort(unique(summary_df$cor_round))) +
  scale_y_continuous(limits = c(-3,3), labels = number_format(accuracy = 1, suffix = "%")) +
  theme_minimal(base_size = 16) +
  theme(panel.spacing = unit(1.2, "lines"),
        strip.text    = element_text(face = "bold"),
        axis.text     = element_text(colour = "black"),
        axis.title    = element_text(colour = "black"))







