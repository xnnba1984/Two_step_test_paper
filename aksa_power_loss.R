library(Rcpp)
library(RcppArmadillo)
library(ggplot2)
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

grid <- expand.grid(pi = seq(0, 0.8, 0.1), N=c(60,90,120), method='AKSA', deltaA=1, 
                    deltaB=c(3,4,5), pattern = c("monotone"), KEEP.OUT.ATTRS = FALSE,
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
  method <- par$method
  deltaA <- par$deltaA
  deltaB <- par$deltaB
  pattern <- par$pattern

  for(r in 1:Rreps){
    # simulate data
    nZero <- round(N * pi); nZero
    B     <- c(rep(0, nZero), runif(N - nZero))
    B     <- sample(B)
    Tt    <- sample(rep(0:1, length.out = N))
    Y     <- rnorm(N, sd = 1)
    
    # simulate treatment effect pattern
    if (pattern == "monotone") {
      rk <- rank(B, ties.method = 'min') / N
      Y[Tt == 1 & B > 0] <- Y[Tt == 1 & B > 0] + deltaB * rk[Tt == 1 & B > 0]
    }
    
    if(method=='AKSA'){
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
  }
  result <- data.frame(pAKSA=pAKSAs)
  key <- paste(as.character(par), collapse = "_")
  result.p[[key]] <- result
}
}); pAKSA <- sapply(result.p, function(df) mean(df$pAKSA <= 0.05)); grid$pAKSA <- pAKSA

#write.csv(grid, 'result/aksa_power_drop.csv', row.names = F)
#saveRDS(result.p, 'result/aksa_result_p_power_drop.rds')

grid <- read.csv('result/aksa_power_drop.csv')

# MC error
grid$se_MC   <- sqrt(grid$pAKSA * (1 - grid$pAKSA) / Rreps)
grid$low95   <- pmax(0, grid$pAKSA - 1.96 * grid$se_MC)
grid$upp95   <- pmin(1, grid$pAKSA + 1.96 * grid$se_MC)

ggplot(grid[grid$N!=30,], aes(x = pi, y = pAKSA,
                 colour = factor(N), group = N)) +
  facet_wrap(~ deltaB, nrow = 1,
             labeller = labeller(deltaB = \(x) paste0("Δ = ", x))) +
  geom_line(size = 1.2) +
  annotate("text",
           x = max(grid$pi) - 0.3,   # a tad to the right of the last x-value
           y = 0.80,                  # same y as the line
           label = "Power = 0.80",
           hjust = 0, vjust = -0.4,   # align just above the line
           size = 5,  colour = "black") +
  geom_point(size = 2.5) +
  ## 95 % MC error bars
  geom_errorbar(aes(ymin = low95, ymax = upp95),
                width = 0.01, linewidth = 0.4) +
  labs(
    x = "\nZero inflation rate",
    y = "Power\n",
    colour = "Sample size"
  ) +
  scale_colour_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "black") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "right",panel.spacing   = unit(2, "lines")
  )




