library(Rcpp)
library(RcppArmadillo)
library(ggplot2)
library(dplyr)
library(tidyr)
library(MASS) 
setwd("~/Library/CloudStorage/Box-Box/Xi/abbvie/AKSA")

grid <- expand.grid(rho=seq(0,0.8,0.1),
                    KEEP.OUT.ATTRS = FALSE,
                    stringsAsFactors = FALSE)
Rreps <- 10000
result.p <- list()

system.time({
  set.seed(2026)
  for (i in seq_len(nrow(grid))) {
    rho_target <- grid[i, ]
    pAs <- c()
    pBs <- c()
    pAKSAs <- c()
    
    Sigma <- matrix(c(1, rho_target, rho_target, 1), 2, 2)
    
    for(r in 1:Rreps){
      z  <- mvrnorm(1, mu = c(0, 0), Sigma = Sigma)
      pAs[r]    <- pnorm(z[1])         # Uniform(0,1) marginals
      pBs[r]    <- pnorm(z[2])
      pAKSAs[r] <- pBs[r]              # dummy (unused for null power)
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
    key <- paste(as.character(rho_target), collapse = "_")
    result.p[[key]] <- result
  }
})
grid$pAKSA <- sapply(result.p, function(df) mean(df$pAKSA <= 0.05, na.rm = T))
grid$pA <- sapply(result.p, function(df) mean(df$pA <= 0.05, na.rm = T))
grid$pB <- sapply(result.p, function(df) mean(df$pB <= 0.05, na.rm = T))
grid$pF <- sapply(result.p, function(df) mean(df$pF <= 0.05, na.rm = T))
grid$pCF <- sapply(result.p, function(df) mean(df$pCF <= 0.05, na.rm = T))

# find new cutoff to control type I error
crit_F <- sapply(result.p, function(df) quantile(df$pF, probs = 0.05, na.rm = TRUE))
grid$critF <- crit_F
write.csv(grid,'result/cor_null.csv', row.names = F)

gg <- read.csv('result/cor_null.csv')

df_long <- gg %>%     
  dplyr::select(rho, Fisher = pF, Brown = pCF) %>% 
  tidyr::pivot_longer(Fisher:Brown, names_to = "Method", values_to = "err")

df_long <- df_long %>% 
  mutate(se  = sqrt(err * (1 - err) / Rreps),
         lo  = pmax(0, err - 1.96 * se),
         hi  = pmin(1, err + 1.96 * se))

df_long$Method <- factor(df_long$Method,
                         levels = c("Fisher", "Brown"))

ggplot(df_long, aes(x = rho, y = err)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25, fill = "#7baac7") +
  geom_line(size = 1.2, colour = "#1f78b4") +
  geom_point(size = 2.5, colour = "#0d4f7d") +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey30") +
  facet_wrap(~ Method, nrow = 1) +
  labs(x = '\nCorrelation',
       y = "Type I error\n") +
  scale_y_continuous(limits = c(0.01,0.1), breaks = seq(0.01, 0.1, by = 0.02)) + 
  theme_minimal(base_size = 16) +
  theme(strip.text = element_text(face = "bold", size = 16),
        panel.spacing.x = unit(1.4, "lines"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"))
