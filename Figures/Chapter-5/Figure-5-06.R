library(mvtnorm)
library(ggplot2)
source(file.path("Figures", "theme.R"))

set.seed(123)

alpha <- 0.10
q1 <- alpha/2
q2 <- 1 - alpha/2
B <- 2e6
n_y_grid <- 151
mu <- 0

df_t <- 5
scale_t <- sqrt((df_t - 2)/df_t)

Sigma_YU <- matrix(c(1, 1, 1, 2), nrow = 2)

z_alpha <- qnorm(1 - alpha/2)

y_all <- mu + scale_t*rt(B, df = df_t)
w_all <- rnorm(B)
u_all <- y_all + w_all
v_all <- y_all - w_all

p_selection_normal <- function(theta, t) {
  1 - pnorm(t, mean = theta, sd = sqrt(2))
}

prob_y_leq_and_selection_normal <- function(y, theta, t) {
  as.numeric(
    pmvnorm(
      lower = c(-Inf, t),
      upper = c(y, Inf),
      mean = c(theta, theta),
      sigma = Sigma_YU
    )
  )
}

F_carve <- function(y, theta, t) {
  p_sel <- p_selection_normal(theta, t)
  
  if (!is.finite(p_sel) || p_sel <= 0) {
    return(NA_real_)
  }
  
  prob_y_leq_and_selection_normal(y, theta, t)/p_sel
}

root_for_q <- function(y, t, q) {
  f <- function(theta) {
    val <- F_carve(y, theta, t) - q
    
    if (!is.finite(val)) {
      return(NA_real_)
    }
    
    val
  }
  
  lower <- y - 8
  upper <- y + 8
  
  f_lower <- f(lower)
  f_upper <- f(upper)
  
  count <- 0
  
  while ((!is.finite(f_lower) || f_lower < 0) && count < 50) {
    lower <- lower - 2
    f_lower <- f(lower)
    count <- count + 1
  }
  
  count <- 0
  
  while ((!is.finite(f_upper) || f_upper > 0) && count < 50) {
    upper <- upper + 2
    f_upper <- f(upper)
    count <- count + 1
  }
  
  if (!is.finite(f_lower) || !is.finite(f_upper)) {
    return(NA_real_)
  }
  
  if (f_lower*f_upper > 0) {
    return(NA_real_)
  }
  
  uniroot(f, interval = c(lower, upper))$root
}

carving_interval <- function(y, t) {
  lower <- root_for_q(y, t, q2)
  upper <- root_for_q(y, t, q1)
  
  c(lower, upper)
}

compute_coverage <- function(t_grid) {
  out <- data.frame(
    t = t_grid,
    carving = NA_real_,
    fission = NA_real_,
    selection = NA_real_
  )
  
  pb <- txtProgressBar(min = 0, max = length(t_grid), style = 3)
  
  for (i in seq_along(t_grid)) {
    t <- t_grid[i]
    
    selected <- u_all >= t
    
    y_sel <- y_all[selected]
    v_sel <- v_all[selected]
    
    out$selection[i] <- mean(selected)
    
    y_grid <- seq(
      quantile(y_sel, 0.001),
      quantile(y_sel, 0.999),
      length.out = n_y_grid
    )
    
    ci_grid <- t(sapply(y_grid, carving_interval, t = t))
    keep <- is.finite(ci_grid[, 1]) & is.finite(ci_grid[, 2])
    
    if (sum(keep) >= 2) {
      lower_fun <- approxfun(y_grid[keep], ci_grid[keep, 1], rule = 2)
      upper_fun <- approxfun(y_grid[keep], ci_grid[keep, 2], rule = 2)
      
      lower_carve <- lower_fun(y_sel)
      upper_carve <- upper_fun(y_sel)
      
      out$carving[i] <- mean(lower_carve <= mu & mu <= upper_carve, na.rm = TRUE)
    }
    
    lower_fission <- v_sel - z_alpha*sqrt(2)
    upper_fission <- v_sel + z_alpha*sqrt(2)
    
    out$fission[i] <- mean(lower_fission <= mu & mu <= upper_fission)
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  out
}

t_grid <- seq(-4, 5, length.out = 75)

res <- compute_coverage(t_grid)

df_long <- rbind(
  data.frame(t = res$t, value = res$carving, method = "Carving"),
  data.frame(t = res$t, value = res$fission, method = "Fission"),
  data.frame(t = res$t, value = res$selection, method = "Selection prob.")
)

p_plot <- ggplot(df_long, aes(x = t, y = value, color = method, linetype = method)) +
  geom_hline(yintercept = 0.9, linewidth = 0.4, linetype = "dashed", color = "grey55") +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  scale_color_manual(values = c("Carving" = "black", "Fission" = "grey45", "Selection prob." = "grey70"), guide = "none") +
  scale_linetype_manual(values = c("Carving" = "solid", "Fission" = "solid", "Selection prob." = "22")) +
  guides(linetype = guide_legend(override.aes = list(color = c("black", "grey45", "grey70"), linewidth = 1.2))) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    name = "Coverage",
    sec.axis = sec_axis(~., name = "Selection prob.", breaks = seq(0, 1, by = 0.2))
  ) +
  coord_cartesian(xlim = c(-4, 5)) +
  labs(x = "t", color = NULL, linetype = NULL) +
  theme_book +
  theme(legend.position = "bottom", legend.key.width = grid::unit(1, "cm"))

ggsave(file.path("Figures", "Outputs", "fig-5-06.pdf"), plot = p_plot, width = 4, height = 3)
