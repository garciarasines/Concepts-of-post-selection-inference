library(patchwork)
source(file.path("Figures", "theme.R"))

set.seed(123)

alpha <- 0.10
c_sel <- 2
pi_signal <- 0.10
prior_var <- 3
prior_sd <- sqrt(prior_var)
y_var_signal <- 1 + prior_var

theta_grid <- seq(-8, 8, by = 0.0025)
dtheta <- theta_grid[2] - theta_grid[1]

y_abs_grid <- seq(c_sel + 0.001, 5, by = 0.01)

marginal_cdf <- function(y) {
  (1 - pi_signal) * pnorm(y) +
    pi_signal * pnorm(y, sd = sqrt(y_var_signal))
}

marginal_pdf <- function(y) {
  (1 - pi_signal) * dnorm(y) +
    pi_signal * dnorm(y, sd = sqrt(y_var_signal))
}

p_selection <- marginal_cdf(-c_sel) + 1 - marginal_cdf(c_sel)

marginal_cdf_selected <- function(y) {
  out <- numeric(length(y))
  
  left <- y < -c_sel
  middle <- y >= -c_sel & y < c_sel
  right <- y >= c_sel
  
  out[left] <- marginal_cdf(y[left]) / p_selection
  out[middle] <- marginal_cdf(-c_sel) / p_selection
  out[right] <- (
    marginal_cdf(-c_sel) +
      marginal_cdf(y[right]) -
      marginal_cdf(c_sel)
  ) / p_selection
  
  out
}

F_theta_selected <- function(y, theta) {
  p_left <- pnorm(-c_sel, mean = theta, sd = 1)
  p_middle_end <- pnorm(c_sel, mean = theta, sd = 1)
  p_E <- p_left + 1 - p_middle_end
  
  ifelse(
    y < -c_sel,
    pnorm(y, mean = theta, sd = 1) / p_E,
    ifelse(
      y < c_sel,
      p_left / p_E,
      (p_left + pnorm(y, mean = theta, sd = 1) - p_middle_end) / p_E
    )
  )
}

q_theta_selected <- function(u, theta) {
  u <- pmin(pmax(u, 1e-12), 1 - 1e-12)
  
  p_left <- pnorm(-c_sel, mean = theta, sd = 1)
  p_middle_end <- pnorm(c_sel, mean = theta, sd = 1)
  p_E <- p_left + 1 - p_middle_end
  
  left_mass <- p_left / p_E
  
  ifelse(
    u <= left_mass,
    qnorm(u * p_E, mean = theta, sd = 1),
    qnorm(p_middle_end + u * p_E - p_left, mean = theta, sd = 1)
  )
}

H <- function(w, theta) {
  y_lower <- q_theta_selected(alpha * w, theta)
  y_upper <- q_theta_selected(alpha * w + 1 - alpha, theta)
  
  marginal_cdf_selected(y_upper) - marginal_cdf_selected(y_lower)
}

w_opt <- sapply(theta_grid, function(theta) {
  optimize(
    f = function(w) H(w, theta),
    interval = c(1e-6, 1 - 1e-6)
  )$minimum
})

acceptance_bounds <- function(w_vec) {
  data.frame(
    theta = theta_grid,
    lower = mapply(
      function(w, theta) q_theta_selected(alpha * w, theta),
      w_vec,
      theta_grid
    ),
    upper = mapply(
      function(w, theta) q_theta_selected(alpha * w + 1 - alpha, theta),
      w_vec,
      theta_grid
    )
  )
}

bounds_Optimal <- acceptance_bounds(w_opt)
bounds_Equal_tailed <- acceptance_bounds(rep(0.5, length(theta_grid)))

confidence_set_size <- function(y, bounds) {
  inside <- bounds$lower <= y & y <= bounds$upper
  inside[is.na(inside)] <- FALSE
  
  if (!any(inside)) {
    return(0)
  }
  
  runs <- rle(inside)
  ends <- cumsum(runs$lengths)
  starts <- ends - runs$lengths + 1
  
  true_runs <- which(runs$values)
  total <- 0
  
  for (k in true_runs) {
    i1 <- starts[k]
    i2 <- ends[k]
    
    theta_left <- bounds$theta[i1]
    theta_right <- bounds$theta[i2]
    
    total <- total + theta_right - theta_left + dtheta
  }
  
  total
}

size_df_wide <- data.frame(
  abs_y = y_abs_grid,
  Optimal = sapply(y_abs_grid, confidence_set_size, bounds = bounds_Optimal),
  Equal_tailed = sapply(y_abs_grid, confidence_set_size, bounds = bounds_Equal_tailed)
)

size_df <- rbind(
  data.frame(
    abs_y = size_df_wide$abs_y,
    size = size_df_wide$Optimal,
    method = "Optimal"
  ),
  data.frame(
    abs_y = size_df_wide$abs_y,
    size = size_df_wide$Equal_tailed,
    method = "Equal_tailed"
  )
)

size_df$method <- factor(size_df$method, levels = c("Optimal", "Equal_tailed"))

density_abs_selected <- function(x) {
  ifelse(
    x >= c_sel,
    2 * marginal_pdf(x) / p_selection,
    0
  )
}

cdf_abs_selected <- function(x) {
  sapply(x, function(z) {
    if (z < c_sel) {
      return(0)
    }
    
    prob <- (marginal_cdf(-c_sel) - marginal_cdf(-z)) +
      (marginal_cdf(z) - marginal_cdf(c_sel))
    
    prob / p_selection
  })
}

density_df <- data.frame(
  abs_y = y_abs_grid,
  density = density_abs_selected(y_abs_grid),
  cdf = cdf_abs_selected(y_abs_grid)
)

N <- 500000

is_signal <- rbinom(N, size = 1, prob = pi_signal)

theta <- ifelse(
  is_signal == 1,
  rnorm(N, mean = 0, sd = prior_sd),
  0
)

y <- rnorm(N, mean = theta, sd = 1)

hist_df <- data.frame(
  abs_y = abs(y[abs(y) > c_sel])
)

hist_df <- hist_df[hist_df$abs_y >= 2 & hist_df$abs_y <= 5, , drop = FALSE]

dy <- y_abs_grid[2] - y_abs_grid[1]

avg_Optimal <- sum(size_df_wide$Optimal * density_df$density) * dy
avg_Equal_tailed <- sum(size_df_wide$Equal_tailed * density_df$density) * dy

cross_x <- 3.6

p_plot_a <- ggplot(size_df, aes(x = abs_y, y = size, color = method)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = cross_x, linetype = "dashed", linewidth = 0.4) +
  scale_color_manual(
    values = c(
      "Optimal" = "black",
      "Equal_tailed" = "grey65"
    )
  ) +
  coord_cartesian(xlim = c(2, 5), ylim = c(1.1, 4.6)) +
  labs(
    x = NULL,
    y = "Size of confidence set",
    color = "Method"
  ) +
  theme_book +
  theme(
    legend.position = "top",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 11),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

p_plot_b <- ggplot() +
  geom_histogram(
    data = hist_df,
    aes(x = abs_y, y = after_stat(density)),
    breaks = seq(2, 5, by = 0.12),
    fill = "white",
    color = "grey45",
    linewidth = 0.25,
    closed = "left"
  ) +
  geom_line(
    data = density_df,
    aes(x = abs_y, y = density),
    linewidth = 0.9,
    color = "black"
  ) +
  geom_line(
    data = density_df,
    aes(x = abs_y, y = cdf),
    linewidth = 0.9,
    color = "black"
  ) +
  geom_vline(xintercept = cross_x, linetype = "dashed", linewidth = 0.4) +
  coord_cartesian(xlim = c(2, 5), ylim = c(0, 1.55)) +
  labs(
    x = expression(abs(y)),
    y = "Density / CDF"
  ) +
  theme_book +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

ggsave(file.path("Figures", "Outputs", "fig-6-12-a.pdf"), plot = p_plot_a, width = 4, height = 3)
ggsave(file.path("Figures", "Outputs", "fig-6-12-b.pdf"), plot = p_plot_b, width = 4, height = 3)

p_plot <- p_plot_a / p_plot_b +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  )

ggsave(file.path("Figures", "Outputs", "fig-6-12.pdf"), plot = p_plot, width = 4, height = 6)
