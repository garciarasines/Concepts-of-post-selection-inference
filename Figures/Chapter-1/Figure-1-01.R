source(file.path("Figures", "theme.R"))

ns <- 1 + 5*(0:20)

df <- data.frame(
  n = ns,
  coverage = pnorm(1.96)^ns - pnorm(-1.96)^ns
)

p_plot <- ggplot(df, aes(x = n, y = coverage)) +
  geom_point(shape = 16, size = 2, color = "black") +
  geom_line(linewidth = 0.4, color = "black") +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "n", y = "Coverage") +
  theme_book

ggsave(file.path("Figures", "Outputs", "fig-1-01.pdf"), plot = p_plot, width = 4, height = 3)
