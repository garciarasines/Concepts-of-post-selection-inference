source(file.path("Figures", "theme.R"))

Kval <- function(delta, k) {
  qnorm((1 + (1 - delta)^(1/k))/2)
}

B <- 5e4
nummeans <- 500
numsel <- 5
alpha <- 0.1
alphaadj <- alpha

set.seed(123)

K1 <- Kval(alpha, nummeans)

es <- matrix(rnorm(nummeans*B), nrow = nummeans, ncol = B)
ws <- matrix(rnorm(nummeans*B), nrow = nummeans, ncol = B)

quant <- qnorm(0.95)
quantadj <- qnorm(1 - alphaadj)

n1 <- 0.9*nummeans
mu <- rnorm(nummeans, sd = 1)
mu[1:n1] <- 0

coveragesUV <- function(gamma) {
  cov <- matrix(0, B, numsel)
  covunadj <- matrix(0, B, numsel)
  covCOS <- matrix(0, B, numsel)
  covUV <- matrix(0, B, numsel)
  len <- matrix(0, B, numsel)
  lenunadj <- matrix(0, B, numsel)
  lenCOS <- matrix(0, B, numsel)
  lenUV <- matrix(0, B, numsel)
  selloss <- matrix(0, B, numsel)
  
  for (b in seq_len(B)) {
    y <- mu + es[, b]
    u <- y + gamma*ws[, b]
    
    rank_mu <- rank(-mu)
    rank_u <- rank(-u)
    
    thresh <- u[rank_u == (numsel + 1)]
    
    for (l in seq_len(numsel)) {
      J <- which(rank_mu == l)
      I <- which(rank_u == l)
      
      selloss[b, l] <- (mu[J] - mu[I])^2
      
      yI <- y[I]
      
      cov[b, l] <- abs(yI - mu[I]) < K1
      len[b, l] <- 2*K1
      
      covunadj[b, l] <- abs(yI - mu[I]) < quant
      lenunadj[b, l] <- 2*quant
      
      yobs <- y[I]
      
      CDF_yobs <- function(theta) {
        pdf_fun <- function(z) {
          exp(dnorm(z, mean = theta, sd = 1, log = TRUE) + pnorm((z - thresh)/gamma, mean = 0, sd = 1, log.p = TRUE) - pnorm((theta - thresh)/sqrt(1 + gamma^2), mean = 0, sd = 1, log.p = TRUE))
        }
        integrate(pdf_fun, lower = -25, upper = yobs)$value
      }
      
      lower <- tryCatch(uniroot(function(x) CDF_yobs(x) - (1 - alphaadj/2), c(yobs - qnorm(0.95) - 18, yobs - qnorm(0.95) + 18))$root, error = function(e) -Inf)
      upper <- tryCatch(uniroot(function(x) CDF_yobs(x) - alphaadj/2, c(yobs + qnorm(0.95) - 18, yobs + qnorm(0.95) + 18))$root, error = function(e) Inf)
      
      covCOS[b, l] <- (mu[I] < upper)*(mu[I] > lower)
      lenCOS[b, l] <- upper - lower
      
      cent <- y[I] - ws[I, b]/gamma
      
      covUV[b, l] <- (cent - sqrt(1 + 1/gamma^2)*qnorm(1 - alphaadj/2) < mu[I])*(cent + sqrt(1 + 1/gamma^2)*qnorm(1 - alphaadj/2) > mu[I])
      lenUV[b, l] <- 2*sqrt(1 + 1/gamma^2)*qnorm(1 - alphaadj/2)
    }
  }
  
  c(gamma, mean(selloss), mean(lenUV), mean(covUV), median(lenCOS), mean(covCOS), mean(len), mean(cov), mean(lenunadj), mean(covunadj))
}

gammas <- seq(0.1, 4, length.out = 60)

pb <- txtProgressBar(min = 0, max = length(gammas), style = 3)

res <- sapply(seq_along(gammas), function(j) {
  out <- coveragesUV(gammas[j])
  setTxtProgressBar(pb, j)
  out
})

close(pb)

df_length <- data.frame(
  x = rep(res[1, ], 4),
  y = c(res[3, ], res[5, ], res[7, ], res[9, ]),
  method = factor(
    rep(c("V", "COS", "PoSI", "Unadjusted"), each = ncol(res)),
    levels = c("V", "COS", "PoSI", "Unadjusted")
  )
)

left_min <- 0
left_max <- 35
right_min <- min(res[2, ], na.rm = TRUE)
right_max <- max(res[2, ], na.rm = TRUE)

to_left <- function(z) {
  left_min + (z - right_min)/(right_max - right_min)*(left_max - left_min)
}

to_right <- function(z) {
  right_min + (z - left_min)/(left_max - left_min)*(right_max - right_min)
}

df_loss <- data.frame(
  x = res[1, ],
  y = to_left(res[2, ])
)

p_plot <- ggplot() +
  geom_line(
    data = df_length,
    aes(x = x, y = y, color = method, linetype = method),
    linewidth = 0.8
  ) +
  geom_line(
    data = df_loss,
    aes(x = x, y = y),
    linewidth = 0.8,
    linetype = "42",
    color = "grey40"
  ) +
  scale_color_manual(values = c(
    "V" = "black",
    "COS" = "grey55",
    "PoSI" = "grey25",
    "Unadjusted" = "black"
  )) +
  scale_linetype_manual(values = c(
    "V" = "solid",
    "COS" = "solid",
    "PoSI" = "22",
    "Unadjusted" = "11"
  )) +
  scale_x_continuous(breaks = seq(0, 4, by = 0.4)) +
  scale_y_continuous(
    name = "Average interval length",
    limits = c(left_min, left_max),
    breaks = seq(0, 35, by = 5),
    sec.axis = sec_axis(
      trans = ~ to_right(.),
      name = "Avg. sel. loss",
      breaks = pretty(res[2, ])
    )
  ) +
  coord_cartesian(xlim = c(0, 4)) +
  labs(
    x = expression(gamma),
    title = "",
    color = NULL,
    linetype = NULL
  ) +
  theme_book +
  theme(
    legend.position = c(0.68, 0.5),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.key.width = grid::unit(1, "cm"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.title.y.right = element_text(color = "grey40"),
    axis.text.y.right = element_text(color = "grey40")
  )

ggsave(file.path("Figures", "Outputs", "fig-5-07.pdf"), plot = p_plot, width = 4, height = 3)
