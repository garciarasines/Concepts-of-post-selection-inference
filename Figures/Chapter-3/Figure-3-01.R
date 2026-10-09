source(file.path("Figures", "theme.R"))

Kval <- function(delta, k) {
  qnorm((1 + (1 - delta)^(1 / k)) / 2)
}

B <- 2e5
nummeans <- 10
alpha <- 0.1
nu <- 0.05 * alpha

set.seed(123)

K1 <- Kval(alpha, nummeans)

const <- numeric(nummeans)
for (j in seq_len(nummeans)) {
  const[j] <- Kval(alpha - nu, j)
}

# Common random numbers across values of spread
es <- matrix(
  rnorm(nummeans * B),
  nrow = nummeans,
  ncol = B
)

mus <- matrix(
  rnorm(nummeans * B),
  nrow = nummeans,
  ncol = B
)

K2 <- Kval(nu, nummeans)

coverages <- function(spread) {
  I <- numeric(B)
  J <- numeric(B)
  agree <- numeric(B)
  
  cov <- numeric(B)
  covunadj <- numeric(B)
  covCOS <- numeric(B)
  covls <- numeric(B)
  
  len <- numeric(B)
  lenunadj <- numeric(B)
  lenCOS <- numeric(B)
  lenls <- numeric(B)
  
  quant <- qnorm(0.95)
  
  for (b in seq_len(B)) {
    
    # Same latent mean vectors are used for every spread
    mu <- spread * mus[, b]
    
    J[b] <- which(mu == max(mu))
    
    y <- mu + es[, b]
    I[b] <- which(y == max(y))
    yI <- y[I[b]]
    
    # Fully simultaneous interval
    cov[b] <- abs(yI - mu[I[b]]) < K1
    len[b] <- 2 * K1
    
    # Unadjusted interval
    covunadj[b] <- abs(yI - mu[I[b]]) < quant
    lenunadj[b] <- 2 * quant
    
    # Conditional selective interval
    thresh <- max(y[-I[b]])
    yobs <- y[I[b]]
    
    CDF_yobs <- function(theta) {
      pdf_fun <- function(z) {
        exp(
          dnorm(z, mean = theta, sd = 1, log = TRUE) -
            pnorm(
              theta - thresh,
              mean = 0,
              sd = 1,
              log.p = TRUE
            )
        )
      }
      
      integrate(
        pdf_fun,
        lower = thresh,
        upper = yobs
      )$value
    }
    
    # Locally simultaneous interval
    index <- sum(y > y[I[b]] - 4 * K2)
    
    temp <- min(
      const[index],
      K1
    )
    
    lowerls <- y[I[b]] - temp
    upperls <- y[I[b]] + temp
    
    covls[b] <- (
      mu[I[b]] < upperls
    ) * (
      mu[I[b]] > lowerls
    )
    
    lenls[b] <- upperls - lowerls
    
    lower <- tryCatch(
      uniroot(
        function(x) CDF_yobs(x) - 0.95,
        c(
          yobs - qnorm(0.95) - 400,
          yobs - qnorm(0.95) + 10
        )
      )$root,
      error = function(e) -Inf
    )
    
    upper <- tryCatch(
      uniroot(
        function(x) CDF_yobs(x) - 0.05,
        c(
          yobs + qnorm(0.95) - 400,
          yobs + qnorm(0.95) + 10
        )
      )$root,
      error = function(e) Inf
    )
    
    covCOS[b] <- (
      mu[I[b]] < upper
    ) * (
      mu[I[b]] > lower
    )
    
    lenCOS[b] <- upper - lower
    
    agree[b] <- as.numeric(I[b] == J[b])
  }
  
  c(
    spread,
    mean(agree),
    median(len),
    median(lenls),
    median(lenCOS),
    mean(cov),
    mean(covls),
    mean(covCOS),
    median(lenunadj),
    mean(covunadj)
  )
}

spreads <- seq(
  1,
  20,
  length.out = 77
)

pb <- txtProgressBar(
  min = 0,
  max = length(spreads),
  style = 3
)

res <- sapply(
  seq_along(spreads),
  function(j) {
    out <- coverages(spreads[j])
    setTxtProgressBar(pb, j)
    out
  }
)

close(pb)

df_length <- data.frame(
  x = rep(res[1, ], 4),
  y = c(
    res[3, ],
    res[4, ],
    res[5, ],
    res[9, ]
  ),
  method = factor(
    rep(
      c(
        "Simultaneous",
        "Local simultaneous",
        "COS",
        "Unadjusted"
      ),
      each = ncol(res)
    ),
    levels = c(
      "Simultaneous",
      "Local simultaneous",
      "COS",
      "Unadjusted"
    )
  )
)

left_min <- 2
left_max <- 8.3
right_min <- 0.75
right_max <- 0.90

to_left <- function(z) {
  left_min +
    (z - right_min) /
    (right_max - right_min) *
    (left_max - left_min)
}

to_right <- function(z) {
  right_min +
    (z - left_min) /
    (left_max - left_min) *
    (right_max - right_min)
}

# Empirical coverage of the unadjusted interval
df_coverage <- data.frame(
  x = res[1, ],
  y = to_left(res[10, ])
)

p_plot <- ggplot() +
  geom_line(
    data = df_length,
    aes(
      x = x,
      y = y,
      color = method,
      linetype = method
    ),
    linewidth = 0.8
  ) +
  geom_line(
    data = df_coverage,
    aes(
      x = x,
      y = y
    ),
    linewidth = 0.8,
    linetype = "11",
    color = "grey40"
  ) +
  scale_color_manual(
    values = c(
      "Simultaneous" = "black",
      "Local simultaneous" = "grey55",
      "COS" = "grey25",
      "Unadjusted" = "black"
    )
  ) +
  scale_linetype_manual(
    values = c(
      "Simultaneous" = "solid",
      "Local simultaneous" = "solid",
      "COS" = "dashed",
      "Unadjusted" = "11"
    )
  ) +
  scale_x_continuous(
    breaks = seq(
      0,
      20,
      by = 4
    )
  ) +
  scale_y_continuous(
    name = "Median interval length",
    limits = c(
      left_min,
      left_max
    ),
    breaks = seq(
      2,
      8,
      by = 1.5
    ),
    sec.axis = sec_axis(
      trans = ~ to_right(.),
      name = "Unadj. cov.",
      breaks = c(
        0.75,
        0.80,
        0.85,
        0.90
      )
    )
  ) +
  coord_cartesian(
    xlim = c(
      0,
      20
    )
  ) +
  labs(
    x = "Std. dev. of theta distribution",
    title = "",
    color = NULL,
    linetype = NULL
  ) +
  theme_book +
  theme(
    legend.position = c(0.7, 0.7),
    legend.background = element_blank(),
    legend.key = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.title.y.right = element_text(
      color = "grey40"
    ),
    axis.text.y.right = element_text(
      color = "grey40"
    )
  )

ggsave(file.path("Figures", "Outputs", "fig-3-01.pdf"), plot = p_plot, width = 4, height = 3)
