# Session 10: contrasing the effects of rubella vaccination
start <- Sys.time()

# To load packages
library(odin)
library(reshape2)
library(ggplot2)
library(RColorBrewer)

# ==== ODIN ====
# To define SEIR model
seir <- odin::odin({
    # To define initial states
    initial(S[]) <- S_init[i]
    initial(E[]) <- E_init[i]
    initial(I[]) <- I_init[i]
    initial(R[]) <- R_init[i]
    initial(ci[]) <- ci_init[i]

    # To vary vaccine start time (convert year to day)
    p <- if (t <= vac_str * 365) 0 else vac_cov * vac_eff

    # To define force of infection
    lambda <- beta * sum(I[])

    # To define derivatives
    # For 0 yrs old (index pos. 1)
    deriv(S[1]) <- ((1 - p) * mu) - lambda * S[i] - eta * S[i]
    deriv(E[1]) <- lambda * S[i] - sigma * E[i] - eta * E[i]
    deriv(I[1]) <- sigma * E[i] - gamma * I[i] - eta * I[i]
    deriv(R[1]) <- (p * mu) + gamma * I[i] - eta * R[i]
    # For 1 to 58 yrs old (index pos. 2 to 59)
    deriv(S[2:(n_grp - 1)]) <- eta * S[i - 1] - lambda * S[i] - eta * S[i]
    deriv(E[2:(n_grp - 1)]) <- eta * E[i - 1] + lambda * S[i] - sigma * E[i] - eta * E[i]
    deriv(I[2:(n_grp - 1)]) <- eta * I[i - 1] + sigma * E[i] - gamma * I[i] - eta * I[i]
    deriv(R[2:(n_grp - 1)]) <- eta * R[i - 1] + gamma * I[i] - eta * R[i]
    # For 59 yrs old (index pos. 60)
    deriv(S[n_grp]) <- eta * S[i - 1] - lambda * S[i]
    deriv(E[n_grp]) <- eta * E[i - 1] + lambda * S[i] - sigma * E[i]
    deriv(I[n_grp]) <- eta * I[i - 1] + sigma * E[i] - gamma * I[i]
    deriv(R[n_grp]) <- eta * R[i - 1] + gamma * I[i]
    # Dummy compartment to hold cumulative infections
    deriv(ci[]) <- lambda * S[i]

    # To define model inputs
    N <- user(60000) # Total population
    R0 <- user(12.19) # Basic reproduction no.
    dur_lat <- user(10) # Avg duration of latency (days)
    sigma <- 1 / dur_lat # Avg rate of being infectious (day^-1)
    dur_inf <- user(11) # Avg duration of being infectious (days)
    gamma <- 1 / dur_inf # Avg rate of recovery (day^-1)
    beta <- R0 / (N * dur_inf)
    eta <- 1 / 365 # Avg rate of aging (day^-1)
    b <- user(1000) # Avg no. of births per year
    mu <- b / 365 # Avg no. of births per day
    vac_cov <- user(0.4) # Vaccine coverage
    vac_eff <- user(1) # Vaccine efficacy
    vac_str <- user(100) # Start of vaccination (year)
    n_grp <- user(60) # No. of age groups
    S_init[] <- user()
    E_init[] <- user()
    I_init[] <- user()
    R_init[] <- user()
    ci_init[] <- user()

    # To define dimensions
    dim(S) <- n_grp
    dim(E) <- n_grp
    dim(I) <- n_grp
    dim(R) <- n_grp
    dim(ci) <- n_grp
    dim(S_init) <- n_grp
    dim(E_init) <- n_grp
    dim(I_init) <- n_grp
    dim(R_init) <- n_grp
    dim(ci_init) <- n_grp
})

# ==== OUTSIDE ODIN ====
# To define initial conditions
n_grp <- 60 # Need to re-assign n_grp outside of odin
data <- c(
    1000, rep(999, (n_grp - 1)), # S[1] = 1000, S[2:60] = 999
    rep(0, n_grp),
    0, rep(1, (n_grp - 1)), # I[1] = 0, I[2:60] = 1
    rep(0, n_grp),
    rep(0, n_grp)
)
ini_con <- matrix(data, nrow = n_grp, ncol = 5, byrow = FALSE)

# To create an instance of the model
mod <- seir$new(
    S_init = ini_con[, 1],
    I_init = ini_con[, 2],
    E_init = ini_con[, 3],
    R_init = ini_con[, 4],
    ci_init = ini_con[, 5]
)

# To assess characteristics of model
print(mod)

# To simulate model over time
t <- seq(from = 0, to = 109500, by = 1)
y <- mod$run(t)

# To convert deSolve matrix to dataframe
y <- as.data.frame(y)

# ==== SUSCEPTIBLES ====
# To create subset of dataframe
## 5, 20, 30 & 40 yrs old (index pos. 6, 21, 31 & 41)
idx <- c(6, 21, 31, 41)
sus <- y[, c("t", paste0("S[", idx, "]"))]

# To label columns
idx_lab <- c("5", "20", "30", "40")
colnames(sus) <- c("t", idx_lab)

# To convert dataframe to long format
sus <- melt(as.data.frame(sus), id = "t")

# To plot model output
plot <-
    ggplot(data = sus, aes(x = t, y = value, color = variable, group = variable)) +
    geom_line(linewidth = 0.5) +
    geom_vline(xintercept = 100 * 365, color = "black", linetype = "dashed", linewidth = 0.3) +
    annotate(geom = "text", label = "Vaccination", color = "black", x = 35000, y = 950, angle = 90, vjust = 0, size = 2) +
    xlab("Time (days)") +
    ylab("No. of susceptible individuals (per 1,000)") +
    labs(color = "Age") +
    scale_color_brewer(palette = "Dark2") +
    theme_bw(base_size = 8) +
    theme(legend.position = "top")
ggsave("rub_sus.png", plot = plot, height = 10, width = 8, units = "cm", dpi = 300)

# ==== NEW INFECTIONS ====
# To create subset of dataframe
inf <- y[, c("t", paste0("ci[", idx, "]"))]

# To calculate no. of new infections per day
for (i in idx) {
    inf[[paste0("ni[", i, "]")]] <- c(0, diff(inf[[paste0("ci[", i, "]")]]))
}

# To select specific columns
new_inf <- inf[, c("t", paste0("ni[", idx, "]"))]

# To scale by 100 (for visualization purposes)
new_inf[, paste0("ni[", idx, "]")] <- new_inf[, paste0("ni[", idx, "]")] * 100

# To label columns
colnames(new_inf) <- c("t", idx_lab)

# To convert dataframe to long format
new_inf <- melt(as.data.frame(new_inf), id = "t")

# To plot model output
plot <-
    ggplot(data = new_inf, aes(x = t, y = value, color = variable, group = variable)) +
    geom_line(linewidth = 0.5) +
    geom_vline(xintercept = 100 * 365, color = "black", linetype = "dashed", linewidth = 0.3) +
    annotate(geom = "text", label = "Vaccination", color = "black", x = 35000, y = 28.5, angle = 90, vjust = 0, size = 2) +
    coord_cartesian(xlim = c(0, 60000), ylim = c(0, 30)) +
    xlab("Time (days)") +
    ylab("No. of new infections (per 100,000)") +
    labs(color = "Age") +
    scale_color_brewer(palette = "Dark2") +
    theme_bw(base_size = 8) +
    theme(legend.position = "top")
ggsave("rub_inf.png", plot = plot, height = 10, width = 8, units = "cm", dpi = 300)

# ==== DAILY FORCE OF INFECTION ====
# To re-assign parameters outside odin
N <- 60000 # Total population
R0 <- 12.19 # Basic reproduction no.
dur_inf <- 11 # Avg duration of being infectious (days)
beta <- R0 / (N * dur_inf)

# To create subset of dataframe
tot_inf <- rowSums(y[, paste0("I[", 1:60, "]")]) # This is a vector
foi <- beta * tot_inf # This is a vector
foi <- data.frame(t = y$t, foi = foi) # This is a dataframe

# To plot model output
plot <-
    ggplot(data = foi, aes(x = t, y = foi)) +
    geom_line(color = "steelblue", linewidth = 0.5) +
    geom_vline(xintercept = 100 * 365, color = "black", linetype = "dashed", linewidth = 0.3) +
    annotate(geom = "text", label = "Vaccination", color = "black", x = 35000, y = 0.0038, angle = 90, vjust = 0, size = 2) +
    coord_cartesian(xlim = c(0, 60000), ylim = c(0, 0.004)) +
    xlab("Time (days)") +
    ylab("Daily force of infection") +
    theme_bw(base_size = 8) +
    theme(legend.position = "top")
ggsave("rub_foi.png", plot = plot, height = 10, width = 8, units = "cm", dpi = 300)

# To calculate time elapsed
end <- Sys.time()
time <- as.integer(difftime(end, start, units = "secs"))
print(paste("Time elapsed: ", time, "s"))
elapsed <- as.numeric(difftime(end, start, units = "secs"))
