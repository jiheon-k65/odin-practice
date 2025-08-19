# Session 11: simulating the effect of non-random mizing on rubella transmission & control
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
    lambda_mat[, ] <- beta[i, j] * I[i]
    lambda[] <- sum(lambda_mat[i, ])

    # To define derivatives
    # For young persons (index pos. 1)
    deriv(S[1]) <- ((1 - p) * mu) - lambda[i] * S[i] - eta[i] * S[i]
    deriv(E[1]) <- lambda[i] * S[i] - sigma * E[i] - eta[i] * E[i]
    deriv(I[1]) <- sigma * E[i] - gamma * I[i] - eta[i] * I[i]
    deriv(R[1]) <- (p * mu) + gamma * I[i] - eta[i] * R[i]
    # For middle-aged persons (index pos. 2)
    deriv(S[2:(n_grp - 1)]) <- eta[i - 1] * S[i - 1] - lambda[i] * S[i] - eta[i] * S[i]
    deriv(E[2:(n_grp - 1)]) <- eta[i - 1] * E[i - 1] + lambda[i] * S[i] - sigma * E[i] - eta[i] * E[i]
    deriv(I[2:(n_grp - 1)]) <- eta[i - 1] * I[i - 1] + sigma * E[i] - gamma * I[i] - eta[i] * I[i]
    deriv(R[2:(n_grp - 1)]) <- eta[i - 1] * R[i - 1] + gamma * I[i] - eta[i] * R[i]
    # For old persons (index pos. 3)
    deriv(S[n_grp]) <- eta[i - 1] * S[i - 1] - lambda[i] * S[i] - delta * S[i]
    deriv(E[n_grp]) <- eta[i - 1] * E[i - 1] + lambda[i] * S[i] - sigma * E[i] - delta * E[i]
    deriv(I[n_grp]) <- eta[i - 1] * I[i - 1] + sigma * E[i] - gamma * I[i] - delta * I[i]
    deriv(R[n_grp]) <- eta[i - 1] * R[i - 1] + gamma * I[i] - delta * R[i]
    # Dummy compartment to hold cumulative infections
    deriv(ci[1:n_grp]) <- lambda[i] * S[i]

    # To define model inputs
    S_init[] <- user()
    E_init[] <- user()
    I_init[] <- user()
    R_init[] <- user()
    ci_init[] <- user()
    n_grp <- user(3) # No. of age groups
    beta[, ] <- user() # Effective contact rate
    dur_lat <- user(10) # Avg duration of latency (days)
    sigma <- 1 / dur_lat # Avg rate of being infectious (day^-1)
    dur_inf <- user(11) # Avg duration of being infectious (days)
    gamma <- 1 / dur_inf # Avg rate of recovery (day^-1)
    dur_grp[] <- user() # Avg time spent in each age group (years)
    eta[] <- 1 / (dur_grp[i] * 365) # Avg rate of aging (day^-1)
    delta <- eta[n_grp] # Avg rate of dying (day^-1)
    pop_o <- 30000 # Total no. of old persons
    mu <- delta * pop_o # Avg no. of births at t
    vac_cov <- user(0) # Vaccine coverage
    vac_eff <- user(1) # Vaccine efficacy
    vac_str <- user(100) # Start of vaccination (year)

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
    dim(beta) <- c(n_grp, n_grp)
    dim(lambda_mat) <- c(n_grp, n_grp)
    dim(lambda) <- n_grp
    dim(dur_grp) <- n_grp
    dim(eta) <- n_grp
})

# ==== OUTSIDE ODIN ====
# To define initial conditions
n_grp <- 3 # Need to re-assign n_grp outside of odin
S0_y <- 5010 # Initial no. of young susceptibles
S0_m <- 3083
S0_o <- 2740
I0_y <- 20
I0_m <- 4
I0_o <- 10
R0_y <- 9970
R0_m <- 11913
R0_o <- 27250
ini_con_data <- c(
    S0_y, S0_m, S0_o,
    rep(0, n_grp), # E
    I0_y, I0_m, I0_o,
    R0_y, R0_m, R0_o,
    rep(0, n_grp) # ci
)
ini_con <- matrix(ini_con_data, nrow = n_grp, ncol = 5, byrow = FALSE)

# To define effective contact rates (b)
b_yy <- 1.81e-5 # Young infect young person
b_ym <- 0 # Middle-aged infects young person
b_yo <- 0 # Old infects young person
b_my <- 0
b_mm <- 2.95e-5
b_mo <- 0
b_oy <- 0
b_om <- 0
b_oo <- 3.32e-5
# To define matrix to hold effective contact rates
b_mat_data <- c(
    b_yy, b_ym, b_yo,
    b_my, b_mm, b_mo,
    b_oy, b_om, b_oo
)
b_mat <- matrix(b_mat_data, nrow = n_grp, ncol = n_grp, byrow = TRUE)

# To define time spent in each age group (years)
dur_y <- 15 # 0 to 14 yrs old
dur_m <- 15 # 15 to 29 yrs old
dur_o <- 30 # 30 to 59 yrs old
# To define vector to hold durations
dur_vec <- c(dur_y, dur_m, dur_o)

# To create an instance of the model
mod <- seir$new(
    S_init = ini_con[, 1],
    I_init = ini_con[, 2],
    E_init = ini_con[, 3],
    R_init = ini_con[, 4],
    ci_init = ini_con[, 5],
    beta = b_mat,
    dur_grp = dur_vec,
    vac_cov = 0.86
)

# To assess characteristics of model
print(mod)

# To simulate model over time
t <- seq(from = 0, to = 73000, by = 1)
y <- mod$run(t)

# To convert deSolve matrix to dataframe
y <- as.data.frame(y)

# ==== PROPORTION SUSCEPTIBLE ====
# To calculate total population for each age group
for (i in 1:n_grp) {
    cols <- paste0(c("S[", "E[", "I[", "R["), i, "]") # This is a string vector
    y[[paste0("N[", i, "]")]] <- rowSums(y[, cols])
}
# To calculate proportion of susceptibles
for (i in 1:n_grp) {
    S_col <- paste0("S[", i, "]") # This is a string
    N_col <- paste0("N[", i, "]") # This is a string
    y[[paste0("Prop_S[", i, "]")]] <- y[[S_col]] / y[[N_col]]
}
# To create a subset of columns
prop_sus <- y[, c("t", paste0("Prop_S[", 1:n_grp, "]"))]

# To label columns
idx_lab <- c("Young", "Middle-aged", "Old")
colnames(prop_sus) <- c("t", idx_lab)

# To convert dataframe to long format
prop_sus <- melt(as.data.frame(prop_sus), id = "t")

# To plot model output
sus_plot <-
    ggplot(data = prop_sus, aes(x = t, y = value, color = variable, group = variable)) +
    geom_line(linewidth = 0.5) +
    geom_vline(xintercept = 100 * 365, color = "black", linetype = "dashed", linewidth = 0.3) +
    annotate(geom = "text", label = "Vaccination", color = "black", x = 35000, y = 0.38, angle = 90, vjust = 0, size = 2) +
    coord_cartesian(ylim = c(0, 0.4)) +
    xlab("Time (days)") +
    ylab("Proportion of susceptible individuals") +
    labs(color = "Age group") +
    scale_color_brewer(palette = "Dark2") +
    theme_bw(base_size = 8) +
    theme(legend.position = "top")
ggsave("rub_mix_sus.png", plot = sus_plot, height = 9, width = 7, units = "cm", dpi = 300)

# ==== NEW INFECTIONS ====
# To calculate no. of new infections per day
for (i in 1:n_grp) {
    ci_col <- paste0("ci[", i, "]")
    y[[paste0("ni[", i, "]")]] <- c(0, diff(y[[ci_col]]))
}
# To scale to 10,000 persons (for visualization purposes)
for (i in 1:n_grp) {
    ni_col <- paste0("ni[", i, "]")
    N_col <- paste0("N[", i, "]")
    y[[paste0("sc_ni[", i, "]")]] <- y[[ni_col]] / y[[N_col]] * 100000
}

# To select specific columns
new_inf <- y[, c("t", paste0("sc_ni[", 1:n_grp, "]"))]

# To label columns
colnames(new_inf) <- c("t", idx_lab)

# To convert dataframe to long format
new_inf <- melt(as.data.frame(new_inf), id = "t")

# To plot model output
inf_plot <-
    ggplot(data = new_inf, aes(x = t, y = value, color = variable, group = variable)) +
    geom_line(linewidth = 0.5) +
    geom_vline(xintercept = 100 * 365, color = "black", linetype = "dashed", linewidth = 0.3) +
    annotate(geom = "text", label = "Vaccination", color = "black", x = 35000, y = 19, angle = 90, vjust = 0, size = 2) +
    coord_cartesian(ylim = c(0, 20)) +
    xlab("Time (days)") +
    ylab("No. of new infections (per 100,000)") +
    labs(color = "Age") +
    scale_color_brewer(palette = "Dark2") +
    theme_bw(base_size = 8) +
    theme(legend.position = "top")
ggsave("rub_mix_inf.png", plot = inf_plot, height = 9, width = 7, units = "cm", dpi = 300)

# To calculate time elapsed
end <- Sys.time()
time <- as.integer(difftime(end, start, units = "secs"))
print(paste("Time elapsed: ", time, "s"))
elapsed <- as.numeric(difftime(end, start, units = "secs"))
