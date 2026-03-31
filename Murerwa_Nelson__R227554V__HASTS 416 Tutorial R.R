# ============================================================
# HASTS416 Tutorial 1 — Markov Chains (FINAL CORRECTED VERSION)
# ============================================================

library(markovchain)
library(diagram)

# ============================================================
# HELPER FUNCTIONS
# ============================================================

simulate_chain <- function(P, n_steps, start = NULL) {
  n_states <- nrow(P)
  if (is.null(start)) start <- sample(1:n_states, 1)
  path <- numeric(n_steps + 1)
  path[1] <- start
  for (t in 1:n_steps) {
    path[t + 1] <- sample(1:n_states, 1, prob = P[path[t], ])
  }
  path
}

mat_power <- function(M, n) {
  result <- diag(nrow(M))
  for (i in seq_len(n)) result <- result %*% M
  result
}

plot_trajectories <- function(trajs, states, title_prefix = "Trajectory") {
  par(mfrow = c(length(trajs), 1), mar = c(3, 4, 2, 1))
  cols <- c("steelblue", "tomato", "forestgreen")
  for (i in seq_along(trajs)) {
    traj <- trajs[[i]]
    plot(0:(length(traj) - 1), traj,
         type = "s", col = cols[i], lwd = 2,
         xlab = "Step", ylab = "State",
         main = sprintf("%s %d (Start: %s)", title_prefix, i, states[traj[1]]),
         yaxt = "n", ylim = c(0.5, length(states) + 0.5))
    axis(2, at = 1:length(states), labels = states)
    abline(h = 1:length(states), lty = 3, col = "gray80")
  }
  par(mfrow = c(1, 1))
}

# ============================================================
# QUESTION A1 — 5-STATE MARKOV CHAIN
# ============================================================

states1 <- paste0("S", 1:5)

P1 <- matrix(c(
  1,   0, 0, 0, 0,
  0.5, 0, 0, 0, 0.5,
  0.2, 0, 0, 0, 0.8,
  0,   0, 1, 0, 0,
  0,   0, 0, 1, 0
), nrow = 5, byrow = TRUE)

rownames(P1) <- colnames(P1) <- states1
mc1 <- new("markovchain", states = states1, transitionMatrix = P1)

cat("Transition Matrix (A1):\n")
print(P1)

# ---- (a) Classification ----
plot(mc1, main = "5-State Markov Chain (A1)")

cat("\n--- A1(a): Classification ---\n")

cat("\nRecurrent Classes:\n")
print(recurrentClasses(mc1))

cat("\nTransient Classes:\n")
print(transientClasses(mc1))

cat("\nAbsorbing States:\n")
print(absorbingStates(mc1))

cat("\nIs Irreducible:", is.irreducible(mc1), "\n")

cat("\nChain Period (based on recurrent class):\n")
print(period(mc1))  # Only meaningful for recurrent class {S1}

# ---- (b) Simulation ----
set.seed(42)
trajs1 <- list(
  simulate_chain(P1, 30),
  simulate_chain(P1, 30),
  simulate_chain(P1, 30)
)
plot_trajectories(trajs1, states1)

cat("\nA1(b) Trajectory Comment:\n")
cat("All trajectories eventually get absorbed into S1.\n")

# ---- (c) Steady state ----
cat("\n--- A1(c): Steady State ---\n")
ss1 <- steadyStates(mc1)
print(ss1)

cat("\nInterpretation:\n")
cat("Chain converges to S1 with probability 1.\n")
cat("Not ergodic (not irreducible).\n")

# ---- (d) Convergence ----
pi0_1 <- rep(1/5, 5)
n_time <- 50
prob1 <- matrix(0, nrow = n_time + 1, ncol = 5)
prob1[1, ] <- pi0_1

Pn1 <- diag(5)
for (t in 1:n_time) {
  Pn1 <- Pn1 %*% P1
  prob1[t + 1, ] <- pi0_1 %*% Pn1
}

matplot(0:n_time, prob1, type = "l", lwd = 2,
        xlab = "Time", ylab = "Probability",
        main = "State Probabilities Over Time (A1)",
        col = 1:5, lty = 1:5)
legend("right", legend = states1, lty = 1:5, col = 1:5)

# ============================================================
# QUESTION A2 — 7-STATE MARKOV CHAIN
# ============================================================

states2 <- paste0("S", 1:7)

P2 <- matrix(c(
  0,   1,   0,   0,   0,   0,   0,
  1,   0,   0,   0,   0,   0,   0,
  0,   0,   0,   0.4, 0.2, 0.2, 0.2,
  0,   0,   0,   0,   0.2, 0.4, 0.4,
  0.3, 0,   0,   0.1, 0.3, 0.1, 0.2,
  0,   0,   0,   0.2, 0.2, 0.3, 0.3,
  0,   0,   0,   0.5, 0.2, 0.2, 0.1
), nrow = 7, byrow = TRUE)

rownames(P2) <- colnames(P2) <- states2
mc2 <- new("markovchain", states = states2, transitionMatrix = P2)

plot(mc2, main = "7-State Markov Chain (A2)")

cat("\n--- A2(b): Classification ---\n")

print(recurrentClasses(mc2))
print(transientClasses(mc2))
print(absorbingStates(mc2))

cat("\nIs Irreducible:", is.irreducible(mc2), "\n")

cat("\nChain Period (recurrent class {S1,S2}):\n")
print(period(mc2))  # Should be 2

# ---- Simulation ----
set.seed(99)
trajs2 <- list(
  simulate_chain(P2, 50),
  simulate_chain(P2, 50)
)
plot_trajectories(trajs2, states2)

# ---- Steady state ----
cat("\n--- A2(d): Stationary Distribution ---\n")
ss2 <- steadyStates(mc2)
print(ss2)

cat("\nInterpretation:\n")
cat("Stationary distribution: (0.5, 0.5, 0, ..., 0)\n")
cat("Chain is periodic ⇒ probabilities oscillate, no true limit.\n")
cat("Not ergodic.\n")

# ============================================================
# QUESTION A3 — TRAFFIC MODEL
# ============================================================

states3 <- c("Light", "Heavy", "Jammed")

P_day <- matrix(c(
  0.4, 0.4, 0.2,
  0.3, 0.4, 0.3,
  0.0, 0.1, 0.9
), nrow = 3, byrow = TRUE)

P_peak <- matrix(c(
  0.1, 0.5, 0.4,
  0.1, 0.3, 0.6,
  0.0, 0.1, 0.9
), nrow = 3, byrow = TRUE)

rownames(P_day) <- colnames(P_day) <- states3
rownames(P_peak) <- colnames(P_peak) <- states3

# ---- Distribution ----
pi0_3 <- c(1, 0, 0)

pi_4pm <- pi0_3 %*% mat_power(P_day, 9)
pi_6pm <- pi_4pm %*% mat_power(P_peak, 6)

cat("\nDistribution at 6PM:\n")
print(round(pi_6pm, 4))

# ---- Simulation ----
set.seed(123)

simulate_traffic <- function() {
  state <- 1
  for (i in 1:9) state <- sample(1:3, 1, prob = P_day[state, ])
  for (i in 1:6) state <- sample(1:3, 1, prob = P_peak[state, ])
  state
}

results <- replicate(10000, simulate_traffic())
sim_probs <- prop.table(table(factor(results, levels = 1:3)))

cat("\nSimulated Distribution:\n")
print(round(sim_probs, 4))

cat("\n--- END ---\n")
