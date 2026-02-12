# ============================================================
# ROLLING 1-STEP VAR-LASSO (BigVAR) - TUTTE LE STAZIONI
# + PROGRESS BAR (progress)
# - Rolling 1-step con IC 95%
# - Baseline pred=0 su differenze
# - 2 plot ggplot sul LIVELLO (solo top 3 per visual)
# ============================================================

# -------------------------
# Pacchetti
# -------------------------
suppressPackageStartupMessages({
  library(BigVAR)
  library(ggplot2)
  library(dplyr)
})

if (!requireNamespace("progress", quietly = TRUE)) {
  install.packages("progress")
}
suppressPackageStartupMessages(library(progress))

# =========================
# CONFIG
# =========================
N_TEST <- 36
GRAN <- c(20, 3)
P_ORDER <- 1

cat("========================================\n")
cat("ROLLING 1-STEP VAR-LASSO (ALL STATIONS)\n")
cat("========================================\n\n")

# =========================
# 1) CARICAMENTO + CLEAN
# =========================
cat("[1/7] Caricamento dati...\n")

pm25_data <- read.csv("TaiwanAirBox032017.csv")
X_raw <- as.matrix(pm25_data[, -1, drop = FALSE])
colnames(X_raw) <- paste0("S", 1:ncol(X_raw))

zero_counts <- colSums(X_raw == 0, na.rm = TRUE)
bad_stations <- which(zero_counts > quantile(zero_counts, 0.95))
X_clean <- if (length(bad_stations) > 0) X_raw[, -bad_stations, drop = FALSE] else X_raw

cat("   Tot stazioni iniziali:", ncol(X_raw), "\n")
cat("   Stazioni rimosse:", length(bad_stations), "\n")
cat("   Stazioni rimaste:", ncol(X_clean), "\n")
cat("   Completato!\n\n")

# =========================
# 2) DIFFERENZE (TUTTE)
# =========================
cat("[2/7] Differenziazione...\n")

X_diff <- diff(X_clean)
colnames(X_diff) <- colnames(X_clean)

T <- nrow(X_diff)
k <- ncol(X_diff)

cat("   Dimensioni X_diff:", T, "x", k, "\n")
cat("   Completato!\n\n")

# =========================
# 3) SETUP CUT + INDICI
# =========================
cat("[3/7] Setup rolling + CUT...\n")

N_TEST <- min(N_TEST, T - 1)

cut_index <- T - N_TEST
time_train <- 1:cut_index
time_test  <- (cut_index + 1):T

stopifnot(length(time_test) == N_TEST, max(time_test) == T)

cat("   T =", T, "\n")
cat("   k =", k, "\n")
cat("   N_TEST =", N_TEST, "\n")
cat("   CUT (fine train) =", cut_index, "\n")
cat("   Test da", min(time_test), "a", max(time_test), "\n")
cat("   Completato!\n\n")

# =========================
# 4) ROLLING 1-STEP con PROGRESS BAR
# =========================
cat("[4/7] Rolling forecast (fit CV + predict + CI)...\n")
cat("   Nota: con k grande può essere lento.\n\n")

pred_mean  <- matrix(NA_real_, N_TEST, k)
pred_lower <- matrix(NA_real_, N_TEST, k)
pred_upper <- matrix(NA_real_, N_TEST, k)
y_true     <- matrix(NA_real_, N_TEST, k)

colnames(pred_mean) <- colnames(pred_lower) <- colnames(pred_upper) <- colnames(y_true) <- colnames(X_diff)

# progress bar
pb <- progress_bar$new(
  format = "  [:bar] :percent | :current/:total | elapsed=:elapsed | eta=:eta | :message",
  total = N_TEST,
  clear = FALSE,
  width = 90
)

start_all <- Sys.time()

for (idx in 1:N_TEST) {
  
  # t = ultimo indice in train, target = osservazione da predire
  t <- cut_index + idx - 1
  target_index <- t + 1
  
  # Training crescente
  Ytrain <- X_diff[1:t, , drop = FALSE]
  y_true[idx, ] <- X_diff[target_index, ]
  
  pb$message(sprintf("idx=%d fit CV (t=%d -> pred %d)", idx, t, target_index))
  
  tryCatch({
    
    TT <- nrow(Ytrain)
    
    # clamp T1/T2 per Rolling CV: T1 < T2 < TT
    T1 <- max(10, floor(0.60 * TT))
    T2 <- max(T1 + 5, floor(0.80 * TT))
    if (T2 >= TT) T2 <- TT - 1
    if (T1 >= T2) T1 <- T2 - 1
    if (T1 < 5) T1 <- 5
    
    Model <- constructModel(
      Y = Ytrain,
      p = P_ORDER,
      struct = "Basic",
      gran = GRAN,
      cv = "Rolling",
      T1 = T1,
      T2 = T2,
      verbose = FALSE,
      model.controls = list(intercept = TRUE)
    )
    
    fit <- cv.BigVAR(Model)
    
    pb$message(sprintf("idx=%d predict+CI", idx))
    
    # predict robusto
    fc <- predict(fit, n.ahead = 1)
    
    if (is.matrix(fc)) {
      if (nrow(fc) == k) {
        yhat <- fc[, 1]
      } else if (ncol(fc) == k) {
        yhat <- fc[1, ]
      } else {
        stop("predict() matrix con dimensioni inattese.")
      }
    } else if (is.array(fc)) {
      yhat <- as.numeric(fc)[1:k]
    } else {
      yhat <- as.numeric(fc)
      if (length(yhat) != k) stop("predict() vettore con lunghezza errata.")
    }
    
    pred_mean[idx, ] <- as.numeric(yhat)
    
    # IC 95% su differenze (residui del fit)
    res <- fit@resids
    Sigma_hat <- cov(res)
    se <- sqrt(diag(Sigma_hat))
    
    pred_lower[idx, ] <- pred_mean[idx, ] - 1.96 * se
    pred_upper[idx, ] <- pred_mean[idx, ] + 1.96 * se
    
  }, error = function(e) {
    pb$message(sprintf("idx=%d ERROR: %s", idx, e$message))
  })
  
  pb$tick()
}

elapsed_total <- as.numeric(difftime(Sys.time(), start_all, units = "mins"))
cat("\nRolling completato in", round(elapsed_total, 2), "minuti.\n")
cat("NA pred_mean =", sum(is.na(pred_mean)), "su", length(pred_mean), "\n\n")

# =========================
# 5) PERFORMANCE + BASELINE
# =========================
cat("[5/7] Performance...\n")

errors <- y_true - pred_mean
mse <- mean(errors^2, na.rm = TRUE)
mae <- mean(abs(errors), na.rm = TRUE)
rmse <- sqrt(mse)

inside_ci <- (y_true >= pred_lower) & (y_true <= pred_upper)
coverage <- mean(inside_ci, na.rm = TRUE) * 100

cor_pred <- cor(as.vector(y_true), as.vector(pred_mean), use = "complete.obs")

cat("   MSE:", round(mse, 4), "\n")
cat("   RMSE:", round(rmse, 4), "\n")
cat("   MAE:", round(mae, 4), "\n")
cat("   Correlazione:", round(cor_pred, 3), "\n")
cat("   Coverage IC 95%:", round(coverage, 1), "%\n\n")

pred0 <- matrix(0, nrow = N_TEST, ncol = k)
mse_model <- colMeans((y_true - pred_mean)^2, na.rm = TRUE)
mse_0     <- colMeans((y_true - pred0)^2, na.rm = TRUE)
improvement_pct <- 100 * (mse_0 - mse_model) / mse_0

df_cmp <- data.frame(
  station = colnames(X_diff),
  mse_model = mse_model,
  mse_baseline0 = mse_0,
  improvement_pct = improvement_pct
)

cat("   Completato!\n\n")

# =========================
# 6) LIVELLO per plotting
# =========================
cat("[6/7] Ricostruzione livello per plotting...\n")

X_level <- X_clean[-1, 1:k, drop = FALSE]
stopifnot(nrow(X_level) == T)

level_true <- X_level[time_test, , drop = FALSE]
level_pred <- X_level[time_test - 1, , drop = FALSE] + pred_mean
level_lower <- X_level[time_test - 1, , drop = FALSE] + pred_lower
level_upper <- X_level[time_test - 1, , drop = FALSE] + pred_upper

cat("   Completato!\n\n")

# =========================
# 7) 2 PLOT ggplot (top 3 per improvement)
# =========================
cat("[7/7] Plot ggplot (top 3 stazioni per improvement)...\n")

top3 <- order(df_cmp$improvement_pct, decreasing = TRUE)[1:3]
station_names <- colnames(X_diff)[top3]
imp_lab <- round(df_cmp$improvement_pct[top3], 1)

# DF FULL
df_full <- data.frame(
  time = rep(1:T, times = length(top3)),
  station = rep(station_names, each = T),
  value = as.vector(X_level[, top3, drop = FALSE])
) %>%
  mutate(
    phase = ifelse(time <= cut_index, "Observed (Train)", "Observed (Test)"),
    panel = paste0(station, " | imp=", rep(imp_lab, each = T), "%")
  )

# DF ZOOM
df_zoom <- data.frame(
  time = rep(time_test, times = length(top3)),
  station = rep(station_names, each = length(time_test)),
  observed = as.vector(level_true[, top3, drop = FALSE]),
  predicted = as.vector(level_pred[, top3, drop = FALSE]),
  lower = as.vector(level_lower[, top3, drop = FALSE]),
  upper = as.vector(level_upper[, top3, drop = FALSE])
) %>%
  mutate(
    panel = paste0(station, " | imp=", rep(imp_lab, each = length(time_test)), "%")
  )

p_full <- ggplot(df_full, aes(x = time, y = value, color = phase)) +
  geom_line(linewidth = 0.6) +
  geom_vline(xintercept = cut_index, linetype = "dashed", linewidth = 0.9) +
  facet_wrap(~ panel, ncol = 1, scales = "free_y") +
  labs(
    title = paste0("FULL (Livello) | CUT=", cut_index, " | Test=", cut_index + 1, "…", T),
    x = "Time index (allineato a X_diff / X_level)",
    y = "PM2.5 (livello)",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

print(p_full)

p_zoom <- ggplot(df_zoom, aes(x = time)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25) +
  geom_line(aes(y = observed, color = "Observed (Test)"), linewidth = 0.9) +
  geom_line(aes(y = predicted, color = "Predicted"), linewidth = 0.9) +
  geom_vline(xintercept = cut_index, linetype = "dashed", linewidth = 0.9) +
  facet_wrap(~ panel, ncol = 1, scales = "free_y") +
  labs(
    title = paste0("ZOOM TEST (Livello) | ultimi N_TEST=", N_TEST, " punti"),
    x = "Time index (test window)",
    y = "PM2.5 (livello)",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

print(p_zoom)

cat("   Plot completati!\n\n")

# =========================
# SALVA RISULTATI
# =========================
results_all <- list(
  T = T, k = k,
  cut_index = cut_index,
  time_test = time_test,
  pred_mean = pred_mean,
  pred_lower = pred_lower,
  pred_upper = pred_upper,
  y_true = y_true,
  performance = list(mse = mse, rmse = rmse, mae = mae, cor = cor_pred, coverage = coverage),
  baseline0 = df_cmp,
  config = list(N_TEST = N_TEST, GRAN = GRAN, P_ORDER = P_ORDER)
)

saveRDS(results_all, "Rolling_VARLASSO_allstations_results.rds")
cat("Risultati salvati: Rolling_VARLASSO_allstations_results.rds\n")
