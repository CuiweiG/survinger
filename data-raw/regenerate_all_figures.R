# ============================================================
# REGENERATE ALL 7 FIGURES FROM SCRATCH
# Every figure uses survinger package functions on real/simulated data
# No manual data, no fake numbers, no post-hoc editing
# ============================================================

.libPaths(Sys.getenv("R_LIBS_USER"))
setwd("C:/Users/openclaw/survinger")
devtools::load_all()

cols <- .surv_colors()
fig_dir <- "man/figures"
dir.create(fig_dir, showWarnings = FALSE)

# ============================================================
# LOAD REAL DATA
# ============================================================
coguk <- readRDS("data-raw/coguk_surveillance.rds")
n_total <- coguk$n_sequences
target_lin <- "XBB.1.5"

# All survinger calls on real data
design <- surv_design(
  data = coguk$sequences, strata = ~ region,
  sequencing_rate = coguk$population[c("region", "seq_rate")],
  population = coguk$population, source_type = "source_type"
)
weighted <- surv_lineage_prevalence(design, target_lin, method = "hajek")
naive    <- surv_naive_prevalence(design, target_lin)
delay    <- surv_estimate_delay(design)
nc       <- surv_nowcast_lineage(design, delay, target_lin)
adj      <- surv_adjusted_prevalence(design, delay, target_lin)
alloc    <- surv_optimize_allocation(design, "min_mse", total_capacity = 2000)
report   <- surv_report(design, target_lineage = target_lin)

diff_mean <- round(mean(abs(weighted$estimates$prevalence -
                            naive$estimates$prevalence), na.rm = TRUE), 4)
max_diff  <- round(max(abs(weighted$estimates$prevalence -
                           naive$estimates$prevalence), na.rm = TRUE), 4)

cat("Design:   ", design$n_obs, "sequences,", design$n_strata, "strata\n")
cat("Gini:     ", round(report$gini, 3), "\n")
cat("Mean diff:", diff_mean * 100, "pp\n")
cat("Max diff: ", max_diff * 100, "pp\n\n")

# Helper: shortened time labels
shorten <- function(x) gsub("2023-", "", x)
thin_labels <- function(vals, every = 3) {
  vals[seq(1, length(vals), by = every)]
}

# ============================================================
# FIG 1: SEQUENCING INEQUALITY — surv_design$strata_info
# ============================================================
info <- design$strata_info
info$label <- gsub("UK-", "", info$region)
info <- info[order(info$seq_rate), ]
info$label <- factor(info$label, levels = info$label)
ratio <- round(max(info$seq_rate) / min(info$seq_rate), 1)

p1 <- ggplot2::ggplot(info, ggplot2::aes(x = .data$label, y = .data$seq_rate)) +
  ggplot2::geom_segment(ggplot2::aes(xend = .data$label, y = 0, yend = .data$seq_rate),
                        color = cols[["primary"]], linewidth = 1.5) +
  ggplot2::geom_point(size = 5, color = cols[["primary"]]) +
  ggplot2::geom_text(ggplot2::aes(label = paste0(round(.data$seq_rate * 100, 1), "%")),
                     hjust = -0.3, size = 4, color = "#333333") +
  ggplot2::geom_hline(yintercept = mean(info$seq_rate), linetype = "dashed",
                      color = cols[["secondary"]], linewidth = 0.5) +
  ggplot2::scale_y_continuous(labels = function(x) paste0(round(x*100,1), "%"),
                              expand = ggplot2::expansion(mult = c(0, 0.25))) +
  ggplot2::coord_flip() +
  ggplot2::labs(title = paste0("Sequencing rates vary ", ratio, "-fold across UK nations"),
                subtitle = paste0("COG-UK data, 2023 H1 (n = ", formatC(n_total, big.mark=","), ")"),
                x = NULL, y = "Sequencing rate", caption = "Data: COG-UK | Dashed: mean") +
  theme_survinger(base_size = 12)

grDevices::png(file.path(fig_dir, "fig1_inequality.png"), width=2400, height=1400, res=300)
print(p1); grDevices::dev.off(); cat("fig1 OK\n")

# ============================================================
# FIG 2: WEIGHTED VS NAIVE — surv_lineage_prevalence output
# ============================================================
w_df <- weighted$estimates[!is.na(weighted$estimates$prevalence), ]
w_df$Method <- "Design-weighted (Hajek)"
n_df <- naive$estimates[!is.na(naive$estimates$prevalence), ]
n_df$Method <- "Naive (unweighted)"
combined <- dplyr::bind_rows(w_df, n_df)
combined$ts <- shorten(combined$time)
labels2 <- thin_labels(sort(unique(combined$ts)))

p2 <- ggplot2::ggplot(combined, ggplot2::aes(
  x = .data$ts, y = .data$prevalence,
  color = .data$Method, fill = .data$Method, group = .data$Method)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper),
                       alpha = 0.10, color = NA) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_point(size = 1.8) +
  ggplot2::scale_color_manual(values = c(cols[["primary"]], cols[["secondary"]])) +
  ggplot2::scale_fill_manual(values = c(cols[["primary"]], cols[["secondary"]])) +
  ggplot2::scale_x_discrete(breaks = labels2) +
  ggplot2::scale_y_continuous(labels = function(x) paste0(round(x*100), "%"), limits = c(0, NA)) +
  ggplot2::annotate("label", x = labels2[2],
                    y = max(w_df$prevalence, na.rm=TRUE) * 0.95,
                    label = paste0("Mean |bias|: ", round(diff_mean*100,1),
                                   " pp\nMax |bias|: ", round(max_diff*100,1), " pp"),
                    hjust = 0, size = 3.2, fill = "#FFFFFFCC", color = "#333333") +
  ggplot2::labs(title = paste0(target_lin, ": design weighting corrects surveillance bias"),
                subtitle = paste0("COG-UK, 4 UK nations, 2023 H1 (n = ", formatC(n_total, big.mark=","), ")"),
                x = "Epiweek", y = "Prevalence", color = NULL, fill = NULL,
                caption = "95% Wald CI | Data: COG-UK") +
  theme_survinger() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1, size=9),
                 legend.position = c(0.70, 0.90),
                 legend.background = ggplot2::element_rect(fill="#FFFFFFDD", color=NA))

grDevices::png(file.path(fig_dir, "fig2_compare.png"), width=3000, height=1800, res=300)
print(p2); grDevices::dev.off(); cat("fig2 OK\n")

# ============================================================
# FIG 3: BIAS HEATMAP — per-region naive vs global weighted
# ============================================================
per_rw <- design$data
per_rw$.target <- as.integer(per_rw$lineage == target_lin)
per_rw <- per_rw |>
  dplyr::group_by(.data$region, .data$epiweek) |>
  dplyr::summarise(naive_prev = mean(.data$.target), .groups = "drop")
gw <- weighted$estimates[c("time","prevalence")]
names(gw) <- c("epiweek","weighted_prev")
per_rw <- dplyr::left_join(per_rw, gw, by = "epiweek")
per_rw$bias <- per_rw$naive_prev - per_rw$weighted_prev
per_rw$rs <- gsub("UK-", "", per_rw$region)
per_rw$es <- shorten(per_rw$epiweek)
ro <- design$strata_info[order(design$strata_info$seq_rate, decreasing=TRUE), ]
per_rw$rs <- factor(per_rw$rs, levels = gsub("UK-","", ro$region))
labels3 <- thin_labels(sort(unique(per_rw$es)))

p3 <- ggplot2::ggplot(per_rw, ggplot2::aes(x=.data$es, y=.data$rs, fill=.data$bias)) +
  ggplot2::geom_tile(color="white", linewidth=0.5) +
  ggplot2::scale_fill_gradient2(low=cols[["primary"]], mid="#FFFFFF", high=cols[["secondary"]],
    midpoint=0, limits=c(-0.15,0.15), oob=scales::squish,
    labels=function(x) paste0(ifelse(x>0,"+",""), round(x*100,1)," pp"), name="Bias") +
  ggplot2::scale_x_discrete(breaks = labels3) +
  ggplot2::labs(title = paste0("Per-nation estimation bias for ", target_lin),
                subtitle = "Ordered by sequencing rate (highest top)",
                x = "Epiweek", y = NULL, caption = "Data: COG-UK") +
  theme_survinger() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1, size=8),
                 legend.position = "right", legend.key.height = ggplot2::unit(40,"pt"))

grDevices::png(file.path(fig_dir, "fig3_bias_heatmap.png"), width=3200, height=1200, res=300)
print(p3); grDevices::dev.off(); cat("fig3 OK\n")

# ============================================================
# FIG 4: NOWCAST — surv_simulate + surv_nowcast_lineage
# (simulated because COG-UK has no report_date)
# ============================================================
set.seed(2024)
sim <- surv_simulate(n_regions=5, n_weeks=20,
  sequencing_rates=c(0.15,0.08,0.03,0.01,0.005),
  delay_params=list(mu=12, size=3), seed=2024)
ref <- max(sim$sequences$collection_date)
trunc_seqs <- sim$sequences[sim$sequences$report_date <= ref, ]
d_trunc <- surv_design(trunc_seqs, ~region,
                       sim$population[c("region","seq_rate")], sim$population)
dfit <- surv_estimate_delay(
  surv_design(sim$sequences, ~region,
              sim$population[c("region","seq_rate")], sim$population))
nc_sim <- surv_nowcast_lineage(d_trunc, dfit, "BA.2.86", ref_date = ref)

full_counts <- sim$sequences[sim$sequences$lineage == "BA.2.86", ] |>
  dplyr::group_by(.data$epiweek) |>
  dplyr::summarise(n_full = dplyr::n(), .groups = "drop")
nc_df <- dplyr::left_join(nc_sim$estimates, full_counts, by=c("time"="epiweek"))
nc_df$n_full[is.na(nc_df$n_full)] <- 0L
labels4 <- thin_labels(sort(unique(nc_df$time)))
nc_start <- nc_df$time[which(nc_df$is_nowcast)[1]]

p4 <- ggplot2::ggplot(nc_df, ggplot2::aes(x=.data$time)) +
  ggplot2::geom_col(ggplot2::aes(y=.data$n_full), fill="#E0E0E0", width=0.7) +
  ggplot2::geom_col(ggplot2::aes(y=.data$n_observed), fill=cols[["neutral"]], alpha=0.8, width=0.5) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin=pmax(0,.data$n_estimated-1.96*.data$se),
                                    ymax=.data$n_estimated+1.96*.data$se, group=1),
                       fill=cols[["quaternary"]], alpha=0.12) +
  ggplot2::geom_line(ggplot2::aes(y=.data$n_estimated, group=1),
                     color=cols[["quaternary"]], linewidth=1.1) +
  ggplot2::geom_point(ggplot2::aes(y=.data$n_estimated,
    shape=ifelse(.data$is_nowcast,"Nowcasted","Complete")),
    color=cols[["quaternary"]], size=2.5) +
  ggplot2::scale_shape_manual(values=c("Complete"=16,"Nowcasted"=17)) +
  ggplot2::scale_x_discrete(breaks=labels4) +
  {if(!is.na(nc_start)) ggplot2::geom_vline(xintercept=nc_start, linetype="dashed", color="#999")} +
  ggplot2::labs(title="Delay-adjusted nowcasting recovers incomplete counts",
                subtitle=paste0("Simulated data | NegBin(\u03bc=",round(dfit$parameters$mu[1],1),
                                ", k=",round(dfit$parameters$size[1],1),") delay"),
                x="Epiweek", y="BA.2.86 count", shape=NULL,
                caption="Light grey: full | Dark grey: reported | Orange: nowcast") +
  theme_survinger() +
  ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1,size=8),
                 legend.position=c(0.15,0.85))

grDevices::png(file.path(fig_dir,"fig4_nowcast.png"), width=3000, height=1800, res=300)
print(p4); grDevices::dev.off(); cat("fig4 OK\n")

# ============================================================
# FIG 5: ALLOCATION — surv_optimize_allocation output
# ============================================================
H <- design$n_strata
strats <- c("equal","proportional","min_mse")
ad <- purrr::map_dfr(strats, function(s) {
  if (s=="equal") {
    n <- rep(2000L%/%H, H); n[1] <- n[1]+2000L-sum(n)
    tibble::tibble(region=gsub("UK-","",design$strata_info$region), n=n, strategy=s)
  } else if (s=="proportional") {
    pi <- design$population$pop_total/sum(design$population$pop_total)
    n <- as.integer(round(pi*2000)); n[1] <- n[1]+2000L-sum(n)
    tibble::tibble(region=gsub("UK-","",design$population$region), n=n, strategy=s)
  } else {
    r <- surv_optimize_allocation(design, s, total_capacity=2000)
    tibble::tibble(region=gsub("UK-","",r$allocation$region), n=r$allocation$n_allocated, strategy=s)
  }
})
ad$strategy <- factor(ad$strategy, levels=c("equal","proportional","min_mse"),
                      labels=c("Equal","Population-proportional","MSE-optimal"))

p5 <- ggplot2::ggplot(ad, ggplot2::aes(x=.data$region, y=.data$n, fill=.data$strategy)) +
  ggplot2::geom_col(position=ggplot2::position_dodge(0.8), width=0.7) +
  ggplot2::scale_fill_manual(values=c(cols[["neutral"]],cols[["light_blue"]],cols[["primary"]])) +
  ggplot2::labs(title="Optimal allocation concentrates resources by population",
                subtitle="2,000 weekly sequences across 4 UK nations",
                x=NULL, y="Sequences/week", fill="Strategy",
                caption="Data: COG-UK") +
  theme_survinger()

grDevices::png(file.path(fig_dir,"fig5_allocation.png"), width=2800, height=1600, res=300)
print(p5); grDevices::dev.off(); cat("fig5 OK\n")

# ============================================================
# FIG 6: BENCHMARK — surv_simulate + surv_lineage_prevalence
# ============================================================
cat("Running benchmark simulation...\n")
set.seed(2024)
n_sim <- 50
gini_levels <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6)

br <- purrr::map_dfr(gini_levels, function(g) {
  purrr::map_dfr(seq_len(n_sim), function(i) {
    rates <- stats::rbeta(5, 1, max(0.5,(1-g)/g))
    rates <- sort(pmax(pmin(rates,0.5),0.005), decreasing=TRUE)
    sim <- surv_simulate(n_regions=5, n_weeks=12, sequencing_rates=rates,
                         seed=i*1000+round(g*100))
    dat <- sim$sequences
    regions <- sort(unique(dat$region))
    prev_by_r <- seq(0.3, 0.05, length.out=5)
    names(prev_by_r) <- regions
    for (r in regions) {
      idx <- which(dat$region == r)
      dat$lineage[idx] <- sample(c("BA.2.86","Other"), length(idx), replace=TRUE,
                                 prob=c(prev_by_r[r], 1-prev_by_r[r]))
    }
    pop_shares <- sim$population$pop_total / sum(sim$population$pop_total)
    true_prev <- sum(pop_shares * prev_by_r[sim$population$region])
    sim$sequences <- dat
    d <- tryCatch(surv_design(dat, ~region, sim$population[c("region","seq_rate")],
                              sim$population), error=function(e) NULL)
    if (is.null(d)) return(NULL)
    w_dat <- .map_weights_to_obs(d)
    w_dat$.t <- as.integer(w_dat$lineage == "BA.2.86")
    tibble::tibble(gini=g, sim=i,
      abs_bias_w = abs(sum(w_dat$weight * w_dat$.t)/sum(w_dat$weight) - true_prev),
      abs_bias_n = abs(mean(w_dat$.t) - true_prev))
  })
})

bs <- br |> dplyr::group_by(.data$gini) |>
  dplyr::summarise(m_w=mean(.data$abs_bias_w), m_n=mean(.data$abs_bias_n),
                   se_w=stats::sd(.data$abs_bias_w)/sqrt(dplyr::n()),
                   se_n=stats::sd(.data$abs_bias_n)/sqrt(dplyr::n()), .groups="drop")

bl <- tidyr::pivot_longer(bs, cols=c("m_w","m_n"), names_to="method", values_to="bias")
bl$se <- ifelse(bl$method=="m_w", bs$se_w, bs$se_n)
bl$Method <- ifelse(bl$method=="m_w", "Hajek (design-weighted)", "Naive (unweighted)")

p6 <- ggplot2::ggplot(bl, ggplot2::aes(x=.data$gini, y=.data$bias,
  color=.data$Method, fill=.data$Method)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin=pmax(0,.data$bias-1.96*.data$se),
                                    ymax=.data$bias+1.96*.data$se), alpha=0.15, color=NA) +
  ggplot2::geom_line(linewidth=1.1) + ggplot2::geom_point(size=3.5) +
  ggplot2::scale_color_manual(values=c(cols[["primary"]],cols[["secondary"]])) +
  ggplot2::scale_fill_manual(values=c(cols[["primary"]],cols[["secondary"]])) +
  ggplot2::scale_x_continuous(breaks=gini_levels) +
  ggplot2::scale_y_continuous(labels=function(x) paste0(round(x*100,1)," pp")) +
  ggplot2::labs(title="Design weighting reduces bias under heterogeneous prevalence",
                subtitle=paste0("5 regions, prevalence 5%-30%, ", n_sim, " reps/level"),
                x="Gini coefficient", y="Mean absolute bias (pp)", color=NULL, fill=NULL,
                caption="Shaded: 95% CI | Prevalence correlated with sequencing rate") +
  theme_survinger() +
  ggplot2::theme(legend.position=c(0.70,0.88),
                 legend.background=ggplot2::element_rect(fill="#FFFFFFDD",color=NA))

grDevices::png(file.path(fig_dir,"fig6_benchmark.png"), width=2800, height=1800, res=300)
print(p6); grDevices::dev.off(); cat("fig6 OK\n")

# ============================================================
# FIG 7: DETECTION — surv_detection_probability on COG-UK
# ============================================================
prev_seq <- c(seq(0.0001,0.001,0.0001), seq(0.001,0.01,0.0005), seq(0.01,0.02,0.001))
prev_seq <- sort(unique(prev_seq))
det_probs <- vapply(prev_seq, function(p) surv_detection_probability(design,p)$overall, numeric(1))
det_df <- tibble::tibble(prevalence=prev_seq*100, detection=det_probs*100)
p95 <- prev_seq[which(det_probs>=0.95)[1]]*100
p50 <- prev_seq[which(det_probs>=0.50)[1]]*100
x_max <- min(max(prev_seq*100), p95*3)
ddc <- det_df[det_df$prevalence <= x_max, ]

p7 <- ggplot2::ggplot(ddc, ggplot2::aes(x=.data$prevalence, y=.data$detection)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin=0, ymax=.data$detection), fill=cols[["light_blue"]], alpha=0.25) +
  ggplot2::geom_line(color=cols[["primary"]], linewidth=1.2) +
  ggplot2::geom_hline(yintercept=95, linetype="dashed", color=cols[["secondary"]], linewidth=0.5) +
  ggplot2::geom_hline(yintercept=50, linetype="dotted", color=cols[["neutral"]], linewidth=0.4) +
  ggplot2::annotate("point", x=p95, y=95, size=4, color=cols[["secondary"]]) +
  ggplot2::annotate("label", x=p95, y=87,
    label=paste0("95% detection\nat ",sprintf("%.2f",p95),"% prevalence"),
    size=3, fill="#FFFFFFDD", color=cols[["secondary"]]) +
  ggplot2::scale_x_continuous(labels=function(x) paste0(x,"%")) +
  ggplot2::scale_y_continuous(labels=function(x) paste0(x,"%"), breaks=c(0,25,50,80,95,100)) +
  ggplot2::labs(title="Detection probability under current COG-UK design",
                subtitle=paste0("n = ",formatC(n_total,big.mark=",")," sequences over 26 weeks"),
                x="True prevalence", y="Weekly detection probability",
                caption="Data: COG-UK | Per-stratum sequencing volumes") +
  theme_survinger()

grDevices::png(file.path(fig_dir,"fig7_detection.png"), width=2800, height=1800, res=300)
print(p7); grDevices::dev.off(); cat("fig7 OK\n")

# Copy all to vignettes/figures
for (f in list.files(fig_dir, "\\.png$", full.names=TRUE))
  file.copy(f, file.path("vignettes/figures", basename(f)), overwrite=TRUE)

cat("\n=== ALL 7 FIGURES REGENERATED FROM SCRATCH ===\n")
figs <- list.files(fig_dir, "\\.png$")
for (f in figs) cat(sprintf("  %s (%d KB)\n", f, round(file.info(file.path(fig_dir,f))$size/1024)))
