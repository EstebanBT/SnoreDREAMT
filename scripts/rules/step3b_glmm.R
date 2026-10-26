# Snakemake rule `step3b_glmm` (idea 4): hierarchical Poisson model of event counts per respiratory-phase bin.
#
#   count ~ cg + sg + cond:cg + cond:sg + (1 | sid:cond) + (0 + cg + sg | sid) + offset(log(occupancy))
#
# cg/sg = cos/sin(bin phase - group phase): phase modulation of the event rate along / orthogonal to the group
# preferred phase; (0 + cg + sg | sid) lets each subject have its own modulation (partial pooling instead of
# unreliable per-subject preferred phases); (1 | sid:cond) absorbs rate differences between conditions.
# beta_cg = modulation depth in clean breathing; gamma = cond:cg = change during the condition.
# Relative change = gamma / beta_cg, with a 95 % interval and P(loss > 50 %) from 20000 draws of the
# (approximately normal) sampling distribution of the fixed effects.
# Run for both exclusion levels: 'none' and 'evoked' (snore-evoked windows removed from events AND occupancy).
# Usage: Rscript step3b_glmm.R <phase_bins_N2.csv> <out.csv>
suppressPackageStartupMessages({library(lme4); library(MASS)})
args <- commandArgs(trailingOnly = TRUE)
d <- read.csv(args[1])
d <- d[d$primary == "True" & d$polarity_ok == "True" & d$occupancy_min > 0, ]
set.seed(1)
res <- list()
for (excl in unique(d$exclusion)) for (es in unique(d$event_set)) {
  for (cond in c("resp_event", "snoring", "hypopnea", "recovery", "clean_odd")) {
    ref <- if (cond == "clean_odd") "clean_even" else "clean"   # clean_odd vs clean_even = control
    x <- d[d$exclusion == excl & d$event_set == es & d$condition %in% c(ref, cond), ]
    # subjects need both conditions and some events in the condition
    tab <- aggregate(count ~ sid + condition, x, sum)
    ok <- intersect(tab$sid[tab$condition == ref & tab$count >= 20], tab$sid[tab$condition == cond & tab$count >= 5])
    x <- x[x$sid %in% ok, ]
    if (length(ok) < 8) next
    x$cond <- as.integer(x$condition == cond)
    x$cg <- cos(x$bin_center - x$mu_group)
    x$sg <- sin(x$bin_center - x$mu_group)
    x$sc <- interaction(x$sid, x$cond)
    x$cgc <- x$cg * x$cond; x$sgc <- x$sg * x$cond
    # random subject slopes for the condition x modulation terms too (otherwise the condition effects ignore
    # between-subject heterogeneity of the change and their SEs are too small); uncorrelated for stability
    fit <- tryCatch(glmer(count ~ cg + sg + cgc + sgc + (1 | sc) + (0 + cg + sg | sid) + (0 + cgc + sgc || sid),
                          offset = log(occupancy_min), family = poisson, data = x,
                          control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))),
                    error = function(e) NULL)
    if (is.null(fit)) next
    fe <- fixef(fit); V <- as.matrix(vcov(fit)); se <- sqrt(diag(V))
    draws <- mvrnorm(20000, fe, V)
    rel <- draws[, "cgc"] / draws[, "cg"]
    rot <- atan2(draws[, "sg"] + draws[, "sgc"], draws[, "cg"] + draws[, "cgc"]) - atan2(draws[, "sg"], draws[, "cg"])
    rot <- atan2(sin(rot), cos(rot)) * 180 / pi
    rot_mean <- atan2(mean(sin(rot * pi / 180)), mean(cos(rot * pi / 180))) * 180 / pi
    rot_c <- atan2(sin((rot - rot_mean) * pi / 180), cos((rot - rot_mean) * pi / 180)) * 180 / pi + rot_mean
    # modulation strength (vector length) in condition relative to reference: separates loss from phase rotation
    len_ratio <- sqrt((draws[, "cg"] + draws[, "cgc"])^2 + (draws[, "sg"] + draws[, "sgc"])^2) /
                 sqrt(draws[, "cg"]^2 + draws[, "sg"]^2)
    res[[length(res) + 1]] <- data.frame(
      exclusion = excl, event_set = es, variant = unique(x$variant), condition = cond, n_subjects = length(ok),
      n_events_condition = sum(x$count[x$cond == 1]), n_events_clean = sum(x$count[x$cond == 0]),
      beta_cg = fe["cg"], beta_cg_se = se["cg"], beta_cg_p = 2 * pnorm(-abs(fe["cg"] / se["cg"])),
      gamma_cg = fe["cgc"], gamma_cg_se = se["cgc"],
      gamma_cg_p_two_sided = 2 * pnorm(-abs(fe["cgc"] / se["cgc"])),
      gamma_cg_p_loss = pnorm(fe["cgc"] / se["cgc"]),
      gamma_sg = fe["sgc"], gamma_sg_p = 2 * pnorm(-abs(fe["sgc"] / se["sgc"])),
      beta_sg = fe["sg"], length_ratio = median(len_ratio), length_ratio_lo = quantile(len_ratio, .025),
      length_ratio_hi = quantile(len_ratio, .975), rotation_deg = rot_mean, rotation_lo = quantile(rot_c, .025), rotation_hi = quantile(rot_c, .975),
      rel_change_median = median(rel), rel_change_lo = quantile(rel, .025), rel_change_hi = quantile(rel, .975),
      p_loss_gt_50pct = mean(rel < -0.5), p_complete_loss = mean(rel < -1),
      singular = isSingular(fit), row.names = NULL)
  }
}
out <- do.call(rbind, res)
write.csv(out, args[2], row.names = FALSE)
print(format(out, digits = 3), row.names = FALSE)
