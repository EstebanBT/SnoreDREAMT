# Snakemake rule `step4_glmm`: does severity moderate the respiratory phase modulation of event rates (clean N2)?
#   count ~ (cg + sg) * (z_log_ahi + z_AGE + male + z_BMI) + (1 | sid) + (0 + cg + sg | sid) + offset(log(occupancy))
# cg:z_log_ahi = change of modulation (along the group phase) per SD of log1p(AHI); relative to the mean modulation cg.
# Usage: Rscript step4_glmm.R <glmm_data.csv> <out.csv>
suppressPackageStartupMessages({library(lme4); library(MASS)})
args <- commandArgs(trailingOnly = TRUE)
d <- read.csv(args[1])
d <- d[d$polarity_ok == "True" & d$occupancy_min > 0, ]
set.seed(1)
res <- list()
for (es in unique(d$event_set)) {
  x <- d[d$event_set == es, ]
  tab <- aggregate(count ~ sid, x, sum)
  x <- x[x$sid %in% tab$sid[tab$count >= 10], ]
  x$cg <- cos(x$bin_center - x$mu_group); x$sg <- sin(x$bin_center - x$mu_group)
  fit <- tryCatch(glmer(count ~ (cg + sg) * (z_log_ahi + z_AGE + male + z_BMI) + (1 | sid) + (0 + cg + sg | sid),
                        offset = log(occupancy_min), family = poisson, data = x,
                        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))),
                  error = function(e) NULL)
  if (is.null(fit)) next
  fe <- fixef(fit); V <- as.matrix(vcov(fit)); se <- sqrt(diag(V))
  dr <- mvrnorm(20000, fe, V)
  rel <- dr[, "cg:z_log_ahi"] / dr[, "cg"]
  res[[length(res) + 1]] <- data.frame(
    event_set = es, variant = unique(x$variant), n_subjects = length(unique(x$sid)), n_events = sum(x$count),
    beta_cg = fe["cg"], beta_cg_p = 2 * pnorm(-abs(fe["cg"] / se["cg"])),
    ahi_x_cg = fe["cg:z_log_ahi"], ahi_x_cg_se = se["cg:z_log_ahi"],
    ahi_x_cg_p = 2 * pnorm(-abs(fe["cg:z_log_ahi"] / se["cg:z_log_ahi"])),
    ahi_x_sg = fe["sg:z_log_ahi"], ahi_x_sg_p = 2 * pnorm(-abs(fe["sg:z_log_ahi"] / se["sg:z_log_ahi"])),
    rel_change_per_sd = median(rel), rel_lo = quantile(rel, .025), rel_hi = quantile(rel, .975),
    age_x_cg = fe["cg:z_AGE"], age_x_cg_p = 2 * pnorm(-abs(fe["cg:z_AGE"] / se["cg:z_AGE"])),
    singular = isSingular(fit), row.names = NULL)
}
out <- do.call(rbind, res)
write.csv(out, args[2], row.names = FALSE)
print(format(out, digits = 3), row.names = FALSE)
