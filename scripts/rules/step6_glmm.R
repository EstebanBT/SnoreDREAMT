# Snakemake rule `step6_glmm`: does severity (z log1p AHI) moderate the phase modulation in clean breathing and its
# change during each breathing condition?
#   count ~ (cg + sg + cgc + sgc) * z_log_ahi + (1 | sid:cond) + (0 + cg + sg | sid) + (0 + cgc + sgc || sid)
#           + offset(log(occupancy))
# cg/sg = cos/sin(phase - group phase); cgc/sgc = condition x cg/sg.
#   sg:z  -> clean preferred phase rotates with AHI (6b, population level)
#   sgc:z -> the condition-related rotation depends on AHI (6a, population level)
#   cgc:z -> the condition-related change in modulation along the clean direction depends on AHI
# Usage: Rscript step6_glmm.R <glmm_data.csv> <out.csv>
suppressPackageStartupMessages(library(lme4))
args <- commandArgs(trailingOnly = TRUE)
d <- read.csv(args[1])
d <- d[d$occupancy_min > 0, ]
res <- list()
for (excl in unique(d$exclusion)) for (es in unique(d$event_set)) for (cond in c("resp_event", "snoring", "hypopnea", "recovery")) {
  if (excl == "evoked" && cond != "snoring") next
  x <- d[d$exclusion == excl & d$event_set == es & d$condition %in% c("clean", cond), ]
  tab <- aggregate(count ~ sid + condition, x, sum)
  ok <- intersect(tab$sid[tab$condition == "clean" & tab$count >= 20], tab$sid[tab$condition == cond & tab$count >= 5])
  x <- x[x$sid %in% ok, ]
  if (length(ok) < 10) next
  x$cond <- as.integer(x$condition == cond)
  x$cg <- cos(x$bin_center - x$mu_group); x$sg <- sin(x$bin_center - x$mu_group)
  x$cgc <- x$cg * x$cond; x$sgc <- x$sg * x$cond; x$sc <- interaction(x$sid, x$cond)
  fit <- tryCatch(glmer(count ~ (cg + sg + cgc + sgc) * z_log_ahi + (1 | sc) + (0 + cg + sg | sid) + (0 + cgc + sgc || sid),
                        offset = log(occupancy_min), family = poisson, data = x,
                        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))),
                  error = function(e) NULL)
  if (is.null(fit)) next
  fe <- fixef(fit); se <- sqrt(diag(as.matrix(vcov(fit))))
  pv <- function(k) 2 * pnorm(-abs(fe[k] / se[k]))
  res[[length(res) + 1]] <- data.frame(
    exclusion = excl, event_set = es, condition = cond, n_subjects = length(ok),
    cg_z = fe["cg:z_log_ahi"], cg_z_p = pv("cg:z_log_ahi"), sg_z = fe["sg:z_log_ahi"], sg_z_p = pv("sg:z_log_ahi"),
    cgc_z = fe["cgc:z_log_ahi"], cgc_z_p = pv("cgc:z_log_ahi"), sgc_z = fe["sgc:z_log_ahi"], sgc_z_p = pv("sgc:z_log_ahi"),
    cgc = fe["cgc"], cgc_p = pv("cgc"), sgc = fe["sgc"], sgc_p = pv("sgc"), singular = isSingular(fit), row.names = NULL)
}
out <- do.call(rbind, res)
write.csv(out, args[2], row.names = FALSE)
print(format(out, digits = 3), row.names = FALSE)
