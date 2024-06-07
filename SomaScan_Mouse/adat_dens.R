## Initialize document
library(tidyverse); library(ggpubr); library(RColorBrewer); library(ggrepel); library(grid)
library(SomaDataIO); library(sm); library(scales)
setwd("~/JDJG/clean/")

# Define plotting functions
replace.na <- function (List, comp, value) {
  arg <- paste(substitute(List), "$", substitute(comp), sep = "")
  arg.value <- eval(parse(text = arg), parent.frame(1))
  if (any(is.na(arg.value))) {
    change <- paste(arg, "<-", deparse(substitute(value)))
    a <- eval(parse(text = change), parent.frame(1))
  }
  invisible()
}

density_compare_function <- function (x, group, h, model = "none",
                                      alpha = 0.7, lty = "mix", leg = group,
                                      top = 0, bot = 0, ...) 
  {
  # Set plot margins
  par(mar = c(5, 4, 4, 2) + 0.1 + c(bot,0,top,0)) # reduce top space
  # low, left, top, right
  
  if (!is.vector(x)) 
    stop("sm.density.compare can handle only 1-d data")
  opt <- sm.options(list(...))
  replace.na(opt, ngrid, 50)
  replace.na(opt, display, "line")
  replace.na(opt, xlab, deparse(substitute(x)))
  replace.na(opt, ylab, "Density")
  replace.na(opt, xlim, c(min(x) - diff(range(x))/4, max(x) + 
                            diff(range(x))/4))
  replace.na(opt, eval.points, seq(opt$xlim[1], opt$xlim[2], 
                                   length = opt$ngrid))
  replace.na(opt, col.band, "cyan")
  if (is.na(opt$band)) {
    if (model == "none") 
      opt$band <- FALSE
    else opt$band <- TRUE
  }
  if ((model == "none") && opt$band) 
    opt$band <- FALSE
  band <- opt$band
  ngrid <- opt$ngrid
  xlim <- opt$xlim
  nboot <- opt$nboot
  y <- x
  if (is.na(opt$test)) {
    if (model == "none") 
      opt$test <- FALSE
    else opt$test <- TRUE
  }
  if ((model == "none") && opt$test) 
    opt$test <- FALSE
  test <- opt$test
  if (opt$display %in% "none") 
    band <- FALSE
  fact <- factor(group)
  fact.levels <- levels(fact)
  nlev <- length(fact.levels)
  ni <- table(fact)
  if (band & (nlev > 2)) {
    cat("Reference band available to compare two groups only.\n")
    band <- FALSE
  }
  # Set lty
  if (length(opt$lty) < nlev) 
    opt$lty <- 1:nlev
  if (!lty == "mix"){
    opt$lty <- rep(lty, nlev)
  }
  # Set colour and alpha level
  if (length(opt$col) < nlev) 
    opt$col <- 2:(nlev + 1)
  opt$col <- scales::alpha(opt$col, alpha)
  if (missing(h)) 
    h <- h.select(x, y = NA, group = group, ...)
  opt$band <- band
  opt$test <- test
  estimate <- matrix(0, ncol = opt$ngrid, nrow = nlev)
  se <- matrix(0, ncol = opt$ngrid, nrow = nlev)
  for (i in 1:nlev) {
    sm <- sm.density(y[fact == fact.levels[i]], h = h, display = "none", 
                     eval.points = opt$eval.points)
    estimate[i, ] <- sm$estimate
    se[i, ] <- sm$se
  }
  eval.points <- sm$eval.points
  if (!(("none" %in% opt$display) | band)) {
    replace.na(opt, yht, 1.1 * max(as.vector(estimate)))
    replace.na(opt, ylim, c(0, opt$yht))
    plot(xlim, opt$ylim, xlab = opt$xlab, ylab = opt$ylab, 
         type = "n")
    for (i in 1:nlev) lines(eval.points, estimate[i, ], lty = opt$lty[i], 
                            col = opt$col[i], lwd = opt$lwd)
  }
  legend("topright", legend = leg, col = unique(opt$col), lty = rep(lty,3), lwd = opt$lwd)
  est <- list(estimate = estimate, eval.points = eval.points, 
              h = h, levels = fact.levels, col = opt$col, lty = opt$lty, 
              lwd = opt$lwd)
  p <- NULL
  invisible(est)
}


density_plot <- function(i = 1, norm = "unnorm", bot = 0, top = 0, title_pos = 0){
  dat <- get(paste0("adat", i, "_", norm))
  dat$SampleId <- as.character(dat$SampleId)
  for (fix_samplenames in "now"){
    dat <- dat[order(dat$SampleId),]
    samples <- dat$SampleId
    cnt <- 0; rootname <- "noname"
    for (j in 1:length(dat$SampleId)){
      # Check for identical sample names
      cur_sample <- samples[j]
      nxt_sample <- samples[j+1]
      if (is.na(nxt_sample)){
        break
      }
      if (cur_sample == nxt_sample){
        rootname <- cur_sample
      }
      # Adjust sample names
      if (cur_sample == rootname){
        cnt <- cnt + 1
        samples[j] <- paste0(rootname, "_", cnt)
      }
      
      if (!nxt_sample == rootname){
        cnt <- 0
      }
    }
    dat$SampleId <- samples
  }
  rownames(dat) <- dat$SampleId
  dat <- dat[order(dat$SampleType),]
  dat$SampleId <- factor(dat$SampleId, levels = dat$SampleId)
  sampletypes <- dat$SampleType
  cols <- rep("#619CFF", length(sampletypes))
  cols[sampletypes == "Sample"] <- "#00BA38"; cols[sampletypes == "Calibrator"] <- "#F8766D"
  dat <- reshape2::melt(dat, id.vars = c("SampleId", "SampleType"))
  dat$value <- log2(dat$value)
  if (i == 1){
    title <- "SomaMouse1"
  }
  if (i == 2){
    title <- "SomaMouse2"
  }
  
  # title <- paste0(title, "Density distribution for each sample")
  norm <- str_split(norm, "_")[[1]][1]
  if (norm == "unnorm"){
    title <- paste0(title, ", unnormalised")
  }
  if (norm == "norm1"){
    title <- paste0(title, ", partially normalised")
  }
  if (norm == "norm2"){
    title <- paste0(title, ", fully normalised")
  }
  density_compare_function(dat$value, dat$SampleId,
                           lwd = 4, lty = 1, col = cols, alpha = 0.7, leg = unique(dat$SampleType),
                           xlab = "Log2-transformed aptamer measurements (RFU)",
                           top = top, bot = bot)
  title(main = title, line = title_pos)
}
# density_plot(i=1, norm="unnorm")
# density_plot(i=2, norm="norm1")


## Load data
load("~/NNEDL/masterprojectxjdjg/curated/Mouse_studies/curated_data_mouse_2.RData")
adat1_unnorm <- read_adat("~/NNEDL/masterprojectxjdjg/raw/Mouse_studies/JYNR210602_GUS2021_390/SS-216967_v4.1_MousePlasma.adat")[,-c(1:6,9:25)]
adat1_norm1 <- read_adat("~/NNEDL/masterprojectxjdjg/raw/Mouse_studies/JYNR210602_GUS2021_390/SS-216967_v4.1_MousePlasma.hybNorm.medNormInt.plateScale.adat")[,-c(1:6,9:30)]
adat1_norm2 <- read_adat("~/NNEDL/masterprojectxjdjg/raw/Mouse_studies/JYNR210602_GUS2021_390/SS-216967_v4.1_MousePlasma.hybNorm.medNormInt.plateScale.medNormSMP.adat")[,-c(1:6,9:30)]
adat2_unnorm <- read_adat("~/NNEDL/masterprojectxjdjg/raw/Mouse_studies/JYNR220601_GUS2022_898/SS-2231349_v4.1_other.adat")[,-c(1:6,9:25)]
adat2_norm1 <- read_adat("~/NNEDL/masterprojectxjdjg/raw/Mouse_studies/JYNR220601_GUS2022_898/SS-2231349_v4.1_other.hybNorm.medNormInt.plateScale.adat")[,-c(1:6,9:30)]
adat2_norm2 <- read_adat("~/NNEDL/masterprojectxjdjg/raw/Mouse_studies/JYNR220601_GUS2022_898/SS-2231349_v4.1_other.hybNorm.medNormInt.plateScale.medNormSMP.adat")[,-c(1:6,9:30)]

# ## Plot and wrangle
# for (letsgo in "go"){
#   png(filename = "~/JDJG/figures/adat_dens.png", width = 3000, height = 1500, units = "px", pointsize = 38)
#   layout(mat = matrix(c(1,2,3,4,5,6), nrow = 2, ncol = 3, byrow=TRUE))
#   for (i in c(1,2)){
#     for (norm in c("unnorm", "norm1", "norm2")){
#       density_plot(i = i, norm = norm)
#     }
#   }
#   dev.off()
# }
# 
# 
# ## Plot and wrangle (thesis version)
# for (letsgo in "go"){
#   png(filename = "~/JDJG/figures/adat_dens_short.png", width = 3000, height = 750, units = "px", pointsize = 26)
#   layout(mat = matrix(c(1,2,3,4,5,6), nrow = 2, ncol = 3, byrow=TRUE))
#   for (i in c(1,2)){
#     for (norm in c("unnorm", "norm1", "norm2")){
#       density_plot(i = i, norm = norm)
#     }
#   }
#   dev.off()
# }


### Outlier detection and removal
## Outlier detection
# JYNR210602
dat <- adat1_norm1[,-c(1,2)]
sampleid <- adat1_norm1$SampleId
sort(sapply(1:nrow(dat), function(x) log2(median(unlist(dat[x,])))))
# Samples with median RFUs 12.403451, 11.470456, 11.376451, and 11.337036
# will be removed.
sapply(1:nrow(dat), function(x) log2(median(unlist(dat[x,]))))
# These values correspond to the samples no. 48, 44, 53, and 18.
sampleid[c(48,44,53,18)]
# Corresponding to samples 4018697918, 4018696639, 4018711275, and 4018700427.
meta_mouse_JYNR210602$Group[meta_mouse_JYNR210602$TubeID %in% c(4018697918, 4018696639, 4018711275, 4018700427)]
# 4018697918 has already been removed from the curated data
outrm1 <- c(48,44,53,18)

# JYNR220601
dat <- adat2_norm1[,-c(1,2)]
sampleid <- adat2_norm1$SampleId
sort(sapply(1:nrow(dat), function(x) log2(median(unlist(dat[x,])))))
# The sample with median RFU 12.303852
# will be removed.
sapply(1:nrow(dat), function(x) log2(median(unlist(dat[x,]))))
# This value corresponds to the sample no. 11
sampleid[c(11)]
# Corresponding to the sample 118606.
meta_mouse_JYNR220601$Group[meta_mouse_JYNR220601$`ANIMAL#` %in% c(118606)]
# 118606 has already been removed from the curated data
outrm2 <- c(11)

## Outlier removal
adat1_unnorm_outrm <- adat1_unnorm[-outrm1,]
adat1_norm1_outrm <- adat1_norm1[-outrm1,]
adat1_norm2_outrm <- adat1_norm2[-outrm1,]

adat2_unnorm_outrm <- adat2_unnorm[-outrm2,]
adat2_norm1_outrm <- adat2_norm1[-outrm2,]
adat2_norm2_outrm <- adat2_norm2[-outrm2,]


# ## Plot and wrangle after outlier removal
# for (letsgo in "go"){
#   png(filename = "~/JDJG/figures/adat_dens_outrm.png", width = 3000, height = 1500, units = "px", pointsize = 38)
#   layout(mat = matrix(c(1,2,3,4,5,6), nrow = 2, ncol = 3, byrow=TRUE))
#   for (i in c(1,2)){
#     for (norm in c("unnorm_outrm", "norm1_outrm", "norm2_outrm")){
#       density_plot(i = i, norm = norm)
#     }
#   }
#   dev.off()
# }



### Plot with samples matching the curated (final) version
## Sample matching
# JYNR210602
sampleid <- adat1_norm1$SampleId
samples <- sampleid %in% meta_mouse_JYNR210602$TubeID
QC <- !adat1_norm1$SampleType == "Sample"
keep_samples1 <- samples | QC

# JYNR220601
sampleid <- adat2_norm1$SampleId
samples <- sampleid %in% meta_mouse_JYNR220601$`ANIMAL#`
QC <- !adat2_norm1$SampleType == "Sample"
keep_samples2 <- samples | QC


## Sample removal
adat1_unnorm_cur <- adat1_unnorm[keep_samples1,]
adat1_norm1_cur <- adat1_norm1[keep_samples1,]
adat1_norm2_cur <- adat1_norm2[keep_samples1,]

adat2_unnorm_cur <- adat2_unnorm[keep_samples2,]
adat2_norm1_cur <- adat2_norm1[keep_samples2,]
adat2_norm2_cur <- adat2_norm2[keep_samples2,]


# ## Plot and wrangle after outlier removal
# for (letsgo in "go"){
#   png(filename = "figures/adat_dens_cur.png", width = 3000, height = 1500, units = "px", pointsize = 38)
#   layout(mat = matrix(c(1,2,3,4,5,6), nrow = 2, ncol = 3, byrow=TRUE))
#   for (i in c(1,2)){
#     for (norm in c("unnorm_cur", "norm1_cur", "norm2_cur")){
#       density_plot(i = i, norm = norm)
#     }
#   }
#   dev.off()
# }





