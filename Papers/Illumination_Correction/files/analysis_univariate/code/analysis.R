#------ load libraries
rm(list=ls())
library(ggplot2)
library(plyr)
library(doMC)
library(yaml)
library(stringr)
library(testthat)
library(reshape2)
library(MASS)
library(grid)
library(gridExtra)
library(knitr)

registerDoMC()
source("profile.data.x.R")
source("zprime.R")

#------ setup
test_type <- 'zprime_os'
THRESH_SCORE <- -1

#------ functions
shorten_feats <- function(l) {
  
  replist <- 
    list(c('_Intensity_', '_'), 
         c('_CorrTub', ''),
         c('Cells_', 'Ce_'),
         c('Cytoplasm_', 'Cy_'),
         c('Nuclei_', 'Nu_'),
         c('Intensity', ''),
         c('Integrated', 'Tot'),
         c('Upper', 'Up'),
         c('Lower', 'Lo'),
         c('Quartile', 'Quar'),
         c('Displacement', 'Disp'))
  
  for (rl in replist) {
    l <- lapply(l, function(x) str_replace(x, rl[1], rl[2]))
  }
  l <- unlist(l)
}

vecdiff <- function(P, func=f_zprime, nbin=5, feat.i=NULL) {
  
  X <- feats(P)
  fac <- factors(P)
  
  Xneg <- X[fac$namecpd == "DMSO",]
  Xpos <- X[fac$namecpd == "Taxol",]
  vec <- foreach(i = seq(NCOL(X)), .combine=c) %do% func(Xneg[,i], Xpos[,i])
  names(vec) <- names(X)
  
  vec
}

color_feats <- function(feat.names) {
  unlist(lapply(feat.names, function(x) {
    s = str_split(x, "_")[[1]][1]
    return(s)
  }))
}

get_week <- function(P) 
  unlist(lapply(factors(P)$Plate, 
                function(x) strsplit(as.character(x), "_")[[1]][[1]]))

trim_comp <- function(feat.names) 
  unlist(lapply(feat.names, function(x) str_split(x, "_", 2)[[1]][2] ))  


# List of plates that failed quality checks based on cell counts
badplates <- c('Week10_40111','Week10_40119','Week2_24381',
               'Week6_31681','Week6_32061',
               'Week7_34381','Week7_34641',
               'Week7_34661','Week7_34681',
               'Week8_38221', 'Week9_39222',
               'Week9_39283','Week9_39301', 
               'Week6_32161')


#------ Load data
for (tag in c("", "_site", "_col", "_row")) { 
  
  # Without correction
  cf0 <- '../input/AZ_posneg_WithoutICFs/well-summary-profile_mean-median-robust_linear-ctrl_norm.yml'
  P0 <- profile.data.x(cf0) 
  
  # With correction
  cf1 <- sprintf('../input/AZ_posneg_WithICFs_Median500M%s/well-summary-profile_mean-median-robust_linear-ctrl_norm.yml', tag)
  P1 <- profile.data.x(cf1) 
  # Blank tag defaults to plate
  if (tag=="") tag <- "_plate" 
  
  P0$data$week <- get_week(P0)
  P1$data$week <- get_week(P1)
  P0$factor_cols <- c(P0$factor_cols, "week")
  P1$factor_cols <- c(P1$factor_cols, "week")
  
  P0$data <- droplevels(subset(P0$data, !(Plate %in% badplates)))
  P1$data <- droplevels(subset(P1$data, !(Plate %in% badplates)))
  P0$data <- P0$data[with(P0$data, order(Plate, Well)), ]
  P1$data <- P1$data[with(P1$data, order(Plate, Well)), ]
  sel <- paste(P1$data$Plate, P1$data$Well, sep="_") %in% 
    paste(P0$data$Plate, P0$data$Well, sep="_")
  P1$data <- P1$data[sel,]
  
  feat.names <- names(feats(P0))
  feat.names <- feat.names[setdiff(
    intersect( grep('_Intensity_', feat.names), grep('CorrTub', feat.names) ),
    grep('Location', feat.names))]
  P0 <- prune.feats(P0, feat.names, keep=T)
  P1 <- prune.feats(P1, feat.names, keep=T)
  expect_true(all(names(feats(P0)) == names(feats(P1))))
  
  #------ Shorten feature names
  P0 <- rename.feats(P0, shorten_feats(P0$feat_cols))
  P1 <- rename.feats(P1, shorten_feats(P1$feat_cols))
  expect_true(all(names(feats(P0)) == names(feats(P1))))
  
  P0$data$namecpd <- mapvalues(P0$data$namecpd, 
                               c("DMSO", "taxol"), c("DMSO", "Taxol"))
  P1$data$namecpd <- mapvalues(P1$data$namecpd, 
                               c("DMSO", "taxol"), c("DMSO", "Taxol"))
  
  #------ Setup_variables
  X0 <- feats(P0)
  X1 <- feats(P1)
  y0 <- factors(P0)$namecpd
  y1 <- factors(P1)$namecpd
  expect_true(all(y0==y1))
  
  #------ Compute differences
  v0 <- vecdiff(P0, f_zprime_os)
  v1 <- vecdiff(P1, f_zprime_os)
  expect_true(all(names(v0)==names(v1)))
  
  #------ Combine data
  Dp <- data.frame(idx=seq(length(v0)), cols=color_feats(names(v0)), 
                   fnames=trim_comp(names(v0)), v0, v1)
  
  #------ Save data for generating plots that require all groupings
  Dpx <- Dp
  Dpx$tag <- tag
  Dpx$test_type <- test_type
  write.csv(Dpx, sprintf("Dp_%s%s.csv", test_type, tag), row.names=F, quote=F)
  
  #------ Create scatter plot of improvements in z'-factor
  
  Dps <- subset(Dp, (v0 > THRESH_SCORE))
  Dps$cols <- mapvalues(Dps$cols, c("Ce", "Cy", "Nu"), 
                        c("Cell", "Cytoplasm", "Nucleus") )
  mx <- 0.25
  p <- ggplot(Dps, aes(v0, v1, color=cols)) + geom_point() + geom_abline()
  p <- p + xlim(THRESH_SCORE, mx) +  ylim(THRESH_SCORE, mx)
  p <- p + ggtitle(sprintf("scatter_1_%s%s", test_type, tag))
  p <- p + xlab("NoIC") + ylab("WithIC")
  p <- p + theme_bw()
  p <- p + guides(color=guide_legend(title=NULL, nrow=1))
  p <- p + theme(legend.position = "bottom")
  ggsave(sprintf("scatter_1_%s%s.pdf", test_type, tag), p, width=6, height=6.5)
  print(p)

}



#------ plot_score_Ce_Tot
rm(list=ls())
fnames <- list.files(path = '.',
                     pattern = 'Dp_zprime_os', 
                     full.names = T) 
D <- do.call("rbind", lapply(fnames, read.csv, header = TRUE)) 
D0 <- subset(D, cols == "Ce" & fnames == "Tot")[1, c("tag", "v0")]
D0$tag <- "Uncorrected"
names(D0)[names(D0)=="v0"] <- "v1"
D <- subset(D, cols == "Ce" & fnames == "Tot")[, c("tag", "v1")]
D <- rbind(D0, D)
D$tag <- mapvalues(D$tag, 
                   c("_col", "_row", "_site", "_plate"),
                   c("Column-wise", "Row-wise", "Site-wise", "Plate-wise"))
names(D)[names(D)=="tag"] <- "Grouping"
names(D)[names(D)=="v1"] <- "Z"
D
D$Grouping <- factor(D$Grouping, levels=c("Plate-wise", 
                                          "Row-wise", "Column-wise", 
                                          "Site-wise",
                                          "Uncorrected"))
p <- ggplot(D, aes(Grouping, Z)) + 
  geom_bar(stat="identity", fill="white", color="black")
p <- p + geom_text(aes(label=sprintf("%.2f", Z)), hjust=1.2)
p <- p + ylab("Z\'-factor") + ylim(-.65, 0)
p <- p + theme_bw() + coord_flip()
p
ggsave("Ce_Tot_zprime_os.pdf", p, width=6, height=6)

#------ plot_score_Ce_Tot
rm(list=ls())
fnames <- list.files(path = '.',
                     pattern = 'Dp_zprime_os', 
                     full.names = T) 
D <- do.call("rbind", lapply(fnames, read.csv, header = TRUE)) 
D$test_type <- NULL
names(D) <- c("Feature_id", "Cell_Compartment", "Feature", "Zprime_without_correction", "Zprime_after_correction", "Grouping")
D$Cell_Compartment <- mapvalues(D$Cell_Compartment, 
                                from = c("Ce", "Cy", "Nu"), 
                                to = c("Cell", "Cytoplasm", "Nucleus"))
D$Grouping <- mapvalues(D$Grouping, 
                                from = c("_col", "_row", "_site", "_plate"), 
                                to = c("column-wise", "row-wise", "site-wise", "plate-wise"))
write.csv(D, "Intensity_features_Zprime_improvement.csv", row.names=F, quote=F)
