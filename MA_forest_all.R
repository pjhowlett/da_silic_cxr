library("meta")
library("readxl")
library("patchwork")
library("gridExtra")
library("grid")


# Set wd 
setwd("INSERT YOUR WORKING DIRECTORY HERE")

# Load files from Documents

HRCT_Normal <- read_xlsx(path = "HRCT_all.xlsx")
CT_Normal <- read_xlsx(path = "CT_all.xlsx")
Autopsy <- read_xlsx(path = "Autopsy.xlsx")

# Check the class of difference variables

HRCT_Normal$sens_tot <- HRCT_Normal$FN + HRCT_Normal$TP

# Metaprop analysis

# Fix this code so it doesnt run an error: HRCT_Sensitivity <- metaprop(HRCT_Normal$TP, $sens_tot, comb.fixed=FALSE,comb.random = TRUE, sm = "plogit", method = "GLMM", studlab=HRCT_Normal$ID) 

HRCT_Sensitivity <- metaprop(HRCT_Normal$TP, HRCT_Normal$sens_tot, comb.fixed=FALSE, comb.random = TRUE, sm = "plogit", studlab=HRCT_Normal$ID)

HRCT_Specificity <- metaprop(HRCT_Normal$TN, HRCT_Normal$TN + HRCT_Normal$FP, comb.fixed=FALSE, comb.random = TRUE, sm = "plogit", studlab=HRCT_Normal$ID)

CT_Sensitivity <- metaprop(CT_Normal$TP, CT_Normal$FN + CT_Normal$TP, comb.fixed=FALSE, comb.random = TRUE, sm = "plogit", studlab=CT_Normal$ID)
CT_Specificity <- metaprop(CT_Normal$TN, CT_Normal$TN + CT_Normal$FP, comb.fixed=FALSE, comb.random = TRUE, sm = "plogit", studlab=CT_Normal$ID)

Autopsy_Sensitivity <- metaprop(Autopsy$TP, Autopsy$FN + Autopsy$TP, comb.fixed=FALSE, comb.random = TRUE, sm = "plogit", method.ci = "CP", studlab=Autopsy$ID)
Autopsy_Specificity <- metaprop(Autopsy$TN, Autopsy$TN + Autopsy$FP, comb.fixed=FALSE, comb.random = TRUE, sm = "plogit", method.ci = "CP", studlab=Autopsy$ID)

# Forest plots (for higher resolution image, export differently)

png(filename = "HRCT_Sensitivity_1_all.png", width = 3600, height = 2262, res = 400)
forest(HRCT_Sensitivity, xlim = c(0,1), pscale =1, digits = 3, fs.xlab=12, fs.study=12, fs.study.lables=12, fs.heading=12, squaresize = 0.5, col.square="navy", 
       col.square.lines="navy", col.diamond="navy", col.diamond.lines="navy", type.random="diamond", ff.fixed="bold.italic", 
       ff.random="bold.italic", print.Q=TRUE, print.pval.Q=TRUE, print.I2=TRUE, print.tau2=FALSE, 
       lty.fixed=0, lty.random=2, col.by="black", comb.fixed=FALSE, xlab="Sensitivity", lwd = 2, cex= 1.2, 
       rightcols=c("effect", "ci"), rightlabs=c("Sensitivity", "95% C.I."), leftcols = c("studlab", "event", "n"), 
       leftlabs = c("Study", "True Positive", "Total"))
dev.off()


png(filename = "HRCT_Specificity_1_all.png", width = 3600, height = 2262, res = 400)
forest(HRCT_Specificity, xlim = c(0,1), pscale =1, digits = 3, fs.xlab=12, fs.study=12, 
       fs.study.lables=12, fs.heading=12, squaresize = 0.5, col.square="navy", 
       col.square.lines="navy", col.diamond="navy", col.diamond.lines="navy", 
       type.random="diamond", ff.fixed="bold.italic", ff.random="bold.italic", 
       print.Q=TRUE, print.pval.Q=TRUE, print.I2=TRUE, print.tau2=FALSE, lty.fixed=0, 
       lty.random=2, col.by="black", comb.fixed=FALSE, xlab="Specificity", lwd = 2, 
       cex= 1.2, rightcols=c("effect", "ci"), rightlabs=c("Specificity", "95% C.I."), 
       leftcols = c("studlab", "event", "n"), leftlabs = c("Study", "True Negative", "Total"))
dev.off()

png(filename = "CT_Sensitivity_1_all.png",  width = 3600, height = 2262, res = 400)
forest(CT_Sensitivity, xlim = c(0,1), pscale =1, digits = 3, fs.xlab=12, fs.study=12, 
       fs.study.lables=12, fs.heading=12, squaresize = 0.5, col.square="navy", 
       col.square.lines="navy", col.diamond="navy", col.diamond.lines="navy", 
       type.random="diamond", ff.fixed="bold.italic", ff.random="bold.italic", 
       print.Q=TRUE, print.pval.Q=TRUE, print.I2=TRUE, print.tau2=FALSE, lty.fixed=0, 
       lty.random=2, col.by="black", comb.fixed=FALSE, xlab="Sensitivity", lwd = 2, 
       cex= 1.2, rightcols=c("effect", "ci"), rightlabs=c("Sensitivity", "95% C.I."), 
       leftcols = c("studlab", "event", "n"), leftlabs = c("Study", "True Positive", "Total"))
dev.off()

png(filename = "CT_Specificity_1_all.png",  width = 3600, height = 2262, res = 400)
forest(CT_Specificity, xlim = c(0,1), pscale =1, digits = 3, fs.xlab=12, fs.study=12, 
       fs.study.lables=12, fs.heading=12, squaresize = 0.5, col.square="navy", 
       col.square.lines="navy", col.diamond="navy", col.diamond.lines="navy", 
       type.random="diamond", ff.fixed="bold.italic", ff.random="bold.italic", 
       print.Q=TRUE, print.pval.Q=TRUE, print.I2=TRUE, print.tau2=FALSE, lty.fixed=0, 
       lty.random=2, col.by="black", comb.fixed=FALSE, xlab="Specificity", lwd = 2, 
      cex= 1.2, rightcols=c("effect", "ci"), rightlabs=c("Specificity", "95% C.I."), leftcols = c("studlab", "event", "n"), leftlabs = c("Study", "True Negative", "Total"))
dev.off()

png(filename = "Autopsy_Sensitivity_1_all.png",  width = 3600, height = 2262, res = 400)
forest(Autopsy_Sensitivity, xlim = c(0,1), pscale =1, digits = 3, fs.xlab=12, fs.study=12, fs.study.lables=12, fs.heading=12, squaresize = 0.5, col.square="navy", col.square.lines="navy", col.diamond="navy", col.diamond.lines="navy", type.random="diamond", ff.fixed="bold.italic", ff.random="bold.italic", print.Q=TRUE, print.pval.Q=TRUE, print.I2=TRUE, print.tau2=FALSE, lty.fixed=0, lty.random=2, col.by="black", comb.fixed=FALSE, xlab="Sensitivity", lwd = 2, main = "Forest Plot for Autopsy Sensitivity", cex= 1.2, rightcols=c("effect", "ci"), rightlabs=c("Sensitivity", "95% C.I."), leftcols = c("studlab", "event", "n"), leftlabs = c("Study", "True Positive", "Total"))
dev.off()

png(filename = "Autopsy_Specificity_1_all.png",  width = 3600, height = 2262, res = 400)
forest(Autopsy_Specificity, xlim = c(0,1), pscale =1, digits = 3, fs.xlab=12, fs.study=12, fs.study.lables=12, fs.heading=12, squaresize = 0.5, col.square="navy", col.square.lines="navy", col.diamond="navy", col.diamond.lines="navy", type.random="diamond", ff.fixed="bold.italic", ff.random="bold.italic", print.Q=TRUE, print.pval.Q=TRUE, print.I2=TRUE, print.tau2=FALSE, lty.fixed=0, lty.random=2, col.by="black", comb.fixed=FALSE, xlab="Specificity", lwd = 2, main = "Forest Plot for Autopsy Specificity", cex= 1.2, rightcols=c("effect", "ci"), rightlabs=c("Specificity", "95% C.I."), leftcols = c("studlab", "event", "n"), leftlabs = c("Study", "True Negative", "Total"))
dev.off()

# Load the saved plots
plot_hrsens <- png::readPNG("HRCT_Sensitivity_1_all.png")
plot_hrspec <- png::readPNG("HRCT_Specificity_1_all.png")
plot_ctsens <- png::readPNG("CT_Sensitivity_1_all.png")
plot_ctspec <- png::readPNG("CT_Specificity_1_all.png")
plot_autsens <- png::readPNG("Autopsy_Sensitivity_1_all.png")
plot_autspec <- png::readPNG("Autopsy_Specificity_1_all.png")

# Arrange the plots side by side
hrct_plot <- grid.arrange(
  rasterGrob(plot_hrsens),
  rasterGrob(plot_hrspec),
  nrow = 1,
  ncol = 2
)

# Save the combined plot
png(filename = "HRCT_Sensitivity_and_Specificity_all.png", width = 7200, height = 2262, res = 400)
grid.draw(hrct_plot)
dev.off()

# Arrange the plots side by side
ct_plot <- grid.arrange(
  rasterGrob(plot_ctsens),
  rasterGrob(plot_ctspec),
  nrow = 1,
  ncol = 2
)

# Save the combined plot
png(filename = "CT_Sensitivity_and_Specificity_all.png", width = 7200, height = 2262, res = 400)
grid.draw(ct_plot)
dev.off()

# Arrange the plots side by side
aut_plot <- grid.arrange(
  rasterGrob(plot_autsens),
  rasterGrob(plot_autspec),
  nrow = 1,
  ncol = 2
)

# Save the combined plot
png(filename = "Autopsy_Sensitivity_and_Specificity_all.png", width = 7200, height = 2262, res = 400)
grid.draw(aut_plot)
dev.off()


