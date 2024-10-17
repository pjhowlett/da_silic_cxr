library(meta)
library(metafor)

dpi <- 450
width_inch <- 3600
height_inch <- 2200

setwd("***INSERT YOUR WORKING DIRECTORY HERE***")

# create forest plot
create_forest_plot <- function(meta_obj, title, filename, measure) {
  png(filename = paste0(filename, ".png"), width = width_inch, height = height_inch, res = dpi)
  forest(meta_obj, 
         xlim = c(0,1), 
         pscale = 1, 
         digits = 3, 
         fs.xlab = 12, 
         fs.study = 12, 
         fs.study.labels = 12, 
         fs.heading = 12, 
         squaresize = 0.5, 
         col.square = "navy", 
         col.square.lines = "navy", 
         col.diamond = "navy", 
         col.diamond.lines = "navy", 
         type.random = "diamond", 
         ff.fixed = "bold.italic", 
         ff.random = "bold.italic", 
         print.Q = TRUE, 
         print.pval.Q = TRUE, 
         print.I2 = TRUE, 
         print.tau2 = FALSE, 
         lty.fixed = 0, 
         lty.random = 2, 
         col.by = "black", 
         comb.fixed = FALSE, 
         xlab = measure, 
         lwd = 2, 
         main = title, 
         cex = 1.2, 
         rightcols = c("effect", "ci"), 
         rightlabs = c(measure, "95% C.I."), 
         leftcols = c("studlab", "event", "n"), 
         leftlabs = c("Study", ifelse(measure == "Sensitivity", "True Positive", "True Negative"), "Total"))
  dev.off()
}

# Data for each cut-off
data_cutoff1 <- data.frame(
  study = c("Antao 2005", "Bergin 1986", "Corbett 1999", "Cowie 1993", "Hnizdo 1993", "Lopes 2008", "Talini 1995"),
  TP = c(19, 17, 40, 46, 163, 40, 14),
  FN = c(3, 1, 40, 2, 163, 4, 6),
  FP = c(3, 0, 12, 9, 25, 0, 5),
  TN = c(16, 5, 149, 13, 206, 1, 2),
  type = c("rad","rad","aut","aut","aut","rad","rad")
)

data_cutoff2 <- data.frame(
  study = c("Antao 2005", "Bergin 1986", "Corbett 1999", "Cowie 1993", "Hnizdo 1993", "Lopes 2008", "Talini 1995"),
  TP = c(17, 17, 26, 18, 100, 27, 7),
  FN = c(1, 0, 9, 1, 62, 0, 0),
  FP = c(5, 0, 26, 38, 88, 13, 12),
  TN = c(18, 6, 180, 14, 307, 4, 8), 
  type = c("rad","rad","aut","aut","aut","rad","rad")
)

data_cutoff3 <- data.frame(
  study = c("Antao 2005", "Bergin 1986", "Corbett 1999", "Cowie 1993", "Hnizdo 1993", "Lopes 2008", "Talini 1995"),
  TP = c(8, 12, 5, 1, 38, 6, 4),
  FN = c(0, 0, 0, 0, 8, 0, 0),
  FP = c(14, 5, 47, 54, 150, 34, 15),
  TN = c(19, 6, 189, 15, 361, 4, 8),
  type = c("rad","rad","aut","aut","aut","rad","rad")
)

# Reorder rows in all dfs so order is now: Bergin 1986, Talini 1995, Antao 2005, Lopes 2008, Hnizdo 1993, Cowie 1993, Corbett 1999
data_cutoff1 <- data_cutoff1[c(2, 7, 1, 6, 5, 4, 3),]
data_cutoff2 <- data_cutoff2[c(2, 7, 1, 6, 5, 4, 3),]
data_cutoff3 <- data_cutoff3[c(2, 7, 1, 6, 5, 4, 3),]

data_cutoff_rad1 <- data_cutoff1[data_cutoff1$type == "rad",]
data_cutoff_rad2 <- data_cutoff2[data_cutoff2$type == "rad",]
data_cutoff_rad3 <- data_cutoff3[data_cutoff3$type == "rad",]

data_cutoff_aut1 <- data_cutoff1[data_cutoff1$type == "aut",]
data_cutoff_aut2 <- data_cutoff2[data_cutoff2$type == "aut",]
data_cutoff_aut3 <- data_cutoff3[data_cutoff3$type == "aut",]


# Create forest plots for each cut-off for radiology
for (i in 1:3) {
  data <- get(paste0("data_cutoff_rad", i))
  
  # Sensitivity
  meta_sens <- metaprop(event = TP, n = TP + FN, studlab = study, data = data, sm = "plogit", comb.fixed = FALSE, comb.random=TRUE)
  create_forest_plot(meta_sens, 
                     paste("Forest Plot for Sensitivity (Cut-off", i, ")"), 
                     paste0("Sensitivity_Cutoff_rad_", i), 
                     "Sensitivity")
  
  # Specificity
  meta_spec <- metaprop(event = TN, n = TN + FP, studlab = study, data = data, sm = "plogit", comb.fixed = FALSE, comb.random=TRUE)
  create_forest_plot(meta_spec, 
                     paste("Forest Plot for Specificity (Cut-off", i, ")"), 
                     paste0("Specificity_Cutoff_rad_", i), 
                     "Specificity")
}

rad_sens1 <- png::readPNG("Sensitivity_Cutoff_rad_1.png")
rad_sens2 <- png::readPNG("Sensitivity_Cutoff_rad_2.png")
rad_sens3 <- png::readPNG("Sensitivity_Cutoff_rad_3.png")
rad_spec1 <- png::readPNG("Specificity_Cutoff_rad_1.png")
rad_spec2 <- png::readPNG("Specificity_Cutoff_rad_2.png")
rad_spec3 <- png::readPNG("Specificity_Cutoff_rad_3.png")

# Arrange the plots side by side
sens_plot <- grid.arrange(
  rasterGrob(rad_sens1),
  rasterGrob(rad_spec1),
  nrow = 1,
  ncol = 2
) 

# Save the combined plot
png(filename = "rad_sens_spec1.png", width = 7200, height = 2000, res = 600)
grid.draw(sens_plot )
dev.off()


# Arrange the plots side by side
sens_plot2 <- grid.arrange(
  rasterGrob(rad_sens2),
  rasterGrob(rad_spec2),
  nrow = 1,
  ncol = 2
) 

# Save the combined plot
png(filename = "rad_sens_spec2.png", width = 7200, height = 2000, res = 600)
grid.draw(sens_plot2 )
dev.off()


# Arrange the plots side by side
sens_plot3 <- grid.arrange(
  rasterGrob(rad_sens3),
  rasterGrob(rad_spec3),
  nrow = 1,
  ncol = 2
) 

# Save the combined plot
png(filename = "rad_sens_spec3.png", width = 7200, height = 2000, res = 600)
grid.draw(sens_plot3 )
dev.off()


# Create forest plots for each cut-off for autopsy
for (i in 1:3) {
  data <- get(paste0("data_cutoff_aut", i))
  
  # Sensitivity
  meta_sens <- metaprop(event = TP, n = TP + FN, studlab = study, data = data, sm = "plogit", comb.fixed = FALSE, comb.random=TRUE)
  create_forest_plot(meta_sens, 
                     paste("Forest Plot for Sensitivity (Cut-off", i, ")"), 
                     paste0("Sensitivity_Cutoff_aut_", i), 
                     "Sensitivity")
  
  # Specificity
  meta_spec <- metaprop(event = TN, n = TN + FP, studlab = study, data = data, sm = "plogit", comb.fixed = FALSE, comb.random=TRUE)
  create_forest_plot(meta_spec, 
                     paste("Forest Plot for Specificity (Cut-off", i, ")"), 
                     paste0("Specificity_Cutoff_aut_", i), 
                     "Specificity")
}

aut_sens1 <- png::readPNG("Sensitivity_Cutoff_aut_1.png")
aut_sens2 <- png::readPNG("Sensitivity_Cutoff_aut_2.png")
aut_sens3 <- png::readPNG("Sensitivity_Cutoff_aut_3.png")
aut_spec1 <- png::readPNG("Specificity_Cutoff_aut_1.png")
aut_spec2 <- png::readPNG("Specificity_Cutoff_aut_2.png")
aut_spec3 <- png::readPNG("Specificity_Cutoff_aut_3.png")

# Arrange the plots side by side
sens_plot <- grid.arrange(
  rasterGrob(aut_sens1),
  rasterGrob(aut_spec1),
  nrow = 1,
  ncol = 2
) 

# Save the combined plot
png(filename = "aut_sens_spec1.png", width = 7200, height = 2000, res = 600)
grid.draw(sens_plot )
dev.off()


# Arrange the plots side by side
sens_plot2 <- grid.arrange(
  rasterGrob(aut_sens2),
  rasterGrob(aut_spec2),
  nrow = 1,
  ncol = 2
) 

# Save the combined plot
png(filename = "aut_sens_spec2.png", width = 7200, height = 2000, res = 600)
grid.draw(sens_plot2 )
dev.off()


# Arrange the plots side by side
sens_plot3 <- grid.arrange(
  rasterGrob(aut_sens3),
  rasterGrob(aut_spec3),
  nrow = 1,
  ncol = 2
) 

# Save the combined plot
png(filename = "aut_sens_spec3.png", width = 7200, height = 2000, res = 600)
grid.draw(sens_plot3 )
dev.off()
