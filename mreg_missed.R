library("readxl")
library("ggplot2")
library("ggpubr")
library("Cairo")
library("ggpmisc")
library("metafor")
library(tidyverse)
library(ggpubr)

# Set wd 
setwd("***INSERT YOUR WORKING DIRECTORY HERE***")

# Load data
props <- read_excel("props.xlsx")
glimpse(props)

# Calculate variance
props$vi <- (props$sens*(1-props$sens)) /props$n

################################
# Linear meta-regression model #
################################

# For proportion of ILO ≥2 cases/ILO ≥2 cases irrespective of reference test
lm_ilo21 <- rma(sens, vi, mods = ~ ilo21, data = props)
summary(lm_ilo21)
coef21 <- coef(lm_ilo21)

# For overall proportion of ILO ≥2 cases confirmed by reference test
lm_ilo2 <- rma(sens, vi, mods = ~ prop2, data = props)
summary(lm_ilo2)
coef2 <- coef(lm_ilo2)

# Graph with confidence intervals using ggplot
a <- ggplot(props, aes(x = ilo21, y = sens, color = Reference)) +
  geom_point(aes(size = n)) +
  geom_smooth(method = "lm", se = T, color = "black") +
  labs(x = "Proportion of ILO ≥2/ILO ≥1 Cases", y = "Sensitivity of CXR") +
  theme_pubr() +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        axis.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 10, face = "bold"),
        legend.title = element_blank(),        
        legend.text = element_text(size = 10)) +
  scale_color_manual(values = c("Autopsy" = "red", "CT" = "blue", "HRCT" = "green")) +
  guides(size = guide_legend(override.aes = list(shape = 19)))
a
ggsave(plot = a, "mr_ilo21.png", h = 4, w = 6, type = "cairo-png")

# Graph with confidence intervals using ggplot
b <- ggplot(props, aes(x = prop2, y = sens, color = Reference)) +
  geom_point(aes(size = n)) +
  geom_smooth(method = "lm", se = T, color = "black") +
  labs(x = "Proportion of ILO ≥2 cases", y = "Sensitivity of CXR") +
  theme_pubr() +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        axis.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 10, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10)) +
  scale_color_manual(values = c("Autopsy" = "red", "CT" = "blue", "HRCT" = "green")) +
  guides(size = guide_legend(override.aes = list(shape = 19)))
b
ggsave(plot = b, "mr_prop2.png", h = 4, w = 6, type = "cairo-png")
b
#####################################
# Associate with exposure for ILO21 #
#####################################

# iom coefficients from Table 4.3b for ISP0-7
coef_ilo2 <- c(intercept_ilo2 = -4.39, ce_lt2_ilo2 = 0.255, se_intercept_ilo2 = 0.526, se_ce_lt2_ilo2 = 0.061)
coef_ilo1 <- c(intercept_ilo1 = -2.09, ce_lt2_ilo1 = 0.202, se_intercept_ilo1 = 0.25, se_ce_lt2_ilo1 = 0.037)

# Function to calculate the probability of ILO 1+ and ILO 2+ using the coefficients
calc_prob <- function(exposure, intercept, ce_lt2) {
  lp <- intercept + ce_lt2 * exposure
  prob <- exp(lp) / (1 + exp(lp))
  return(prob)
}

# Create exposure range for coefficients above
df <- data.frame(exposure_range = seq(0, 41, by = 0.1))

# Apply exposure ranges to calculate probabilities of ILO 1+ and ILO 2+
df$prob_ilo1plus <- sapply(df$exposure_range, calc_prob, intercept = coef_ilo1["intercept_ilo1"], ce_lt2 = coef_ilo1["ce_lt2_ilo1"])
df$prob_ilo2plus <- sapply(df$exposure_range, calc_prob, intercept = coef_ilo2["intercept_ilo2"], ce_lt2 = coef_ilo2["ce_lt2_ilo2"])

# Calculating the sensitivity of CXR for ILO2/ILO1 CXR cases
df$sens_ilo21 <- coef21["intrcpt"] + coef21["ilo21"] * (df$prob_ilo2plus / df$prob_ilo1plus)
df$sens_ilo21 <- pmin(df$sens_ilo21, 1)

# Calculate upper and lower bounds using standard errors
# For IOM Coef for >= ILO 1
df$upper_prob_ilo1plus <- sapply(df$exposure_range, calc_prob, intercept = coef_ilo1["intercept_ilo1"] + 1.96 * coef_ilo1["se_intercept_ilo1"],
                              ce_lt2 = coef_ilo1["ce_lt2_ilo1"] + 1.96 * coef_ilo1["se_ce_lt2_ilo1"])

df$lower_prob_ilo1plus <- sapply(df$exposure_range, calc_prob, intercept = coef_ilo1["intercept_ilo1"] - 1.96 * coef_ilo1["se_intercept_ilo1"],
                              ce_lt2 = coef_ilo1["ce_lt2_ilo1"] - 1.96 * coef_ilo1["se_ce_lt2_ilo1"])

# For IOM Coef for >= ILO 2
df$upper_prob_ilo2plus <- sapply(df$exposure_range, calc_prob, intercept = coef_ilo2["intercept_ilo2"] + 1.96 * coef_ilo2["se_intercept_ilo2"],
                              ce_lt2 = coef_ilo2["ce_lt2_ilo2"] + 1.96 * coef_ilo2["se_ce_lt2_ilo2"])
df$lower_prob_ilo2plus <- sapply(df$exposure_range, calc_prob, intercept = coef_ilo2["intercept_ilo2"] - 1.96 * coef_ilo2["se_intercept_ilo2"],
                              ce_lt2 = coef_ilo2["ce_lt2_ilo2"] - 1.96 * coef_ilo2["se_ce_lt2_ilo2"])

# Apply upper and lower bounds to main coefficients of Meta-regression model
df$upper_sensitivity_ilo1plus <- coef21["intrcpt"] + coef21["ilo21"] * (df$upper_prob_ilo2plus / df$upper_prob_ilo1plus)
df$lower_sensitivity_ilo1plus <- coef21["intrcpt"] + coef21["ilo21"] * (df$lower_prob_ilo2plus / df$lower_prob_ilo1plus)

df$upper_sensitivity_ilo1plus <- pmin(df$upper_sensitivity_ilo1plus, 1)
df$lower_sensitivity_ilo1plus <- pmin(df$lower_sensitivity_ilo1plus, 1)

# Transform exposure data to correct units 
df$exposure <- df$exposure_range * 1000 / 2160
head(df)

# Graph with confidence intervals using ggplot
c <- ggplot(df, aes(x = exposure, y = sens_ilo21)) +
  geom_ribbon(aes(ymin = lower_sensitivity_ilo1plus, ymax = upper_sensitivity_ilo1plus), fill = "darkgrey", alpha = 0.4) +
  geom_line(color = "black", size = 1.4) +
  scale_x_continuous(breaks = seq(0, max(df$exposure), by = 2), limits = c(0, max(df$exposure))) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
  labs(x = "Cumulative Silica Exposure (mg/m3-years)",
       y = "Sensitivity of CXR") +
  theme_pubr() +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        axis.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 10, face = "bold"))
c

# Arrange d and e with ggpubr, no titles 
fig_5 <- ggpubr::ggarrange(a, c, # list of plots
                           labels = "AUTO", # labels
                           align = "hv", # Align them both, horizontal and vertical
                           ncol = 2)  # number of rows

fig_5

# Save d and e plots as DR21_plots in 1 line code 
ggsave("Fig_5_MR.png", fig_5, width = 10, height = 6)

#############################
# Table for paper - miners #
#############################

exp_values <- c(1,2,4,6,8,10)

# Create function for all of the below with inputs of exp_values and coefs
full_df <- function(exp_range, coefs_ilo1, coefs_ilo2, coef21) {
  
df_tab <- data.frame(exposure = exp_values)
df_tab$exposure_range <- df_tab$exposure / (1000 / 2160)

# Apply exposure ranges to calculate probabilities of ILO 1+ and ILO 2+
df_tab$prob_ilo1plus <- sapply(df_tab$exposure_range, calc_prob, intercept = coef_ilo1["intercept_ilo1"], ce_lt2 = coef_ilo1["ce_lt2_ilo1"])
df_tab$prob_ilo2plus <- sapply(df_tab$exposure_range, calc_prob, intercept = coef_ilo2["intercept_ilo2"], ce_lt2 = coef_ilo2["ce_lt2_ilo2"])
df_tab$ilo2_1 <- df_tab$prob_ilo2plus / df_tab$prob_ilo1plus

# Calculating the sensitivity of CXR for ILO2/ILO1 CXR cases
df_tab$sens_ilo21 <- coef21["intrcpt"] + coef21["ilo21"] * (df_tab$prob_ilo2plus / df_tab$prob_ilo1plus)
df_tab$sens_ilo21 <- pmin(df_tab$sens_ilo21, 1)

# Calculate upper and lower bounds using standard errors
# For IOM Coef for >= ILO 1
df_tab$upper_prob_ilo1plus <- sapply(df_tab$exposure_range, calc_prob, intercept = coef_ilo1["intercept_ilo1"] + 1.96 * coef_ilo1["se_intercept_ilo1"],
                                 ce_lt2 = coef_ilo1["ce_lt2_ilo1"] + 1.96 * coef_ilo1["se_ce_lt2_ilo1"])

df_tab$lower_prob_ilo1plus <- sapply(df_tab$exposure_range, calc_prob, intercept = coef_ilo1["intercept_ilo1"] - 1.96 * coef_ilo1["se_intercept_ilo1"],
                                 ce_lt2 = coef_ilo1["ce_lt2_ilo1"] - 1.96 * coef_ilo1["se_ce_lt2_ilo1"])

# For IOM Coef for >= ILO 2
df_tab$upper_prob_ilo2plus <- sapply(df_tab$exposure_range, calc_prob, intercept = coef_ilo2["intercept_ilo2"] + 1.96 * coef_ilo2["se_intercept_ilo2"],
                                 ce_lt2 = coef_ilo2["ce_lt2_ilo2"] + 1.96 * coef_ilo2["se_ce_lt2_ilo2"])
df_tab$lower_prob_ilo2plus <- sapply(df_tab$exposure_range, calc_prob, intercept = coef_ilo2["intercept_ilo2"] - 1.96 * coef_ilo2["se_intercept_ilo2"],
                                 ce_lt2 = coef_ilo2["ce_lt2_ilo2"] - 1.96 * coef_ilo2["se_ce_lt2_ilo2"])
df_tab$ilo2_1_upper <- df_tab$upper_prob_ilo2plus / df_tab$upper_prob_ilo1plus
df_tab$ilo2_1_lower <- df_tab$lower_prob_ilo2plus / df_tab$lower_prob_ilo1plus

# Apply upper and lower bounds to main coefficients of Meta-regression model
df_tab$upper_sensitivity_ilo1plus <- coef21["intrcpt"] + coef21["ilo21"] * (df_tab$upper_prob_ilo2plus / df_tab$upper_prob_ilo1plus)
df_tab$lower_sensitivity_ilo1plus <- coef21["intrcpt"] + coef21["ilo21"] * (df_tab$lower_prob_ilo2plus / df_tab$lower_prob_ilo1plus)

df_tab$upper_sensitivity_ilo1plus <- pmin(df_tab$upper_sensitivity_ilo1plus, 1)
df_tab$lower_sensitivity_ilo1plus <- pmin(df_tab$lower_sensitivity_ilo1plus, 1)

print(df_tab)

}

df_tab <- full_df(exp_values, coef_ilo1, coef_ilo2, coef21)

# Round to 2 decimal places
df_tab <- round(df_tab, 2)

# Create new column Sensitivity (relative); should read sens_ilo21(95% CI lower_sensitivity_ilo1plus, upper_sensitivity_ilo1plus)
df_tab$sensitivity <- paste0(df_tab$sens_ilo21, " (", df_tab$lower_sensitivity_ilo1plus, ", ", df_tab$upper_sensitivity_ilo1plus, ")")

# Supplement table 
glimpse(df_tab)

# Load csv file pred_min_tab.csv from wd FOR MINERS 
pred_mine <- read_csv("pred_mine_tab.csv")
pred_mine$cases <- pred_mine$pred * 419.9
pred_mine$cases_upper <- pred_mine$ci.ub * 419.9
pred_mine$cases_lower <- pred_mine$ci.lb * 419.9
# Rename dose to exposure
pred_mine <- pred_mine %>% rename(exposure = dose)
glimpse(pred_mine)

# Join dt_tab and pred_mine by exposure 
df_tab_mine <- left_join(df_tab, pred_mine, by = "exposure")

# ADd column of sens_fixed all with value "0.76 (0.63, 0.86)"
sens_fixed <- 0.76
fixed_lower <- 0.63 
fixed_upper <- 0.86
df_tab_mine$sens_fixed <- "0.76 (0.63, 0.86)"
glimpse(df_tab_mine)
# Turn below into a function to create a df 

missed_df <- function(df_tab_mine){

# New column missed_rel cases scaled by the loss of sensitivity due to sens_ilo21 
df_tab_mine$missed_rel <- round(((df_tab_mine$cases/df_tab_mine$sens_ilo21) - df_tab_mine$cases),0)
# Lower CI
df_tab_mine$rel_lower <- round(((df_tab_mine$cases/df_tab_mine$lower_sensitivity_ilo1plus) - df_tab_mine$cases),0)
# Upper CI 
df_tab_mine$rel_upper <- round(((df_tab_mine$cases/df_tab_mine$upper_sensitivity_ilo1plus) - df_tab_mine$cases),0)

# Paste into 95% column 
df_tab_mine$relative_sens <- paste0(df_tab_mine$missed_rel, " (", df_tab_mine$rel_upper, ", ", df_tab_mine$rel_lower, ")")
glimpse(df_tab_mine)

# Number of cases total 
df_tab_mine$cases_rel <- df_tab_mine$cases/df_tab_mine$sens_ilo21

df_tab_mine$cases_fixed <- df_tab_mine$cases/sens_fixed

# New column missed_abs cases scaled by the loss of sensitivity due to fixed sens with round, 0
df_tab_mine$missed_abs <- round((df_tab_mine$cases/sens_fixed) - df_tab_mine$cases,0)
# Lower CI with round dp 0 
df_tab_mine$abs_lower <- round((df_tab_mine$cases/fixed_lower) - df_tab_mine$cases,0)
# Upper CI with round dp 0 
df_tab_mine$abs_upper <- round((df_tab_mine$cases/fixed_upper) - df_tab_mine$cases,0)

# Paste into 95% column
df_tab_mine$absolute_sens <- paste0(df_tab_mine$missed_abs, " (", df_tab_mine$abs_upper, ", ", df_tab_mine$abs_lower, ")")

# NNS 
df_tab_mine$nns_rel_x <- round(1000/df_tab_mine$missed_rel, 0)
df_tab_mine$nns_rel_ui <- round(1000/df_tab_mine$rel_upper, 0)
df_tab_mine$nns_rel_li <- round(1000/df_tab_mine$rel_lower, 0)
df_tab_mine$nns_rel <- paste0(df_tab_mine$nns_rel_x, " (", df_tab_mine$nns_rel_li, ", ", df_tab_mine$nns_rel_ui, ")")

df_tab_mine$nns_abs_x <- round(1000/df_tab_mine$missed_abs, 0)
df_tab_mine$nns_abs_ui <- round(1000/df_tab_mine$abs_upper, 0)
df_tab_mine$nns_abs_li <- round(1000/df_tab_mine$abs_lower, 0)
df_tab_mine$nns_abs <- paste0(df_tab_mine$nns_abs_x, " (", df_tab_mine$nns_abs_li, ", ", df_tab_mine$nns_abs_ui, ")")

# Difference 
df_tab_mine$diff <- df_tab_mine$missed_rel - df_tab_mine$missed_abs
df_tab_mine
}
df_tab_mine

# Apply function to df_tab_mine
df_tab_mine <- missed_df(df_tab_mine)

# Extract table for paper 
tab_mine_pub <- df_tab_mine %>% select(exposure, cases, sensitivity, relative_sens, nns_rel, sens_fixed, absolute_sens, nns_abs)
tab_mine_pub$cases <- round(tab_mine_pub$cases, 0)



tab_mine_pub 
# Write to csv
write_csv(tab_mine_pub, "tab_mine_pub.csv")

# Example table to show working 
glimpse(df_tab_mine)
df_tab_mine$ilo1_full <- paste0(df_tab_mine$prob_ilo1plus, " (", df_tab_mine$lower_prob_ilo1plus, ", ", df_tab_mine$upper_prob_ilo1plus, ")")
df_tab_mine$ilo2_full <- paste0(df_tab_mine$prob_ilo2plus, " (", df_tab_mine$lower_prob_ilo2plus, ", ", df_tab_mine$upper_prob_ilo2plus, ")")
df_tab_mine$ilo2_1_full <- paste0(df_tab_mine$ilo2_1, " (", df_tab_mine$ilo2_1_lower, ", ", df_tab_mine$ilo2_1_upper, ")")
example_tab <- df_tab_mine %>% 
  select(exposure, exposure_range,ilo1_full, ilo2_full,  ilo2_1_full, sensitivity)
glimpse(example_tab)

# Write to csv
write_csv(example_tab, "example_tab.csv")

# Dose-response curve 
# Load csv pred_mine.csv from wd
pred_mine_dr <- read_csv("pred_mine.csv")
pred_mine_dr$cases <- pred_mine_dr$pred * 419.9
# Rename dose to exposure
pred_mine_dr <- pred_mine_dr %>% rename(exposure = dose)

# Create df for dose response curve 
exp_values <- seq(0, 20, by = 0.01)
df_dr <- full_df(exp_values, coef_ilo1, coef_ilo2, coef21)

glimpse(df_dr)
# Join df_dr and pred_mine_dr by exposure
df_dr_mine <- left_join(df_dr, pred_mine_dr, by = "exposure")

# Apply missed_df function to df_dr_mine
df_dr_mine <- missed_df(df_dr_mine)

df_dr_mine$sens_fixed <- "0.76 (0.63, 0.86)"

# Data for plot
df_dr_mine_plot <- df_dr_mine %>% select(exposure, cases, cases_rel, cases_fixed)
glimpse(df_dr_mine_plot)

# Make pmin of 1000 for cases, cases_rel, cases_fixed use tidyverse to apply to all 3 cols in one piece of code
df_dr_mine_plot <- df_dr_mine_plot %>% mutate_all(~pmin(1000, .))

# GGplot of dose-response curve df_dr_mine_plot with line for each cases, cases_rel, cases_fixed
d <- ggplot(df_dr_mine_plot, aes(x = exposure)) +
  geom_smooth(aes(y = cases, colour = "cases"), method = "loess", span = 0.1) +  # Adjust span for smoother curve
  geom_smooth(aes(y = cases_rel, colour = "cases_rel"), method = "loess", span = 0.1) +
  geom_smooth(aes(y = cases_fixed, colour = "cases_fixed"), method = "loess", span = 0.1) +
  scale_x_continuous(breaks = seq(0, max(df$exposure), by = 2), limits = c(0, 11)) +
  labs(subtitle = "Miners",
       x = "Cumulalative silica exposure (mg/m^3)",
       y = "Cases per 1000") +
  scale_colour_manual(values = c("cases" = "black", "cases_rel" = "red", "cases_fixed" = "blue"), 
                      labels = c("cases"= "Unadjusted", "cases_rel" ="Relative sensitivity", "cases_fixed" = "Fixed sensitivity"), 
                      name = NULL) +
  theme_pubr()

d
################################
# Rpt for NONMINERS with ILO21 #
################################

# Load csv file pred_non_mine_tab.csv from wd FOR NON-MINERS
pred_non_mine <- read_csv("pred_non_mine_tab.csv")
pred_non_mine$cases <- pred_non_mine$pred * 50.9
pred_non_mine$cases_upper <- pred_non_mine$ci.ub * 50.9
pred_non_mine$cases_lower <- pred_non_mine$ci.lb * 50.9
# Rename dose to exposure
pred_non_mine <- pred_non_mine %>% rename(exposure = dose)
glimpse(pred_non_mine)

# Left join df_dr and pred_non_mine by exposure
df_tab_non_mine <- left_join(df_tab, pred_non_mine, by = "exposure")


# Apply missed_df function to df_dr_non_mine
df_tab_non_mine <- missed_df(df_tab_non_mine)
df_tab_non_mine$sens_fixed <- "0.76 (0.63, 0.86)"

# Extract table for paper
tab_non_mine_pub <- df_tab_non_mine %>% select(exposure, cases, sensitivity, relative_sens, nns_rel,
                                              sens_fixed, absolute_sens, nns_abs)
tab_non_mine_pub$cases <- round(tab_non_mine_pub$cases, 0)

tab_non_mine_pub

# Write to csv
write_csv(tab_non_mine_pub, "tab_non_mine_pub.csv")

# Now for the plot 

# Load csv file pred_non_mine.csv from wd
pred_nonmine_dr <- read_csv("pred_nonmine.csv")
pred_nonmine_dr$cases <- pred_nonmine_dr$pred * 50.9
# Rename dose to exposure
pred_nonmine_dr <- pred_nonmine_dr %>% rename(exposure = dose)

# Create df for dose response curve 
exp_values <- seq(0, 20, by = 0.01)

# Apply exposure ranges to calculate probabilities of ILO 1+ and ILO 2+
df_dr_nonmine <- full_df(exp_values, coef_ilo1, coef_ilo2, coef21)
glimpse(df_dr_nonmine)
glimpse(pred_nonmine_dr)
# Join df_dr_nonmine and pred_non_mine by exposure
df_dr_nonmine <- left_join(df_dr_nonmine, pred_nonmine_dr, by = "exposure")
glimpse(df_dr_nonmine)

# Apply missed_df function to df_dr_nonmine
df_dr_nonmine <- missed_df(df_dr_nonmine)

# Data for plot
df_dr_nonmine_plot <- df_dr_nonmine %>% select(exposure, cases, cases_rel, cases_fixed)
head(df_dr_nonmine_plot)

# GGplot of dose-response curve df_dr_nonmine_plot with line for each cases, cases_rel, cases_fixed
e <- ggplot(df_dr_nonmine_plot, aes(x = exposure)) +
  geom_smooth(aes(y = cases, colour = "cases", method = "loess", span = 0.1)) +
  geom_smooth(aes(y = cases_rel, colour = "cases_rel",method = "loess", span = 0.1)) +
  geom_smooth(aes(y = cases_fixed, colour = "cases_fixed", method = "loess", span = 0.1)) +
  scale_x_continuous(breaks = seq(0, max(df$exposure), by = 2), limits = c(0, 11)) +
  labs(subtitle = "Non-miners",
       x = "Cumulalative silica exposure (mg/m^3)",
       y = "Cases per 1000") +
  scale_colour_manual(values = c("cases" = "black", "cases_rel" = "red", "cases_fixed" = "blue"), 
                      labels = c("cases"= "Unadjusted", "cases_rel" ="Relative sensitivity", "cases_fixed" = "Fixed sensitivity"), 
                      name = NULL) +
  theme_pubr()
e
# Arrange d and e with ggpubr, no titles 
fig_6 <- ggpubr::ggarrange(d, e, # list of plots
                  labels = "AUTO", # labels
                  align = "hv", # Align them both, horizontal and vertical
                  ncol = 2)  # number of rows
fig_6


# Save d and e plots as DR21_plots in 1 line code 
ggsave("Fig_6_DR.png", fig_6, width = 10, height = 6)


#########################################
# Associate with exposure for prop ILO2 #
#########################################

# Create new exposure range for coefficients above

df2 <- data.frame(exposure = seq(0, 20, by = 0.01))
df2$exposure_range <- df2$exposure / (1000 / 2160)

# Apply exposure ranges to calculate probabilities of ILO 1+ and ILO 2+
df2$prob_ilo1plus <- sapply(df2$exposure_range, calc_prob, intercept = coef_ilo1["intercept_ilo1"], ce_lt2 = coef_ilo1["ce_lt2_ilo1"])
df2$prob_ilo2plus <- sapply(df2$exposure_range, calc_prob, intercept = coef_ilo2["intercept_ilo2"], ce_lt2 = coef_ilo2["ce_lt2_ilo2"])

# Calculating the sensitivity of CXR for ILO2/ILO1 CXR cases
df2$sens_ilo2 <- coef2["intrcpt"] + coef2["prop2"] * df2$prob_ilo2plus
df2$sens_ilo2 <- pmin(df2$sens_ilo2, 1)

# For IOM Coef for >= ILO 2
df2$upper_prob_ilo2plus <- sapply(df2$exposure_range, calc_prob, intercept = coef_ilo2["intercept_ilo2"] + 1.96 * coef_ilo2["se_intercept_ilo2"],
                                 ce_lt2 = coef_ilo2["ce_lt2_ilo2"] + 1.96 * coef_ilo2["se_ce_lt2_ilo2"])
df2$lower_prob_ilo2plus <- sapply(df2$exposure_range, calc_prob, intercept = coef_ilo2["intercept_ilo2"] - 1.96 * coef_ilo2["se_intercept_ilo2"],
                                 ce_lt2 = coef_ilo2["ce_lt2_ilo2"] - 1.96 * coef_ilo2["se_ce_lt2_ilo2"])

# Apply upper and lower bounds to main coefficients of Meta-regression model
df2$upper_sensitivity_ilo2plus <- coef2["intrcpt"] + coef2["prop2"] * df2$upper_prob_ilo2plus 
df2$lower_sensitivity_ilo2plus <- coef2["intrcpt"] + coef2["prop2"] * df2$lower_prob_ilo2plus 

df2$upper_sensitivity_ilo2plus <- pmin(df2$upper_sensitivity_ilo2plus, 1)
df2$lower_sensitivity_ilo2plus <- pmin(df2$lower_sensitivity_ilo2plus, 1)

# Graph with confidence intervals using ggplot
g <- ggplot(df2, aes(x = exposure, y = sens_ilo2)) +
  geom_ribbon(aes(ymin = lower_sensitivity_ilo2plus, ymax = upper_sensitivity_ilo2plus), fill = "darkgrey", alpha = 0.4) +
  geom_line(color = "black", size = 1.4) +
  scale_x_continuous(breaks = seq(0, max(df2$exposure), by = 2), limits = c(0, max(df2$exposure))) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
  labs(x = "Cumulative Silica Exposure (mg/m3-years)",
       y = "Sensitivity of CXR") +
  theme_pubr() +
  theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        axis.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 10, face = "bold"))
g

supp_8 <- ggpubr::ggarrange(b, g, # list of plots
                           labels = "AUTO", # labels
                           align = "hv", # Align them both, horizontal and vertical
                           ncol = 2)  # number of rows

supp_8

# Save d and e plots as DR21_plots in 1 line code 
ggsave("Supp_8_MR.png", supp_8, width = 10, height = 6)

# Table of missed cases

glimpse(df2)

# Extract rows with exposures of 1,2,4,6,8,10
df2_tab <- df2 %>% filter(exposure %in% c(1,2,4,6,8,10))

###############
# ILO2 + MINE #
###############

# Join dt_tab and pred_mine by exposure 
df2_tab_mine <- left_join(df2_tab, pred_mine, by = "exposure")
glimpse(df2_tab_mine)
# Round to 2 decimal places
df2_tab_mine <- round(df2_tab_mine, 2)

df2_tab_mine$sens_fixed <- "0.76 (0.63, 0.86)"

# Make new missed cases function 

missed_df2 <- function(df_tab_mine){
  
  # New column missed_rel cases scaled by the loss of sensitivity due to sens_ilo21 
  df_tab_mine$missed_rel <- round(((df_tab_mine$cases/df_tab_mine$sens_ilo2) - df_tab_mine$cases),0)
  # Lower CI
  df_tab_mine$rel_lower <- round(((df_tab_mine$cases/df_tab_mine$lower_sensitivity_ilo2plus) - df_tab_mine$cases),0)
  # Upper CI 
  df_tab_mine$rel_upper <- round(((df_tab_mine$cases/df_tab_mine$upper_sensitivity_ilo2plus) - df_tab_mine$cases),0)
  
  # Paste into 95% column 
  df_tab_mine$relative_sens <- paste0(df_tab_mine$missed_rel, " (", df_tab_mine$rel_upper, ", ", df_tab_mine$rel_lower, ")")
  glimpse(df_tab_mine)
  
  # Number of cases total 
  df_tab_mine$cases_rel <- df_tab_mine$cases/df_tab_mine$sens_ilo2
  
  df_tab_mine$cases_fixed <- df_tab_mine$cases/sens_fixed
  
  # New column missed_abs cases scaled by the loss of sensitivity due to fixed sens with round, 0
  df_tab_mine$missed_abs <- round((df_tab_mine$cases/sens_fixed) - df_tab_mine$cases,0)
  # Lower CI with round dp 0 
  df_tab_mine$abs_lower <- round((df_tab_mine$cases/fixed_lower) - df_tab_mine$cases,0)
  # Upper CI with round dp 0 
  df_tab_mine$abs_upper <- round((df_tab_mine$cases/fixed_upper) - df_tab_mine$cases,0)
  
  # Paste into 95% column
  df_tab_mine$absolute_sens <- paste0(df_tab_mine$missed_abs, " (", df_tab_mine$abs_upper, ", ", df_tab_mine$abs_lower, ")")
  glimpse(df_tab_mine)
  # Difference 
  df_tab_mine$diff <- df_tab_mine$missed_rel - df_tab_mine$missed_abs
  
  # Create new column Sensitivity (relative); should read sens_ilo2(95% CI lower_sensitivity_ilo2plus, upper_sensitivity_ilo2plus)
  df2_tab_mine$sensitivity <- paste0(df2_tab_mine$sens_ilo2, " (", df2_tab_mine$lower_sensitivity_ilo2plus, ", ", df2_tab_mine$upper_sensitivity_ilo2plus, ")")
  
  # NNS 
  df_tab_mine$nns_rel_x <- round(1000/df_tab_mine$missed_rel, 0)
  df_tab_mine$nns_rel_ui <- round(1000/df_tab_mine$rel_upper, 0)
  df_tab_mine$nns_rel_li <- round(1000/df_tab_mine$rel_lower, 0)
  df_tab_mine$nns_rel <- paste0(df_tab_mine$nns_rel_x, " (", df_tab_mine$nns_rel_li, ", ", df_tab_mine$nns_rel_ui, ")")
  
  df_tab_mine$nns_abs_x <- round(1000/df_tab_mine$missed_abs, 0)
  df_tab_mine$nns_abs_ui <- round(1000/df_tab_mine$abs_upper, 0)
  df_tab_mine$nns_abs_li <- round(1000/df_tab_mine$abs_lower, 0)
  df_tab_mine$nns_abs <- paste0(df_tab_mine$nns_abs_x, " (", df_tab_mine$nns_abs_li, ", ", df_tab_mine$nns_abs_ui, ")")
  
  
  # Round sensitivity to 2 dp
  df_tab_mine
}

df2_tab_mine <- missed_df2(df2_tab_mine)

# Create new column Sensitivity (relative); should read sens_ilo21(95% CI lower_sensitivity_ilo1plus, upper_sensitivity_ilo1plus)
df2_tab_mine$sensitivity <- paste0(df2_tab_mine$sens_ilo2, " (", df2_tab_mine$lower_sensitivity_ilo2plus, 
                             ", ", df2_tab_mine$upper_sensitivity_ilo2plus, ")")


# Select columns for table
tab2_mine_pub <- df2_tab_mine %>% select(exposure, cases, sensitivity, relative_sens, nns_rel,
                                         sens_fixed, absolute_sens, nns_abs)

# Round cases to 0 
tab2_mine_pub$cases <- round(tab2_mine_pub$cases,0)
tab2_mine_pub
# Write to csv
write.csv(tab2_mine_pub, "tab2_mine_pub.csv", row.names = FALSE)

# DR for mine and ILO2 
df2_dr_mine <- left_join(df2, pred_mine_dr, by = "exposure")

# Apply missed_df function to df_dr_nonmine
df2_dr_mine <- missed_df2(df2_dr_mine)

# Data for plot
df2_dr_mine <- df2_dr_mine %>% select(exposure, cases, cases_rel, cases_fixed)

# Make any value in the df have a pmin of 1000
df2_dr_mine <- df2_dr_mine %>% mutate_all(~pmin(1000, .))

# Remove rows with NA in cases
df2_dr_mine <- df2_dr_mine %>% filter(!is.na(cases))

# DR Plot
f <- ggplot(df2_dr_mine, aes(x = exposure)) +
  geom_smooth(aes(y = cases, colour = "cases"), method = "loess", span = 0.1) +  # method and span outside aes()
  geom_smooth(aes(y = cases_rel, colour = "cases_rel"), method = "loess", span = 0.1) + 
  geom_smooth(aes(y = cases_fixed, colour = "cases_fixed"), method = "loess", span = 0.1) +
  scale_x_continuous(breaks = seq(0, max(df2_dr_mine$exposure), by = 2), limits = c(0, 11)) +  # Corrected variable df
  labs(title = "Miners",
       x = "Cumulalative silica exposure (mg/m^3)",
       y = "Cases per 1000") +
  scale_colour_manual(values = c("cases" = "black", "cases_rel" = "red", "cases_fixed" = "blue"), 
                      labels = c("cases"= "Unadjusted", "cases_rel" ="Relative sensitivity", "cases_fixed" = "Fixed sensitivity"), 
                      name = NULL) +
  theme_pubr()
f
##################
# ILO2 + NONMINE #
##################

# Join dt_tab and pred_nonmine by exposure
df2_tab_nonmine <- left_join(df2_tab, pred_non_mine, by = "exposure")

# Round to 2 decimal places
df2_tab_nonmine <- round(df2_tab_nonmine, 2)
df2_tab_nonmine$sens_fixed <- "0.76 (0.63, 0.86)"
glimpse(df2_tab_nonmine)

# Apply missed_df function to df_dr_non_mine
df2_tab_nonmine <- missed_df2(df2_tab_nonmine)

# Create new column Sensitivity (relative); should read sens_ilo21(95% CI lower_sensitivity_ilo1plus, upper_sensitivity_ilo1plus)
df2_tab_nonmine$sensitivity <- paste0(df2_tab_nonmine$sens_ilo2, " (", df2_tab_nonmine$lower_sensitivity_ilo2plus, 
                                   ", ", df2_tab_nonmine$upper_sensitivity_ilo2plus, ")")

# Select columns for table
tab2_nonmine_pub <- df2_tab_nonmine %>% select(exposure, cases, sensitivity, relative_sens, nns_rel, 
                                               sens_fixed, absolute_sens, nns_abs)

# Round cases to 0
tab2_nonmine_pub$cases <- round(tab2_nonmine_pub$cases,0)
tab2_nonmine_pub

# Write to csv
write.csv(tab2_nonmine_pub, "tab2_nonmine_pub.csv", row.names = FALSE)

# DR for nonmine and ILO2
# Join df2 and pred_nonmine_dr by exposure
df2_dr_nonmine <- left_join(df2, pred_nonmine_dr, by = "exposure")

# Apply missed_df function to df_dr_nonmine
df2_dr_nonmine <- missed_df2(df2_dr_nonmine)

# Data for plot
df2_dr_nonmine <- df2_dr_nonmine %>% select(exposure, cases, cases_rel, cases_fixed)

# DR Plot
g <- ggplot(df2_dr_nonmine, aes(x = exposure)) +
  geom_smooth(aes(y = cases, colour = "cases", method = "loess", span = 0.1)) +
  geom_smooth(aes(y = cases_rel, colour = "cases_rel", method = "loess", span = 0.1)) + 
  geom_smooth(aes(y = cases_fixed, colour = "cases_fixed", method = "loess", span = 0.1)) +
  scale_x_continuous(breaks = seq(0, max(df$exposure), by = 2), limits = c(0, 11)) +
  labs(title = "Non-miners",
       x = "Cumulalative silica exposure (mg/m^3)",
       y = "Cases per 1000") +
  scale_colour_manual(values = c("cases" = "black", "cases_rel" = "red", "cases_fixed" = "blue"), 
                      labels = c("cases"= "Unadjusted", "cases_rel" ="Relative sensitivity", "cases_fixed" = "Fixed sensitivity"), 
                      name = NULL) +
  theme_pubr()

g

# Arrange d and e with ggpubr, no titles 
supp_9 <- ggpubr::ggarrange(f, g, # list of plots
                           labels = "AUTO", # labels
                           align = "hv", # Align them both, horizontal and vertical
                           ncol = 2)  # number of rows


# Save d and e plots as DR21_plots in 1 line code 
ggsave("Supp_9_DR.png", supp_9, width = 10, height = 6)
