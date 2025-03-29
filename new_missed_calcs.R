# libraries if not loaded in prior
library(tidyverse)
library(purrr)
library(ggpubr)

# set wd 
setwd("~/Library/CloudStorage/Dropbox/PhD/BSc 2024/Silicosis - Systematic Review/Missed cases")

##############
# Set up DFs #
##############
#Hypothetical populations (cases per 1000)
pop <- data.frame(
  prevalence = c("Low", "Medium", "High"),
  cases_cxr = c(50, 150, 300)  # 5%, 15%, 30% of 1000
)

# Proportions (ILO >=2/1 relative to >=1/0 ref-unconfirmed)
severity <- c(0.20, 0.40)  # 20% and 40%

# Fixed sensitivity from HRCT meta-analysis
fixed_sens <- 0.76
fixed_lower <- 0.63
fixed_upper <- 0.86

# Meta-regression coefficients and SEs from lm_ilo21 (ref-unconfirmed)
intercept <- 0.1595
slope <- 1.1468
se_intercept <- 0.1991
se_slope <- 0.3622

# Calculate missed cases & NNS
calc_metrics <- function(cases_cxr, sens) {
  true_cases <- cases_cxr / sens
  missed <- round(true_cases - cases_cxr, 0)
  nns <- ifelse(missed == 0, Inf, round(1000 / missed, 0))
  return(list(missed = missed, nns = nns))
}

# Create results data frame
results <- expand.grid(
  prevalence = pop$prevalence,
  severity = severity,
  scenario = c("Fixed", "Relative")
) %>%
  left_join(pop, by = "prevalence") %>%
  mutate(
    ilo21_prop = severity,
    sens = case_when(
      scenario == "Fixed" ~ fixed_sens,
      scenario == "Relative" ~ pmin(intercept + slope * ilo21_prop, 1)
    ),
    sens_se = case_when(
      scenario == "Fixed" ~ NA_real_,
      scenario == "Relative" ~ sqrt(se_intercept^2 + (ilo21_prop * se_slope)^2)
    ),
    sens_lower = case_when(
      scenario == "Fixed" ~ fixed_lower,
      scenario == "Relative" ~ pmax(intercept + slope * ilo21_prop - 1.96 * sens_se, 0)
    ),
    sens_upper = case_when(
      scenario == "Fixed" ~ fixed_upper,
      scenario == "Relative" ~ pmin(intercept + slope * ilo21_prop + 1.96 * sens_se, 1)
    )
  ) %>%
  mutate(
    metrics = pmap(list(cases_cxr = cases_cxr, sens = sens), calc_metrics),
    missed = map_dbl(metrics, "missed"),
    nns = map_dbl(metrics, "nns"),
    missed_lower = pmap_dbl(list(cases_cxr, sens_upper), ~ round(..1 / ..2 - ..1, 0)),
    missed_upper = pmap_dbl(list(cases_cxr, sens_lower), ~ round(..1 / ..2 - ..1, 0)),
    nns_lower = pmap_dbl(list(missed_upper), ~ ifelse(..1 == 0, Inf, round(1000 / ..1, 0))),
    nns_upper = pmap_dbl(list(missed_lower), ~ ifelse(..1 == 0, Inf, round(1000 / ..1, 0)))
  ) %>%
  mutate(
    sens_ci = sprintf("%.2f (%.2f-%.2f)", sens, sens_lower, sens_upper),  # Hyphen instead of en-dash
    missed_ci = pmap_chr(list(missed, missed_lower, missed_upper), ~ {
      m <- ifelse(is.finite(..1), sprintf("%d", ..1), "Inf")
      l <- ifelse(is.finite(..2), sprintf("%d", ..2), "Inf")
      u <- ifelse(is.finite(..3), sprintf("%d", ..3), "Inf")
      paste0(m, " (", l, "-", u, ")")
    }),
    nns_ci = pmap_chr(list(nns, nns_lower, nns_upper), ~ {
      n <- ifelse(is.finite(..1), sprintf("%d", ..1), "Inf")
      l <- ifelse(is.finite(..2), sprintf("%d", ..2), "Inf")
      u <- ifelse(is.finite(..3), sprintf("%d", ..3), "Inf")
      paste0(n, " (", l, "-", u, ")")
    })
  )

# First table improved formatting
table_output <- results %>%
  mutate(severity_label = ifelse(severity == 0.20, "Less severe (20%)", "More severe (40%)")) %>%
  select(prevalence, scenario, severity_label, sens_ci, missed_ci, nns_ci) %>%
  pivot_wider(
    names_from = scenario,
    values_from = c(sens_ci, missed_ci, nns_ci),
    names_glue = "{scenario}_{.value}"  # No prefix here, we'll rename explicitly
  ) %>%
  rename(
    "Prevalence" = prevalence,
    "Severity" = severity_label,
    "CXR Fixed Sensitivity (95% CI)" = Fixed_sens_ci,
    "CXR Relative Sensitivity (95% CI)" = Relative_sens_ci,
    "CXR Fixed Missed Cases (95% CI)" = Fixed_missed_ci,
    "CXR Relative Missed Cases (95% CI)" = Relative_missed_ci,
    "CXR Fixed NNS (95% CI)" = Fixed_nns_ci,
    "CXR Relative NNS (95% CI)" = Relative_nns_ci
  ) %>%
  arrange(Prevalence, Severity)

# Second table sameres
table_publication <- table_output %>%
  mutate(
    Prevalence = case_when(
      Prevalence == "Low" ~ "Low (5%, 50/1000)",
      Prevalence == "Medium" ~ "Medium (15%, 150/1000)",
      Prevalence == "High" ~ "High (30%, 300/1000)"
    )
  ) %>%
  arrange(factor(Prevalence, levels = c("Low (5%, 50/1000)", "Medium (15%, 150/1000)", "High (30%, 300/1000)")))

# Print sample
print(table_output)
print(table_publication)

# CSV
write_csv(table_output, "table_analysis_improved.csv")
write_csv(table_publication, "table_analysis_2.csv")

# Create continuous prevalence values
prevalence_seq <- seq(1, 40, by = 0.5)  # 1% to 40% by 0.5% increments

# Expand grid with severity and scenario
continuous_results <- expand.grid(
  prevalence = prevalence_seq,
  severity = c(0.20, 0.40),  # Corrected to match your original severity values (20% and 40%)
  scenario = c("Fixed", "Relative")
) %>%
  mutate(
    cases_cxr = prevalence * 10,  # prevalence (%) * 10 to get cases per 1000
    ilo21_prop = severity,
    sens = case_when(
      scenario == "Fixed" ~ fixed_sens,
      scenario == "Relative" ~ pmin(intercept + slope * ilo21_prop, 1)
    )
  ) %>%
  mutate(
    metrics = pmap(list(cases_cxr = cases_cxr, sens = sens), calc_metrics),
    missed = map_dbl(metrics, "missed")
  ) %>%
  mutate(
    severity_label = case_when(
      scenario == "Fixed" ~ "Fixed Sensitivity",
      scenario == "Relative" & severity == 0.20 ~ "Relative - Less severe (20%)",
      scenario == "Relative" & severity == 0.40 ~ "Relative - More severe (40%)"
    )
  )

# Plot
ggplot(continuous_results, aes(x = prevalence, y = missed, color = severity_label)) +
  geom_smooth(se = FALSE, method = "loess", span = 0.8, size = 0.9) +
  labs(
    x = "ILO ≥1/0 Silicosis Prevalence (%)",
    y = "Missed Cases per 1000 persons",
    colour = "Scenario"
  ) +
  theme_pubr(base_size = 14) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  scale_color_brewer(palette = "Set1")


####################
# ILO 2/1 over popn #
####################

pop <- data.frame(
  cases_cxr = c(50, 150, 300)  # 5%, 15%, 30% of 1000
)

# Meta-regression coefficients and SEs from lm_ilo2 (ref-unconfirmed)
intercept <- 0.4106803
slope <- 1.0392226
se_intercept <- 0.1431
se_slope <- 0.0559

pop$fixed <- 0.76
(pop$cases_cxr/1000)*.4
pop$rel_sens_low <- intercept + slope * (0.20*(pop$cases_cxr/1000))
pop$rel_missed_low <- round(pop$cases_cxr/pop$rel_sens_low,0) - pop$cases_cxr
pop$rel_nns_low <-  round(1000/pop$rel_missed_low,0)
pop$rel_sens_high <- intercept + slope * (0.40*(pop$cases_cxr/1000))
pop$rel_missed_high <- round(pop$cases_cxr/pop$rel_sens_high,0) - pop$cases_cxr
pop$rel_nns_high <-  round(1000/pop$rel_missed_high,0)
pop$fixed_missed <- round(pop$cases_cxr/pop$fixed,0) - pop$cases_cxr
pop$fixed_nns <-  round(1000/pop$fixed_missed,0)
pop

# Assuming cases_cxr represents the number of cases per 1000 (which is a prevalence rate)
pop <- data.frame(cases_cxr = seq(0, 400, by = 1))  # cases per 1000, adjusted from your example

fixed <- 0.76
rel_low <- pop
rel_low$sens <- pmin(intercept + slope * (0.20 * (rel_low$cases_cxr / 1000)),1)
rel_low$missed <- round(rel_low$cases_cxr / rel_low$sens, 0) - rel_low$cases_cxr

rel_high <- pop
rel_high$sens <- pmin(intercept + slope * (0.40 * (rel_high$cases_cxr / 1000)),1)
rel_high$missed <- round(rel_high$cases_cxr / rel_high$sens, 0) - rel_high$cases_cxr
fix <- pop
fix$sens <- 0.76  
fix$missed <- round(fix$cases_cxr / fix$sens, 0) - fix$cases_cxr

pop_long <- rbind(
  cbind(fix, scenario = "Fixed"),
  cbind(rel_low, scenario = "Relative: 20% severe"),
  cbind(rel_high, scenario = "Relative: 40% severe")
)
pop_long$prevalence <- pop_long$cases_cxr / 10  # Convert cases per 1000 to prevalence
tail(pop_long)

# Plot
ilo2_impact <- ggplot(pop_long, aes(x = prevalence, y = missed, color = scenario)) +
  geom_smooth(se = FALSE, method = "loess", span = 0.8, size = 0.9) +
  labs(
    x = "ILO ≥1/0 Silicosis Prevalence (%)",  # Adjusted label
    y = "Missed Cases per 1000 persons",
    colour = "Scenario"
  ) +
  theme_pubr(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank(), 
  ) +
  scale_color_brewer(palette = "Set2")
ilo2_impact

# Assuming cases_cxr represents the number of cases per 1000 (which is a prevalence rate)
pop2 <- data.frame(cases_cxr = c(50, 150, 300))  # cases per 1000, adjusted from your example

# Meta-regression coefficients and SEs from lm_ilo21 (ref-unconfirmed)
intercept <- 0.1595
slope <- 1.1468
se_intercept <- 0.1991
se_slope <- 0.3622

pop2$fixed <- 0.76
pop2$rel_sens_low <- pmin(intercept + (slope*0.20), 1)
pop2$rel_missed_low <- round(pop2$cases_cxr/pop2$rel_sens_low,0) - pop2$cases_cxr
pop2$rel_nns_low <-  round(1000/pop2$rel_missed_low,0)
pop2$rel_sens_high <- pmin(intercept + (slope*0.40), 1)
pop2$rel_missed_high <- round(pop2$cases_cxr/pop2$rel_sens_high,0) - pop2$cases_cxr
pop2$rel_nns_high <-  round(1000/pop2$rel_missed_high,0)
pop2$fixed_missed <- round(pop2$cases_cxr/pop2$fixed, 0) - pop2$cases_cxr
pop2$fixed_nns <-  round(1000/pop2$fixed_missed,0)
pop2


# Assuming cases_cxr represents the number of cases per 1000 (which is a prevalence rate)
pop2 <- data.frame(cases_cxr = seq(0, 400, by = 1))  # cases per 1000, adjusted from your example

fixed <- 0.76
rel_low2 <- pop2
rel_low2$sens <- pmin(intercept + (slope * 0.20),1)
rel_low2$missed <- round(rel_low2$cases_cxr / rel_low2$sens, 0) - rel_low2$cases_cxr

rel_high2 <- pop2
rel_high2$sens <- pmin(intercept + (slope *0.40),1)
rel_high2$missed <- round(rel_high2$cases_cxr / rel_high2$sens, 0) - rel_high2$cases_cxr

fix2 <- pop2
fix2$sens <- 0.76  
fix2$missed <- round(fix2$cases_cxr / fix2$sens, 0) - fix2$cases_cxr

pop_long2 <- rbind(
  cbind(fix2, scenario = "Fixed"),
  cbind(rel_low2, scenario = "Relative - Less severe (20%)"),
  cbind(rel_high2, scenario = "Relative - More severe (40%)")
)
pop_long2$prevalence <- pop_long2$cases_cxr / 10  # Convert cases per 1000 to prevalence
tail(pop_long2)

# Plot
ilo21_impact <- ggplot(pop_long2, aes(x = prevalence, y = missed, color = scenario)) +
  geom_smooth(se = FALSE, method = "loess", span = 0.8, size = 0.9) +
  labs(
    x = "ILO ≥1/0 Silicosis Prevalence (%)",  # Adjusted label
    y = "Missed Cases per 1000 persons",
    colour = "Scenario"
  ) +
  theme_pubr(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank(), 
  ) +
  scale_color_brewer(palette = "Set2")

ilo21_impact
 
ggsave("Supp_ilo2_impact.png", ilo2_impact, width = 6, height = 6)
ggsave("Supp_ilo21_impact.png", ilo21_impact, width = 6, height = 6)

# Arrange  and e with ggpubr, no titles 
impact <- ggpubr::ggarrange(ilo21_impact,ilo2_impact, # list of plots
                           labels = "AUTO", # labels
                           align = "hv", # Align them both, horizontal and vertical
                           ncol = 2, 
                           common.legend = TRUE,               # Create a shared legend
                           legend = "bottom" )  # number of rows
impact
ggsave("Supp_impact.jpg", impact, width = 10, height = 6)
