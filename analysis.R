library(tidyverse)
library(readxl)
library(parallel)
library(lmerTest)
library(parameters)
library(brms)
library(bayestestR)
library(ggplot2)
library(ggdist)

set.seed(157)

##Data preparation

vickers <- read_xls("Vickers2006.xls", sheet = 3)

vickers_score <- vickers %>%
  mutate(pk_change = pk5 - pk1) %>% 
  select(id, pk1, pk5, group, pk_change) %>% 
  drop_na() %>% 
  rename("pre" = "pk1")

vickers_score_long <- vickers_score %>% 
  pivot_longer(cols = c(pre, pk5), names_to = "time", values_to = "score") %>% 
  mutate(
    time = case_when(
      time == "pre" ~ 0,
      .default = 1
    ))

## Frequentist analyses

## Group means

vickers_means <- vickers_score %>% 
  group_by(group) %>% 
  summarise(mean_pre = mean(pre), mean_post = mean(pk5)) %>% 
  mutate(mean_pre = round(mean_pre, 1), mean_post = round(mean_post, 1))

## Frequentist models

anova_change <- lm(pk_change ~ group, data = vickers_score)

fancova1 <- lm(pk5 ~ pre + group, data = vickers_score)

fancova_inter <- lm(pk5 ~ pre + group + pre:group, data = vickers_score)

lmm <- lmer(score ~ time * group + (1 | id), data = vickers_score_long)

## Frequentist results tables

sum_anova <- summary(anova_change)

sum_fancova1 <- summary(fancova1)
sum_fancova_inter <- summary(fancova_inter)

sum_lmm <- summary(lmm, ddf = "Kenward-Roger")

ci_anova <- confint(anova_change)

ci_fancova1 <- confint(fancova1)
ci_fancova_inter <- confint(fancova_inter)

ci_lmm <- model_parameters(lmm, ci_method = "kenward")

table_anova1 <- bind_cols(sum_anova$coefficients, ci_anova) %>% 
  tibble() %>% 
  rename("P-Value" = "Pr(>|t|)") %>%
  mutate(Coefficient = c("Intercept", "Group"),
         `95% CI` = paste0("(", round(`2.5 %`, 1), ", ",
                           round(`97.5 %`, 1), ")"),
         Estimate = round(Estimate, 1),
         `P-Value` = round(`P-Value`, 4)) %>%
  select(Coefficient, Estimate, `P-Value`, `95% CI`)

table_ancova1 <- bind_cols(sum_fancova1$coefficients, ci_fancova1) %>% 
  tibble() %>% 
  rename("P-Value" = "Pr(>|t|)") %>%
  mutate(Coefficient = c("Intercept", "pre", "Group"),
         `95% CI` = paste0("(", round(`2.5 %`, 1), ", ",
                           round(`97.5 %`, 1), ")"),
         Estimate = round(Estimate, 1),
         `P-Value` = round(`P-Value`, 4)) %>% 
  select(Coefficient, Estimate, `P-Value`, `95% CI`)

table_ancova_inter <- bind_cols(sum_fancova_inter$coefficients,
                                ci_fancova_inter) %>% 
  tibble() %>% 
  rename("P-Value" = "Pr(>|t|)") %>%
  mutate(Coefficient = c("Intercept", "pre", "Group", "pk1:Group"),
         `95% CI` = paste0("(", round(`2.5 %`, 1), ", ",
                           round(`97.5 %`, 1), ")"),
         Estimate = round(Estimate, 1),
         `P-Value` = round(`P-Value`, 4)) %>% 
  select(Coefficient, Estimate, `P-Value`, `95% CI`)

table_lmm <- bind_cols(sum_lmm$coefficients, drop_na(as_tibble(ci_lmm$CI_low)),
                       drop_na(as_tibble(ci_lmm$CI_high))) %>% 
  rename("P-Value" = "Pr(>|t|)") %>%
  mutate(Coefficient = c("Intercept", "Time", "Group", "Time:Group"),
         `95% CI` = paste0("(", round(`value...6`, 1), ", ",
                           round(`value...7`, 1), ")"),
         Estimate = round(Estimate, 1),
         `P-Value` = round(`P-Value`, 4)) %>% 
  select(Coefficient, Estimate, `P-Value`, `95% CI`)

##Plotting with interaction

plot(pk5 ~ pre, data = vickers_score,
     col = ifelse(vickers_score$group == "1", "red", "blue"),
     pch = 16, xlab = "Baseline (Pre)", ylab = "Outcome (Post)")

xseq <- seq(min(vickers_score$pre), max(vickers_score$pk5), length.out = 100)

lines(xseq, 0.4 + 0.8 * xseq, col = "blue", lwd = 2)

lines(xseq, (0.4 + 1.9) + (0.8 + -0.3) * xseq, col = "red", lwd = 2)

lines(xseq, 3.5 + 0.7 * xseq - 4.6, col = "orange", lwd = 2)

## Bayesian analyses

## sd(y) for models

sd(vickers_score$pk5)

## Prior plots

x_vals_neu <- seq(0 - 4 * 100, 0 + 4 * 100, length.out = 1000)
df_neu <- tibble(x = x_vals_neu, y = dnorm(x_vals_neu, mean = 0, sd = 100))

x_vals_scep <- seq(0 - 4 * 0.5, 0 + 4 * 0.5, length.out = 1000)
df_scep <- tibble(x = x_vals_scep, y = dnorm(x_vals_scep, mean = 0, sd = 0.5))

x_vals_enth <- seq(-10 - 4 * 0.5, -10 + 4 * 0.5, length.out = 1000)
df_enth <- tibble(x = x_vals_enth, y = dnorm(x_vals_enth, mean = -10, sd = 0.5))

x_vals_sigma <- seq(0, 8 * 10, length.out = 1000)
df_sigma <- tibble(x = x_vals_sigma, y = dcauchy(x_vals_sigma, 0, 8))

priorplot_neu <- ggplot(df_neu, aes(x = x, y = y)) +
  geom_line(linewidth = 1.2) + 
  labs(x = "Neutral",
       y = "Probability density")

priorplot_scep <- ggplot(df_scep, aes(x = x, y = y)) +
  geom_line(linewidth = 1.2) + 
  labs(x = "Sceptical",
       y = "Probability density")

priorplot_enth <- ggplot(df_enth, aes(x = x, y = y)) +
  geom_line(linewidth = 1.2) + 
  labs(x = "Enthusiastic",
       y = "Probability density")

priorplot_sigma <- ggplot(df_sigma, aes(x = x, y = y)) +
  geom_line(linewidth = 1.2) + 
  labs(x = "SD of model outcome",
       y = "Probability density")

priorplots <- gridExtra::arrangeGrob(priorplot_neu, priorplot_scep,
                                     priorplot_sigma, priorplot_enth)

ggsave(
  "priorsplot.png",
  priorplots,
  height = 4,
  width = 6.5,
  dpi = 600
)

## brms model fitting

b_ancova_neu <- brm(pk5 ~ pre + group, data = vickers_score, silent = 1,
                    thin = 1, cores = 8, iter = 20000, warmup = 10000,
                    prior = c(prior(normal(0, 100), class = "Intercept"),
                              prior(normal(0, 100), class = "b", coef = "pre"),
                              prior(normal(0, 100), class = "b",
                                    coef = "group"),
                              prior(cauchy(0, 8), class = "sigma")),
                    seed = 157, chains = 4)

b_ancova_scep <- brm(pk5 ~ pre + group, data = vickers_score, silent = 1,
                     thin = 1, cores = 8, iter = 20000, warmup = 10000,
                     prior = c(prior(normal(0, 100), class = "Intercept"),
                               prior(normal(0, 100), class = "b", coef = "pre"),
                               prior(normal(0, 0.5), class = "b",
                                     coef = "group"),
                               prior(cauchy(0, 8), class = "sigma")),
                     seed = 157, chains = 4)

b_ancova_enth <- brm(pk5 ~ pre + group, data = vickers_score, silent = 1,
                     thin = 1, cores = 8, iter = 20000, warmup = 10000,
                     prior = c(prior(normal(0, 100), class = "Intercept"),
                               prior(normal(0, 100), class = "b", coef = "pre"),
                               prior(normal(-10, 0.5), class = "b",
                                     coef = "group"),
                               prior(cauchy(0, 8), class = "sigma")),
                     seed = 157, chains = 4)

## Doubled iters

b_ancova_neu_x2 <- brm(pk5 ~ pre + group, data = vickers_score, silent = 1,
                       thin = 1, cores = 8, iter = 40000, warmup = 20000,
                       prior = c(prior(normal(0, 100), class = "Intercept"),
                                 prior(normal(0, 100), class = "b",
                                       coef = "pre"),
                                 prior(normal(0, 100), class = "b",
                                       coef = "group"),
                                 prior(cauchy(0, 8), class = "sigma")),
                       seed = 157, chains = 4)

b_ancova_scep_x2 <- brm(pk5 ~ pre + group, data = vickers_score, silent = 1,
                        thin = 1, cores = 8, iter = 40000, warmup = 20000,
                        prior = c(prior(normal(0, 100), class = "Intercept"),
                                  prior(normal(0, 100), class = "b",
                                        coef = "pre"),
                                  prior(normal(0, 0.5), class = "b",
                                        coef = "group"),
                                  prior(cauchy(0, 8), class = "sigma")),
                        seed = 157, chains = 4)

b_ancova_enth_x2 <- brm(pk5 ~ pre + group, data = vickers_score, silent = 1,
                        thin = 1, cores = 8, iter = 40000, warmup = 20000,
                        prior = c(prior(normal(0, 100), class = "Intercept"),
                                  prior(normal(0, 100), class = "b",
                                        coef = "pre"),
                                  prior(normal(-10, 0.5), class = "b",
                                        coef = "group"),
                                  prior(cauchy(0, 8), class = "sigma")),
                        seed = 157, chains = 4)

## Diagnostics

## Traces

mcmc_plot(b_ancova_neu, type = "trace")

ggsave(
  "traceplot_neu.png",
  height = 4,
  width = 6.5,
  dpi = 600
)

mcmc_plot(b_ancova_scep, type = "trace")

ggsave(
  "traceplot_scep.png",
  height = 4,
  width = 6.5,
  dpi = 600
)

mcmc_plot(b_ancova_enth, type = "trace")

ggsave(
  "traceplot_enth.png",
  height = 4,
  width = 6.5,
  dpi = 600
)

mcmc_plot(b_ancova_neu_x2, type = "trace")

ggsave(
  "traceplot_neu_x2.png",
  height = 4,
  width = 6.5,
  dpi = 600
)

mcmc_plot(b_ancova_scep_x2, type = "trace")

ggsave(
  "traceplot_scep_x2.png",
  height = 4,
  width = 6.5,
  dpi = 600
)

mcmc_plot(b_ancova_enth_x2, type = "trace")

ggsave(
  "traceplot_enth_x2.png",
  height = 4,
  width = 6.5,
  dpi = 600
)

## Histograms

mcmc_plot(b_ancova_neu, type = "hist")

ggsave(
  "histplot_neu.png",
  height = 4,
  width = 6.5,
  dpi = 600
)

mcmc_plot(b_ancova_scep, type = "hist")

ggsave(
  "histplot_scep.png",
  height = 4,
  width = 6.5,
  dpi = 600
)

mcmc_plot(b_ancova_enth, type = "hist")

ggsave(
  "histplot_enth.png",
  height = 4,
  width = 6.5,
  dpi = 600
)

## Autocorrelation

mcmc_plot(b_ancova_neu, type = "acf")

ggsave(
  "acfplot_neu.png",
  height = 4,
  width = 6.5,
  dpi = 600
)

mcmc_plot(b_ancova_scep, type = "acf")

ggsave(
  "acfplot_scep.png",
  height = 4,
  width = 6.5,
  dpi = 600
)

mcmc_plot(b_ancova_enth, type = "acf")

ggsave(
  "acfplot_enth.png",
  height = 4,
  width = 6.5,
  dpi = 600
)

## Bayesian results tables

draws_bancova_neu <- as_draws_df(b_ancova_neu) %>% 
  mutate(model = "Neutral")
draws_bancova_scep <- as_draws_df(b_ancova_scep) %>% 
  mutate(model = "Sceptical")
draws_bancova_enth <- as_draws_df(b_ancova_enth) %>% 
  mutate(model = "Enthusiastic")

sum_bancova_neu <- summary(b_ancova_neu)
sum_bancova_scep <- summary(b_ancova_scep)
sum_bancova_enth <- summary(b_ancova_enth)

sd_bancova_neu <- draws_bancova_neu %>% 
  select(b_Intercept:sigma) %>% 
  sapply(sd) %>% 
  as_tibble() %>% 
  rename("SD" = "value")
sd_bancova_scep <- draws_bancova_scep %>% 
  select(b_Intercept:sigma) %>% 
  sapply(sd) %>% 
  as_tibble() %>% 
  rename("SD" = "value")
sd_bancova_enth <- draws_bancova_enth %>% 
  select(b_Intercept:sigma) %>%
  sapply(sd) %>%
  as_tibble() %>% 
  rename("SD" = "value")

sigma_ci_neu <- bayestestR::hdi(draws_bancova_neu$sigma) %>% 
  mutate(Parameter = "sigma") %>% 
  select(Parameter, CI, CI_low, CI_high)
sigma_ci_scep <- bayestestR::hdi(draws_bancova_scep$sigma) %>% 
  mutate(Parameter = "sigma") %>% 
  select(Parameter, CI, CI_low, CI_high)
sigma_ci_enth <- bayestestR::hdi(draws_bancova_enth$sigma) %>% 
  mutate(Parameter = "sigma") %>% 
  select(Parameter, CI, CI_low, CI_high)

ci_bancova_neu <- bayestestR::hdi(b_ancova_neu) %>% 
  select(Parameter, CI, CI_low, CI_high) %>% 
  rbind(sigma_ci_neu)
ci_bancova_scep <- bayestestR::hdi(b_ancova_scep) %>% 
  select(Parameter, CI, CI_low, CI_high) %>% 
  rbind(sigma_ci_scep)
ci_bancova_enth <- bayestestR::hdi(b_ancova_enth) %>% 
  select(Parameter, CI, CI_low, CI_high) %>% 
  rbind(sigma_ci_enth)

table_neu <- rbind(sum_bancova_neu$fixed, sum_bancova_neu$spec_pars) %>% 
  cbind(ci_bancova_neu) %>%
  cbind(sd_bancova_neu) %>% 
  tibble() %>% 
  mutate(Coefficient = c("Intercept", "pk1", "Group", "sigma"),
         Estimate = round(Estimate, 1),
         Rhat = round(Rhat, 2),
         `95% HDI` = paste0("(", round(CI_low, 1), ", ",
                             round(CI_high, 1), ")"),
         SD = round(SD, 1))%>% 
  select(Coefficient, Estimate, `95% HDI`, SD, Rhat)

table_scep <- rbind(sum_bancova_scep$fixed, sum_bancova_scep$spec_pars) %>% 
  cbind(ci_bancova_scep) %>%
  cbind(sd_bancova_scep) %>% 
  tibble() %>% 
  mutate(Coefficient = c("Intercept", "pk1", "Group", "sigma"),
         Estimate = round(Estimate, 1),
         Rhat = round(Rhat, 2),
         `95% HDI` = paste0("(", round(CI_low, 1), ", ",
                             round(CI_high, 1), ")"),
         SD = round(SD, 1))%>% 
  select(Coefficient, Estimate, `95% HDI`, SD, Rhat)

table_enth <- rbind(sum_bancova_enth$fixed, sum_bancova_enth$spec_pars) %>% 
  cbind(ci_bancova_enth) %>%
  cbind(sd_bancova_enth) %>% 
  tibble() %>% 
  mutate(Coefficient = c("Intercept", "pk1", "Group", "sigma"),
         Estimate = round(Estimate, 1),
         Rhat = round(Rhat, 2),
         `95% HDI` = paste0("(", round(CI_low, 1), ", ",
                             round(CI_high, 1), ")"),
         SD = round(SD, 1))%>% 
  select(Coefficient, Estimate, `95% HDI`, SD, Rhat)

## Summaries just for the Rhat value (double iters)

summary(b_ancova_neu_x2)
summary(b_ancova_scep_x2)
summary(b_ancova_enth_x2)

## Posterior probabilities

mean(draws_bancova_neu$b_group < 0)
mean(draws_bancova_neu$b_group < -3)

mean(draws_bancova_scep$b_group < 0)
mean(draws_bancova_scep$b_group < -3)

mean(draws_bancova_enth$b_group < 0)
mean(draws_bancova_enth$b_group < -3)

## Preparation for posterior visualisation

draws_all <- bind_rows(draws_bancova_neu, draws_bancova_scep,
                       draws_bancova_enth)

draws_neu_long <- draws_bancova_neu %>% 
  select(b_Intercept, b_pre, b_group, sigma) %>% 
  pivot_longer(cols = everything(),
               names_to = "Parameter", values_to = "Value")

draws_scep_long <- draws_bancova_scep %>% 
  select(b_Intercept, b_pre, b_group, sigma) %>% 
  pivot_longer(cols = everything(),
               names_to = "Parameter", values_to = "Value")

draws_enth_long <- draws_bancova_enth %>% 
  select(b_Intercept, b_pre, b_group, sigma) %>% 
  pivot_longer(cols = everything(),
               names_to = "Parameter", values_to = "Value")

draws_all_long <- draws_all %>% 
  select(b_Intercept, b_pre, b_group, sigma, model) %>% 
  pivot_longer(cols = b_Intercept:sigma,
               names_to = "Parameter", values_to = "Value")

## Group posterior plots

plotPosterior <- function(df, param) {
  
  df <- df %>%
    filter(Parameter == param)
  
  ggplot(data = df, aes(x = Value, y = Parameter)) +
    stat_halfeye(aes(fill = after_stat(level)), .width = c(.95, 1),
                 point_interval = mean_hdi) +
    scale_fill_brewer(na.translate = F, palette = "YlOrRd") +
    labs(title = "Posterior PDF with HDI",
         x = paste0("Estimated difference in headache score at 1 year"),
         y = NULL,
         fill = "HDI")
}

plotPosteriorTEline <- function(df, param) {
  
  df <- df %>%
    filter(Parameter == param)
  
  ggplot(data = df, aes(x = Value, y = Parameter)) +
    stat_halfeye(aes(fill = after_stat(level)), .width = c(.95, 1),
                 point_interval = mean_hdi) +
    geom_vline(xintercept = 0) +
    scale_fill_brewer(na.translate = F, palette = "YlOrRd") +
    labs(title = "Posterior PDF with HDI",
         x = paste0("Estimated difference in headache score at 1 year"),
         y = NULL,
         fill = "HDI")
}

plotPosteriorCDF <- function(df, param) {
  
  df <- df %>%
    filter(Parameter == param)
  
  ggplot(data = df, aes(x = Parameter, y = Value)) +
    stat_ccdfinterval(aes(fill = after_stat(level)), .width = c(.95, 1),
                      point_interval = mean_hdi, justification = 1) +
    geom_vline(xintercept = 0) +
    scale_fill_brewer(na.translate = F, palette = "YlOrRd") +
    labs(title = "Posterior CDF with HDI",
         x = paste0(param, " Estimate"),
         y = NULL,
         fill = "HDI")
}

plotPosteriorsTogether <- function(df, param) {
  
  df <- df %>%
    filter(Parameter == param)
  
  ggplot(data = df, aes(x = Value, y = model)) +
    stat_halfeye(aes(fill = after_stat(level)), .width = c(.95, 1),
                 point_interval = mean_hdi) +
    geom_vline(xintercept = 0) +
    scale_fill_brewer(na.translate = F, palette = "YlOrRd") +
    labs(title = "Posterior PDF with HDI",
         x = paste0("Estimated difference in headache score at 1 year"),
         y = NULL,
         fill = "HDI")
}

plotPosteriorTEline(draws_neu_long, "b_group")

ggsave(
  "teplot_neu.png",
  height = 4,
  width = 6.5,
  dpi = 600
)

plotPosteriorTEline(draws_scep_long, "b_group")

ggsave(
  "teplot_scep.png",
  height = 4,
  width = 6.5,
  dpi = 600
)

plotPosterior(draws_enth_long, "b_group")

ggsave(
  "teplot_enth.png",
  height = 4,
  width = 6.5,
  dpi = 600
)

plotPosteriorsTogether(draws_all_long, "b_group")

ggsave(
  "teplot_all.png",
  height = 4,
  width = 6.5,
  dpi = 600
)

plotPosteriorCDF(draws_neu_long, "b_group")
plotPosteriorCDF(draws_scep_long, "b_group")
plotPosteriorCDF(draws_enth_long, "b_group")