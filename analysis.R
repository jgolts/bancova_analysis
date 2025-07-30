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
  drop_na()

vickers_score_long <- vickers_score %>% 
  pivot_longer(cols = c(pk1, pk5), names_to = "time", values_to = "score") %>% 
  mutate(
    time = case_when(
      time == "pk1" ~ 0,
      .default = 1
    ))

vickers_meds <- vickers %>% 
  mutate(painmedspk1 = painmedspk1/4, painmedspk5 = painmedspk5/4) %>%
  select(id, painmedspk1, painmedspk5, group, age, sex, migraine,
         chronicity) %>% 
  drop_na()

# Frequentist analyses

anova_change <- lm(pk_change ~ group, data = vickers_score)

fancova1 <- lm(pk5 ~ pk1 + group, data = vickers_score)
fancova2 <- lm(painmedspk5 ~ painmedspk1 + group + age + sex + migraine
               + chronicity, data = vickers_meds)

fancova_inter <- lm(pk5 ~ pk1 + group + pk1:group, data = vickers_score)

lmm <- lmer(score ~ time * group + (1 | id), data = vickers_score_long)

## Frequentist results tables

sum_anova <- summary(anova_change)

sum_fancova1 <- summary(fancova1)
sum_fancova2 <- summary(fancova2)
sum_fancova_inter <- summary(fancova_inter)

sum_lmm <- summary(lmm, ddf = "Kenward-Roger")

ci_anova <- confint(anova_change)

ci_fancova1 <- confint(fancova1)
ci_fancova2 <- confint(fancova2)
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
  mutate(Coefficient = c("Intercept", "pk1", "Group"),
         `95% CI` = paste0("(", round(`2.5 %`, 1), ", ",
                           round(`97.5 %`, 1), ")"),
         Estimate = round(Estimate, 1),
         `P-Value` = round(`P-Value`, 4)) %>% 
  select(Coefficient, Estimate, `P-Value`, `95% CI`)

table_ancova2 <- bind_cols(sum_fancova2$coefficients, ci_fancova2) %>% 
  tibble() %>% 
  rename("P-Value" = "Pr(>|t|)") %>%
  mutate(Coefficient = c("Intercept", "pk1", "Group", "Age", "Sex", "Migraine",
                         "Chronicity"),
         `95% CI` = paste0("(", round(`2.5 %`, 1), ", ",
                           round(`97.5 %`, 1), ")"),
         Estimate = round(Estimate, 1),
         `P-Value` = round(`P-Value`, 4)) %>% 
  select(Coefficient, Estimate, `P-Value`, `95% CI`)

table_ancova_inter <- bind_cols(sum_fancova_inter$coefficients,
                                ci_fancova_inter) %>% 
  tibble() %>% 
  rename("P-Value" = "Pr(>|t|)") %>%
  mutate(Coefficient = c("Intercept", "pk1", "Group", "pk1:Group"),
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

## Bayesian analyses

## sd(y) for models

sd(vickers_score$pk5)
sd(vickers_meds$painmedspk5)

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
  labs(x = "Coefficient",
       y = "Probability density")

priorplot_scep <- ggplot(df_scep, aes(x = x, y = y)) +
  geom_line(linewidth = 1.2) + 
  labs(x = "Treatment effect",
       y = "Probability density")

priorplot_enth <- ggplot(df_enth, aes(x = x, y = y)) +
  geom_line(linewidth = 1.2) + 
  labs(x = "Treatment effect",
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

b_ancova_neu <- brm(pk5 ~ pk1 + group, data = vickers_score, silent = 1,
                    thin = 1, cores = 8, iter = 20000, warmup = 10000,
                    prior = c(prior(normal(0, 100), class = "Intercept"),
                              prior(normal(0, 100), class = "b", coef = "pk1"),
                              prior(normal(0, 100), class = "b", coef = "group"),
                              prior(cauchy(0, 8), class = "sigma")),
                    seed = 157, chains = 4)

b_ancova_scep <- brm(pk5 ~ pk1 + group, data = vickers_score, silent = 1,
                     thin = 1, cores = 8, iter = 20000, warmup = 10000,
                     prior = c(prior(normal(0, 100), class = "Intercept"),
                               prior(normal(0, 100), class = "b", coef = "pk1"),
                               prior(normal(0, 0.5), class = "b",
                                     coef = "group"),
                               prior(cauchy(0, 8), class = "sigma")),
                     seed = 157, chains = 4)

b_ancova_enth <- brm(pk5 ~ pk1 + group, data = vickers_score, silent = 1,
                     thin = 1, cores = 8, iter = 20000, warmup = 10000,
                     prior = c(prior(normal(0, 100), class = "Intercept"),
                               prior(normal(0, 100), class = "b", coef = "pk1"),
                               prior(normal(-10, 0.5), class = "b",
                                     coef = "group"),
                               prior(cauchy(0, 8), class = "sigma")),
                     seed = 157, chains = 4)

b_ancova_meds <- brm(painmedspk5 ~ painmedspk1 + group + age + sex + migraine
                     + chronicity, data = vickers_meds,
                     silent = 1, thin = 1, cores = 8, iter = 20000,
                     warmup = 10000,
                     prior = c(prior(normal(0, 100), class = "Intercept"),
                               prior(normal(0, 100), class = "b",
                                     coef = "painmedspk1"),
                               prior(normal(0, 100), class = "b", coef = "group"),
                               prior(cauchy(0, 8), class = "sigma")),
                     seed = 157, chains = 4)

## Doubled iters

b_ancova_neu_x2 <- brm(pk5 ~ pk1 + group, data = vickers_score, silent = 1,
                       thin = 1, cores = 8, iter = 40000, warmup = 20000,
                       prior = c(prior(normal(0, 100), class = "Intercept"),
                                 prior(normal(0, 100), class = "b", coef = "pk1"),
                                 prior(normal(0, 100), class = "b",
                                       coef = "group"),
                                 prior(cauchy(0, 8), class = "sigma")),
                       seed = 157, chains = 4)

b_ancova_scep_x2 <- brm(pk5 ~ pk1 + group, data = vickers_score, silent = 1,
                        thin = 1, cores = 8, iter = 40000, warmup = 20000,
                        prior = c(prior(normal(0, 100), class = "Intercept"),
                                  prior(normal(0, 100), class = "b", coef = "pk1"),
                                  prior(normal(0, 0.5), class = "b",
                                        coef = "group"),
                                  prior(cauchy(0, 8), class = "sigma")),
                        seed = 157, chains = 4)

b_ancova_enth_x2 <- brm(pk5 ~ pk1 + group, data = vickers_score, silent = 1,
                        thin = 1, cores = 8, iter = 40000, warmup = 20000,
                        prior = c(prior(normal(0, 100), class = "Intercept"),
                                  prior(normal(0, 100), class = "b", coef = "pk1"),
                                  prior(normal(-10, 0.5), class = "b",
                                        coef = "group"),
                                  prior(cauchy(0, 8), class = "sigma")),
                        seed = 157, chains = 4)

b_ancova_meds_x2 <- brm(painmedspk5 ~ painmedspk1 + group + age + sex + migraine
                        + chronicity, data = vickers_meds,
                        silent = 1, thin = 1, cores = 8, iter = 40000,
                        warmup = 20000,
                        prior = c(prior(normal(0, 100), class = "Intercept"),
                                  prior(normal(0, 100), class = "b",
                                        coef = "painmedspk1"),
                                  prior(normal(0, 100), class = "b", coef = "group"),
                                  prior(cauchy(0, 8), class = "sigma")),
                        seed = 157, chains = 4)

## Diagnostics

## Traces

mcmc_plot(b_ancova_neu, type = "trace")
mcmc_plot(b_ancova_scep, type = "trace")
mcmc_plot(b_ancova_enth, type = "trace")
mcmc_plot(b_ancova_meds, type = "trace")

mcmc_plot(b_ancova_neu_x2, type = "trace")
mcmc_plot(b_ancova_scep_x2, type = "trace")
mcmc_plot(b_ancova_enth_x2, type = "trace")
mcmc_plot(b_ancova_meds_x2, type = "trace")

## Histograms

mcmc_plot(b_ancova_neu, type = "hist")
mcmc_plot(b_ancova_scep, type = "hist")
mcmc_plot(b_ancova_enth, type = "hist")
mcmc_plot(b_ancova_meds, type = "hist")

## Autocorrelation

mcmc_plot(b_ancova_neu, type = "acf")
mcmc_plot(b_ancova_scep, type = "acf")
mcmc_plot(b_ancova_enth, type = "acf")
mcmc_plot(b_ancova_meds, type = "acf")

## Bayesian results tables

draws_bancova_neu <- as_draws_df(b_ancova_neu)
draws_bancova_scep <- as_draws_df(b_ancova_scep)
draws_bancova_enth <- as_draws_df(b_ancova_enth)
draws_bancova_meds <- as_draws_df(b_ancova_meds)

sum_bancova_neu <- summary(b_ancova_neu)
sum_bancova_scep <- summary(b_ancova_scep)
sum_bancova_enth <- summary(b_ancova_enth)
sum_bancova_meds <- summary(b_ancova_meds)

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
sd_bancova_meds <- draws_bancova_meds %>% 
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
sigma_ci_meds <- bayestestR::hdi(draws_bancova_meds$sigma) %>% 
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
ci_bancova_meds <- bayestestR::hdi(b_ancova_meds) %>% 
  select(Parameter, CI, CI_low, CI_high) %>% 
  rbind(sigma_ci_meds)

table_neu <- rbind(sum_bancova_neu$fixed, sum_bancova_neu$spec_pars) %>% 
  cbind(ci_bancova_neu) %>%
  cbind(sd_bancova_neu) %>% 
  tibble() %>% 
  mutate(Coefficient = c("Intercept", "pk1", "Group", "sigma"),
         Estimate = round(Estimate, 1),
         Rhat = round(Rhat, 2),
         `95% HDPI` = paste0("(", round(CI_low, 1), ", ",
                             round(CI_high, 1), ")"),
         SD = round(SD, 1))%>% 
  select(Coefficient, Estimate, `95% HDPI`, SD, Rhat)

table_scep <- rbind(sum_bancova_scep$fixed, sum_bancova_scep$spec_pars) %>% 
  cbind(ci_bancova_scep) %>%
  cbind(sd_bancova_scep) %>% 
  tibble() %>% 
  mutate(Coefficient = c("Intercept", "pk1", "Group", "sigma"),
         Estimate = round(Estimate, 1),
         Rhat = round(Rhat, 2),
         `95% HDPI` = paste0("(", round(CI_low, 1), ", ",
                             round(CI_high, 1), ")"),
         SD = round(SD, 1))%>% 
  select(Coefficient, Estimate, `95% HDPI`, SD, Rhat)

table_enth <- rbind(sum_bancova_enth$fixed, sum_bancova_enth$spec_pars) %>% 
  cbind(ci_bancova_enth) %>%
  cbind(sd_bancova_enth) %>% 
  tibble() %>% 
  mutate(Coefficient = c("Intercept", "pk1", "Group", "sigma"),
         Estimate = round(Estimate, 1),
         Rhat = round(Rhat, 2),
         `95% HDPI` = paste0("(", round(CI_low, 1), ", ",
                             round(CI_high, 1), ")"),
         SD = round(SD, 1))%>% 
  select(Coefficient, Estimate, `95% HDPI`, SD, Rhat)

table_meds <- rbind(sum_bancova_meds$fixed, sum_bancova_meds$spec_pars) %>% 
  cbind(ci_bancova_meds) %>% 
  cbind(sd_bancova_meds) %>% 
  tibble() %>% 
  mutate(Coefficient = c("Intercept", "pk1", "Group", "Age", "Sex", "Migraine",
                         "Chronicity", "sigma"),
         Estimate = round(Estimate, 1),
         Rhat = round(Rhat, 2),
         `95% HDPI` = paste0("(", round(CI_low, 1), ", ",
                             round(CI_high, 1), ")"),
         SD = round(SD, 1))%>% 
  select(Coefficient, Estimate, `95% HDPI`, SD, Rhat)

## Summaries just for the Rhat value (double iters)

summary(b_ancova_neu_x2)
summary(b_ancova_scep_x2)
summary(b_ancova_enth_x2)
summary(b_ancova_meds_x2)

## Posterior probabilities

mean(draws_bancova_neu$b_group < 0)

mean(draws_bancova_scep$b_group < 0)

mean(draws_bancova_enth$b_group < 0)

mean(draws_bancova_meds$b_group < 0)

## Preparation for posterior visualisation

draws_neu_long <- draws_bancova_neu %>% 
  select(b_Intercept, b_pk1, b_group, sigma) %>% 
  pivot_longer(cols = everything(),
               names_to = "Parameter", values_to = "Value")

draws_scep_long <- draws_bancova_scep %>% 
  select(b_Intercept, b_pk1, b_group, sigma) %>% 
  pivot_longer(cols = everything(),
               names_to = "Parameter", values_to = "Value")

draws_enth_long <- draws_bancova_enth %>% 
  select(b_Intercept, b_pk1, b_group, sigma) %>% 
  pivot_longer(cols = everything(),
               names_to = "Parameter", values_to = "Value")

draws_meds_long <- draws_bancova_meds %>% 
  select(b_Intercept, b_painmedspk1, b_group, sigma) %>% 
  pivot_longer(cols = everything(),
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
         x = paste0(param, " Estimate"),
         y = NULL,
         fill = "HDI")
}

plotPosteriorCDF <- function(df, param) {
  
  df <- df %>%
    filter(Parameter == param)
  
  ggplot(data = df, aes(x = Parameter, y = Value)) +
    stat_ccdfinterval(aes(fill = after_stat(level)), .width = c(.95, 1),
                      point_interval = mean_hdi, justification = 1) +
    scale_fill_brewer(na.translate = F, palette = "YlOrRd") +
    labs(title = "Posterior CDF with HDI",
         x = paste0(param, " Estimate"),
         y = NULL,
         fill = "HDI")
}

plotPosterior(draws_neu_long, "b_group")

ggsave(
  "teplot_neu.png",
  height = 4,
  width = 6.5,
  dpi = 600
)

plotPosterior(draws_scep_long, "b_group")

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

plotPosterior(draws_meds_long, "b_group")

ggsave(
  "teplot_meds.png",
  height = 4,
  width = 6.5,
  dpi = 600
)

plotPosteriorCDF(draws_neu_long, "b_group")
plotPosteriorCDF(draws_scep_long, "b_group")
plotPosteriorCDF(draws_enth_long, "b_group")
plotPosteriorCDF(draws_meds_long, "b_group")