library(rema)
library(patchwork)
library(tidyverse)
library(here)

#---Fit smoother to capensis ----
read_csv(here("mods","hcr1","capensis.csv"))
?prepare_rema_input
cap<-prepare_rema_input(
  model_name = "Capensis",
  multi_survey = 0,
  admb_re = NULL,
  biomass_dat = read_csv(here("mods","hcr1","capensis.csv")) |> mutate(strata="Capensis"),
  cpue_dat = NULL,
  sum_cpue_index = FALSE,
  start_year = NULL,
  end_year = 2024,
  wt_biomass = NULL,
  wt_cpue = NULL,
  PE_options = NULL,
  q_options = NULL,
  zeros = NULL,
  extra_biomass_cv = NULL,
  extra_cpue_cv = NULL
)

m <- fit_rema(cap)

#---Cape-hakes, fit smoother to both species----
caphakes<-prepare_rema_input(
  model_name = "Same PE",
  multi_survey = 0,
  admb_re = NULL,
  biomass_dat = read_csv(here("mods","hcr1","CapeHakes.csv")),
  cpue_dat = NULL,
  sum_cpue_index = FALSE,
  start_year = NULL,
  end_year = 2024,
  wt_biomass = NULL,
  wt_cpue = NULL,
  PE_options = NULL,
  q_options = NULL,
  zeros = NULL,
  extra_biomass_cv = NULL,
  extra_cpue_cv = NULL
)

df_obs<- tibble(year=caphakes$data$model_yrs,
                capensis=caphakes$data$biomass_obs[,1],
                paradoxus=caphakes$data$biomass_obs[,2],
                p=paradoxus/(capensis+paradoxus))
df_obs

m <- fit_rema(caphakes)

caphakes2 <-prepare_rema_input(
  model_name = "Variable PE",
  biomass_dat = read_csv(here("mods","hcr1","CapeHakes.csv")),
  end_year = 2025,
)
caphakes2$data$pointer_PE_biomass
caphakes2$data$pointer_PE_biomass<-rep(0,2)
caphakes2$map$log_PE
caphakes2$map$log_PE <-  as.factor(caphakes$par$log_PE <- c(1,1))

m2 <- fit_rema(caphakes2)
names(m)
output <- tidy_rema(rema_model = m)
output <- tidy_rema(rema_model = m2)
names(output)
output$parameter_estimates # estimated fixed effects parameters
output$biomass_by_strata # data.frame of predicted and observed biomass by stratum
output$total_predicted_biomass # total predicted biomass (same as biomass_by_strata for univariate models)
unique(output$biomass_by_strata$strata) # total predicted biomass (same as biomass_by_strata for univariate models)
names(output$biomass_by_strata )
output$biomass_by_strata |>
  ggplot(aes(x=year,ymin=pred_lci,ymax=pred_uci,fill=strata,y=obs)) +
  geom_ribbon(alpha=.6) +
  scale_y_continuous(labels = scales::comma, expand = c(0, 0), limits = c(0, NA)) +
  geom_errorbar(aes(x = year, ymin = obs_lci, ymax = obs_uci),width=.2,alpha=.4) +
  geom_line(aes(x = year, y = pred)) + ylab("Biomass (kt)") +
  ggthemes::theme_few() + geom_point(alpha=.7)

# Sum over the species
tail(output$proportion_biomass_by_strata)
summary(output$proportion_biomass_by_strata$paradoxus)
min(output$proportion_biomass_by_strata$paradoxus)
summary((output$proportion_biomass_by_strata$paradoxus)/
mean(output$proportion_biomass_by_strata$paradoxus))
p1<-output$proportion_biomass_by_strata |>
  tidyr::pivot_longer(cols = -c(model_name, year)) |>
  group_by(name) |> mutate(mn=mean(value)) |> ungroup() |>
  ggplot(aes(x = factor(year), y = value,group=1, fill = reorder(name, (value)))) +
  #fill = reorder(area, desc(p)))) +
  geom_bar(position="stack", stat="identity") +
  #geom_point(stat='summary', fun.y=mn) +
  #stat_summary(fun.y=mn, geom="line")
  scale_fill_brewer(palette = 'Greys') +
  labs(x = NULL, y = NULL, fill = NULL,
       title = 'Proportion biomass by species') +
  scale_y_continuous(expand = expansion(mult = c(0, 0)));p1
# Compute proportions in data
p1$data |>
  filter(name=="capensis") |> transmute(year,p=(1-value)/(1-mn)) |> print(n=Inf)
p1$data |>
  filter(name=="capensis") |>
  ggplot(aes(x=year,y=1-value),color="blue") + geom_line(size=1,color="blue")+
  geom_line(aes(x=year,y=1-mn)) +
  ylab("Proportion paradoxus")+ ylim(c(0,1)) +
  ggthemes::theme_few() +
  geom_point(data=df_obs,aes(x=year,y=p),size=2,shape=3)

plots$proportion_biomass_by_strata$data

# Make a shade plot of control rul

df_hcr <- tibble(B_Btarget=seq(0,2,by=.01                    ))
df_hcr
# (5) Generate model plots
?plot_rema
plots <- plot_rema(tidy_rema = output,
                   # optional y-axis label
                   biomass_ylab = 'Biomass (kt)')
plots$biomass_by_strata
plots$total_predicted_biomass + ylim(c(0,1900) ) + ggthemes::theme_few()

?get_osa_residuals
osa <- get_osa_residuals(m, options = list(method = "cdf"))
cowplot::plot_grid(osa$plots$biomass_resids,
                   osa$plots$biomass_qqplot,
                   osa$plots$biomass_hist,
                   osa$plots$biomass_fitted,
                   ncol = 1)

# (7) Compare with ADMB RE model results
compare <- compare_rema_models(rema_models = list(m,m2),
                               biomass_ylab = 'Biomass (t)')
compare$plots$biomass_by_strata
compare$aic
names(compare)
names(compare$output)
