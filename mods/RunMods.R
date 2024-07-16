## Quick file to explore some things we can do with the TMB
## version of the model.
library(here)
library(tidyverse)
theme_set(ggthemes::theme_few())
# Read in nh_r.rep file
source("R/read-admb2.R")
ctl_dat <- read_dat(here("mods", "test", "log_input.rep"))
View(ctl_dat)

# Create a new package -------------------------------------------------
library(usethis)
path <- file.path(here())
create_package(path)
proj_activate(path)
use_description()
use_mit_license("Jim Ianelli")

use_package("ggplot2", "Suggests")
#> ✔ Adding 'ggplot2' to Suggests field in DESCRIPTION
#> • Use `requireNamespace("ggplot2", quietly = TRUE)` to test if package is installed
#> • Then directly refer to functions with `ggplot2::fun()`

# Set up other files -------------------------------------------------
use_readme_md()
#> ✔ Writing 'README.md'
#> • Update 'README.md' to include installation instructions.

use_news_md()
#> ✔ Writing 'NEWS.md'

use_test("my-test")
#> ✔ Adding 'testthat' to Suggests field in DESCRIPTION
#> ✔ Adding '3' to Config/testthat/edition
#> ✔ Creating 'tests/testthat/'
#> ✔ Writing 'tests/testthat.R'
#> ✔ Writing 'tests/testthat/test-my-test.R'
#> • Edit 'tests/testthat/test-my-test.R'

x <- 1
y <- 2
use_data(x, y)
#> ✔ Adding 'R' to Depends field in DESCRIPTION
#> ✔ Creating 'data/'
#> ✔ Setting LazyData to 'true' in 'DESCRIPTION'
#> ✔ Saving 'x', 'y' to 'data/x.rda', 'data/y.rda'
#> • Document your data (see 'https://r-pkgs.org/data.html')

# Use git ------------------------------------------------------------
use_git()
#> ✔ Initialising Git repo
#> ✔ Adding '.Rproj.user', '.Rhis

#--Read in age comparison datas-----
df_age <- read_csv(here("mods", "data", "AgeCompare.csv"))
glimpse(df_age)
df_age |> pivot_longer(cols=2:11, names_to="Age", values_to="Value") |>
  mutate(Age=as.numeric(Age)) |>
  ggplot(aes(x = Age, y = Value, color = Source, shape=Source )) +
  ggtitle("Survey 1") +
  geom_point() +
  geom_line(stat='Identity') +
  facet_wrap(. ~ Year)
ggsave(here("mods","figs","AgeCompareSurv.png"),width=9,height=8)
df_age <- read_csv(here("mods", "data", "AgeCompareFish.csv"))
glimpse(df_age)
df_age |> pivot_longer(cols=2:11, names_to="Age", values_to="Value") |>
  mutate(Age=as.numeric(Age)) |>
  ggplot(aes(x = Age, y = Value, color = Source, shape=Source )) +
  ggtitle("Fishery ") +
  geom_point() +
  geom_line(stat='Identity') +
  facet_wrap(. ~ Year)
ggsave(here("mods","figs","AgeCompareFish.png"),width=9,height=8)

#--Read in model output files-----

  obc <- read_rep(here("mods", "Mod1", "nh_R.rep"))
  bc <- read_rep(here("mods", "bc", "nh_R.rep"))
  m1 <- read_rep(here("mods", "m1", "nh_R.rep"))

  m2 <- read_rep(here("mods", "m2", "nh_R.rep"))
  m3 <- read_rep(here("mods", "m3", "nh_R.rep"))
  m4 <- read_rep(here("mods", "m4", "nh_R.rep"))
  m5 <- read_rep(here("mods", "m5", "nh_R.rep"))
  m6 <- read_rep(here("mods", "m6", "nh_R.rep"))
#  h4a<- read_rep(here("mods", "h4asymp", "nh_R.rep"))
  h4a<- read_rep(here("mods", "h4a", "nh_R.rep"))
  h4 <- read_rep(here("mods", "h4", "nh_R.rep"))
  h5 <- read_rep(here("mods", "h5", "nh_R.rep"))
  h7 <- read_rep(here("mods", "h7", "nh_R.rep"))
  h7tvs <- read_rep(here("mods", "h7tvs", "nh_R.rep"))
  minus1 <- read_rep(here("mods", "minus1", "nh_R.rep"))
  minus0 <- read_rep(here("mods", "minus0", "nh_R.rep"))
  h9 <- read_rep(here("mods", "h9", "nh_R.rep"))

  df <- rbind(
    data.frame(SSB = bc$SSB, R = bc$Pred_Rec, Model = "Base Case"),
    data.frame(SSB = m1$SSB, R = m1$Pred_Rec, Model = "Model 1"),
    data.frame(SSB = m2$SSB, R = m2$Pred_Rec, Model = "Model 2"),
    data.frame(SSB = m3$SSB, R = m3$Pred_Rec, Model = "Model 3"),
    data.frame(SSB = m4$SSB, R = m4$Pred_Rec, Model = "Model 4"),
    data.frame(SSB = m5$SSB, R = m5$Pred_Rec, Model = "Model 5"),
    data.frame(SSB = m6$SSB, R = m6$Pred_Rec, Model = "Model 6")
  )

  df <- rbind(
    #data.frame(SSB = bc$SSB, R = bc$Pred_Rec, Model = "Base Case"),
    #data.frame(SSB = h4$SSB, R = h4$Pred_Rec, Model = "Steepness=0.4 "),
    data.frame(SSB = h5$SSB, R = h5$Pred_Rec, Model = "Steepness=0.5 "),
    data.frame(SSB = h7$SSB, R = h7$Pred_Rec, Model = "Steepness=0.7 "),
    data.frame(SSB = h7tvs$SSB, R = h7tvs$Pred_Rec, Model = "Steepness=0.7, TVS"),
    data.frame(SSB = minus1$SSB, R = minus1$Pred_Rec, Model = "Steepness=0.7, minus1"),
    data.frame(SSB = h9$SSB, R = h9$Pred_Rec, Model = "Steepness=0.9 ")
  )

  #-- Read in the output files for diagnostics and error bars---
{
  mods <-  rbind(
    read_csv(here("mods", "bc", "nh_out.csv")) |> mutate(Model = "Base Case"),
    read_csv(here("mods", "m1", "nh_out.csv")) |> mutate(Model = "Model 1"),
    read_csv(here("mods", "m2", "nh_out.csv")) |> mutate(Model = "Model 2"),
    read_csv(here("mods", "m3", "nh_out.csv")) |> mutate(Model = "Model 3"),
    #read_csv(here("mods", "m4", "nh_out.csv")) |> mutate(Model = "Model 4"),
    #read_csv(here("mods", "m5", "nh_out.csv")) |> mutate(Model = "Model 5"),
    read_csv(here("mods", "m6", "nh_out.csv")) |> mutate(Model = "Model 6")
    #read_csv(here("mods", "Mod1", "nh_out.csv")) |> mutate(Model = "Base Case"),
    #read_csv(here("mods", "h4", "nh_out.csv")) |> mutate(Model = "Steepness=0.4 "),
    #read_csv(here("mods", "h4a", "nh_out.csv")) |> mutate(Model = "Steepness=0.4a "),
    #read_csv(here("mods", "h5", "nh_out.csv")) |> mutate(Model = "Steepness=0.5 "),
    #read_csv(here("mods", "h7", "nh_out.csv")) |> mutate(Model = "Steepness=0.7 "),
    #read_csv(here("mods", "h7tvs", "nh_out.csv")) |> mutate(Model = "Steepness=0.7, TVS "),
    #read_csv(here("mods", "minus0", "nh_out.csv")) |> mutate(Model = "Steepness=0.7, minus0 "),
    #read_csv(here("mods", "minus1", "nh_out.csv")) |> mutate(Model = "Steepness=0.7, minus1 ")
    #read_csv(here("mods", "h9", "nh_out.csv")) |> mutate(Model = "Steepness=0.9 ")
  )
}

#--Plot the SRR curves------------------------
p1 <- df %>% ggplot(aes(x = SSB, y = R, color = Model, shape=Model)) +
  geom_point() +
  geom_line();p1
plotly::ggplotly(p1)

#--Plot the SSB and R curves------------------------

#--SSB---
tail(mods)
mods %>%
  filter(Year > 1960, Variable == "SSB") |>
  ggplot(aes(x = Year, y = value, ymin = ymin, ymax = ymax, type = Model, fill = Model)) +
#  facet_wrap(. ~ Model, scales="fixed") +
  geom_ribbon(alpha = .24) +
  ggthemes::theme_few() +
  coord_cartesian(ylim=c(0, 3e3)) +
  geom_line(aes(color=Model )) + ylim(0,NA) +
  geom_point(aes(color=Model, shape=Model),size=1.) +
  ylab("SSB") +
  xlab("Year")

#--Recruits---
mods %>%
  filter(Year > 2010, Variable == "R") |>
  ggplot(aes(x = Year, y = value, ymin = ymin, ymax = ymax, color = Model, fill = Model)) +
  geom_errorbar(width=0.5,position="dodge",alpha = .3) +
  ggthemes::theme_few() +
  geom_bar(width=0.5,stat = "Identity",position="dodge") +
  ylab("R") +
  xlab("Year") #+
  facet_wrap(. ~ Model)

#--SRR
mods %>%
  filter(Variable == "R" | Variable == "SSB") |>
  select(c(1:3, 6)) |>
  pivot_wider(names_from = Variable, values_from = value) |>
  ggplot(aes(x = SSB, y = R, label = Year, color = Model, fill = Model)) +
  geom_point() +
  geom_text(alpha = .3) +
  ggthemes::theme_few() +
  facet_wrap(. ~ Model) +
  geom_line(data = df, aes(y = R, x = SSB, color = Model), inherit.aes = FALSE)

#--Plot the reference points-------------------
unique(mods$Variable)
mods |>
  filter(Variable != "R" & Variable != "SSB") |>
  ggplot(aes(x = Model, y = value, ymin = ymin, ymax = ymax, color = Model, fill = Model)) +
  geom_point(size = 3) +
  ggthemes::theme_few() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_errorbar() +
  ylab("SSB") +
  xlab("") +
  facet_wrap(. ~ Variable, scales = "free")
# q: how do I exclude legend
# a: theme(legend.position="none")

#--Age fits using tidyverse----

PlotAgeFit(x=minus0, title="Minus age group=0",type="survey1",fage=0)
PlotAgeFit(x=minus1, title="Minus age group=1",type="survey1",fage=1,lage=7)
PlotAgeFit(x=minus1, title="Minus age group=1",type="fishery",fage=2,lage=7)
PlotAgeFit(x=minus1, title="Minus age group=1",type="fishery",fage=2,lage=7)
PlotAgeFit(x=h7tvs, title="Minus age group=2",type="fishery",fage=2,lage=7)
PlotAgeFit(x=bc, title="old base case",type="fishery",fage=2,lage=7)
PlotAgeFit(x=h7tvs, title="Minus age group=2",type="survey1",fage=2,lage=7)
PlotAgeFit()
PlotAgeFit(x=h9,title="Steepness=0.9")
PlotAgeFit(x=minus1, type="survey1")
PlotAgeFit(x=h4,title="Steepness=0.4",type="fishery")
PlotAgeFit(x=h4a,title="Steepness=0.4a",type="fishery")
PlotAgeFit(x=h4,title="Steepness=0.4",type="survey1")
PlotAgeFit(x=h4a,title="Steepness=0.4a",type="survey1")
PlotAgeFit(x=h5,title="Steepness=0.5",type="survey1")
PlotAgeFit(x=h7,title="Steepness=0.7",type="survey1")
PlotAgeFit(x=h9,title="Steepness=0.9",type="survey1")

#--Index fits----
  dfidx <- rbind(
    data.frame(Year = 1964:2023,Obs = bc$Obs_Survey_1, predicted = bc$Pre_Survey_1, Model = "Base Case"),
    #data.frame(Year = 1964:2023,Obs = h4$Obs_Survey_1, predicted = h4$Pre_Survey_1, Model = "h=0.4"),
    #data.frame(Year = 1964:2023,Obs = h4a$Obs_Survey_1, predicted = h4a$Pre_Survey_1, Model = "h=0.4a"),
    data.frame(Year = 1964:2023,Obs = minus0$Obs_Survey_1, predicted = minus0$Pre_Survey_1, Model = "h=0.7, minus0"),
    data.frame(Year = 1964:2023,Obs = minus1$Obs_Survey_1, predicted = minus1$Pre_Survey_1, Model = "h=0.7, minus1"),
    #data.frame(Year = 1964:2023,Obs = h5$Obs_Survey_1, predicted = h5$Pre_Survey_1, Model = "h=0.5"),
    data.frame(Year = 1964:2023,Obs = h7tvs$Obs_Survey_1, predicted = h7tvs$Pre_Survey_1, Model = "h=0.7, TVS")
    #data.frame(Year = 1964:2023,Obs = h9$Obs_Survey_1, predicted = h9$Pre_Survey_1, Model = "h=0.9")
  )
  dfidx |> filter(Obs>-0) |> ggplot(aes(x = Year, y = Obs, color = Model)) +
    geom_point(color="black") + ylim(0,NA)+
    geom_line(aes(y = predicted)) + geom_point(aes(y=predicted))

  dfcpue <- rbind(
    data.frame(Year = 1964:2023,Obs = bc$Obs_CPUE_3, predicted = bc$e_CPUE_3, Model = "Base Case"),
    data.frame(Year = 1964:2023,Obs = m1$Obs_CPUE_3, predicted = m1$e_CPUE_3, Model = "Model 1"),
    data.frame(Year = 1964:2023,Obs = m2$Obs_CPUE_3, predicted = m2$e_CPUE_3, Model = "Model 2"),
    data.frame(Year = 1964:2023,Obs = m3$Obs_CPUE_3, predicted = m3$e_CPUE_3, Model = "Model 3"),
    data.frame(Year = 1964:2023,Obs = m4$Obs_CPUE_3, predicted = m4$e_CPUE_3, Model = "Model 4"),
    data.frame(Year = 1964:2023,Obs = m6$Obs_CPUE_3, predicted = m6$e_CPUE_3, Model = "Model 6")
  #  data.frame(Year = 1964:2023,Obs = h4$Obs_CPUE_3, predicted = h4$e_CPUE_3, Model = "h=0.4"),
    #data.frame(Year = 1964:2023,Obs = h4a$Obs_CPUE_3, predicted = h4a$e_CPUE_3, Model = "h=0.4a"),
    #data.frame(Year = 1964:2023,Obs = h7tvs$Obs_CPUE_3, predicted = h7tvs$e_CPUE_3, Model = "h=0.7, minus2"),
    #data.frame(Year = 1964:2023,Obs = minus0$Obs_CPUE_3, predicted = minus0$e_CPUE_3, Model = "h=0.7, minus0"),
    #data.frame(Year = 1964:2023,Obs = minus1$Obs_CPUE_3, predicted = minus1$e_CPUE_3, Model = "h=0.7, minus1")
    #data.frame(Year = 1964:2023,Obs = h5$Obs_CPUE_3, predicted = h5$e_CPUE_3, Model = "h=0.5"),
    #data.frame(Year = 1964:2023,Obs = h7$Obs_CPUE_3, predicted = h7$e_CPUE_3, Model = "h=0.7"),
    #data.frame(Year = 1964:2023,Obs = h9$Obs_CPUE_3, predicted = h9$e_CPUE_3, Model = "h=0.9")
  )
  dfcpue |> filter(Obs>0) |> ggplot(aes(x = Year, y = Obs, color = Model)) +
    geom_point(color="black") +
    geom_line(aes(y = predicted)) + geom_point(aes(y=predicted))+
    ggtitle("CPUE 3") + ylim(0,NA)

  dfcpue <- rbind(
    data.frame(Year = 1964:2023,Obs = bc$Obs_CPUE_6, predicted = bc$e_CPUE_6, Model = "Basecase     "),
    data.frame(Year = 1964:2023,Obs = h7tvs$Obs_CPUE_6, predicted = h7tvs$e_CPUE_6, Model = "h=0.7, minus2"),
    data.frame(Year = 1964:2023,Obs = minus0$Obs_CPUE_6, predicted = minus0$e_CPUE_6, Model = "h=0.7, minus0"),
    data.frame(Year = 1964:2023,Obs = minus1$Obs_CPUE_6, predicted = minus1$e_CPUE_6, Model = "h=0.7, minus1")
  )
  dfcpue |> filter(Obs>0) |> ggplot(aes(x = Year, y = Obs, color = Model)) +
    geom_point(color="black") +
    geom_line(aes(y = predicted)) + geom_point(aes(y=predicted))+
    ggtitle("CPUE 6") + ylim(0,NA)

  dfcpue <- rbind(
    data.frame(Year = 1964:2023,Obs = bc$Obs_CPUE_1, predicted = bc$e_CPUE_1, Model = "Basecase     "),
    data.frame(Year = 1964:2023,Obs = h7tvs$Obs_CPUE_1, predicted = h7tvs$e_CPUE_1, Model = "h=0.7, minus2"),
    data.frame(Year = 1964:2023,Obs = minus0$Obs_CPUE_1, predicted = minus0$e_CPUE_1, Model = "h=0.7, minus0"),
    data.frame(Year = 1964:2023,Obs = minus1$Obs_CPUE_1, predicted = minus1$e_CPUE_1, Model = "h=0.7, minus1")
  )
  dfcpue |> filter(Obs>0) |> ggplot(aes(x = Year, y = Obs, color = Model)) +
    geom_point(color="black") +
    geom_line(aes(y = predicted)) + geom_point(aes(y=predicted))+
    ggtitle("CPUE 1") + ylim(0,NA)

#---Selectivities----------------
dfsel <- rbind(
    data.frame(Year = 1964:2023,Sel = bc$S[1:60,],  Model = "Base Case"),
    #data.frame(Year = 1964:2023,Sel = h4$S[1:60,],  Model = "h=0.4"),
    data.frame(Year = 1964:2023,Sel = h4a$S[1:60,], Model = "h=0.4a")
    data.frame(Year = 1964:2023,Sel = h7tvs$S[1:60,], Model = "h=0.7, TVS")
  )
names(dfsel) <- c("Year", 0:8, "Model")
glimpse(dfsel)
M<- data.frame(Year = 1964:2023,Sel = h4a$S[1:60,], Model = "h=0.4a")
M<- data.frame(Year = 1964:2023,Sel = bc$S[1:60,], Model = "Basecase")
M<- data.frame(Year = 1964:2023,Sel = m4$S[1:60,], Model = "Model 4, dome-shaped fishery")
M<- data.frame(Year = 1964:2023,Sel = h7tvs$S[1:60,], Model = "h=0.7, TVS")
plot_sel()
sel
M[,2:10]
  Year=M$Year;sel=M[,2:10];styr=1964; fage=NULL; lage=NULL; alpha=0.2;scale=3.8;fill="purple"
  plot_sel <- function(Year=M$Year,sel=M[,2:10],styr=1964, fage=NULL, lage=NULL, alpha=0.2,scale=3.8,fill="purple")
  {
    df        <- data.frame(Year=Year,sel=sel );
    if (is.null(fage)) fage      <- 0
    if (is.null(lage)) lage      <- length(sel[1,])-1
    df <- df |> select(1:(lage-fage+2))
    names(df) <- c("Year",fage:lage)
    nages     <- length(fage:lage)
    names(df)
    sdf       <-  pivot_longer(df,names_to="age",values_to="sel",cols=2:(nages+1)) %>% filter(Year>=styr) %>% mutate(age=as.numeric(age)) #+ arrange(age,yr)
    p1  <- ggplot(sdf,aes(x=age,y=as.factor(Year),height = sel)) + geom_density_ridges(stat = "identity",scale = scale, alpha = alpha,
                                                                                       fill=fill,color="black") + ggthemes::theme_few() +
      ylab("Year") + xlab("Age (years)") +
      scale_x_continuous(limits=c(fage,lage),breaks=fage:lage) +
      scale_y_discrete(limits=rev(levels(as.factor(sdf$Year))))
    return(p1)
  }
  p1
  plot_sel()

#--Try Bayesian posterior integration using adnuts--------
  library(adnuts)                         # needs to be 1.0.9000
  library(snowfall)
  library(rstan)
  library(shinystan)
  reps <- parallel::detectCores()-2 # chains to run in parallel
  reps
  ## Reproducible seeds are passed to ADMB
  set.seed(352)
  seeds <- sample(1:1e4, size=reps)
  seeds

  m <- 'nh -steepness 0.7'
  d <- 'mcmc'
  setwd(d)
  system(paste(m, '-nox -binp nh.bar -hbf 1 -iprint 200 -mcmc 10'))
  setwd('../..')
  setwd('..')
  getwd()

  get.inits <- function(fit, chains){
    post <- extract_samples(fit)
    ind <- sample(1:nrow(post), size=chains)
    lapply(ind, function(i) as.numeric(post[i,]))
  }

  ## Here we assume the pm.exe model is in a folder called 'pm'
  ## as well. This folder gets copied during parallel runs.
  d<- 'mods/mcmc'
  m <- 'nh' # the model name, folder is also assumed to be called this runs/base/mcmc
  ## First optimize the model to make sure the Hessian is good.
  setwd(d);
  system('nh -nox -mcmc 15 -hbf 1 -binp nh.bar -phase 50');
  setwd('../..')


  ## Run--
  iter <- 500 # maybe too many...depends are number cores...I used 8...
  chains=8
  #iter <- 4000*thin; warmup <- iter/#8

  inits <- NULL ## start chains from MLE
  fit.mle <- sample_nuts(model=m, path=d, iter=iter, warmup=iter/4,
                         chains=chains, cores=chains, control=list(control=list(max_treedepth=14,
                                                                                metric='mle')))

  summary(fit.mle)
  plot_uncertainties(fit.mle)
  pairs_admb(fit.mle, pars=1:6, order='slow')
  pairs_admb(fit.mle, pars=1:6, order='fast')

  pdf("pairs_rdev.pdf")
  pairs_admb(fit.mle, pars=68:78)
  pairs_admb(fit.mle2, pars=68:78)
  dev.off()
  pdf("marginals.pdf")
  plot_marginals(fit.mle)
  dev.off()
  print(fit.mle)
  plot_sampler_params(fit.mle)
  launch_shinyadmb(fit.mle)

  ## It doesn't really need any fixes so rerun with NUTS. Reoptimize to get
  ## the correct mass matrix for NUTS. Note the -hbf 1 argument. This is a
  ## technical requirement b/c NUTS uses a different set of bounding
  ## functions and thus the mass matrix will be different.
  ## Use default MLE covariance (mass matrix) and short parallel NUTS chains
  ## started from the MLE.

  ## If good, run again for inference using updated mass matrix. Increase
  ## adapt_delta toward 1 if you have divergences (runs will take longer).
  mass <- fit.mle$covar.est # note this is in unbounded parameter space
  inits <- get.inits(fit.mle, reps) ## use inits from pilot run
  reps
  fit.mle2 <- sample_nuts(model=m, path=d, iter=2000, warmup=iter/4,
                          chains=chains, cores=chains, control=list(control=list(max_treedepth=14,
                                                                                 metric=mass,adapt_delta=0.95)))
  plot_sampler_params(fit.mle2)
  launch_shinyadmb(fit.mle)
  launch_shinyadmb(fit.mle2)
  pairs_admb(fit.mle2, pars=1:6, order='slow')
  summary(fit.mle)
  summary(fit.mle2)
  saveRDS(fit.mle, file='fit.mle.RDS')
  saveRDS(fit.mle2, file='fit.mle2.RDS')

  ## Again check for issues of nonconvergence and other standard checks. Then
  ## use for inference.
  ess <- monitor(fit.mle2$samples, warmup=nuts.updated$warmup, print=FALSE)[,'n_eff']
  nuts.updated$ess <- ess
  saveRDS(nuts.updated, file='nuts.updated.RDS')
  launch_shinyadmb(nuts.updated)
  slow <- names(sort(ess))[1:8]
  png('pairs_slow_nuts2.png', width=7, height=5, units='in', res=500)
  pairs_admb(fit=nuts.updated, pars=slow)
  dev.off()

  save.image() # save to file for later


#--Read in model output files parallel-----
dir_list  <- c("Mod1", "h4", "h4asymp", "h5", "h7", "h9")
fn <- paste0("mods/", dir_list, "/nh_R.rep")
read_rep <- read_admb
fn
mod_list  <- c("bc", "h4", "h4asymp", "h5", "h7", "h9")
mod_names <- c("Base case", "h=0.4", "h=0.4, asymptotic", "h=0.5", "h=0.7", "h=0.9")

#---FishLife steepness--------
library(FishLife)
Taxa = Search_species( Genus = "Merluccius", Species = "capensis")$match_taxonomy
Taxa = Search_species( Genus = "Merluccius", Species = "paradoxus")$match_taxonomy
Predict = Plot_taxa( Taxa,  mfrow=c(3,2) )
