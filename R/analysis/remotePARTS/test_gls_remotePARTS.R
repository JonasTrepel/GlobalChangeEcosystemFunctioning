#### reproducible example remotePARTS 
library(remotePARTS)
source("R/functions/extract_gls_estimates.R")

dat <- data.frame(
  X = runif(10000, min = -180, max = 180),
  Y = runif(10000, min = -90, max = 90),
  t_coef = rnorm(10000), 
  var_1 = rnorm(10000, 2, 3), 
  var_2 = rnorm(10000, 4, 3), 
  var_3 = rnorm(10000, 1, 3),
  var_4 = rnorm(10000, 5, 3),
  var_5 = rnorm(10000, 1.4, 1)
)


#partition the data 

pm <- sample_partitions(npix = nrow(dat), partsize = 1500, npart = NA)
dim(pm)

## 
set.seed(1)
#Intercept only model works fine, same p value and estimates when running it repeatedly 
gls_fit1 <- fitGLS_partition(t_coef ~ 1,
                           partmat = pm,
                           covar_FUN = "covar_exp",
                           covar.pars = list(range = 0.01),
                           data = dat,
                           nugget = NA,
                           ncores = 4,
                           progressbar = TRUE, 
                           parallel = T, 
                           coord.names = c("X", "Y")
)
gls_fit1

set.seed(161)
gls_fit2 <- fitGLS_partition(t_coef ~ 1 +
                               var_1 + var_2 + var_3 + var_4 + var_5,
                             partmat = pm,
                             covar_FUN = "covar_exp",
                             covar.pars = list(range = 0.1),
                             data = dat,
                             nugget = NA,
                             ncores = 4,
                             progressbar = TRUE, 
                             parallel = T, 
                             coord.names = c("X", "Y")
)
gls_fit2

dt_est <- extract_gls_estimates(gls_fit2, part = TRUE)
dt_est %>% 
  ggplot() + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(y = term, x = estimate, xmin = ci_lb, xmax = ci_ub, color = sig),
                  alpha = 0.9, linewidth = 1.2) +
  scale_color_manual(values = c("significant" = "forestgreen", "non-significant" = "grey")) +
  theme_classic() +
  theme(legend.position = "none")
