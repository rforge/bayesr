library(bamlss)
library(funData)
library(tidyverse)
source("simJM.R")
sim_dat <- simJM(nsub = 10, long_setting = "fpc", alpha_setting = "constant",
                 long_df = 3, sigma = 0.05, seed = 120, eval_scale = 2)

# The data are generated so that in Yao (2007) notation
# Formula 6: h_i(t) = h_0(t)*exp(\gamma X_i(t) + V_i(t)'\zeta)
# 
# with
# t = sim_dat$obstime
# h_0(t) = 1.4*log((time + 10)/1000) - 1.5 ---- function lambda()
# \gamma = 1 ---- function alpha()
# V_i(t) = sim_dat$x1 ---- function gamma()
# \zeta = 0.3 ---- function gamma()
# 
# Formula 1: X_i(t_ij) = \mu(t_ij) + Z_i(t_ij)'\beta + 
#               \sum_{k=1}^K \xi_ik\phi_k(t_ij) + \epsilon_ij ---- function mu()
# \mu_ij = 1.25
# \Z_i(t_ij)' = (sin(sim_dat$x2), obstime)
# \beta = (0.6, -0.01)
# K = 3
# \xi_ik ~ N(0, \lambda_k) with 
# \lambda_1 = 2*1, \lambda_2 = 2*2/3, \lambda_3 = 2*1/3 ---- function gen_fpc()
# ---- eVal()
# \phi_1, \phi_2, \phi_3 --- function mu() ---- eFun() ---- funData:::efPoly()
# \epsilon_ij ~ N(0, 0.05)
 
# Short data frame for events
sim_dat_short <- sim_dat %>%
  group_by(id) %>%
  slice(tail(row_number(), 1)) %>%
  ungroup()

# Plot the data
ggplot(data = sim_dat, aes(x = obstime, y = y, color = id)) + 
  geom_line() +
  geom_point(data = sim_dat_short, 
             mapping = aes(x = survtime, y = y, shape = factor(event)),
             size = 2) +
  scale_shape_manual(values = c(1, 4))
