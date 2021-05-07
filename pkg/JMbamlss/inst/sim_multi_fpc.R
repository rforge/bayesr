library(tidyverse)
source("simMultiJM.R")
sim_dat <- simMultiJM(nsub = 10, M = 4,
                      sigma = function(t, x) {-3+I(x$marker == "m2")*1},
                      mfpc_args = list(type = "split", eFunType = "PolyHigh",
                                       ignoreDeg = c(1:3), eValType = "linear",
                                       eValScale = 5),
                      seed = 1,
                      full = TRUE)

# Short data frame for events
sim_dat_short <- sim_dat$data %>%
  group_by(id, marker) %>%
  slice(tail(row_number(), 1)) %>%
  ungroup()


# Plot the data
ggplot(data = sim_dat$data, aes(x = obstime, y = y, color = id)) + 
  geom_line() +
  facet_grid(marker ~ ., scales = "free_y") +
  geom_point(data = sim_dat_short, 
             mapping = aes(x = survtime, y = y, shape = factor(event)),
             size = 2) +
  scale_shape_manual(values = c(1, 4))
  
