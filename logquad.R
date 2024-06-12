library(logquad5q0)
library(tidyverse)
library(SUMMER)

input <- format_data(
  rate      = 0.00804138,
  lower_age = 28,
  upper_age = 365.25*5,
  type      = "qx",
  sex       = "total")

lagrange5q0(data = input)