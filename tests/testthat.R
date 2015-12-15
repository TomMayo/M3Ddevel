Sys.setenv("R_TESTS" = "")
library(testthat)
library(M3Ddevel)

test_check("M3Ddevel")
