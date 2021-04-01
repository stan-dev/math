
data(mcycle, package="MASS")
setwd("~/Code/laplace_approximation/math/test/unit/math/laplace")

write.table(t(mcycle$times), file = "motorcycle_gp/x_vec.csv",
  row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ")

write.table(t(mcycle$accel), file = "motorcycle_gp/y_vec.csv",
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ")

# write.csv(as.vector(mcycle$times), "motorcycle_gp/x_vec.csv")
# write.csv(mcycle$accel, "motorcyle_gp/y_vec.csv")
