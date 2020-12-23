library(lmerTest)
library(variancePartition)

tmp <- dream(t(mtcars$mpg), ~ cyl, data = mtcars, L = c(0 ,1))

#test <- contest(tmp, 1)
