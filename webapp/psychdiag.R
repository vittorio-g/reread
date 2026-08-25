suppressMessages(library(careless))
m <- as.matrix(read.csv("sim_bench/300_0_M.csv", header=FALSE))
R <- cor(m); ut <- R[upper.tri(R)]
cat("max|r|", round(max(abs(ut)),3), " pos>=.4:", sum(ut>=.4), " pos>=.3:", sum(ut>=.3), "\n")
for(cv in c(.4,.3,.25)){ v<-careless::psychsyn(m,critval=cv); cat("psychsyn@",cv," finite",sum(is.finite(v))," mean",round(mean(v,na.rm=TRUE),3),"\n") }
