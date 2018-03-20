library(doParallel)
cl <- makeCluster(2)
registerDoParallel(cl)

# standard %dopar% loops are _not_ reproducible
set.seed(1234)
s1 <- foreach(i=1:4) %dopar% { runif(1) }
set.seed(1234)
s2 <- foreach(i=1:4) %dopar% { runif(1) }
identical(s1, s2)

# single %dorng% loops are reproducible
r1 <- foreach(i=1:4, .options.RNG=1234, .packages = 'doRNG') %dorng% { runif(1) }
r2 <- foreach(i=1:4, .options.RNG=1234, .packages = 'doRNG') %dorng% { runif(1) }
identical(r1, r2)
# the sequence os RNG seed is stored as an attribute
attr(r1, 'rng')

# stop cluster
stopCluster(cl)

# More examples can be found in demo `doRNG`
## Not run: ------------------------------------
# demo('doRNG')
## ---------------------------------------------