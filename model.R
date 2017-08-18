require(TMB)
compile("betabinomial.cpp")
dyn.load(dynlib("betabinomial"))

N <- 1000
x <- rnorm(N)
eta <- .3 + 1*x
p <- plogis(eta)
gamma <- 1.1
n <- 5
alpha <- p*(n - 1)/(gamma - 1)
beta <- (1-p)*(n - 1)/(gamma - 1)
y <- rbinom(N, size=n, prob = rbeta(N, alpha,beta))

data <- data.frame(x,y,n)
parameters <- list(a=0,b=0,log_gamma=0)
obj <- MakeADFun(data,parameters, DLL="betabinomial")
opt <- do.call(optim,obj)
summary(sdreport(obj))
obj$fn(opt$par)

summary(glm(cbind(y, n - y) ~ x, quasibinomial, data))



warnings()
