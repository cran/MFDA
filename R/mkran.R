mkrandom<-function (formula, data) 
{
    attach(data)
    form.wk <- terms.formula(formula)[[2]]
    if (!("|" %in% strsplit(deparse(form.wk), "")[[1]])) 
        stop("gss error in mkran: missing | in grouping formula")
    term.wk <- strsplit(deparse(form.wk), " \\| ")[[1]]
    z2.wk <- eval(parse(text = term.wk[2]))  
    z <- NULL
    lvl.z2 <- levels(z2.wk)
    for (i in lvl.z2) z <- cbind(z, as.numeric(z2.wk == i))
    init <- 0
    env <- length(levels(z2.wk))
    fun <- function(zeta, env) diag(10^(-zeta), env)
    sigma <- list(fun = fun, env = env)    
    detach(data)	
    list(z = z, sigma = sigma, init = init)
}