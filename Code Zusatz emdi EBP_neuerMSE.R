ebp_neue = function(formel, universe, sample, smp_domains, uni_domains, L, threshold, 
                    B, MSE, par_bound=c(-1,2), sel = 4, method)
{
  
set.seed(100) 
 
# Generalities ------------------------------------------------------------
  library(nlme)#for the mixed model
  library(FNN) #for Kulback Leibler
  library(MuMIn)#for rsq
  library(nloptr) # for johnson optimisation
  
# PovertyFunctions ---------------------------------------------------------
  
  hcr_function = function(y){
    mean(y<threshold)
  }
  pgap_function = function(y){
    mean((y<threshold)*(threshold-y)/threshold)
  }
  qsr_function = function(y){
    sum(y[(y>quantile(y,0.8))])/sum(y[(y<quantile(y,0.2))])
  }
  fast_mean = function(X)
  {
    .Internal(mean(X))
  }
  gini_function = function(X)
  {
    n <- length(X)
    X <- sort(X)
    G <- sum(X*1L:n)
    G <- 2*G/sum(X) - (n + 1L)
    return(G/n)
  }
  
# Transformation related functions ----------------------------------------
  
  geometric.mean<-function (x) #for RMLE in the parameter estimation
  {
    exp(fast_mean(log(x)))
  }
  
  select_john <- function(z, y){
    Pxis <- pnorm(q = z * c(3, 1, -1, -3))  
    Xs <- quantile(x = y , probs = Pxis)
    dfs <- - diff(Xs)
    m <- dfs[1]
    p <- dfs[2]
    n <- dfs[3]
    critical <- m*n/p^2
    cut(critical,breaks =  c(-Inf,1 - 1e-2, 1 + 1e-2, Inf), 
        labels =c("Sb", "Sl", "Su" ))
  }
  
  skewness = function (x, na.rm = FALSE) #from moments package
  {
    if (is.matrix(x))
      apply(x, 2, skewness, na.rm = na.rm)
    else if (is.vector(x)) {
      if (na.rm)
        x <- x[!is.na(x)]
      n <- length(x)
      (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
    }
    else if (is.data.frame(x))
      sapply(x, skewness, na.rm = na.rm)
    else skewness(as.vector(x), na.rm = na.rm)
  }
  kurtosis <- function (x, na.rm = FALSE) 
  {
    if (is.matrix(x)) 
      apply(x, 2, kurtosis, na.rm = na.rm)
    else if (is.vector(x)) {
      if (na.rm) 
        x <- x[!is.na(x)]
      n <- length(x)
      n * sum((x - mean(x))^4)/(sum((x - mean(x))^2)^2)
    }
    else if (is.data.frame(x)) 
      sapply(x, kurtosis, na.rm = na.rm)
    else kurtosis(as.vector(x), na.rm = na.rm)
  }
  
  bound_estimation = function(y , max_range_bounds = NULL, m=NULL) # bounds for the estimation parameter
  {
    
    if(sel == 1) #for log shift enshure that the data are positive. The other transformations use the given interval
    {
      span=range(y)
      if( (span[1]+1) <= 1)
      {
        lower = abs(span[1])+1
      }
      else
      {
        lower = 0
      }
      upper = diff(span)
      
      return(c(lower,upper))
    }
    # else if(sel == 10)
    # {
    #   if(m == "Su")
    #   {
    #     return(max_range_bounds)
    #   }
    #   else if(m == "Sl"
    #   {
    #     
    #   }
    #   else if(m == "Sb")
    #   {
    #     
    #   }
    #   else
    #   {
    #     stop("This should not have happened Johnny")
    #   }
    # }
    else
    {
      return(max_range_bounds)
    }
  }
  
# Loading the selected transformation method for y -------------------------------------------------------------
  if(sel==1)
  {
    transformation = function(l,y,inv=F, m=NULL) #Log-shift transformation
    {
      m=0
      if(!inv)
      {
        y=log(y+l)
      }
      else
      {
        y=exp(y)-l
      }
      return(list(y=y, m=m))
    }
  }
  else if(sel==2)
  {
    transformation = function(l, y, inv=FALSE, m=NULL) #no transformation
    {
      return(list(y=y, m=m))
    }
  }
  else if(sel==3)
  {
    transformation = function(l,y, inv=F,m=NULL){ #convex-to-concave transformation
      m=0
      if(!inv){
        sm=y<0
        if(any(!sm))
        {
          if(abs(l)<=1e-12)
          {
            y[!sm]=log(y[!sm]+1)
          }
          else
          {
            y[!sm]=((y[!sm]+1)^l-1)/l
          }
        }
        if(any(sm))
        {
          if(abs(l-2)<=1e-12)
          {
            y[sm]=-log(-y[sm]+1)
          }
          else
          {
            y[sm]=-((-y[sm]+1)^(2-l)-1)/(2-l)
          }
        }
      }
      
      else{
        if(abs(l)<=1e-12)
        {
          sm = y >= 0-
            if(any(sm))
            {
              y[sm]=exp(y[sm]) - 1
            }
          if(any(!sm))
          {
            y[!sm] = 1 - exp(-y[!sm])
          }
        }
        else
        {
          if(l<0)
          {
            sm = y < -1/l
          }
          else
          {
            sm  = y > -1/l
          }
          if(l<2)
          {
            sm2 = y < -1/(l-2) & !sm
          }
          else
          {
            sm2 = y > -1/(l-2)  & !sm
          }
          if(any(sm2))
          {
            y[sm2] = -((y[sm2]*(l-2)+1)^(1/(2-l)) - 1)
          }
          if(any(sm))
          {
            y[sm] = (y[sm] * l + 1)^(1/l) - 1
          }
        }
      }
      return(list(y=y, m=m))
    }
  }
  else if(sel==4)
  {
    transformation = function(l, y, inv=FALSE, m=NULL) #Box-Cox transformation (lambda=l)
    {
      if(!inv)
      {
        if(is.null(m))
        {
          m = 0
        }
        if((s=min(y))<=0) #with shift(=m) parameter for making data positive (>0)
        {
          s = abs(s)+1
        }
        else
        {
          s=0
        }
        m=m+s
        
        if(abs(l)<=1e-12) #case lambda=0
        {
          y = log(y+m)
        }
        else
        {
          y = ((y+m)^l-1)/l
        }
      }
      else
      {
        if(is.null(m)) #inverse transformation
        {
          m = 0
        }
        if(abs(l)<=1e-12) #case lambda=0
        {
          y = exp(y) - m
        }
        else
        {
          y = (l*y+1)^(1/l)-m
        }
      }
      return(list(y = y, m = m)) #return of transformed data and shift (overwriten y)
    }
  }
  else if(sel==5)
  {
    transformation = function(l, y, inv=FALSE, m) #Log transformation
    {
      if(!inv)
      {
        m = 0
        if(min(y)<=0)
        {
          m = abs(min(y)) + 1
          y = y + m
        }
        y=log(y)
      }
      else
      {
        y=exp(y)-m
      }
      return(list(y=y, m=m))
    }
  }
  else if(sel==6)
  {
    transformation = function(l, y, inv=FALSE, m) #signed Log transformation
    {
      m <- 0
      if(!inv)
      {
        y <- log(abs(y)^sign(y))
      }
      else
      {
        y <- sign(y) * exp(sign(y)*y)
      }
      return(list(y=y, m=m))
    }
  }
  else if(sel==7)
  {
    transformation = function(l, y, inv=FALSE, m) #exponetial
    {
      m <- 0
      if(!inv)
      {
        if(abs(l)<=1e-12)
        {
          ;
        }
        else
        {
          y <- (exp(l*y)-1)/l
        }
      }
      else
      {
        if(abs(l)<=1e-12)
        {
          ;
        }
        else
        {
          y <- log(l*y+1)/l
        }
      }
      return(list(y=y, m=m))
    }
  }
  else if(sel == 8)#modulus transformation
  {
    transformation = function(l, y, inv=FALSE, m){
      m <- 0
      if(!inv)
      {
        if(abs(l)<= 1e-12)
        {
          y <- sign(y) * log(sign(y)*y + 1)
        }
        else
        {
          y <- sign(y) * (((abs(y) + 1 )^l ) -1 ) / l
        }
      }
      else
      {
        if(abs(l)<= 1e-12)
        {
          y <- sign(y) * (exp(abs(y)) - 1)
        }
        else
        {
          y <- sign(y) * ((abs(y)*l + 1)^(1/l) - 1) - 1
        }
      }
      return(list(y=y, m=m))
    }
  }
  else if (sel == 9)
  {
    transformation = function(l, y, inv=FALSE, m)
    {
      if(!inv)
      {
        m <- 0
        if(min(y) <= 0)
        {
          m = abs(min(y)) + 1
          y = y + m
        }
        if(abs(l)<= 1e-12)
        {
          y <- log(y)
        }
        else
        {
          y <- (y^l - y^(-l)) / (2 * l)
        }
      }
      else
      {
        if(abs(l)<= 1e-12)
        {
          y <- exp(y) - m
        }
        else
        {
          y <- (l * y + sqrt(l^2 * y^2 + 1))^(1/l) - m 
        }
      }
      return(list(y=y, m=m))
    }
  }
  else if(sel == 10)#johnson 2parameter
  {
    transformation <- function(l,y,inv=F, m="Su"){
      if(m=="Su")
      {
        if(!inv){
          y <- asinh((y - l[1])/l[2])
          if(any(is.na(y)))
          {
            browser()
          }
        }
        else
        {
          y <- l[1] + (l[2] * sinh(y))
        }
      }
      else if(m=="Sb")
      {
        if(!inv){
          y <- log((y-l[1])/(l[2]+l[1]-y))
        }
        else
        {
          y <-  l[1]+l[2] - l[2]/(exp(y)+1) 
        }
        return(list(y=y, m=m))
      }
      else if(m=="Sl")
      {
        if(!inv){
          y <- log((y - l[1])/l[2])
        }
        else
        {
          y <- (exp(y) * l[2]) + l[1] 
        }
      }
      return(list(y=y, m=m))
    }
  }
  else if(sel == 11)
  {
    transformation = function(l, y, inv=FALSE, m)
    {
      m <-0
      if(!inv)
      {
        y <- (sign(y) * abs(y)^l -1)/l
      }
      else
      {
        tmp <- y / l + 1
        y <- sign(tmp) * abs(tmp)^(1/l)
      }
      return(list(y=y, m=m))
    }
  }
  else if (sel == 12)#arcsinh
  {
    transformation = function(l, y, inv=FALSE, m)
    {
      m <- NULL
      if (!inv)
      {
        y <- sinh(l[1] * asinh(y) - l[2])
      }
      else
      {
        y <- sinh((asinh(y)+l[2])/l[1])
      } 
      return(list(y=y, m=m))
    }
  }
  else if (sel == 13) # johnson 4 parameter
  {
    transformation <- function(l,y,inv=F, m="Su"){
      if(m=="Su")
      {
        if(!inv){
          y <- l[1] + l[2] * asinh((y - l[3])/l[4])
        }
        else
        {
          y <- l[3] + (l[4] * sinh((y - l[1])/l[2]))
        }
      }
      else if(m=="Sb")
      {
        if(!inv){
          y <- l[1]  +  l[2]  * log((y-l[3])/(l[4]+l[3]-y))
        }
        else
        {
          y <-  l[3]+l[4] - (l[4]/(exp((y - l[1])/l[2])+1)) 
        }
        return(list(y=y, m=m))
      }
      else if(m=="Sl")
      {
        if(!inv){
          y <- l[1] + l[2] * log((y - l[3])/l[4])
        }
        else
        {
          y <- (exp((y - l[1])/l[2]) * l[4]) + l[3] 
        }
      }
      return(list(y=y, m=m))
    }
  }
  
# Loading the selected standardized transformation method -----------------------------------
  
  
  if(sel == 1)
  {
    sd_trans = function(y, l, m) # standardized Log-shift
    {
      y = y + l
      gm = geometric.mean(y)
      y = gm*log(y)
    }
  }
  else if(sel == 2)  # non transformation
  {
    sd_trans = function(y, l, m)
    {
      return(NULL)
    }
  }
  else if(sel == 3)  # standardized convex-to-concave transformation
  {
    sd_trans = function(y, l, m)
    {
      sm=y<0
      if(any(!sm))
      {
        if(abs(l)<=1e-12)
        {
          y[!sm] = y[!sm] + 1
          gm = geometric.mean(y[!sm])
          y[!sm] = gm*log(y[!sm])
        }
        else
        {
          gm = geometric.mean( y[!sm] + 1)
          y[!sm]=((y[!sm]+1)^l-1)/(l*((gm)^(l-1)))
        }
      }
      if(any(sm))
      {
        if(abs(l-2)<=1e-12)
        {
          gm = geometric.mean( 1 - y[sm])
          y[sm]=-log(-y[sm]+1)*gm
        }
        else
        {
          gm = geometric.mean( 1 - y[sm])
          y[sm]=(-((-y[sm]+1)^(2-l)-1)/(2-l))*gm^(l-1)
        }
      }
      return(y)
    }
  }
  else if(sel==4) # standardized Box-Cox
  {
    sd_trans = function(y, l, m)
    {
      if((m=min(y))<=0)
      {
        y=y-m
        y=y+1
      }
      
      gm=geometric.mean(y)
      if(abs(l)>1e-12)
      {
        y=(y^l-1)/(l*((gm)^(l-1)))
      }
      else
      {
        y=gm*log(y)
      }
      return(y)
    }
    
  }
  else if(sel == 5)  # for log transformation
  {
    sd_trans = function(y, l, m)
    {
      return(NULL)
    }
  }
  else if(sel == 6)
  {
    sd_trans = function(y, l, m)
    {
      geo_p <- geometric.mean(y[y>1])
      geo_n <- geometric.mean(abs(y[y<-1]))
      geo <- rep(0, length(y))
      geo[y>1] <- geo_p
      geo[y<-1] <- geo_n
      y<- log(abs(y)^sign(y))*geo
      return(y)
    }
  }
  else if(sel == 7)
  {
    sd_trans = function(y, l, m)
    {
      y <- exp(-l * fast_mean(y)) * exp(l * y - 1) / l
      return(y)
    }
  }
  else if(sel == 8)
  {
    sd_trans <- function(y, l, m)
    {
      if(abs(l) <= 1e-12)
      {
        y <- sign(y) * geometric.mean(y) * (abs(y)+1) * log(abs(y)+1)
      }
      else
      {
        y <- sign(y) * ((abs(y) +1)^l -1)/(l * geometric.mean(y)^(l-1) * (abs(y)+1))
      }
    }
  }
  else if(sel == 9)
  {
    sd_trans <- function(y, l, m)
    {
      if(min(y) <= 0)
      {
        m = abs(min(y)) + 1
        y = y + m
      }
      if(abs(l) <= 1e-12)
      {
        y <- geometric.mean(y) * log(y)
      }
      else
      {
        geo <- geometric.mean(y^(l-1)+y^(-l-1))
        y <- ((y^l -y^(-l))/(2*l))*2/geo
      }
    }
  }
  else if(sel == 10)#johnson 2 parameter
  {
    sd_trans <- function(y, l, m)
    {
      if(m == "Su")
      {
        y <- geometric.mean(sqrt(1+((y - l[1])/l[2])^2))/l[2] * asinh((y-l[1])/l[2]) 
      }
      else if(m == "Sb")
      {
        y <- geometric.mean((y - l[1]) * (l[2] + l[1] - y)) / l[2] *
          log( (y - l[1]) / (l[2] + l[1] - y) )
      }
      else if(m == "Sl")
      {
        y <-  geometric.mean(y-l[1])  *
          log((y - l[1])/l[2])
      }
      return(y)
    }
  }
  else if(sel == 11)
  {
    sd_trans = function(y, l, m)
    {
      geo_p <- geometric.mean(y[y>1])
      geo_n <- geometric.mean(abs(y[y<-1]))
      geo <- rep(0, length(y))
      geo[y>1] <- geo_p
      geo[y<-1] <- geo_n
      y<- sign(y) * (abs(y)^l -1) / (l * geo^(l-1))
      return(y)
    }
  }
  else if(sel == 12)#asinh trafo
  {
    sd_trans <- function(y, l, m = NULL)
    {
      
      hlp <- cosh(l[1] * asinh(y) - l[2]) / sqrt(y^2 + 1)
      y <-  geometric.mean(hlp) / l[1] * sinh( l[1] * asinh(y) - l[2])
      return(y)
    }
  }
  else if (sel == 13) # johnson 4 parameter
  {
    sd_trans <- function(y, l, m)
    {
      if(m == "Su")
      {
        y <- geometric.mean(sqrt(1+((y - l[1])/l[2])^2)) * (l[2]/l[4]) * (l[1] + l[2]*asinh((y-l[3])/l[4])) 
      }
      else if(m == "Sb")
      {
        y <- geometric.mean((y - l[3]) * (l[4] + l[3] - y)) / ( l[4] * l[2]) *
          (l[1] + l[2]*log((y - l[3]) / (l[4] + l[3] - y)))
      }
      else if(m == "Sl")
      {
        y <-  geometric.mean(y-l[3])/l[2]  *
          ( l[1]+l[2] * log((y - l[3])/l[4]))
      }
      return(y)
    }
  }
  
# Loading the generict opt function -----------------------------------
  
  generic_opt=function(l, y, dat, form, idD, meth,m = NULL) # Estimation methos for lambda
  {
    method <- meth
    if(method != 6) #ML (method 6) uses another approach
    {
      y = transformation(l=l,y=y,inv=F, m=m)$y
      mod <- NULL
      try(mod <- lme(y~-1+dat,random=~1|as.factor(idD),method="REML"))
      if(is.null(mod))
      {
        return(999999)
      }
      res <- residuals(mod,level=0,type="pearson")
    }
    if(method==1)
    {
      skres=skewness(res) #pooled skewness minimization (Natalia)
      re=as.matrix(random.effects(mod))[,1]
      skre=skewness(re)
      sigmae2est=mod$sigma^2
      sigmau2est=as.numeric(VarCorr(mod)[1,1])
      w=sigmae2est/(sigmae2est+sigmau2est)
      
      result=w*abs(skres)+(1-w)*abs(skre)
    }
    else if(method == 2) #skewness minimization (Molina)
    {
      skres=skewness(res)
      result = abs(skres)
    }
    else if(method == 3) #Divergence minimization (KS: Kolmogorov Smirnoff)
    {
      step.length=10^-4
      eval.probs=seq(0,1, by=step.length)
      eval.points = qnorm(eval.probs, mean=mean(res), sd=sd(res))
      test.probs = ecdf(res)(eval.points)
      difs=eval.probs - test.probs
      result = max(abs(difs)) #minimization of the supremum of the diferences
    }
    else if(method == 4) #Divergence minimization (Craemer von Mises)
    {
      step.length=10^-4
      eval.probs=seq(0,1, by=step.length)
      eval.points = qnorm(eval.probs, mean=mean(res), sd=sd(res))
      test.probs = ecdf(res)(eval.points)
      difs=eval.probs - test.probs
      result = sum((difs)^2) #minimization of the sum of the squared diferences
    }
    else if(method == 5) #Divergence minimization (KL: Kullback-Leibler from FNN package)
    {
      step.length=10^-4
      eval.probs=seq(0,1, by=step.length)
      eval.probs=eval.probs[-c(1,length(eval.probs))]
      eval.points = qnorm(eval.probs, mean=mean(res), sd=sd(res))
      test.probs = quantile(res, probs = eval.probs)
      result=KL.divergence(eval.probs,test.probs,k=5)
      result=median(result)
    }
    else if(method == 6) #RMLE method
    {
      y = sd_trans(y = y, l = l, m = m)
      if(any(is.na(y))){
        return(999999)
      }
      mod = lme(y ~ -1 + dat, random=~1|as.factor(idD), method="REML")
      if(logLik(mod) > 0)
      {
        if(sel == 10){
          return(Inf)
        }
        else{
          return(999999)
        }
      }
      return(- logLik(mod))
    }
    else if(method == 7) #kurtosis to 3
    {
      result <- abs(kurtosis(res)-3)
    }
    else if(method == 8)
    {
      result <- 1 - r.squaredGLMM(mod)[["R2c"]]
    }
    return(result)
  }
  
# Point Estimation Function ------------------------------------
  point_estim <- function(Y_smp)
  {
    
    if(sel != 2 && sel !=5 && sel != 12 && sel != 10 && sel != 13 && sel !=6) # if no estimation, no optmimization
    {
      par_bound = bound_estimation(Y_smp, max_range_bounds = par_bound)
      
      # Estimation of Transformation parameters
      optpar = optimize(generic_opt, y = Y_smp , dat = X_smp, 
                      form = formel, idD = smp_domains, meth = method
                      , interval = par_bound, maximum = F
      )$minimum
    }
    else if(sel == 12) #arcsinh
    {
      optpar <- optim(par = c(1,1), fn = generic_opt, dat = X_smp, form=formel, 
                      idD=smp_domains, meth=method,y=Y_smp)$par
    }
    else if(sel == 10 | sel == 13 ) #johnson
    {
      
      m <- select_john(z = 1, 
                       y = resid(lme(Y_smp~-1+X_smp,random=~1|as.factor(smp_domains),
                                     method="REML")))
      if(sel == 10)
      {
        par_bound = bound_estimation(Y_smp, max_range_bounds = par_bound, m = m)
        
        optpar <- sbplx(  x0 = c(1,1), # am besten  bisher
                          fn = generic_opt, 
                          upper = c(Inf,Inf), 
                          lower = c(-Inf,0.001),
                          control = list(xtol_rel = 1e-16),
                          nl.info = FALSE,  
                          m = m, dat = X_smp, 
                          form = formel, 
                          idD = smp_domains, 
                          meth = method,
                          y = Y_smp 
        )$par
      }
      else
      {
        optpar <- optim (   par = c(1,1,1,1), fn = generic_opt, 
                            upper = c(Inf,Inf,Inf,Inf), 
                            lower = c(-Inf, 0.01, -Inf, 0.00001),
                            method = "L-BFGS-B",
                            m = m, dat = X_smp, 
                            form = formel, 
                            idD = smp_domains, 
                            meth = method,
                            y = Y_smp 
        )$par
      }
    }
    else
    {
      optpar = NULL
    }
    tmp = transformation(y=Y_smp, l=optpar, inv=F)
    Y_smp = tmp$y
    par_m = tmp$m
    rm(tmp)
    # estimation of the mixed linear model
    nlm_smp = lme(Y_smp~-1+X_smp,random=~1|as.factor(smp_domains),method="REML")
    betas=fixed.effects(nlm_smp)
    
    rand_eff=rep(0, length(unique(uni_domains)))
    rand_eff[dist_obs_dom]=(random.effects(nlm_smp)[[1]])
    sigmae2est=nlm_smp$sigma^2    # Estimated error variance
    sigmau2est=as.numeric(VarCorr(nlm_smp)[1,1]) # VarCorr(fit2) is
    
    gamma = sigmau2est / (sigmau2est + sigmae2est / n_smp)
    sigmav2est = sigmau2est * (1-gamma)
    
    rand_eff_uni = rep(rand_eff, n_uni)
    
    mu = X_uni %*% betas + rand_eff_uni
    
    quant10s = quant25s = quant75s = quant90s = ginis = qsrs = pgaps = hcrs = means = matrix(nrow=N_dom_uni, ncol=L)
    
    # for-loop calculating the pseudo populations
    for(l in 1 : L)
    {
      eps = rnorm(N_uni, 0, sqrt(sigmae2est))
      
      vu=vector(length=N_uni)
      #browser()
      vu[!obs_dom] = rep(rnorm(N_dom_unobs, 0, sqrt(sigmau2est)),n_uni[!dist_obs_dom])
      vu[obs_dom]  = rep(rnorm(rep(1, N_dom_smp), 0, sqrt(sigmav2est)), n_uni[dist_obs_dom])
      y_pred = mu + eps + vu
      y_pred = transformation(y=y_pred, l=optpar,inv=T, m=par_m)$y
      y_pred[!is.finite(y_pred)] = 0
      
      quant10s[,l] = tapply(y_pred, uni_domains, quantile, probs=0.1)
      quant25s[,l] = tapply(y_pred, uni_domains, quantile, probs=0.25)
      quant75s[,l] = tapply(y_pred, uni_domains, quantile, probs=0.75)
      quant90s[,l] = tapply(y_pred, uni_domains, quantile, probs=0.9)
      ginis[,l] = tapply(y_pred, uni_domains, gini_function)
      means[,l] = tapply(y_pred, uni_domains, fast_mean)
      hcrs[,l] = tapply(y_pred, uni_domains, hcr_function)
      qsrs[,l] =  tapply(y_pred, uni_domains, qsr_function)
      pgaps[,l] =  tapply(y_pred, uni_domains, pgap_function)
    }
    
    # point estimations
    results = data.frame(
      dom = unique(uni_domains),
      quant10 = rowMeans(quant10s),
      quant25 = rowMeans(quant25s),
      quant75 = rowMeans(quant75s),
      quant90 = rowMeans(quant90s),
      gini = rowMeans(ginis),
      mean = rowMeans(means),
      hcr = rowMeans(hcrs),
      qsr = rowMeans(qsrs),
      pgap = rowMeans(pgaps)
    )
    
    return(list(POV = results,
                optpar = optpar,
                par_m = par_m,
                mu = mu,
                betas = betas,
                rand_eff = rand_eff,
                sigmae2est = sigmae2est,
                sigmau2est = sigmau2est,
                sigmav2est = sigmav2est,
                mix_Mod = nlm_smp)
    )
  }
  
# MSE Estim Function--------------------------------------------
  mse_estim = function(optpar_,par_m_,mu_uni_, mu_smp_,sigmae2est_,sigmau2est_,sigmav2est_)
  {
    eps   =  vector(length = N_uni )
    eps[obs_dom]  =  rnorm(sum(obs_dom), 0, sqrt(sigmae2est_))
    
    eps[!obs_dom] =  rnorm(sum(!obs_dom), 0, sqrt(sigmae2est_+sigmau2est_))
    
    vu_tmp=rnorm(N_dom_uni, 0, sqrt(sigmau2est_))
    
    vu_uni  = rep(vu_tmp, n_uni)
    vu_smp  = rep(vu_tmp[dist_obs_dom], n_smp)
    
    rm(vu_tmp)
    Y_uni_b = mu_uni_ + eps + vu_uni
    eps = rnorm(N_smp, 0, sqrt(sigmae2est_))
    
    
    Y_smp_b = mu_smp_ + eps + vu_smp
    
    Y_smp_b_del = Y_smp_b
    Y_smp_b = transformation(y=Y_smp_b, l=optpar_, inv=T, m=par_m_)$y
    
    Y_smp_b[!is.finite(Y_smp_b)]=0
    
    alles = point_estim(Y_smp = Y_smp_b)
    estims = as.matrix(alles[[1]][,-1])
    optpar_b = alles[[2]]
    
    Y_uni_b = transformation(y=Y_uni_b, l=optpar_, inv=T, m=par_m_)$y
    Y_uni_b[!is.finite(Y_uni_b)] = 0
    truth = as.matrix(data.frame(
      quant10 = tapply(Y_uni_b, uni_domains, quantile, probs=0.1),
      quant25 = tapply(Y_uni_b, uni_domains, quantile, probs=0.25),
      quant75 = tapply(Y_uni_b, uni_domains, quantile, probs=0.75),
      quant90 = tapply(Y_uni_b, uni_domains, quantile, probs=0.9),
      gini = tapply(Y_uni_b, uni_domains, gini_function),
      mean = tapply(Y_uni_b, uni_domains, fast_mean),
      hcr =  tapply(Y_uni_b, uni_domains, hcr_function),
      qsr = tapply(Y_uni_b, uni_domains, qsr_function),
      pgap = tapply(Y_uni_b, uni_domains, pgap_function)
    ))
    
    return(list((estims-truth)^2,optpar_b))
  }
  
  
# General Variable Definitions ---------------------------------
  rearrange_universe = order(uni_domains)
  rearrange_sample = order(smp_domains)
  
  universe = universe[rearrange_universe, ,drop = FALSE]
  uni_domains = uni_domains[rearrange_universe]
  #sampling_weight = uni_domains[rearrange_universe]
  #rururbv = rururbv[rearrange_universe]
  
  sample = sample[rearrange_sample,,drop = FALSE]
  smp_domains = smp_domains[rearrange_sample]
  
  
  N_uni = length(uni_domains)
  N_smp = length(smp_domains)
  N_unobs = N_uni-N_smp
  N_dom_smp = length(unique(smp_domains))
  N_dom_uni = length(unique(uni_domains))
  N_dom_unobs = N_dom_uni - N_dom_smp
  n_smp = as.vector(table(smp_domains))
  n_uni = as.vector(table(uni_domains))
  
  obs_dom = uni_domains %in% unique(smp_domains)
  dist_obs_dom = sort(unique(uni_domains)) %in% unique(smp_domains)
  tmp = formel
  tmp[2] = formula(uni_domains~1)[2]
  X_uni = model.matrix(tmp, data.frame(uni_domains, universe))
  rm(tmp)
  X_smp = model.matrix(formel,sample)
  Y_smp = as.matrix(sample[paste(formel[2])])
  
# Calling the functions -------------------------------------------------------------
  Pov_ests = point_estim(Y_smp = Y_smp)
  
  mu_smp = X_smp %*% Pov_ests$betas
  mu_uni = X_uni %*% Pov_ests$betas
  
  if(MSE)
  {
    
    mses = replicate(B, mse_estim(
      optpar_         = Pov_ests$optpar ,
      par_m_          = Pov_ests$par_m ,
      mu_uni_         = mu_uni ,
      mu_smp_         = mu_smp , 
      sigmae2est_     = Pov_ests$sigmae2est ,
      sigmau2est_     = Pov_ests$sigmau2est ,
      sigmav2est_     = Pov_ests$sigmav2est
    ),simplify = F)
    
    optpar_bs = sapply(mses, function(X){X[[2]]})
    
    mses = sapply(mses, function(X){X[[1]]}, simplify="array")
    
    mses = apply(mses, c(1,2), fast_mean)
    mses = data.frame(Domains = unique(uni_domains),mses)
    result = list(POV = Pov_ests[[1]], mix_Mod = Pov_ests$mix_Mod, optpar = Pov_ests[[2]], 
                  MSE = mses, Pov_ests[-1], optpar_bs = optpar_bs)
  }
  else
  {
    result = Pov_ests
  }
##### ------
  
  return(result)
}

