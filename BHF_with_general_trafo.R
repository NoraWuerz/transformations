# 1. Normal Seznario--------
# Simulationen ------------------

Sim_Anzahl <- 500

source("function_BHF_general_trafo.R")

domains  <- "idD"
formel <- y ~ x

# ohne trafo ------------
Est_Tau_ebp_Matrix<- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_bc_Matrix <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_n_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
TrueMeans_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
ResEstPar         <- list()

ohne_trafo <- list()

load(file = "H:/Transformationen/Pop-Data/s1_normal_500.RData")

for(i in 1:Sim_Anzahl){
  print(paste("ohne Trafo", i))
  
  data_smp <- Pop[[i]]
  data_pop <- attr(Pop[[i]], "pop")
    
  Ergebnisse <- bhf_trafo(data_pop = data_pop, data_smp = data_smp,
                          domains = domains, formel = formel, sel = 2)
  
  parameter <- list()
  parameter[[1]] <- Ergebnisse$ResEstPar$optpar
  parameter[[2]] <- Ergebnisse$ResEstPar$betas
  names(parameter) <- c("optpar", "beta")
  
  ResEstPar[[i]]        <- parameter
  Est_Tau_bc_Matrix[,i] <- Ergebnisse$mean_est_bc
  Est_Tau_n_Matrix[,i]  <- Ergebnisse$mean_est_naive
  Est_Tau_ebp_Matrix[,i]<- Ergebnisse$mean_est_ebp
  TrueMeans_Matrix[,i]  <- Ergebnisse$mean_true
}

ohne_trafo[[1]] <- ResEstPar
ohne_trafo[[2]] <- Est_Tau_bc_Matrix
ohne_trafo[[3]] <- Est_Tau_n_Matrix
ohne_trafo[[4]] <- Est_Tau_ebp_Matrix
ohne_trafo[[5]] <- TrueMeans_Matrix
names(ohne_trafo) <- c("ResEstPar", "mean_bc", "mean_naive", "mean_ebp", "mean_true" )

rm(ResEstPar, Est_Tau_bc_Matrix, Est_Tau_n_Matrix, Est_Tau_ebp_Matrix, TrueMeans_Matrix, 
   parameter, Ergebnisse, data_pop, data_smp, Pop)

# log ------------
Est_Tau_ebp_Matrix<- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_bc_Matrix <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_n_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
TrueMeans_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
ResEstPar         <- list()

log <- list()

load(file = "H:/Transformationen/Pop-Data/s1_normal_500.RData")

for(i in 1:Sim_Anzahl){
  print(paste("log", i))
  
  data_smp <- Pop[[i]]
  data_pop <- attr(Pop[[i]], "pop")
  
  Ergebnisse <- bhf_trafo(data_pop = data_pop, data_smp = data_smp,
                          domains = domains, formel = formel, sel = 5)
  
  parameter <- list()
  parameter[[1]] <- Ergebnisse$ResEstPar$optpar
  parameter[[2]] <- Ergebnisse$ResEstPar$betas
  names(parameter) <- c("optpar", "beta")
  
  print(Ergebnisse$mean_est_bc)
  
  ResEstPar[[i]]        <- parameter
  Est_Tau_bc_Matrix[,i] <- Ergebnisse$mean_est_bc
  Est_Tau_n_Matrix[,i]  <- Ergebnisse$mean_est_naive
  Est_Tau_ebp_Matrix[,i]<- Ergebnisse$mean_est_ebp
  TrueMeans_Matrix[,i]  <- Ergebnisse$mean_true
}

log[[1]] <- ResEstPar
log[[2]] <- Est_Tau_bc_Matrix
log[[3]] <- Est_Tau_n_Matrix
log[[4]] <- Est_Tau_ebp_Matrix
log[[5]] <- TrueMeans_Matrix
names(log) <- c("ResEstPar", "mean_bc", "mean_naive", "mean_ebp", "mean_true" )

rm(ResEstPar, Est_Tau_bc_Matrix, Est_Tau_n_Matrix, Est_Tau_ebp_Matrix, TrueMeans_Matrix, 
   parameter, Ergebnisse, data_pop, data_smp, Pop)

# log-shift ------------
Est_Tau_ebp_Matrix<- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_bc_Matrix <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_n_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
TrueMeans_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
ResEstPar         <- list()

log_shift <- list()

load(file = "H:/Transformationen/Pop-Data/s1_normal_500.RData")

for(i in 1:Sim_Anzahl){
  print(paste("log-shift", i))
  
  data_smp <- Pop[[i]]
  data_pop <- attr(Pop[[i]], "pop")
  
  Ergebnisse <- bhf_trafo(data_pop = data_pop, data_smp = data_smp,
                          domains = domains, formel = formel, sel = 1)

  print(Ergebnisse$mean_est_bc)
  
  parameter <- list()
  parameter[[1]] <- Ergebnisse$ResEstPar$optpar
  parameter[[2]] <- Ergebnisse$ResEstPar$betas
  names(parameter) <- c("optpar", "beta")
  
  ResEstPar[[i]] <- parameter
  Est_Tau_bc_Matrix[,i] <- Ergebnisse$mean_est_bc
  Est_Tau_n_Matrix[,i]  <- Ergebnisse$mean_est_naive
  Est_Tau_ebp_Matrix[,i]<- Ergebnisse$mean_est_ebp
  TrueMeans_Matrix[,i]  <- Ergebnisse$mean_true
}

log_shift[[1]] <- ResEstPar
log_shift[[2]] <- Est_Tau_bc_Matrix
log_shift[[3]] <- Est_Tau_n_Matrix
log_shift[[4]] <- Est_Tau_ebp_Matrix
log_shift[[5]] <- TrueMeans_Matrix
names(log_shift) <- c("ResEstPar", "mean_bc", "mean_naive", "mean_ebp", "mean_true" )

rm(ResEstPar, Est_Tau_bc_Matrix, Est_Tau_n_Matrix, Est_Tau_ebp_Matrix, TrueMeans_Matrix, 
   parameter, Ergebnisse, data_pop, data_smp, Pop)


# Box-Cox ------------
Est_Tau_ebp_Matrix<- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_bc_Matrix <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_n_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
TrueMeans_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
ResEstPar         <- list()

box_cox <- list()

load(file = "H:/Transformationen/Pop-Data/s1_normal_500.RData")

for(i in 1:Sim_Anzahl){
  print(paste("box-cox", i))
  
  data_smp <- Pop[[i]]
  data_pop <- attr(Pop[[i]], "pop")
  
  Ergebnisse <- bhf_trafo(data_pop = data_pop, data_smp = data_smp,
                          domains = domains, formel = formel, sel = 4)
  
  parameter <- list()
  parameter[[1]] <- Ergebnisse$ResEstPar$optpar
  parameter[[2]] <- Ergebnisse$ResEstPar$betas
  names(parameter) <- c("optpar", "beta")
  
  ResEstPar[[i]] <- parameter
  Est_Tau_bc_Matrix[,i] <- Ergebnisse$mean_est_bc
  Est_Tau_n_Matrix[,i]  <- Ergebnisse$mean_est_naive
  Est_Tau_ebp_Matrix[,i]<- Ergebnisse$mean_est_ebp
  TrueMeans_Matrix[,i]  <- Ergebnisse$mean_true
}

box_cox[[1]] <- ResEstPar
box_cox[[2]] <- Est_Tau_bc_Matrix
box_cox[[3]] <- Est_Tau_n_Matrix
box_cox[[4]] <- Est_Tau_ebp_Matrix
box_cox[[5]] <- TrueMeans_Matrix
names(box_cox) <- c("ResEstPar", "mean_bc", "mean_naive", "mean_ebp", "mean_true" )

rm(ResEstPar, Est_Tau_bc_Matrix, Est_Tau_n_Matrix, Est_Tau_ebp_Matrix, TrueMeans_Matrix, 
   parameter, Ergebnisse, data_pop, data_smp, Pop)


# Dual-Power ------------
Est_Tau_ebp_Matrix<- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_bc_Matrix <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_n_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
TrueMeans_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
ResEstPar         <- list()

dual_power <- list()

load(file = "H:/Transformationen/Pop-Data/s1_normal_500.RData")

for(i in 1:Sim_Anzahl){
  print(paste("dual-power", i))
  
  data_smp <- Pop[[i]]
  data_pop <- attr(Pop[[i]], "pop")
  
  Ergebnisse <- bhf_trafo(data_pop = data_pop, data_smp = data_smp,
                          domains = domains, formel = formel, sel = 9)
  
  parameter <- list()
  parameter[[1]] <- Ergebnisse$ResEstPar$optpar
  parameter[[2]] <- Ergebnisse$ResEstPar$betas
  names(parameter) <- c("optpar", "beta")
  
  ResEstPar[[i]] <- parameter
  Est_Tau_bc_Matrix[,i] <- Ergebnisse$mean_est_bc
  Est_Tau_n_Matrix[,i]  <- Ergebnisse$mean_est_naive
  Est_Tau_ebp_Matrix[,i]<- Ergebnisse$mean_est_ebp
  TrueMeans_Matrix[,i]  <- Ergebnisse$mean_true
}

dual_power[[1]] <- ResEstPar
dual_power[[2]] <- Est_Tau_bc_Matrix
dual_power[[3]] <- Est_Tau_n_Matrix
dual_power[[4]] <- Est_Tau_ebp_Matrix
dual_power[[5]] <- TrueMeans_Matrix
names(dual_power) <- c("ResEstPar", "mean_bc", "mean_naive", "mean_ebp", "mean_true" )

rm(ResEstPar, Est_Tau_bc_Matrix, Est_Tau_n_Matrix, Est_Tau_ebp_Matrix, TrueMeans_Matrix, 
   parameter, Ergebnisse, data_pop, data_smp, Pop)


# Ergebnisse speicher -----------
setwd("Sim Ergebnisse/normal/")

save(ohne_trafo, file = "ohne_trafo.Rdata")
save(log, file = "log.Rdata")
save(log_shift, file = "log_shift.Rdata")
save(box_cox, file = "box_cox.Rdata")
save(dual_power, file = "dual_power.Rdata")


# 2. Log-Scale Seznario--------
rm(normal, log, log_shift, box_cox, dual_power)
# Simulationen ------------------

Sim_Anzahl <- 500

domains  <- "idD"
formel <- y ~ x

# ohne trafo ------------
Est_Tau_ebp_Matrix<- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_bc_Matrix <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_n_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
TrueMeans_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
ResEstPar         <- list()

ohne_trafo <- list()

load(file = "H:/Transformationen/Pop-Data/s5_logscale_500.RData")

for(i in 1:Sim_Anzahl){
  print(paste("ohne Trafo", i))
  
  data_smp <- Pop[[i]]
  data_pop <- attr(Pop[[i]], "pop")
  
  Ergebnisse <- bhf_trafo(data_pop = data_pop, data_smp = data_smp,
                          domains = domains, formel = formel, sel = 2)
  
  parameter <- list()
  parameter[[1]] <- Ergebnisse$ResEstPar$optpar
  parameter[[2]] <- Ergebnisse$ResEstPar$betas
  names(parameter) <- c("optpar", "beta")
  
  ResEstPar[[i]]        <- parameter
  Est_Tau_bc_Matrix[,i] <- Ergebnisse$mean_est_bc
  Est_Tau_n_Matrix[,i]  <- Ergebnisse$mean_est_naive
  Est_Tau_ebp_Matrix[,i]<- Ergebnisse$mean_est_ebp
  TrueMeans_Matrix[,i]  <- Ergebnisse$mean_true
}

ohne_trafo[[1]] <- ResEstPar
ohne_trafo[[2]] <- Est_Tau_bc_Matrix
ohne_trafo[[3]] <- Est_Tau_n_Matrix
ohne_trafo[[4]] <- Est_Tau_ebp_Matrix
ohne_trafo[[5]] <- TrueMeans_Matrix
names(ohne_trafo) <- c("ResEstPar", "mean_bc", "mean_naive", "mean_ebp", "mean_true" )

rm(ResEstPar, Est_Tau_bc_Matrix, Est_Tau_n_Matrix, Est_Tau_ebp_Matrix, TrueMeans_Matrix, 
   parameter, Ergebnisse, data_pop, data_smp, Pop)

# log ------------
Est_Tau_ebp_Matrix<- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_bc_Matrix <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_n_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
TrueMeans_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
ResEstPar         <- list()

log <- list()

load(file = "H:/Transformationen/Pop-Data/s5_logscale_500.RData")

for(i in 1:Sim_Anzahl){
  print(paste("log", i))
  
  data_smp <- Pop[[i]]
  data_pop <- attr(Pop[[i]], "pop")
  
  Ergebnisse <- bhf_trafo(data_pop = data_pop, data_smp = data_smp,
                          domains = domains, formel = formel, sel = 5)
  
  parameter <- list()
  parameter[[1]] <- Ergebnisse$ResEstPar$optpar
  parameter[[2]] <- Ergebnisse$ResEstPar$betas
  names(parameter) <- c("optpar", "beta")
  
  ResEstPar[[i]]        <- parameter
  Est_Tau_bc_Matrix[,i] <- Ergebnisse$mean_est_bc
  Est_Tau_n_Matrix[,i]  <- Ergebnisse$mean_est_naive
  Est_Tau_ebp_Matrix[,i]<- Ergebnisse$mean_est_ebp
  TrueMeans_Matrix[,i]  <- Ergebnisse$mean_true
}

log[[1]] <- ResEstPar
log[[2]] <- Est_Tau_bc_Matrix
log[[3]] <- Est_Tau_n_Matrix
log[[4]] <- Est_Tau_ebp_Matrix
log[[5]] <- TrueMeans_Matrix
names(log) <- c("ResEstPar", "mean_bc", "mean_naive", "mean_ebp", "mean_true" )

rm(ResEstPar, Est_Tau_bc_Matrix, Est_Tau_n_Matrix, Est_Tau_ebp_Matrix, TrueMeans_Matrix, 
   parameter, Ergebnisse, data_pop, data_smp, Pop)

# log-shift ------------
Est_Tau_ebp_Matrix<- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_bc_Matrix <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_n_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
TrueMeans_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
ResEstPar         <- list()

log_shift <- list()

load(file = "H:/Transformationen/Pop-Data/s5_logscale_500.RData")

for(i in 1:Sim_Anzahl){
  print(paste("log-shift", i))
  
  data_smp <- Pop[[i]]
  data_pop <- attr(Pop[[i]], "pop")
  
  Ergebnisse <- bhf_trafo(data_pop = data_pop, data_smp = data_smp,
                          domains = domains, formel = formel, sel = 1)
  
  parameter <- list()
  parameter[[1]] <- Ergebnisse$ResEstPar$optpar
  parameter[[2]] <- Ergebnisse$ResEstPar$betas
  names(parameter) <- c("optpar", "beta")
  
  ResEstPar[[i]] <- parameter
  Est_Tau_bc_Matrix[,i] <- Ergebnisse$mean_est_bc
  Est_Tau_n_Matrix[,i]  <- Ergebnisse$mean_est_naive
  Est_Tau_ebp_Matrix[,i]<- Ergebnisse$mean_est_ebp
  TrueMeans_Matrix[,i]  <- Ergebnisse$mean_true
}

log_shift[[1]] <- ResEstPar
log_shift[[2]] <- Est_Tau_bc_Matrix
log_shift[[3]] <- Est_Tau_n_Matrix
log_shift[[4]] <- Est_Tau_ebp_Matrix
log_shift[[5]] <- TrueMeans_Matrix
names(log_shift) <- c("ResEstPar", "mean_bc", "mean_naive", "mean_ebp", "mean_true" )

rm(ResEstPar, Est_Tau_bc_Matrix, Est_Tau_n_Matrix, Est_Tau_ebp_Matrix, TrueMeans_Matrix, 
   parameter, Ergebnisse, data_pop, data_smp, Pop)


# Box-Cox ------------
Est_Tau_ebp_Matrix<- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_bc_Matrix <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_n_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
TrueMeans_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
ResEstPar         <- list()

box_cox <- list()

load(file = "H:/Transformationen/Pop-Data/s5_logscale_500.RData")

for(i in 1:Sim_Anzahl){
  print(paste("box-cox", i))
  
  data_smp <- Pop[[i]]
  data_pop <- attr(Pop[[i]], "pop")
  
  Ergebnisse <- bhf_trafo(data_pop = data_pop, data_smp = data_smp,
                          domains = domains, formel = formel, sel = 4)
  
  parameter <- list()
  parameter[[1]] <- Ergebnisse$ResEstPar$optpar
  parameter[[2]] <- Ergebnisse$ResEstPar$betas
  names(parameter) <- c("optpar", "beta")
  
  print(which(is.na(Ergebnisse$mean_est_naive)))
  
  ResEstPar[[i]] <- parameter
  Est_Tau_bc_Matrix[,i] <- Ergebnisse$mean_est_bc
  Est_Tau_n_Matrix[,i]  <- Ergebnisse$mean_est_naive
  Est_Tau_ebp_Matrix[,i]<- Ergebnisse$mean_est_ebp
  TrueMeans_Matrix[,i]  <- Ergebnisse$mean_true
}

box_cox[[1]] <- ResEstPar
box_cox[[2]] <- Est_Tau_bc_Matrix
box_cox[[3]] <- Est_Tau_n_Matrix
box_cox[[4]] <- Est_Tau_ebp_Matrix
box_cox[[5]] <- TrueMeans_Matrix
names(box_cox) <- c("ResEstPar", "mean_bc", "mean_naive", "mean_ebp", "mean_true" )

rm(ResEstPar, Est_Tau_bc_Matrix, Est_Tau_n_Matrix, Est_Tau_ebp_Matrix, TrueMeans_Matrix, 
   parameter, Ergebnisse, data_pop, data_smp, Pop)


# Dual-Power ------------
Est_Tau_ebp_Matrix<- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_bc_Matrix <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_n_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
TrueMeans_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
ResEstPar         <- list()

dual_power <- list()

load(file = "H:/Transformationen/Pop-Data/s5_logscale_500.RData")

for(i in 1:Sim_Anzahl){
  print(paste("dual-power", i))
  
  data_smp <- Pop[[i]]
  data_pop <- attr(Pop[[i]], "pop")
  
  Ergebnisse <- bhf_trafo(data_pop = data_pop, data_smp = data_smp,
                          domains = domains, formel = formel, sel = 9)
  
  parameter <- list()
  parameter[[1]] <- Ergebnisse$ResEstPar$optpar
  parameter[[2]] <- Ergebnisse$ResEstPar$betas
  names(parameter) <- c("optpar", "beta")
  
  ResEstPar[[i]] <- parameter
  Est_Tau_bc_Matrix[,i] <- Ergebnisse$mean_est_bc
  Est_Tau_n_Matrix[,i]  <- Ergebnisse$mean_est_naive
  Est_Tau_ebp_Matrix[,i]<- Ergebnisse$mean_est_ebp
  TrueMeans_Matrix[,i]  <- Ergebnisse$mean_true
}

dual_power[[1]] <- ResEstPar
dual_power[[2]] <- Est_Tau_bc_Matrix
dual_power[[3]] <- Est_Tau_n_Matrix
dual_power[[4]] <- Est_Tau_ebp_Matrix
dual_power[[5]] <- TrueMeans_Matrix
names(dual_power) <- c("ResEstPar", "mean_bc", "mean_naive", "mean_ebp", "mean_true" )

rm(ResEstPar, Est_Tau_bc_Matrix, Est_Tau_n_Matrix, Est_Tau_ebp_Matrix, TrueMeans_Matrix, 
   parameter, Ergebnisse, data_pop, data_smp, Pop)


# Ergebnisse speicher -----------
setwd("../logscale/")

save(ohne_trafo, file = "ohne_trafo.Rdata")
save(log, file = "log.Rdata")
save(log_shift, file = "log_shift.Rdata")
save(box_cox, file = "box_cox.Rdata")
save(dual_power, file = "dual_power.Rdata")

# 3. GB2 Seznario--------
rm(ohne_trafo, log, log_shift, box_cox, dual_power)
# Simulationen ------------------

Sim_Anzahl <- 500

domains  <- "idD"
formel <- y ~ x

# ohne trafo ------------
Est_Tau_ebp_Matrix<- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_bc_Matrix <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_n_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
TrueMeans_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
ResEstPar         <- list()

ohne_trafo <- list()

load(file = "H:/Transformationen/Pop-Data/s9_GB2_500.Rdata")

for(i in 1:Sim_Anzahl){
  print(paste("ohne Trafo", i))
  
  data_smp <- Pop[[i]]
  data_pop <- attr(Pop[[i]], "pop")
  
  Ergebnisse <- bhf_trafo(data_pop = data_pop, data_smp = data_smp,
                          domains = domains, formel = formel, sel = 2)
  
  parameter <- list()
  parameter[[1]] <- Ergebnisse$ResEstPar$optpar
  parameter[[2]] <- Ergebnisse$ResEstPar$betas
  names(parameter) <- c("optpar", "beta")
  
  ResEstPar[[i]]        <- parameter
  Est_Tau_bc_Matrix[,i] <- Ergebnisse$mean_est_bc
  Est_Tau_n_Matrix[,i]  <- Ergebnisse$mean_est_naive
  Est_Tau_ebp_Matrix[,i]<- Ergebnisse$mean_est_ebp
  TrueMeans_Matrix[,i]  <- Ergebnisse$mean_true
}

ohne_trafo[[1]] <- ResEstPar
ohne_trafo[[2]] <- Est_Tau_bc_Matrix
ohne_trafo[[3]] <- Est_Tau_n_Matrix
ohne_trafo[[4]] <- Est_Tau_ebp_Matrix
ohne_trafo[[5]] <- TrueMeans_Matrix
names(ohne_trafo) <- c("ResEstPar", "mean_bc", "mean_naive", "mean_ebp", "mean_true" )

rm(ResEstPar, Est_Tau_bc_Matrix, Est_Tau_n_Matrix, Est_Tau_ebp_Matrix, TrueMeans_Matrix, 
   parameter, Ergebnisse, data_pop, data_smp, Pop)

# log ------------
Est_Tau_ebp_Matrix<- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_bc_Matrix <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_n_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
TrueMeans_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
ResEstPar         <- list()

log <- list()

load(file = "H:/Transformationen/Pop-Data/s9_GB2_500.Rdata")

for(i in 1:Sim_Anzahl){
  print(paste("log", i))
  
  data_smp <- Pop[[i]]
  data_pop <- attr(Pop[[i]], "pop")
  
  Ergebnisse <- bhf_trafo(data_pop = data_pop, data_smp = data_smp,
                          domains = domains, formel = formel, sel = 5)
  
  parameter <- list()
  parameter[[1]] <- Ergebnisse$ResEstPar$optpar
  parameter[[2]] <- Ergebnisse$ResEstPar$betas
  names(parameter) <- c("optpar", "beta")
  
  ResEstPar[[i]]        <- parameter
  Est_Tau_bc_Matrix[,i] <- Ergebnisse$mean_est_bc
  Est_Tau_n_Matrix[,i]  <- Ergebnisse$mean_est_naive
  Est_Tau_ebp_Matrix[,i]<- Ergebnisse$mean_est_ebp
  TrueMeans_Matrix[,i]  <- Ergebnisse$mean_true
}

log[[1]] <- ResEstPar
log[[2]] <- Est_Tau_bc_Matrix
log[[3]] <- Est_Tau_n_Matrix
log[[4]] <- Est_Tau_ebp_Matrix
log[[5]] <- TrueMeans_Matrix
names(log) <- c("ResEstPar", "mean_bc", "mean_naive", "mean_ebp", "mean_true" )

rm(ResEstPar, Est_Tau_bc_Matrix, Est_Tau_n_Matrix, Est_Tau_ebp_Matrix, TrueMeans_Matrix, 
   parameter, Ergebnisse, data_pop, data_smp, Pop)

# log-shift ------------
Est_Tau_ebp_Matrix<- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_bc_Matrix <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_n_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
TrueMeans_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
ResEstPar         <- list()

log_shift <- list()

load(file = "H:/Transformationen/Pop-Data/s9_GB2_500.Rdata")

for(i in 1:Sim_Anzahl){
  print(paste("log-shift", i))
  
  data_smp <- Pop[[i]]
  data_pop <- attr(Pop[[i]], "pop")
  
  Ergebnisse <- bhf_trafo(data_pop = data_pop, data_smp = data_smp,
                          domains = domains, formel = formel, sel = 1)
  
  parameter <- list()
  parameter[[1]] <- Ergebnisse$ResEstPar$optpar
  parameter[[2]] <- Ergebnisse$ResEstPar$betas
  names(parameter) <- c("optpar", "beta")
  
  ResEstPar[[i]] <- parameter
  Est_Tau_bc_Matrix[,i] <- Ergebnisse$mean_est_bc
  Est_Tau_n_Matrix[,i]  <- Ergebnisse$mean_est_naive
  Est_Tau_ebp_Matrix[,i]<- Ergebnisse$mean_est_ebp
  TrueMeans_Matrix[,i]  <- Ergebnisse$mean_true
}

log_shift[[1]] <- ResEstPar
log_shift[[2]] <- Est_Tau_bc_Matrix
log_shift[[3]] <- Est_Tau_n_Matrix
log_shift[[4]] <- Est_Tau_ebp_Matrix
log_shift[[5]] <- TrueMeans_Matrix
names(log_shift) <- c("ResEstPar", "mean_bc", "mean_naive", "mean_ebp", "mean_true" )

rm(ResEstPar, Est_Tau_bc_Matrix, Est_Tau_n_Matrix, Est_Tau_ebp_Matrix, TrueMeans_Matrix, 
   parameter, Ergebnisse, data_pop, data_smp, Pop)


# Box-Cox ------------
Est_Tau_ebp_Matrix<- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_bc_Matrix <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_n_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
TrueMeans_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
ResEstPar         <- list()

box_cox <- list()

load(file = "H:/Transformationen/Pop-Data/s9_GB2_500.Rdata")

for(i in 1:Sim_Anzahl){
  print(paste("box-cox", i))
  
  data_smp <- Pop[[i]]
  data_pop <- attr(Pop[[i]], "pop")
  
  Ergebnisse <- bhf_trafo(data_pop = data_pop, data_smp = data_smp,
                          domains = domains, formel = formel, sel = 4)
  
  parameter <- list()
  parameter[[1]] <- Ergebnisse$ResEstPar$optpar
  parameter[[2]] <- Ergebnisse$ResEstPar$betas
  names(parameter) <- c("optpar", "beta")
  
  ResEstPar[[i]] <- parameter
  Est_Tau_bc_Matrix[,i] <- Ergebnisse$mean_est_bc
  Est_Tau_n_Matrix[,i]  <- Ergebnisse$mean_est_naive
  Est_Tau_ebp_Matrix[,i]<- Ergebnisse$mean_est_ebp
  TrueMeans_Matrix[,i]  <- Ergebnisse$mean_true
}

box_cox[[1]] <- ResEstPar
box_cox[[2]] <- Est_Tau_bc_Matrix
box_cox[[3]] <- Est_Tau_n_Matrix
box_cox[[4]] <- Est_Tau_ebp_Matrix
box_cox[[5]] <- TrueMeans_Matrix
names(box_cox) <- c("ResEstPar", "mean_bc", "mean_naive", "mean_ebp", "mean_true" )

rm(ResEstPar, Est_Tau_bc_Matrix, Est_Tau_n_Matrix, Est_Tau_ebp_Matrix, TrueMeans_Matrix, 
   parameter, Ergebnisse, data_pop, data_smp, Pop)


# Dual-Power ------------
Est_Tau_ebp_Matrix<- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_bc_Matrix <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_n_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
TrueMeans_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
ResEstPar         <- list()

dual_power <- list()

load(file = "H:/Transformationen/Pop-Data/s9_GB2_500.Rdata")

for(i in 1:Sim_Anzahl){
  print(paste("dual-power", i))
  
  data_smp <- Pop[[i]]
  data_pop <- attr(Pop[[i]], "pop")
  
  Ergebnisse <- bhf_trafo(data_pop = data_pop, data_smp = data_smp,
                          domains = domains, formel = formel, sel = 9)
  
  parameter <- list()
  parameter[[1]] <- Ergebnisse$ResEstPar$optpar
  parameter[[2]] <- Ergebnisse$ResEstPar$betas
  names(parameter) <- c("optpar", "beta")
  
  ResEstPar[[i]] <- parameter
  Est_Tau_bc_Matrix[,i] <- Ergebnisse$mean_est_bc
  Est_Tau_n_Matrix[,i]  <- Ergebnisse$mean_est_naive
  Est_Tau_ebp_Matrix[,i]<- Ergebnisse$mean_est_ebp
  TrueMeans_Matrix[,i]  <- Ergebnisse$mean_true
}

dual_power[[1]] <- ResEstPar
dual_power[[2]] <- Est_Tau_bc_Matrix
dual_power[[3]] <- Est_Tau_n_Matrix
dual_power[[4]] <- Est_Tau_ebp_Matrix
dual_power[[5]] <- TrueMeans_Matrix
names(dual_power) <- c("ResEstPar", "mean_bc", "mean_naive", "mean_ebp", "mean_true" )

rm(ResEstPar, Est_Tau_bc_Matrix, Est_Tau_n_Matrix, Est_Tau_ebp_Matrix, TrueMeans_Matrix, 
   parameter, Ergebnisse, data_pop, data_smp, Pop)


# Ergebnisse speicher -----------
setwd("../GB2/")

save(ohne_trafo, file = "ohne_trafo.Rdata")
save(log, file = "log.Rdata")
save(log_shift, file = "log_shift.Rdata")
save(box_cox, file = "box_cox.Rdata")
save(dual_power, file = "dual_power.Rdata")

# 4. Paretro Seznario--------
rm(ohne_trafo, log, log_shift, box_cox, dual_power)
# Simulationen ------------------

Sim_Anzahl <- 500

domains  <- "idD"
formel <- y ~ x

# ohne trafo ------------
Est_Tau_ebp_Matrix<- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_bc_Matrix <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_n_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
TrueMeans_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
ResEstPar         <- list()

ohne_trafo <- list()

load(file = "H:/Transformationen/Pop-Data/s7_pareto_500.RData")

for(i in 1:Sim_Anzahl){
  print(paste("ohne Trafo", i))
  
  data_smp <- Pop[[i]]
  data_pop <- attr(Pop[[i]], "pop")
  
  Ergebnisse <- bhf_trafo(data_pop = data_pop, data_smp = data_smp,
                          domains = domains, formel = formel, sel = 2)
  
  parameter <- list()
  parameter[[1]] <- Ergebnisse$ResEstPar$optpar
  parameter[[2]] <- Ergebnisse$ResEstPar$betas
  names(parameter) <- c("optpar", "beta")
  
  ResEstPar[[i]]        <- parameter
  Est_Tau_bc_Matrix[,i] <- Ergebnisse$mean_est_bc
  Est_Tau_n_Matrix[,i]  <- Ergebnisse$mean_est_naive
  Est_Tau_ebp_Matrix[,i]<- Ergebnisse$mean_est_ebp
  TrueMeans_Matrix[,i]  <- Ergebnisse$mean_true
}

ohne_trafo[[1]] <- ResEstPar
ohne_trafo[[2]] <- Est_Tau_bc_Matrix
ohne_trafo[[3]] <- Est_Tau_n_Matrix
ohne_trafo[[4]] <- Est_Tau_ebp_Matrix
ohne_trafo[[5]] <- TrueMeans_Matrix
names(ohne_trafo) <- c("ResEstPar", "mean_bc", "mean_naive", "mean_ebp", "mean_true" )

rm(ResEstPar, Est_Tau_bc_Matrix, Est_Tau_n_Matrix, Est_Tau_ebp_Matrix, TrueMeans_Matrix, 
   parameter, Ergebnisse, data_pop, data_smp, Pop)

# log ------------
Est_Tau_ebp_Matrix<- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_bc_Matrix <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_n_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
TrueMeans_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
ResEstPar         <- list()

log <- list()

load(file = "H:/Transformationen/Pop-Data/s7_pareto_500.RData")

for(i in 1:Sim_Anzahl){
  print(paste("log", i))
  
  data_smp <- Pop[[i]]
  data_pop <- attr(Pop[[i]], "pop")
  
  Ergebnisse <- bhf_trafo(data_pop = data_pop, data_smp = data_smp,
                          domains = domains, formel = formel, sel = 5)
  
  parameter <- list()
  parameter[[1]] <- Ergebnisse$ResEstPar$optpar
  parameter[[2]] <- Ergebnisse$ResEstPar$betas
  names(parameter) <- c("optpar", "beta")
  
  ResEstPar[[i]]        <- parameter
  Est_Tau_bc_Matrix[,i] <- Ergebnisse$mean_est_bc
  Est_Tau_n_Matrix[,i]  <- Ergebnisse$mean_est_naive
  Est_Tau_ebp_Matrix[,i]<- Ergebnisse$mean_est_ebp
  TrueMeans_Matrix[,i]  <- Ergebnisse$mean_true
}

log[[1]] <- ResEstPar
log[[2]] <- Est_Tau_bc_Matrix
log[[3]] <- Est_Tau_n_Matrix
log[[4]] <- Est_Tau_ebp_Matrix
log[[5]] <- TrueMeans_Matrix
names(log) <- c("ResEstPar", "mean_bc", "mean_naive", "mean_ebp", "mean_true" )

rm(ResEstPar, Est_Tau_bc_Matrix, Est_Tau_n_Matrix, Est_Tau_ebp_Matrix, TrueMeans_Matrix, 
   parameter, Ergebnisse, data_pop, data_smp, Pop)

# log-shift ------------
Est_Tau_ebp_Matrix<- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_bc_Matrix <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_n_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
TrueMeans_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
ResEstPar         <- list()

log_shift <- list()

load(file = "H:/Transformationen/Pop-Data/s7_pareto_500.RData")

for(i in 1:Sim_Anzahl){
  print(paste("log-shift", i))
  
  data_smp <- Pop[[i]]
  data_pop <- attr(Pop[[i]], "pop")
  
  Ergebnisse <- bhf_trafo(data_pop = data_pop, data_smp = data_smp,
                          domains = domains, formel = formel, sel = 1)
  
  parameter <- list()
  parameter[[1]] <- Ergebnisse$ResEstPar$optpar
  parameter[[2]] <- Ergebnisse$ResEstPar$betas
  names(parameter) <- c("optpar", "beta")
  
  ResEstPar[[i]] <- parameter
  Est_Tau_bc_Matrix[,i] <- Ergebnisse$mean_est_bc
  Est_Tau_n_Matrix[,i]  <- Ergebnisse$mean_est_naive
  Est_Tau_ebp_Matrix[,i]<- Ergebnisse$mean_est_ebp
  TrueMeans_Matrix[,i]  <- Ergebnisse$mean_true
}

log_shift[[1]] <- ResEstPar
log_shift[[2]] <- Est_Tau_bc_Matrix
log_shift[[3]] <- Est_Tau_n_Matrix
log_shift[[4]] <- Est_Tau_ebp_Matrix
log_shift[[5]] <- TrueMeans_Matrix
names(log_shift) <- c("ResEstPar", "mean_bc", "mean_naive", "mean_ebp", "mean_true" )

rm(ResEstPar, Est_Tau_bc_Matrix, Est_Tau_n_Matrix, Est_Tau_ebp_Matrix, TrueMeans_Matrix, 
   parameter, Ergebnisse, data_pop, data_smp, Pop)


# Box-Cox ------------
Est_Tau_ebp_Matrix<- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_bc_Matrix <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_n_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
TrueMeans_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
ResEstPar         <- list()

box_cox <- list()

load(file = "H:/Transformationen/Pop-Data/s7_pareto_500.RData")

for(i in 1:Sim_Anzahl){
  print(paste("box-cox", i))
  
  data_smp <- Pop[[i]]
  data_pop <- attr(Pop[[i]], "pop")
  
  Ergebnisse <- bhf_trafo(data_pop = data_pop, data_smp = data_smp,
                          domains = domains, formel = formel, sel = 4)
  
  parameter <- list()
  parameter[[1]] <- Ergebnisse$ResEstPar$optpar
  parameter[[2]] <- Ergebnisse$ResEstPar$betas
  names(parameter) <- c("optpar", "beta")
  
  ResEstPar[[i]] <- parameter
  Est_Tau_bc_Matrix[,i] <- Ergebnisse$mean_est_bc
  Est_Tau_n_Matrix[,i]  <- Ergebnisse$mean_est_naive
  Est_Tau_ebp_Matrix[,i]<- Ergebnisse$mean_est_ebp
  TrueMeans_Matrix[,i]  <- Ergebnisse$mean_true
}

box_cox[[1]] <- ResEstPar
box_cox[[2]] <- Est_Tau_bc_Matrix
box_cox[[3]] <- Est_Tau_n_Matrix
box_cox[[4]] <- Est_Tau_ebp_Matrix
box_cox[[5]] <- TrueMeans_Matrix
names(box_cox) <- c("ResEstPar", "mean_bc", "mean_naive", "mean_ebp", "mean_true" )

rm(ResEstPar, Est_Tau_bc_Matrix, Est_Tau_n_Matrix, Est_Tau_ebp_Matrix, TrueMeans_Matrix, 
   parameter, Ergebnisse, data_pop, data_smp, Pop)


# Dual-Power ------------
Est_Tau_ebp_Matrix<- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_bc_Matrix <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
Est_Tau_n_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
TrueMeans_Matrix  <- matrix(NA, ncol = Sim_Anzahl, nrow = 50)
ResEstPar         <- list()

dual_power <- list()

load(file = "H:/Transformationen/Pop-Data/s7_pareto_500.RData")

for(i in 1:Sim_Anzahl){
  print(paste("dual-power", i))
  
  data_smp <- Pop[[i]]
  data_pop <- attr(Pop[[i]], "pop")
  
  Ergebnisse <- bhf_trafo(data_pop = data_pop, data_smp = data_smp,
                          domains = domains, formel = formel, sel = 9)
  
  parameter <- list()
  parameter[[1]] <- Ergebnisse$ResEstPar$optpar
  parameter[[2]] <- Ergebnisse$ResEstPar$betas
  names(parameter) <- c("optpar", "beta")
  
  ResEstPar[[i]] <- parameter
  Est_Tau_bc_Matrix[,i] <- Ergebnisse$mean_est_bc
  Est_Tau_n_Matrix[,i]  <- Ergebnisse$mean_est_naive
  Est_Tau_ebp_Matrix[,i]<- Ergebnisse$mean_est_ebp
  TrueMeans_Matrix[,i]  <- Ergebnisse$mean_true
}

dual_power[[1]] <- ResEstPar
dual_power[[2]] <- Est_Tau_bc_Matrix
dual_power[[3]] <- Est_Tau_n_Matrix
dual_power[[4]] <- Est_Tau_ebp_Matrix
dual_power[[5]] <- TrueMeans_Matrix
names(dual_power) <- c("ResEstPar", "mean_bc", "mean_naive", "mean_ebp", "mean_true" )

rm(ResEstPar, Est_Tau_bc_Matrix, Est_Tau_n_Matrix, Est_Tau_ebp_Matrix, TrueMeans_Matrix, 
   parameter, Ergebnisse, data_pop, data_smp, Pop)


# Ergebnisse speicher -----------
setwd("../Paretro/")

save(ohne_trafo, file = "ohne_trafo.Rdata")
save(log, file = "log.Rdata")
save(log_shift, file = "log_shift.Rdata")
save(box_cox, file = "box_cox.Rdata")
save(dual_power, file = "dual_power.Rdata")
