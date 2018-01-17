rm(list=ls())

main_dir <- "C:\\merrill\\MS-COM\\mscom_sim"
R_dir <- file.path(main_dir, "R")

funs <- list.files(R_dir)
ignore <- sapply(1:length(funs), function(x) source(file.path(R_dir, funs[x])))

input_mat <- data.frame("SpeciesName"=c("tuna","billfish","shark"),
						"Fdynamics"="One-way",
						"InitialDepl"=c(0.9, 0.8, 0.5),
						"SigmaR"=0.737, "rho"=0.436,
						"SigmaF"=0.2,
						"linf"=c(50, 80, 100), "vbk"=c(0.2, 0.2, 0.1), "t0"=-0.01,
						"lwa"=0.005, "lwb"=3,
						"M"=c(0.2, 0.18, 0.1),
						"M50"=c(3, 4, 8),
						"S50"=c(3, 4, 3), "S95"=c(4, 5, 4),
						"h"=1,
						"R0"=1000)

sim <- sim_pops(input_mat=input_mat,
				nyears=20,
				nburn=100,
				seed=143)

par(mfrow=c(4,3))
for(x in 1:nrow(input_mat)){
	plot(sim[[x]]$F_t, type="l", lwd=2, ylim=c(0, max(sim[[x]]$F_t)*1.1))
}
for(x in 1:nrow(input_mat)){
	plot(sim[[x]]$R_t, type="l", lwd=2, ylim=c(0, max(sim[[x]]$R_t)*1.1))
}
for(x in 1:nrow(input_mat)){
	plot(sim[[x]]$D_t, type="l", lwd=2, ylim=c(0, max(sim[[x]]$D_t)*1.1))
}
for(x in 1:nrow(input_mat)){
	plot(sim[[x]]$N_t, type="l", lwd=2, ylim=c(0, max(sim[[x]]$N_t)*1.1))
}