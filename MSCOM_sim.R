rm(list=ls())

library(dplyr)
library(reshape2)
library(ggplot2)

main_dir <- "C:\\merrill\\MS-COM\\mscom_sim"

R_dir <- file.path(main_dir, "R")

funs <- list.files(R_dir)
ignore <- sapply(1:length(funs), function(x) source(file.path(R_dir, funs[x])))

## tuna = yellowfin Thunnus albacares
## billfish = striped marlin Kajikia audax
## shark = mako Isurus oxyrinchus
species <- c("tuna", "billfish", "shark")
oneway_mat <- data.frame("SpeciesName"=species,
						"Fdynamics"="One-way",
						"InitialDepl"=1,
						"PercentFcrash"=c(0.4, 0.6, 0.8),
						"SigmaR"=0, "rho"=0,
						"SigmaF"=0,
						"r" = c(0.6, 0.67, 0.11),
						"K"=c(304, 227, 738))

sim_oneway <- sim_pops(input_mat=oneway_mat,
				nyears=50,
				seed=123, 
				model="biomass-dynamic")

p_oneway <- ggplot(sim_oneway %>% dplyr::filter(variable %in% c("RelativeCatch","ExploitRate","RelativeEffort","Depletion"))) +
	geom_line(aes(x=Year, y=value, colour=Species), lwd=2) +
	facet_wrap(~variable, scales='free_y') +
	theme_lsd() +
	coord_cartesian(ylim=c(0,1.01))
	
## get catch data
catch_oneway <- sim_oneway %>% filter(variable=="Catch")
input_oneway <- sapply(1:length(species), function(x){
	sub <- catch_oneway %>% filter(Species==species[x])
	df <- sub$value
	return(df)
})
colnames(input_oneway) <- species




#set of species for which catches are pooled
# one species goes up, another species goes down
# in general --> more stable than single species -- not as much contrast




constant_mat <- data.frame("SpeciesName"=c("tuna","billfish","shark"),
						"Fdynamics"="Constant",
						"InitialDepl"=1,
						"PercentFcrash"=c(0.4, 0.6, 0.8),
						"SigmaR"=0, "rho"=0,
						"SigmaF"=0,
						"r" = c(0.6, 0.67, 0.11),
						"K"=c(304, 227, 738))

sim_constant <- sim_pops(input_mat=constant_mat,
				nyears=50,
				seed=123, 
				model="biomass-dynamic")


p_const <- ggplot(sim_constant %>% dplyr::filter(variable %in% c("RelativeCatch","ExploitRate","RelativeEffort","Depletion"))) +
	geom_line(aes(x=Year, y=value, colour=Species), lwd=2) +
	facet_wrap(~variable, scales='free_y') +
	theme_lsd() + 
	coord_cartesian(ylim=c(0,1.01))

