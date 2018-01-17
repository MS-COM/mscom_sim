#' Operating model
#'
#' \code{sim_pops} Simulate population dynamics for multiple species in the same fishery
#'
#' @author M.B. Rudd, K. Kleisner, et al.
#' @param species vector of names of life history types -- programmed include tuna, billfish, shark
#' @importFrom FishLife
#' @return named list of attributes of true population/data
#' @export

sim_pops <- function(input_mat,
					nyears,
					nburn, 
					seed){

#######################
## life history types
#######################
	## **** calc r and K from life history types required for age-structure

	species <- unique(input_mat$SpeciesName)
	nspecies <- length(input_mat$SpeciesName)

	lh <- lapply(1:nspecies, function(x){
		sp <- input_mat %>% dplyr::filter(species==species[x])
		lh <- with(sp, lh_list(linf=linf, vbk=vbk, t0=t0,
						 lwa=lwa, lwb=lwb,
						 M=M, M50=M50,
						 S50=S50, S95=S95,
						 h=h, R0=R0))
		return(lh)
	})
	names(lh) <- species
 	

#######################
## Random variables
#######################
	tyears <- nburn + nyears

	set.seed(seed)

	## Recruitment deviations
	RecDev <- lapply(1:nspecies, function(x){
		sp <- input_mat %>% dplyr::filter(species==species[x])
		raw <- with(sp, rnorm(tyears, mean= -(SigmaR ^ 2)/2, sd=SigmaR))

		out <- rep(NA, tyears)
		out[1] <- raw[1]
		for(t in 2:length(raw)){
			out[t] <- with(sp, out[t-1] * rho + sqrt(1 - rho ^ 2) * raw[t])
		}
		return(out)
	})

	FishDev <- lapply(1:nspecies, function(x){
		sp <- input_mat %>% dplyr::filter(species==species[x])
		out <- with(sp, rnorm(tyears, mean= -(SigmaF ^ 2)/2, sd=SigmaF))
		return(out)
	})

#######################################
## initiate fishing mortality dynamics
#######################################

	## initial fishing mortality rate corresponding to init_depl
	## where init_depl == SPR not depletion
	Finit <- lapply(1:nspecies, function(x){
		choose_lh <- lh[[which(species == species[x])]]
		sp <- input_mat %>% dplyr::filter(species == species[x])

		out_Finit <- with(choose_lh, tryCatch(uniroot(
					            calc_ref,
					            lower = 0,
					            upper = 200,
					            ages = ages,
					            Mat_a = Mat_a,
					            W_a = W_a,
					            M = M,
					            S_a = S_a,
					            ref = sp$InitialDepl
					          )$root,
					          error = function(e)
					            NA))
	})

	## fishing mortality rate where SPR = 0.05
	Fcrash <- lapply(1:nspecies, function(x){
		choose_lh <- lh[[which(species == species[x])]]
		out_Fcrash <- with(choose_lh, tryCatch(uniroot(
					            calc_ref,
					            lower = 0,
					            upper = 200,
					            ages =ages,
					            Mat_a =Mat_a,
					            W_a =W_a,
					            M =M,
					            S_a =S_a,
					            ref = 0.05
					          )$root,
					          error = function(e)
					            NA))
	})

	## fishing mortality pattern over time
	F_t <- lapply(1:nspecies, function(x){
		sp <- input_mat %>% dplyr::filter(species == species[x])

		## initialize
		out <- rep(NA, tyears)

		if(sp$Fdynamics=="Constant") out <- rep(Finit[[x]], tyears) * exp(FishDev[[x]])
		if(sp$Fdynamics=="One-way"){
			for(t in 1:tyears){
				if(t <= (nburn+1)) out[t] <- Finit[[x]]
				if(t > (nburn+1)){
					out[t] <- ifelse(out[t-1] <= Fcrash[[x]], out[t-1]*1.05, Fcrash[[x]]) * exp(FishDev[[x]][t])
				}
			}
		}
		if(sp$Fdynamics=="EffortDyn"){
			out[1] <- Finit[[x]]
		}

		return(out)
	})

#######################################
## initiate age structure - year 1
#######################################

	## recruitment
	R_t <- lapply(1:nspecies, function(x){
		choose_lh <- lh[[which(species == species[x])]]
		out <- rep(NA, tyears)
		out[1] <- choose_lh$R0 * exp(RecDev[[x]][1])
		return(out)
	})

	## values at age over time
	N_at <- N_at0 <- Cn_at <- Cw_at <- lapply(1:nspecies, function(x){
		choose_lh <- lh[[which(species == species[x])]]
		return(matrix(NA, nrow=length(choose_lh$ages), ncol=tyears))
	})

	## values over time
	SB_t <- VB_t <- TB_t <- lapply(1:nspecies, function(x){
		return(rep(NA, tyears))
	})

	## loop over species
	for(x in 1:nspecies){
		choose_lh <- lh[[which(species == species[x])]]

		## abundance at age over time
		for(a in 1:length(choose_lh$ages)){
			if(a==1){
				N_at[[x]][a,1] <- R_t[[x]][1]
				N_at0[[x]][a,1] <- R_t[[x]][1]
			}
			if(a > 1 & a < length(choose_lh$ages)){
				N_at[[x]][a,1] <- N_at[[x]][a-1,1] * exp(-choose_lh$M - F_t[[x]][1]*choose_lh$S_a[a-1])
				N_at0[[x]][a,1] <- N_at0[[x]][a-1,1] * exp(-choose_lh$M)
			}
			if(a == length(choose_lh$ages)){
				N_at[[x]][a,1] <- N_at[[x]][a-1,1] * exp(-choose_lh$M - F_t[[x]][1]*choose_lh$S_a[a]) / (1 - exp(-choose_lh$M - F_t[[x]][1]*choose_lh$S_a[a]))
				N_at0[[x]][a,1] <- N_at0[[x]][a-1,1] * exp(-choose_lh$M) / (1 - exp(-choose_lh$M))
			}
		}

		## biomass over time
		VB_t[[x]][1] <- sum(N_at[[x]][,1] * choose_lh$W_a * choose_lh$S_a)
		TB_t[[x]][1] <- sum(N_at[[x]][,1] * choose_lh$W_a)
		SB_t[[x]][1] <- sum(N_at[[x]][,1] * choose_lh$W_a * choose_lh$Mat_a)

		
		## catchMSY
		Cn_at[[x]][,1] <- N_at[[x]][,1] * (1 - exp(-choose_lh$M - F_t[[x]][1] * choose_lh$S_a)) * (F_t[[x]][1] * choose_lh$S_a) / (choose_lh$M + F_t[[x]][1] * choose_lh$S_a)

		Cw_at[[x]][,1] <- choose_lh$W_a * Cn_at[[x]][,1]

	}

		SB0 <- lapply(1:nspecies, function(x){
			choose_lh <- lh[[which(species == species[x])]]
			SB0 <- with(choose_lh, sum(calc_equil_abund(ages = ages, S_a = S_a, M = M, F = 0, R0 = R0) * Mat_a * W_a))
		})

#######################################
## project age structure
#######################################

	for(x in 1:nspecies){
		choose_lh <- lh[[which(species == species[x])]]
		sp <- input_mat %>% dplyr::filter(species == species[x])
		
		for(y in 2:tyears){
			if(sp$Fdynamics == "EffortDyn"){
				if(y <= (nburn+1)) F_t[[x]][y] <- Finit[[x]] * exp(FishDev[[x]][y])
				if(y > (nburn+1)) F_t[[x]][y] <- F_t[[x]][y-1] * (SB_t[[x]][y-1] / (0.2 * SB0[[x]]/2)) ^ 0.2 * exp(FishDev[[x]][y])
			}
			
			R_t[[x]][y] <- (4 * choose_lh$h * choose_lh$R0 * SB_t[[x]][y-1] / (SB0[[x]] * (1 - choose_lh$h) + SB_t[[x]][y-1] * (5 - choose_lh$h - 1))) * exp(RecDev[[x]][y])

			for(a in 1:length(choose_lh$ages)){
				if(a == 1){
					N_at[[x]][a,y] <- R_t[[x]][y]
					N_at0[[x]][a,y] <- R_t[[x]][y]
				}
				if(a > 1 & a < length(choose_lh$ages)){
					N_at[[x]][a,y] <- N_at[[x]][a-1, y-1] * exp(-choose_lh$M - F_t[[x]][y-1] * choose_lh$S_a[a-1])
					N_at0[[x]][a,y] <- N_at0[[x]][a-1,y-1] * exp(-choose_lh$M)
				}
				if(a == length(choose_lh$ages)){
					N_at[[x]][a,y] <- (N_at[[x]][a-1,y-1] * exp(-choose_lh$M - F_t[[x]][y-1] * choose_lh$S_a[a-1])) + (N_at[[x]][a, y-1] * exp(-choose_lh$M - F_t[[x]][y-1] * choose_lh$S_a[a]))
					N_at0[[x]][a,y] <- (N_at0[[x]][a-1,y-1] * exp(-choose_lh$M)) + (N_at0[[x]][a, y-1] * exp(-choose_lh$M))
				}
			}		

			SB_t[[x]][y] <- sum(N_at[[x]][,y] * choose_lh$W_a * choose_lh$Mat_a)
			VB_t[[x]][y] <- sum(N_at[[x]][,y] * choose_lh$W_a * choose_lh$S_a)
			TB_t[[x]][y] <- sum(N_at[[x]][,y] * choose_lh$W_a)		

			Cn_at[[x]][,y] <- N_at[[x]][,y] * (1 - exp(-choose_lh$M - F_t[[x]][y] * choose_lh$S_a)) * (F_t[[x]][y] * choose_lh$S_a) / (choose_lh$M + F_t[[x]][y] * choose_lh$S_a)
			Cw_at[[x]][,y] <- Cn_at[[x]][,y] * choose_lh$W_a
		}
	}

	trim <- lapply(1:nspecies, function(x){
		out <- NULL
		out$N_at <- N_at[[x]][,-c(1:nburn)]
		out$N_t <- colSums(out$N_at)
		out$Cn_at <- Cn_at[[x]][,-c(1:nburn)]
		out$Cw_at <- Cw_at[[x]][,-c(1:nburn)]
		out$SB_t <- SB_t[[x]][-c(1:nburn)]
		out$VB_t <- VB_t[[x]][-c(1:nburn)]
		out$TB_t <- TB_t[[x]][-c(1:nburn)]
		out$D_t <- SB_t[[x]][-c(1:nburn)]/SB0[[x]]
		out$F_t <- F_t[[x]][-c(1:nburn)]
		out$R_t <- R_t[[x]][-c(1:nburn)]
		return(out)
	})
	return(trim)

}