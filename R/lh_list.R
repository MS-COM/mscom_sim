lh_list <- function(linf, vbk, t0, lwa, lwb, M, M50, S50, S95, h, R0){

	Amax <- ceiling(-log(0.01)/M)
	ages <- 0:Amax

	L_a <- linf*(1 - exp(-vbk*(ages - t0)))
	W_a <- lwa * L_a ^ lwb

	Mat_a <- 1 / (1 + exp(M50 - ages))

	S_a <- (1 / (1 + exp(-log(19) * (ages - S50) / (S95 - S50))))

	out <- NULL
	out$Amax <- Amax
	out$ages <- ages
	out$L_a <- L_a
	out$W_a <- W_a
	out$Mat_a <- Mat_a
	out$M <- M
	out$h <- h
	out$S_a <- S_a
	out$R0 <- R0
	return(out)

}