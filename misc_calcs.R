	# L_a_bill <- 275*(1 - exp(-0.25*(0 - -0.01)))
	# R0_bill <- 227 / (0.006 * L_a_bill ^ 2.9) ## 102321

	# ages_tuna <- 0:ceiling(-log(0.01)/0.5)
	# L_a_tuna <- 183*(1 - exp(-0.5*(ages_tuna - -0.01)))
	# N0_tuna <- (606000 * exp(-0.5 * ages_tuna))
	# B0_tuna <- sum(N0_tuna * 0.0224 * L_a_tuna ^ 2.94)

	# ages_shark <- 0:ceiling(-log(0.01)/0.17)
	# L_a_shark <- 321*(1 - exp(-0.1*(ages_shark - -0.01)))
	# N0_shark <- 100000 * exp(-0.1 * ages_shark)
	# B0_shark <- sum(N0_shark * 0.00524 * L_a_shark ^ 3.141)

	# L_a_shark <- 321*(1 - exp(-0.1 * (0 - -0.01)))
	# K_shark <- 100000 * 0.00524 * L_a_shark ^ 3.141



						#"M"=c(0.5, 0.47, 0.17),

						#"R0"=c(606000, 102321, 100000),
						# "linf"=c(183, 275, 321), "vbk"=c(0.5, 0.25, 0.1), "t0"=-0.01,
						# "lwa"=c(0.0224, 0.0066, 0.00524), "lwb"=c(2.94, 2.9, 3.141),
						# "M50"=c(2, 2, 17), ## billfish = 150 mm
						# "S1"=c(3, 3, 6), "S2"=c(4, 10, 11), ## mako -- peak age 8, 50% at age 6 and 11, billfish = knife edge at age 10
						# "Selex_type"=c('logistic','logistic','logistic'),
						# "h"=c(0.75, 0.9, 0.4),
