#tabla 2x2 del algoritmo

cm <- data.frame(a = 16, b = 120, c = 9, d = 180)

#numero total de enfermos y sanos
enf <- cm$a + cm$c
san <- cm$d + cm$b

#conteos boostrap
a_vec <- rbinom(1000, enf , cm$a/enf)
c_vec <- enf - a_vec

d_vec <- rbinom(1000, san, cm$d/san)
b_vec <- san - d_vec

#sensibilidades y especificidades boostrap
se <- a_vec/enf
es <- d_vec/san

#region de credibilidad para una prevalencia de 0.1
prev <- 0.1
vpp <- (prev*se)/(prev*se+(1-prev)*(1-es))

vpn <- ((1-prev)*es)/((1-prev)*es + (1-se)*prev)

lr_pos <- se/(1-es)
lr_neg <- (1-se)/es
#intervalos de confianza valores predictivos
ic_vpp <- quantile(vpp, probs = c(0.025, 0.975));ic_vpp
ic_vpn <- quantile(vpn, probs = c(0.025, 0.975));ic_vpn

#intervalos de confianza razones de verosimilitud

ic_lr_pos <- quantile(lr_pos, probs = c(0.025, 0.975));ic_lr_pos
ic_lr_neg <- quantile(lr_neg, probs = c(0.025, 0.975));ic_lr_neg



#### simulación de la prevalencia #######

prev <- seq(0, 1, by = 0.001)
se <- cm$a/(cm$a+cm$c)
es <- cm$d/(cm$d+cm$b)

vpp <- (prev*se)/(prev*se+(1-prev)*(1-es))
vpn <- ((1-prev)*es)/((1-prev)*es + (1-se)*prev)

#grafico
plot(prev, vpp, col = "red", xlim = c(min(prev), max(prev)), ylim = c(0,1), type = "l", xlab = "Prevalencia", ylab = "Valor predictivo")
lines(prev, vpn, col = "blue")
legend("right", legend = c("VPP", "VPN"), lty = c(1,1),col = c("red", "blue"), box.lty = 1, cex = 0.55, text.font=4, xjust = 0)

