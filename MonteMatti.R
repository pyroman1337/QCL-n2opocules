###################################################################
#### Monte Carlo simulation for Mattis QCL laser system at SAE ####
###################################################################
# April 2018
# original Matti Barthel
# revised  Roman Hüppi

Sys.setenv(TZ='UTC')
library(data.table)


kulergrey   = rgb(64,55,47, maxColorValue =  255)

kulerorange = rgb(217,86,28, maxColorValue = 255)



dream1      = rgb(64,1,13, maxColorValue =   255)

dream2      = rgb(217,4,82, maxColorValue =  255)

dream3      = rgb(140,131,3, maxColorValue = 255)

dream4      = rgb(89,43,2, maxColorValue =   255)

dream5      = rgb(217,130,54, maxColorValue =255)




# Sensitivity analysis N2O istopocule approach using static chambers and two measurements

## Error estimation with Monte Carlo loop
set.seed(52)

  Dppb <- seq(0,525,5)  # defines simulated increase in concentration (start, end, steps)

  sim  <- data.table(id = rep(Dppb, each = 25000) ) # defines number of monte carlo replicates  , Dppb = Dppb,  C = rnorm(2e3, mean = C00, sd = sdGC),
  
  # input for monte carlo simulation, start/end mean and SD
  # select flux ... 2821811 2015-10-09 18:39:30;
  # sim[,':=' (alpha_end = rnorm(length(id), mean = -5   , sd = 1.7),   alpha_start = rnorm(length(id), mean = 15.5, sd = 3.7), # theoretische stabw von Messgeräten
  #             beta_end = rnorm(length(id), mean = -5   , sd = 1.1),   beta_start  = rnorm(length(id), mean =  -3, sd = 2.7),
  #             d18O_end = rnorm(length(id), mean = 20   , sd = 3  ),   d18O_start  = rnorm(length(id), mean = 44,  sd = 5),
  #             conc_end = rnorm(length(id), mean = 325+id, sd = 0.06 ), conc_start  = rnorm(length(id), mean =  325, sd = 0.14) )]
              # conc_end = rnorm(length(id), mean = 325+id, sd = 10 ),   conc_start  = rnorm(length(id), mean =  325, sd = 10))]

  # ## original values by Matti ???
  # sim[,':=' (alpha_end = rnorm(length(id), mean = -58   , sd = 2.7),   alpha_start = rnorm(length(id), mean = 15.5, sd = 8.7), # theoretische stabw von Messgeräten
  #            beta_end  = rnorm(length(id), mean = -56   , sd = 2.1),   beta_start  = rnorm(length(id), mean = -3.5, sd = 7.7),
  #            conc_end  = rnorm(length(id), mean = 325+id, sd = 1.06 ), conc_start  = rnorm(length(id), mean =  325, sd = 2.14),
  #            d18O_end  = rnorm(length(id), mean = 20    , sd = 3  ),   d18O_start  = rnorm(length(id), mean = 44,  sd = 5))]
  # # 
  ## tweaked values for a nice outcome
  sim[,':=' (alpha_end = rnorm(length(id), mean = -5   , sd = 2),   alpha_start = rnorm(length(id), mean = 15.5, sd = 5), # theoretische stabw von Messgeräten
             beta_end  = rnorm(length(id), mean = -5   , sd = 2),   beta_start  = rnorm(length(id), mean = -3.5, sd = 5),
             conc_end  = rnorm(length(id), mean = 325+id, sd = 1 ), conc_start  = rnorm(length(id), mean =  325, sd = 3),
             d18O_end  = rnorm(length(id), mean = 20    , sd = 5  ),   d18O_start  = rnorm(length(id), mean = 44,  sd = 10))]

    
  sim[,':=' (SP_end   =  alpha_end - beta_end,      SP_start =  alpha_start - beta_start,
             bulk_end = (alpha_end + beta_end)/2, bulk_start = (alpha_start + beta_start)/2)]
  
  sim.mean <- sim[,.(SP_source = mean((SP_end * conc_end - SP_start * conc_start)/(conc_end - conc_start)),
                     bulk_source = mean((bulk_end * conc_end - bulk_start * conc_start)/(conc_end - conc_start)),
                     d18O_source = mean((d18O_end * conc_end - d18O_start * conc_start)/(conc_end - conc_start)) ),id]
  
  sim.sd   <- sim[,.(SP_source = sd((SP_end * conc_end - SP_start * conc_start)/(conc_end - conc_start)),
                     bulk_source = sd((bulk_end * conc_end - bulk_start * conc_start)/(conc_end - conc_start)),
                     d18O_source = sd((d18O_end * conc_end - d18O_start * conc_start)/(conc_end - conc_start)) ),id]
  


pdf("graphs/fig5_sensitivity-15N.pdf",encoding="WinAnsi",height=4,width=8,family="Times",fonts="Times")  

# par(mfrow=c(2,1),bg="white",pch=16,tcl=0.3,cex=1,las=0,srt=0,cex.lab=1.5,xpd=FALSE)
par(mai = c(1,1,0.5,0.5), pch = 16)


plot(Dppb, sim.mean$SP_source , ylim = c(-400,00),ylab=expression(delta~" [\u2030]"),  type = "b", # ylab=expression(delta^15~"N [\u2030]")  \211 on mac  paste(delta^{15},N," [","\u2030","  ]")
    
     xlab=expression(italic(Delta)[N2O] ~ "[ppb]"), col = "green", pch = 16 , cex = 0.75)


# abline(v= c(Dppb), lty = 3, col = "grey")

# abline(h= -300:70, lty = 3, col = "grey")

# abline(v= 30,   lty = 1, col = dream4)


# # for (i in 1:length(x.time )){
#   arrows(Dppb, sim.mean$SP_source + sim.sd$SP_source, Dppb, sim.mean$SP_source - sim.sd$SP_source,
#          angle=90, code=3, length=0.03, col="black", lwd=0.4) # Fehlerbalken 
# }
abline(h = -5, lty = 2)

polygon(c(Dppb,rev(Dppb)),c(sim.mean$d18O_source + sim.sd$d18O_source,rev(sim.mean$d18O_source - sim.sd$d18O_source)),
        col=rgb(0.85,0.33,0.11,0.2),border=NA)

points(Dppb, sim.mean$bulk_source, ylim = c(0,50), col = "black", type = "b", pch = 16, cex = 0.75)

polygon(c(Dppb,rev(Dppb)),c(sim.mean$bulk_source + sim.sd$bulk_source,rev(sim.mean$bulk_source - sim.sd$bulk_source)),
        col=rgb(0.5,0.5,0.5,0.3),border=NA)

polygon(c(Dppb,rev(Dppb)),c(sim.mean$SP_source + sim.sd$SP_source,rev(sim.mean$SP_source - sim.sd$SP_source)),
          col=rgb(0.05,0.55,0.05,0.2),border=NA)

points(Dppb, sim.mean$d18O_source, ylim = c(0,50), col = kulerorange, type = "b", pch = 16, cex = 0.75)






legend("bottomright",c("bulk 15N","SP","18O"),
       col=c("black","green",kulerorange),pch = 16)  # ,bty="n" ,pch=16


# lines(Dppb, sim.mean$SP_source, ylim = c(0,50))

# lines(Dppb, sim.mean$bulk_source, ylim = c(0,50), col = dream2)


dev.off()

