# *******************************************
# R SCRIPT FOR NITROGEN ISOTOPOMER ANALYSIS #
# Graphical Part                            #
# M. Barthel |  24 Sep 2015                 #
# version: 0.1 - 24 Jan 2018                #
# adopted for GitHub XXX 2018 (R. Hüppi)    #
#############################################
library(nlme)   # package for function lmList
library(data.table)  # data table functions (largely compatible with R-base dataframes)
library(lubridate)   # round date vector in order to average the data vector
library(corrplot)
library(plotly)
# introduce Mattis cool colour codes
kulergrey   = rgb(64,55,47, maxColorValue =  255)
kulerorange = rgb(217,86,28, maxColorValue = 255)

dream1      = rgb(64,1,13, maxColorValue =   255)
dream2      = rgb(217,4,82, maxColorValue =  255)
dream3      = rgb(140,131,3, maxColorValue = 255)
dream4      = rgb(89,43,2, maxColorValue =   255)
dream5      = rgb(217,130,54, maxColorValue =255)

## read in data from QCL-main.R script
# FLUXDATA <-  fread("output/FLUXDATA.csv", header=T, sep = ",")
# QCL.DATA     <-  fread("output/QCL-DATA.csv", header=T, sep = ",")
# QCL <- QCL.DATA
# KEELING.RES     <-  fread("output/keeling-results.csv", header=T, sep = ",")

NA.RM.1  <- !is.na(QCL.DATA$calcode) & !is.na(QCL$spec.446a) & # !is.na(QCL$spec.546a) & !is.na(QCL$spec.456a) & #  & !is.na(QCL$d15Na)    & !is.na(QCL$d15Nb) & !is.na(QCL$d18O) &
  QCL$calcode != c(2740) & QCL$calcode != 2750 # take out certain dates where the dilution cal didn't work


# OFFSET CORRECTION USING ETHZSAEHIGH1 at VICI port #3
OFFSET <- NA.RM.1 & QCL$VICI == 3 & QCL$hourinseconds > 250 &  QCL$hourinseconds < 310 & QCL$hour == 3 |  # last minute of cal cycle at VICI#3 before dilution starts, 
  NA.RM.1 & QCL$VICI == 3 & QCL$hourinseconds > 250 &  QCL$hourinseconds < 310 & QCL$hour == 15

#selection vectors
SPAN   <- QCL.DATA$VICI == 4 & QCL.DATA$hourinseconds > 480 & QCL.DATA$TIMESTAMP > as.POSIXct("2015-09-30 00:00:00, tz=UTC")  # 2015-12-05 is the actual starting date
ANCHOR <- NA.RM.1 & QCL.DATA$VICI == 3 & QCL.DATA$hourinseconds > 250 &  QCL.DATA$hourinseconds < 310 & QCL.DATA$hour == 3  |  # & QCL.DATA$TIMESTAMP > as.POSIXct("2015-12-05 00:00:00, tz=UTC") |  # last minute of cal cycle at VICI#3 before dilution starts, 
          NA.RM.1 & QCL.DATA$VICI == 3 & QCL.DATA$hourinseconds > 250 &  QCL.DATA$hourinseconds < 310 & QCL.DATA$hour == 15  # & QCL.DATA$TIMESTAMP > as.POSIXct("2015-12-05 00:00:00, tz=UTC")

# NA.RM.1  <- !is.na(QCL.DATA$calcode) & !is.na(QCL.DATA$spec.446a)   & !is.na(QCL.DATA$d15Na)    & !is.na(QCL.DATA$d15Nb) & !is.na(QCL.DATA$d18O)

DILUTION <- NA.RM.1 & QCL.DATA$VICI == 3 & QCL.DATA$hourinseconds > 310 & QCL.DATA$hour == 3 | 
            NA.RM.1 & QCL.DATA$VICI == 3 &
##############################################################################################################

#########################
### publication Fig 1 ###
#########################
## dilution, anchor and other data validation plots
## Dilution validation plot ##   

pdf("graphs/fig1_dilution.pdf",encoding="MacRoman",height=10,width=12,family="Times",fonts="Times")   
par(mfrow=c(3,3),bg="white",pch=16,tcl=0.3,cex=1,las=0,srt=0,cex.lab=1*1.3,xpd=FALSE)
par(omi = c(0.5,0,0,0.5))

par(mai = c(0.5,1,0.25,0))

plot(QCL.DATA$spec.446a[DILUTION], QCL.DATA$d15Na[DILUTION], xlim =c (0,3000), ylim = c(-60,10),
     xlab = expression("[446] in ppb"), ylab = expression(delta~"456"), col = kulergrey )  # expression(delta^15*"N"^beta)
abline(lm(QCL.DATA$d15Na[DILUTION]~QCL.DATA$spec.446a[DILUTION]))
points(QCL.DATA$spec.446a[DILUTION],QCL.DATA$d15Na.corr2[DILUTION],col=dream2)
abline(-1.1,0)
legend("bottomright",c("raw data","dilution correction"),
       col=c(kulergrey,dream2),pch=16,bty="n")

plot(QCL.DATA$spec.446a[DILUTION], QCL.DATA$d15Nb[DILUTION], xlim =c (0,3000),  ylim = c(-40,10),
     xlab = expression("[446] in ppb"), ylab = expression(delta~"546"), col = kulergrey )  # expression(delta^15*"N"^beta)
abline(lm(QCL.DATA$d15Nb[DILUTION]~QCL.DATA$spec.446a[DILUTION]))
points(QCL.DATA$spec.446a[DILUTION],QCL.DATA$d15Nb.corr2[DILUTION],col=dream3)
abline(2.2,0)
legend("bottomright",c("raw data","dilution correction"),
       col=c(kulergrey,dream3),pch=16,bty="n")

mtext(expression("[446] in ppb"),1,line=2.5, cex = 1.2)

plot(QCL.DATA$spec.446a[DILUTION], QCL.DATA$d18O[DILUTION], xlim =c (0,3000),  ylim = c(0,60),
     xlab = expression("[446] in ppb"), ylab = expression(delta~"448"), col = kulergrey )  # expression(delta^15*"N"^beta)
abline(lm(QCL.DATA$d18O[DILUTION]~QCL.DATA$spec.446a[DILUTION]))
points(QCL.DATA$spec.446a[DILUTION],QCL.DATA$d18O.corr2[DILUTION],col=dream5)
abline(40.4,0)
legend("bottomright",c("raw data","dilution correction"),
       col=c(kulergrey,dream5),pch=16,bty="n")


par(mai = c(0,1,0.5,0))


plot(calcode/10,offset.d15Na,  col = dream2, ylab = expression(delta~"456 intercept"),xaxt="n", xlab = "", 
     ylim = c(max(offset.d15Na)*0.8,min(offset.d15Na)*1.2) )  # check the variability of the intercept  parse_date_time(floor(calcode/10),j)  ylim = c(-70,-30),
# plot(calcode,slope.d15Na, ylim = c(0,0.02))  # check the variability of the slope
plot(calcode/10,offset.d15Nb, col = dream3, ylab = expression(delta~"546 intercept"),xaxt="n", xlab = "",
     ylim = c(max(offset.d15Nb)*0.8,min(offset.d15Nb)*1.2) )  # check the variability of the slope 

plot(calcode/10,offset.d18O, col = dream5, ylab = expression(delta~"448 intercept"),xaxt="n", xlab = "",  # check the variability of the intercept   ylim = c(0,60), 
     ylim = c(max(offset.d18O)*1.2,min(offset.d18O)*0.8) )  

par(mai = c(0.5,1,0,0))

plot(calcode/10,slope.d15Na,  col = dream2, ylab = expression(delta~"456 slope"),  # check the variability of the intercept  ylim = c(-70,-30),
     ylim = c(max(slope.d15Na)*1.2,min(slope.d15Na)*0.8) )  
plot(calcode/10,slope.d15Nb, col = dream3, ylab = expression(delta~"546 slope"),  # check the variability of the intercept   ylim = c(-35,-15),
     ylim = c(max(slope.d15Nb)*1.2,min(slope.d15Nb)*0.8) )  

     mtext("Date [DOY]",1,line=2.2, cex = 1.2)

plot(calcode/10,slope.d18O, col = dream5, ylab = expression(delta~"448 slope"),  # check the variability of the intercept   ylim = c(0,60), 
     ylim = c(max(slope.d18O)*0.8,min(slope.d18O)*1.2) )  


dev.off()

# # pdf("graphs/anchor.pdf",encoding="MacRoman",height=6,width=12,family="Times",fonts="Times")   
# par(mfrow=c(1,1),bg="white",pch=16,tcl=0.3,cex=1,las=0,srt=0,cex.lab=1*1.3,xpd=FALSE)
# par(omi = c(0.7,0,0.09,0.09))
# par(mai = c(1,1,0.1,0))
# 
# plot(QCL.DATA$TIMESTAMP[ANCHOR],QCL.DATA$d15Nb.corr2[ANCHOR],ylim=c(-20,60),col="white",xlab='Time [months]', ylab=expression(paste(delta^{15},N[beta])),
#      xlim=c(as.POSIXct("2015-10-30 00:00:00, tz=UTC"),as.POSIXct("2016-01-28 00:00:00, tz=UTC")))  # 2015-12-05 is the actual starting date
# 
# points(QCL.DATA$TIMESTAMP[ANCHOR],QCL.DATA$d15Nb.corr2[ANCHOR],col=dream1,cex=0.3,pch=16)
# abline(h=ETHZSAEHIGH1.d15Nb,col=dream1)
# 
# points(QCL.DATA$TIMESTAMP[SPAN],QCL.DATA$d15Nb.corr2[SPAN],col=dream2,cex=0.3,pch=16)
# abline(h=ETHZSAEHIGH7.d15Nb,col=dream2)
# # dev.off()


#########################
### publication Fig 2 ###
#########################
# ANCHOR (ETHZSAEHIGH-1) & SPAN (ETHZSAEHIGH-7) TARGET CHECK

pdf("graphs/fig2_anchor-span.pdf",encoding="WinAnsi", height=6,width=12,family="Times",fonts="Times")  # ,encoding="MacRoman"
par(mfrow=c(4,1),bg="white",pch=16,tcl=0.3,cex=1,las=0,srt=0,cex.lab=1,xpd=F)
par(omi = c(0.75,0,0.25,0.25)) # graph margin
par(mai = c(0,1,0,0))  # plot margin

# x.time  <- 288:310   # 274:304  
# x.time  <- min(QCL.DATA$DOY[SPAN]):max(QCL.DATA$DOY[SPAN])
x.time  <- na.omit(unique(QCL.DATA$DOY[SPAN]))

# d15Na
plot(x.time,tapply(QCL.DATA$d15Na.corr2[SPAN],QCL.DATA$DOY[SPAN],mean),col="white",xlab="", xaxt = "n", ylim=c(32.5,38.5),
     ylab= expression(delta~"456 [\u2030]")) #  sorry, this is different on other plattforms, i.e. \211 for windows?!

abline(h = ETHZSAEHIGH7.d15Na, col = dream2)
rect( min(x.time-10), ETHZSAEHIGH7.d15Na-0.71 ,max(x.time+10),ETHZSAEHIGH7.d15Na+0.71, col = rgb(217,4,82, alpha = 50, maxColorValue =  255), border = F)

sd.d15Na   <- tapply(QCL.DATA$d15Na.corr2[SPAN],QCL.DATA$DOY[SPAN],sd)
mean.d15Na <- tapply(QCL.DATA$d15Na.corr2[SPAN],QCL.DATA$DOY[SPAN],mean)
# n.d15Na    <- tapply(QCL.DATA$d15Na.corr2[SPAN],QCL.DATA$DOY[SPAN],length)

points(x.time ,mean.d15Na,col=dream2)

for (i in 1:length(x.time )){
  arrows(x.time[i], mean.d15Na[i] + sd.d15Na[i], x.time[i], mean.d15Na[i] - sd.d15Na[i],
         angle=90, code=3, length=0.03, col=kulergrey, lwd=0.4) # Fehlerbalken 
}
# legend("topright",c("456"), cex = 1, 
#        col=c(dream2),pch=16,bty="n")

# d15Nb
plot(x.time,tapply(QCL.DATA$d15Nb.corr2[SPAN],QCL.DATA$DOY[SPAN],mean),col="white",xlab="", xaxt = "n", yaxt = "n", ylim=c(32,38),
     ylab = expression(delta~"546 [\u2030]")) #  sorry, this is different on other plattforms, i.e. \211 for windows?!

abline(h = ETHZSAEHIGH7.d15Nb, col = dream3)
rect( min(x.time-10), ETHZSAEHIGH7.d15Nb-0.5 ,max(x.time+10),ETHZSAEHIGH7.d15Nb+0.5, col = rgb(140,131,3, alpha = 50, maxColorValue =  255), border = F)
axis(2, at = seq(33,38,1))

sd.d15Nb   = tapply(QCL.DATA$d15Nb.corr2[SPAN],QCL.DATA$DOY[SPAN],sd)
mean.d15Nb = tapply(QCL.DATA$d15Nb.corr2[SPAN],QCL.DATA$DOY[SPAN],mean)

points(x.time ,mean.d15Nb,col=dream3)

for (i in 1:length(x.time )){
  arrows(x.time[i], mean.d15Nb[i] + sd.d15Nb[i], x.time[i], mean.d15Nb[i] - sd.d15Nb[i],
         angle=90, code=3, length=0.03, col=kulergrey, lwd=0.4) # Fehlerbalken 
}
# legend("topright",c("546"), cex = 1, 
#        col=c(dream3),pch=16,bty="n")

# SP
QCL.DATA$SP.corr2 <- QCL.DATA$d15Na.corr2 - QCL.DATA$d15Nb.corr2
plot(x.time,tapply(QCL.DATA$SP.corr2[SPAN],QCL.DATA$DOY[SPAN],mean,na.omit = TRUE),col="white",xlab="", xaxt = "n", yaxt = "n", ylim=c(-3.5,2.75),
     ylab = expression("SP [\u2030]")) #  sorry, this is different on other plattforms, i.e. \211 for windows?!

abline(h = ETHZSAEHIGH7.d15Na - ETHZSAEHIGH7.d15Nb, col = "green")
rect( min(x.time-10), ETHZSAEHIGH7.d15Nb-0.5 ,max(x.time+10),ETHZSAEHIGH7.d15Nb+0.5, col = rgb(140,131,3, alpha = 50, maxColorValue =  255), border = F)
axis(2, at = seq(-3,2,1))

sd.d15Nb   = tapply(QCL.DATA$SP.corr2[SPAN],QCL.DATA$DOY[SPAN],sd)
mean.d15Nb = tapply(QCL.DATA$SP.corr2[SPAN],QCL.DATA$DOY[SPAN],mean)

points(x.time ,mean.d15Nb,col="green")

for (i in 1:length(x.time )){
  arrows(x.time[i], mean.d15Nb[i] + sd.d15Nb[i], x.time[i], mean.d15Nb[i] - sd.d15Nb[i],
         angle=90, code=3, length=0.03, col="grey50", lwd=0.4) # Fehlerbalken 
}
# legend("topright",c("546"), cex = 1, 
#        col=c(dream3),pch=16,bty="n")

# d18O
plot(x.time,tapply(QCL.DATA$d18O.corr2[SPAN],QCL.DATA$DOY[SPAN],mean),col="white",xlab='time [DOY]',  ylim=c(31,43),
     ylab= expression(delta~"448 [\u2030]")) #  sorry, this is different on other plattforms, i.e. \211 for windows?!

abline(h = ETHZSAEHIGH7.d18O, col = dream5)
rect( min(x.time-10), ETHZSAEHIGH7.d18O-0.35 ,max(x.time+10),ETHZSAEHIGH7.d18O+0.35, col = rgb(217,130,54, alpha = 50, maxColorValue =  255), border = F)

sd.d18O   = tapply(QCL.DATA$d18O.corr2[SPAN],QCL.DATA$DOY[SPAN],sd)
mean.d18O = tapply(QCL.DATA$d18O.corr2[SPAN],QCL.DATA$DOY[SPAN],mean)

points(x.time ,mean.d18O,col=dream5)

for (i in 1:length(x.time )){
  arrows(x.time[i], mean.d18O[i] + sd.d18O[i], x.time[i], mean.d18O[i] - sd.d18O[i],
         angle=90, code=3, length=0.03, col="grey50", lwd=0.4) # Fehlerbalken 
}
# legend("topright",c("448"), cex = 1, 
#        col=c(dream5),pch=16,bty="n")

mtext("Date [DOY]",1,line=2.5, cex = 1.2)

dev.off()

## fig 2.2
# 
# qcl.filter <- QCL %>%  filter(OFFSET)    # %>%  filter(ampmnumeric == 0 )
# 
# var.x <- qcl.filter$TIMESTAMP
# var.y <- qcl.filter$d15Na
# var.y2 <- qcl.filter$Tcell
# fit <-  lm(var.y ~ var.x)
# ay <- list(
#   tickfont = list(color = "red"),
#   overlaying = "y",
#   side = "right",
#   title = "Cell temperature"
# )
# 
# qcl.filter %>% 
#   plot_ly() %>% 
#   # add_markers(x = ~var.x)
#   add_markers(x = ~var.x, y = ~var.y, name = "d15Na" ) %>% 
#   add_markers(x = ~var.x, y = ~var.y2, name = "tcell", yaxis = "y2") %>% 
#   add_lines(x = ~var.x, y = fitted(fit), name = "linear regression")  %>% 
#   layout(showlegend = TRUE,
#     xaxis = list(title = "Time"),
#     yaxis = list(title = "d15Na",
#     title = "Double Y Axis", 
#     yaxis2 = ay)
#   )



#########################
### publication Fig 3 ###
#########################   
pdf("graphs/fig3_keeling.pdf", width=9, height=8, encoding="WinAnsi")  # WinAnsi   MacRoman

par(bg="white",pch=16,tcl=0.3,cex=0.9,cex.lab=1.5,las=0,srt=0,xpd=FALSE)
par(mfrow=c(3,2))
par(omi = c(0.7,0,0.5,0.25))
par(mai = c(0,1,0,0))

chamber2.1 = QCL.DATA$hourinseconds > 900 & QCL.DATA$hourinseconds < 3600 & QCL.DATA$VICI > 4 & QCL.DATA$VICI.time == 2830611 # 2821811 #  October 10 chamber 7
# chamber2.2 = QCL.DATA$hourinseconds > 900 & QCL.DATA$hourinseconds < 3600 & QCL.DATA$VICI > 4 & QCL.DATA$VICI.time == 2830611
# concentrations alpha
###########################################################################################################

plot(QCL.DATA$hourinseconds[chamber2.1],QCL.DATA$spec.446a[chamber2.1],col="white",cex=.5,ylim=c(300,420),xaxt="n",
     ylab = "456 & 446 [ppb]", xlab="Time [sec]")

points(QCL.DATA$hourinseconds[chamber2.1],QCL.DATA$spec.446a[chamber2.1],col=dream1,cex=.5,ylim=c(300,420))
points(QCL.DATA$hourinseconds[chamber2.1],QCL.DATA$spec.456a[chamber2.1],col=dream2,cex=.5,ylim=c(300,420))

legend("bottomright",c("456","446"), cex = 1.2, 
       col=c(dream2,dream1),pch=16,bty="n")
#points(QCL.DATA$hourinseconds[chamber2.2],QCL.DATA$spec.446a[chamber2.2],col=dream1,cex=.5,ylim=c(300,420))
#points(QCL.DATA$hourinseconds[chamber2.2],QCL.DATA$spec.456a[chamber2.2],col=dream2,cex=.5,ylim=c(300,420))

### calculate SD for a single flux ###
# detrend.species <- QCL.DATA$spec.446a[chamber2.1]
# detrend.stats <- summary( lm(detrend.species ~ QCL.DATA$hourinseconds[chamber2.1]) )
# detrend.stats$coefficients[2,1]
# 
# detrended <- detrend.species - detrend.stats$coefficients[2,1] * QCL.DATA$hourinseconds[chamber2.1]
# plot(detrended)
# sd.start <- sd(detrended[1:100])
# sd.end   <- sd(detrended[150:250])

# keeling plots alpha
plot(QCL.DATA$spec.446a.inv[chamber2.1],QCL.DATA$d15Na.corr2[chamber2.1],col="white",cex=.5,ylim=c(-10,40),xlim=c(0.00245,0.00295),xaxt="n",
     ylab= expression(delta^15~"N [\u2030]"))    # expression(paste(delta^{15},N[alpha])))

points(QCL.DATA$spec.446a.inv[chamber2.1],QCL.DATA$d15Na.corr2[chamber2.1],col=dream1,cex=.5)
abline(lm(QCL.DATA$d15Na.corr2[chamber2.1]~QCL.DATA$spec.446a.inv[chamber2.1]))

#points(QCL.DATA$spec.446a.inv[chamber2.2],QCL.DATA$d15Na.corr2[chamber2.2],col=dream4,cex=.5)
#abline(lm(QCL.DATA$d15Na.corr2[chamber2.2]~QCL.DATA$spec.446a.inv[chamber2.2]))


# concentrations beta
###########################################################################################################

plot(QCL.DATA$hourinseconds[chamber2.1],QCL.DATA$spec.446a[chamber2.1],col="white",cex=.5,ylim=c(300,420),xaxt="n", yaxt = "n",
     xlab="Time [sec]", ylab = "546 & 446 [ppb]" )  #ylab=expression(""^{14}*"N"^{14}*"N"^{16}*"O & "^{15}*"N"^{14}*"N"^{16}*"O [ppb]"))

points(QCL.DATA$hourinseconds[chamber2.1],QCL.DATA$spec.446a[chamber2.1],col=dream1,cex=.5,ylim=c(300,420))
points(QCL.DATA$hourinseconds[chamber2.1],QCL.DATA$spec.546a[chamber2.1],col=dream3,cex=.5,ylim=c(300,420))

legend("bottomright",c("546","446"),cex = 1.2, 
       col=c(dream3,dream1),pch=16,bty="n")
axis(2, at = seq(320,400,20))
#points(QCL.DATA$hourinseconds[chamber2.2],QCL.DATA$spec.446a[chamber2.2],col=dream1,cex=.5,ylim=c(300,420))
#points(QCL.DATA$hourinseconds[chamber2.2],QCL.DATA$spec.546a[chamber2.2],col=dream3,cex=.5,ylim=c(300,420))


# keeling plots beta
plot(QCL.DATA$spec.446a.inv[chamber2.1],QCL.DATA$d15Nb.corr2[chamber2.1],col="white",cex=.5,ylim=c(-30,10),xlim=c(0.00245,0.00295),xaxt="n", yaxt = "n",
     ylab = expression(delta^15~"N [\u2030]"))

points(QCL.DATA$spec.446a.inv[chamber2.1],QCL.DATA$d15Nb.corr2[chamber2.1],col=dream1,cex=.5)
abline(lm(QCL.DATA$d15Nb.corr2[chamber2.1]~QCL.DATA$spec.446a.inv[chamber2.1]))
axis(2, at = seq(-30,10,10))

#points(QCL.DATA$spec.446a.inv[chamber2.2],QCL.DATA$d15Nb.corr2[chamber2.2],col=dream4,cex=.5)
#abline(lm(QCL.DATA$d15Nb.corr2[chamber2.2]~QCL.DATA$spec.446a.inv[chamber2.2]))


# concentrations 18O
###########################################################################################################

plot(QCL.DATA$hourinseconds[chamber2.1],QCL.DATA$spec.446a[chamber2.1],col="white",cex=.5,ylim=c(300,420),
     ylab="448 and 446 [ppb]", xlab="Time [sec]")

points(QCL.DATA$hourinseconds[chamber2.1],QCL.DATA$spec.446a[chamber2.1],col=dream1,cex=.5,ylim=c(300,420))
points(QCL.DATA$hourinseconds[chamber2.1],QCL.DATA$spec.448a[chamber2.1],col=dream5,cex=.5,ylim=c(300,420))

legend("bottomright",c("448","446"), cex = 1.2, 
       col=c(dream5,dream1),pch=16,bty="n")
#points(QCL.DATA$hourinseconds[chamber2.2],QCL.DATA$spec.446a[chamber2.2],col=dream1,cex=.5,ylim=c(300,420))
#points(QCL.DATA$hourinseconds[chamber2.2],QCL.DATA$spec.448a[chamber2.2],col=dream5,cex=.5,ylim=c(300,420))


mtext("chamber closure time [s]",1,line=3)

# keeling plots 18O
plot(QCL.DATA$spec.446a.inv[chamber2.1],QCL.DATA$d18O.corr2[chamber2.1],col="white",cex=.5,ylim=c(-50,100),xlim=c(0.00245,0.00295),
     ylab = expression(delta^15~"N [\u2030]"))

points(QCL.DATA$spec.446a.inv[chamber2.1],QCL.DATA$d18O.corr2[chamber2.1],col=dream1,cex=.5)
abline(lm(QCL.DATA$d18O.corr2[chamber2.1]~QCL.DATA$spec.446a.inv[chamber2.1]))

mtext("1/[446]",1,line=3)

#points(QCL.DATA$spec.446a.inv[chamber2.2],QCL.DATA$d18O.corr2[chamber2.2],col=dream4,cex=.5)
#abline(lm(QCL.DATA$d18O.corr2[chamber2.2]~QCL.DATA$spec.446a.inv[chamber2.2]))

dev.off()


#########################
### publication Fig 4 ###
#########################

pdf("graphs/fig4-1_flux.pdf",encoding="MacRoman",height=3,width=7,family="Times",fonts="Times")
par(mfrow=c(1,1),bg="white",pch=16,tcl=0.3,cex=0.8,las=0,srt=0,cex.lab=1,xpd=F)
par(omi = c(0.7,0,0.1,0.1))
par(mai = c(0,0.75,0.1,0))

plot(ch.time,n2o.flux,
    ylab=expression("flux N"[2]*"O [nmol m"^{-2}*" s"^{-1}*"]"),xlab='Time',#ylim = c(-0.1,1.2),
   col=dream1,cex=.5)
abline(h=0,lty=3)

dev.off()

pdf("graphs/fig4-2_fluxvsdelta.pdf",encoding="WinAnsi",height=4,width=7,family="Times",fonts="Times")
par(mfrow=c(1,3),bg="white",pch=16,tcl=0.3,cex=1,las=0,srt=0,cex.lab=1,xpd=F)
par(omi = c(0.7,0.5,0.09,0.1))
par(mai = c(0,0.5,0.1,0))

plot(n2o.flux,coef(keeling.d15Na)[,1],
     xlab=expression("N"[2]*"O flux [nmol m"^{-2}*" s"^{-1}*"]"),
     ylab = "",
     # ylab=expression(delta~"[per mille]"),
     cex=.5,xlim=c(0,1.1),ylim=c(-1000,1000),col=dream2)
# rect( 0, -50 ,1.1, 50, col = rgb(140,131,3, alpha = 50, maxColorValue =  255), border = F)
legend("topright",c("456"),
       col=c(dream2),pch=16,bty="n")
abline(h=0)

mtext(side = 2, expression(delta^15~"N [\u2030]"), 2.5)


plot(n2o.flux,coef(keeling.d15Nb)[,1],
    xlab=expression("N"[2]*"O flux [nmol m"^{-2}*" s"^{-1}*"]"),
    ylab = "",
    # ylab=expression(delta~"[per mille]"),
    cex=.5,xlim=c(0,1.1),ylim=c(-1000,1000),col=dream3)
# rect( 0, -50 ,1.1, 50, col = rgb(140,131,3, alpha = 50, maxColorValue =  255), border = F)
legend("topright",c("546"),
       col=c(dream3),pch=16,bty="n")
mtext(side = 1, expression("N"[2]*"O flux [nmol m"^{-2}*" s"^{-1}*"]"), 2.5)

abline(h=0)

plot(n2o.flux,coef(keeling.d18O)[,1],
     xlab=expression("N"[2]*"O flux [nmol m"^{-2}*" s"^{-1}*"]"),
     ylab = "",
     # ylab=expression(delta~"[per mille]"),
     cex=.5,xlim=c(0,1.1),ylim=c(-1000,1000),col=dream5)
# rect( 0, -50 ,1.1, 50, col = rgb(140,131,3, alpha = 50, maxColorValue =  255), border = F)
legend("topright",c("448"),
       col=c(dream5),pch=16,bty="n")
# mtext(side = 2, expression(delta~"[\u2030]"), 2)

abline(h=0)

dev.off()





# pdf("/Users/mbarthel/Desktop/QCL1/graphs/cumflux.pdf",encoding="MacRoman",height=4,width=7,family="Times",fonts="Times")   
par(bg="white",pch=16,tcl=0.3,cex=1,las=0,srt=0,cex.lab=1*1.3,xpd=FALSE)
par(omi = c(0.7,0,0.09,0.09))
par(mai = c(1,1,0.1,0))

plot(ch.time[ch12],cumsum(n2o.flux[ch12]),col="white",ylim=c(0,40),xlab="Time",ylab=expression("flux N"[2]*"O [nmol m"^{-2}*" s"^{-1}*"]"))
    
 
    lines(ch.time[ch01][0:200],cumsum(n2o.flux[ch01][0:200]),col="turquoise1",cex=0.50) # Zinal
    lines(ch.time[ch02][0:200],cumsum(n2o.flux[ch02][0:200]),col="yellow1"   ,cex=0.50) # Probus
    lines(ch.time[ch03][0:200],cumsum(n2o.flux[ch03][0:200]),col="purple1"   ,cex=0.50) # CH Claro
    lines(ch.time[ch04][0:200],cumsum(n2o.flux[ch04][0:200]),col="orange1"   ,cex=0.50) # Monte Calme
    lines(ch.time[ch05][0:200],cumsum(n2o.flux[ch05][0:200]),col="yellow1"   ,cex=0.75) # Probus
    lines(ch.time[ch06][0:200],cumsum(n2o.flux[ch06][0:200]),col="purple1"   ,cex=0.75) # CH Claro
    lines(ch.time[ch07][0:200],cumsum(n2o.flux[ch07][0:200]),col="orange1"   ,cex=0.75) # Monte Calme
    lines(ch.time[ch08][0:200],cumsum(n2o.flux[ch08][0:200]),col="purple1"   ,cex=1.00) # CH Claro
    lines(ch.time[ch09][0:200],cumsum(n2o.flux[ch09][0:200]),col="yellow1"   ,cex=1.00) # Probus
    lines(ch.time[ch10][0:200],cumsum(n2o.flux[ch10][0:200]),col="turquoise1",cex=0.75) # Zinal
    lines(ch.time[ch11][0:200],cumsum(n2o.flux[ch11][0:200]),col="turquoise1",cex=1.00) # Zinal
    lines(ch.time[ch12][0:200],cumsum(n2o.flux[ch12][0:200]),col="orange1"   ,cex=1.00) # Monte Calme

# dev.off()





# # graphical illustration of R2 obtained from regression for flux calculations
# 
# plot(ch.time,n2o.R2,col="white",ylim=c(0,1),xlab="Time",ylab=expression("N"[2]*"O R"^{2}*" from slope"))
#     
#     points(ch.time[ch01][0:200],n2o.R2[ch01][0:200],col="turquoise1",cex=0.50) # Zinal
#     points(ch.time[ch02][0:200],n2o.R2[ch02][0:200],col="yellow1"   ,cex=0.50) # Probus
#     points(ch.time[ch03][0:200],n2o.R2[ch03][0:200],col="purple1"   ,cex=0.50) # CH Claro
#     points(ch.time[ch04][0:200],n2o.R2[ch04][0:200],col="orange1"   ,cex=0.50) # Monte Calme
#     points(ch.time[ch05][0:200],n2o.R2[ch05][0:200],col="yellow1"   ,cex=0.75) # Probus
#     points(ch.time[ch06][0:200],n2o.R2[ch06][0:200],col="purple1"   ,cex=0.75) # CH Claro
#     points(ch.time[ch07][0:200],n2o.R2[ch07][0:200],col="orange1"   ,cex=0.75) # Monte Calme
#     points(ch.time[ch08][0:200],n2o.R2[ch08][0:200],col="purple1"   ,cex=1.00) # CH Claro
#     points(ch.time[ch09][0:200],n2o.R2[ch09][0:200],col="yellow1"   ,cex=1.00) # Probus
#     points(ch.time[ch10][0:200],n2o.R2[ch10][0:200],col="turquoise1",cex=0.75) # Zinal
#     points(ch.time[ch11][0:200],n2o.R2[ch11][0:200],col="turquoise1",cex=1.00) # Zinal
#     points(ch.time[ch12][0:200],n2o.R2[ch12][0:200],col="orange1"   ,cex=1.00) # Monte Calme
#  
# 
#  plot(ch.time,co2.R2,col="white",ylim=c(0,1),xlab="Time",ylab=expression("CO"[2]*" R"^{2}*" from slope"))
#     
#     points(ch.time[ch01][0:200],co2.R2[ch01][0:200],col="turquoise1",cex=0.50) # Zinal
#     points(ch.time[ch02][0:200],co2.R2[ch02][0:200],col="yellow1"   ,cex=0.50) # Probus
#     points(ch.time[ch03][0:200],co2.R2[ch03][0:200],col="purple1"   ,cex=0.50) # CH Claro
#     points(ch.time[ch04][0:200],co2.R2[ch04][0:200],col="orange1"   ,cex=0.50) # Monte Calme
#     points(ch.time[ch05][0:200],co2.R2[ch05][0:200],col="yellow1"   ,cex=0.75) # Probus
#     points(ch.time[ch06][0:200],co2.R2[ch06][0:200],col="purple1"   ,cex=0.75) # CH Claro
#     points(ch.time[ch07][0:200],co2.R2[ch07][0:200],col="orange1"   ,cex=0.75) # Monte Calme
#     points(ch.time[ch08][0:200],co2.R2[ch08][0:200],col="purple1"   ,cex=1.00) # CH Claro
#     points(ch.time[ch09][0:200],co2.R2[ch09][0:200],col="yellow1"   ,cex=1.00) # Probus
#     points(ch.time[ch10][0:200],co2.R2[ch10][0:200],col="turquoise1",cex=0.75) # Zinal
#     points(ch.time[ch11][0:200],co2.R2[ch11][0:200],col="turquoise1",cex=1.00) # Zinal
#     points(ch.time[ch12][0:200],co2.R2[ch12][0:200],col="orange1"   ,cex=1.00) # Monte Calme
#     


##############################################################################################################
# cumulative sum of all treatments CO2

plot(ch.time[ch01][0:200],ZINAL.co2.flux/1000,col="white",ylim=c(-150,50),xlab="Time",ylab=expression("flux CO"[2]*" ["*mu*"mol m"^{-2}*" s"^{-1}*"]"))

lines(ch.time[ch01][0:200],cumsum(ZINAL.co2.flux)/1000, col="turquoise1",cex=0.50) # Zinal
lines(ch.time[ch01][0:200],cumsum(PROBUS.co2.flux)/1000,col="yellow1"   ,cex=0.50) # Probus
lines(ch.time[ch01][0:200],cumsum(CLARO.co2.flux)/1000, col="purple1"   ,cex=0.50) # CH Claro
lines(ch.time[ch01][0:200],cumsum(MONTE.co2.flux)/1000, col="orange1"   ,cex=0.50) # Monte Calme
abline(h=0,lty=2)


# cumulative N2O sum of single chamber measurements, no difference between using all data or data with good R2 only.

plot(ch.time[ch12][0:200],cumsum(n2o.flux[ch12][0:200]),col="white",ylim=c(0,40),xlab="Time",ylab=expression("flux N"[2]*"O [nmol m"^{-2}*" s"^{-1}*"]"))

lines(ch.time[ch01 & n2o.R2 > 0.7][0:200],cumsum(n2o.flux[ch01 & n2o.R2 > 0.7][0:200]),col="turquoise1",cex=0.50) # Zinal
lines(ch.time[ch02 & n2o.R2 > 0.7][0:200],cumsum(n2o.flux[ch02 & n2o.R2 > 0.7][0:200]),col="yellow1"   ,cex=0.50) # Probus
lines(ch.time[ch03 & n2o.R2 > 0.7][0:200],cumsum(n2o.flux[ch03 & n2o.R2 > 0.7][0:200]),col="purple1"   ,cex=0.50) # CH Claro
lines(ch.time[ch04 & n2o.R2 > 0.7][0:200],cumsum(n2o.flux[ch04 & n2o.R2 > 0.7][0:200]),col="orange1"   ,cex=0.50) # Monte Calme
lines(ch.time[ch05 & n2o.R2 > 0.7][0:200],cumsum(n2o.flux[ch05 & n2o.R2 > 0.7][0:200]),col="yellow1"   ,cex=0.75) # Probus
lines(ch.time[ch06 & n2o.R2 > 0.7][0:200],cumsum(n2o.flux[ch06 & n2o.R2 > 0.7][0:200]),col="purple1"   ,cex=0.75) # CH Claro
lines(ch.time[ch07 & n2o.R2 > 0.7][0:200],cumsum(n2o.flux[ch07 & n2o.R2 > 0.7][0:200]),col="orange1"   ,cex=0.75) # Monte Calme
lines(ch.time[ch08 & n2o.R2 > 0.7][0:200],cumsum(n2o.flux[ch08 & n2o.R2 > 0.7][0:200]),col="purple1"   ,cex=1.00) # CH Claro
lines(ch.time[ch09 & n2o.R2 > 0.7][0:200],cumsum(n2o.flux[ch09 & n2o.R2 > 0.7][0:200]),col="yellow1"   ,cex=1.00) # Probus
lines(ch.time[ch10 & n2o.R2 > 0.7][0:200],cumsum(n2o.flux[ch10 & n2o.R2 > 0.7][0:200]),col="turquoise1",cex=0.75) # Zinal
lines(ch.time[ch11 & n2o.R2 > 0.7][0:200],cumsum(n2o.flux[ch11 & n2o.R2 > 0.7][0:200]),col="turquoise1",cex=1.00) # Zinal
lines(ch.time[ch12 & n2o.R2 > 0.7][0:200],cumsum(n2o.flux[ch12 & n2o.R2 > 0.7][0:200]),col="orange1"   ,cex=1.00) # Monte Calme



# R2 in relation to flux. Fluxes close to zero having a bad R2    
plot(n2o.R2,n2o.flux,ylim=c(-.5,2))
plot(co2.R2,co2.flux)   






