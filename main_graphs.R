# *******************************************
# R SCRIPT FOR NITROGEN ISOTOPOMER ANALYSIS #
# Graphical Part                            #
# M. Barthel |  24 Sep 2015                 #
# version: 0.1 - 24 Jan 2018                #
# adopted for GitHub XXX 2018 (R. Hüppi)    #
#############################################

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
# KEELING.RES     <-  fread("output/keeling-results.csv", header=T, sep = ",")

#selection vectors
SPAN   <- QCL.DATA$VICI == 4 & QCL.DATA$hourinseconds > 480 & QCL.DATA$TIMESTAMP > as.POSIXct("2015-10-30 00:00:00, tz=UTC")  # 2015-12-05 is the actual starting date
ANCHOR <- NA.RM.1 & QCL.DATA$VICI == 3 & QCL.DATA$hourinseconds > 250 &  QCL.DATA$hourinseconds < 310 & QCL.DATA$hour == 3  &  # & QCL.DATA$TIMESTAMP > as.POSIXct("2015-12-05 00:00:00, tz=UTC") |  # last minute of cal cycle at VICI#3 before dilution starts, 
          NA.RM.1 & QCL.DATA$VICI == 3 & QCL.DATA$hourinseconds > 250 &  QCL.DATA$hourinseconds < 310 & QCL.DATA$hour == 15  # & QCL.DATA$TIMESTAMP > as.POSIXct("2015-12-05 00:00:00, tz=UTC")

NA.RM.1  <- !is.na(QCL.DATA$calcode) & !is.na(QCL.DATA$spec.446a)   & !is.na(QCL.DATA$d15Na)    & !is.na(QCL.DATA$d15Nb) & !is.na(QCL.DATA$d18O)

DILUTION <- NA.RM.1 & QCL.DATA$VICI == 3 & QCL.DATA$hourinseconds > 310 & QCL.DATA$hour == 3 | 
            NA.RM.1 & QCL.DATA$VICI == 3 & QCL.DATA$hourinseconds > 310 & QCL.DATA$hour == 15

##############################################################################################################

#########################
### publication Fig 1 ###
#########################
## dilution, anchor and other data validation plots
## Dilution validation plot ##   

plot(calcode,offset.d15Nb, ylim = c(-32,-20))  # check the variability of the intercept



plot(QCL.DATA$spec.446a[DILUTION], QCL.DATA$d15Nb[DILUTION], xlim =c (0,3000), ylim = c(-60,10),
    xlab = expression("total "^{14}*"N"^{14}*"N"^{16}*"O signal [?]"), ylab = expression(delta^15*"N"^beta) )
abline(lm(QCL.DATA$d15Nb[DILUTION]~QCL.DATA$spec.446a[DILUTION]))
points(QCL.DATA$spec.446a[DILUTION],QCL.DATA$d15Nb.corr2[DILUTION],col="red")
abline(0,0)
legend("bottomright",c("raw data","dilution correction"),
    col=c("black","red"),pch=21,bty="n")


# ANCHOR (ETHZSAEHIGH-1) & SPAN (ETHZSAEHIGH-7) TARGET CHECK

# pdf("graphs/anchor.pdf",encoding="MacRoman",height=6,width=12,family="Times",fonts="Times")   
par(mfrow=c(1,1),bg="white",pch=16,tcl=0.3,cex=1,las=0,srt=0,cex.lab=1*1.3,xpd=FALSE)
par(omi = c(0.7,0,0.09,0.09))
par(mai = c(1,1,0.1,0))

plot(QCL.DATA$TIMESTAMP[ANCHOR],QCL.DATA$d15Nb.corr2[ANCHOR],ylim=c(-20,60),col="white",xlab='Time [months]', ylab=expression(paste(delta^{15},N[beta])),
     xlim=c(as.POSIXct("2015-10-30 00:00:00, tz=UTC"),as.POSIXct("2016-01-28 00:00:00, tz=UTC")))  # 2015-12-05 is the actual starting date

points(QCL.DATA$TIMESTAMP[ANCHOR],QCL.DATA$d15Nb.corr2[ANCHOR],col=dream1,cex=0.3,pch=16)
abline(h=ETHZSAEHIGH1.d15Nb,col=dream1)

points(QCL.DATA$TIMESTAMP[SPAN],QCL.DATA$d15Nb.corr2[SPAN],col=dream2,cex=0.3,pch=16)
abline(h=ETHZSAEHIGH7.d15Nb,col=dream2)
# dev.off()

#########################
### publication Fig 2 ###
#########################




#########################
### publication Fig 3 ###
#########################   
# pdf("/Users/mbarthel/Desktop/QCL1/graphs/keeling.pdf", width=9, height=8, encoding="WinAnsi")  

par(bg="white",pch=16,tcl=0.3,cex=0.6*1,cex.lab=1*1.5,las=0,srt=0,xpd=FALSE)
par(mfrow=c(3,2))
par(omi = c(0.7,0,0.5,0.05))
par(mai = c(0,1,0,0))

chamber2.1 = QCL.DATA$hourinseconds > 900 & QCL.DATA$hourinseconds < 3600 & QCL.DATA$VICI > 4 & QCL.DATA$VICI.time == 2830611 # October 10 chamber 7
chamber2.2 = QCL.DATA$hourinseconds > 900 & QCL.DATA$hourinseconds < 3600 & QCL.DATA$VICI > 4 & QCL.DATA$VICI.time == 2830611
# concentrations alpha
###########################################################################################################

plot(QCL.DATA$hourinseconds[chamber2.1],QCL.DATA$spec.446a[chamber2.1],col="white",cex=.5,ylim=c(300,420),xaxt="n",
     ylab=expression(""^{14}*"N"^{14}*"N"^{16}*"O & "^{14}*"N"^{15}*"N"^{16}*"O [ppb]"),xlab="Time [sec]")

points(QCL.DATA$hourinseconds[chamber2.1],QCL.DATA$spec.446a[chamber2.1],col=dream1,cex=.5,ylim=c(300,420))
points(QCL.DATA$hourinseconds[chamber2.1],QCL.DATA$spec.456a[chamber2.1],col=dream2,cex=.5,ylim=c(300,420))

#points(QCL.DATA$hourinseconds[chamber2.2],QCL.DATA$spec.446a[chamber2.2],col=dream1,cex=.5,ylim=c(300,420))
#points(QCL.DATA$hourinseconds[chamber2.2],QCL.DATA$spec.456a[chamber2.2],col=dream2,cex=.5,ylim=c(300,420))


# keeling plots alpha
plot(QCL.DATA$spec.446a.inv[chamber2.1],QCL.DATA$d15Na.corr2[chamber2.1],col="white",cex=.5,ylim=c(-20,50),xlim=c(0.0025,0.00295),xaxt="n",
     ylab=expression(paste(delta^{15},N[alpha])))

points(QCL.DATA$spec.446a.inv[chamber2.1],QCL.DATA$d15Na.corr2[chamber2.1],col=dream1,cex=.5)
abline(lm(QCL.DATA$d15Na.corr2[chamber2.1]~QCL.DATA$spec.446a.inv[chamber2.1]))

#points(QCL.DATA$spec.446a.inv[chamber2.2],QCL.DATA$d15Na.corr2[chamber2.2],col=dream4,cex=.5)
#abline(lm(QCL.DATA$d15Na.corr2[chamber2.2]~QCL.DATA$spec.446a.inv[chamber2.2]))


# concentrations beta
###########################################################################################################

plot(QCL.DATA$hourinseconds[chamber2.1],QCL.DATA$spec.446a[chamber2.1],col="white",cex=.5,ylim=c(300,420),xaxt="n",
     ylab=expression(""^{14}*"N"^{14}*"N"^{16}*"O & "^{15}*"N"^{14}*"N"^{16}*"O [ppb]"),xlab="Time [sec]")

points(QCL.DATA$hourinseconds[chamber2.1],QCL.DATA$spec.446a[chamber2.1],col=dream1,cex=.5,ylim=c(300,420))
points(QCL.DATA$hourinseconds[chamber2.1],QCL.DATA$spec.546a[chamber2.1],col=dream3,cex=.5,ylim=c(300,420))

#points(QCL.DATA$hourinseconds[chamber2.2],QCL.DATA$spec.446a[chamber2.2],col=dream1,cex=.5,ylim=c(300,420))
#points(QCL.DATA$hourinseconds[chamber2.2],QCL.DATA$spec.546a[chamber2.2],col=dream3,cex=.5,ylim=c(300,420))


# keeling plots beta
plot(QCL.DATA$spec.446a.inv[chamber2.1],QCL.DATA$d15Nb.corr2[chamber2.1],col="white",cex=.5,ylim=c(-40,15),xlim=c(0.0025,0.00295),xaxt="n",
     ylab=expression(paste(delta^{15},N[beta])))

points(QCL.DATA$spec.446a.inv[chamber2.1],QCL.DATA$d15Nb.corr2[chamber2.1],col=dream1,cex=.5)
abline(lm(QCL.DATA$d15Nb.corr2[chamber2.1]~QCL.DATA$spec.446a.inv[chamber2.1]))

#points(QCL.DATA$spec.446a.inv[chamber2.2],QCL.DATA$d15Nb.corr2[chamber2.2],col=dream4,cex=.5)
#abline(lm(QCL.DATA$d15Nb.corr2[chamber2.2]~QCL.DATA$spec.446a.inv[chamber2.2]))


# concentrations 18O
###########################################################################################################

plot(QCL.DATA$hourinseconds[chamber2.1],QCL.DATA$spec.446a[chamber2.1],col="white",cex=.5,ylim=c(300,420),
     ylab=expression(""^{14}*"N"^{14}*"N"^{16}*"O & "^{14}*"N"^{14}*"N"^{18}*"O [ppb]"),xlab="Time [sec]")

points(QCL.DATA$hourinseconds[chamber2.1],QCL.DATA$spec.446a[chamber2.1],col=dream1,cex=.5,ylim=c(300,420))
points(QCL.DATA$hourinseconds[chamber2.1],QCL.DATA$spec.448a[chamber2.1],col=dream5,cex=.5,ylim=c(300,420))

#points(QCL.DATA$hourinseconds[chamber2.2],QCL.DATA$spec.446a[chamber2.2],col=dream1,cex=.5,ylim=c(300,420))
#points(QCL.DATA$hourinseconds[chamber2.2],QCL.DATA$spec.448a[chamber2.2],col=dream5,cex=.5,ylim=c(300,420))


mtext("Time [secs]",1,line=3)

# keeling plots 18O
plot(QCL.DATA$spec.446a.inv[chamber2.1],QCL.DATA$d18O.corr2[chamber2.1],col="white",cex=.5,ylim=c(-100,100),xlim=c(0.0025,0.00295),
     ylab=expression(paste(delta^{18},O)))

points(QCL.DATA$spec.446a.inv[chamber2.1],QCL.DATA$d18O.corr2[chamber2.1],col=dream1,cex=.5)
abline(lm(QCL.DATA$d18O.corr2[chamber2.1]~QCL.DATA$spec.446a.inv[chamber2.1]))

#points(QCL.DATA$spec.446a.inv[chamber2.2],QCL.DATA$d18O.corr2[chamber2.2],col=dream4,cex=.5)
#abline(lm(QCL.DATA$d18O.corr2[chamber2.2]~QCL.DATA$spec.446a.inv[chamber2.2]))

# dev.off()    


#########################
### publication Fig 4 ###
#########################

# pdf("/Users/mbarthel/Desktop/QCL1/graphs/fluxvsdelta.pdf",encoding="MacRoman",height=8,width=7,family="Times",fonts="Times")   
par(mfrow=c(2,1),bg="white",pch=16,tcl=0.3,cex=1,las=0,srt=0,cex.lab=1*1.3,xpd=FALSE)
par(omi = c(0.7,0,0.09,0.09))
par(mai = c(1,1,0.1,0))

plot(ch.time,n2o.flux,
    ylab=expression("flux N"[2]*"O [nmol m"^{-2}*" s"^{-1}*"]"),xlab='Time',ylim = c(-0.5,2),
   col=dream1,cex=.5)
abline(h=0,lty=3)



plot(n2o.flux,coef(keeling.d15Nb)[,1],
    xlab=expression("flux N"[2]*"O [nmol m"^{-2}*" s"^{-1}*"]"),
    ylab=expression(paste(delta^{15},N[beta])),
    cex=.5,xlim=c(0,1.5),ylim=c(-10000,10000),col=dream1)

# dev.off()





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




# pdf("/Users/mbarthel/Desktop/QCL1/graphs/anchor-span.pdf",encoding="MacRoman",height=6,width=12,family="Times",fonts="Times")   
par(mfrow=c(1,1),bg="white",pch=16,tcl=0.3,cex=1,las=0,srt=0,cex.lab=1*1.3,xpd=FALSE)
par(omi = c(0.7,0,0.09,0.09))
par(mai = c(1,1,0.1,0))

plot(339:390,tapply(QCL.DATA$d15Nb.corr2[SPAN],QCL.DATA$DOY[SPAN],mean),ylim=c(30,50),col="white",xlab='Time [months]', ylab=expression(paste(delta^{15},N[beta])))

# d15Na
x.time     = 339:390
sd.d15Na   = tapply(QCL.DATA$d15Na.corr2[SPAN],QCL.DATA$DOY[SPAN],sd)
mean.d15Na = tapply(QCL.DATA$d15Na.corr2[SPAN],QCL.DATA$DOY[SPAN],mean)

points(x.time ,mean.d15Na,col=dream2)

    for (i in 1:length(x.time )){
        arrows(x.time[i], mean.d15Na[i] + sd.d15Na[i], x.time[i], mean.d15Na[i] - sd.d15Na[i],
        angle=90, code=3, length=0.03, col="grey50", lwd=0.4) # Fehlerbalken 
    }


# d15Nb
sd.d15Nb   = tapply(QCL.DATA$d15Nb.corr2[SPAN],QCL.DATA$DOY[SPAN],sd)
mean.d15Nb = tapply(QCL.DATA$d15Nb.corr2[SPAN],QCL.DATA$DOY[SPAN],mean)

points(x.time ,mean.d15Nb,col=dream3)

    for (i in 1:length(x.time )){
        arrows(x.time[i], mean.d15Nb[i] + sd.d15Nb[i], x.time[i], mean.d15Nb[i] - sd.d15Nb[i],
        angle=90, code=3, length=0.03, col="grey50", lwd=0.4) # Fehlerbalken 
    }


# d18O
sd.d18O   = tapply(QCL.DATA$d18O.corr2[SPAN],QCL.DATA$DOY[SPAN],sd)
mean.d18O = tapply(QCL.DATA$d18O.corr2[SPAN],QCL.DATA$DOY[SPAN],mean)

points(x.time ,mean.d18O,col=dream5)

    for (i in 1:length(x.time )){
        arrows(x.time[i], mean.d18O[i] + sd.d18O[i], x.time[i], mean.d18O[i] - sd.d18O[i],
        angle=90, code=3, length=0.03, col="grey50", lwd=0.4) # Fehlerbalken 
    }




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






