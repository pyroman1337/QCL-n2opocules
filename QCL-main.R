# *****************************************
# R SCRIPT FOR NITROGEN ISOTOPOMER ANALYSIS
# M. Barthel |  24 Sep 2015 
# version: 0.1 
# adopted for GitHub 10. Jan 2017 (R. Hüppi)

# libraries to include
library(nlme)   # package for function lmList

# GLOBAL PARAMs
#  ==============================================================================================================================================================================================================================================

  rm(list=ls())
  Sys.setenv(TZ='UTC')

data.folder <- "data/"  # set data folder for QCL rawdata (containing .stc and .str files)
  
  
kulergrey   <- rgb(64,55,47, maxColorValue =  255)
kulerorange <- rgb(217,86,28, maxColorValue = 255)

dream1      <- rgb(64,1,13, maxColorValue =   255)
dream2      <- rgb(217,4,82, maxColorValue =  255)
dream3      <- rgb(140,131,3, maxColorValue = 255)
dream4      <- rgb(89,43,2, maxColorValue =   255)
dream5      <- rgb(217,130,54, maxColorValue =255)

# FUNCTIONS TO REMOVE DATA AT START AND END AND THEN APPLY MEAN OR SD         

cutMean <-  function(x, start, end, ...) {
return(mean(x[start:(length(x)-end)], ...))
}
cutSD <- function(x, start, end, ...) {
return(sd(x[start:(length(x)-end)], ...))
}
cutSum <- function(x, start, end, ...) {
return(sum(x[start:(length(x)-end)], ...))
}


# TDL Wintel output *.stc and *.str time vector is given in seconds from 01 January 1904; that had been inherited from the IGOR programming language 
# species 1-5 are measured by Laser 1 (2188 cm-1)
# species 6-11 are measured by Laser 2 (2202 cm-1)        
  
# HITRAN ABUNDANCIES FOR N2O Species
# 446 9.903e-1
# 456 3.641e-3
# 546 3.641e-3
# 448 1.986e-3
# 447 3.693e-4

# all 9.903e-1+3.641e-3+3.641e-3+1.986e-3+3.693e-4
# (3.641e-3/9.903e-1)/0.0036782


# INTERNATIONAL STANDARDS
# AIR-N2 0.0036782 Isotope ratio of AIR-N2 from Werner & Brandt 2001
# V-SMOW 0.0020052 Isotope ratio of V-SMOW from Werner & Brandt 2001

 AIR.N2 <- 0.0036782 # Isotope ratio of AIR-N2 from Werner & Brandt 2001
 V.SMOW <- 0.0020052 # Isotope ratio of V-SMOW from Werner & Brandt 2001



 
 ETHZSAEHIGH1.d15Na <-  0.00 #plusminus 0.32
 ETHZSAEHIGH1.d15Nb <-  2.16 #plusminus 0.33
 ETHZSAEHIGH1.d18O  <- 38.98 #plusminus 0.33

 ETHZSAEHIGH7.d15Na <- 34.01 #plusminus 0.71
 ETHZSAEHIGH7.d15Nb <- 35.98 #plusminus 0.50
 ETHZSAEHIGH7.d18O  <- 37.99 #plusminus 0.35



 ETHZSAEHIGH1.R456  <- (ETHZSAEHIGH1.d15Na/1000 + 1) * AIR.N2 #conversion from d-value to ratio
 ETHZSAEHIGH1.R546  <- (ETHZSAEHIGH1.d15Nb/1000 + 1) * AIR.N2 #conversion from d-value to ratio
 ETHZSAEHIGH1.R448  <- (ETHZSAEHIGH1.d18O/ 1000 + 1) * V.SMOW #conversion from d-value to ratio






# location.graphs                              <-  "/Users/mbarthel/Desktop/QCL1/graphs/"  
# graphname.concentrations.ETHZSAEHIGH1.high   <-  "/Users/mbarthel/Desktop/QCL1/graphs/cal/2015-03-13-ETHZSAEHIGH1-mixhigh-concentrations.png"
# graphname.deltavalues.ETHZSAEHIGH1.high      <-  "/Users/mbarthel/Desktop/QCL1/graphs/cal/2015-03-13-ETHZSAEHIGH1-mixhigh-delta-values.png"
#  
# graphname.concentrations.ethzsaehigh <-  "/Users/mbarthel/Desktop/QCL1/graphs/cal/2015-03-13-ETHZSAEHIGH-concentrations.png"
# graphname.deltavalues.ethzsaehigh    <-  "/Users/mbarthel/Desktop/QCL1/graphs/cal/2015-03-13-ETHZSAEHIGH-delta-values.png"
# 
# graphname.concentrations.dilution <-  "/Users/mbarthel/Desktop/QCL1/graphs/cal/2015-03-13-DILUTION-concentrations.png"
# graphname.deltavalues.dilution    <-  "/Users/mbarthel/Desktop/QCL1/graphs/cal/2015-03-13-DILUTION-delta-values.png"


#  Reads in temperature data, obtained from Rscript 'ibutton.R'
#  ==============================================================================================================================================================================================================================================      

#TEMPERATURE <- read.table("/Users/mbarthel/Desktop/QCL1/data/GH-temperature-mod.txt",sep=",")
#TEMPERATURE$TIMESTAMP  <-  as.POSIXct(strptime(TEMPERATURE$TIMESTAMP,format="%Y-%m-%d %H:%M:%S",tz="UTC"))

#SUBSET        <- TEMPERATURE$TIMESTAMP > as.POSIXct("2015-09-15 00:00:00, tz=UTC") &  TEMPERATURE$TIMESTAMP < as.POSIXct("2016-01-25 23:59:59, tz=UTC") 
#TEMPERATURE   <- TEMPERATURE[SUBSET,]


#  Reads in all QCL *.str files
#  ==============================================================================================================================================================================================================================================      
#  file header: "time","spec.546a","spec.456a","spec.446a","spec.h2oa","spec.co2a","spec.448a","spec.n2oa","spec.446b","spec.co2b","spec.COa","spec.h2o", 
#  OCL$TIMESTAMP <- as.POSIXct(strptime(QCL$TIMESTAMP,format="%Y-%m-%d %H:%M:%S",))
#  the single headers line in the *.str files starts with time stream file was openend in four different formats: (1) local time as timestamp (2) local time in seconds since preceding midnight (3) local time in IGOR time, ie seconds since 1904 (4) UTC time


 # data.location.str              <-  "/DATA"                # data is placed on desktop
 # setwd(data.location.str) 
 filename.str                   <-  list.files(path = data.folder, pattern="*.str")   # relative Data path from package


        STR    <- data.frame(files=NULL, Month=NULL, Gain=NULL , Loss=NULL, sum=NULL, mean=NULL) # creates empty data frame container
        
        for(i in seq(along=filename.str)) {
                x <- read.table(paste0(data.folder,filename.str[i]), header=FALSE, skip=1,sep="",fill=TRUE,
                col.names = c("time","spec.546a","spec.456a","spec.446a","spec.h2oa","spec.co2a","spec.448a","spec.n2o","spec.446b","spec.co2b","spec.CO","spec.h2ob"))
                STR <- rbind(STR, x) # appends data from new files to data frame 'STR'
        }
        
       STR$TIMESTAMP     <- ISOdatetime(1904,1,1,0,0,0,tz="UTC") + STR[,"time"]
       STR$TIMESTAMP     <- as.POSIXct(strptime(STR$TIMESTAMP,format="%Y-%m-%d %H:%M:%S",tz="UTC"))           

#  Reads in all QCL *.stc files and creates a subset data frame
#  ==============================================================================================================================================================================================================================================      

 # data.location.stc              <-  "/data"                # data is placed on desktop
 # setwd(data.location.stc) 
 filename.stc                   <-  list.files(path = data.folder ,pattern="*.stc") # relative Data path from package
    
  
  
#  Reads in all QCL *.stc files
#  ==============================================================================================================================================================================================================================================        
        QCL.stc    <- data.frame(files=NULL, Month=NULL, Gain=NULL , Loss=NULL, sum=NULL, mean=NULL) # creates empty data frame container
        
        for(i in seq(along=filename.stc)) {
                x <- read.table(paste0(data.folder,filename.stc[i]), header=FALSE, skip=2,sep=",",fill=TRUE,
                col.names = c("time","rangeF1L1","rangeF1L2","rangeF2L1","rangeF2L2","Pcell","Tcell","Pref","Tref","AD8","AD9","AD10","AD11","AD12","AD13","AD14","AD15",
                              "statusW","VICI_W","USBByte","NI6024Byte","SPEFile","Tlaser1","Vlaser1","LWlaser1","Tlaser2","Vlaser2","LWlaser2","dT1","dT2","X1","pos1","X2","pos2"))
                QCL.stc <- rbind(QCL.stc, x) # appends data from new files to data frame 'QCL.stc'
        }
        
        # creates data subset from QCL and removes QCL 
        STC <- QCL.stc[,c("time","rangeF1L1","rangeF1L2","rangeF2L1","rangeF2L2","Pcell","Tcell","Pref","Tref","AD8","AD9","AD10","AD11","AD12","AD13","AD14","AD15",
                              "statusW","VICI_W","Tlaser1","Vlaser1","LWlaser1","Tlaser2","Vlaser2","LWlaser2","dT1","dT2","X1","pos1","X2","pos2")]

        STC$TIMESTAMP     <- ISOdatetime(1904,1,1,0,0,0,tz="UTC") + STC[,"time"] 
        STC$TIMESTAMP     <- as.POSIXct(strptime(STC$TIMESTAMP,format="%Y-%m-%d %H:%M:%S",tz="UTC"))        
        rm(QCL.stc)
        







        # MERGE DATA FRAMES


        #FOO  <- merge(STC,TEMPERATURE,by ="TIMESTAMP",all=TRUE, sort=TRUE, incomparables=TRUE)
        QCL  <- merge(STC,STR        ,by ="TIMESTAMP",all=TRUE, sort=TRUE, incomparables=TRUE)
        
        rm(TEMPERATURE)
        rm(STR)
        rm(STC)
        rm(FOO)
       
        # Aggregate the Data in flexible integration time (ideally suggested by the lowest Allan variance)  
        str(QCL)
        library(data.table)
        QCL.dt <- as.data.table(QCL)  
        QCL.dt[ ]
        
        # IMPLEMENT TIME VECTORS
    
        QCL$DOY            <-  as.numeric(format(QCL$TIMESTAMP,format="%j"))
        FOO1               <-  QCL$TIMESTAMP > as.POSIXct("2016-01-01 00:59:59, tz=UTC")
        QCL$DOY[FOO1]      <-  QCL$DOY[FOO1]  + 365           # create continuous DOY value       

        
        QCL$doy.frac       <-  as.numeric(format(QCL$TIMESTAMP,format="%j")) + as.numeric(format(QCL$TIMESTAMP,format="%H"))/24 + as.numeric(format(QCL$TIMESTAMP,format="%M"))/24/60 + as.numeric(format(QCL$TIMESTAMP,format="%S"))/24/60/60         
        QCL$weekday        <-  as.numeric(format(QCL$TIMESTAMP,format="%u")) 
        QCL$week           <-  as.numeric(format(QCL$TIMESTAMP,format="%W")) # week of year, UK convention (week starts Monday)        
        QCL$hour           <-  as.numeric(format(QCL$TIMESTAMP,format="%H"))
        QCL$HOUR           <-  format(QCL$TIMESTAMP,format="%H")
        QCL$hourinseconds  <-  as.numeric(format(QCL$TIMESTAMP,format= "%M"))*60 + as.numeric(format(QCL$TIMESTAMP,format="%S"))
        QCL$hourinminutes  <-  as.numeric(format(QCL$TIMESTAMP,format= "%M")) 
        
        QCL$timecode       <-  as.numeric(paste(QCL$DOY,QCL$hour,sep = ""))

        QCL$ampm           <-  format(QCL$TIMESTAMP,format="%p") # AM/PM indicator
        QCL$ampmnumeric    <-  numeric(length(QCL[,2]))
        QCL$ampmnumeric[QCL$ampm == "PM"] =  1
        QCL$calcode        <-  as.numeric(paste(QCL$DOY,QCL$ampmnumeric,sep=""))





QCL$VICI  <- QCL$VICI_W + 1 # TDL Wintel is writing VICI multivalve positions from 0 to 15, however, in software used valve numbers are from # 1 to 16        
QCL$N     <- !is.na(QCL$spec.446a)  # CREATE LOGICAL VECTOR TO KNOW n of array after aggregation          

# calculate d-values based on QCL mixing ratios and using HITRAN abundancies
QCL$d15Na <- ((((0.003641 * QCL$spec.456a) / (0.9903 * QCL$spec.446a)) / AIR.N2) - 1) * 1000  # referenced against AIR-N2
QCL$d15Nb <- ((((0.003641 * QCL$spec.546a) / (0.9903 * QCL$spec.446a)) / AIR.N2) - 1) * 1000  # referenced against AIR-N2
QCL$d18O  <- ((((0.001986 * QCL$spec.448a) / (0.9903 * QCL$spec.446a)) / V.SMOW) - 1) * 1000  # referenced against V-SMOW
QCL$bulk  <- (QCL$d15Na + QCL$d15Nb) / 2
QCL$SP    <-  QCL$d15Na - QCL$d15Nb


# calculate Ratios 
QCL$R456  <- (0.003641 * QCL$spec.456a)/(0.9903 * QCL$spec.446a)
QCL$R546  <- (0.003641 * QCL$spec.546a)/(0.9903 * QCL$spec.446a)
QCL$R448  <- (0.001986 * QCL$spec.448a)/(0.9903 * QCL$spec.446a)




# CALIBRATION CALCULATION AND EQUATIONs
# =======================================================================================================================================================================
NA.RM.1  <- !is.na(QCL$calcode) & !is.na(QCL$spec.446a)   & !is.na(QCL$d15Na)    & !is.na(QCL$d15Nb) & !is.na(QCL$d18O)
  

# OFFSET CORRECTION USING ETHZSAEHIGH1 at VICI port #3
OFFSET <- NA.RM.1 & QCL$VICI == 3 & QCL$hourinseconds > 250 &  QCL$hourinseconds < 310 & QCL$hour == 3 |  # last minute of cal cycle at VICI#3 before dilution starts, 
         NA.RM.1 & QCL$VICI == 3 & QCL$hourinseconds > 250 &  QCL$hourinseconds < 310 & QCL$hour == 15
# average of ETHZSAEHIGH1 measurement for each hour

QCLETHZSAEHIGH1.446    <- tapply(QCL$spec.446a[OFFSET],QCL$calcode[OFFSET],mean,na.rm=TRUE) 
QCLETHZSAEHIGH1.d15Na  <- tapply(QCL$d15Na[OFFSET],    QCL$calcode[OFFSET],mean,na.rm=TRUE) 
QCLETHZSAEHIGH1.d15Nb  <- tapply(QCL$d15Nb[OFFSET],    QCL$calcode[OFFSET],mean,na.rm=TRUE) 
QCLETHZSAEHIGH1.d18O   <- tapply(QCL$d18O[OFFSET],     QCL$calcode[OFFSET],mean,na.rm=TRUE) 

# standard deviation of ETHZSAEHIGH1 measurement for each hour
QCLETHZSAEHIGH1.d15Na.error  <- sd(QCL$d15Na[OFFSET])
QCLETHZSAEHIGH1.d15Nb.error  <- sd(QCL$d15Nb[OFFSET])
QCLETHZSAEHIGH1.d18O.error   <- sd(QCL$d18O[OFFSET])

# difference between target and measured value
diff.d15Na <- ((ETHZSAEHIGH1.d15Na - QCLETHZSAEHIGH1.d15Na)/(QCLETHZSAEHIGH1.d15Na+1000))*1000
diff.d15Nb <- ((ETHZSAEHIGH1.d15Nb - QCLETHZSAEHIGH1.d15Nb)/(QCLETHZSAEHIGH1.d15Nb+1000))*1000
diff.d18O  <- ((ETHZSAEHIGH1.d18O  - QCLETHZSAEHIGH1.d18O)/(QCLETHZSAEHIGH1.d18O+1000))*1000





# DILUTION CORRECTION USING ETHZSAEHIGH1 at VICI port #3

 # PERIODS USEABLE FOR CORRECTION
 # VICI #3 ETHZSAEHIGH-1

     DILUTION <- NA.RM.1 & QCL$VICI == 3 & QCL$hourinseconds > 310 & QCL$hour == 3 | 
                NA.RM.1 & QCL$VICI == 3 & QCL$hourinseconds > 310 & QCL$hour == 15


    # library(nlme)   # package for function lmList
    QCL.DILUTION.SUBSET          <-      subset(QCL[DILUTION,], select=c("calcode","spec.446a","d15Na","d15Nb","d18O"))
    # DILUTION fit parameters  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     

    fit.d15Na            <-      lmList(d15Na ~ spec.446a | calcode,data = QCL.DILUTION.SUBSET)
    fit.d15Nb            <-      lmList(d15Nb ~ spec.446a | calcode,data = QCL.DILUTION.SUBSET)
    fit.d18O             <-      lmList(d18O ~ spec.446a  | calcode,data = QCL.DILUTION.SUBSET)
    
    # d15Na
    R2.d15Na             <-      summary(fit.d15Na)$r.squared
    df.d15Na             <-      summary(fit.d15Na)$df[,2]
    offset.d15Na         <-      coef(fit.d15Na)[,1]
    slope.d15Na          <-      coef(fit.d15Na)[,2]

    # d15Nb
    R2.d15Nb             <-      summary(fit.d15Nb)$r.squared
    df.d15Nb             <-      summary(fit.d15Nb)$df[,2]
    offset.d15Nb         <-      coef(fit.d15Nb)[,1]
    slope.d15Nb          <-      coef(fit.d15Nb)[,2]


    # d18O
    R2.d18O             <-      summary(fit.d18O)$r.squared
    df.d18O             <-      summary(fit.d18O)$df[,2]
    offset.d18O         <-      coef(fit.d18O)[,1]
    slope.d18O          <-      coef(fit.d18O)[,2]


    calcode                <- tapply(QCL$calcode[DILUTION],QCL$calcode[DILUTION],mean,na.rm=TRUE) 

    CALIBRATION <-  data.frame(calcode,QCLETHZSAEHIGH1.446,diff.d15Na,R2.d15Na,df.d15Na,offset.d15Na,slope.d15Na,
                              diff.d15Nb,R2.d15Nb,df.d15Nb,offset.d15Nb,slope.d15Nb,
                              diff.d18O,R2.d18O,df.d18O,offset.d18O,slope.d18O)



DATA    <-    merge(QCL,CALIBRATION     ,by="calcode",all=TRUE, sort=FALSE, incomparables=NA)
DATA    <-    DATA[order(DATA$TIMESTAMP),] # reorder data.frame as merge in previous line fails to sort correctly 



# APPLYING CORRECTIONS

# OFFSET CORRECTION
DATA$d15Na.corr1 <- DATA$d15Na + DATA$diff.d15Na + (DATA$d15Na * DATA$diff.d15Na * 0.001)
DATA$d15Nb.corr1 <- DATA$d15Nb + DATA$diff.d15Nb + (DATA$d15Nb * DATA$diff.d15Nb * 0.001)
DATA$d18O.corr1  <- DATA$d18O  + DATA$diff.d18O  + (DATA$d18O  * DATA$diff.d18O  * 0.001)


# DILUTION CORRECTION
# based on the below assessed linear function with *.corr2 from the dilution calibration the data is corrected for with y2 = m(x2-x1) + y1 from m = (y2 - y1) / (x2 - x1)
 QCLETHZSAEHIGH1.446 # equals true concentration of tank for which values are corrected for

#dilution correction    = Anstieg von corr2 * concentration tank       - momentane conc   +  gemessene delta.corr1
 DATA$d15Na.corr2       <- DATA$slope.d15Na * (DATA$QCLETHZSAEHIGH1.446 - DATA$spec.446a)  +  DATA$d15Na.corr1
 DATA$d15Nb.corr2       <- DATA$slope.d15Nb * (DATA$QCLETHZSAEHIGH1.446 - DATA$spec.446a)  +  DATA$d15Nb.corr1
 DATA$d18O.corr2        <- DATA$slope.d18O  * (DATA$QCLETHZSAEHIGH1.446 - DATA$spec.446a)  +  DATA$d18O.corr1
 
 # validation plot 
 # plot(QCL$TIMESTAMP[DILUTION],QCL$spec.446a[DILUTION])
 plot(QCL$spec.446a[DILUTION],QCL$d15Nb[DILUTION],xlim=c(0,3000),ylim=c(-70,10))
 abline(lm(QCL$d15Nb[DILUTION]~QCL$spec.446a[DILUTION]))
 points(DATA$spec.446a[DILUTION],DATA$d15Nb.corr2[DILUTION],col="red")
 abline(0,0)
 
 plot(calcode,offset.d15Nb, ylim = c(-32,-20))  # check the variability of the intercept
 plot(calcode,slope.d15Nb, ylim = c(0,0.02))    # check the variability of the slope
 
# ANCHOR (ETHZSAEHIGH-1) & SPAN (ETHZSAEHIGH-7) TARGET CHECK
SPAN   <- DATA$VICI == 4 & DATA$hourinseconds > 480 & DATA$TIMESTAMP > as.POSIXct("2015-10-30 00:00:00, tz=UTC")  # 2015-12-05 is the actual starting date
ANCHOR <- NA.RM.1 & DATA$VICI == 3 & DATA$hourinseconds > 250 &  DATA$hourinseconds < 310 & DATA$hour == 3   # & DATA$TIMESTAMP > as.POSIXct("2015-12-05 00:00:00, tz=UTC") |  # last minute of cal cycle at VICI#3 before dilution starts, 
         NA.RM.1 & DATA$VICI == 3 & DATA$hourinseconds > 250 &  DATA$hourinseconds < 310 & DATA$hour == 15  # & DATA$TIMESTAMP > as.POSIXct("2015-12-05 00:00:00, tz=UTC")


pdf("graphs/anchor.pdf",encoding="MacRoman",height=6,width=12,family="Times",fonts="Times")   
par(mfrow=c(1,1),bg="white",pch=16,tcl=0.3,cex=1,las=0,srt=0,cex.lab=1*1.3,xpd=FALSE)
par(omi = c(0.7,0,0.09,0.09))
par(mai = c(1,1,0.1,0))

plot(DATA$TIMESTAMP[ANCHOR],DATA$d15Nb.corr2[ANCHOR],ylim=c(-20,60),col="white",xlab='Time [months]', ylab=expression(paste(delta^{15},N[beta])),
  xlim=c(as.POSIXct("2015-10-30 00:00:00, tz=UTC"),as.POSIXct("2016-01-28 00:00:00, tz=UTC")))  # 2015-12-05 is the actual starting date

points(DATA$TIMESTAMP[ANCHOR],DATA$d15Nb.corr2[ANCHOR],col=dream1,cex=0.3,pch=16)
abline(h=ETHZSAEHIGH1.d15Nb,col=dream1)

points(DATA$TIMESTAMP[SPAN],DATA$d15Nb.corr2[SPAN],col=dream2,cex=0.3,pch=16)
abline(h=ETHZSAEHIGH7.d15Nb,col=dream2)
dev.off()


# CHAMBER VOLUME CALCULATIONS

# radius chamber cylinder 0.235 m
# volume remaining cylinder airspace in liter (0.235m^2*pi*rimtosoildistance.in.meter) * 1000 
# volume chamber lid in liter (0.5m*0.5m*0.15m) * 1000
# volume extensions in liter (0.5m*0.5m*0.5m) * 1000 

# started experiment with chamber lids and one extension thus total volume = ((0.235^2*pi*0.11) + (0.5*0.5*0.15) + (0.5*0.5*0.5)) * 1000
ch.volumes <- numeric(length(DATA[,2]))

ch.volumes[DATA$VICI == 5]  <- 162.5 + 26.5
ch.volumes[DATA$VICI == 6]  <- 162.5 + 28.7
ch.volumes[DATA$VICI == 7]  <- 162.5 + 32.6
ch.volumes[DATA$VICI == 8]  <- 162.5 + 33.2
ch.volumes[DATA$VICI == 9]  <- 162.5 + 38.0
ch.volumes[DATA$VICI == 10] <- 162.5 + 33.2
ch.volumes[DATA$VICI == 11] <- 162.5 + 31.8
ch.volumes[DATA$VICI == 12] <- 162.5 + 38.4
ch.volumes[DATA$VICI == 13] <- 162.5 + 35.9
ch.volumes[DATA$VICI == 14] <- 162.5 + 35.0
ch.volumes[DATA$VICI == 15] <- 162.5 + 33.0
ch.volumes[DATA$VICI == 16] <- 162.5 + 35.1


ROUND1  <- as.POSIXct("2015-11-06 00:00:00, tz=UTC")
ROUND2  <- as.POSIXct("2015-11-13 00:00:00, tz=UTC") 
ROUND3  <- as.POSIXct("2015-12-11 00:00:00, tz=UTC") 
NA.RM.2  <-  !is.na(DATA$TIMESTAMP) & !is.na(DATA$VICI) 


ch.volumes[DATA$TIMESTAMP > ROUND1 & DATA$VICI == 8  & NA.RM.2] <- ch.volumes[DATA$TIMESTAMP > ROUND1 & DATA$VICI == 8 & NA.RM.2]  + (0.5*0.5*0.5)*1000 # extension on chamber 4 on 6 Nov 2015
ch.volumes[DATA$TIMESTAMP > ROUND1 & DATA$VICI == 16 & NA.RM.2] <- ch.volumes[DATA$TIMESTAMP > ROUND1 & DATA$VICI == 16 & NA.RM.2] + (0.5*0.5*0.5)*1000 # extension on chamber 12 on 6 Nov 2015

ch.volumes[DATA$TIMESTAMP > ROUND2 & DATA$VICI == 11 & NA.RM.2] <- ch.volumes[DATA$TIMESTAMP > ROUND2 & DATA$VICI == 11 & NA.RM.2] + (0.5*0.5*0.5)*1000 # extension on chamber 7 on 13 Nov 2015
ch.volumes[DATA$TIMESTAMP > ROUND2 & DATA$VICI == 13 & NA.RM.2] <- ch.volumes[DATA$TIMESTAMP > ROUND2 & DATA$VICI == 13 & NA.RM.2] + (0.5*0.5*0.5)*1000 # extension on chamber 9 on 13 Nov 2015

ch.volumes[DATA$TIMESTAMP > ROUND3 & DATA$VICI == 6 & NA.RM.2]  <- ch.volumes[DATA$TIMESTAMP > ROUND3 & DATA$VICI == 6 & NA.RM.2]  + (0.5*0.5*0.5)*1000 # extension on chamber 2  on 11 Dec 2015 
ch.volumes[DATA$TIMESTAMP > ROUND3 & DATA$VICI == 9 & NA.RM.2]  <- ch.volumes[DATA$TIMESTAMP > ROUND3 & DATA$VICI == 9 & NA.RM.2]  + (0.5*0.5*0.5)*1000 # extension on chamber 5 on 11 Dec 2015 

DATA$ch.volumes <- ch.volumes


# ===============================================================================================================================================================================================================================================
# FLUX CALCULATIONS based on Betsys Excel spreadsheet, but linear regression only
#NA.RM.3  <-  !is.na(DATA$spec.446a) & !is.na(DATA$VICI) 

# DATA$N2O.nmol[NA.RM.3 & DATA$VICI == 5]       <- 1.00035 * DATA$spec.446a[NA.RM.3 & DATA$VICI == 5]  / (0.08206 * (273.15 + DATA$temp_C1_mod[NA.RM.3 & DATA$VICI  == 5]))  # nmol/L
# DATA$N2O.nmol[NA.RM.3 & DATA$VICI == 6]       <- 1.00035 * DATA$spec.446a[NA.RM.3 & DATA$VICI == 6]  / (0.08206 * (273.15 + DATA$temp_C2_mod[NA.RM.3 & DATA$VICI  == 6]))  # nmol/L
# DATA$N2O.nmol[NA.RM.3 & DATA$VICI == 7]       <- 1.00035 * DATA$spec.446a[NA.RM.3 & DATA$VICI == 7]  / (0.08206 * (273.15 + DATA$temp_C3_mod[NA.RM.3 & DATA$VICI  == 7]))  # nmol/L
# DATA$N2O.nmol[NA.RM.3 & DATA$VICI == 8]       <- 1.00035 * DATA$spec.446a[NA.RM.3 & DATA$VICI == 8]  / (0.08206 * (273.15 + DATA$temp_C4_mod[NA.RM.3 & DATA$VICI  == 8]))  # nmol/L
# DATA$N2O.nmol[NA.RM.3 & DATA$VICI == 9]       <- 1.00035 * DATA$spec.446a[NA.RM.3 & DATA$VICI == 9]  / (0.08206 * (273.15 + DATA$temp_C5_mod[NA.RM.3 & DATA$VICI  == 9]))  # nmol/L
# DATA$N2O.nmol[NA.RM.3 & DATA$VICI == 10]      <- 1.00035 * DATA$spec.446a[NA.RM.3 & DATA$VICI == 10] / (0.08206 * (273.15 + DATA$temp_C6_mod[NA.RM.3 & DATA$VICI  == 10])) # nmol/L
# DATA$N2O.nmol[NA.RM.3 & DATA$VICI == 11]      <- 1.00035 * DATA$spec.446a[NA.RM.3 & DATA$VICI == 11] / (0.08206 * (273.15 + DATA$temp_C7_mod[NA.RM.3 & DATA$VICI  == 11])) # nmol/L
# DATA$N2O.nmol[NA.RM.3 & DATA$VICI == 12]      <- 1.00035 * DATA$spec.446a[NA.RM.3 & DATA$VICI == 12] / (0.08206 * (273.15 + DATA$temp_C8_mod[NA.RM.3 & DATA$VICI  == 12])) # nmol/L
# DATA$N2O.nmol[NA.RM.3 & DATA$VICI == 13]      <- 1.00035 * DATA$spec.446a[NA.RM.3 & DATA$VICI == 13] / (0.08206 * (273.15 + DATA$temp_C9_mod[NA.RM.3 & DATA$VICI  == 13])) # nmol/L
# DATA$N2O.nmol[NA.RM.3 & DATA$VICI == 14]      <- 1.00035 * DATA$spec.446a[NA.RM.3 & DATA$VICI == 14] / (0.08206 * (273.15 + DATA$temp_C10_mod[NA.RM.3 & DATA$VICI == 14])) # nmol/L
# DATA$N2O.nmol[NA.RM.3 & DATA$VICI == 15]      <- 1.00035 * DATA$spec.446a[NA.RM.3 & DATA$VICI == 15] / (0.08206 * (273.15 + DATA$temp_C11_mod[NA.RM.3 & DATA$VICI == 15])) # nmol/L
# DATA$N2O.nmol[NA.RM.3 & DATA$VICI == 16]      <- 1.00035 * DATA$spec.446a[NA.RM.3 & DATA$VICI == 16] / (0.08206 * (273.15 + DATA$temp_C12_mod[NA.RM.3 & DATA$VICI == 16])) # nmol/L
DATA$N2O.nmol[NA.RM.3 ]      <- 1.00035 * DATA$spec.446a[NA.RM.3] / (0.08206 * (273.15 + 20)) # nmol/L  - no variable temperature data



# DATA$CO2.nmol[NA.RM.3 & DATA$VICI == 5]       <- 1.00035 * DATA$spec.co2b[NA.RM.3 & DATA$VICI == 5]  / (0.08206 * (273.15 + DATA$temp_C1_mod[NA.RM.3 & DATA$VICI  == 5]))  # nmol/L
# DATA$CO2.nmol[NA.RM.3 & DATA$VICI == 6]       <- 1.00035 * DATA$spec.co2b[NA.RM.3 & DATA$VICI == 6]  / (0.08206 * (273.15 + DATA$temp_C2_mod[NA.RM.3 & DATA$VICI  == 6]))  # nmol/L
# DATA$CO2.nmol[NA.RM.3 & DATA$VICI == 7]       <- 1.00035 * DATA$spec.co2b[NA.RM.3 & DATA$VICI == 7]  / (0.08206 * (273.15 + DATA$temp_C3_mod[NA.RM.3 & DATA$VICI  == 7]))  # nmol/L
# DATA$CO2.nmol[NA.RM.3 & DATA$VICI == 8]       <- 1.00035 * DATA$spec.co2b[NA.RM.3 & DATA$VICI == 8]  / (0.08206 * (273.15 + DATA$temp_C4_mod[NA.RM.3 & DATA$VICI  == 8]))  # nmol/L
# DATA$CO2.nmol[NA.RM.3 & DATA$VICI == 9]       <- 1.00035 * DATA$spec.co2b[NA.RM.3 & DATA$VICI == 9]  / (0.08206 * (273.15 + DATA$temp_C5_mod[NA.RM.3 & DATA$VICI  == 9]))  # nmol/L
# DATA$CO2.nmol[NA.RM.3 & DATA$VICI == 10]      <- 1.00035 * DATA$spec.co2b[NA.RM.3 & DATA$VICI == 10] / (0.08206 * (273.15 + DATA$temp_C6_mod[NA.RM.3 & DATA$VICI  == 10])) # nmol/L
# DATA$CO2.nmol[NA.RM.3 & DATA$VICI == 11]      <- 1.00035 * DATA$spec.co2b[NA.RM.3 & DATA$VICI == 11] / (0.08206 * (273.15 + DATA$temp_C7_mod[NA.RM.3 & DATA$VICI  == 11])) # nmol/L
# DATA$CO2.nmol[NA.RM.3 & DATA$VICI == 12]      <- 1.00035 * DATA$spec.co2b[NA.RM.3 & DATA$VICI == 12] / (0.08206 * (273.15 + DATA$temp_C8_mod[NA.RM.3 & DATA$VICI  == 12])) # nmol/L
# DATA$CO2.nmol[NA.RM.3 & DATA$VICI == 13]      <- 1.00035 * DATA$spec.co2b[NA.RM.3 & DATA$VICI == 13] / (0.08206 * (273.15 + DATA$temp_C9_mod[NA.RM.3 & DATA$VICI  == 13])) # nmol/L
# DATA$CO2.nmol[NA.RM.3 & DATA$VICI == 14]      <- 1.00035 * DATA$spec.co2b[NA.RM.3 & DATA$VICI == 14] / (0.08206 * (273.15 + DATA$temp_C10_mod[NA.RM.3 & DATA$VICI == 14])) # nmol/L
# DATA$CO2.nmol[NA.RM.3 & DATA$VICI == 15]      <- 1.00035 * DATA$spec.co2b[NA.RM.3 & DATA$VICI == 15] / (0.08206 * (273.15 + DATA$temp_C11_mod[NA.RM.3 & DATA$VICI == 15])) # nmol/L
# DATA$CO2.nmol[NA.RM.3 & DATA$VICI == 16]      <- 1.00035 * DATA$spec.co2b[NA.RM.3 & DATA$VICI == 16] / (0.08206 * (273.15 + DATA$temp_C12_mod[NA.RM.3 & DATA$VICI == 16])) # nmol/L
DATA$CO2.nmol[NA.RM.3 ]      <- 1.00035 * DATA$spec.co2b[NA.RM.3] / (0.08206 * (273.15 + 20)) # nmol/L   - no variable temperature data


DATA$VICI.time     <- paste(DATA$DOY,DATA$HOUR,DATA$VICI,sep = "") # use "HOUR" instead of "hour" as former is including 0 when counting from 1-9 
DATA$LABEL         <- paste(DATA$weekday ,DATA$ampm,sep = "_")

DATA$spec.446a.inv <- 1/ DATA$spec.446a




NA.RM.4 <- !is.na(DATA$hourinseconds) & !is.na(DATA$VICI.time)   & !is.na(DATA$N2O.nmol)    & !is.na(DATA$CO2.nmol) & 
           !is.na(DATA$spec.446a.inv) & !is.na(DATA$d15Na.corr2) & !is.na(DATA$d15Nb.corr2) & !is.na(DATA$d18O.corr2) & !is.na(DATA$LABEL)
NA.RM.4 <- as.vector(NA.RM.4)



# In the following respective Keeling plots are done for all chambers
KEELING <- # NA.RM.4 & 
           DATA$hourinseconds > 1200 & DATA$hourinseconds < 3540 & DATA$VICI > 4 & # start keeling plots 20min into the hour for all VICI ports > 4 (all chambers)
           DATA$LABEL != "4_PM" & DATA$LABEL != "5_AM"         # exclude 13C label periods Friday morning

KEELING <- as.vector(KEELING)

# library(nlme)   # package for function lmList  moved to top

    DATA.KEELING.SUBSET      <-      na.omit(subset(DATA[KEELING,], select=c("hourinseconds","VICI.time","N2O.nmol","CO2.nmol","spec.446a.inv","d15Na.corr2","d15Nb.corr2","d18O.corr2")))
    

    # DILUTION fit parameters  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
    volume                   <-      tapply(DATA$ch.volumes[KEELING],DATA$VICI.time[KEELING],mean)
    flux.time                <-      tapply(DATA$TIMESTAMP[KEELING],DATA$VICI.time[KEELING],mean)
    ch.time                  <-      ISOdatetime(1970,1,1,0,0,0) + flux.time    
    port                     <-      tapply(DATA$VICI[KEELING],DATA$VICI.time[KEELING],mean)
    
    
    # Keeling plot for d15Na
    keeling.d15Na            <-      lmList(d15Na.corr2 ~ spec.446a.inv | VICI.time,data = DATA.KEELING.SUBSET)    
    keeling.R2.d15Na         <-      summary(keeling.d15Na)$r.squared
    keeling.df.d15Na         <-      summary(keeling.d15Na)$df[,2]
    keeling.offset.d15Na     <-      coef(keeling.d15Na)[,1]
    keeling.slope.d15Na      <-      coef(keeling.d15Na)[,2]
    
    plot(ch.time,keeling.R2.d15Na)
    
    # Keeling plot for d15Nb
    keeling.d15Nb            <-      lmList(d15Nb.corr2 ~ spec.446a.inv | VICI.time,data = DATA.KEELING.SUBSET)    
    keeling.R2.d15Nb         <-      summary(keeling.d15Nb)$r.squared
    keeling.df.d15Nb         <-      summary(keeling.d15Nb)$df[,2]
    keeling.offset.d15Nb     <-      coef(keeling.d15Nb)[,1]
    keeling.slope.d15Nb      <-      coef(keeling.d15Nb)[,2]  

    # Keeling plot for d18O
    keeling.d18O            <-      lmList(d18O.corr2 ~ spec.446a.inv | VICI.time,data = DATA.KEELING.SUBSET)    
    keeling.R2.d18O         <-      summary(keeling.d18O)$r.squared
    keeling.df.d18O         <-      summary(keeling.d18O)$df[,2]
    keeling.offset.d18O     <-      coef(keeling.d18O)[,1]
    keeling.slope.d18O      <-      coef(keeling.d18O)[,2]    
    

    
    
# In the following all fluxes are computed    
FLUX.N2O    <- NA.RM.4 & 
              DATA$hourinseconds > 2000 & DATA$hourinseconds < 3540 & DATA$VICI > 4 & # start keeling plots 20min into the hour for all VICI ports > 4 (all chambers)
              DATA$LABEL != "4_PM" & DATA$LABEL != "5_AM"         # exclude 13C label periods Friday morning



FLUX.CO2    <- NA.RM.4 & 
              DATA$hourinseconds >  900 & DATA$hourinseconds < 2000 & DATA$VICI > 4 & # start keeling plots 20min into the hour for all VICI ports > 4 (all chambers)
              DATA$LABEL != "4_PM" & DATA$LABEL != "5_AM"         # exclude 13C label periods Friday morning


FLUX.N2O    <- as.vector(FLUX.N2O)
FLUX.CO2    <- as.vector(FLUX.CO2)

# library(nlme)   # package for function lmList  -> moved to top

    DATA.FLUX.SUBSET.N2O      <-      na.omit(subset(DATA[FLUX.N2O,], select=c("hourinseconds","VICI.time","N2O.nmol","CO2.nmol","spec.446a.inv","d15Na.corr2","d15Nb.corr2","d18O.corr2","spec.446a")))
    DATA.FLUX.SUBSET.CO2      <-      na.omit(subset(DATA[FLUX.CO2,], select=c("hourinseconds","VICI.time","N2O.nmol","CO2.nmol","spec.446a.inv","d15Na.corr2","d15Nb.corr2","d18O.corr2","spec.446a")))
        
    
    # CO2 Flux calculation
    co2.flux.regression      <-      lmList(CO2.nmol ~ hourinseconds | VICI.time,data = DATA.FLUX.SUBSET.CO2)
    co2.R2                   <-      summary(co2.flux.regression)$r.squared
    co2.flux                 <-      coef(co2.flux.regression)[,2]  * (volume/(0.235^2*pi))

    # N2O Flux calculation          
    n2o.flux.regression      <-      lmList(N2O.nmol ~ hourinseconds | VICI.time,data = DATA.FLUX.SUBSET.N2O)
    n2o.R2                   <-      summary(n2o.flux.regression)$r.squared    
    n2o.flux                 <-      coef(n2o.flux.regression)[,2]  * (volume/(0.235^2*pi))

    
    
   plot(n2o.flux) # check output
   plot(co2.flux) # check output
   
    
    
    