# *****************************************
# R SCRIPT FOR NITROGEN ISOTOPOMER ANALYSIS
# M. Barthel |  24 Sep 2015 
# version: 0.2
# adopted for GitHub 10. Jan 2017 (R. Hüppi)

# libraries to include
library(nlme)   # package for function lmList
library(data.table)  # data table functions (largely compatible with R-base dataframes)
library(lubridate)   # round date vector in order to average the data vector
library(corrplot)
library(plotly)

# GLOBAL PARAMs
#  ==============================================================================================================================================================================================================================================

  rm(list=ls())
  Sys.setenv(TZ='UTC')
  lct <- Sys.getlocale("LC_TIME"); Sys.setlocale("LC_TIME", "C")      # sets time settings to english (no confusion with local settings)
  
data.folder <- "data/2015-10/"  # set data folder for QCL rawdata (containing .stc and .str files)
# data.folder <- "/home/hueppir/DATA-QCL/"  # set data folder for QCL rawdata (containing .stc and .str files)

# kulergrey   <- rgb(64,55,47, maxColorValue =  255) # remaining example of kuler colours

# TDL Wintel output *.stc and *.str time vector is given in seconds from 01 January 1904; that had been inherited from the IGOR programming language 
# species 1-5 are measured by Laser 1 (2188 cm-1)
# species 6-11 are measured by Laser 2 (2202 cm-1)        
  

# INTERNATIONAL STANDARDS
 AIR.N2 <- 0.0036782 # Isotope ratio of AIR-N2 from Werner & Brandt 2001
 V.SMOW <- 0.0020052 # Isotope ratio of V-SMOW from Werner & Brandt 2001

# SAE ETH Lab standards
 ETHZSAEHIGH1.d15Na <-  0.00 #plusminus 0.32
 ETHZSAEHIGH1.d15Nb <-  2.16 #plusminus 0.33
 ETHZSAEHIGH1.d18O  <- 38.98 #plusminus 0.33

 ETHZSAEHIGH7.d15Na <- 34.01 #plusminus 0.71
 ETHZSAEHIGH7.d15Nb <- 35.98 #plusminus 0.50
 ETHZSAEHIGH7.d18O  <- 37.99 #plusminus 0.35


 ETHZSAEHIGH1.R456  <- (ETHZSAEHIGH1.d15Na/1000 + 1) * AIR.N2 #conversion from d-value to ratio
 ETHZSAEHIGH1.R546  <- (ETHZSAEHIGH1.d15Nb/1000 + 1) * AIR.N2 #conversion from d-value to ratio
 ETHZSAEHIGH1.R448  <- (ETHZSAEHIGH1.d18O/ 1000 + 1) * V.SMOW #conversion from d-value to ratio


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


 filename.str                   <-  list.files(path = data.folder, pattern="*.str")   # relative Data path from package

      STR     <- fread(paste0(data.folder,filename.str[1]),skip=1,fill=TRUE)
      
      for(i in seq(along=filename.str)[-1]) {
        x <- fread(paste0(data.folder,filename.str[i]), header=F, skip=1,fill=TRUE)
        STR <- rbind(STR, x) # appends data from new files to data frame 'QCL.stc'
      }
      
      names(STR) <-  c("time","spec.546a","spec.456a","spec.446a","spec.h2oa","spec.co2a","spec.448a","spec.n2o","spec.446b","spec.co2b","spec.CO","spec.h2ob")      
      
      STR$TIMESTAMP    <- as_datetime(STR$time, origin = "1904-01-01 UTC")
                  
#  Reads in all QCL *.stc files and creates a subset data frame
#  ==============================================================================================================================================================================================================================================      

 filename.stc                   <-  list.files(path = data.folder ,pattern="*.stc") # relative Data path from package
    
#  Reads in all QCL *.stc files
#  ==============================================================================================================================================================================================================================================        
          QCL.stc     <- fread(paste0(data.folder,filename.stc[1]),skip=1,fill=TRUE)

          for(i in seq(along=filename.stc)[-1]) {
              x <- fread(paste0(data.folder,filename.stc[i]), header=T, skip=1,fill=TRUE)
              QCL.stc <- rbind(QCL.stc, x) # appends data from new files to data frame 'QCL.stc'
              }

          names(QCL.stc) <-  c("time","rangeF1L1","rangeF1L2","rangeF2L1","rangeF2L2","Pcell","Tcell","Pref","Tref","AD8","AD9","AD10","AD11","AD12","AD13","AD14","AD15",
                            "statusW","VICI_W","USBByte","NI6024Byte","SPEFile","Tlaser1","Vlaser1","LWlaser1","Tlaser2","Vlaser2","LWlaser2","dT1","dT2","X1","pos1","X2","pos2")


        # creates data subset from QCL and removes QCL
        STC <- QCL.stc[,c("time","rangeF1L1","rangeF1L2","rangeF2L1","rangeF2L2","Pcell","Tcell","Pref","Tref","AD8","AD9","AD10","AD11","AD12","AD13","AD14","AD15",
                              "statusW","VICI_W","Tlaser1","Vlaser1","LWlaser1","Tlaser2","Vlaser2","LWlaser2","dT1","dT2","X1","pos1","X2","pos2")]
        
        STC$TIMESTAMP     <- as_datetime(STC$time, origin = "1904-01-01 UTC")
        

        # MERGE DATA FRAMES
        #FOO  <- merge(STC,TEMPERATURE,by ="TIMESTAMP",all=TRUE, sort=TRUE, incomparables=TRUE)
        QCL  <- merge(STC,STR        ,by ="TIMESTAMP",all=TRUE, sort=TRUE, incomparables=TRUE)
        
        # rm(TEMPERATURE)
        # rm(STR)
        # rm(STC)
        # rm(FOO)
        # rm(TEMPERATURE,STR,STC,F00,QCL.stc)
        
   
  # Aggregate the Data in flexible integration time (ideally suggested by the lowest Allan variance)
  #  ===========================================================================================================================================================================================================
        # str(QCL)
        QCL.dt <- as.data.table(QCL)  
        # QCL.dt$agg.unit <- format(QCL.dt$TIMESTAMP,format="%Y-%m-%d %H:%M")
        QCL.dt$agg.unit   <- round_date(QCL.dt$TIMESTAMP, "10 seconds")   # set integration time here
        QCL.agg           <- QCL.dt[, lapply(.SD,mean), by = agg.unit]    # data is averaged over the defined integration time above
        
        setnames(QCL.agg, c("agg.unit","TIMESTAMP"), c("TIMESTAMP","TIMESTAMP.mean"))  # replace the aggregated time frame with the original TIMESTAMP vector
        # system.time(QCL.agg[,TIMESTAMP := as.POSIXct(strptime(TIMESTAMP,format="%Y-%m-%d %H:%M:%S",tz="UTC"))])
        
        QCL <- QCL.agg # lets get back to the script
      
       
        
  # IMPLEMENT TIME VECTORS
  #  ===========================================================================================================================================================================================================
    
        QCL$DOY            <-  as.numeric(format(QCL$TIMESTAMP,format="%j"))
        # QCL[,DOY :=  as.numeric(format(QCL$TIMESTAMP,format="%j"))]  # datatable pendant
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
        QCL$ampmnumeric    <-  numeric(length(QCL[,TIMESTAMP]))
        QCL$ampmnumeric[QCL$ampm == "PM"] =  1
        QCL$calcode        <-  as.numeric(paste(QCL$DOY,QCL$ampmnumeric,sep=""))

        QCL$VICI  <- QCL$VICI_W + 1 # TDL Wintel is writing VICI multivalve positions from 0 to 15, however, in software used valve numbers are from # 1 to 16        
        QCL$N     <- !is.na(QCL$spec.446a)  # CREATE LOGICAL VECTOR TO KNOW n of array after aggregation          

        
   #  Create Selections and Filter from the complete dataset 
   #  ===========================================================================================================================================================================================================
        
        NA.RM.1  <- !is.na(QCL$calcode) & !is.na(QCL$spec.446a) & # !is.na(QCL$spec.546a) & !is.na(QCL$spec.456a) & #  & !is.na(QCL$d15Na)    & !is.na(QCL$d15Nb) & !is.na(QCL$d18O) &
          QCL$calcode != c(2740) & QCL$calcode != 2750 # take out certain dates where the dilution cal didn't work
        
        
        # OFFSET CORRECTION USING ETHZSAEHIGH1 at VICI port #3
        OFFSET <- NA.RM.1 & QCL$VICI == 3 & QCL$hourinseconds > 250 &  QCL$hourinseconds < 310 & QCL$hour == 3 |  # last minute of cal cycle at VICI#3 before dilution starts, 
          NA.RM.1 & QCL$VICI == 3 & QCL$hourinseconds > 250 &  QCL$hourinseconds < 310 & QCL$hour == 15
        
        
        # Parameter diagnostics
        
        # QCL.backup    <- QCL
        # QCL   <-  QCL.backup
        # head(QCL)
        # cor(QCL[,3:15])
        
        
        qcl.filter <- QCL %>%  filter(OFFSET)
      
        var.x <- qcl.filter$Tcell
        var.y <- qcl.filter$d18O
        fit <-  lm(var.y~var.x)
        qcl.filter %>% 
          plot_ly(x = ~var.x) %>% 
          add_markers(y = ~var.y ) %>% 
          add_lines(x = ~var.x, y = fitted(fit))  %>% 
          layout(showlegend = FALSE)
        
        # %>%
        #   add_text(x =  mean( var.x),text = round(summary(fit)$r.squared,4))
  

        # fit2 <- lm( QCL$spec.456a.mtt[OFFSET] ~ QCL$Tcell[OFFSET])
        # corrplot
        # write.table(erkan.data, "output/erkan_output.csv", sep = ",")
        
        ###  Corelation plot ###
        # cor.subset <- QCL[OFFSET,.(time =as.numeric(TIMESTAMP), rangeF1L1,rangeF1L2,Pcell,Tcell,Pref,Tref,#AD8,AD9,AD10,AD12,AD13,AD14,AD15,
        #              Tlaser1,Tlaser2,X1,X2,  # ,Vlaser1 ,Vlaser2
        #              spec.446a,spec.456a,spec.546a,spec.448a,spec.h2oa,spec.co2a,spec.n2o,
        #              # d15Na,d15Nb,d18O,  # ratios (can be omited if checked before corection)
        #              spec.446b,spec.co2b,spec.CO,spec.h2ob)]
        # cor.matrix <- cor(cor.subset)
        # 
        #  cor.mtest <- function(mat, ...) {
        #   mat <- as.matrix(mat)
        #   n <- ncol(mat)
        #   p.mat<- matrix(NA, n, n)
        #   diag(p.mat) <- 0
        #   for (i in 1:(n - 1)) {
        #     for (j in (i + 1):n) {
        #       tmp <- cor.test(mat[, i], mat[, j], ...)
        #       p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
        #     }
        #   }
        #   colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
        #   p.mat
        # }
        #  p.mat <- cor.mtest(cor.subset)
        #  
        # corrplot(cor.matrix, method = "number", order="hclust",type="upper",p.mat = p.mat, sig.level = 0.01,  number.cex = 0.75) # , insig = "blank"
        # 
        # plot_ly(QCL[OFFSET], x = ~ TIMESTAMP, y = ~ d15Na,type = "scatter", mode = "markers")  #  [OFFSET,.(AD10,spec.456a)]  x = ~ AD10,
        # fit <-  lm(QCL$AD10[OFFSET]~QCL$spec.456a[OFFSET])
        # 
        # qcl.filter <- QCL %>%  filter(OFFSET)
        # qcl.filter %>% 
        #   plot_ly(x = ~spec.456a) %>% 
        #   add_markers(y = ~AD10) %>% 
        #   add_lines(x = ~spec.456a, y = fitted(fit)) %>%
        # 
        # 
# calculate d-values and ratio based on QCL mixing ratios and using HITRAN abundancies
# =======================================================================================================================================================================
# calculate Ratios 
QCL$R456  <- (0.003641 * QCL$spec.456a)/(0.9903 * QCL$spec.446a)
QCL$R546  <- (0.003641 * QCL$spec.546a)/(0.9903 * QCL$spec.446a)
QCL$R448  <- (0.001986 * QCL$spec.448a)/(0.9903 * QCL$spec.446a)
    
# delta-values    
QCL$d15Na <- (QCL$R456 / AIR.N2 - 1) * 1000  # referenced against AIR-N2
QCL$d15Nb <- (QCL$R546 / AIR.N2 - 1) * 1000  # referenced against AIR-N2
QCL$d18O  <- (QCL$R448 / V.SMOW - 1) * 1000  # referenced against V-SMOW

# Parameter corrections
# =======================================================================================================================================================================
## Tcell correction on 456
QCL[OFFSET,d15Na.orig := d15Na]
QCL[OFFSET,d15Na := d15Na.orig - (Tcell- mean(Tcell)) * coef(lm(d15Na.orig ~ Tcell))[2]]
QCL[OFFSET,d15Na := d15Na - (rangeF1L1- mean(rangeF1L1)) * coef(lm(d15Na ~ rangeF1L1))[2]]
## X1 correction on 546
QCL[OFFSET,d15Nb.orig := d15Nb]
QCL[OFFSET,d15Nb := d15Nb.orig - (X1- mean(X1)) * coef(lm(d15Nb.orig ~ X1))[2]]
## spec.co2b correction on 448
QCL[OFFSET,d18O.orig := d18O]
QCL[OFFSET,d18O := d18O.orig - (spec.co2b- mean(spec.co2b)) * coef(lm(d18O.orig ~ spec.co2b))[2]]       


# bluk and SP
QCL$bulk  <- (QCL$d15Na + QCL$d15Nb) / 2         # bulk 15N is the average of alpha and beta
QCL$SP    <-  QCL$d15Na - QCL$d15Nb              # calculate site preference value




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
diff.d18O  <- ((ETHZSAEHIGH1.d18O  - QCLETHZSAEHIGH1.d18O) /(QCLETHZSAEHIGH1.d18O +1000))*1000



# CALIBRATION CALCULATION AND EQUATIONs
# =======================================================================================================================================================================

# DILUTION CORRECTION USING ETHZSAEHIGH1 at VICI port #3

 # PERIODS USEABLE FOR CORRECTION
 # VICI #3 ETHZSAEHIGH-1

     DILUTION <- NA.RM.1 & QCL$VICI == 3 & QCL$hourinseconds > 310 & QCL$hour == 3 | 
                NA.RM.1 & QCL$VICI == 3 & QCL$hourinseconds > 310 & QCL$hour == 15


    # library(nlme)   # package for function lmList
    QCL.DILUTION.SUBSET          <-      subset(QCL[DILUTION,], select=c("calcode","spec.446a","d15Na","d15Nb","d18O","Tcell","spec.h2oa","Pcell"))
    # DILUTION fit parameters  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     

    fit.d15Na            <-      lmList(d15Na ~ spec.446a | calcode,data = QCL.DILUTION.SUBSET)
    fit.d15Nb            <-      lmList(d15Nb ~ spec.446a | calcode,data = QCL.DILUTION.SUBSET)
    fit.d18O             <-      lmList(d18O ~  spec.446a  | calcode,data = QCL.DILUTION.SUBSET)

    # # Tcell / H2O   => there is no temperature, water or pressure dependency of the dillution calibration factor
    # mean.h2o             <- QCL.DILUTION.SUBSET[,mean(spec.h2oa),calcode]
    # mean.Tcell           <- QCL.DILUTION.SUBSET[,mean(Tcell),calcode]
    # mean.Pcell           <- QCL.DILUTION.SUBSET[,mean(Pcell),calcode]
    
    
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



QCL.DATA    <-    merge(QCL,CALIBRATION     ,by="calcode",all=TRUE, sort=FALSE, incomparables=NA)
QCL.DATA    <-    QCL.DATA[order(QCL.DATA$TIMESTAMP),] # reorder data.frame as merge in previous line fails to sort correctly 



# APPLYING CORRECTIONS

# OFFSET CORRECTION
QCL.DATA$d15Na.corr1 <- QCL.DATA$d15Na + QCL.DATA$diff.d15Na + (QCL.DATA$d15Na * QCL.DATA$diff.d15Na * 0.001)
QCL.DATA$d15Nb.corr1 <- QCL.DATA$d15Nb + QCL.DATA$diff.d15Nb + (QCL.DATA$d15Nb * QCL.DATA$diff.d15Nb * 0.001)
QCL.DATA$d18O.corr1  <- QCL.DATA$d18O  + QCL.DATA$diff.d18O  + (QCL.DATA$d18O  * QCL.DATA$diff.d18O  * 0.001)


# DILUTION CORRECTION
# based on the below assessed linear function with *.corr2 from the dilution calibration the data is corrected for with y2 = m(x2-x1) + y1 from m = (y2 - y1) / (x2 - x1)
 QCLETHZSAEHIGH1.446 # equals true concentration of tank for which values are corrected for

#dilution correction    = Anstieg von corr2 * concentration tank       - momentane conc   +  gemessene delta.corr1
 QCL.DATA$d15Na.corr2       <- QCL.DATA$slope.d15Na * (QCL.DATA$QCLETHZSAEHIGH1.446 - QCL.DATA$spec.446a)  +  QCL.DATA$d15Na.corr1
 QCL.DATA$d15Nb.corr2       <- QCL.DATA$slope.d15Nb * (QCL.DATA$QCLETHZSAEHIGH1.446 - QCL.DATA$spec.446a)  +  QCL.DATA$d15Nb.corr1
 QCL.DATA$d18O.corr2        <- QCL.DATA$slope.d18O  * (QCL.DATA$QCLETHZSAEHIGH1.446 - QCL.DATA$spec.446a)  +  QCL.DATA$d18O.corr1
 
 # validation plot 

 # plot(calcode,offset.d15Nb, ylim = c(-32,-20))  # check the variability of the intercept
 # plot(mean.Tcell$V1 , offset.d15Nb, ylim= c(-30,-20))   # is there a temperature dependency?
 # plot(mean.Tcell$V1 , slope.d15Nb, ylim= c(0,0.02))   # is there a temperature dependency?
# sd(offset.d15Nb)  # 1 sec => 1.578 ;   1.382   for 10 sec time integration; 30 sec => 2.006       ; 
 # sd(slope.d15Nb)   # 1 sec => 0.001211; 0.001303 for 10 sec time integration 30 sec => 0.001822    ; 
 


# CHAMBER VOLUME CALCULATIONS

# radius chamber cylinder 0.235 m
# volume remaining cylinder airspace in liter (0.235m^2*pi*rimtosoildistance.in.meter) * 1000 
# volume chamber lid in liter (0.5m*0.5m*0.15m) * 1000
# volume extensions in liter (0.5m*0.5m*0.5m) * 1000 

# started experiment with chamber lids and one extension thus total volume = ((0.235^2*pi*0.11) + (0.5*0.5*0.15) + (0.5*0.5*0.5)) * 1000
ch.volumes <- numeric(length(QCL.DATA[,2]))

ch.volumes[QCL.DATA$VICI == 5]  <- 162.5 + 26.5
ch.volumes[QCL.DATA$VICI == 6]  <- 162.5 + 28.7
ch.volumes[QCL.DATA$VICI == 7]  <- 162.5 + 32.6
ch.volumes[QCL.DATA$VICI == 8]  <- 162.5 + 33.2
ch.volumes[QCL.DATA$VICI == 9]  <- 162.5 + 38.0
ch.volumes[QCL.DATA$VICI == 10] <- 162.5 + 33.2
ch.volumes[QCL.DATA$VICI == 11] <- 162.5 + 31.8
ch.volumes[QCL.DATA$VICI == 12] <- 162.5 + 38.4
ch.volumes[QCL.DATA$VICI == 13] <- 162.5 + 35.9
ch.volumes[QCL.DATA$VICI == 14] <- 162.5 + 35.0
ch.volumes[QCL.DATA$VICI == 15] <- 162.5 + 33.0
ch.volumes[QCL.DATA$VICI == 16] <- 162.5 + 35.1


ROUND1  <- as.POSIXct("2015-11-06 00:00:00, tz=UTC")
ROUND2  <- as.POSIXct("2015-11-13 00:00:00, tz=UTC") 
ROUND3  <- as.POSIXct("2015-12-11 00:00:00, tz=UTC") 
NA.RM.2  <-  !is.na(QCL.DATA$TIMESTAMP) & !is.na(QCL.DATA$VICI) 


ch.volumes[QCL.DATA$TIMESTAMP > ROUND1 & QCL.DATA$VICI == 8  & NA.RM.2] <- ch.volumes[QCL.DATA$TIMESTAMP > ROUND1 & QCL.DATA$VICI == 8 & NA.RM.2]  + (0.5*0.5*0.5)*1000 # extension on chamber 4 on 6 Nov 2015
ch.volumes[QCL.DATA$TIMESTAMP > ROUND1 & QCL.DATA$VICI == 16 & NA.RM.2] <- ch.volumes[QCL.DATA$TIMESTAMP > ROUND1 & QCL.DATA$VICI == 16 & NA.RM.2] + (0.5*0.5*0.5)*1000 # extension on chamber 12 on 6 Nov 2015

ch.volumes[QCL.DATA$TIMESTAMP > ROUND2 & QCL.DATA$VICI == 11 & NA.RM.2] <- ch.volumes[QCL.DATA$TIMESTAMP > ROUND2 & QCL.DATA$VICI == 11 & NA.RM.2] + (0.5*0.5*0.5)*1000 # extension on chamber 7 on 13 Nov 2015
ch.volumes[QCL.DATA$TIMESTAMP > ROUND2 & QCL.DATA$VICI == 13 & NA.RM.2] <- ch.volumes[QCL.DATA$TIMESTAMP > ROUND2 & QCL.DATA$VICI == 13 & NA.RM.2] + (0.5*0.5*0.5)*1000 # extension on chamber 9 on 13 Nov 2015

ch.volumes[QCL.DATA$TIMESTAMP > ROUND3 & QCL.DATA$VICI == 6 & NA.RM.2]  <- ch.volumes[QCL.DATA$TIMESTAMP > ROUND3 & QCL.DATA$VICI == 6 & NA.RM.2]  + (0.5*0.5*0.5)*1000 # extension on chamber 2  on 11 Dec 2015 
ch.volumes[QCL.DATA$TIMESTAMP > ROUND3 & QCL.DATA$VICI == 9 & NA.RM.2]  <- ch.volumes[QCL.DATA$TIMESTAMP > ROUND3 & QCL.DATA$VICI == 9 & NA.RM.2]  + (0.5*0.5*0.5)*1000 # extension on chamber 5 on 11 Dec 2015 

QCL.DATA$ch.volumes <- ch.volumes


# ===============================================================================================================================================================================================================================================
# FLUX CALCULATIONS (linear regression)
NA.RM.3  <-  !is.na(QCL.DATA$spec.446a) & !is.na(QCL.DATA$VICI)

# QCL.DATA$N2O.nmol[NA.RM.3 & QCL.DATA$VICI == 5]       <- 1.00035 * QCL.DATA$spec.446a[NA.RM.3 & QCL.DATA$VICI == 5]  / (0.08206 * (273.15 + QCL.DATA$temp_C1_mod[NA.RM.3 & QCL.DATA$VICI  == 5]))  # nmol/L
# QCL.DATA$N2O.nmol[NA.RM.3 & QCL.DATA$VICI == 6]       <- 1.00035 * QCL.DATA$spec.446a[NA.RM.3 & QCL.DATA$VICI == 6]  / (0.08206 * (273.15 + QCL.DATA$temp_C2_mod[NA.RM.3 & QCL.DATA$VICI  == 6]))  # nmol/L
# QCL.DATA$N2O.nmol[NA.RM.3 & QCL.DATA$VICI == 7]       <- 1.00035 * QCL.DATA$spec.446a[NA.RM.3 & QCL.DATA$VICI == 7]  / (0.08206 * (273.15 + QCL.DATA$temp_C3_mod[NA.RM.3 & QCL.DATA$VICI  == 7]))  # nmol/L
# QCL.DATA$N2O.nmol[NA.RM.3 & QCL.DATA$VICI == 8]       <- 1.00035 * QCL.DATA$spec.446a[NA.RM.3 & QCL.DATA$VICI == 8]  / (0.08206 * (273.15 + QCL.DATA$temp_C4_mod[NA.RM.3 & QCL.DATA$VICI  == 8]))  # nmol/L
# QCL.DATA$N2O.nmol[NA.RM.3 & QCL.DATA$VICI == 9]       <- 1.00035 * QCL.DATA$spec.446a[NA.RM.3 & QCL.DATA$VICI == 9]  / (0.08206 * (273.15 + QCL.DATA$temp_C5_mod[NA.RM.3 & QCL.DATA$VICI  == 9]))  # nmol/L
# QCL.DATA$N2O.nmol[NA.RM.3 & QCL.DATA$VICI == 10]      <- 1.00035 * QCL.DATA$spec.446a[NA.RM.3 & QCL.DATA$VICI == 10] / (0.08206 * (273.15 + QCL.DATA$temp_C6_mod[NA.RM.3 & QCL.DATA$VICI  == 10])) # nmol/L
# QCL.DATA$N2O.nmol[NA.RM.3 & QCL.DATA$VICI == 11]      <- 1.00035 * QCL.DATA$spec.446a[NA.RM.3 & QCL.DATA$VICI == 11] / (0.08206 * (273.15 + QCL.DATA$temp_C7_mod[NA.RM.3 & QCL.DATA$VICI  == 11])) # nmol/L
# QCL.DATA$N2O.nmol[NA.RM.3 & QCL.DATA$VICI == 12]      <- 1.00035 * QCL.DATA$spec.446a[NA.RM.3 & QCL.DATA$VICI == 12] / (0.08206 * (273.15 + QCL.DATA$temp_C8_mod[NA.RM.3 & QCL.DATA$VICI  == 12])) # nmol/L
# QCL.DATA$N2O.nmol[NA.RM.3 & QCL.DATA$VICI == 13]      <- 1.00035 * QCL.DATA$spec.446a[NA.RM.3 & QCL.DATA$VICI == 13] / (0.08206 * (273.15 + QCL.DATA$temp_C9_mod[NA.RM.3 & QCL.DATA$VICI  == 13])) # nmol/L
# QCL.DATA$N2O.nmol[NA.RM.3 & QCL.DATA$VICI == 14]      <- 1.00035 * QCL.DATA$spec.446a[NA.RM.3 & QCL.DATA$VICI == 14] / (0.08206 * (273.15 + QCL.DATA$temp_C10_mod[NA.RM.3 & QCL.DATA$VICI == 14])) # nmol/L
# QCL.DATA$N2O.nmol[NA.RM.3 & QCL.DATA$VICI == 15]      <- 1.00035 * QCL.DATA$spec.446a[NA.RM.3 & QCL.DATA$VICI == 15] / (0.08206 * (273.15 + QCL.DATA$temp_C11_mod[NA.RM.3 & QCL.DATA$VICI == 15])) # nmol/L
# QCL.DATA$N2O.nmol[NA.RM.3 & QCL.DATA$VICI == 16]      <- 1.00035 * QCL.DATA$spec.446a[NA.RM.3 & QCL.DATA$VICI == 16] / (0.08206 * (273.15 + QCL.DATA$temp_C12_mod[NA.RM.3 & QCL.DATA$VICI == 16])) # nmol/L
QCL.DATA$N2O.nmol[NA.RM.3 ]      <- 1.00035 * QCL.DATA$spec.446a[NA.RM.3] / (0.08206 * (273.15 + 20)) # nmol/L  - no variable temperature data



# QCL.DATA$CO2.nmol[NA.RM.3 & QCL.DATA$VICI == 5]       <- 1.00035 * QCL.DATA$spec.co2b[NA.RM.3 & QCL.DATA$VICI == 5]  / (0.08206 * (273.15 + QCL.DATA$temp_C1_mod[NA.RM.3 & QCL.DATA$VICI  == 5]))  # nmol/L
# QCL.DATA$CO2.nmol[NA.RM.3 & QCL.DATA$VICI == 6]       <- 1.00035 * QCL.DATA$spec.co2b[NA.RM.3 & QCL.DATA$VICI == 6]  / (0.08206 * (273.15 + QCL.DATA$temp_C2_mod[NA.RM.3 & QCL.DATA$VICI  == 6]))  # nmol/L
# QCL.DATA$CO2.nmol[NA.RM.3 & QCL.DATA$VICI == 7]       <- 1.00035 * QCL.DATA$spec.co2b[NA.RM.3 & QCL.DATA$VICI == 7]  / (0.08206 * (273.15 + QCL.DATA$temp_C3_mod[NA.RM.3 & QCL.DATA$VICI  == 7]))  # nmol/L
# QCL.DATA$CO2.nmol[NA.RM.3 & QCL.DATA$VICI == 8]       <- 1.00035 * QCL.DATA$spec.co2b[NA.RM.3 & QCL.DATA$VICI == 8]  / (0.08206 * (273.15 + QCL.DATA$temp_C4_mod[NA.RM.3 & QCL.DATA$VICI  == 8]))  # nmol/L
# QCL.DATA$CO2.nmol[NA.RM.3 & QCL.DATA$VICI == 9]       <- 1.00035 * QCL.DATA$spec.co2b[NA.RM.3 & QCL.DATA$VICI == 9]  / (0.08206 * (273.15 + QCL.DATA$temp_C5_mod[NA.RM.3 & QCL.DATA$VICI  == 9]))  # nmol/L
# QCL.DATA$CO2.nmol[NA.RM.3 & QCL.DATA$VICI == 10]      <- 1.00035 * QCL.DATA$spec.co2b[NA.RM.3 & QCL.DATA$VICI == 10] / (0.08206 * (273.15 + QCL.DATA$temp_C6_mod[NA.RM.3 & QCL.DATA$VICI  == 10])) # nmol/L
# QCL.DATA$CO2.nmol[NA.RM.3 & QCL.DATA$VICI == 11]      <- 1.00035 * QCL.DATA$spec.co2b[NA.RM.3 & QCL.DATA$VICI == 11] / (0.08206 * (273.15 + QCL.DATA$temp_C7_mod[NA.RM.3 & QCL.DATA$VICI  == 11])) # nmol/L
# QCL.DATA$CO2.nmol[NA.RM.3 & QCL.DATA$VICI == 12]      <- 1.00035 * QCL.DATA$spec.co2b[NA.RM.3 & QCL.DATA$VICI == 12] / (0.08206 * (273.15 + QCL.DATA$temp_C8_mod[NA.RM.3 & QCL.DATA$VICI  == 12])) # nmol/L
# QCL.DATA$CO2.nmol[NA.RM.3 & QCL.DATA$VICI == 13]      <- 1.00035 * QCL.DATA$spec.co2b[NA.RM.3 & QCL.DATA$VICI == 13] / (0.08206 * (273.15 + QCL.DATA$temp_C9_mod[NA.RM.3 & QCL.DATA$VICI  == 13])) # nmol/L
# QCL.DATA$CO2.nmol[NA.RM.3 & QCL.DATA$VICI == 14]      <- 1.00035 * QCL.DATA$spec.co2b[NA.RM.3 & QCL.DATA$VICI == 14] / (0.08206 * (273.15 + QCL.DATA$temp_C10_mod[NA.RM.3 & QCL.DATA$VICI == 14])) # nmol/L
# QCL.DATA$CO2.nmol[NA.RM.3 & QCL.DATA$VICI == 15]      <- 1.00035 * QCL.DATA$spec.co2b[NA.RM.3 & QCL.DATA$VICI == 15] / (0.08206 * (273.15 + QCL.DATA$temp_C11_mod[NA.RM.3 & QCL.DATA$VICI == 15])) # nmol/L
# QCL.DATA$CO2.nmol[NA.RM.3 & QCL.DATA$VICI == 16]      <- 1.00035 * QCL.DATA$spec.co2b[NA.RM.3 & QCL.DATA$VICI == 16] / (0.08206 * (273.15 + QCL.DATA$temp_C12_mod[NA.RM.3 & QCL.DATA$VICI == 16])) # nmol/L
QCL.DATA$CO2.nmol[NA.RM.3 ]      <- 1.00035 * QCL.DATA$spec.co2b[NA.RM.3] / (0.08206 * (273.15 + 20)) # nmol/L   - no variable temperature data


QCL.DATA$VICI.time     <- paste(QCL.DATA$DOY,QCL.DATA$HOUR,QCL.DATA$VICI,sep = "") # use "HOUR" instead of "hour" as former is including 0 when counting from 1-9 
QCL.DATA$LABEL         <- paste(QCL.DATA$weekday ,QCL.DATA$ampm,sep = "_")

QCL.DATA$spec.446a.inv <- 1/ QCL.DATA$spec.446a


## outtakes NA
NA.RM.4 <- !is.na(QCL.DATA$hourinseconds) & !is.na(QCL.DATA$VICI.time)   & !is.na(QCL.DATA$N2O.nmol)    & !is.na(QCL.DATA$CO2.nmol) & 
           !is.na(QCL.DATA$spec.446a.inv) & !is.na(QCL.DATA$d15Na.corr2) & !is.na(QCL.DATA$d15Nb.corr2) & !is.na(QCL.DATA$d18O.corr2) & !is.na(QCL.DATA$LABEL)
NA.RM.4 <- as.vector(NA.RM.4)


# In the following respective Keeling plots are done for all chambers
KEELING <- NA.RM.4 & 
           QCL.DATA$hourinseconds > 1200 & QCL.DATA$hourinseconds < 3540 & QCL.DATA$VICI > 4 & # start keeling plots 20min into the hour for all VICI ports > 4 (all chambers)
           QCL.DATA$LABEL != "4_PM" & QCL.DATA$LABEL != "5_AM"         # exclude 13C label periods Friday morning

KEELING <- as.vector(KEELING)

# library(nlme)   # package for function lmList  moved to top

    DATA.KEELING.SUBSET      <-      na.omit(subset(QCL.DATA[KEELING,], select=c("hourinseconds","VICI.time","N2O.nmol","CO2.nmol","spec.446a.inv","d15Na.corr2","d15Nb.corr2","d18O.corr2")))
   
    # DILUTION fit parameters  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++     
    volume                   <-      tapply(QCL.DATA$ch.volumes[KEELING],QCL.DATA$VICI.time[KEELING],mean)
    flux.time                <-      tapply(QCL.DATA$TIMESTAMP[KEELING],QCL.DATA$VICI.time[KEELING],mean)
    ch.time                  <-      ISOdatetime(1970,1,1,0,0,0) + flux.time    
    port                     <-      tapply(QCL.DATA$VICI[KEELING],QCL.DATA$VICI.time[KEELING],mean)
    
    # Keeling plot for d15Na
    keeling.d15Na            <-      lmList(d15Na.corr2 ~ spec.446a.inv | VICI.time,data = DATA.KEELING.SUBSET)    
    keeling.R2.d15Na         <-      summary(keeling.d15Na)$r.squared
    keeling.df.d15Na         <-      summary(keeling.d15Na)$df[,2]
    keeling.offset.d15Na     <-      coef(keeling.d15Na)[,1]
    keeling.slope.d15Na      <-      coef(keeling.d15Na)[,2]

    # plot(ch.time,keeling.R2.d15Na)
    
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
              QCL.DATA$hourinseconds > 2000 & QCL.DATA$hourinseconds < 3540 & QCL.DATA$VICI > 4 & # start keeling plots 20min into the hour for all VICI ports > 4 (all chambers)
              QCL.DATA$LABEL != "4_PM" & QCL.DATA$LABEL != "5_AM"         # exclude 13C label periods Friday morning


FLUX.CO2    <- NA.RM.4 & 
              QCL.DATA$hourinseconds >  900 & QCL.DATA$hourinseconds < 2000 & QCL.DATA$VICI > 4 & # start keeling plots 20min into the hour for all VICI ports > 4 (all chambers)
              QCL.DATA$LABEL != "4_PM" & QCL.DATA$LABEL != "5_AM"         # exclude 13C label periods Friday morning


FLUX.N2O    <- as.vector(FLUX.N2O)
FLUX.CO2    <- as.vector(FLUX.CO2)

# library(nlme)   # package for function lmList  -> moved to top

    DATA.FLUX.SUBSET.N2O      <-      na.omit(subset(QCL.DATA[FLUX.N2O,], select=c("hourinseconds","VICI.time","N2O.nmol","CO2.nmol","spec.446a.inv","d15Na.corr2","d15Nb.corr2","d18O.corr2","spec.446a")))
    DATA.FLUX.SUBSET.CO2      <-      na.omit(subset(QCL.DATA[FLUX.CO2,], select=c("hourinseconds","VICI.time","N2O.nmol","CO2.nmol","spec.446a.inv","d15Na.corr2","d15Nb.corr2","d18O.corr2","spec.446a")))
        
    
    # CO2 Flux calculation
    co2.flux.regression      <-      lmList(CO2.nmol ~ hourinseconds | VICI.time,data = DATA.FLUX.SUBSET.CO2)
    co2.R2                   <-      summary(co2.flux.regression)$r.squared
    co2.flux                 <-      coef(co2.flux.regression)[,2]  * (volume/(0.235^2*pi))

    # N2O Flux calculation          
    n2o.flux.regression      <-      lmList(N2O.nmol ~ hourinseconds | VICI.time,data = DATA.FLUX.SUBSET.N2O)
    n2o.R2                   <-      summary(n2o.flux.regression)$r.squared    
    n2o.flux                 <-      coef(n2o.flux.regression)[,2]  * (volume/(0.235^2*pi))

    
   # plot(n2o.flux) # check output
   # plot(co2.flux) # check output
    
    #summary(keeling.d15Nb)
    KEELING.RES <- data.frame(ch.time[n2o.flux > 0.9],
                              port[n2o.flux > 0.9],
                              n2o.flux[n2o.flux > 0.9],
                              coef(keeling.d15Na)[,1][n2o.flux > 0.9],
                              coef(keeling.d15Nb)[,1][n2o.flux > 0.9],
                              coef(keeling.d18O) [,1][n2o.flux > 0.9])
    
    colnames(KEELING.RES) <- c("time","VICI","flux","d15Na","d15Nb","d18O")
    KEELING.RES$d15Nbulk  <- (KEELING.RES$d15Na + KEELING.RES$d15Nb)/2
    KEELING.RES$SP        <-  KEELING.RES$d15Na - KEELING.RES$d15Nb
    
   
   # FILTERING OF unreasonable CO2 values as general criteria 
   ch01 <- port == 5  & co2.flux <  15000 & co2.flux >  -15000 
   ch02 <- port == 6  & co2.flux <  15000 & co2.flux >  -15000
   ch03 <- port == 7  & co2.flux <  15000 & co2.flux >  -15000
   ch04 <- port == 8  & co2.flux <  15000 & co2.flux >  -15000 
   ch05 <- port == 9  & co2.flux <  15000 & co2.flux >  -15000 
   ch06 <- port == 10 & co2.flux <  15000 & co2.flux >  -15000 
   ch07 <- port == 11 & co2.flux <  15000 & co2.flux >  -15000 
   ch08 <- port == 12 & co2.flux <  15000 & co2.flux >  -15000 
   ch09 <- port == 13 & co2.flux <  15000 & co2.flux >  -15000 
   ch10 <- port == 14 & co2.flux <  15000 & co2.flux >  -15000 
   ch11 <- port == 15 & co2.flux <  15000 & co2.flux >  -15000 
   ch12 <- port == 16 & co2.flux <  15000 & co2.flux >  -15000
   
  
   
   # CREATE MEAN of all treatments, 200 is common length of each chamber measurement    
   
   ZINAL  <- data.frame(co2.flux[ch01][0:200],co2.flux[ch10][0:200],co2.flux[ch11][0:200]);ZINAL.co2.flux <- rowMeans(ZINAL)
   CLARO  <- data.frame(co2.flux[ch03][0:200],co2.flux[ch06][0:200],co2.flux[ch08][0:200]);CLARO.co2.flux <- rowMeans(CLARO)
   PROBUS <- data.frame(co2.flux[ch02][0:200],co2.flux[ch05][0:200],co2.flux[ch09][0:200]);PROBUS.co2.flux<- rowMeans(PROBUS)
   MONTE  <- data.frame(co2.flux[ch04][0:200],co2.flux[ch07][0:200],co2.flux[ch12][0:200]);MONTE.co2.flux <- rowMeans(MONTE)
   
   
   FLUXDATA <- data.frame(
     ch.time[ch01][0:200],n2o.flux[ch01][0:200],
     ch.time[ch02][0:200],n2o.flux[ch02][0:200],
     ch.time[ch03][0:200],n2o.flux[ch03][0:200],
     ch.time[ch04][0:200],n2o.flux[ch04][0:200],
     ch.time[ch05][0:200],n2o.flux[ch05][0:200],
     ch.time[ch06][0:200],n2o.flux[ch06][0:200],
     ch.time[ch07][0:200],n2o.flux[ch07][0:200],
     ch.time[ch08][0:200],n2o.flux[ch08][0:200],
     ch.time[ch09][0:200],n2o.flux[ch09][0:200],
     ch.time[ch10][0:200],n2o.flux[ch10][0:200],
     ch.time[ch11][0:200],n2o.flux[ch11][0:200],
     ch.time[ch12][0:200],n2o.flux[ch12][0:200]
   )
   # head(FLUXDATA)
   
   ### write results to output folder
   # fwrite(QCL.DATA    ,"output/QCL-DATA.csv"    ,sep=",")       # write processed rawdata into output file
   # 
   # fwrite(FLUXDATA,"output/FLUXDATA.csv",sep=",")               # write calculated fluxes into output file
   # 
   # fwrite(KEELING.RES,"output/KEELING-RESULTS.csv",sep=",")     # write keeling plot results
    
   