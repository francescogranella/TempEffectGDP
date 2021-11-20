# # Persistent Effect of Temperature on GDP Identified from Lower Frequency Temperature Variability
# # Bastien-Olvera, B. A. and Moore, F. C.

# ## Index
#     # 1. Setup
#     # 2. Simulation (Figure 1 and 2)
#         # 2.1. Temperature fluctuations and filters - Figure 1
#         # 2.2. Global temperature profile
#         # 2.3. Perform MAIN simulation - Figure 2
#         # 2.4. Additional simulations used for Response to Reviewers
#     # 3. Empirical analysis (Figure 3 and Table 1)
#         # 3.1 Country-level regressions
#         # 3.2. Panel Regression
#         # 3.2. Categorizing - Figure 3
#         # 3.3. FELM abs estimate by filter - Table 1, columns 1 and 2
#         # 3.4. Mean effect accross filters and socioeconomic variables; Supp Fig 1 and 2
#         # 3.5. Other datasets - Table 1, columns 3 to 6; Supp Fig 3

## 1. Setup (start)
    #set working directory
    dir <- "C:/Users/bastien/Documents/GitHub/TempEffectGDP/"
    library(here)
    dir <- here()
    setwd(dir)
    x<-c("sandwich","TTR", "ggpubr","mFilter", "ncdf4","TSA","tidyverse","ggplot2",
        "ggpubr","dplR","reshape2","raster","dplyr",
        "RColorBrewer","colorspace","spData","sf","countrycode", "ggstatsplot",
        "ggsignif","dlnm","lfe","directlabels","splines","timeSeries",
        "MASS","stargazer","jtools","wbstats","dplR","lmtest","sandwich")
    lapply(x, require, character.only = TRUE)
    # some functions and color pallettes
    tsyears <- function(ts) as.numeric(trunc(time(ts)))
    cols <- c("#38170B","#BF1B0B", "#FFC465", "#66ADE5", "#252A52")
    cols=rev(c("#76d3ae","#0dcfca","#055692"))

    #Load data
    barro <- read.csv("Data/Barro_UDel.csv") #Barro-Ursua and University of Delaware
    mad <- read.csv("Data/Maddison_UDel.csv") #Maddison Project and University of Delaware
    wb <- read.csv("Data/WB_UDel.csv")   #World Bank and University of Delaware
    bhm <- read.csv("Data/BHM.csv") #BHM dataset (World Bank and University of Delaware)
    

    #Code to merge World Bank and Temperature Data data (to generate the file WB_UDel.csv; currently commented)
        #wbg <- read.csv(paste(dir,"Data/GDPpc_2015usd_WDI.csv",sep=""))
        #udel_t <- read.csv(paste(dir,"Data/UDel_T_popaverage.csv",sep=""))
        #udel_p <- read.csv(paste(dir,"Data/UDel_P_popaverage.csv",sep=""))
        # wbg <- wbg[,-c(1)]
        # names(wbg) <- c("Country.Code",   "Series.Name"   , "Series.Code", sprintf("%s",seq(1960:2020)))
        # years <- sprintf("%s",seq(1960:2020))
        # #arranging data
        # for (i in 1:length(years)){
        #     yeari <- years[i]
        #     w2 <- wbg[,which(colnames(wbg) %in% 
        #         c("Country.Code","Series.Name",yeari))]
        #     w3 <- reshape(w2, idvar="Country.Code", 
        #         timevar="Series.Name", direction="wide")
        #     w3<-w3[,c(1:2)]
        #     w3[,(dim(w3)[2]+1)] <- yeari
        #     colnames(w3) <- c("countrycode","gdppc","year")
        #     if (i==1){Wbg <- w3} else {
        #         Wbg <- rbind(Wbg,w3)
        #     }
        # }
        # glimpse(Wbg)
        # Wbg$year <- as.double(Wbg$year)+1959
        
        # Wbg <- Wbg[order( Wbg[,1], Wbg[,3] ),] # Sort by column index [1] then [3]
        # Wbg$gdppc <- as.double(Wbg$gdppc)
        # Wbg <- Wbg[-which(Wbg$countrycode==""),]

        # growth <- Wbg$gdppc*NA
        # for(i in 2:(nrow(Wbg)-1)) {
        #         if (Wbg$countrycode[i]==Wbg$countrycode[i-1]){
        #             growth[i] <- log(Wbg$gdppc[i]) - log(Wbg$gdppc[i-1])
        #         }
        #         }
        # Wbg$growth <- growth
        # wb <- merge(Wbg,udel_p, by=c("countrycode","year"))
        # wb <- merge(wb,udel_t, by=c("countrycode","year"))
        # glimpse(wb)
        # wb <- wb[,-which(names(wb)%in%c("X.x","X.y"))]
        # names(wb)<-c( "countrycode", "year" ,"gdppc","growth" ,"UDel_pop_preci" , "UDel_pop_temp")
        # write.csv(wb,file=paste(dir,"Data/WB_UDel.csv",sep=""))
        # #Compare with BHM
        
        #     bhm_data <- bhm[,which(names(bhm) %in% c("iso","year","growthWDI","rgdpCAPgr","UDel_temp_popweight","UDel_precip_popweight"))]
        #     names(bhm_data) <- c("countrycode","year","growthWDI","rgdpCAPgr","t_bhm_udel","p_bhm_udel")
        #     merged_bhm <- merge(wb,bhm_data,by=c("countrycode","year"))   

        #     summary(felm(growth~growthWDI|countrycode|0|0,data=merged_bhm))

        #     ggplot(data=merged_bhm)+
        #         geom_point(aes(x=growth,y=growthWDI))
            
        #     summary(felm(t_bhm_udel~UDel_pop_temp|countrycode|0|0,data=merged_bhm))

        #     ggplot(data=merged_bhm)+
        #         geom_point(aes(x=t_bhm_udel,y=UDel_pop_temp))
        #Compare data
    #Code to merge World Bank and Temperature Data data (to generate the file WB_UDel.csv currently commented)
   
## 1. Setup (end)

## 2. Simulation - Figures 1 and 2 (start)

    ## 2.1. Temperature fluctuations and filters - Figure 1 (start)
        
        wbUS <- wb[wb$countrycode=="USA",]
        mt <- lm(UDel_pop_temp~year, data = wbUS)
        t <- resid(mt)    
        t <- interpNA(t, method = "linear")
        t <- removeNA(t)
        tdf <- data.frame(year = mt$model[,2], temp= t+1.8, Filter = "Unfiltered", mean = 1.8)
        t3 <- pass.filt(t, W=3, type="low", method="Butterworth")
        tdf <- rbind(tdf,data.frame(year = mt$model[,2], temp= t3+1, Filter = "2-3 years", mean = 1))
        t5 <- pass.filt(t, W=5, type="low", method="Butterworth")
        tdf <- rbind(tdf,data.frame(year = mt$model[,2], temp= t5+0.5, Filter = "2-5 years", mean = 0.5))
        t10 <- pass.filt(t, W=10, type="low", method="Butterworth")
        tdf <- rbind(tdf,data.frame(year = mt$model[,2], temp= t10+0.25, Filter = "2-10 years", mean = 0.25))
        t15 <- pass.filt(t, W=15, type="low", method="Butterworth")
        tdf <- rbind(tdf,data.frame(year = mt$model[,2], temp= t15+0.1, Filter = "2-15 years", mean = 0.1))
        
        tdf$Filter <- factor(tdf$Filter, levels = c("Unfiltered", "2-3 years", "2-5 years","2-10 years","2-15 years"))
        

        ggplot(data=tdf, aes(x=year, y=temp, color = Filter))+geom_line()+theme_bw()+
            theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border = element_blank(),plot.title = element_text(hjust = 0.5))+
            geom_errorbar(aes(x=1958,ymin=1,ymax=2), color = "black")+
            geom_text(aes(x=1960,y=1.5),label="1C", color="black")+
            #geom_dl(aes(label = Filter), method = list(dl.combine("last.points")), cex = 0.9) +
            xlim(1958,2018)+ ggtitle("United States temperature fluctuation")+ labs(color= "Filtered oscillations \n (periods)")+
            geom_line(data = tdf, aes(x=year, y =mean, color = Filter))
    ## 2.1. Temperature fluctuations and filters - Figure 1 (end)

    ## 2.2. Global temperature profile (start)
        
        #Read and save global temperature file. NetCDF too large to be in GitHub

        #file=nc_open(paste(dir,"Data/gmt_MCruns_ensemble_full.nc",sep=""))
        #gtemp=ncvar_get(file,"gmt")
        #gtemp=apply(gtemp,MARGIN=3,FUN=mean)
        #saveRDS(gtemp, file = "Data/gtemp.rds")
        
        # Uncomment for US temp profile (end)
            #file=nc_open(paste(dir,"Data/air_MCruns_ensemble_mean_LMRv2.0.nc",sep=""))
            #gtemp=ncvar_get(file,"air",start=c(1,1,1,1),count=c(-1,-1,1,-1))
            #lon <- ncvar_get(file,"lon")
            #lat <- ncvar_get(file,"lat")
            #time <- ncvar_get(file,"time")
            #tunits <- ncatt_get(file,"time","units")
            #data(world)
            #geom_iso <- world$geom[world$iso_a2=="US"]
            #geom_iso <- st_shift_longitude(geom_iso)
            #geom_iso <- st_cast(geom_iso, "POLYGON")
            #geom_iso <-as_Spatial(geom_iso)
            
            #gtemp <- brick(gtemp, xmn=min(lat), xmx=max(lat), ymn=min(lon), ymx=max(lon), crs=CRS("+proj=longlat +datum=WGS84 +no_defs"))
            #ustemp_lmr <- rep(NA,1500)
            #for (i in 1:1500){
            #    gtemp_year <- raster::subset(gtemp,i)
            #   #gtemp_year <- rotate(gtemp_year)
            #  gtemp_year <- t(gtemp_year)
            # ustemp <- extract(gtemp_year, geom_iso,method='simple')
                #ustemp_lmr[i] <- mean(ustemp[[1]])

            #}
            
            #print(file)
            #ustemp_lmr_orig <- ustemp_lmr
            #ustemp_lmr=lm(ustemp_lmr_orig~x)$residuals
            #p=periodogram(ustemp_lmr)
        
        # Uncomment for US temp profile (end)
        
        # Load gtemp
        gtemp <- readRDS(file = "Data/gtemp.rds")
        
        #take first 1500 years, prior to anthropogenic influence
        gtemp=gtemp[1:1500]

        #take out linear time trend
        x=1:1500
        gtemp=lm(gtemp~x)$residuals
        p=periodogram(gtemp)
        dd=data.frame(freq=p$freq,spec=p$spec)
        dd=dd[order(-dd$spec),]

        #look at top 50 frequencies
        dd=dd[1:50,]
        dd$period=1/dd$freq

        #create randomized time series 
        randomts=function(timeseries){
            #returns a time series with the same spectral profile as the argument, but with randomly chosen phases
            ft=fft(timeseries)
            N=length(timeseries)
            rphase=runif(N,min=0,max=360)
            newft=complex(real=abs(ft)*cos(rphase),imaginary=abs(ft)*sin(rphase))
            return(fft(newft,inverse=T)/N)
        }
    ## 2.2. Global temperature profile (end)

    ## 2.3. Perform MAIN simulation - WITH DRIFT ADJUSTMENT (start)


            time=350 #or 350 years
            basegr=0.01 #2% per year baseline growth rate
            start=100
            baseline=rep(basegr,time)
            coef=-0.05 #effect size - change in growth per degree warming()
            coef2=-0.05
            growthsd=0.005 #standard deviation of growth variability unexplained by temperature
            periods <- c(0,3,5,10,15)
            nsims=10000
            sims=array(dim=c(nsims,length(periods),2))
            sims_adjusted=array(dim=c(nsims,length(periods),2))
            sims_ar=array(dim=c(nsims,length(periods),2))
            for(i in 1:nsims){
                randomtemp=Re(randomts(gtemp))[1:time]
                #randomtemp=Re(randomts(ustemp_lmr))[1:time]
                randomgrowth_g=basegr #growth impact model
                randomgrowth_l=basegr #levels impact model
                for(j in 2:time){
                #linear
                randomgrowth_g=c(randomgrowth_g, basegr+(randomtemp[j]*coef)+rnorm(1,mean=0,sd=growthsd))
                randomgrowth_l=c(randomgrowth_l, basegr+(randomtemp[j]-randomtemp[j-1])*coef+rnorm(1,mean=0,sd=growthsd))
                }
                dataset=data.frame(years=1:time,temp=randomtemp,g=randomgrowth_g,l=randomgrowth_l)
                dataset$templag =c(NA,dataset$temp[1:(dim(dataset)[1]-1)])
                for(k in 1:length(periods)){
                mod_g=lm(g~temp+templag,data=dataset)$coef[2]
                mod_l=lm(l~temp+templag,data=dataset)$coef[2]

                if (k==1){
                        sims[i,k,]=c(mod_g,mod_l)
                        sims_adjusted[i,k,]=c(mod_g,mod_l)
                        
                    }else{
                    x <- seq(1:length(randomtemp))
                    randomtemp <- lm(randomtemp~x)$residuals
                    tempts <- pass.filt(randomtemp, W=periods[k], type="low", method="Butterworth")
                    
                    ratio <- median(randomtemp/tempts)
                    growth <- randomgrowth_g
                    level <- randomgrowth_l
                    filt <- data.frame(temp = unclass(tempts), growth = unclass(growth),
                        level = unclass(level))
                    filt$templag=c(NA,filt$temp[1:(dim(filt)[1]-1)])
                    names(filt) <- c("temp","growth","levels","templag")
                    mod_gfilt=lm(growth~temp,data=filt)$coef[2]
                    mod_lfilt=lm(levels~temp,data=filt)$coef[2]
                    dw <- dwtest(lm(growth~temp,data=filt)) #test for autocorrelation
                    sims[i,k,]=c(mod_gfilt,mod_lfilt)
                    sims_adjusted[i,k,]=c(mod_gfilt/ratio,mod_lfilt/ratio)
                    sims_ar[i,k,]=c(dw$statistic[1],dw$p.value[1])
                }
                

                }
                if(i%%50==0) print(i)
            }
            #
          
                # Test - Figure 2
                    
                    simsfiltmean=apply(sims_ar,MARGIN=c(2,3),FUN=mean)
                    simsfiltsd=apply(sims_ar,MARGIN=c(2,3),FUN=sd)
                    colnames(simsfiltmean)=c("Growth","Level");rownames(simsfiltmean)=periods
                    simsfiltdat=cbind(melt(simsfiltmean),melt(simsfiltsd)[,3])
                    colnames(simsfiltdat)=c("frequencies","ImpactType","dw_estimate","dw_pvalue")
                    theme_set(theme_bw(base_size = 20))
                    
                    simsfiltdat$periodsregationPeriod <- c("Unfiltered", "10 years", "20 years", "50 years", "100 years","Unfiltered", "10 years", "20 years", "50 years", "100 years")
                    simsfiltdat$periodsregationPeriod <- factor(simsfiltdat$periodsregationPeriod,
                        levels = c("Unfiltered", "10 years", "20 years", "50 years", "100 years"))
                    simsfiltdat_ar <- simsfiltdat

                    
                    simsfiltmean=apply(sims_adjusted,MARGIN=c(2,3),FUN=mean)
                    simsfiltsd=apply(sims_adjusted,MARGIN=c(2,3),FUN=sd)
                    colnames(simsfiltmean)=c("Growth","Level");rownames(simsfiltmean)=periods
                    simsfiltdat=cbind(melt(simsfiltmean),melt(simsfiltsd)[,3])
                    colnames(simsfiltdat)=c("periodsregationPeriod","ImpactType","MeanEffect","SDEffect")
                    theme_set(theme_bw(base_size = 20))
                    
                    simsfiltdat$periodsregationPeriod <- c("Unfiltered", "3 years", "5 years", "10 years", "15 years","Unfiltered", "3 years", "5 years", "10 years", "15 years")
                    simsfiltdat$periodsregationPeriod <- factor(simsfiltdat$periodsregationPeriod,
                    levels = c("Unfiltered", "3 years", "5 years", "10 years", "15 years"))
                    simsfiltdat_adjusted <- simsfiltdat

                    simsfiltmean=apply(sims,MARGIN=c(2,3),FUN=mean)
                    simsfiltsd=apply(sims,MARGIN=c(2,3),FUN=sd)
                    colnames(simsfiltmean)=c("Growth","Level");rownames(simsfiltmean)=periods
                    simsfiltdat=cbind(melt(simsfiltmean),melt(simsfiltsd)[,3])
                    colnames(simsfiltdat)=c("periodsregationPeriod","ImpactType","MeanEffect","SDEffect")
                    theme_set(theme_bw(base_size = 20))
                    
                    simsfiltdat$periodsregationPeriod <- c("Unfiltered", "3 years", "5 years", "10 years", "15 years","Unfiltered", "3 years", "5 years", "10 years", "15 years")
                    #simsfiltdat$periodsregationPeriod <- c("Unfiltered", "10 years", "20 years", "50 years", "100 years","Unfiltered", "10 years", "20 years", "50 years", "100 years")
                    simsfiltdat$periodsregationPeriod <- factor(simsfiltdat$periodsregationPeriod,
                        levels = c("Unfiltered", "3 years", "5 years", "10 years", "15 years"))
                    
                    
                    a=ggplot(simsfiltdat_adjusted,aes(x=factor(periodsregationPeriod),y=MeanEffect*100,col=ImpactType,group=ImpactType))+geom_line(lwd=1.25)
                    a=a+geom_point(size=4)+geom_errorbar(aes(ymin=(MeanEffect-1.96*SDEffect)*100,ymax=(MeanEffect+1.96*SDEffect)*100),width=0.1,lwd=1.25)
                    a=a+geom_hline(yintercept=0,lwd=1.5,lty=3)
                    a=a+labs(x="Filters (minimum periodicity)",y="Estimated Growth Impact\n(% per Degree)",col="Impact Model")
                    a=a+scale_color_manual(values=c("#d3818c","#7375a4")) + 
                    ggtitle("Impacts estimation using low-pass filters") #+ ylim(-10,5)
                    a
                    
                    #ggsave("Fig2.png",dpi=600)
                    
    ## 2.3. Perform MAIN simulation - WITH DRIFT ADJUSTMENT (end)

    ## 2.4. Additional simulations and figures (Used for Response to Reviewers)
        ## 2.4a. Perform simulation - WITH QUADRATIC TERM (start)

            time=350 #or 350 years
            basegr=0.01 #2% per year baseline growth rate
            start=100
            baseline=rep(basegr,time)
            coef=-0.05 #effect size - change in growth per degree warming()
            coef2=-0.05
            growthsd=0.005 #standard deviation of growth variability unexplained by temperature
            periods <- c(0,3,5,10,15)
            nsims=500
            sims=array(dim=c(nsims,length(periods),2))
            for(i in 1:nsims){
                randomtemp=Re(randomts(gtemp))[1:time]
                #randomtemp=Re(randomts(ustemp_lmr))[1:time]
                randomgrowth_g=basegr #growth impact model
                randomgrowth_l=basegr #levels impact model
                for(j in 2:time){
                #quadratic
                randomgrowth_g=c(randomgrowth_g, basegr+(randomtemp[j]*coef)+(coef2*randomtemp[j]^2)+rnorm(1,mean=0,sd=growthsd))
                randomgrowth_l=c(randomgrowth_l, basegr+(randomtemp[j]-randomtemp[j-1])*coef+coef2*(randomtemp[j]^2-randomtemp[j-1]^2)+rnorm(1,mean=0,sd=growthsd))
                #linear
                #randomgrowth_g=c(randomgrowth_g, basegr+(randomtemp[j]*coef)+rnorm(1,mean=0,sd=growthsd))
                #randomgrowth_l=c(randomgrowth_l, basegr+(randomtemp[j]-randomtemp[j-1])*coef+rnorm(1,mean=0,sd=growthsd))
                }
                dataset=data.frame(years=1:time,temp=randomtemp,g=randomgrowth_g,l=randomgrowth_l)
                dataset$templag =c(NA,dataset$temp[1:(dim(dataset)[1]-1)])
                for(k in 1:length(periods)){
                mod_g=lm(g~temp+templag,data=dataset)$coef[2]
                mod_l=lm(l~temp+templag,data=dataset)$coef[2]

                if (k==1){sims[i,k,]=c(mod_g,mod_l)
                    }else{
                    tempts <- pass.filt(randomtemp, W=periods[k], type="low", method="Butterworth")
                    growth <- randomgrowth_g
                    level <- randomgrowth_l
                    filt <- data.frame(temp = unclass(tempts), growth = unclass(growth),
                        level = unclass(level))
                    filt$templag=c(NA,filt$temp[1:(dim(filt)[1]-1)])
                    names(filt) <- c("temp","growth","levels","templag")
                    mod_gfilt=lm(growth~temp,data=filt)$coef[2]
                    mod_lfilt=lm(levels~temp,data=filt)$coef[2]
                    sims[i,k,]=c(mod_gfilt,mod_lfilt)
                }
                

                }
                if(i%%50==0) print(i)
            }
            # Comparing quadratic
                #Unfiltered Growth (start)
                    linreg <- lm(g~temp,data=dataset)
                    nonlinreg <- lm(g~temp+I(temp^2),data=dataset)
                    summary(linreg)            
                    summary(nonlinreg)

                    model <- nonlinreg
                    Sigma <- vcov(model)
                    coefT <- "temp"
                    start1 <- which(names(coef(model))==coefT)
                    end1 <- which(names(coef(model))==paste("I(",coefT,"^2)",sep=""))
                    sigma = Sigma[c(1:end1),c(1:end1)]
                    beta.hat <- coef(model)[c(1:end1)]
                    x <- seq(from=min(dataset$temp),to=max(dataset$temp), length=length(dataset$temp))
                    xmat <- cbind(1,x, x^2)
                    gestimated <- colSums(beta.hat*t(xmat)) 
                    ci12 <- gestimated + 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
                    ci22 <- gestimated -  1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))

                    glimpse(dataset)
                    dataset$gestimated2 <- gestimated
                    dataset$ci12 <- ci12
                    dataset$ci22 <- ci22
                    dataset$x <- xmat[,2]

                    model <- linreg
                    sigma <- vcov(model)
                    beta.hat <- coef(model)
                    x <- seq(from=min(dataset$temp),to=max(dataset$temp), length=length(dataset$temp))
                    xmat <- cbind(1,x)
                    gestimated <- colSums(beta.hat*t(xmat)) 
                    ci1 <- gestimated + 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
                    ci2 <- gestimated -  1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))

                    glimpse(dataset)
                    dataset$gestimated <- gestimated
                    dataset$ci1<- ci1
                    dataset$ci2 <- ci2

                    cols=c("#d3818c","#7375a4")
                    u_g <- ggplot(data=dataset, aes(x=temp,y=g))+
                    geom_point(alpha=0.5)+
                    geom_line(aes(x=x,y=gestimated),col=cols[1])+
                    geom_ribbon(aes(ymin=ci1,ymax=ci2,x=x),alpha=0.2,fill=cols[1])+
                    geom_line(aes(x=x,y=gestimated2),col=cols[2])+
                    geom_ribbon(aes(ymin=ci12,ymax=ci22,x=x),alpha=0.2,fill=cols[2])+
                    xlab("Temperature change")+
                    ylab("Growth") + 
                    ggtitle("Unfiltered Growth World")
                #Unfiltered Growth (end)
                
                #Unfiltered levels (start)
                    linreg <- lm(l~temp+templag,data=dataset)
                    nonlinreg <- lm(l~temp+I(temp^2)+templag,data=dataset)
                    summary(linreg)            
                    summary(nonlinreg)

                    model <- nonlinreg
                    Sigma <- vcov(model)
                    coefT <- "temp"
                    start1 <- which(names(coef(model))==coefT)
                    end1 <- which(names(coef(model))==paste("I(",coefT,"^2)",sep=""))
                    sigma = Sigma[c(1:end1),c(1:end1)]
                    beta.hat <- coef(model)[c(1:end1)]
                    x <- seq(from=min(dataset$temp),to=max(dataset$temp), length=length(dataset$temp))
                    xmat <- cbind(1,x, x^2)
                    gestimated <- colSums(beta.hat*t(xmat)) 
                    ci12 <- gestimated + 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
                    ci22 <- gestimated -  1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))

                    glimpse(dataset)
                    dataset$gestimated2 <- gestimated
                    dataset$ci12 <- ci12
                    dataset$ci22 <- ci22
                    dataset$x <- xmat[,2]

                    model <- linreg
                    sigma <- vcov(model)
                    start1 <- which(names(coef(model))==coefT)
                    end1 <- which(names(coef(model))==coefT)
                    sigma = Sigma[c(1:end1),c(1:end1)]
                    beta.hat <- coef(model)[c(1:end1)]
                    x <- seq(from=min(dataset$temp),to=max(dataset$temp), length=length(dataset$temp))
                    xmat <- cbind(1,x)
                    gestimated <- colSums(beta.hat*t(xmat)) 
                    ci1 <- gestimated + 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
                    ci2 <- gestimated -  1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))

                    glimpse(dataset)
                    dataset$gestimated <- gestimated
                    dataset$ci12 <- ci1
                    dataset$ci22 <- ci2

                    cols=c("#d3818c","#7375a4")
                    u_l <- ggplot(data=dataset, aes(x=temp,y=l))+
                    geom_point(alpha=0.5)+
                    geom_line(aes(x=x,y=gestimated),col=cols[1])+
                    geom_ribbon(aes(ymin=ci1,ymax=ci2,x=x),alpha=0.2,fill=cols[1])+
                    geom_line(aes(x=x,y=gestimated2),col=cols[2])+
                    geom_ribbon(aes(ymin=ci12,ymax=ci22,x=x),alpha=0.2,fill=cols[2])+
                    xlab("Temperature change")+
                    ylab("Growth") + 
                    ggtitle("Unfiltered Levels World")
                #Unfiltered levels (end)

                #15y-filtered Growth (start)
                    linreg <- lm(growth~temp,data=filt)
                    nonlinreg <- lm(growth~temp+I(temp^2),data=filt)
                    summary(linreg)            
                    summary(nonlinreg)

                    model <- nonlinreg
                    Sigma <- vcov(model)
                    coefT <- "temp"
                    start1 <- which(names(coef(model))==coefT)
                    end1 <- which(names(coef(model))==paste("I(",coefT,"^2)",sep=""))
                    sigma = Sigma[c(1:end1),c(1:end1)]
                    beta.hat <- coef(model)[c(1:end1)]
                    x <- seq(from=min(filt$temp),to=max(filt$temp), length=length(filt$temp))
                    xmat <- cbind(1,x, x^2)
                    gestimated <- colSums(beta.hat*t(xmat)) 
                    ci12 <- gestimated + 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
                    ci22 <- gestimated -  1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))

                    glimpse(filt)
                    filt$gestimated2 <- gestimated
                    filt$ci12 <- ci12
                    filt$ci22 <- ci22
                    filt$x <- xmat[,2]

                    model <- linreg
                    sigma <- vcov(model)
                    beta.hat <- coef(model)
                    x <- seq(from=min(filt$temp),to=max(filt$temp), length=length(filt$temp))
                    xmat <- cbind(1,x)
                    gestimated <- colSums(beta.hat*t(xmat)) 
                    ci1 <- gestimated + 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
                    ci2 <- gestimated -  1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))

                    glimpse(filt)
                    filt$gestimated <- gestimated
                    filt$ci1 <- ci1
                    filt$ci2 <- ci2

                    cols=c("#d3818c","#7375a4")
                    g_15 <- ggplot(data=filt, aes(x=temp,y=growth))+
                    geom_point(alpha=0.5)+
                    geom_line(aes(x=x,y=gestimated),col=cols[1])+
                    geom_ribbon(aes(ymin=ci1,ymax=ci2,x=x),alpha=0.2,fill=cols[1])+
                    geom_line(aes(x=x,y=gestimated2),col=cols[2])+
                    geom_ribbon(aes(ymin=ci12,ymax=ci22,x=x),alpha=0.2,fill=cols[2])+
                    xlab("Temperature change")+
                    ylab("Growth") + 
                    ggtitle("15-y filtered Growth World")
                #15y-filtered Growth (end)
                
                #15y-filtered levels (start)
                    linreg <- lm(levels~temp+templag,data=filt)
                    nonlinreg <- lm(levels~temp+I(temp^2)+templag,data=filt)
                    summary(linreg)            
                    summary(nonlinreg)

                    model <- nonlinreg
                    Sigma <- vcov(model)
                    coefT <- "temp"
                    start1 <- which(names(coef(model))==coefT)
                    end1 <- which(names(coef(model))==paste("I(",coefT,"^2)",sep=""))
                    sigma = Sigma[c(1:end1),c(1:end1)]
                    beta.hat <- coef(model)[c(1:end1)]
                    x <- seq(from=min(filt$temp),to=max(filt$temp), length=length(filt$temp))
                    xmat <- cbind(1,x, x^2)
                    gestimated <- colSums(beta.hat*t(xmat)) 
                    ci12 <- gestimated + 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
                    ci22 <- gestimated -  1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))

                    glimpse(filt)
                    filt$gestimated2 <- gestimated
                    filt$ci12 <- ci12
                    filt$ci22 <- ci22
                    filt$x <- xmat[,2]

                    model <- linreg
                    sigma <- vcov(model)
                    start1 <- which(names(coef(model))==coefT)
                    end1 <- which(names(coef(model))==coefT)
                    sigma = Sigma[c(1:end1),c(1:end1)]
                    beta.hat <- coef(model)[c(1:end1)]
                    x <- seq(from=min(filt$temp),to=max(filt$temp), length=length(filt$temp))
                    xmat <- cbind(1,x)
                    gestimated <- colSums(beta.hat*t(xmat)) 
                    ci1 <- gestimated + 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
                    ci2 <- gestimated -  1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))

                    glimpse(filt)
                    filt$gestimated <- gestimated
                    filt$ci1 <- ci1
                    filt$ci2 <- ci2

                    cols=c("#d3818c","#7375a4")
                    l_15 <- ggplot(data=filt, aes(x=temp,y=levels))+
                    geom_point(alpha=0.5)+
                    geom_line(aes(x=x,y=gestimated),col=cols[1])+
                    geom_ribbon(aes(ymin=ci1,ymax=ci2,x=x),alpha=0.2,fill=cols[1])+
                    geom_line(aes(x=x,y=gestimated2),col=cols[2])+
                    geom_ribbon(aes(ymin=ci12,ymax=ci22,x=x),alpha=0.2,fill=cols[2])+
                    xlab("Temperature change")+
                    ylab("Growth") + 
                    ggtitle("15-y filtered Levels World")
                #15y-filteredd levels (end)
                
                # Test 
                    simsfiltmean=apply(sims,MARGIN=c(2,3),FUN=mean)
                    simsfiltsd=apply(sims,MARGIN=c(2,3),FUN=sd)
                    colnames(simsfiltmean)=c("Growth","Level");rownames(simsfiltmean)=periods
                    simsfiltdat=cbind(melt(simsfiltmean),melt(simsfiltsd)[,3])
                    colnames(simsfiltdat)=c("periodsregationPeriod","ImpactType","MeanEffect","SDEffect")
                    theme_set(theme_bw(base_size = 20))
                    
                    simsfiltdat$periodsregationPeriod <- c("Unfiltered", "3 years", "5 years", "10 years", "15 years","Unfiltered", "3 years", "5 years", "10 years", "15 years")
                    simsfiltdat$periodsregationPeriod <- factor(simsfiltdat$periodsregationPeriod,
                        levels = c("Unfiltered", "3 years", "5 years", "10 years", "15 years"))
                    
                    glimpse(simsfiltdat$periodsregationPeriod)
                    
                    
                    a=ggplot(simsfiltdat,aes(x=factor(periodsregationPeriod),y=MeanEffect*100,col=ImpactType,group=ImpactType))+geom_line(lwd=1.25)
                    a=a+geom_point(size=4)+geom_errorbar(aes(ymin=(MeanEffect-1.96*SDEffect)*100,ymax=(MeanEffect+1.96*SDEffect)*100),width=0.1,lwd=1.25)
                    a=a+geom_hline(yintercept=0,lwd=1.5,lty=3)
                    a=a+labs(x="Filters (minimum periodicity)",y="Estimated Growth Impact\n(% per Degree)",col="Impact Model")
                    a=a+scale_color_manual(values=c("#d3818c","#7375a4")) + ggtitle("Impacts estimation using low pass filters") #+ ylim(-10,5)
                    a
                # Comparing quadratic (end)
                
                ggarrange(ggarrange(u_l,l_15,u_g,g_15),a,heights = c(2, 0.7),nrow=2,ncol=1)
                #ggsave("FigResponse_Quadratic.png",dpi=600)
        ## 2.4a. Perform simulation - WITH QUADRATIC TERM (end)

        ## 2.4b. Perform simulation - Combination Gammas and Betas (start)


            time=350 #or 350 years
            basegr=0.01 #2% per year baseline growth rate
            start=100
            baseline=rep(basegr,time)
            coef_gamma=-0.05 #effect size - change in growth per degree warming()
            coef_beta=-0.05
            coef_beta2=0.05
            coef_beta3=0.07
            growthsd=0.005 #standard deviation of growth variability unexplained by temperature
            periods <- c(0,3,5,10,15)
            nsims=500
            sims=array(dim=c(nsims,length(periods),2))
            sims_adjusted=array(dim=c(nsims,length(periods),3))
            for(i in 1:nsims){
                randomtemp=Re(randomts(gtemp))[1:time]
                #randomtemp=Re(randomts(ustemp_lmr))[1:time]
                randomgrowth_g=basegr #growth impact model
                randomgrowth_l=basegr #levels impact model
                randomgrowth_l2=basegr #levels impact model
                for(j in 2:time){
                randomgrowth_g=c(randomgrowth_g, basegr+(randomtemp[j]*coef_gamma)+(randomtemp[j]-randomtemp[j-1])*coef_beta+rnorm(1,mean=0,sd=growthsd))
                randomgrowth_l=c(randomgrowth_l, basegr+(randomtemp[j]*coef_gamma)+(randomtemp[j]-randomtemp[j-1])*coef_beta2+rnorm(1,mean=0,sd=growthsd))
                randomgrowth_l2=c(randomgrowth_l2, basegr+(randomtemp[j]*coef_gamma)+(randomtemp[j]-randomtemp[j-1])*coef_beta3+rnorm(1,mean=0,sd=growthsd))
                
                }
                dataset=data.frame(years=1:time,temp=randomtemp,g=randomgrowth_g,l=randomgrowth_l,l2=randomgrowth_l2)
                dataset$templag =c(NA,dataset$temp[1:(dim(dataset)[1]-1)])
                for(k in 1:length(periods)){
                mod_g=lm(g~temp+templag,data=dataset)$coef[2]
                mod_l=lm(l~temp+templag,data=dataset)$coef[2]
                mod_l2=lm(l2~temp+templag,data=dataset)$coef[2]

                if (k==1){
                        sims[i,k,]=c(mod_g,mod_l)
                        sims_adjusted[i,k,]=c(mod_g,mod_l,mod_l2)
                        
                    }else{
                    x <- seq(1:length(randomtemp))
                    randomtemp <- lm(randomtemp~x)$residuals
                    tempts <- pass.filt(randomtemp, W=periods[k], type="low", method="Butterworth")
                    
                    ratio <- median(randomtemp/tempts)
                    #ratio <- median(randomtemp / tempts)
                    #ratio <- median( abs(randomtemp) / abs(tempts))
                    growth <- randomgrowth_g
                    level <- randomgrowth_l
                    level2 <- randomgrowth_l2
                    filt <- data.frame(temp = unclass(tempts), growth = unclass(growth),
                        level = unclass(level),level2 = unclass(level2))
                    filt$templag=c(NA,filt$temp[1:(dim(filt)[1]-1)])
                    names(filt) <- c("temp","growth","levels","levels2","templag")
                    mod_gfilt=lm(growth~temp,data=filt)$coef[2]
                    mod_lfilt=lm(levels~temp,data=filt)$coef[2]
                    mod_lfilt2=lm(levels2~temp,data=filt)$coef[2]
                    sims[i,k,]=c(mod_gfilt,mod_lfilt)
                    sims_adjusted[i,k,]=c(mod_gfilt/ratio,mod_lfilt/ratio,mod_lfilt2/ratio)
                }
                

                }
                if(i%%50==0) print(i)
            }
            # Comparing quadratic
          
                # Test 
                    simsfiltmean=apply(sims_adjusted,MARGIN=c(2,3),FUN=mean)
                    simsfiltsd=apply(sims_adjusted,MARGIN=c(2,3),FUN=sd)
                    colnames(simsfiltmean)=c("Attenuating","Intensifying","Changing Sign");rownames(simsfiltmean)=periods
                    simsfiltdat=cbind(melt(simsfiltmean),melt(simsfiltsd)[,3])
                    colnames(simsfiltdat)=c("periodsregationPeriod","ImpactType","MeanEffect","SDEffect")
                    theme_set(theme_bw(base_size = 20))
                    
                    simsfiltdat$periodsregationPeriod <- c("Unfiltered", "3 years", "5 years", "10 years", "15 years","Unfiltered", "3 years", "5 years", "10 years", "15 years","Unfiltered", "3 years", "5 years", "10 years", "15 years")
                    simsfiltdat$periodsregationPeriod <- factor(simsfiltdat$periodsregationPeriod,
                        levels = c("Unfiltered", "3 years", "5 years", "10 years", "15 years"))
                    simsfiltdat_adjusted <- simsfiltdat

                    # simsfiltmean=apply(sims,MARGIN=c(2,3),FUN=mean)
                    # simsfiltsd=apply(sims,MARGIN=c(2,3),FUN=sd)
                    # colnames(simsfiltmean)=c("Growth","Level");rownames(simsfiltmean)=periods
                    # simsfiltdat=cbind(melt(simsfiltmean),melt(simsfiltsd)[,3])
                    # colnames(simsfiltdat)=c("periodsregationPeriod","ImpactType","MeanEffect","SDEffect")
                    # theme_set(theme_bw(base_size = 20))
                    
                    # #simsfiltdat$periodsregationPeriod <- c("Unfiltered", "3 years", "5 years", "10 years", "15 years","Unfiltered", "3 years", "5 years", "10 years", "15 years")
                    # simsfiltdat$periodsregationPeriod <- c("Unfiltered", "10 years", "20 years", "50 years", "100 years","Unfiltered", "10 years", "20 years", "50 years", "100 years")
                    # simsfiltdat$periodsregationPeriod <- factor(simsfiltdat$periodsregationPeriod,
                    #     levels = c("Unfiltered", "10 years", "20 years", "50 years", "100 years"))
                    
                    
                    a=ggplot(simsfiltdat_adjusted,aes(x=factor(periodsregationPeriod),y=MeanEffect*100,col=ImpactType,group=ImpactType))+geom_line(lwd=1.25)
                    a=a+geom_point(size=4)+geom_errorbar(aes(ymin=(MeanEffect-1.96*SDEffect)*100,ymax=(MeanEffect+1.96*SDEffect)*100),width=0.1,lwd=1.25)
                    a=a+geom_hline(yintercept=0,lwd=1.5,lty=3)
                    a=a+labs(x="Filters (minimum periodicity)",y="Estimated Growth Impact\n(% per Degree)",col="Impact Model")
                    a=a+scale_color_manual(values=c("#7fc97f","#beaed4","#fdc086"),
                        labels=c(expression(paste(beta,"=-0.05; ",gamma,"=-0.05")),expression(paste(beta,"=0.05; ",gamma,"=-0.05")),expression(paste(beta,"=0.07; ",gamma,"=-0.05")))) + 
                    ggtitle("Combined effects") #+ ylim(-10,5)
                    a


                    b=ggplot(simsfiltdat,aes(x=factor(periodsregationPeriod),y=MeanEffect*100,col=ImpactType,group=ImpactType))+geom_line(lwd=1.25)
                    b=b+geom_point(size=4)+geom_errorbar(aes(ymin=(MeanEffect-1.96*SDEffect)*100,ymax=(MeanEffect+1.96*SDEffect)*100),width=0.1,lwd=1.25)
                    b=b+geom_hline(yintercept=0,lwd=1.5,lty=3)
                    b=b+labs(x="Filters (minimum periodicity)",y="Estimated Growth Impact\n(% per Degree)",col="Impact Model")
                    b=b+scale_color_manual(values=c("#7fc97f","#beaed4","#fdc086")) + ggtitle("Impacts estimation using low pass filters") #+ ylim(-10,5) 
                    b
                # Comparing quadratic (end)
                
                ggarrange(b,a,nrow=2,ncol=1)
        ## 2.4b. Perform simulation - Combination Gammas and Betas (end)

    ## 2.4. Additional simulations and figures (Used for Response to Reviewers)


## 2. Simulation - Figures 1 and 2 (end)

## 3. Empirical analysis - Figure 3 and Table 1 (start)

    ## 3.1. Country-level regressions (start)
        dataset <- c("wb","barro","mad") 
        datasetweather <- c("UDel")
        periods <- c(0,3,5,10,15)
        countries <- unique(wb$countrycode)
        #periods <- c(0,10,20,25,27) Uncomment to get Supp Fig 2
        fullmods_filter=array(dim=c(length(countries),2,length(periods),length(dataset),length(datasetweather)))
        fullmods_filter_var=array(dim=c(length(countries),2,length(periods),length(dataset),length(datasetweather)))
        fullmods_filter_ar=array(dim=c(length(countries),2,length(periods),length(dataset),length(datasetweather)))
        fullmods_filter_p=array(dim=c(length(countries),2,length(periods),length(dataset),length(datasetweather)))
        panel_data <- data.frame(years = integer(),temp = double(), growth = double(), 
            preci = double(), countrycode = factor(), climdata = character(),
            econdata = character(), filter = character(), meant = double(), meanp=double(), growth_wdi=double())
        for (mm in 1:length(datasetweather)){
            tempname <- paste(datasetweather[mm],"_pop_temp",sep="")
            preciname <- paste(datasetweather[mm],"_pop_preci",sep="")
            for (jj in (1:length(dataset))){
                DATA <- get(dataset[jj])
                countries=unique(factor(DATA$countrycode))
                #fullmods=array(dim=c(length(countries),2,length(periods)))
                for(k in 1:length(periods)){
                for(i in 1:length(countries)){
                    dat=DATA[which(DATA$countrycode==countries[i]),
                        which(colnames(DATA)%in%c("countrycode","year","growth",tempname,preciname))] 
                    dat <- dat[is.finite(dat$growth),]

                        gonext <- 0
                    for (nn in (1:(dim(dat)[2]))){
                        if(sum(is.na(dat[,nn]))==dim(dat)[1]){
                            print(paste(countries[i],i,"missing complete column"))
                            gonext <- 1 }
                        
                    }
                    if(gonext==1){
                        gonext <- 0
                        next}
                    if(sum(complete.cases(dat))<2){
                            print(paste(countries[i],i,"not complete cases"))
                            next
                    }
                    dat <- dat %>% rename("temp"=tempname,"preci"=preciname)
                    dat <- dat[complete.cases(dat),]
                    if(max(dat$year)-min(dat$year) != (dim(dat)[1]-1)){
                        print(paste("Missing observation. Gap. To Do",dataset[jj]))
                        next

                    }
                    #dat <- dat[-c(1),] #first observation is missing
                    if(dim(dat)[1]<5){next}
                    mt <- lm(temp~year+I(year^2), data = dat)
                    meanT <- mean(dat$temp, na.rm = TRUE)
                    t <- resid(mt)    
                    mp <- lm(preci~year+I(year^2), data = dat)
                    p <- resid(mp)
                    meanP <- mean(dat$preci, na.rm = TRUE)
                    mg <- lm(growth~year+I(year^2), data = dat)
                    g <- resid(mg)
                    growth_wdi <- dat$growth
                    if(length(t)<(2*(periods[k]))){
                            print(paste(countries[i],i,"not enough data for this filter","length=",length(t),"periods=",periods[k]))
                            next}
                    if(k==1){
                        tempts <- t
                        precits <- p
                        ratio <- 1
                    } else{
                        tempts <- pass.filt(t, W=periods[k], type="low", method="Butterworth")
                        precits <- pass.filt(p, W=periods[k], type="low", method="Butterworth")
                        x <- seq(1:length(tempts))
                        demeaned_t <- t
                        demeaned_tempts <- tempts
                        ratio <- abs(median(demeaned_t/demeaned_tempts))
                        if(median(demeaned_t/demeaned_tempts)<0){
                            print(paste("ratio was negative for ", countries[i], "at period",toString(periods[k])))
                        }
                     
                    }
                     
                    temp <- data.frame(year = dat$year, temp= unclass(tempts))
                    preci <- data.frame(year = dat$year, preci= unclass(precits))
                    growth <- data.frame(year = dat$year, growth = unclass(g))
                    #growth_wdi <- data.frame(year = dat$year, growth = unclass(dat$growth))
                    filterdata <- merge(temp,growth, by = "year")
                    filterdata <- merge(filterdata,preci, by = "year")
                    #filterdata <- merge(filterdata,growth_wdi, by = "year")
                    names(filterdata) <- c("years","temp","growth","preci")

                    if(k==1){
                        filterdata2 <- filterdata
                        filterdata2$templag =c(NA,filterdata2$temp[1:(dim(filterdata2)[1]-1)])
                        mod_gfilterdata=lm(growth~temp+preci+templag,data=filterdata2)
                    } else{
                        mod_gfilterdata=lm(growth~temp+preci,data=filterdata)
                    }
                    dw <- dwtest(mod_gfilterdata) #test for autocorrelation
                    if (dw$p.value[1]<0.1){
                        vcov_ <- sandwich::vcovHAC(mod_gfilterdata)
                        hacse <- unname(sqrt(diag(vcov_)))[2] #s.e. corrected for autocorrelation
                    } else {
                        hacse <- summary(mod_gfilterdata)$coefficients[2,2]
                        vcov_ <-vcov(mod_gfilterdata)
                    } #We could try Newey West
                    
                    #filterdata$years <- dat$year
                    filterdata$countrycode <- rep(dat$countrycode[1],dim(filterdata)[1])
                    filterdata$climdata <- rep(datasetweather[mm],dim(filterdata)[1])
                    filterdata$econdata <- rep(dataset[jj],dim(filterdata)[1])
                    filterdata$filter <- rep(paste(periods[k],sep="-"),dim(filterdata)[1])
                    filterdata$meant <- rep(meanT,dim(filterdata)[1])
                    filterdata$meanp <- rep(meanP,dim(filterdata)[1])
                    filterdata$growth_wdi <- dat$growth 
                    panel_data <- bind_rows(panel_data,filterdata)
                    fullmods_filter[i,,k,jj,mm]=c(summary(mod_gfilterdata)$coefficients[2,1]/ratio,hacse/ratio)
                    fullmods_filter_var[i,,k,jj,mm]=c(summary(mod_gfilterdata)$coefficients[2,1]/ratio,vcov_[2,2]/(ratio^2))
                    fullmods_filter_ar[i,,k,jj,mm]=c(dw$statistic[1],dw$p.value[1])
                    fullmods_filter_p[i,,k,jj,mm]=summary(mod_gfilterdata)$coefficients[3,1:2] #Not HAC s.e. If usd should be divided by the ratio of precipietation timeseries
                    }
                    
                
                }
            }
        }
        original_fullmods_filter <- fullmods_filter
        original_fullmods_filter_var <- fullmods_filter_var
        glimpse(fullmods_filter)
        ranges <- paste(periods,sep='-')
        countries <- unique(factor(wb$countrycode))
        countries_wb <- unlist(lapply(countries, as.character))
        dimnames(fullmods_filter)=list(countries,c("Estimate","StandardError"),ranges,dataset,datasetweather)
        fullmods_filterm <- melt(fullmods_filter)
        glimpse(fullmods_filterm)
        #
        names(fullmods_filterm) <- c("countrycode","variable","frequencies","econdata","climdata","value")
        new_fft <- fullmods_filterm[order(fullmods_filterm[,2], fullmods_filterm[,3] ),] # Sort by column index [1] then [3]
        glimpse(new_fft)
        estimates <- new_fft$value[which(new_fft$variable=="Estimate")]        
        ses <- new_fft$value[which(new_fft$variable=="StandardError")]
        new_fft <- new_fft[c(1:length(estimates)),-c(2,6)]    
        new_fft$Estimate <- estimates
        new_fft$StandardError <- ses   
        glimpse(new_fft) 
        new_fft <- new_fft[order(new_fft[,1]),] # Sort by column index [1] then [3]     
        fullmods_filter <- new_fft
        #
        # countriesbarro <- unique(factor(barro$countrycode))
        # countries_barro <- unlist(lapply(countriesbarro , as.character))
        # countries_barro_extended <- c(countries_barro,countries_wb[(length(countries_barro)+1):length(countries_wb)*NA])
        # countriesmad <- unique(factor(mad$countrycode))
        # countries_mad <- unlist(lapply(countriesmad , as.character))
        # countries_mad_extended <- c(countries_mad,countries_wb[(length(countries_mad)+1):length(countries_wb)*NA])
        # countries_wbrep <- rep(countries_wb,2*5)
        # countries_barrorep <- rep(countries_barro_extended,2*5)
        # countries_madrep <- rep(countries_mad_extended,2*5) #2*5
        # names1 <- c(countries_wbrep,countries_barrorep,countries_madrep)
        # #names2 <- c(names1,names1)
        # names2 <- names1 #only usiong UDel and not LMR
        # fullmods_filterm$Var1 <- names2
        # fullmods_filter=dcast(fullmods_filterm,Var1+Var3+Var4+Var5~Var2, fun = mean)
        # fullmods_filter$Estimate[which(is.infinite(fullmods_filter$StandardError))]=NA
        # names(fullmods_filter) <- c("countrycode","frequencies","econdata","climdata","Estimate","StandardError")
        
        ranges <- paste(periods,sep='-')
        countries <- unique(factor(wb$countrycode))
        countries_wb <- unlist(lapply(countries, as.character))
        dimnames(fullmods_filter_var)=list(countries,c("Estimate","StandardError"),ranges,dataset,datasetweather)
        fullmods_filter_varm <- melt(fullmods_filter_var)
        glimpse(fullmods_filter_var)
        #
        names(fullmods_filter_varm) <- c("countrycode","variable","frequencies","econdata","climdata","value")
        new_fft_var <- fullmods_filter_varm[order(fullmods_filter_varm[,2], fullmods_filter_varm[,3] ),] # Sort by column index [1] then [3]
        estimates <- new_fft_var$value[which(new_fft_var$variable=="Estimate")]        
        ses <- new_fft_var$value[which(new_fft_var$variable=="StandardError")]
        new_fft_var <- new_fft_var[c(1:length(estimates)),-c(2,6)]    
        new_fft_var$Estimate <- estimates
        new_fft_var$Variance <- ses   
        glimpse(new_fft_var) 
        new_fft_var  <- new_fft_var[order(new_fft_var[,1]),] # Sort by column index [1] then [3]     
        fullmods_filter_var <- new_fft_var    


        # dimnames(fullmods_filter_var)=list(countries,c("Estimate","Var"),ranges,dataset,datasetweather)
        # fullmods_filterm <- melt(fullmods_filter_var)
        # fullmods_filterm$Var1 <- names2
        # fullmods_filter_v=dcast(fullmods_filterm,Var1+Var3+Var4+Var5~Var2, fun = mean)
        # fullmods_filter_v$Estimate[which(is.infinite(fullmods_filter_v$Var))]=NA
        # names(fullmods_filter_v) <- c("countrycode","frequencies","econdata","climdata","Estimate","Variance")

        # dimnames(fullmods_filter_ar)=list(countries,c("dw_test","dw_pvalue"),ranges,dataset,datasetweather)
        # fullmods_filterm <- melt(fullmods_filter_ar)
        # fullmods_filterm$Var1 <- names2
        # fullmods_filter_ar=dcast(fullmods_filterm,Var1+Var3+Var4+Var5~Var2, fun = mean)
        # fullmods_filter_ar$dw_test[which(is.infinite(fullmods_filter_ar$Var))]=NA
        # names(fullmods_filter_ar) <- c("countrycode","frequencies","econdata","climdata","dw_test","dw_pvalue")

        # Is it autocorrelated?
            
            # aa <- ggplot(data=fullmods_filter_ar, aes(dw_pvalue,fill=frequencies,group=frequencies))+theme_bw()+  # `data` must be a data frame, or other object coercible by `fortify()`, not a numeric vector.
            # geom_histogram(breaks=c(0,0.001,0.005,0.01,1))+ scale_x_continuous(trans='log',breaks=c(0,0.001,0.005,0.01,0.2,0.5,1))+
            # ggtitle("P-values of the DW test for autocorrelation")
            #
            #
            # f_ar <- fullmods_filter_ar[fullmods_filter_ar$econdata=="wb",]
            # f_ar <- f_ar[f_ar$climdata=="UDel",]
            # f_ar$dw_pvalue_ad <- p.adjust(f_ar$dw_pvalue,"bonferroni")
            #
            # aa2 <- ggplot(data=f_ar, aes(dw_pvalue,fill=frequencies,group=frequencies))+theme_bw()+
            # geom_histogram(breaks=c(2e-5,0.05,1))+ scale_x_continuous(trans='log',breaks=c(2e-5,0.05,1))+
            # ggtitle("P-values of the DW test for autocorrelation")
            # aa2
            #
            # aa3 <- ggplot(data=f_ar, aes(dw_pvalue_ad,fill=frequencies,group=frequencies))+theme_bw()+
            # geom_histogram(breaks=c(2e-5,0.05,1))+ scale_x_continuous(trans='log',breaks=c(2e-5,0.05,1))+
            # ggtitle("P-values of the DW test for autocorrelation\nBonferroni adjusted")
            # aa3
            #
            # ggarrange(aa2,aa3,common.legend=TRUE,legend="bottom")
            # #ggsave("dwtest_autocorrelation.png",dpi=600)
            #
            # #ab <- ggplot(data=fullmods_filter_ar, aes(dw_test,fill=frequencies,group=frequencies))+theme_bw()+
            # #geom_histogram()+
            # #ggtitle("Autocorrelation Estimate")
            #
            # #ggarrange(ab, aa,common.legend=TRUE,ncol=1,nrow=2,legend="bottom")


        #Is it autocorrelated?
    ## 3.1. Country-level regressions (end)
        
    ## 3.2 Panel Regression (start)
        
        coefT <- "temp"
        coefP <- "preci"
        panel_data$temp <- panel_data$temp + panel_data$meant
        panel_data$preci <- panel_data$preci + panel_data$meanp
        
        panel_results <- data.frame(temp=double(),dgdt=double(),ci1=double(),ci2=double(),filter=double())
        
        for(period in periods){
        p1 <- panel_data[which(panel_data$climdata=="UDel" & panel_data$econdata=="wb"  & panel_data$filter==toString(period)), ] #take only unfiltered data
        felm_panel1 <- felm(growth_wdi ~ temp+I(temp^2)+preci+I(preci^2)+countrycode:years+countrycode:I(years^2)|countrycode+years|0|countrycode, data =p1)
            #felm_panel1 <- felm(growth ~ temp+I(temp^2)+preci+I(preci^2)|countrycode+years|0|countrycode, data =p1)
            assign(paste("felm_panel",toString(period),sep=""),felm_panel1)
        # Getting marginal effect of an additional degree (start)
            #summary(felm_panel1)
            model <- felm_panel1
            Sigma <- vcov(model)
            start1 <- which(names(coef(model))==coefT)
            end1 <- which(names(coef(model))==paste("I(",coefT,"^2)",sep=""))
            sigma = Sigma[c(start1,end1),c(start1,end1)]
            beta.hat <- coef(model)[c(start1,end1)]
            x <- seq(from=min(p1$temp),to=max(p1$temp))
            xmat <- cbind(1,2*x)
            dgdt <- colSums(beta.hat*t(xmat)) 
            ci1 <- dgdt + 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
            ci2 <- dgdt- 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
            
             panel_results <-  rbind(panel_results,data.frame(temp=x,dgdt=dgdt,ci1=ci1,ci2=ci2,filter=rep(period,length(x))))

            
        }
            #stargazer(felm_panel0,felm_panel3,felm_panel5,felm_panel10,felm_panel15,type="html", out="felm_panel.html")
            stargazer(felm_panel0,felm_panel3,felm_panel5,felm_panel10,felm_panel15,type="text")

            ggplot(panel_results,aes(x=temp,y=dgdt,ymin=ci2,ymax=ci1,fill=factor(filter)))+geom_line(aes(color=factor(filter)))+
            scale_fill_discrete(name="Filters")+scale_color_discrete(name="Filters")+
            theme_bw()+xlab("Temperature (C)") + ylab("Change in Growth")+
            geom_hline(yintercept=0, linetype="dashed", color = "black")+
            geom_ribbon(alpha=0.3)+facet_wrap(~filter)
            #ggsave("Panel_reg_marginal.png",dpi=600)
        # Getting marginal effect of an additional degree (start)
        

    ## 3.2 Panel Regression (end)

    ## 3.3. categorizing - Plotting Figure 3  (start)

        # Categorizing statistically different estimates (start)
         
            fmod_fft <- new_fft
            fmod_fft$Variance <- fullmods_filter_var$Variance
            
            #fmod_fft <- fullmods_filter
            fmod_fft <- fmod_fft[fmod_fft$econdata=="wb",]
            fmod_fft <- fmod_fft[fmod_fft$climdata=="UDel",]
            fmod_fft$high95 <- fmod_fft$Estimate + fmod_fft$StandardError*1.64
            fmod_fft$low95 <- fmod_fft$Estimate - fmod_fft$StandardError*1.64
            numcountries <- length(unique(fmod_fft$countrycode))
            #fmod_fft$absestimate <- abs(fmod_fft$Estimate)
            fmod_fft$lowsignificant <- 0
            fmod_fft$unfilteredsignificant <- 0
            fmod_fft$signlofreq <- 0
                            
            uncertain <- c(which(fmod_fft$high95>0 & fmod_fft$low95 <0))
            `%notin%` <- Negate(`%in%`)
            fmod_fft$sign <- rep(0,dim(fmod_fft)[1])
            
            numcountries <- length(unique(fmod_fft$countrycode))
            glimpse(fmod_fft)
            for (i in 1:numcountries){
                if(is.null(fmod_fft$Estimate[1+5*(i-1)] )){next}
                if(!is.na(fmod_fft$Estimate[5+5*(i-1)])){
                    lastfreq <- 5
                    }else if(!is.na(fmod_fft$Estimate[4+5*(i-1)])){
                        lastfreq <- 4
                    }else if(!is.na(fmod_fft$Estimate[3+5*(i-1)])){
                        lastfreq <- 3
                    }else if(!is.na(fmod_fft$Estimate[2+5*(i-1)])){
                        lastfreq <- 2
                    }else{next}
                if(fmod_fft$Estimate[lastfreq+5*(i-1)] >0){fmod_fft$signlofreq[(1+5*(i-1)):(5+5*(i-1))]=1}
                m <- fmod_fft$Estimate[1+5*(i-1)] * fmod_fft$Estimate[lastfreq+5*(i-1)]
                if (is.na(m)){next}
                if(m>0 ){
                    fmod_fft$sign[(1+5*(i-1)):(5+5*(i-1))]=1
                    }
            }

                converging_stat <- 0
                converging <- 0
                intensifying <- 0
                intensifying_stat <- 0
                for (i in 1:numcountries){
                    if(!is.na(fmod_fft$Estimate[5+5*(i-1)])){
                            lastfreq <- 5
                            }else if(!is.na(fmod_fft$Estimate[4+5*(i-1)])){
                                lastfreq <- 4
                            }else if(!is.na(fmod_fft$Estimate[3+5*(i-1)])){
                                lastfreq <- 3
                            }else if(!is.na(fmod_fft$Estimate[2+5*(i-1)])){
                                lastfreq <- 2
                            }else{next}
                            
                            if(is.na(fmod_fft$Estimate[1+5*(i-1)])){next}
                    theta <- abs(fmod_fft$Estimate[1+5*(i-1)]) - abs(fmod_fft$Estimate[lastfreq+5*(i-1)])
                    var <- fmod_fft$Variance[lastfreq+5*(i-1)] + fmod_fft$Variance[1+5*(i-1)]
                    conf95 <-  (var^0.5)*1.65 #one-tail  95%
                    
                    if(abs(fmod_fft$Estimate[lastfreq+5*(i-1)]) > abs(fmod_fft$Estimate[1+5*(i-1)])){
                        intensifying <- c(intensifying,(1+5*(i-1)):(5+5*(i-1))) 
                        if(theta+conf95 < 0 ){
                            intensifying_stat <- c(intensifying_stat,(1+5*(i-1)):(5+5*(i-1))) 
                        }
                    }
                    
                    if(abs(fmod_fft$Estimate[lastfreq+5*(i-1)]) <= abs(fmod_fft$Estimate[1+5*(i-1)])){
                            converging <- c(converging,(1+5*(i-1)):(5+5*(i-1))) 
                            if(theta-conf95 >= 0 ){
                            converging_stat <- c(converging_stat,(1+5*(i-1)):(5+5*(i-1))) 
                        }
                    }
                    
                    
                        
                        if((lastfreq+5*(i-1)) %notin% uncertain){
                        fmod_fft$lowsignificant[(1+5*(i-1)):(5+5*(i-1))] <- 1     
                                       }
                        if((1+5*(i-1)) %notin% uncertain){
                            fmod_fft$unfilteredsignificant[(1+5*(i-1)):(5+5*(i-1))] <- 1     
                                       }
                    }
                    
                    
                filt_names <- c("Unfiltered","3 years","5 years", "10 years", "15 years")
                
                fmod_fft$filters <- rep(filt_names,length(unique(wb$countrycode)))
                fmod_fft$filters <- factor(fmod_fft$filters, levels = filt_names)            
                fmod_fft$category <- "Undefined"
                fmod_fft$category[converging] <- "converging"
                fmod_fft$category[intensifying] <- "intensifying"
                fmod_fft$category[converging_stat] <- "converging_stat"
                fmod_fft$category[intensifying_stat] <- "intensifying_stat"
                table(fmod_fft$category)/5
                
                fmod_fft$Category <- "Undefined"
                fmod_fft$Category[converging] <- "converging"
                fmod_fft$Category[intensifying] <- "intensifying"
                table(fmod_fft$Category)/5
                
                fmod_fft$CoeffDifferents <- 0
                fmod_fft$CoeffDifferents[intensifying_stat] <- 1
                fmod_fft$CoeffDifferents[converging_stat] <- 1
                table(fmod_fft$CoeffDifferents)/5
                
                table(fmod_fft$lowsignificant)/5

                fmod_fft_na<-fmod_fft[fmod_fft$category=="Undefined" | 
                    fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate<quantile(fmod_fft$Estimate, na.rm=TRUE,0.01))] |
                    fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate>quantile(fmod_fft$Estimate, na.rm=TRUE,0.99))],]
                fmod_fft<-fmod_fft[fmod_fft$category!="Undefined",]
                library('rnaturalearth')
                world <- ne_countries(scale = "medium", returnclass = "sf")

                
                fmod_fft<- fmod_fft[fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate>quantile(fmod_fft$Estimate, na.rm=TRUE,0.99))],]
                fmod_fft<- fmod_fft[fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate<quantile(fmod_fft$Estimate, na.rm=TRUE,0.01))],]
               
                fmod_fft$iso_a3 <- fmod_fft$countrycode
                fmod_fft_map <- merge(world,fmod_fft,by="iso_a3")
                fmod_fft_na$iso_a3 <- fmod_fft_na$countrycode
                fmod_fft_map_na <- merge(world,fmod_fft_na,by="iso_a3")
                

                #factor(fmod_fft$category)
                fmod_fft$category <- factor(fmod_fft$category, levels = c("intensifying_stat", "intensifying", "converging","converging_stat"))
                df1 <- fmod_fft[which(fmod_fft$econdata=="wb" &fmod_fft$climdata=="UDel" & fmod_fft$filter=="Unfiltered"),]
                pal1 <- sequential_hcl(7, palette = "Oranges")
                pal2 <- sequential_hcl(7, palette = "Teal")
                df1$unfilteredsignificant <- factor(df1$unfilteredsignificant)
                df1$lowsignificant <- factor(df1$lowsignificant)
                    
                hist_bf <- ggplot() + theme_bw() + 
                    geom_histogram(data=df1, na.rm= TRUE, mapping = aes(x=lowsignificant, fill=factor(category)), 
                    #stat='count')+ scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],pal2[1]), 
                    #name="",labels=c("*Intensifying","Intensifying","Converging","*Converging"))+
                    stat='count')+ scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3]), 
                    name="",labels=c("*Intensifying","Intensifying","Converging"))+ #No converging statistically significant found
                    scale_x_discrete(name="",
                        labels=c(expression(paste("|",theta[f],"| =0",sep="")),
                        expression(paste("|",theta[f],"| >0",sep=""))))+
                    scale_y_continuous(name="Number of countries")+theme(legend.position="none")

                df2 <- fmod_fft[which(fmod_fft$econdata=="wb" &fmod_fft$climdata=="UDel" & fmod_fft$filter=="Unfiltered" & fmod_fft$unfilteredsignificant==1),]
                df2$unfilteredsignificant <- factor(df2$unfilteredsignificant)
                glimpse(df2)     
                hist_bu <- ggplot() + theme_bw()+
                    geom_histogram(data=df2, na.rm= TRUE, mapping = aes(x=unfilteredsignificant, fill=factor(category)),stat='count')+ 
                    #scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],pal2[1]), 
                    #name="",labels=c("*Intensifying","Intensifying","Converging","*Converging"))+
                    scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3]), 
                    name="",labels=c("*Intensifying","Intensifying","Converging"))+ #No converging statistically significant found
                    scale_x_discrete(name="",labels=c(expression(paste("|",theta[U],"| >0",sep=""))))+
                    scale_y_continuous(name="Number of countries") + theme(legend.position="none")
                 
                 plot_fg <- ggplot(data=fmod_fft,aes(x=filters,y=Estimate*100, group = countrycode,color=factor(category)))+
                    geom_line()+
                    scale_colour_manual(name="Categories", values=c(pal1[1],pal1[3],pal2[3],pal2[1]),
                    labels=c("Converging","Intensifying","Statistically Intensifying","Statistically Converging")) +
                    theme_bw() + xlab("Minimum Periodicity after Filtering")+
                    geom_hline(yintercept=0,lty=2)+
                    #geom_dl(data=fg,aes(label = countrycode), method = list(dl.combine("last.points")), cex = 0.9)+
                    ylab("Estimated Effect of \n 1 Degree Warming on Growth (pp)") +
                    theme(axis.line = element_line(colour = "black"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(),legend.position="none") +
                    #panel.background = element_blank()) +
                    ggtitle("")
                    plot_fg      
                
                
                plot_categoriescolors <- ggplot() + geom_histogram(data=data.frame(x=c("*Intensifying","Intensifying","Converging","*Converging")), na.rm= TRUE, mapping = aes(x=x, fill=factor(x)),stat='count')+ 
                scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],pal2[1]),name="",labels=c("*Intensifying","Intensifying","Converging","*Converging"))
                leg_cat <- get_legend(plot_categoriescolors)                    
                
                fig_3 <- ggarrange(plot_fg,ggarrange(hist_bu,leg_cat,hist_bf,ncol=3,nrow=1),ncol=1,nrow=2)    
                fig_3
                #ggsave("Figures/Fig3.png",dpi=600)
                
                map_categories <- ggplot(data=fmod_fft_map) +
                            geom_sf(data = fmod_fft_map_na,aes(fill=NA))+
                            theme_minimal()+
                            geom_sf(aes(fill = category))+
                            scale_fill_manual(name = "",
                                            #labels=c("Converging","Intensifying","*Intensifying"),
                                            values=c(pal2[3],pal1[3],pal1[1]))+
                            theme(legend.position="none")+
                            ggtitle("Behavior of coefficients")
                    ggarrange(map_categories,leg_cat,ncol=2,nrow=1,widths=c(5,1))
                    
                    #ggsave("Fig3_map.png",dpi=600)

                fmod_fft_map_low <- fmod_fft_map[which(fmod_fft_map$filters=="15 years"),]
                missingcountries15 <- fmod_fft_map_low$countrycode[is.na(fmod_fft_map_low$Estimate)]
                fmod_fft_map_10 <- fmod_fft_map[which(fmod_fft_map$filters=="10 years"),]
                
                fmod_fft_map_low$Estimate[which(fmod_fft_map_low$countrycode %in% missingcountries15)] <- fmod_fft_map_10$Estimate[which(fmod_fft_map_10$countrycode %in% missingcountries15)]
                ggplot(data = fmod_fft_map_low) +
                                geom_sf(data=fmod_fft_map_na,fill=NA)+theme_minimal()+
                                geom_sf(data=fmod_fft_map_low,aes(fill = Estimate*100))+
                                scale_fill_gradient2(
                                    name = "Estimated impact \n (% per Degree)",
                                    low = "red",
                                    mid = "white",
                                    high = "#00BFC4",
                                    midpoint = 0,
                                    space = "Lab",
                                    na.value = "gray",
                                    guide = "colourbar",
                                    aesthetics = "fill")+
                                    ggtitle("Growth effects")+
                                theme(legend.position="bottom")
                    #ggsave("Map_correcteddata.png",dpi=600)
                           
        # Categorizing statistically different estimates (end)
    ## 3.3. categorizing - Plotting Figure 3 (end) 

    ## 3.4. FELM abs estimate by filter - Table 1, columns 1 and 2 (start)
        fmod_fft$absestimate <- abs(fmod_fft$Estimate)
        fmod_fft$continent <- countrycode(sourcevar = fmod_fft$countrycode,
                                origin = "iso3c",
                                destination = "continent")
        fmod_fft$invse <- 1/fmod_fft$StandardError
        fmod_fft <- fmod_fft[!is.na(fmod_fft$invse),]
        felm_1 <- felm(absestimate ~ factor(frequencies)|0|0|continent, data = fmod_fft,weights = fmod_fft$invse)
        felm_2 <- felm(absestimate ~ factor(frequencies)|0|0|continent, data = fmod_fft)
        
        #stargazer(felm_2,felm_1, type="html", out="felm_1_2.html")
        stargazer(felm_2,felm_1,type="text")

        #felm_1_within <- felm(absestimate ~ frequencies|countrycode|0|continent, data = fmod_fft,weights = fmod_fft$invse)
        #felm_2_within <- felm(absestimate ~ frequencies|countrycode|0|continent, data = fmod_fft)
        #summary(felm_1_within)
        #stargazer(felm_2_within,felm_1_within,type="text")
        
        #stargazer(felm_2,felm_1,felm_2_within,felm_1_within, type="html", out="felm_1_2_within.html")
        
    ## 3.4. FELM abs estimate by filter - Table 1, columns 1 and 2 (end)
    
    ## 3.5. Mean effect accross filters and socioeconomic variables; Supp Fig 1 and 2 (start)
        # Adding socioeconomic variables (start)
            wb2 <- wb[!is.na(wb$gdppc),]
            glimpse(wb2)
            minyear <- aggregate(wb2$year, by = list(wb2$countrycode), FUN = min, na.rm = TRUE )
            names(minyear) <- c("countrycode","year")
            minyear$extra <- paste(minyear$countrycode,minyear$year,sep="")
            wb2$extra <- paste(wb2$countrycode,wb2$year,sep="")
            wb2 <- wb2[which(wb2$extra %in% minyear$extra),]
            wb2$loggdppc <- log(wb2$gdppc)
            wb2 <- wb2[,which(names(wb2) %in% c("countrycode","loggdppc"))]
            fmod_fft2 <- merge(fmod_fft, wb2,by = "countrycode")
            meangdppc <- aggregate(wb$gdppc, by = list(wb$countrycode), FUN = mean, na.rm = TRUE )
            names(meangdppc ) <-  c("countrycode", "meangdppc ")
            fmod_fft2 <- merge(fmod_fft2, meangdppc , by = "countrycode")
            pop_data <- wb_data("SP.POP.TOTL", start_date = 2019, end_date = 2019)
            names(pop_data)[2] <- "countrycode"
            fmod_fft2 <- merge(fmod_fft2, pop_data,by = "countrycode") 
            a <- wb_search("GDP.*PPP")
            a <- a[6,1]
            gdp_data <- wb_data(a, start_date = 2019, end_date = 2019)
            names(gdp_data)[2] <- "countrycode"
            fmod_fft2 <- merge(fmod_fft2, gdp_data,by = "countrycode")         
            a <- wb_search("GDP per capita")
            a <- a[10,1]
            gdp_data <- wb_data(a, start_date = 2019, end_date = 2019)
            names(gdp_data)[2] <- "countrycode"
            fmod_fft2 <- merge(fmod_fft2, gdp_data,by = "countrycode")
            meanT <- aggregate(wb$UDel_pop_temp, by = list(wb$countrycode), FUN = mean, na.rm = TRUE )
            names(meanT ) <-  c("countrycode", "meanT")
            fmod_fft2 <- merge(fmod_fft2, meanT,by = "countrycode")
            fmod_fft2 <- fmod_fft2[, !duplicated(colnames(fmod_fft))]
            #glimpse(fmod_fft2)
            #ggplot(data=fmod_fft2, aes(x=meanT,y= log(NY.GDP.PCAP.PP.KD), color=category))+geom_point()+theme_bw()
        # Adding socioeconomic variables (end)

        # weighted average (start)

            meanestimate <- data.frame(mean = rep(0,15),filter = rep(0,15), weight = rep(0,15), var = rep(0,15))
            frequency <- c("Unfiltered","3 years","5 years","10 years","15 years")
            #fmod_fft15 <- fmod_fft2[(fmod_fft2$filters=="Unfiltered"),]
            #fmod_fft15 <- fmod_fft15[!is.na(fmod_fft15$Estimate),]
            c15 <- levels(factor(fmod_fft2$countrycode))
            fmod_fft15 <- fmod_fft2[(fmod_fft2$countrycode %in% c15),]
            weighted.var <- function(x, w, na.rm = FALSE) {
                if (na.rm) {
                    w <- w[i <- !is.na(x)]
                    x <- x[i]
                }
                sum.w <- sum(w)
                sum.w2 <- sum(w^2)
                mean.w <- sum(x * w) / sum(w)
                (sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2, na.rm =
                na.rm)
                }
            for (i in 1:5){
                meanestimate$filter[(1 + ((i-1)*3)):(3 + ((i-1)*3))] <- frequency[i]
                meanestimate$mean[1 + ((i-1)*3)]=mean(fmod_fft15$Estimate[fmod_fft15$filters==frequency[i]],na.rm=TRUE)
                meanestimate$weight[1 + ((i-1)*3)]="Unweighted"
                meanestimate$var[1 + ((i-1)*3)] <- var(fmod_fft15$Estimate[fmod_fft15$filters==frequency[i]],na.rm=TRUE)


                x1 <- fmod_fft15[!is.na(fmod_fft15$SP.POP.TOTL),]
                x1 <- x1[!is.na(x1$Estimate),]
                est <- x1$Estimate[x1$filters==frequency[i]]
                popw <- x1$SP.POP.TOTL[x1$filters==frequency[i]]
                
                meanestimate$mean[2 + ((i-1)*3)] <-sum((est*popw)/sum(popw),na.rm=TRUE)
                meanestimate$weight[2 + ((i-1)*3)] <- "Population in 2019"
                meanestimate$var[2 + ((i-1)*3)] <- weighted.var(x=est, w = popw,na.rm=TRUE)

                #x1 <- fmod_fft15[!is.na(fmod_fft15$Estimate) & !is.na(fmod_fft15$NY.GDP.PCAP.PP.KD),]
                x1 <- fmod_fft15[!is.na(fmod_fft15$Estimate) & !is.na(fmod_fft15$NY.GDP.MKTP.PP.KD),]
                gdpw <- x1$NY.GDP.MKTP.PP.KD[x1$filters==frequency[i]]
                #gdpw <- x1$NY.GDP.PCAP.PP.KD[x1$filters==frequency[i]]
                est <- x1$Estimate[x1$filters==frequency[i]]
                meanestimate$mean[3 + ((i-1)*3)] <- sum((est*gdpw)/sum(gdpw),na.rm=TRUE)
                meanestimate$weight[3 + ((i-1)*3)] <- "GDP, PPP constant 2017 USD"
                meanestimate$var[3 + ((i-1)*3)]  <- weighted.var(x=x1$Estimate[x1$filters==frequency[i]], w = x1$NY.GDP.MKTP.PP.KD[x1$filters==frequency[i]],na.rm=TRUE)

            }
            #See weighted means
            meanestimate
        # weighted average (end)

        # Socioeconomic characteristics across categories (start)
            
            fmod_fft2$catsign <- paste(fmod_fft2$category,fmod_fft2$signlofreq,sep=".")
            fmod_fft2 <- fmod_fft2[fmod_fft2$category!="other",]
            fmod_fft2$logGDPpc <- log(fmod_fft2$NY.GDP.PCAP.PP.KD)
            densityGDP6 <-  ggstatsplot::ggbetweenstats(
                data = fmod_fft2,
                x = catsign,
                y = logGDPpc
                )

            densitymeanT6 <-  ggstatsplot::ggbetweenstats(
            data = fmod_fft2,
            x = catsign,
            y = meanT
            )

            densityGDP3 <-  ggstatsplot::ggbetweenstats(
            data = fmod_fft2,
            x = category,
            y = logGDPpc
            )

            densitymeanT3 <-  ggstatsplot::ggbetweenstats(
            data = fmod_fft2,
            x = category,
            y = meanT
            )

            ggarrange(densityGDP3,densitymeanT3)
            ggsave("Category_byTandGDP.png",dpi=600)
        
            cat <- data.frame(T = rep(0,6),minT=rep(0,6),maxT=rep(0,6),G = rep(0,6),minG=rep(0,6),maxG=rep(0,6),categories=c("Constant.negative","Constant.positive","Converging.negative","Converging.positive","Intensifying.negative","Intensifying.positive"))
            cat$T <- aggregate(fmod_fft2$meanT, list(fmod_fft2$catsign), mean, na.rm=TRUE)[1:6,2]
            cat$G <- aggregate(fmod_fft2$logGDPpc, list(fmod_fft2$catsign), mean, na.rm=TRUE)[1:6,2]
            cat$minT <- aggregate(fmod_fft2$meanT, list(fmod_fft2$catsign), min, na.rm=TRUE)[1:6,2]
            cat$maxT <- aggregate(fmod_fft2$meanT, list(fmod_fft2$catsign),max, na.rm=TRUE)[1:6,2]
            cat$minG <- aggregate(fmod_fft2$logGDPpc, list(fmod_fft2$catsign), min, na.rm=TRUE)[1:6,2]
            cat$maxG <- aggregate(fmod_fft2$logGDPpc, list(fmod_fft2$catsign),max, na.rm=TRUE)[1:6,2]
            cat$category <- c("Constant","Constant","Converging","Converging","Intensifying","Intensifying")
            cat$signlofreq <- c("Negative","Positive","Negative","Positive","Negative","Positive")
            
            cat6spread <- ggplot(data=cat, aes(x=T,y=G,xmin=minT,xmax=maxT,ymin=minG,ymax=maxG, shape=category, color = signlofreq,group=categories))+
            geom_point()+geom_errorbar()+geom_errorbarh()+
            scale_color_discrete(name="Direction of the effect")+
            scale_color_brewer(name="Direction of the effect",palette="Dark2")+
            scale_shape_discrete(name="Category")+
            labs(x="Country's population-weighted temperature", y="log(GDP per capita)")+theme_bw()

            cat6datapoints <- ggplot(data=fmod_fft2, aes(x=meanT,y=logGDPpc, shape=category, color = factor(signlofreq)))+
            geom_point()+
            scale_color_discrete(name="Direction of the effect")+
            scale_color_brewer(name="Direction of the effect",palette="Dark2")+
            scale_shape_discrete(name="Category")+
            labs(x="Country's population-weighted temperature", y="log(GDP per capita)")+theme_bw()

            ggarrange(cat6spread,cat6datapoints,ncol=1,common.legend=TRUE)
            ggsave("SupFig3_cat6.png",dpi=600)
            cat <- data.frame(T = rep(0,3),minT=rep(0,3),maxT=rep(0,3),G = rep(0,3),minG=rep(0,3),maxG=rep(0,3),category=c("Constant","Converging","Intensifying"))
            cat$T <- aggregate(fmod_fft2$meanT, list(fmod_fft2$category), mean, na.rm=TRUE)[,2]
            cat$G <- aggregate(fmod_fft2$logGDPpc, list(fmod_fft2$category), mean, na.rm=TRUE)[,2]
            cat$minT <- aggregate(fmod_fft2$meanT, list(fmod_fft2$category), min, na.rm=TRUE)[,2]
            cat$maxT <- aggregate(fmod_fft2$meanT, list(fmod_fft2$category),max, na.rm=TRUE)[,2]
            cat$minG <- aggregate(fmod_fft2$logGDPpc, list(fmod_fft2$category), min, na.rm=TRUE)[,2]
            cat$maxG <- aggregate(fmod_fft2$logGDPpc, list(fmod_fft2$category),max, na.rm=TRUE)[,2]
            cat3spread <- ggplot(data=cat, aes(x=T,y=G,xmin=minT,xmax=maxT,ymin=minG,ymax=maxG,color=category))+
            geom_point()+geom_errorbar()+geom_errorbarh()+
            labs(x="Mean Temperature", y="logGDPpc")+theme_bw()

            cat3datapoints <- ggplot(data=fmod_fft2, aes(x=meanT,y=logGDPpc, color=category))+
            geom_point()+
            labs(x="Mean Temperature", y="logGDPpc")+theme_bw()

            ggarrange(densityGDP6,densitymeanT6,cat6spread,cat6datapoints,ncol=1,nrow=4)
            #ggsave("soecioeconomic_6cat.png",dpi=600)
            ggarrange(densityGDP3,densitymeanT3,cat3spread,cat3datapoints,ncol=1,nrow=4)
            #ggsave("soecioeconomic_3cat.png",dpi=600)
            head(fmod_fft2[fmod_fft2$filters=="Unfiltered",])
            a1=ggplot(fmod_fft2[fmod_fft2$filters=="Unfiltered",], 
            aes(x=meanT,y=Estimate))+geom_point()+
            geom_smooth(method = "lm",
                    color = "#00BFC4", show.legend = FALSE) +  theme_bw()+
                    geom_hline(yintercept=0,lty=2) +
                    labs(title="Temperature effect Filter = Unfiltered")
            

            b1=ggplot(fmod_fft2[fmod_fft2$filters=="Unfiltered",],aes(x=log(NY.GDP.PCAP.PP.KD),
            y=Estimate))+geom_point()+
            #geom_errorbar()+
            geom_smooth(method = "lm", 
                    color = "#00BFC4", show.legend = FALSE) +  theme_bw()+
                    geom_hline(yintercept=0,lty=2) +
                    xlab("log GDP per capita")+
                    labs(title="Temperature effect Filter = Unfiltered")
            a=ggplot(fmod_fft2[fmod_fft2$filters=="15 years",], 
            aes(x=meanT,y=Estimate))+geom_point()+
            geom_smooth(method = "lm",
                    color = "#00BFC4", show.legend = FALSE) +  theme_bw()+
                    geom_hline(yintercept=0,lty=2) +
                    labs(title="Temperature effect Filter = 15 year")
            

            b=ggplot(fmod_fft2[fmod_fft2$filters=="15 years",],aes(x=log(NY.GDP.PCAP.PP.KD),
            y=Estimate))+geom_point()+
            #geom_errorbar()+
            geom_smooth(method = "lm", 
                    color = "#00BFC4", show.legend = FALSE) +  theme_bw()+
                    geom_hline(yintercept=0,lty=2) +
                    xlab("log GDP per capita")+
                    labs(title="Temperature effect Filter = 15 year")
            ggarrange(a1,a, b1, b)
            
             #ggsave("lowfreq_meanTlogGDP.png",dpi=600)
        # Socioeconomic characteristics across categories (end)
    ## 3.5. Mean effect accross filters and socioeconomic variables; Supp Fig 1 and 2 (end)

    ## 3.6. Other datasets - Table 1, columns 3 to 6; Supp Fig 3 (start)
        # Barro
            ## 3.2. categorizing - Plotting Figure 3  (start)

                    # Categorizing statistically different estimates (start)
                        fmod_fft <- new_fft
                        fmod_fft$Variance <- fullmods_filter_var$Variance
                        
                        #fmod_fft <- fullmods_filter
                        fmod_fft <- fmod_fft[fmod_fft$econdata=="barro",]
                        fmod_fft <- fmod_fft[fmod_fft$climdata=="UDel",]
                        fmod_fft$high95 <- fmod_fft$Estimate + fmod_fft$StandardError*1.64
                        fmod_fft$low95 <- fmod_fft$Estimate - fmod_fft$StandardError*1.64
                        numcountries <- length(unique(fmod_fft$countrycode))
                        #fmod_fft$absestimate <- abs(fmod_fft$Estimate)
                        fmod_fft$lowsignificant <- 0
                        fmod_fft$unfilteredsignificant <- 0
                                        
                        uncertain <- c(which(fmod_fft$high95>0 & fmod_fft$low95 <0))
                        fmod_fft$sign <- rep(0,dim(fmod_fft)[1])
                        
                        numcountries <- length(unique(fmod_fft$countrycode))
                        for (i in 1:numcountries){
                            if(is.null(fmod_fft$Estimate[1+5*(i-1)] )){next}
                            if(!is.na(fmod_fft$Estimate[5+5*(i-1)])){
                                lastfreq <- 5
                                }else if(!is.na(fmod_fft$Estimate[4+5*(i-1)])){
                                    lastfreq <- 4
                                }else if(!is.na(fmod_fft$Estimate[3+5*(i-1)])){
                                    lastfreq <- 3
                                }else if(!is.na(fmod_fft$Estimate[2+5*(i-1)])){
                                    lastfreq <- 2
                                }else{next}
                            m <- fmod_fft$Estimate[1+5*(i-1)] * fmod_fft$Estimate[lastfreq+5*(i-1)]
                            if (is.na(m)){next}
                            if(m>0 ){
                                fmod_fft$sign[(1+5*(i-1)):(5+5*(i-1))]=1
                                }
                        }

                            converging_stat <- 0
                            converging <- 0
                            intensifying <- 0
                            intensifying_stat <- 0
                            for (i in 1:numcountries){
                                if(!is.na(fmod_fft$Estimate[5+5*(i-1)])){
                                        lastfreq <- 5
                                        }else if(!is.na(fmod_fft$Estimate[4+5*(i-1)])){
                                            lastfreq <- 4
                                        }else if(!is.na(fmod_fft$Estimate[3+5*(i-1)])){
                                            lastfreq <- 3
                                        }else if(!is.na(fmod_fft$Estimate[2+5*(i-1)])){
                                            lastfreq <- 2
                                        }else{next}
                                        
                                        if(is.na(fmod_fft$Estimate[1+5*(i-1)])){next}
                                theta <- abs(fmod_fft$Estimate[1+5*(i-1)]) - abs(fmod_fft$Estimate[lastfreq+5*(i-1)])
                                var <- fmod_fft$Variance[lastfreq+5*(i-1)] + fmod_fft$Variance[1+5*(i-1)]
                                conf95 <-  (var^0.5)*1.65 #one-tail  95%
                                
                                if(abs(fmod_fft$Estimate[lastfreq+5*(i-1)]) > abs(fmod_fft$Estimate[1+5*(i-1)])){
                                    intensifying <- c(intensifying,(1+5*(i-1)):(5+5*(i-1))) 
                                    if(theta+conf95 < 0 ){
                                        intensifying_stat <- c(intensifying_stat,(1+5*(i-1)):(5+5*(i-1))) 
                                    }
                                }
                                
                                if(abs(fmod_fft$Estimate[lastfreq+5*(i-1)]) <= abs(fmod_fft$Estimate[1+5*(i-1)])){
                                        converging <- c(converging,(1+5*(i-1)):(5+5*(i-1))) 
                                        if(theta-conf95 >= 0 ){
                                        converging_stat <- c(converging_stat,(1+5*(i-1)):(5+5*(i-1))) 
                                    }
                                }
                                
                                
                                    
                                    if((lastfreq+5*(i-1)) %notin% uncertain){
                                    fmod_fft$lowsignificant[(1+5*(i-1)):(5+5*(i-1))] <- 1     
                                                }
                                    if((1+5*(i-1)) %notin% uncertain){
                                        fmod_fft$unfilteredsignificant[(1+5*(i-1)):(5+5*(i-1))] <- 1     
                                                }
                                }
                                
                                
                            filt_names <- c("Unfiltered","3 years","5 years", "10 years", "15 years")
                            
                            fmod_fft$filters <- rep(filt_names,length(unique(wb$countrycode)))
                            fmod_fft$filters <- factor(fmod_fft$filters, levels = filt_names)            
                            fmod_fft$category <- "Undefined"
                            fmod_fft$category[converging] <- "converging"
                            fmod_fft$category[intensifying] <- "intensifying"
                            fmod_fft$category[converging_stat] <- "converging_stat"
                            fmod_fft$category[intensifying_stat] <- "intensifying_stat"
                            table(fmod_fft$category)/5
                            
                            fmod_fft$Category <- "Undefined"
                            fmod_fft$Category[converging] <- "converging"
                            fmod_fft$Category[intensifying] <- "intensifying"
                            table(fmod_fft$Category)/5
                            
                            fmod_fft$CoeffDifferents <- 0
                            fmod_fft$CoeffDifferents[intensifying_stat] <- 1
                            fmod_fft$CoeffDifferents[converging_stat] <- 1
                            table(fmod_fft$CoeffDifferents)/5
                            
                            table(fmod_fft$lowsignificant)/5

                            fmod_fft_na<-fmod_fft[fmod_fft$category=="Undefined" | 
                                fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate<quantile(fmod_fft$Estimate, na.rm=TRUE,0.01))] |
                                fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate>quantile(fmod_fft$Estimate, na.rm=TRUE,0.99))],]
                            fmod_fft<-fmod_fft[fmod_fft$category!="Undefined",]
                            world <- ne_countries(scale = "medium", returnclass = "sf")

                            
                            fmod_fft<- fmod_fft[fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate>quantile(fmod_fft$Estimate, na.rm=TRUE,0.99))],]
                            fmod_fft<- fmod_fft[fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate<quantile(fmod_fft$Estimate, na.rm=TRUE,0.01))],]
                        
                            fmod_fft$iso_a3 <- fmod_fft$countrycode
                            fmod_fft_map <- merge(world,fmod_fft,by="iso_a3")
                            fmod_fft_na$iso_a3 <- fmod_fft_na$countrycode
                            fmod_fft_map_na <- merge(world,fmod_fft_na,by="iso_a3")
                            

                            #factor(fmod_fft$category)
                            fmod_fft$category <- factor(fmod_fft$category, levels = c("intensifying_stat", "intensifying", "converging","converging_stat"))
                            df1 <- fmod_fft[which(fmod_fft$econdata=="barro" &fmod_fft$climdata=="UDel" & fmod_fft$filter=="Unfiltered"),]
                            df1$unfilteredsignificant <- factor(df1$unfilteredsignificant)
                            df1$lowsignificant <- factor(df1$lowsignificant)
                                
                            hist_bf <- ggplot() + theme_bw() + 
                                geom_histogram(data=df1, na.rm= TRUE, mapping = aes(x=lowsignificant, fill=factor(category)), 
                                #stat='count')+ scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],pal2[1]), 
                                #name="",labels=c("*Intensifying","Intensifying","Converging","*Converging"))+
                                stat='count')+ scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3]), 
                                name="",labels=c("*Intensifying","Intensifying","Converging"))+ #No converging statistically significant found
                                scale_x_discrete(name="",
                                    labels=c(expression(paste("|",theta[f],"| =0",sep="")),expression(paste("|",theta[f],"| >0",sep=""))))+
                                scale_y_continuous(name="Number of countries")+theme(legend.position="none")

                            df2 <- fmod_fft[which(fmod_fft$econdata=="barro" &fmod_fft$climdata=="UDel" & fmod_fft$filter=="Unfiltered" & fmod_fft$unfilteredsignificant==1),]
                            df2$unfilteredsignificant <- factor(df2$unfilteredsignificant)
                            glimpse(df2)     
                            hist_bu <- ggplot() + theme_bw()+
                                geom_histogram(data=df2, na.rm= TRUE, mapping = aes(x=unfilteredsignificant, fill=factor(category)),stat='count')+ 
                                #scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],pal2[1]), 
                                #name="",labels=c("*Intensifying","Intensifying","Converging","*Converging"))+
                                scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3]), 
                                name="",labels=c("*Intensifying","Intensifying","Converging"))+ #No converging statistically significant found
                                scale_x_discrete(name="",labels=c(expression(paste("|",beta[U],"| >0",sep=""))))+
                                scale_y_continuous(name="Number of countries") + theme(legend.position="none")
                            
                            plot_fg <- ggplot(data=fmod_fft,aes(x=filters,y=Estimate*100, group = countrycode,color=factor(category)))+
                                geom_line()+
                                scale_colour_manual(name="Categories", values=c(pal1[1],pal1[3],pal2[3],pal2[1]),labels=c("Converging","Intensifying","Statistically Intensifying","Statistically Converging")) +
                                theme_bw() + xlab("Minimum Periodicity after Filtering")+
                                geom_hline(yintercept=0,lty=2)+
                                #geom_dl(data=fg,aes(label = countrycode), method = list(dl.combine("last.points")), cex = 0.9)+
                                ylab("Estimated Effect of \n 1 Degree Warming on Growth (pp)") +
                                theme(axis.line = element_line(colour = "black"),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                panel.border = element_blank(),
                                panel.background = element_blank(),legend.position="none") +
                                #panel.background = element_blank()) +
                                ggtitle("")
                                plot_fg      
                            
                            
                            plot_categoriescolors <- ggplot() + geom_histogram(data=data.frame(x=c("*Intensifying","Intensifying","Converging","*Converging")), na.rm= TRUE, mapping = aes(x=x, fill=factor(x)),stat='count')+ 
                            scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],pal2[1]),name="",labels=c("*Intensifying","Intensifying","Converging","*Converging"))
                            leg_cat <- get_legend(plot_categoriescolors)                    
                            
                            fig_3 <- ggarrange(plot_fg,ggarrange(hist_bu,leg_cat,hist_bf,ncol=3,nrow=1),ncol=1,nrow=2)    
                            fig_3
                            #ggsave("Fig3_barro.png",dpi=600)
                                    
                    # Categorizing statistically different estimates (end)
            ## 3.2. categorizing - Plotting Figure 3 (end) 

            ## 3.3. FELM abs estimate by filter - Table 1, columns 1 and 2 (start)
                fmod_fft$absestimate <- abs(fmod_fft$Estimate)
                fmod_fft$continent <- countrycode(sourcevar = fmod_fft$countrycode,
                                        origin = "iso3c",
                                        destination = "continent")
                fmod_fft$invse <- 1/fmod_fft$StandardError
                fmod_fft <- fmod_fft[!is.na(fmod_fft$invse),]
                felm_3 <- felm(absestimate ~ factor(frequencies)|0|0|continent, data = fmod_fft,weights = fmod_fft$invse)
                felm_4 <- felm(absestimate ~ factor(frequencies)|0|0|continent, data = fmod_fft)
            ## 3.3. FELM abs estimate by filter - Table 1, columns 1 and 2 (end)
        # Barro  

        # Maddison
            ## 3.2. categorizing - Plotting Figure 3  (start)

                    # Categorizing statistically different estimates (start)
                        fmod_fft <- new_fft
                        fmod_fft$Variance <- fullmods_filter_var$Variance
                        
                        #fmod_fft <- fullmods_filter
                        fmod_fft <- fmod_fft[fmod_fft$econdata=="mad",]
                        fmod_fft <- fmod_fft[fmod_fft$climdata=="UDel",]
                        fmod_fft$high95 <- fmod_fft$Estimate + fmod_fft$StandardError*1.64
                        fmod_fft$low95 <- fmod_fft$Estimate - fmod_fft$StandardError*1.64
                        numcountries <- length(unique(fmod_fft$countrycode))
                        #fmod_fft$absestimate <- abs(fmod_fft$Estimate)
                        fmod_fft$lowsignificant <- 0
                        fmod_fft$unfilteredsignificant <- 0
                                        
                        uncertain <- c(which(fmod_fft$high95>0 & fmod_fft$low95 <0))
                        fmod_fft$sign <- rep(0,dim(fmod_fft)[1])
                        
                        numcountries <- length(unique(fmod_fft$countrycode))
                        for (i in 1:numcountries){
                            if(is.null(fmod_fft$Estimate[1+5*(i-1)] )){next}
                            if(!is.na(fmod_fft$Estimate[5+5*(i-1)])){
                                lastfreq <- 5
                                }else if(!is.na(fmod_fft$Estimate[4+5*(i-1)])){
                                    lastfreq <- 4
                                }else if(!is.na(fmod_fft$Estimate[3+5*(i-1)])){
                                    lastfreq <- 3
                                }else if(!is.na(fmod_fft$Estimate[2+5*(i-1)])){
                                    lastfreq <- 2
                                }else{next}
                            m <- fmod_fft$Estimate[1+5*(i-1)] * fmod_fft$Estimate[lastfreq+5*(i-1)]
                            if (is.na(m)){next}
                            if(m>0 ){
                                fmod_fft$sign[(1+5*(i-1)):(5+5*(i-1))]=1
                                }
                        }

                            converging_stat <- 0
                            converging <- 0
                            intensifying <- 0
                            intensifying_stat <- 0
                            for (i in 1:numcountries){
                                if(!is.na(fmod_fft$Estimate[5+5*(i-1)])){
                                        lastfreq <- 5
                                        }else if(!is.na(fmod_fft$Estimate[4+5*(i-1)])){
                                            lastfreq <- 4
                                        }else if(!is.na(fmod_fft$Estimate[3+5*(i-1)])){
                                            lastfreq <- 3
                                        }else if(!is.na(fmod_fft$Estimate[2+5*(i-1)])){
                                            lastfreq <- 2
                                        }else{next}
                                        
                                        if(is.na(fmod_fft$Estimate[1+5*(i-1)])){next}
                                theta <- abs(fmod_fft$Estimate[1+5*(i-1)]) - abs(fmod_fft$Estimate[lastfreq+5*(i-1)])
                                var <- fmod_fft$Variance[lastfreq+5*(i-1)] + fmod_fft$Variance[1+5*(i-1)]
                                conf95 <-  (var^0.5)*1.65 #one-tail  95%
                                
                                if(abs(fmod_fft$Estimate[lastfreq+5*(i-1)]) > abs(fmod_fft$Estimate[1+5*(i-1)])){
                                    intensifying <- c(intensifying,(1+5*(i-1)):(5+5*(i-1))) 
                                    if(theta+conf95 < 0 ){
                                        intensifying_stat <- c(intensifying_stat,(1+5*(i-1)):(5+5*(i-1))) 
                                    }
                                }
                                
                                if(abs(fmod_fft$Estimate[lastfreq+5*(i-1)]) <= abs(fmod_fft$Estimate[1+5*(i-1)])){
                                        converging <- c(converging,(1+5*(i-1)):(5+5*(i-1))) 
                                        if(theta-conf95 >= 0 ){
                                        converging_stat <- c(converging_stat,(1+5*(i-1)):(5+5*(i-1))) 
                                    }
                                }
                                
                                
                                    
                                    if((lastfreq+5*(i-1)) %notin% uncertain){
                                    fmod_fft$lowsignificant[(1+5*(i-1)):(5+5*(i-1))] <- 1     
                                                }
                                    if((1+5*(i-1)) %notin% uncertain){
                                        fmod_fft$unfilteredsignificant[(1+5*(i-1)):(5+5*(i-1))] <- 1     
                                                }
                                }
                                
                                
                            filt_names <- c("Unfiltered","3 years","5 years", "10 years", "15 years")
                            
                            fmod_fft$filters <- rep(filt_names,length(unique(wb$countrycode)))
                            fmod_fft$filters <- factor(fmod_fft$filters, levels = filt_names)            
                            fmod_fft$category <- "Undefined"
                            fmod_fft$category[converging] <- "converging"
                            fmod_fft$category[intensifying] <- "intensifying"
                            fmod_fft$category[converging_stat] <- "converging_stat"
                            fmod_fft$category[intensifying_stat] <- "intensifying_stat"
                            table(fmod_fft$category)/5
                            
                            fmod_fft$Category <- "Undefined"
                            fmod_fft$Category[converging] <- "converging"
                            fmod_fft$Category[intensifying] <- "intensifying"
                            table(fmod_fft$Category)/5
                            
                            fmod_fft$CoeffDifferents <- 0
                            fmod_fft$CoeffDifferents[intensifying_stat] <- 1
                            fmod_fft$CoeffDifferents[converging_stat] <- 1
                            table(fmod_fft$CoeffDifferents)/5
                            
                            table(fmod_fft$lowsignificant)/5

                            fmod_fft_na<-fmod_fft[fmod_fft$category=="Undefined" | 
                                fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate<quantile(fmod_fft$Estimate, na.rm=TRUE,0.01))] |
                                fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate>quantile(fmod_fft$Estimate, na.rm=TRUE,0.99))],]
                            fmod_fft<-fmod_fft[fmod_fft$category!="Undefined",]
                            world <- ne_countries(scale = "medium", returnclass = "sf")

                            
                            fmod_fft<- fmod_fft[fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate>quantile(fmod_fft$Estimate, na.rm=TRUE,0.99))],]
                            fmod_fft<- fmod_fft[fmod_fft$countrycode %notin% fmod_fft$countrycode[which(fmod_fft$Estimate<quantile(fmod_fft$Estimate, na.rm=TRUE,0.01))],]
                        
                            fmod_fft$iso_a3 <- fmod_fft$countrycode
                            fmod_fft_map <- merge(world,fmod_fft,by="iso_a3")
                            fmod_fft_na$iso_a3 <- fmod_fft_na$countrycode
                            fmod_fft_map_na <- merge(world,fmod_fft_na,by="iso_a3")
                            

                            #factor(fmod_fft$category)
                            fmod_fft$category <- factor(fmod_fft$category, levels = c("intensifying_stat", "intensifying", "converging","converging_stat"))
                            df1 <- fmod_fft[which(fmod_fft$econdata=="mad" &fmod_fft$climdata=="UDel" & fmod_fft$filter=="Unfiltered"),]
                            df1$unfilteredsignificant <- factor(df1$unfilteredsignificant)
                            df1$lowsignificant <- factor(df1$lowsignificant)
                                
                            hist_bf <- ggplot() + theme_bw() + 
                                geom_histogram(data=df1, na.rm= TRUE, mapping = aes(x=lowsignificant, fill=factor(category)), 
                                #stat='count')+ scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],pal2[1]), 
                                #name="",labels=c("*Intensifying","Intensifying","Converging","*Converging"))+
                                stat='count')+ scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3]), 
                                name="",labels=c("*Intensifying","Intensifying","Converging"))+ #No converging statistically significant found
                                scale_x_discrete(name="",
                                    labels=c(expression(paste("|",theta[f],"| =0",sep="")),expression(paste("|",theta[f],"| >0",sep=""))))+
                                scale_y_continuous(name="Number of countries")+theme(legend.position="none")

                            df2 <- fmod_fft[which(fmod_fft$econdata=="mad" &fmod_fft$climdata=="UDel" & fmod_fft$filter=="Unfiltered" & fmod_fft$unfilteredsignificant==1),]
                            df2$unfilteredsignificant <- factor(df2$unfilteredsignificant)
                            glimpse(df2)     
                            hist_bu <- ggplot() + theme_bw()+
                                geom_histogram(data=df2, na.rm= TRUE, mapping = aes(x=unfilteredsignificant, fill=factor(category)),stat='count')+ 
                                #scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],pal2[1]), 
                                #name="",labels=c("*Intensifying","Intensifying","Converging","*Converging"))+
                                scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3]), 
                                name="",labels=c("*Intensifying","Intensifying","Converging"))+ #No converging statistically significant found
                                scale_x_discrete(name="",labels=c(expression(paste("|",beta[U],"| >0",sep=""))))+
                                scale_y_continuous(name="Number of countries") + theme(legend.position="none")
                            
                            plot_fg <- ggplot(data=fmod_fft,aes(x=filters,y=Estimate*100, group = countrycode,color=factor(category)))+
                                geom_line()+
                                scale_colour_manual(name="Categories", values=c(pal1[1],pal1[3],pal2[3],pal2[1]),labels=c("Converging","Intensifying","Statistically Intensifying","Statistically Converging")) +
                                theme_bw() + xlab("Minimum Periodicity after Filtering")+
                                geom_hline(yintercept=0,lty=2)+
                                #geom_dl(data=fg,aes(label = countrycode), method = list(dl.combine("last.points")), cex = 0.9)+
                                ylab("Estimated Effect of \n 1 Degree Warming on Growth (pp)") +
                                theme(axis.line = element_line(colour = "black"),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                panel.border = element_blank(),
                                panel.background = element_blank(),legend.position="none") +
                                #panel.background = element_blank()) +
                                ggtitle("")
                                plot_fg      
                            
                            
                            plot_categoriescolors <- ggplot() + geom_histogram(data=data.frame(x=c("*Intensifying","Intensifying","Converging","*Converging")), na.rm= TRUE, mapping = aes(x=x, fill=factor(x)),stat='count')+ 
                            scale_fill_manual(values = c(pal1[1],pal1[3],pal2[3],pal2[1]),name="",labels=c("*Intensifying","Intensifying","Converging","*Converging"))
                            leg_cat <- get_legend(plot_categoriescolors)                    
                            
                            fig_3 <- ggarrange(plot_fg,ggarrange(hist_bu,leg_cat,hist_bf,ncol=3,nrow=1),ncol=1,nrow=2)    
                            fig_3
                            #ggsave("Fig3_mad.png",dpi=600)
                                    
                    # Categorizing statistically different estimates (end)
            ## 3.2. categorizing - Plotting Figure 3 (end) 

            ## 3.3. FELM abs estimate by filter - Table 1, columns 1 and 2 (start)
                fmod_fft$absestimate <- abs(fmod_fft$Estimate)
                fmod_fft$continent <- countrycode(sourcevar = fmod_fft$countrycode,
                                        origin = "iso3c",
                                        destination = "continent")
                fmod_fft$invse <- 1/fmod_fft$StandardError
                fmod_fft <- fmod_fft[!is.na(fmod_fft$invse),]
                felm_5 <- felm(absestimate ~ factor(frequencies)|0|0|continent, data = fmod_fft,weights = fmod_fft$invse)
                felm_6 <- felm(absestimate ~ factor(frequencies)|0|0|continent, data = fmod_fft)
                #stargazer(felm_2,felm_1,felm_4,felm_3,felm_6,felm_5,type="html",out="Regression_on_estimates.html")
            ## 3.3. FELM abs estimate by filter - Table 1, columns 1 and 2 (end)
        # Maddison
    ## 3.5. Other datasets - Table 1, columns 3 to 6; Supp Fig 3 (end)

## 3. Empricial analysis - Figure 3 and Table 1  (end)

