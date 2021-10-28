# # Persistent Effect of Temperature on GDP Identified from Lower Frequency Temperature Variability
# # Bastien-Olvera, B. A. and Moore, F. C.

# ## Index
#     # 1. Setup
#     # 2. Simulation (Figure 1 and 2)
#         # 2.1 Temperature fluctuations and filters - Figure 1
#         # 2.2 Global temperature simulation
#         # 2.3. Perform simulation - Figure 2
#     # 3. Empirical analysis (Figure 3 and Table 1)
#         # 3.1 Country-level regressions
#         # 3.2. Categorizing - Figure 3
#         # 3.3. FELM abs estimate by filter - Table 1, columns 1 and 2
#         # 3.4. Mean effect accross filters and socioeconomic variables; Supp Fig 1 and 2
#         # 3.5. Other datasets - Table 1, columns 3 to 6; Supp Fig 3

## 1. Setup (start)
    #set working directory
    dir <- "C:/Users/bastien/Documents/GitHub/TempEffectGDP/"
    setwd(dir)
    x<-c("TTR", "ggpubr","mFilter", "ncdf4","TSA","tidyverse","ggplot2",
        "ggpubr","dplR","reshape2","raster","dplyr",
        "RColorBrewer","colorspace","spData","sf","countrycode", "ggstatsplot",
        "ggsignif","dlnm","lfe","directlabels","splines","timeSeries",
        "MASS","stargazer","jtools","wbstats","dplR")
    lapply(x, require, character.only = TRUE)
    # some functions
    tsyears <- function(ts) as.numeric(trunc(time(ts)))
    cols <- c("#38170B","#BF1B0B", "#FFC465", "#66ADE5", "#252A52")
    cols=rev(c("#76d3ae","#0dcfca","#055692"))

    #Load data
    barro <- read.csv(paste(dir,"Data/Barro_UDel.csv",sep=""))
    mad <- read.csv(paste(dir,"Data/Maddison_UDel.csv",sep="")) 
    wb <- read.csv(paste(dir,"Data/WB_UDel.csv",sep=""))   
    countries <- unique(factor(wb$countrycode))


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
    ## 2.1 Temperature fluctuations and filters - Figure 1 (end)

    ## 2.2. Global temperature simulation (start)
        file=nc_open(paste(dir,"Data/gmt_MCruns_ensemble_full.nc",sep=""))
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
        gtemp=ncvar_get(file,"gmt")
        gtemp=apply(gtemp,MARGIN=3,FUN=mean)
        
        #take first 1500 years, prior to anthropogenic influence
        #gtemp=gtemp[1:1500]

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
    ## 2.2. Global temperature (end)

    ## 2.3. Perform simulation - WITH QUADRATIC TERM (start)


        bhm <- function(T){
                    G<- 0.0127*T-0.0005*(T^2)
                    return(G)
                }
                T <- seq(13,16,length=100)
                response <- bhm(T)
                ustemprange <- range(wbUS$UDel_pop_temp,na.rm=TRUE)
                ust <- seq(from=min(wbUS$UDel_pop_temp,na.rm=TRUE),to=max(wbUS$UDel_pop_temp,na.rm=TRUE),length=100)         
                curveT <- data.frame(T,response)

                tan <- function(T){
                    slope <-  0.0127-2*0.0005*mean(ust)
                    intercept <- (0.0127*mean(ust)-0.0005*(mean(ust)^2) )-(slope*mean(ust))
                    Tan <- intercept + slope * (T)
                    return(Tan)
                }
                
                Tan<-tan(ust)
                plot(T,response)
                lines(ust,Tan)

                abs(ust[1]-mean(ust)) #1 degree of difference
                plot(ust,tan(ust)-bhm(ust))

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
    ## 2.3. Perform simulation - WITH QUADRATIC TERM (end)


    ## 2.4. Perform simulation - WITH DRIFT ADJUSTMENT (start)


            time=350 #or 350 years
            basegr=0.01 #2% per year baseline growth rate
            start=100
            baseline=rep(basegr,time)
            coef=-0.05 #effect size - change in growth per degree warming()
            coef2=-0.05
            growthsd=0.005 #standard deviation of growth variability unexplained by temperature
            periods <- c(0,10,20,50,100)
            nsims=500
            sims=array(dim=c(nsims,length(periods),2))
            sims_adjusted=array(dim=c(nsims,length(periods),2))
            for(i in 1:nsims){
                randomtemp=Re(randomts(gtemp))[1:time]
                #randomtemp=Re(randomts(ustemp_lmr))[1:time]
                randomgrowth_g=basegr #growth impact model
                randomgrowth_l=basegr #levels impact model
                for(j in 2:time){
                #quadratic
                #randomgrowth_g=c(randomgrowth_g, basegr+(randomtemp[j]*coef)+(coef2*randomtemp[j]^2)+rnorm(1,mean=0,sd=growthsd))
                #randomgrowth_l=c(randomgrowth_l, basegr+(randomtemp[j]-randomtemp[j-1])*coef+coef2*(randomtemp[j]^2-randomtemp[j-1]^2)+rnorm(1,mean=0,sd=growthsd))
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
                    #ratio <- median(randomtemp / tempts)
                    #ratio <- median( abs(randomtemp) / abs(tempts))
                    growth <- randomgrowth_g
                    level <- randomgrowth_l
                    filt <- data.frame(temp = unclass(tempts), growth = unclass(growth),
                        level = unclass(level))
                    filt$templag=c(NA,filt$temp[1:(dim(filt)[1]-1)])
                    names(filt) <- c("temp","growth","levels","templag")
                    mod_gfilt=lm(growth~temp,data=filt)$coef[2]
                    mod_lfilt=lm(levels~temp,data=filt)$coef[2]
                    sims[i,k,]=c(mod_gfilt,mod_lfilt)
                    sims_adjusted[i,k,]=c(mod_gfilt/ratio,mod_lfilt/ratio)
                }
                

                }
                if(i%%50==0) print(i)
            }
            # Comparing quadratic
          
                # Test 
                    simsfiltmean=apply(sims_adjusted,MARGIN=c(2,3),FUN=mean)
                    simsfiltsd=apply(sims_adjusted,MARGIN=c(2,3),FUN=sd)
                    colnames(simsfiltmean)=c("Growth","Level");rownames(simsfiltmean)=periods
                    simsfiltdat=cbind(melt(simsfiltmean),melt(simsfiltsd)[,3])
                    colnames(simsfiltdat)=c("periodsregationPeriod","ImpactType","MeanEffect","SDEffect")
                    theme_set(theme_bw(base_size = 20))
                    
                    simsfiltdat$periodsregationPeriod <- c("Unfiltered", "10 years", "20 years", "50 years", "100 years","Unfiltered", "10 years", "20 years", "50 years", "100 years")
                    simsfiltdat$periodsregationPeriod <- factor(simsfiltdat$periodsregationPeriod,
                        levels = c("Unfiltered", "10 years", "20 years", "50 years", "100 years"))
                    simsfiltdat_adjusted <- simsfiltdat

                    simsfiltmean=apply(sims,MARGIN=c(2,3),FUN=mean)
                    simsfiltsd=apply(sims,MARGIN=c(2,3),FUN=sd)
                    colnames(simsfiltmean)=c("Growth","Level");rownames(simsfiltmean)=periods
                    simsfiltdat=cbind(melt(simsfiltmean),melt(simsfiltsd)[,3])
                    colnames(simsfiltdat)=c("periodsregationPeriod","ImpactType","MeanEffect","SDEffect")
                    theme_set(theme_bw(base_size = 20))
                    
                    #simsfiltdat$periodsregationPeriod <- c("Unfiltered", "3 years", "5 years", "10 years", "15 years","Unfiltered", "3 years", "5 years", "10 years", "15 years")
                    simsfiltdat$periodsregationPeriod <- c("Unfiltered", "10 years", "20 years", "50 years", "100 years","Unfiltered", "10 years", "20 years", "50 years", "100 years")
                    simsfiltdat$periodsregationPeriod <- factor(simsfiltdat$periodsregationPeriod,
                        levels = c("Unfiltered", "10 years", "20 years", "50 years", "100 years"))
                    
                    
                    a=ggplot(simsfiltdat_adjusted,aes(x=factor(periodsregationPeriod),y=MeanEffect*100,col=ImpactType,group=ImpactType))+geom_line(lwd=1.25)
                    a=a+geom_point(size=4)+geom_errorbar(aes(ymin=(MeanEffect-1.96*SDEffect)*100,ymax=(MeanEffect+1.96*SDEffect)*100),width=0.1,lwd=1.25)
                    a=a+geom_hline(yintercept=0,lwd=1.5,lty=3)
                    a=a+labs(x="Filters (minimum periodicity)",y="Estimated Growth Impact\n(% per Degree)",col="Impact Model")
                    a=a+scale_color_manual(values=c("#d3818c","#7375a4")) + 
                    ggtitle("Drift-Adjusted") #+ ylim(-10,5)
                    a


                    b=ggplot(simsfiltdat,aes(x=factor(periodsregationPeriod),y=MeanEffect*100,col=ImpactType,group=ImpactType))+geom_line(lwd=1.25)
                    b=b+geom_point(size=4)+geom_errorbar(aes(ymin=(MeanEffect-1.96*SDEffect)*100,ymax=(MeanEffect+1.96*SDEffect)*100),width=0.1,lwd=1.25)
                    b=b+geom_hline(yintercept=0,lwd=1.5,lty=3)
                    b=b+labs(x="Filters (minimum periodicity)",y="Estimated Growth Impact\n(% per Degree)",col="Impact Model")
                    b=b+scale_color_manual(values=c("#d3818c","#7375a4")) + ggtitle("Impacts estimation using low pass filters") #+ ylim(-10,5)
                    b
                # Comparing quadratic (end)
                
                ggarrange(b,a,nrow=2,ncol=1)
    ## 2.4. Perform simulation - WITH DRIFT ADJUSTMENT (end)

    ## 2.4. Perform simulation - Combaination Gammas and Betas (start)


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
                    b=b+scale_color_manual(values=c("#d3818c","#7375a4")) + ggtitle("Impacts estimation using low pass filters") #+ ylim(-10,5)
                    b
                # Comparing quadratic (end)
                
                ggarrange(b,a,nrow=2,ncol=1)
    ## 2.4. Perform simulation - WITH DRIFT ADJUSTMENT (end)

    
    #2.5. Perform Simulation using FFT (start)
        # Useful functions (start)
            convert.fft <- function(cs, sample.rate=1) {
            cs <- cs / length(cs) # normalize

            distance.center <- function(c)signif( Mod(c),        4)
            angle           <- function(c)signif( 180*Arg(c)/pi, 3)
            
            df <- data.frame(cycle    = 0:(length(cs)-1),
                            freq     = 0:(length(cs)-1) * sample.rate / length(cs),
                            strength = sapply(cs, distance.center),
                            delay    = sapply(cs, angle))
            df
            }

            get.trajectory <- function(X.k,ts,acq.freq) {
                # returns the x.n time series for a given time sequence (ts) and
                # a vector with the amount of frequencies k in the signal (X.k)
                
            N   <- length(ts)
            i   <- complex(real = 0, imaginary = 1)
            x.n <- rep(0,N)           # create vector to keep the trajectory
            ks  <- 0:(length(X.k)-1)
            
            for(n in 0:(N-1)) {       # compute each time point x_n based on freqs X.k
                x.n[n+1] <- sum(X.k * exp(i*2*pi*ks*n/N)) / N
            }
            
            x.n * acq.freq 
            }
        
            plot.harmonic <- function(Xk, i, ts, acq.freq, color="red") {
            Xk.h <- rep(0,length(Xk))
            Xk.h[i+1] <- Xk[i+1] # i-th harmonic
            harmonic.trajectory <- get.trajectory(Xk.h, ts, acq.freq=acq.freq)
            points(ts, harmonic.trajectory, type="l", col=color)
            }

            plot.frequency.spectrum <- function(X.k, xlimits=c(0,length(X.k))) {
            plot.data  <- cbind(0:(length(X.k)-1), Mod(X.k))


            plot.data[2:length(X.k),2] <- 2*plot.data[2:length(X.k),2] 
            
            plot(plot.data, t="h", lwd=2, main="", 
                xlab="Frequency (Hz)", ylab="Strength", 
                xlim=xlimits, ylim=c(0,max(Mod(plot.data[,2]))))
            }

            plot.show <- function(trajectory, time=1, harmonics=-1, plot.freq=FALSE) {

            acq.freq <- length(trajectory)/time      # data acquisition frequency (Hz)
            ts  <- seq(0,time-1/acq.freq,1/acq.freq) # vector of sampling time-points (s) 
            
            X.k <- fft(trajectory)
            x.n <- get.trajectory(X.k,ts, acq.freq=acq.freq) / acq.freq
            
            if (plot.freq)
                plot.frequency.spectrum(X.k)
            
            max.y <- ceiling(1.5*max(Mod(x.n)))
            
            if (harmonics[1]==-1) {
                min.y <- floor(min(Mod(x.n)))-1
            } else {
                min.y <- ceiling(-1.5*max(Mod(x.n)))
            }
            
            plot(ts,x.n, type="l",ylim=c(min.y,max.y))
            abline(h=min.y:max.y,v=0:time,lty=3)
            points(ts,trajectory,pch=19,col="red")  # the data points we know
            
            if (harmonics[1]>-1) {
                for(i in 0:length(harmonics)) {
                plot.harmonic(X.k, harmonics[i], ts, acq.freq, color=i+1)
                }
            }
            }

            convert.fft <- function(cs, sample.rate=1) {
            cs <- cs / length(cs) # normalize

            distance.center <- function(c)signif( Mod(c),        4)
            angle           <- function(c)signif( 180*Arg(c)/pi, 3)
            
            df <- data.frame(cycle    = 0:(length(cs)-1),
                            freq     = 0:(length(cs)-1) * sample.rate / length(cs),
                            strength = sapply(cs, distance.center),
                            delay    = sapply(cs, angle))
            df
            }

        # Useful Functions (end)

        #Simulation FFT (start)
            time=350 #or 350 years
            basegr=0.01 #2% per year baseline growth rate
            start=100
            baseline=rep(basegr,time)
            coef=-0.05 #effect size - change in growth per degree warming()
            coef2=-0.05
            growthsd=0.005 #standard deviation of growth variability unexplained by temperature
            periods <- c(0,345,346,347,348)
            nsims=500
            sims=array(dim=c(nsims,length(periods),2))
            sims_adjusted=array(dim=c(nsims,length(periods),2))

            acq.freq <- 1    # data acquisition frequency (Hz)
            ts  <- seq(0,time-1/acq.freq,1/acq.freq) # vector of sampling time-points (s) 
      
            for(i in 1:nsims){
                randomtemp=Re(randomts(gtemp))[1:time]
                #randomtemp=Re(randomts(ustemp_lmr))[1:time]
                randomgrowth_g=basegr #growth impact model
                randomgrowth_l=basegr #levels impact model
                for(j in 2:time){
                #quadratic
                #randomgrowth_g=c(randomgrowth_g, basegr+(randomtemp[j]*coef)+(coef2*randomtemp[j]^2)+rnorm(1,mean=0,sd=growthsd))
                #randomgrowth_l=c(randomgrowth_l, basegr+(randomtemp[j]-randomtemp[j-1])*coef+coef2*(randomtemp[j]^2-randomtemp[j-1]^2)+rnorm(1,mean=0,sd=growthsd))
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
                        range <- length(randomtemp) - periods[k]
                        X.temp <- fft(randomtemp)
                        Xtemp.h <- rep(0,length(X.temp))
                        #Xtemp.h[(k+1)] <- X.temp[(k+1)] #  i-th harmonic
                        Xtemp.h[1:range] <- X.temp[1:range] # up to range-th harmonic
                        tempts <- Re(get.trajectory(Xtemp.h, ts, acq.freq))

                        X.g <- fft(randomgrowth_g)
                        Xg.h <- rep(0,length(X.g))
                        Xg.h[1:range] <- X.g[1:range] # up to i-th harmonic
                        growth <- Re(get.trajectory(Xg.h, ts, acq.freq))

                        X.l <- fft(randomgrowth_l)
                        Xl.h <- rep(0,length(X.l))
                        Xl.h[1:range] <- X.l[1:range] # up to i-th harmonic
                        level <- Re(get.trajectory(Xl.h, ts, acq.freq))

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
        # Simulaiton FFT (end)
        
        #Plot Test
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
        #Plot Test
    #2.5 Perform simulation using FFT (end)

    
## 2. Simulation - Figures 1 and 2 (end)

## 3. Empirical analysis - Figure 3 and Table 1 (start)

    ## 3.1. Country-level regressions (start)
        dataset <- c("wb","barro","mad") 
        datasetweather <- c("LMR","UDel")
        periods <- c(0,3,5,10,15)
        #periods <- c(0,10,20,25,27) Uncomment to get Supp Fig 2
        fullmods_filter=array(dim=c(length(countries),2,length(periods),length(dataset),length(datasetweather)))
        fullmods_filter_var=array(dim=c(length(countries),2,length(periods),length(dataset),length(datasetweather)))
        fullmods_filter_p=array(dim=c(length(countries),2,length(periods),length(dataset),length(datasetweather)))
        panel_data <- data.frame(years = integer(),temp = double(), growth = double(), 
            preci = double(), countrycode = factor(), climdata = character(),
            econdata = character(), filter = character(), meant = double())
        for (mm in 1:length(datasetweather)){
            tempname <- paste(datasetweather[mm],"_pop_temp",sep="")
            preciname <- paste(datasetweather[mm],"_pop_preci",sep="")
            for (jj in (1:length(dataset))){
                DATA <- get(dataset[jj])
                countries=unique(factor(DATA$countrycode))
                fullmods=array(dim=c(length(countries),2,length(periods)))
                for(k in 1:length(periods)){
                for(i in 1:length(countries)){
                    dat=DATA[which(DATA$countrycode==countries[i]),
                        which(colnames(DATA)%in%c("countrycode","year","growth",tempname,preciname))] 
                    dat <- dat[is.finite(dat$growth),]
                        gonext <- 0
                    for (nn in (1:(dim(dat)[2]))){
                        if(sum(is.na(dat[,nn]))==dim(dat)[1]){
                            #print(paste(countries[i],i,"missing complete column"))
                            gonext <- 1 }
                        
                    }
                    if(gonext==1){
                        gonext <- 0
                        next}
                    if(sum(complete.cases(dat))<2){
                            #print(paste(countries[i],i,"not complete cases"))
                            next
                    }
                    names(dat) <- c("countrycode","year","growth","temp","preci")
                    mt <- lm(temp~year+I(year^2), data = dat)
                    meanT <- mean(dat$temp, na.rm = TRUE)
                    t <- resid(mt)    
                    t <- timeSeries::interpNA(t, method = "linear")
                    t <- timeSeries::removeNA(t)
                    if(sum(complete.cases(t))<(2*(periods[k]+3))){
                            #print(paste(countries[i],i,"not enough data for this filter"))
                            next}
                    if(k==1){
                        tempts <- t
                        ratio <- 1
                    } else{
                    tempts <- pass.filt(t, W=periods[k], type="low", method="Butterworth")
                    x <- seq(1:length(tempts))
                    #demeaned_t <- lm(t~x)$residuals
                    #demeaned_tempts <- lm(tempts~x)$residuals
                    ratio <- median(demeaned_t/demeaned_tempts)
                     
                    }
                     
                    temp <- data.frame(year = 1:length(tempts), temp= unclass(tempts))
                    if(sum(complete.cases(temp))<4){next}
                    mp <- lm(preci~year+I(year^2), data = dat)
                    p <- resid(mp)
                    p <- interpNA(p, method = "linear")
                    p <- removeNA(p)
                    if(k==1){
                        precits <- p
                    } else{
                    precits <- pass.filt(p, W=periods[k], type="low", method="Butterworth")
                    }
                    preci <- data.frame(year = 1:length(precits), preci= unclass(precits))
                    mg <- lm(growth~year+I(year^2), data = dat)
                    g <- resid(mg)
                    g <- interpNA(g, method = "linear")
                    g <- removeNA(g)
                    growth <- data.frame(year = 1:length(g), growth = unclass(g))
                    filterdata <- merge(temp,growth, by = "year")
                    filterdata <- merge(filterdata,preci, by = "year")
                    names(filterdata) <- c("years","temp","growth","preci")
                    if(k==1){
                        filterdata$templag =c(NA,filterdata$temp[1:(dim(filterdata)[1]-1)])
                        mod_gfilterdata=lm(growth~temp+preci+templag,data=filterdata)
                    } else{
                        mod_gfilterdata=lm(growth~temp+preci,data=filterdata)
                        # Uncomment to include templag
                    }

                    filterdata$years <- mt$model[[2]]
                    filterdata$countrycode <- rep(dat$countrycode[1],dim(filterdata)[1])
                    filterdata$climdata <- rep(datasetweather[mm],dim(filterdata)[1])
                    filterdata$econdata <- rep(dataset[jj],dim(filterdata)[1])
                    filterdata$filter <- rep(paste(periods[k],sep="-"),dim(filterdata)[1])
                    filterdata$meant <- rep(meanT,dim(filterdata)[1])
                    panel_data <- bind_rows(panel_data,filterdata)
                    fullmods_filter[i,,k,jj,mm]=c(summary(mod_gfilterdata)$coefficients[2,1:2]/ratio)
                    fullmods_filter_var[i,,k,jj,mm]=c(summary(mod_gfilterdata)$coefficients[2,1]/ratio,vcov(mod_gfilterdata)[2,2]/(ratio^2))
                    fullmods_filter_p[i,,k,jj,mm]=summary(mod_gfilterdata)$coefficients[3,1:2] #If usd should be divided by the ratio of precipietation timeseries
                    }
                    
                
                }
            }
        }

        ranges <- paste(periods,sep='-')
        length(levels(wb$countrycode))
        countries <- unique(factor(wb$countrycode))
        countries_wb <- unlist(lapply(countries, as.character))
        dimnames(fullmods_filter)=list(countries,c("Estimate","StandardError"),ranges,dataset,datasetweather)
        fullmods_filterm <- melt(fullmods_filter)
        countriesbarro <- unique(factor(barro$countrycode))
        countries_barro <- unlist(lapply(countriesbarro , as.character))
        countries_barro_extended <- c(countries_barro,countries_wb[(length(countries_barro)+1):length(countries_wb)*NA])
        countriesmad <- unique(factor(mad$countrycode))
        countries_mad <- unlist(lapply(countriesmad , as.character))
        countries_mad_extended <- c(countries_mad,countries_wb[(length(countries_mad)+1):length(countries_wb)*NA])
        countries_wbrep <- rep(countries_wb,2*5)
        countries_barrorep <- rep(countries_barro_extended,2*5)
        countries_madrep <- rep(countries_mad_extended,2*5)
        names1 <- c(countries_wbrep,countries_barrorep,countries_madrep)
        names2 <- c(names1,names1)
        fullmods_filterm$Var1 <- names2
        fullmods_filter=dcast(fullmods_filterm,Var1+Var3+Var4+Var5~Var2, fun = mean)
        fullmods_filter$Estimate[which(is.infinite(fullmods_filter$StandardError))]=NA
        names(fullmods_filter) <- c("countrycode","frequencies","econdata","climdata","Estimate","StandardError")


        dimnames(fullmods_filter_var)=list(countries,c("Estimate","Var"),ranges,dataset,datasetweather)
        fullmods_filterm <- melt(fullmods_filter_var)
        fullmods_filterm$Var1 <- names2
        fullmods_filter_v=dcast(fullmods_filterm,Var1+Var3+Var4+Var5~Var2, fun = mean)
        fullmods_filter_v$Estimate[which(is.infinite(fullmods_filter_v$Var))]=NA
        names(fullmods_filter_v) <- c("countrycode","frequencies","econdata","climdata","Estimate","Variance")


    ## 3.1. Country-level regressions (end)
        
    ## 3.1b Panel Regression (start)
        
        coefT <- "temp"
        x <- c(0,1,2,3)
        xmat <- cbind(1, 2*x)
            
        p1 <- panel_data[which(panel_data$climdata=="UDel" & panel_data$econdata=="wb"  & panel_data$filter=="0"),]
        felm_panel1 <- felm(growth ~ temp+I(temp^2)+preci+I(preci^2)+years+I(years^2)|countrycode|0|countrycode, data =p1)
        # Getting marginal effect of an additional degree (start)
            #summary(felm_panel1)
            model <- felm_panel1
            Sigma <- vcov(model)
            sigma = Sigma[c(start1:end1),c(start1:end1)]
            start1 <- which(names(coef(model))==coefT)
            end1 <- which(names(coef(model))==paste("I(",coefT,"^2)",sep=""))
            beta.hat <- coef(model)[c(start1:end1)]
            dgdt <- colSums(beta.hat*t(xmat)) 
            ci1 <- dgdt + 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
            ci2 <- dgdt- 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
        # Getting marginal effect of an additional degree (start)
            
        p2 <- panel_data[which(panel_data$climdata=="UDel" & panel_data$econdata=="wb"  & panel_data$filter=="3"),]
        felm_panel2 <- felm(growth ~ temp+I(temp^2)+preci+I(preci^2)+years+I(years^2)|countrycode|0|countrycode, data =p2)
        # Getting marginal effect of an additional degree (start)
            model <- felm_panel2
            Sigma <- vcov(model)
            sigma = Sigma[c(start1:end1),c(start1:end1)]
            start1 <- which(names(coef(model))==coefT)
            end1 <- which(names(coef(model))==paste("I(",coefT,"^2)",sep=""))
            beta.hat <- coef(model)[c(start1:end1)]
            dgdt <- colSums(beta.hat*t(xmat)) 
            ci1 <- dgdt + 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
            ci2 <- dgdt- 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
        # Getting marginal effect of an additional degree (start)
         
        p3 <- panel_data[which(panel_data$climdata=="UDel" & panel_data$econdata=="wb"  & panel_data$filter=="5"),]
        felm_panel3 <- felm(growth ~ temp+I(temp^2)+preci+I(preci^2)+years+I(years^2)|countrycode|0|countrycode, data =p3)
        # Getting marginal effect of an additional degree (start)
            model <- felm_panel3
            Sigma <- vcov(model)
            sigma = Sigma[c(start1:end1),c(start1:end1)]
            start1 <- which(names(coef(model))==coefT)
            end1 <- which(names(coef(model))==paste("I(",coefT,"^2)",sep=""))
            beta.hat <- coef(model)[c(start1:end1)]
            dgdt <- colSums(beta.hat*t(xmat)) 
            ci1 <- dgdt + 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
            ci2 <- dgdt- 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
        # Getting marginal effect of an additional degree (start)
         
        p4 <- panel_data[which(panel_data$climdata=="UDel" & panel_data$econdata=="wb"  & panel_data$filter=="10"),]
        felm_panel4 <- felm(growth ~ temp+I(temp^2)+preci+I(preci^2)+years+I(years^2)|countrycode|0|countrycode, data =p4)
        # Getting marginal effect of an additional degree (start)
            model <- felm_panel4
            Sigma <- vcov(model)
            sigma = Sigma[c(start1:end1),c(start1:end1)]
            start1 <- which(names(coef(model))==coefT)
            end1 <- which(names(coef(model))==paste("I(",coefT,"^2)",sep=""))
            beta.hat <- coef(model)[c(start1:end1)]
            dgdt <- colSums(beta.hat*t(xmat)) 
            ci1 <- dgdt + 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
            ci2 <- dgdt- 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
        # Getting marginal effect of an additional degree (start)
         
        p5 <- panel_data[which(panel_data$climdata=="UDel" & panel_data$econdata=="wb"  & panel_data$filter=="15"),]
        felm_panel5 <- felm(growth ~ temp+I(temp^2)+preci+I(preci^2)+years+I(years^2)|countrycode|0|countrycode, data =p5)
        # Getting marginal effect of an additional degree (start)
            model <- felm_panel5
            Sigma <- vcov(model)
            sigma = Sigma[c(start1:end1),c(start1:end1)]
            start1 <- which(names(coef(model))==coefT)
            end1 <- which(names(coef(model))==paste("I(",coefT,"^2)",sep=""))
            beta.hat <- coef(model)[c(start1:end1)]
            dgdt <- colSums(beta.hat*t(xmat)) 
            ci1 <- dgdt + 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
            ci2 <- dgdt- 1.96*sqrt(diag((xmat %*% sigma) %*% t(xmat)))
        # Getting marginal effect of an additional degree (start)
         
        
        stargazer(felm_panel1,felm_panel2,felm_panel3,felm_panel4,felm_panel5,type="text")

    ## 3.1b Panel Regression (end)

    ## 3.2. categorizing - Plotting Figure 3  (start)

        

        #   Original Categorization (start)
            fmod_fft$high95 <- fmod_fft$Estimate + fmod_fft$StandardError*1.64
            fmod_fft$low95 <- fmod_fft$Estimate - fmod_fft$StandardError*1.64
            uncertain <- c(which(fmod_fft$high95>0 & fmod_fft$low95 <0))
            positive_Constant <- 0
            positive_Intensifying <- 0
            positive_Converging <- 0
            fmod_fft$absestimate <- abs(fmod_fft$Estimate)
            for (i in 1:numcountries){
                        if(is.null(fmod_fft$absestimate[1+5*(i-1)] )){next}
                        if(!is.na(fmod_fft$absestimate[5+5*(i-1)])){
                            lastfreq <- 5
                            }else if(!is.na(fmod_fft$absestimate[4+5*(i-1)])){
                                lastfreq <- 4
                            }else if(!is.na(fmod_fft$absestimate[3+5*(i-1)])){
                                lastfreq <- 3
                            }else if(!is.na(fmod_fft$absestimate[2+5*(i-1)])){
                                lastfreq <- 2
                            }else{next}
                        if(fmod_fft$absestimate[1+5*(i-1)]>0){ #all because absolute value
                            if(fmod_fft$absestimate[lastfreq+5*(i-1)]>0){ #all because using abs value
                                if((fmod_fft$absestimate[lastfreq+5*(i-1)] )>(fmod_fft$absestimate[1+5*(i-1)] /2)){ #larger than half the first estimate
                                    if((fmod_fft$absestimate[lastfreq+5*(i-1)] )>(fmod_fft$absestimate[1+5*(i-1)] *1.5)){
                                        positive_Intensifying <- c(positive_Intensifying,(1+5*(i-1)):(5+5*(i-1)))
                                    } else{
                                        if((fmod_fft$Estimate[lastfreq+5*(i-1)]*fmod_fft$Estimate[1+5*(i-1)])<0){
                                            positive_Converging <- c(positive_Converging,(1+5*(i-1)):(5+5*(i-1)))
                                        } else{
                                        #same sign? if not, is decreasinfg, else:
                                        positive_Constant <- c(positive_Constant,(1+5*(i-1)):(5+5*(i-1)))}} }
                                        else{
                                            positive_Converging <- c(positive_Converging,(1+5*(i-1)):(5+5*(i-1)))
                                        }}}}
                fmod_fft$significant <- 1
                fmod_fft$significant[uncertain]<- 0
                #is the lowest frequency estimate significant?
                for (i in 1:numcountries){
                    if(!is.na(fmod_fft$absestimate[5+5*(i-1)])){
                            lastfreq <- 5
                            }else if(!is.na(fmod_fft$absestimate[4+5*(i-1)])){
                                lastfreq <- 4
                            }else if(!is.na(fmod_fft$absestimate[3+5*(i-1)])){
                                lastfreq <- 3
                            }else if(!is.na(fmod_fft$absestimate[2+5*(i-1)])){
                                lastfreq <- 2
                            }else{next}
                    if((lastfreq+5*(i-1)) %in% uncertain){ 
                        fmod_fft$significant[1+5*(i-1)] <- 0
                        fmod_fft$significant[2+5*(i-1)] <- 0
                        fmod_fft$significant[3+5*(i-1)] <- 0
                        fmod_fft$significant[4+5*(i-1)] <- 0
                        fmod_fft$significant[5+5*(i-1)] <- 0
                    } else {
                        fmod_fft$significant[1+5*(i-1)] <- 1
                        fmod_fft$significant[2+5*(i-1)] <- 1
                        fmod_fft$significant[3+5*(i-1)] <- 1
                        fmod_fft$significant[4+5*(i-1)] <- 1
                        fmod_fft$significant[5+5*(i-1)] <- 1
                    }
                }
                filt_names <- c("Unfiltered","3 years","5 years", "10 years", "15 years")
                #filt_names <- c("Unfiltered","10 years","15 years", "20 years", "25 years") #uncomment to get supp Fig 2
                fmod_fft$filters <- rep(filt_names,217)
                fmod_fft$filters <- factor(fmod_fft$filters, levels = filt_names)            
                fmod_fft$category <-"other"
                fmod_fft$category[positive_Constant] <- "Constant"
                fmod_fft$category[positive_Intensifying] <- "Intensifying"
                fmod_fft$category[positive_Converging] <- "Converging"
                

                fpc <- fmod_fft[positive_Constant,]
                fpi <- fmod_fft[positive_Intensifying,]
                fpd <- fmod_fft[positive_Converging,]
                fpo <- fmod_fft[fmod_fft$category=="other",] #no estimates
                    
                #Removing outliers
                fpi <- fpi[fpi$countrycode!="SSD",]
                fpi <- fpi[fpi$countrycode!="SLE",]
                fpi <- fpi[fpi$countrycode!="GNQ",]
                #fpc <- fpc[fpc$countrycode!="LBY",]
                fpi <- fpi[fpi$countrycode!="KOR",]
                fpi <- fpi[fpi$countrycode!=fpi$countrycode[which(fpi$Estimate==min(fpi$Estimate,na.rm=TRUE))],]
                
                # Figure 3 (start)
                    plot_fpi <- ggplot(data=fpi,aes(x=filters,y=Estimate*100, group = countrycode,color=factor(significant)))+
                    geom_line()+
                    scale_colour_discrete(guide = 'none') +
                    theme_bw()+
                    geom_hline(yintercept=0,lty=2)+
                    geom_dl(data=fpi[fpi$significant==1,],aes(label = countrycode), method = list(dl.combine("last.points")), cex = 0.9)+
                    labs(title="Intensifying")  + xlab("") + ylab("Estimated Effect of \n 1 Degree Warming on Growth (pp)") +
                    theme(axis.line = element_line(colour = "black"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank()) 



                    
                    plot_fpc <- ggplot(data=fpc,aes(x=filters,y=Estimate*100, group = countrycode,color=factor(significant)))+
                    geom_line()+
                    scale_colour_discrete(name="Low-frequency estimate significantly different from zero") +
                    theme_bw()+
                    geom_hline(yintercept=0,lty=2)+
                    geom_dl(data=fpc[fpc$significant==1,],aes(label = countrycode), method = list(dl.combine("last.points")), cex = 0.9)+
                    labs(title="Constant")  + xlab("Minimum Periodicity after Filtering") + ylab("") +
                    theme(axis.line = element_line(colour = "black"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank()) 



                    #fpd <- fpd[fpd$countrycode!="LBY",]
                    plot_fpd <-ggplot(data=fpd,aes(x=filters,y=Estimate*100, group = countrycode,color=factor(significant)))+
                    geom_line()+
                    scale_colour_discrete(guide = 'none') +
                    theme_bw()+
                    geom_hline(yintercept=0,lty=2)+
                    geom_dl(data=fpd[fpd$significant==1,],aes(label = countrycode), method = list(dl.combine("last.points")), cex = 0.9)+
                    labs(title="Converging")  + xlab("") + ylab("")+
                    theme(axis.line = element_line(colour = "black"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank()) 

                    d <- ggarrange(plot_fpi,plot_fpc,plot_fpd,ncol=3,common.legend=TRUE) 
                    d
                    table(fmod_fft$category)/5
                    #ggsave('Fig3_cat.png',dpi=500) 
                    
                # Figure 3 (start)

                #Plotting some outliers 
                    lby <- ggplot(data=wb[wb$c
                    ountrycode=="LBY",], aes(x=year, y=growth))+geom_line()+theme_bw()+ggtitle("Economic growth in Libya") #ilitary intervention in 2011
                    gnq <- ggplot(data=wb[wb$countrycode=="GNQ",], aes(x=year, y=growth))+geom_line()+theme_bw()+ggtitle("Economic growth in Equatorial Guinea") #In 1995 Mobil discovered oil in Guinea
                    ssd <- ggplot(data=wb[wb$countrycode=="SSD",], aes(x=year, y=growth))+
                    geom_line()+theme_bw()+ggtitle("Economic growth in South Sudan") #In 1995 Mobil discovered oil in Guinea              
                    ukr <- ggplot(data=wb[wb$countrycode=="UKR",], aes(x=year, y=growth))+
                    geom_line()+theme_bw()+ggtitle("Economic growth in Ukraine") #In 1995 Mobil discovered oil in Guinea
                    ggarrange(d,ggarrange(gnq,lby,ncol=2),nrow=2)
                    #ggsave('Fig_pathways_outliers.png',dpi=500) #Figure 3
        #   Original Categorization (end)

        # Categorizing significant intensifying (start)
            fmod_fft <- fullmods_filter
            fmod_fft <- fmod_fft[fmod_fft$econdata=="wb",]
            fmod_fft <- fmod_fft[fmod_fft$climdata=="UDel",]
                fmod_fft$high95 <- fmod_fft$Estimate + fmod_fft$StandardError*1.64
                fmod_fft$low95 <- fmod_fft$Estimate - fmod_fft$StandardError*1.64
                fmod_fft$absestimate <- abs(fmod_fft$Estimate)
                fmod_fft$significant <- 0
                levels <- 0
                growth <- 0              
                uncertain <- c(which(fmod_fft$high95>0 & fmod_fft$low95 <0))
                `%notin%` <- Negate(`%in%`)
                for (i in 1:numcountries){
                    if(!is.na(fmod_fft$absestimate[5+5*(i-1)])){
                            lastfreq <- 5
                            }else if(!is.na(fmod_fft$absestimate[4+5*(i-1)])){
                                lastfreq <- 4
                            }else if(!is.na(fmod_fft$absestimate[3+5*(i-1)])){
                                lastfreq <- 3
                            }else if(!is.na(fmod_fft$absestimate[2+5*(i-1)])){
                                lastfreq <- 2
                            }else{next}
                     if((lastfreq+5*(i-1)) %in% uncertain){

                        levels <- c(levels,(1+5*(i-1)):(5+5*(i-1))) 
                        
                    } else if((lastfreq+5*(i-1)) %notin% uncertain) {
                        growth <- c(growth,(1+5*(i-1)):(5+5*(i-1))) 
                        
                    }
                    if((1+5*(i-1)) %notin% uncertain){
                        fmod_fft$significant[(1+5*(i-1)):(5+5*(i-1))] <- 1     
                                       }}
                filt_names <- c("Unfiltered","3 years","5 years", "10 years", "15 years")
                #filt_names <- c("Unfiltered","10 years","15 years", "20 years", "25 years") #uncomment to get supp Fig 2
                fmod_fft$filters <- rep(filt_names,217)
                fmod_fft$filters <- factor(fmod_fft$filters, levels = filt_names)            
                fmod_fft$category <-"other"
                fmod_fft$category[growth] <- "growth"
                fmod_fft$category[levels] <- "levels"
                
                

                fg <- fmod_fft[growth,]
                fl <- fmod_fft[levels,]
                    
                #Removing outliers
                fg <- fg[fg$countrycode!="SSD",]
                fg <- fg[fg$countrycode!="SLE",]
                fg <- fg[fg$countrycode!="GNQ",]
                #fpc <- fpc[fpc$countrycode!="LBY",]
                fg <- fg[fg$countrycode!="KOR",]
                fg <- fg[fg$countrycode!=fg$countrycode[which(fg$Estimate==min(fg$Estimate,na.rm=TRUE))],]
                fl <- fl[fl$countrycode!=fl$countrycode[which(fl$Estimate==min(fl$Estimate,na.rm=TRUE))],]
                fl <- fl[fl$countrycode!=fl$countrycode[which(fl$Estimate==min(fl$Estimate,na.rm=TRUE))],]
                
                # Figure 3 (start)
                    plot_fg <- ggplot(data=fg,aes(x=filters,y=Estimate*100, group = countrycode,color=factor(significant)))+
                    geom_line()+
                    scale_colour_discrete(name="Unfiltered estimate significantly different from zero") +
                    theme_bw() + xlab("Minimum Periodicity after Filtering")+
                    geom_hline(yintercept=0,lty=2)+
                    #geom_dl(data=fg[fg$significant==1,],aes(label = countrycode), method = list(dl.combine("last.points")), cex = 0.9)+
                    labs(title="Evidence of Growth Effects \n (Low-freq estimate statistically different from zero)")  + ylab("Estimated Effect of \n 1 Degree Warming on Growth (pp)") +
                    theme(axis.line = element_line(colour = "black"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(),legend.position="bottom") 



                    
                    plot_fl <- ggplot(data=fl,aes(x=filters,y=Estimate*100, group = countrycode,color=factor(significant)))+
                    geom_line()+
                    scale_colour_discrete(name="Unfiltered estimate significantly different from zero") +
                    theme_bw()+
                    geom_hline(yintercept=0,lty=2)+
                    geom_dl(data=fpc[fpc$significant==1,],aes(label = countrycode), method = list(dl.combine("last.points")), cex = 0.9)+
                    labs(title="No evidence of growth effects")  + xlab("Minimum Periodicity after Filtering") + ylab("") +
                    theme(axis.line = element_line(colour = "black"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank()) 


                    d <- ggarrange(plot_fg,plot_fl,ncol=2,common.legend=TRUE,legend="bottom") 
                    d
                    table(fmod_fft$category)/5
                    #ggsave('Fig3_cat.png',dpi=500) 
                    
                    # Mapping estimates
                        glimpse(fmod_fft)
                        unique(fmod_fft$countrycode)
                        install.packages('rnaturalearthdata')
                        library('rnaturalearthdata')
                        world <- ne_countries(scale = "medium", returnclass = "sf")
                        
                        
                        fmod_fft$iso_a3 <- fmod_fft$countrycode
                        fmod_fft$test1 <- NA
                        fmod_fft$test1[fmod_fft$category=="growth"] <- 1
                        fmod_fft_map <- merge(world,fmod_fft,by="iso_a3")
                        library(RColorBrewer)
                            cols <- c("1" = "#8dd3c7", "0" = "#ffffb3")
                            map_growth <- ggplot(data = fmod_fft_map) +
                                geom_sf(aes(fill = factor(test1)))+
                                scale_fill_manual(name = "Growth Effects",
                                                    labels = c("Detected","Not detected"),
                                values=cols)+
                                theme(legend.position="bottom")+
                                ggtitle("Location of Growth Effects")

                            
                           
                            growtheff <- fmod_fft
                            cg <- unique(growtheff$countrycode)
                            for (i in 1:length(cg)){
                                if(is.na(growtheff$Estimate[5+5*(i-1)])){
                                    growtheff$Estimate[5+5*(i-1)] <- growtheff$Estimate[4+5*(i-1)]
                                    
                                if(is.na(growtheff$Estimate[4+5*(i-1)])){
                                    growtheff$Estimate[5+5*(i-1)] <- growtheff$Estimate[3+5*(i-1)]
                                    
                                if(is.na(growtheff$Estimate[3+5*(i-1)])){
                                    growtheff$Estimate[5+5*(i-1)] <- growtheff$Estimate[2+5*(i-1)]
                                    }}}
                            }
                            growtheff$Estimate
                            growtheff <- growtheff[which(growtheff$test1==1),]

                            growtheff <- growtheff[which(growtheff$frequencies==max(growtheff$frequencies)),]
                            growtheff <- merge(world,growtheff,by="iso_a3")

                            map_significantgrowth <- ggplot(data = fmod_fft_map) +
                                geom_sf(fill=NA)+
                                geom_sf(data=growtheff,aes(fill = Estimate*100))+
                                scale_fill_gradient2(
                                    name = "Estimated impact \n (% per Degree)",
                                    low = "red",
                                    mid = "white",
                                    high = "#00BFC4",
                                    midpoint = 0,
                                    space = "Lab",
                                    na.value = "grey50",
                                    guide = "colourbar",
                                    aesthetics = "fill")+
                                    ggtitle("Detected growth effects")+
                                theme(legend.position="bottom")
                            map_significantgrowth    
                                
                                
                             

                        ggarrange(plot_fg,map_growth,ncol=1,nrows=2)

                    # Mapping estimates
                # Figure 3 (start)

                #Plotting some outliers 
                    lby <- ggplot(data=wb[wb$c
                    ountrycode=="LBY",], aes(x=year, y=growth))+geom_line()+theme_bw()+ggtitle("Economic growth in Libya") #ilitary intervention in 2011
                    gnq <- ggplot(data=wb[wb$countrycode=="GNQ",], aes(x=year, y=growth))+geom_line()+theme_bw()+ggtitle("Economic growth in Equatorial Guinea") #In 1995 Mobil discovered oil in Guinea
                    ssd <- ggplot(data=wb[wb$countrycode=="SSD",], aes(x=year, y=growth))+
                    geom_line()+theme_bw()+ggtitle("Economic growth in South Sudan") #In 1995 Mobil discovered oil in Guinea              
                    ukr <- ggplot(data=wb[wb$countrycode=="UKR",], aes(x=year, y=growth))+
                    geom_line()+theme_bw()+ggtitle("Economic growth in Ukraine") #In 1995 Mobil discovered oil in Guinea
                    ggarrange(d,ggarrange(gnq,lby,ncol=2),nrow=2)
                    #ggsave('Fig_pathways_outliers.png',dpi=500) #Figure 3
        # Categorizing significant intensifying (end)

        # Categorizing statistically different estimates (start)
            glimpse(fullmods_filter_v)
            glimpse(fullmods_filter)
            fullmods_filter$Variance <- fullmods_filter_v$Variance
            fmod_fft <- fullmods_filter
            fmod_fft <- fmod_fft[fmod_fft$econdata=="wb",]
            fmod_fft <- fmod_fft[fmod_fft$climdata=="UDel",]
                fmod_fft$high95 <- fmod_fft$Estimate + fmod_fft$StandardError*1.64
                fmod_fft$low95 <- fmod_fft$Estimate - fmod_fft$StandardError*1.64
                #fmod_fft$absestimate <- abs(fmod_fft$Estimate)
                fmod_fft$lowsignificant <- 0
                fmod_fft$unfilteredsignificant <- 0
                              
                uncertain <- c(which(fmod_fft$high95>0 & fmod_fft$low95 <0))
                `%notin%` <- Negate(`%in%`)
                fmod_fft$sign <- rep(0,dim(fmod_fft)[1])
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

                converging <- 0
                intensifying <- 0
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
                    theta <- abs(fmod_fft$Estimate[1+5*(i-1)]) - abs(fmod_fft$Estimate[lastfreq+5*(i-1)])
                    var <- fmod_fft$Variance[lastfreq+5*(i-1)] + fmod_fft$Variance[1+5*(i-1)]
                    conf95 <-  (var^0.5)*1.65 #one-tail  95%
                    if(fmod_fft$sign[lastfreq+5*(i-1)]==1){
                        if(theta-conf95 >0 ){
                            converging <- c(converging,(1+5*(i-1)):(5+5*(i-1))) 
                        }
                        if(theta+conf95 < 0 ){
                            intensifying <- c(intensifying,(1+5*(i-1)):(5+5*(i-1))) 
                        }
                        if((lastfreq+5*(i-1)) %notin% uncertain){
                        fmod_fft$lowsignificant[(1+5*(i-1)):(5+5*(i-1))] <- 1     
                                       }
                        if((1+5*(i-1)) %notin% uncertain){
                            fmod_fft$unfilteredsignificant[(1+5*(i-1)):(5+5*(i-1))] <- 1     
                                       }
                    }
                    }
                    
                filt_names <- c("Unfiltered","3 years","5 years", "10 years", "15 years")
                #filt_names <- c("Unfiltered","10 years","15 years", "20 years", "25 years") #uncomment to get supp Fig 2
                fmod_fft$filters <- rep(filt_names,217)
                fmod_fft$filters <- factor(fmod_fft$filters, levels = filt_names)            
                fmod_fft$category <- "Undefined"
                fmod_fft$category[converging] <- "converging"
                fmod_fft$category[intensifying] <- "intensifying"
                table(fmod_fft$significant)/5
                table(fmod_fft$category)/5
                

                fl <- fmod_fft[converging,]
                fg <- fmod_fft[intensifying,]
                fl <- fmod_fft[-c(intensifying),]
                    
                #Removing outliers
                #fg <- fg[fg$countrycode!="SSD",]
                #fg <- fg[fg$countrycode!="SLE",]
                #fg <- fg[fg$countrycode!="GNQ",]
                #fpc <- fpc[fpc$countrycode!="LBY",]
                #fg <- fg[fg$countrycode!="KOR",]
                #fg <- fg[fg$countrycode!=fg$countrycode[which(fg$Estimate==min(fg$Estimate,na.rm=TRUE))],]
                #fl <- fl[fl$countrycode!=fl$countrycode[which(fl$Estimate==min(fl$Estimate,na.rm=TRUE))],]
                #fl <- fl[fl$countrycode!=fl$countrycode[which(fl$Estimate==min(fl$Estimate,na.rm=TRUE))],]
                
                # Figure 3 (start)
                    plot_fg <- ggplot(data=fg,aes(x=filters,y=Estimate*100, group = countrycode,color=factor(lowsignificant)))+
                    geom_line()+
                    scale_colour_discrete(name="Low-freq estimate significantly different from zero") +
                    theme_bw() + xlab("Minimum Periodicity after Filtering")+
                    geom_hline(yintercept=0,lty=2)+
                    #geom_dl(data=fg[fg$significant==1,],aes(label = countrycode), method = list(dl.combine("last.points")), cex = 0.9)+
                    ylab("Estimated Effect of \n 1 Degree Warming on Growth (pp)") +
                    theme(axis.line = element_line(colour = "black"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank(),legend.position="bottom") +
                    ggtitle("Evidence of Intensifying Effects \n (Low-freq estimate significantly larger than unfiltered estimate)")



                    
                    plot_fl <- ggplot(data=fl[fl$unfilteredsignificant==1,],aes(x=filters,y=Estimate*100, group = countrycode,color=factor(lowsignificant)))+
                    geom_line()+
                    scale_colour_discrete(name="Loefreq estimate significantly different from zero") +
                    theme_bw()+
                    geom_hline(yintercept=0,lty=2)+
                    #geom_dl(data=fpc[fpc$significant==1,],aes(label = countrycode), method = list(dl.combine("last.points")), cex = 0.9)+
                    xlab("Minimum Periodicity after Filtering")+ ylab("")+
                    theme(axis.line = element_line(colour = "black"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank()) +
                    ggtitle("Evidence of Converging towards zero \n (Low-freq estimate significantly smaller than unfiltered estimate)")



                    d <- ggarrange(map_growth,ggarrange(plot_fg,plot_fl,ncol=2,nrow=1,common.legend=FALSE,legend="bottom"),ncol=1,nrow=2) 
                    d
                    table(fmod_fft$category)/5
                    table(fmod_fft$significant)/5

                    #ggsave('Fig3_cat.png',dpi=500) 
                    
                    # Mapping estimates
                        glimpse(fmod_fft)
                        unique(fmod_fft$countrycode)
                        install.packages('rnaturalearthdata')
                        library('rnaturalearthdata')
                        world <- ne_countries(scale = "medium", returnclass = "sf")
                        
                        
                        fmod_fft$iso_a3 <- fmod_fft$countrycode
                        fmod_fft$test1 <- NA
                        fmod_fft$test1[fmod_fft$category=="growth"] <- 1
                        fmod_fft_map <- merge(world,fmod_fft,by="iso_a3")
                        library(RColorBrewer)
                            cols <- c("1" = "#8dd3c7", "0" = "#ffffb3")
                            map_growth <- ggplot(data = fmod_fft_map) +
                                geom_sf(aes(fill = factor(test1)))+
                                scale_fill_manual(name = "Growth Effects",
                                                    labels = c("Detected","Not detected"),
                                values=cols)+
                                theme(legend.position="bottom")+
                                ggtitle("Location of Growth Effects")

                            
                           
                            growtheff <- fmod_fft
                            cg <- unique(growtheff$countrycode)
                            for (i in 1:length(cg)){
                                if(is.na(growtheff$Estimate[5+5*(i-1)])){
                                    growtheff$Estimate[5+5*(i-1)] <- growtheff$Estimate[4+5*(i-1)]
                                    
                                if(is.na(growtheff$Estimate[4+5*(i-1)])){
                                    growtheff$Estimate[5+5*(i-1)] <- growtheff$Estimate[3+5*(i-1)]
                                    
                                if(is.na(growtheff$Estimate[3+5*(i-1)])){
                                    growtheff$Estimate[5+5*(i-1)] <- growtheff$Estimate[2+5*(i-1)]
                                    }}}
                            }
                            growtheff$Estimate
                            growtheff <- growtheff[which(growtheff$test1==1),]

                            growtheff <- growtheff[which(growtheff$frequencies==15),]
                            growtheff <- merge(world,growtheff,by="iso_a3")

                            ggplot(data = fmod_fft_map) +
                                geom_sf(fill=NA)+
                                geom_sf(data=growtheff,aes(fill = Estimate*100))+
                                scale_fill_gradient2(
                                    name = "Estimated impact \n (% per Degree)",
                                    low = "red",
                                    mid = "white",
                                    high = "#00BFC4",
                                    midpoint = 0,
                                    space = "Lab",
                                    na.value = "grey50",
                                    guide = "colourbar",
                                    aesthetics = "fill")+
                                    ggtitle("Detected growth effects")+
                                theme(legend.position="bottom")
                                
                                
                                
                                scale_fill_manual(name = "Growth Effects",
                                                    labels = c("Detected","Not detected"),
                                values=cols)+
                                theme(legend.position="bottom")+
                                ggtitle("Location of Growth Effects")


                        ggarrange(plot_fg,map_growth,ncol=1,nrows=2)


                            ggplot(data = fmod_fft_map) +
                                geom_sf()+
                                facet_wrap(~frequencies, labeller = labeller(frequencies = 
                                    c("0"="Unfiltered",
                                    "3" = "f3",
                                    "5" = "f5",
                                    "10" = "f10",
                                    "15" = "f15"))) +
                                geom_sf(data = fmod_fft_map, aes(fill = Estimate)) +
                                scale_fill_gradient2(low = "red",
                                    mid = "white",
                                    high = "#00BFC4",
                                    midpoint = 0,
                                    space = "Lab",
                                    na.value = "grey50",
                                    guide = "colourbar",
                                    aesthetics = "fill")
                        
                        ggplot(data = fmod_fft_map) +
                            geom_sf()+
                            geom_sf(data = fmod_fft_map[fmod_fft_map$test1==1,], aes(fill = Estimate)) +
                            scale_fill_gradient2(low = "red",
                                mid = "white",
                                high = "#00BFC4",
                                midpoint = 0,
                                space = "Lab",
                                na.value = "grey50",
                                guide = "colourbar",
                                aesthetics = "fill")

                    # Mapping estimates
                # Figure 3 (start)

                #Plotting some outliers 
                    lby <- ggplot(data=wb[wb$c
                    ountrycode=="LBY",], aes(x=year, y=growth))+geom_line()+theme_bw()+ggtitle("Economic growth in Libya") #ilitary intervention in 2011
                    gnq <- ggplot(data=wb[wb$countrycode=="GNQ",], aes(x=year, y=growth))+geom_line()+theme_bw()+ggtitle("Economic growth in Equatorial Guinea") #In 1995 Mobil discovered oil in Guinea
                    ssd <- ggplot(data=wb[wb$countrycode=="SSD",], aes(x=year, y=growth))+
                    geom_line()+theme_bw()+ggtitle("Economic growth in South Sudan") #In 1995 Mobil discovered oil in Guinea              
                    ukr <- ggplot(data=wb[wb$countrycode=="UKR",], aes(x=year, y=growth))+
                    geom_line()+theme_bw()+ggtitle("Economic growth in Ukraine") #In 1995 Mobil discovered oil in Guinea
                    ggarrange(d,ggarrange(gnq,lby,ncol=2),nrow=2)
                    #ggsave('Fig_pathways_outliers.png',dpi=500) #Figure 3
        # Categorizing significant intensifying (end)
    ## 3.2. categorizing - Plotting Figure 3 (end) 

    ## 3.3. FELM abs estimate by filter - Table 1, columns 1 and 2 (start)
        fmod_fft$continent <- countrycode(sourcevar = fmod_fft$countrycode,
                                origin = "iso3c",
                                destination = "continent")
        fmod_fft$invse <- 1/fmod_fft$StandardError
        fmod_fft <- fmod_fft[!is.na(fmod_fft$invse),]
        felm_1 <- felm(absestimate ~ factor(frequencies)|0|0|continent, data = fmod_fft,weights = fmod_fft$invse)
        felm_2 <- felm(absestimate ~ factor(frequencies)|0|0|continent, data = fmod_fft)
        
        #stargazer(felm_2,felm_1, type="html", out="felm_1_2.html")
        stargazer(felm_2,felm_1,type="text")
    ## 3.3. FELM abs estimate by filter - Table 1, columns 1 and 2 (end)
    
    ## 3.4. Mean effect accross filters and socioeconomic variables; Supp Fig 1 and 2 (start)
        # Adding socioeconomic variables (start)
            wb2 <- wb[!is.na(wb$gdppc),]
            minyear <- aggregate(wb2$year, by = list(wb2$countrycode), FUN = min, na.rm = TRUE )
            names(minyear) <- c("countrycode","year")
            minyear$extra <- paste(minyear$countrycode,minyear$year,sep="")
            wb2 <- wb2[which(wb2$extra %in% minyear$extra),]
            wb2$loggdppc <- log(wb2$gdppc)
            wb2 <- wb2[,which(names(wb2) %in% c("countrycode","loggdppc"))]
            fmod_fft <- merge(fmod_fft, wb2,by = "countrycode")
            meangdppc <- aggregate(wb$gdppc, by = list(wb$countrycode), FUN = mean, na.rm = TRUE )
            names(meangdppc ) <-  c("countrycode", "meangdppc ")
            fmod_fft <- merge(fmod_fft, meangdppc , by = "countrycode")
            pop_data <- wb_data("SP.POP.TOTL", start_date = 2019, end_date = 2019)
            names(pop_data)[2] <- "countrycode"
            fmod_fft <- merge(fmod_fft, pop_data,by = "countrycode") 
            a <- wb_search("GDP.*PPP")
            a <- a[6,1]
            gdp_data <- wb_data(a, start_date = 2019, end_date = 2019)
            glimpse(gdp_data)
            names(gdp_data)[2] <- "countrycode"
            fmod_fft <- merge(fmod_fft, gdp_data,by = "countrycode")         
            a <- wb_search("GDP per capita")
            a <- a[10,1]
            gdp_data <- wb_data(a, start_date = 2019, end_date = 2019)
            names(gdp_data)[2] <- "countrycode"
            fmod_fft <- merge(fmod_fft, gdp_data,by = "countrycode")
            meanT <- aggregate(wb$UDel_pop_temp, by = list(wb$countrycode), FUN = mean, na.rm = TRUE )
            names(meanT ) <-  c("countrycode", "meanT")
            fmod_fft <- merge(fmod_fft, meanT,by = "countrycode")
            fmod_fft <- fmod_fft[, !duplicated(colnames(fmod_fft))]
            #ggplot(data=fmod_fft, aes(x=meanT,y= log(NY.GDP.PCAP.PP.CD), color=category))+geom_point()+theme_bw()
        # Adding socioeconomic variables (end)

        # weighted average (start)

            meanestimate <- data.frame(mean = rep(0,15),filter = rep(0,15), weight = rep(0,15), var = rep(0,15))
            frequency <- c(0,3,5,10,15)
            fmod_fft15 <- fmod_fft[(fmod_fft$frequencies==15),]
            fmod_fft15 <- fmod_fft15[!is.na(fmod_fft15$Estimate),]
            c15 <- levels(factor(fmod_fft15$countrycode))
            fmod_fft15 <- fmod_fft[(fmod_fft$countrycode %in% c15),]
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
                meanestimate$mean[1 + ((i-1)*3)]=mean(fmod_fft15$Estimate[fmod_fft15$frequencies==frequency[i]],na.rm=TRUE)
                meanestimate$weight[1 + ((i-1)*3)]="Unweighted"
                meanestimate$var[1 + ((i-1)*3)] <- var(fmod_fft15$Estimate[fmod_fft15$frequencies==frequency[i]])


                x1 <- fmod_fft15[!is.na(fmod_fft15$SP.POP.TOTL),]
                meanestimate$mean[2 + ((i-1)*3)] <- weighted.mean(x=x1$Estimate[x1$frequencies==frequency[i]], x1$SP.POP.TOTL[x1$frequencies==frequency[i]],na.rm=TRUE)
                meanestimate$weight[2 + ((i-1)*3)] <- "Population in 2019"
                meanestimate$var[2 + ((i-1)*3)] <- weighted.var(x=x1$Estimate[x1$frequencies==frequency[i]], w = x1$SP.POP.TOTL[x1$frequencies==frequency[i]])

                x1 <- fmod_fft15[!is.na(fmod_fft15$NY.GDP.MKTP.PP.KD),]
                meanestimate$mean[3 + ((i-1)*3)] <- weighted.mean(x=x1$Estimate[x1$frequencies==frequency[i]], x1$NY.GDP.MKTP.PP.KD[x1$frequencies==frequency[i]],na.rm=TRUE)
                meanestimate$weight[3 + ((i-1)*3)] <- "GDP, PPP constant 2017 USD"
                meanestimate$var[3 + ((i-1)*3)]  <- weighted.var(x=x1$Estimate[x1$frequencies==frequency[i]], w = x1$NY.GDP.MKTP.PP.KD[x1$frequencies==frequency[i]])

            }
            #See weighted means
            meanestimate
        # weighted average (end)

        # Socioeconomic characteristics across categories (start)
            
            fmod_fft$catsign <- paste(fmod_fft$category,fmod_fft$signlofreq,sep=".")
            fmod_fft <- fmod_fft[fmod_fft$category!="other",]
            fmod_fft$logGDPpc <- log(fmod_fft$NY.GDP.PCAP.PP.KD)
            densityGDP6 <-  ggstatsplot::ggbetweenstats(
                data = fmod_fft,
                x = catsign,
                y = logGDPpc
                )

            densitymeanT6 <-  ggstatsplot::ggbetweenstats(
            data = fmod_fft,
            x = catsign,
            y = meanT
            )

            densityGDP3 <-  ggstatsplot::ggbetweenstats(
            data = fmod_fft,
            x = category,
            y = logGDPpc
            )

            densitymeanT3 <-  ggstatsplot::ggbetweenstats(
            data = fmod_fft,
            x = category,
            y = meanT
            )
        
            cat <- data.frame(T = rep(0,6),minT=rep(0,6),maxT=rep(0,6),G = rep(0,6),minG=rep(0,6),maxG=rep(0,6),categories=c("Constant.negative","Constant.positive","Converging.negative","Converging.positive","Intensifying.negative","Intensifying.positive"))
            cat$T <- aggregate(fmod_fft$meanT, list(fmod_fft$catsign), mean, na.rm=TRUE)[1:6,2]
            cat$G <- aggregate(fmod_fft$logGDPpc, list(fmod_fft$catsign), mean, na.rm=TRUE)[1:6,2]
            cat$minT <- aggregate(fmod_fft$meanT, list(fmod_fft$catsign), min, na.rm=TRUE)[1:6,2]
            cat$maxT <- aggregate(fmod_fft$meanT, list(fmod_fft$catsign),max, na.rm=TRUE)[1:6,2]
            cat$minG <- aggregate(fmod_fft$logGDPpc, list(fmod_fft$catsign), min, na.rm=TRUE)[1:6,2]
            cat$maxG <- aggregate(fmod_fft$logGDPpc, list(fmod_fft$catsign),max, na.rm=TRUE)[1:6,2]
            cat$category <- c("Constant","Constant","Converging","Converging","Intensifying","Intensifying")
            cat$signlofreq <- c("Negative","Positive","Negative","Positive","Negative","Positive")
            
            cat6spread <- ggplot(data=cat, aes(x=T,y=G,xmin=minT,xmax=maxT,ymin=minG,ymax=maxG, shape=category, color = signlofreq,group=categories))+
            geom_point()+geom_errorbar()+geom_errorbarh()+
            scale_color_discrete(name="Direction of the effect")+
            scale_color_brewer(name="Direction of the effect",palette="Dark2")+
            scale_shape_discrete(name="Category")+
            labs(x="Country's population-weighted temperature", y="log(GDP per capita)")+theme_bw()

            cat6datapoints <- ggplot(data=fmod_fft, aes(x=meanT,y=logGDPpc, shape=category, color = signlofreq))+
            geom_point()+
            scale_color_discrete(name="Direction of the effect")+
            scale_color_brewer(name="Direction of the effect",palette="Dark2")+
            scale_shape_discrete(name="Category")+
            labs(x="Country's population-weighted temperature", y="log(GDP per capita)")+theme_bw()

            ggarrange(cat6spread,cat6datapoints,ncol=1,common.legend=TRUE)
            #ggsave("SupFig3_cat6.png",dpi=600)
            cat <- data.frame(T = rep(0,3),minT=rep(0,3),maxT=rep(0,3),G = rep(0,3),minG=rep(0,3),maxG=rep(0,3),category=c("Constant","Converging","Intensifying"))
            cat$T <- aggregate(fmod_fft$meanT, list(fmod_fft$category), mean, na.rm=TRUE)[,2]
            cat$G <- aggregate(fmod_fft$logGDPpc, list(fmod_fft$category), mean, na.rm=TRUE)[,2]
            cat$minT <- aggregate(fmod_fft$meanT, list(fmod_fft$category), min, na.rm=TRUE)[,2]
            cat$maxT <- aggregate(fmod_fft$meanT, list(fmod_fft$category),max, na.rm=TRUE)[,2]
            cat$minG <- aggregate(fmod_fft$logGDPpc, list(fmod_fft$category), min, na.rm=TRUE)[,2]
            cat$maxG <- aggregate(fmod_fft$logGDPpc, list(fmod_fft$category),max, na.rm=TRUE)[,2]
            cat3spread <- ggplot(data=cat, aes(x=T,y=G,xmin=minT,xmax=maxT,ymin=minG,ymax=maxG,color=category))+
            geom_point()+geom_errorbar()+geom_errorbarh()+
            labs(x="Mean Temperature", y="logGDPpc")+theme_bw()

            cat3datapoints <- ggplot(data=fmod_fft, aes(x=meanT,y=logGDPpc, color=category))+
            geom_point()+
            labs(x="Mean Temperature", y="logGDPpc")+theme_bw()

            ggarrange(densityGDP6,densitymeanT6,cat6spread,cat6datapoints,nrow=4)
            #ggsave("soecioeconomic_6cat.png",dpi=600)
            ggarrange(densityGDP3,densitymeanT3,cat3spread,cat3datapoints,nrow=4)
            #ggsave("soecioeconomic_3cat.png",dpi=600)

            a=ggplot(fmod_fft[fmod_fft$frequencies=="15",],
            aes(x=meanT,y=Estimate))
            a = a+geom_point()+
            geom_smooth(method = "lm",
                    color = "#00BFC4", show.legend = FALSE) +  theme_bw()+
                    geom_hline(yintercept=0,lty=2) +
                    labs(title="Temperature effect Filter = 15")
            a

            b=ggplot(fmod_fft[fmod_fft$frequencies=="15",],aes(x=logGDPpc,
            y=Estimate))
            b = b+geom_point()+
            #geom_errorbar()+
            geom_smooth(method = "lm", 
                    color = "#00BFC4", show.legend = FALSE) +  theme_bw()+
                    geom_hline(yintercept=0,lty=2) +
                    labs(title="Temperature effect Filter = 15")
            ggarrange(a,b)
             #ggsave("lowfreq_meanTlogGDP.png",dpi=600)
        # Socioeconomic characteristics across categories (end)
    ## 3.4. Mean effect accross filters and socioeconomic variables; Supp Fig 1 and 2 (end)

    ## 3.5. Other datasets - Table 1, columns 3 to 6; Supp Fig 3 (start)
        # Analizing Barro-Ursua dataset (start)
            # Supp Fig. 3 - bottom panel (start)
                fmod_fft <- fullmods_filter
                fmod_fft <- fmod_fft[fmod_fft$econdata=="barro",]
                fmod_fft <- fmod_fft[fmod_fft$climdata=="UDel",]
                fmod_fft$sign <- rep(0,dim(fmod_fft)[1])
                glimpse(fmod_fft)
                    for (i in 1:217){
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
                    table(fmod_fft$sign)


                    fmod_fft$signlofreq <- rep("negative",dim(fmod_fft)[1])

                    for (i in 1:217){
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
                        m <- fmod_fft$Estimate[lastfreq+5*(i-1)]
                        if (is.na(m)){next}
                        if(m>0 ){
                            fmod_fft$signlofreq[(1+5*(i-1)):(5+5*(i-1))]="positive"
                            }
                    }

                # Categorizing  (start)
                    fmod_fft$high95 <- fmod_fft$Estimate + fmod_fft$StandardError*1.64
                    fmod_fft$low95 <- fmod_fft$Estimate - fmod_fft$StandardError*1.64
                    uncertain <- c(which(fmod_fft$high95>0 & fmod_fft$low95 <0))
                    positive_Constant <- 0
                    positive_Intensifying <- 0
                    positive_Converging <- 0
                    fmod_fft$absestimate <- abs(fmod_fft$Estimate)
                    for (i in 1:217){
                                if(is.null(fmod_fft$absestimate[1+5*(i-1)] )){next}
                                if(!is.na(fmod_fft$absestimate[5+5*(i-1)])){
                                    lastfreq <- 5
                                    }else if(!is.na(fmod_fft$absestimate[4+5*(i-1)])){
                                        lastfreq <- 4
                                    }else if(!is.na(fmod_fft$absestimate[3+5*(i-1)])){
                                        lastfreq <- 3
                                    }else if(!is.na(fmod_fft$absestimate[2+5*(i-1)])){
                                        lastfreq <- 2
                                    }else{next}
                                if(fmod_fft$absestimate[1+5*(i-1)]>0){ #all because absolute value
                                    if(fmod_fft$absestimate[lastfreq+5*(i-1)]>0){ #all because using abs value
                                        if((fmod_fft$absestimate[lastfreq+5*(i-1)] )>(fmod_fft$absestimate[1+5*(i-1)] /2)){ #larger than half the first estimate
                                            if((fmod_fft$absestimate[lastfreq+5*(i-1)] )>(fmod_fft$absestimate[1+5*(i-1)] *1.5)){
                                                positive_Intensifying <- c(positive_Intensifying,(1+5*(i-1)):(5+5*(i-1)))
                                            } else{
                                                if((fmod_fft$Estimate[lastfreq+5*(i-1)]*fmod_fft$Estimate[1+5*(i-1)])<0){
                                                    positive_Converging <- c(positive_Converging,(1+5*(i-1)):(5+5*(i-1)))
                                                } else{
                                                #same sign? if not, is decreasinfg, else:
                                                positive_Constant <- c(positive_Constant,(1+5*(i-1)):(5+5*(i-1)))}} }
                                                else{
                                                    positive_Converging <- c(positive_Converging,(1+5*(i-1)):(5+5*(i-1)))
                                                }}}}
                        fmod_fft$significant <- 1
                        fmod_fft$significant[uncertain]<- 0
                        for (i in 1:217){
                            if(!is.na(fmod_fft$absestimate[5+5*(i-1)])){
                                    lastfreq <- 5
                                    }else if(!is.na(fmod_fft$absestimate[4+5*(i-1)])){
                                        lastfreq <- 4
                                    }else if(!is.na(fmod_fft$absestimate[3+5*(i-1)])){
                                        lastfreq <- 3
                                    }else if(!is.na(fmod_fft$absestimate[2+5*(i-1)])){
                                        lastfreq <- 2
                                    }else{next}
                            if((lastfreq+5*(i-1)) %in% uncertain){ 
                                fmod_fft$significant[1+5*(i-1)] <- 0
                                fmod_fft$significant[2+5*(i-1)] <- 0
                                fmod_fft$significant[3+5*(i-1)] <- 0
                                fmod_fft$significant[4+5*(i-1)] <- 0
                                fmod_fft$significant[5+5*(i-1)] <- 0
                            } else {
                                fmod_fft$significant[1+5*(i-1)] <- 1
                                fmod_fft$significant[2+5*(i-1)] <- 1
                                fmod_fft$significant[3+5*(i-1)] <- 1
                                fmod_fft$significant[4+5*(i-1)] <- 1
                                fmod_fft$significant[5+5*(i-1)] <- 1
                            }
                        }
                        glimpse(fmod_fft)
                        
                        fmod_fft$filters <- rep(c("Unfiltered","3 years","5 years", "10 years", "15 years"),43)
                        fmod_fft$filters <- factor(fmod_fft$filters, levels = c("Unfiltered", "3 years", "5 years", "10 years", "15 years"))
                    
                        fmod_fft$category <-"other"
                        fmod_fft$category[positive_Constant] <- "Constant"
                        fmod_fft$category[positive_Intensifying] <- "Intensifying"
                        fmod_fft$category[positive_Converging] <- "Converging"
                        

                        fpc <- fmod_fft[positive_Constant,]
                        fpi <- fmod_fft[positive_Intensifying,]
                        fpd <- fmod_fft[positive_Converging,]
                        
                        fpi <- fpi[fpi$countrycode!="SSD",]
                        fpi <- fpi[fpi$countrycode!="GNQ",]

                    
                        plot_fpi <- ggplot(data=fpi,aes(x=filters,y=Estimate*100, group = countrycode,color=factor(significant)))+
                        geom_line()+
                        # scale_x_discrete(breaks=c("0","5","10","15"), labels=c("Unfiltered","5 years", "10 years", "15 years"))+  
                        #scale_x_discrete( labels=c("Unfiltered","5 years", "10 years", "15 years")) +
                        scale_colour_discrete(guide = 'none') +
                        theme_bw()+
                        geom_hline(yintercept=0,lty=2)+
                        geom_dl(data=fpi[fpi$significant==1,],aes(label = countrycode), method = list(dl.combine("last.points")), cex = 0.9)+
                        labs(title="Intensifying")  + xlab("") + ylab("Estimated Effect of \n 1 Degree Warming on Growth (pp)") +
                        theme(axis.line = element_line(colour = "black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank()) 




                        plot_fpc <- ggplot(data=fpc,aes(x=filters,y=Estimate*100, group = countrycode,color=factor(significant)))+
                        geom_line()+
                        # scale_x_discrete(breaks=c("0","5","10","15"), labels=c("Unfiltered","5 years", "10 years", "15 years"))+  
                        #scale_x_discrete( labels=c("Unfiltered","5 years", "10 years", "15 years")) +
                        scale_colour_discrete(name="Low-frequency estimate significantly different from zero") +
                        theme_bw()+
                        geom_hline(yintercept=0,lty=2)+
                        geom_dl(data=fpc[fpc$significant==1,],aes(label = countrycode), method = list(dl.combine("last.points")), cex = 0.9)+
                        labs(title="Constant")  + xlab("Minimum Periodicity after Filtering") + ylab("") +
                        theme(axis.line = element_line(colour = "black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank()) 



                        fpd <- fpd[fpd$countrycode!="LBY",]
                        plot_fpd <-ggplot(data=fpd,aes(x=filters,y=Estimate*100, group = countrycode,color=factor(significant)))+
                        geom_line()+
                        # scale_x_discrete(breaks=c("0","5","10","15"), labels=c("Unfiltered","5 years", "10 years", "15 years"))+  
                        #scale_x_discrete( labels=c("Unfiltered","5 years", "10 years", "15 years")) +
                        scale_colour_discrete(guide = 'none') +
                        theme_bw()+
                        geom_hline(yintercept=0,lty=2)+
                        geom_dl(data=fpd[fpd$significant==1,],aes(label = countrycode), method = list(dl.combine("last.points")), cex = 0.9)+
                        labs(title="Converging")  + xlab("") + ylab("")+
                        theme(axis.line = element_line(colour = "black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank()) 

                        ggarrange(plot_fpd,plot_fpc,plot_fpi,ncol=3,common.legend=TRUE)      
                        dbarro <- ggarrange(plot_fpi,plot_fpc,plot_fpd,ncol=3,common.legend=TRUE) 
                        dbarro
                        table(fmod_fft$category)/5
                # Categorizing  (start) 

    
            # Supp Fig. 3 - bottom panel (start)
    
            # FELM abs estimate by filter - Table 1, columns 3 and 4 (start)
                library(countrycode)
                fmod_fft$continent <- countrycode(sourcevar = fmod_fft$countrycode,
                                        origin = "iso3c",
                                        destination = "continent")
                fmod_fft$invse <- 1/fmod_fft$StandardError
                fmod_fft <- fmod_fft[!is.na(fmod_fft$invse),]
                felm_1barro <- felm(absestimate ~ factor(frequencies)|0|0|continent, data = fmod_fft,weights = fmod_fft$invse)
                summary(felm_1barro)
                felm_2barro <- felm(absestimate ~ factor(frequencies)|0|0|continent, data = fmod_fft)
                summary(felm_2barro)
                
                stargazer(felm_2,felm_1,felm_2barro,felm_1barro, type="html", out="felm_1_2_barro.html")
            # FELM abs estimate by filter - Table 1, columns 3 and 4 (end)

        # Analizing Barro-Ursua dataset (end)

        # Analizing Maddison dataset (start)
            # Supp Fig. 3 - top panel (start)
                fmod_fft <- fullmods_filter
                fmod_fft <- fmod_fft[fmod_fft$econdata=="mad",]
                fmod_fft <- fmod_fft[fmod_fft$climdata=="UDel",]
                #fmod_fft <- fmod_fft[!is.na(fmod_fft$Estimate),]
                fmod_fft$sign <- rep(0,dim(fmod_fft)[1])
                glimpse(fmod_fft)
                    for (i in 1:217){
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
                    table(fmod_fft$sign)


                    fmod_fft$signlofreq <- rep("negative",dim(fmod_fft)[1])

                    for (i in 1:217){
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
                            m <- fmod_fft$Estimate[lastfreq+5*(i-1)]
                            if (is.na(m)){next}
                            if(m>0 ){
                                fmod_fft$signlofreq[(1+5*(i-1)):(5+5*(i-1))]="positive"
                                }
                    }
                    table(fmod_fft$signlofreq)
                #sign of low freq (start)

                #categorizing (start)
                    fmod_fft$high95 <- fmod_fft$Estimate + fmod_fft$StandardError*1.64
                    fmod_fft$low95 <- fmod_fft$Estimate - fmod_fft$StandardError*1.64
                    uncertain <- c(which(fmod_fft$high95>0 & fmod_fft$low95 <0))
                    positive_Constant <- 0
                    positive_Intensifying <- 0
                    positive_Converging <- 0
                    fmod_fft$absestimate <- abs(fmod_fft$Estimate)
                    for (i in 1:217){
                                if(is.null(fmod_fft$absestimate[1+5*(i-1)] )){next}
                                if(!is.na(fmod_fft$absestimate[5+5*(i-1)])){
                                    lastfreq <- 5
                                    }else if(!is.na(fmod_fft$absestimate[4+5*(i-1)])){
                                        lastfreq <- 4
                                    }else if(!is.na(fmod_fft$absestimate[3+5*(i-1)])){
                                        lastfreq <- 3
                                    }else if(!is.na(fmod_fft$absestimate[2+5*(i-1)])){
                                        lastfreq <- 2
                                    }else{next}
                                if(fmod_fft$absestimate[1+5*(i-1)]>0){ #all because absolute value
                                    if(fmod_fft$absestimate[lastfreq+5*(i-1)]>0){ #all because using abs value
                                        if((fmod_fft$absestimate[lastfreq+5*(i-1)] )>(fmod_fft$absestimate[1+5*(i-1)] /2)){ #larger than half the first estimate
                                            if((fmod_fft$absestimate[lastfreq+5*(i-1)] )>(fmod_fft$absestimate[1+5*(i-1)] *1.5)){
                                                positive_Intensifying <- c(positive_Intensifying,(1+5*(i-1)):(5+5*(i-1)))
                                            } else{
                                                if((fmod_fft$Estimate[lastfreq+5*(i-1)]*fmod_fft$Estimate[1+5*(i-1)])<0){
                                                    positive_Converging <- c(positive_Converging,(1+5*(i-1)):(5+5*(i-1)))
                                                } else{
                                                #same sign? if not, is decreasinfg, else:
                                                positive_Constant <- c(positive_Constant,(1+5*(i-1)):(5+5*(i-1)))}} }
                                                else{
                                                    positive_Converging <- c(positive_Converging,(1+5*(i-1)):(5+5*(i-1)))
                                                }}}}
                        fmod_fft$significant <- 1
                        fmod_fft$significant[uncertain]<- 0
                        for (i in 1:217){
                            if(!is.na(fmod_fft$absestimate[5+5*(i-1)])){
                                    lastfreq <- 5
                                    }else if(!is.na(fmod_fft$absestimate[4+5*(i-1)])){
                                        lastfreq <- 4
                                    }else if(!is.na(fmod_fft$absestimate[3+5*(i-1)])){
                                        lastfreq <- 3
                                    }else if(!is.na(fmod_fft$absestimate[2+5*(i-1)])){
                                        lastfreq <- 2
                                    }else{next}
                            if((lastfreq+5*(i-1)) %in% uncertain){ 
                                fmod_fft$significant[1+5*(i-1)] <- 0
                                fmod_fft$significant[2+5*(i-1)] <- 0
                                fmod_fft$significant[3+5*(i-1)] <- 0
                                fmod_fft$significant[4+5*(i-1)] <- 0
                                fmod_fft$significant[5+5*(i-1)] <- 0
                            } else {
                                fmod_fft$significant[1+5*(i-1)] <- 1
                                fmod_fft$significant[2+5*(i-1)] <- 1
                                fmod_fft$significant[3+5*(i-1)] <- 1
                                fmod_fft$significant[4+5*(i-1)] <- 1
                                fmod_fft$significant[5+5*(i-1)] <- 1
                            }
                        }
                        fmod_fft$filters <- rep(c("Unfiltered","3 years","5 years", "10 years", "15 years"),170)
                        fmod_fft$filters <- factor(fmod_fft$filters, levels = c("Unfiltered", "3 years", "5 years", "10 years", "15 years"))
                    
                        fmod_fft$category <-"other"
                        fmod_fft$category[positive_Constant] <- "Constant"
                        fmod_fft$category[positive_Intensifying] <- "Intensifying"
                        fmod_fft$category[positive_Converging] <- "Converging"
                        

                        fpc <- fmod_fft[positive_Constant,]
                        fpi <- fmod_fft[positive_Intensifying,]
                        fpd <- fmod_fft[positive_Converging,]
                        
                        fpi <- fpi[fpi$countrycode!="SSD",]
                        fpi <- fpi[fpi$countrycode!="GNQ",]

                        plot_fpi <- ggplot(data=fpi,aes(x=filters,y=Estimate*100, group = countrycode,color=factor(significant)))+
                        geom_line()+
                        # scale_x_discrete(breaks=c("0","5","10","15"), labels=c("Unfiltered","5 years", "10 years", "15 years"))+  
                        #scale_x_discrete( labels=c("Unfiltered","5 years", "10 years", "15 years")) +
                        scale_colour_discrete(guide = 'none') +
                        theme_bw()+
                        geom_hline(yintercept=0,lty=2)+
                        geom_dl(data=fpi[fpi$significant==1,],aes(label = countrycode), method = list(dl.combine("last.points")), cex = 0.9)+
                        labs(title="Intensifying")  + xlab("") + ylab("Estimated Effect of \n 1 Degree Warming on Growth (pp)") +
                        theme(axis.line = element_line(colour = "black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank()) 




                        plot_fpc <- ggplot(data=fpc,aes(x=filters,y=Estimate*100, group = countrycode,color=factor(significant)))+
                        geom_line()+
                        # scale_x_discrete(breaks=c("0","5","10","15"), labels=c("Unfiltered","5 years", "10 years", "15 years"))+  
                        #scale_x_discrete( labels=c("Unfiltered","5 years", "10 years", "15 years")) +
                        scale_colour_discrete(name="Low-frequency estimate significantly different from zero") +
                        theme_bw()+
                        geom_hline(yintercept=0,lty=2)+
                        geom_dl(data=fpc[fpc$significant==1,],aes(label = countrycode), method = list(dl.combine("last.points")), cex = 0.9)+
                        labs(title="Constant")  + xlab("Minimum Periodicity after Filtering") + ylab("") +
                        theme(axis.line = element_line(colour = "black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank()) 



                        fpd <- fpd[fpd$countrycode!="LBY",]
                        plot_fpd <-ggplot(data=fpd,aes(x=filters,y=Estimate*100, group = countrycode,color=factor(significant)))+
                        geom_line()+
                        # scale_x_discrete(breaks=c("0","5","10","15"), labels=c("Unfiltered","5 years", "10 years", "15 years"))+  
                        #scale_x_discrete( labels=c("Unfiltered","5 years", "10 years", "15 years")) +
                        scale_colour_discrete(guide = 'none') +
                        theme_bw()+
                        geom_hline(yintercept=0,lty=2)+
                        geom_dl(data=fpd[fpd$significant==1,],aes(label = countrycode), method = list(dl.combine("last.points")), cex = 0.9)+
                        labs(title="Converging")  + xlab("") + ylab("")+
                        theme(axis.line = element_line(colour = "black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank()) 

                        ggarrange(plot_fpd,plot_fpc,plot_fpi,ncol=3,common.legend=TRUE)      
                        dmadd <- ggarrange(plot_fpi,plot_fpc,plot_fpd,ncol=3,common.legend=TRUE) 
                        dmadd
                        table(fmod_fft$category)/5
                        ggarrange(dmadd,dbarro,ncol=1)
                        #ggsave('barro_mad_categories.png',dpi=600)
                                    text <- "Barro Ursua Dataset"

                        # Create a text grob
                        tgrob <- text_grob(text,size = 20)
                        # Draw the text
                        plot_barro <- as_ggplot(tgrob)


                        text <- "Maddison Porject Dataset"

                        # Create a text grob
                        tgrob <- text_grob(text,size = 20)
                        # Draw the text
                        plot_mad <- as_ggplot(tgrob) 

                        ggarrange(plot_barro,dbarro,plot_mad, dmadd,
                                ncol = 1,nrow = 4)
                #categorizing (start) 
            # Supp Fig. 3 - top panel (end)
           

            # FELM abs estimate by filter - Table 1, columns 5 and 6 (start)
                    fmod_fft$continent <- countrycode(sourcevar = fmod_fft$countrycode,
                                            origin = "iso3c",
                                            destination = "continent")
                    fmod_fft$invse <- 1/fmod_fft$StandardError
                    fmod_fft <- fmod_fft[!is.na(fmod_fft$invse),]
                    felm_1mad <- felm(absestimate ~ factor(frequencies)|0|0|continent, data = fmod_fft,weights = fmod_fft$invse)
                    felm_2mad <- felm(absestimate ~ factor(frequencies)|0|0|continent, data = fmod_fft)
                    
                    stargazer(felm_2,felm_1,felm_2barro,felm_1barro, felm_2mad,felm_1mad,type="html", out="felm_wb_barro_mad.html")
            # FELM abs estimate by filter - Table 1, columns 5 and 6 (end)
        
        # Analizing Maddison dataset (end)
    ## 3.5. Other datasets - Table 1, columns 3 to 6; Supp Fig 3 (end)

## 3. Empricial analysis - Figure 3 and Table 1  (end)



#Create Temperature Shocks
    year <- seq(from=1,to=10)
    duration_y <- c("1-year","2-year") #duration of the shock
    filters <- c("Unfiltered",">3",">4",">5") #filtering the shock

    for (i in 1:length(duration_y)){
        tempshock <- rep(0,10)
        for (j in 1:length(filters)){
            tempshock[2:(1+i)] <- 1
            duration <- rep(duration_y[i],10)
            filter <- rep(filters[j],10)
            if(j==1){
                growthshock <- 2 - 0.5 * tempshock
            }
            if (i == 1 & j==1){
                shock <- data.frame(year, tempshock,duration,filter,growthshock)
            }else{
                if(j!=1) tempshock <- pass.filt(tempshock, W=1+j, type="low", method="Butterworth")        
                shock <- rbind(shock,data.frame(year, tempshock,duration,filter,growthshock))
            }  
        }
        
    }




#Create Temperature Shocks

#Plot temperature Shocks and filters
    ggplot(data=shock, aes(x=year, y=tempshock, group=interaction(duration,filter)))+
    theme_bw()+
    geom_line(aes(color=filter,lty=duration))
#Plot temperature shocks and filters

# Simulate GDP growth shock
    lm(growthshock~tempshock,data=shock[shock$filter=="Unfiltered",])$coef[2]
    lm(growthshock~tempshock,data=shock[shock$filter==">3",])$coef[2]
    lm(growthshock~tempshock,data=shock[shock$filter==">4",])$coef[2]
    lm(growthshock~tempshock,data=shock[shock$filter==">5",])$coef[2]
# Simulate GDP growth shock
