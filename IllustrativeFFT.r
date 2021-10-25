#libraries (start)
  install.packages('forecast')
    library(TTR)
    tsyears <- function(ts) as.numeric(trunc(time(ts)))
    library(mFilter)
    library(ncdf4)
    library('TSA')
    library(tidyverse)
    library(ggplot2)
    library('spectral')
    library('forecast')
#libraries (end)

# Define functions for FFT analysis (start)

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

# Define functions for FFT analysis (end)



#read in global temperature record for Fourier decomposition
  file=nc_open("C:/Users/bastien/Box/Long Run GDP Growth/Data/LMR Data/gmt_MCruns_ensemble_full.nc")
  gtemp=ncvar_get(file,"gmt")
  gtemp=apply(gtemp,MARGIN=3,FUN=mean)


  years <- 1:51
  y     <- gtemp[years]
  #y <- rep(0,50)
  #y[20] <- 0.1
  tempdf <- data.frame(y[3:length(y)],y[2:(length(y)-1)],y[1:(length(y)-2)])
  names(tempdf) <- c("temp","templag","years")
  tempdf$g <- -0.5 * tempdf$temp 
  tempdf$l <- -0.5 * tempdf$temp- (-0.5)*tempdf$templag
  #tempdf$gl <- -0.2 * tempdf$temp - (-0.3)*tempdf$templag  - 0.3 * tempdf$temp 
  tempdf$gl <- -0.2 * tempdf$temp - (-0.3)*tempdf$templag  - 0.3 * tempdf$temp 
  tempdf$harmonic <- "original"
  tempdf$coeff_g <- 0
  tempdf$coeff_l <- 0
  tempdf$coeff_gl <- 0

  tempdf$se_g <- 0
  tempdf$se_l <- 0
  tempdf$se_gl <- 0

  tempdf$tempsum <- 0
  tempdf$gsum <- 0
  tempdf$lsum <- 0
  tempdf$glsum <- 0
  tempdf$coeff_gsum <- 0
  tempdf$coeff_lsum <- 0
  tempdf$coeff_glsum <- 0

  X.temp <- fft(tempdf$temp)                   # get amount of each frequency k
  X.l <- fft(tempdf$l)                   # get amount of each frequency k
  X.g <- fft(tempdf$g)                   # get amount of each frequency k
  X.gl <- fft(tempdf$gl)                   # get amount of each frequency k

  par(mfrow=c(3,1))
  plot.frequency.spectrum(X.temp, xlimits=c(0,100))
  plot.frequency.spectrum(X.l, xlimits=c(0,100))
  plot.frequency.spectrum(X.g, xlimits=c(0,100))



  time     <- 49                            # measuring time interval (seconds)
  acq.freq <- 1                          # data acquisition frequency (Hz)
  ts  <- seq(0,time-1/acq.freq,1/acq.freq) # vector of sampling time-points (s) 
  plot.harmonic(X.temp, 1, ts, acq.freq, color="red") 
  years = 1:49
  i <- 1
  for (i in 1:24){
  Xtemp.h <- rep(0,length(X.temp))
  Xtemp.h[(i+1)] <- X.temp[(i+1)] # i-th harmonic
  harmonic.trajectory.temp <- get.trajectory(Xtemp.h, ts, acq.freq=acq.freq)
  Xg.h <- rep(0,length(X.g))
  Xg.h[(i+1)] <- X.g[(i+1)] # i-th harmonic
  harmonic.trajectory.g <- get.trajectory(Xg.h, ts, acq.freq=acq.freq)
  Xl.h <- rep(0,length(X.l))
  Xl.h[(i+1)] <- X.l[(i+1)] # i-th harmonic
  convert.fft(Xl.h)
  harmonic.trajectory.l <- get.trajectory(Xl.h, ts, acq.freq=acq.freq)
  #plot(Re(harmonic.trajectory.l))
  Xgl.h <- rep(0,length(X.gl))
  Xgl.h[(i+1)] <- X.gl[(i+1)] # i-th harmonic
  harmonic.trajectory.gl <- get.trajectory(Xgl.h, ts, acq.freq=acq.freq)
  le <- length(Re(harmonic.trajectory.l))
  #lm_l <- lm(Re(harmonic.trajectory.l)[2:le] ~ Re(harmonic.trajectory.temp)[2:le] + Re(harmonic.trajectory.temp)[1:le-1]  )
  lm_l <- lm(Re(harmonic.trajectory.l) ~ Re(harmonic.trajectory.temp))
  coeff_l <- summary(lm_l)$coefficients[2,1]
  se_l <- summary(lm_l)$coefficients[2,2]
  #lm_g <- lm(Re(harmonic.trajectory.g)[2:le] ~ Re(harmonic.trajectory.temp)[2:le] +  Re(harmonic.trajectory.temp)[1:le-1])
  lm_g <- lm(Re(harmonic.trajectory.g) ~ Re(harmonic.trajectory.temp))
  coeff_g <- summary(lm_g)$coefficients[2,1]
  se_g <-  summary(lm_g)$coefficients[2,2]
  #
  lm_gl <- lm(Re(harmonic.trajectory.gl) ~ Re(harmonic.trajectory.temp))
  coeff_gl <- summary(lm_gl)$coefficients[2,1]
  se_gl <-  summary(lm_gl)$coefficients[2,2]
  #sum of harmonics so far
  tempsum <- tempdf$tempsum[(((i-1)*49)+1):(((i-1)*49)+49)] + Re(harmonic.trajectory.temp)
  gsum <- tempdf$gsum[(((i-1)*49)+1):(((i-1)*49)+49)] + Re(harmonic.trajectory.g)
  lsum <- tempdf$lsum[(((i-1)*49)+1):(((i-1)*49)+49)] + Re(harmonic.trajectory.l)
  glsum <- tempdf$glsum[(((i-1)*49)+1):(((i-1)*49)+49)] + Re(harmonic.trajectory.gl)
  coeff_gsum <- summary(lm(gsum ~ tempsum))$coefficients[2,1]
  coeff_lsum <- summary(lm(lsum ~ tempsum))$coefficients[2,1]
  coeff_glsum <- summary(lm(glsum ~ tempsum))$coefficients[2,1]
  hdf <- data.frame(temp=Re(harmonic.trajectory.temp),templag=Re(harmonic.trajectory.temp),years=years,g=Re(harmonic.trajectory.g),l=Re(harmonic.trajectory.l),gl=Re(harmonic.trajectory.gl),
    harmonic=rep(toString(i),49),coeff_g = rep(coeff_g,49),coeff_l = rep(coeff_l,49),coeff_gl = rep(coeff_gl,49),se_g = rep(se_g,49),
    se_l = rep(se_l,49),se_gl = rep(se_gl,49),tempsum=tempsum, gsum = gsum, lsum = lsum, glsum = glsum,coeff_gsum = rep(coeff_gsum,49),coeff_lsum=rep(coeff_lsum,49),coeff_glsum=rep(coeff_glsum,49))
  tempdf <- rbind(tempdf,hdf)
  }

  ggplot(data = tempdf, aes(x = years, y = tempsum, color = harmonic))+geom_line()

  cbp1 <- c( "#0072B2","#009E73","#E69F00")
  tempdf$harmonic <- factor(tempdf$harmonic)
  tempdf <- tempdf[which(tempdf$harmonic != "original"),]
  #sorted_labels <- paste(sort(as.integer(levels(tempdf$harmonic))))
  #tempdf$harmonic <- factor(tempdf$harmonic, levels = sorted_labels)

  tempdf$harmonic <- as.numeric(as.character(tempdf$harmonic))
  glimpse(tempdf)
  tempdf$period <- (time / tempdf$harmonic)*2 #period in years
  ggplot(data=tempdf, aes(x =years, y = temp, color = harmonic, group = harmonic)) + 
  geom_line()
#Getting harmonics


#plots by period (start)
  tempdf$temp <- tempdf$temp +0.015 
  waves_t <- ggplot(data=tempdf[which(tempdf$harmonic != "original"),], aes(x =years, y = temp, group = period)) + 
  #geom_line()
  geom_line(aes(color=period),position = position_stack(reverse = FALSE)) +
  theme_minimal() + 
  scale_color_gradientn(colours = cbp1, values = c(1.0, 0.3, 0.05, 0)) + 
  labs(title="Sum of Temperature armonics",
          y ="", x = "Years")+
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+ 
      geom_errorbar(aes(x=0,ymin=0.2,ymax=0.4)) + 
      geom_text(aes(x=2,y=0.25,label="0.1C")) + 
      labs(color = "Period (years)")




  tempdf$g <- tempdf$g +0.015 
  waves_g <- ggplot(data=tempdf[which(tempdf$harmonic != "original"),], aes(x =years, y = g, group = period)) + 
  #geom_line()
  geom_line(aes(color=period),position = position_stack(reverse = FALSE))+
  theme_minimal() + 
  scale_color_gradientn(colours = cbp1, values = c(1.0, 0.3, 0.05, 0)) + 
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+ 
  labs(title="Sum of GDP harmonics \n (growth effect)",
          y ="", x = "Years") + 
      labs(color = "Period (years)")



  tempdf$l <- tempdf$l +0.015 
  waves_l <- ggplot(data=tempdf[which(tempdf$period != "original"),], aes(x =years, y = l, group = period)) + 
  #geom_line()
  geom_line(aes(color=period),position = position_stack(reverse = FALSE))+
  theme_minimal() + 
  scale_color_gradientn(colours = cbp1, values = c(1.0, 0.3, 0.05, 0)) + 
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+ 
  labs(title="Sum of GDP harmonics \n (levels effect)",
          x ="Years", y = "") + 
      labs(color = "Period (years)")


      

  tempdf$gl <- tempdf$gl +0.015 
  waves_gl <- ggplot(data=tempdf[which(tempdf$period != "original"),], aes(x =years, y = gl, group = period)) + 
  #geom_line()
  geom_line(aes(color=period),position = position_stack(reverse = FALSE))+
  theme_minimal() + 
  scale_color_gradientn(colours = cbp1, values = c(1.0, 0.3, 0.05, 0)) + 
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+ 
  labs(title="Sum of GDP harmonics \n (combined effect)",
          x ="Years", y = "") + 
      labs(color = "Period (years)")

  library(ggpubr)
  ggarrange(waves_t,ggarrange(waves_l,waves_g),nrow=2)



  waves_coeff_l <- ggplot(data=tempdf[which(tempdf$harmonic != "original"),], aes(x =period, y = coeff_l)) + 
  #geom_line(aes(color=harmonic),position = position_stack(reverse = FALSE))+theme_bw()
  geom_point(aes(color=period)) + theme_bw()+ 
  scale_color_gradientn(colours = cbp1, values = c(1.0, 0.3, 0.05, 0)) + 
  labs(title="Regression on individual harmonics \n (levels effect)",
          y ="Estimated coefficient", x = "Period (years)") +
      geom_errorbar(aes(x=period,ymin=coeff_l-se_l,ymax=coeff_l+se_l, colour = period)) + 
      labs(color = "Period (years)") 
  #geom_point()

  waves_coeff_g <- ggplot(data=tempdf[which(tempdf$harmonic != "original"),], aes(x =period, y = coeff_g)) + 
  #geom_line(aes(color=harmonic),position = position_stack(reverse = FALSE))+theme_bw()
  geom_point(aes(color=period)) + theme_bw()+ 
  scale_color_gradientn(colours = cbp1, values = c(1.0, 0.3, 0.05, 0)) + 
  labs(title="Regression on individual harmonics \n (growth effect)",
          y ="Estimated coefficient", x = "Period (years)")+
      geom_errorbar(aes(x=period,ymin=coeff_g-se_g,ymax=coeff_g+se_g, colour = period)) + 
      labs(color = "Period (years)")
  #geom_point()


  waves_coeff_gl <- ggplot(data=tempdf[which(tempdf$harmonic != "original"),], aes(x =period, y = coeff_gl)) + 
  #geom_line(aes(color=harmonic),position = position_stack(reverse = FALSE))+theme_bw()
  geom_point(aes(color=period)) + theme_bw()+ 
  scale_color_gradientn(colours = cbp1, values = c(1.0, 0.3, 0.05, 0)) + 
  labs(title="Regression on individual harmonics \n (combined effect)",
          y ="Estimated coefficient", x = "Period (years)")+
      geom_errorbar(aes(x=period,ymin=coeff_gl-se_gl,ymax=coeff_gl+se_gl, colour = period)) + 
      labs(color = "Period (years)")
  #geom_point()

  ggarrange(waves_t,ggarrange(waves_l,waves_g,legend ="none"),ggarrange(waves_coeff_l,waves_coeff_g,common.legend =TRUE, legend = "bottom"),nrow=3, legend ="none")

  ggarrange(waves_t,ggarrange(waves_l,waves_g,waves_gl,legend ="none",ncol = 3),ggarrange(waves_coeff_l,waves_coeff_g,waves_coeff_gl,common.legend =TRUE, legend = "bottom", ncol = 3),nrow=3, legend ="none")
#plots by period (end)

#plots by harmonic (start)
  waves_t <- ggplot(data=tempdf[which(tempdf$harmonic != "original"),], aes(x =years, y = temp, group = harmonic)) + 
  #geom_line()
  geom_line(aes(color=harmonic),position = position_stack(reverse = TRUE)) +
  theme_minimal() + 
  scale_color_gradientn(colours = cbp1,breaks = c(3,22),labels = c("low","high"))+
  #, values = c(1.0, 0.3, 0.05, 0)) + 
  labs(title="Sum of Temperature ahrmonics",
          y ="", x = "Years")+
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+ 
      geom_errorbar(aes(x=0,ymin=0.2,ymax=0.4)) + 
      geom_text(aes(x=2,y=0.25,label="0.1C")) + 
      labs(color = "Frequency")




  waves_g <- ggplot(data=tempdf[which(tempdf$harmonic != "original"),], aes(x =years, y = g, group = harmonic)) + 
  #geom_line()
  geom_line(aes(color=harmonic),position = position_stack(reverse = TRUE))+
  theme_minimal() + 
  scale_color_gradientn(colours = cbp1,breaks = c(3,22),labels = c("low","high"))+
    #, values = c(1.0, 0.3, 0.05, 0)) + 
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+ 
  labs(title="Sum of GDP harmonics \n (growth effect)",
          y ="", x = "Years") + 
      labs(color = "Frequency")



  waves_l <- ggplot(data=tempdf[which(tempdf$period != "original"),], aes(x =years, y = l, group = harmonic)) + 
  #geom_line()
  geom_line(aes(color=harmonic),position = position_stack(reverse = TRUE))+
  theme_minimal() + 
  scale_color_gradientn(colours = cbp1,breaks = c(3,22),labels = c("low","high"))+
    #, values = c(1.0, 0.3, 0.05, 0)) + 
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+ 
  labs(title="Sum of GDP harmonics \n (levels effect)",
          x ="Years", y = "") + 
      labs(color = "Frequency")


      

  waves_gl <- ggplot(data=tempdf[which(tempdf$period != "original"),], aes(x =years, y = gl, group = harmonic)) + 
  #geom_line()
  geom_line(aes(color =harmonic),position = position_stack(reverse = TRUE))+
  theme_minimal() + 
  scale_color_gradientn(colours = cbp1,breaks = c(3,22),labels = c("low","high"))+  
    #, values = c(1.0, 0.3, 0.05, 0)) + 
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+ 
  labs(title="Sum of GDP harmonics \n (combined effect)",
          x ="Years", y = "") + 
      labs(color = "Frequency")

  library(ggpubr)
  ggarrange(waves_t,ggarrange(waves_l,waves_g),nrow=2)



  waves_coeff_l <- ggplot(data=tempdf[which(tempdf$harmonic != "original"),], aes(x = harmonic, y = coeff_l)) + 
  #geom_line(aes(color=harmonic),position = position_stack(reverse = TRUE))+theme_bw()
  geom_point(aes(color =harmonic)) + theme_bw()+ 
  scale_color_gradientn(colours = cbp1,breaks = c(3,22),labels = c("low","high"))+
  labs(title="Regression on individual harmonics \n (levels effect)",
          y ="Estimated coefficient", x = "Harmonics") +
      geom_errorbar(aes(x=harmonic,ymin=coeff_l-se_l,ymax=coeff_l+se_l, colour = harmonic)) + 
      labs(color = "Frequency") 
  #geom_point()

  waves_coeff_g <- ggplot(data=tempdf[which(tempdf$harmonic != "original"),], aes(x = harmonic, y = coeff_g)) + 
  #geom_line(aes(color=harmonic),position = position_stack(reverse = TRUE))+theme_bw()
  geom_point(aes(color =harmonic)) + theme_bw()+ 
  scale_color_gradientn(colours = cbp1,breaks = c(3,22),labels = c("low","high"))+
  labs(title="Regression on individual harmonics \n (growth effect)",
          y ="Estimated coefficient", x = "Harmonics")+
      geom_errorbar(aes(x=harmonic,ymin=coeff_g-se_g,ymax=coeff_g+se_g, colour = harmonic)) + 
      labs(color = "Frequency")
  #geom_point()


  waves_coeff_gl <- ggplot(data=tempdf[which(tempdf$harmonic != "original"),], aes(x = harmonic, y = coeff_gl)) + 
  #geom_line(aes(color=harmonic),position = position_stack(reverse = TRUE))+theme_bw()
  geom_point(aes(color =harmonic)) + theme_bw()+ 
  scale_color_gradientn(colours = cbp1,breaks = c(3,22),labels = c("low","high"))+
  labs(title="Regression on individual harmonics \n (combined effect)",
          y ="Estimated coefficient", x = "Harmonics")+
      geom_errorbar(aes(x=harmonic,ymin=coeff_gl-se_gl,ymax=coeff_gl+se_gl, colour = harmonic)) + 
      labs(color = "Frequency")
  #geom_point()

  ggarrange(waves_t,ggarrange(waves_l,waves_g,legend ="none"),ggarrange(waves_coeff_l,waves_coeff_g,common.legend =TRUE, legend = "bottom"),nrow=3, legend ="none")

  ggarrange(waves_t,ggarrange(waves_l,waves_g,waves_gl,legend ="none",ncol = 3),ggarrange(waves_coeff_l,waves_coeff_g,waves_coeff_gl,common.legend =TRUE, legend = "bottom", ncol = 3),nrow=3, legend ="none")
#plots by harmonic (end)

glimpse(tempdf)
estimatedf <- aggregate(tempdf, 
          by = list(unique.values = tempdf$harmonic), 
                               FUN = mean)

par(mfrow=c(3,1))
plot(
 spec.fft(estimatedf$coeff_l,estimatedf$period),
ylab = "Amplitude",
xlab = "Frequency",
#type = "l",
#xlim = c(-1,1),
main = "Spectrum"
)
plot(
 spec.fft(estimatedf$coeff_g,estimatedf$period),
ylab = "Amplitude",
xlab = "Frequency",
#type = "l",
#xlim = c(-1,1),
main = "Spectrum"
)
plot(
 spec.fft(estimatedf$coeff_gl,estimatedf$period),
ylab = "Amplitude",
xlab = "Frequency",
#type = "l",
#xlim = c(-1,1),
main = "Spectrum"
)
plot(estimatedf$coeff_l)
fftl <- fft(estimatedf$coeff_l)
plot.show(estimatedf$coeff_l)
plot.show(estimatedf$coeff_l, time=1, harmonics=0:7, plot.freq=TRUE)
plot.frequency.spectrum(fftl, xlimits=c(0,12))
plot(fftl)
fftg <- fft(estimatedf$coeff_g)
plot.frequency.spectrum(fftg, xlimits=c(0,50))
fftgl <- fft(estimatedf$coeff_gl)
plot.frequency.spectrum(fftgl, xlimits=c(0,50))
cfl <- convert.fft(fftl)
cfl[1:20,]
cfg <- convert.fft(fftg)
cfg[1:20,]
cfgl <-convert.fft(fftgl)
cfgl[1:20,]

time     <- length(estimatedf$coeff_l)                    # measuring time interval (seconds)
#time     <- 1                # measuring time interval (seconds)
acq.freq <- 1                        # data acquisition frequency (Hz)
ts  <- seq(0,time-(1/acq.freq),1/acq.freq) # vector of sampling time-points (s) 
X.k <- fftl
x.n <- get.trajectory(X.k,ts,acq.freq) 
plot(ts,x.n,type="l",ylim=c(-2,1),lwd=2)
plot(ts,x.n,type="l",lwd=2)
abline(v=0:time,h=-2:4,lty=3); abline(h=0)
plot.harmonic(X.k,1,ts,acq.freq,"red")
plot.harmonic(X.k,2,ts,acq.freq,"green")
plot.harmonic(X.k,3,ts,acq.freq,"blue")
plot.harmonic(X.k,15,ts,acq.freq,"blue")
plot.harmonic(X.k,22,ts,acq.freq,"green")
plot.harmonic(X.k,42,ts,acq.freq,"red")

    plot.datal  <- data.frame(0:(length(fftl)-1), Mod(fftl) / length(fftl))
    plot.datal[2:length(fftl),2] <- 2*plot.datal[2:length(fftl),2]
    freqspecl <- ggplot(data = plot.datal)+
    geom_segment(aes(x=plot.datal[,1],xend=plot.datal[,1], y = 0, yend = plot.datal[,2])) + 
    xlim(0,12) + theme_bw() +
    labs(title="Frequency spectrum of \n impacts on levels",
        y ="Strength", x = "Frequency (cycles per year)")

    plot.datag  <- data.frame(0:(length(fftg)-1), Mod(fftg) / length(fftg))
    plot.datag[2:length(fftg),2] <- 2*plot.datag[2:length(fftg),2]
    freqspecg <- ggplot(data = plot.datag)+
    geom_segment(aes(x=plot.datag[,1],xend=plot.datag[,1], y = 0, yend = plot.datag[,2])) + 
    xlim(0,12) + theme_bw() +
    labs(title="Frequency spectrum of \n  impacts on growth",
        y ="Strength", x = "Frequency (cycles per year)")

    plot.datagl  <- data.frame(0:(length(fftgl)-1), Mod(fftgl) / length(fftgl))
    plot.datagl[2:length(fftgl),2] <- 2*plot.datagl[2:length(fftgl),2]
    freqspecgl <- ggplot(data = plot.datagl)+
    geom_segment(aes(x=plot.datagl[,1],xend=plot.datagl[,1], y = 0, yend = plot.datagl[,2])) + 
    xlim(0,12) + theme_bw() +
    #xlim(0,10) + theme_bw() +
    labs(title="Frequency spectrum of \n impacts on growth and levels",
        y ="Strength", x = "Frequency (cycles per year)")
   

   ggarrange(ggarrange(waves_coeff_l,waves_coeff_g,waves_coeff_gl,common.legend =TRUE, legend = "bottom", ncol = 3),
   ggarrange(freqspecl,freqspecg,freqspecgl,common.legend =TRUE, legend = "bottom", ncol = 3),nrow=2, legend ="none")

time     <- length(estimatedf$coeff_l)                    # measuring time interval (seconds)
acq.freq <- 1                          # data acquisition frequency (Hz)
ts  <- seq(0,time-1/acq.freq,1/acq.freq) # vector of sampling time-points (s) 

coef_harmonics = data.frame(coef_l_harmonic = double(), coef_g_harmonic = double(), 
  coef_gl_harmonic=double(), harmonic = integer(), previous_harmonic = integer())
for (i in 1:((length(fftl)/2))){
  Xl.h <- rep(0,length(fftl))
  Xl.h[(i)] <- fftl[(i)] # i-th harmonic  Mod(fftl[i]) / length(fftgl)
  Xg.h <- rep(0,length(fftg))
  Xg.h[(i)] <- fftg[(i)] # i-th harmonic
  Xgl.h <- rep(0,length(fftgl))
  Xgl.h[(i)] <- fftgl[(i)] # i-th harmonic

  harmonic.trajectory.l <- get.trajectory(Xl.h, ts, acq.freq=acq.freq)
  harmonic.trajectory.g <- get.trajectory(Xg.h, ts, acq.freq=acq.freq)
  harmonic.trajectory.gl <- get.trajectory(Xgl.h, ts, acq.freq=acq.freq)
  

  df <- cbind(coef_l_harmonic = Re(harmonic.trajectory.l),coef_g_harmonic = Re(harmonic.trajectory.g),
    coef_gl_harmonic = Re(harmonic.trajectory.gl),  harmonic = i, 
    previous_harmonic = (seq(0:length( harmonic.trajectory.gl))-1))
  coef_harmonics <- rbind(coef_harmonics,df)
}

glimpse(coef_harmonics)
#coef_harmonics$coef_l_harmonic <- coef_harmonics$coef_l_harmonic+0.15 
waves_l_coef <- ggplot(data=coef_harmonics, aes(x =previous_harmonic, y = coef_l_harmonic, group = harmonic)) + 
#geom_line()
#geom_line(aes(color=harmonic),position = position_stack(reverse = FALSE))+
geom_line(aes(color=harmonic))+
theme_minimal() + 
  scale_color_gradientn(colours = cbp1,breaks = c(3,22),labels = c("low","high"))+
 
labs(title="Sum of GDP harmonics \n (levels effect)",
        x ="harmonic", y = "") + 
    labs(color = "Frequency")

waves_g_coef <- ggplot(data=coef_harmonics, aes(x =previous_harmonic, y = coef_g_harmonic, group = harmonic)) + 
#geom_line()
#geom_line(aes(color=harmonic),position = position_stack(reverse = FALSE))+
geom_line(aes(color=harmonic))+
theme_minimal() + 
  scale_color_gradientn(colours = cbp1,breaks = c(3,22),labels = c("low","high"))+
labs(title="Sum of GDP harmonics \n (growth effect)",
        x ="harmonic", y = "") + 
    labs(color = "Frequency")

waves_gl_coef <- ggplot(data=coef_harmonics, aes(x =previous_harmonic, y = coef_gl_harmonic, group = harmonic)) + 
#geom_line()
#geom_line(aes(color=harmonic),position = position_stack(reverse = FALSE))+
geom_line(aes(color=harmonic))+
theme_minimal() + 
  scale_color_gradientn(colours = cbp1,breaks = c(3,22),labels = c("low","high"))+
 
labs(title="Sum of GDP harmonics \n (combined effect)",
        x ="harmonic", y = "") + 
    labs(color = "Frequency")

ggarrange(waves_l_coef,waves_g_coef,waves_gl_coef, nrow = 1)


  ggarrange(ggarrange(waves_coeff_l,waves_coeff_g,waves_coeff_gl,common.legend =TRUE, legend = "bottom", ncol = 3),
   ggarrange(freqspecl,freqspecg,freqspecgl,common.legend =TRUE, legend = "bottom", ncol = 3),
   ggarrange(waves_l_coef,waves_g_coef,waves_gl_coef,common.legend =TRUE, legend = "bottom", ncol = 3),nrow=3, legend ="none")






waves_coeff_l <- ggplot(data=tempdf[which(tempdf$harmonic != "original"),], aes(x =harmonic, y = coeff_l)) + 
#geom_line(aes(color=harmonic),position = position_stack(reverse = FALSE))+theme_bw()
geom_point(aes(color=harmonic)) + theme_bw()+ 
scale_color_gradientn(colours = cbp1, values = c(1.0, 0.3, 0.05, 0)) + 
labs(title="Regression on individual harmonics \n (levels effect)",
        y ="Estimated coefficient", x = "harmonic (cycles per 50 years)") +
    geom_errorbar(aes(x=harmonic,ymin=coeff_l-se_l,ymax=coeff_l+se_l, colour = harmonic)) + 
    labs(color = "harmonic (cycles per 50 years)") 
#geom_point()

waves_coeff_g <- ggplot(data=tempdf[which(tempdf$harmonic != "original"),], aes(x =harmonic, y = coeff_g)) + 
#geom_line(aes(color=harmonic),position = position_stack(reverse = FALSE))+theme_bw()
geom_point(aes(color=harmonic)) + theme_bw()+ 
scale_color_gradientn(colours = cbp1, values = c(1.0, 0.3, 0.05, 0)) + 
labs(title="Regression on individual harmonics \n (growth effect)",
        y ="Estimated coefficient", x = "harmonic (cycles per 50 years)")+
    geom_errorbar(aes(x=harmonic,ymin=coeff_g-se_g,ymax=coeff_g+se_g, colour = harmonic)) + 
    labs(color = "harmonic (cycles per 50 years)")
#geom_point()


waves_coeff_gl <- ggplot(data=tempdf[which(tempdf$harmonic != "original"),], aes(x =harmonic, y = coeff_gl)) + 
#geom_line(aes(color=harmonic),position = position_stack(reverse = FALSE))+theme_bw()
geom_point(aes(color=harmonic)) + theme_bw()+ 
scale_color_gradientn(colours = cbp1, values = c(1.0, 0.3, 0.05, 0)) + 
labs(title="Regression on individual harmonics \n (combined effect)",
        y ="Estimated coefficient", x = "harmonic (cycles per 50 years)")+
    geom_errorbar(aes(x=harmonic,ymin=coeff_gl-se_gl,ymax=coeff_gl+se_gl, colour = harmonic)) + 
    labs(color = "harmonic (cycles per 50 years)")
#geom_point()

  ggarrange(ggarrange(waves_coeff_l,waves_coeff_g,waves_coeff_gl,common.legend =TRUE, legend = "bottom", ncol = 3),
   ggarrange(freqspecl,freqspecg,freqspecgl,common.legend =TRUE, legend = "bottom", ncol = 3),
   ggarrange(waves_l_coef,waves_g_coef,waves_gl_coef,common.legend =TRUE, legend = "bottom", ncol = 3),nrow=3, legend ="none")




















#########################################################
###########################
################################












waves_coeffsum_l <- ggplot(data=tempdf[which(tempdf$harmonic != "original"),], aes(x =period, y = coeff_lsum)) + 
#geom_line(aes(color=harmonic),position = position_stack(reverse = FALSE))+theme_bw()
geom_point(aes(color=period)) + theme_bw()+ 
scale_color_gradientn(colours = cbp1, values = c(1.0, 0.3, 0.05, 0)) + 
labs(title="Regression on individual harmonics \n (levels effect)",
        y ="Estimated coefficient", x = "Period (years)") +
    #geom_errorbar(aes(x=period,ymin=coeff_l-se_l,ymax=coeff_l+se_l, colour = period)) + 
    labs(color = "Period (years)") 
#geom_point()

waves_coeffsum_g <- ggplot(data=tempdf[which(tempdf$harmonic != "original"),], aes(x =period, y = coeff_gsum)) + 
#geom_line(aes(color=harmonic),position = position_stack(reverse = FALSE))+theme_bw()
geom_point(aes(color=period)) + theme_bw()+ 
scale_color_gradientn(colours = cbp1, values = c(1.0, 0.3, 0.05, 0)) + 
labs(title="Regression on individual harmonics \n (growth effect)",
        y ="Estimated coefficient", x = "Period (years)")+
    #geom_errorbar(aes(x=period,ymin=coeff_g-se_g,ymax=coeff_g+se_g, colour = period)) + 
    labs(color = "Period (years)")
#geom_point()

ggarrange(waves_t,ggarrange(waves_l,waves_g,legend ="none"),ggarrange(waves_coeffsum_l,waves_coeffsum_g,common.legend =TRUE, legend = "bottom"),nrow=3, legend ="none")





glimpse(tempdf)
`%notin%` <- Negate(`%in%`)

tempdf_3filter <- tempdf[which(tempdf$harmonic <10),]
glimpse(tempdf_3filter)
tempdf_3filtersum <- aggregate(tempdf_3filter, 
                               by = list(unique.values = tempdf_3filter$years), 
                               FUN = mean)
                               
ggplot(data = tempdf_3filtersum, aes(x = years, y = temp))+ geom_line()


glimpse(tempdf_3filtersum)
a <- lm(tempdf_3filtersum$g ~ tempdf_3filtersum$temp)
b <- lm(tempdf_3filtersum$l ~ tempdf_3filtersum$temp)
summary(a)
summary(b)



plot.frequency.spectrum(X.temp, xlimits=c(0,50))


table(tempdf$harmonic)






#####END

lml <- lm(tempdf$l ~ tempdf$temp)
lmg <- lm(tempdf$g ~ tempdf$temp)
summary(lmg)
summary(lml)

lml <- lm(l ~ temp, data = tempdf[which(tempdf$harmonic>40),])
lmg <- lm(g ~ temp, data = tempdf[which(tempdf$harmonic>40),])
summary(lmg)
summary(lml)


lml <- lm(l ~ temp, data = tempdf[which(tempdf$harmonic<20),])
lmg <- lm(g ~ temp, data = tempdf[which(tempdf$harmonic<20),])
summary(lmg)
summary(lml)
2*pi*49/50

glimpse(tempdf)

plot.harmonic(X.k,2,ts,acq.freq,"green")
plot.harmonic(X.k,3,ts,acq.freq,"blue")
plot.harmonic(X.k,14,ts,acq.freq,"black")




x.n <- get.trajectory(X.k,ts,acq.freq)   # create time wave

plot(tempdf$temp,type="l")
plot(ts,x.n,type="l",ylim=c(-0.1,0.1),lwd=2)
abline(v=0:time,h=(seq(1:20)-10)*0.01,lty=3); abline(h=0)

plot.harmonic(X.k,1,ts,acq.freq,"red")

X.k <- fft(tempdf$g)                   # get amount of each frequency k

time     <- 50                            # measuring time interval (seconds)
acq.freq <- 1                          # data acquisition frequency (Hz)
ts  <- seq(0,time-1/acq.freq,1/acq.freq) # vector of sampling time-points (s) 
x.n <- get.trajectory(X.k,ts,acq.freq)   # create time wave

plot.harmonic(X.k,1,ts,acq.freq,"blue")



plot.harmonic(X.k,2,ts,acq.freq,"green")
plot.harmonic(X.k,3,ts,acq.freq,"blue")
plot.harmonic(X.k,14,ts,acq.freq,"black")
50/15

plot(Mod(X.k),type="l")
plot(Arg(X.k), type="l")
i <- 1
j <- 2
Xk.h <- rep(0,length(X.k))
Xk.h[(i+1):(i+j)] <- X.k[(i+1):(i+j)] # i-th harmonic
harmonic.trajectory <- get.trajectory(Xk.h, ts, acq.freq=acq.freq)

points(ts, Re(harmonic.trajectory), type="l", col='blue')



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


x <- 1:4
class(x)
fx <- fft(x)
class(fx)
plot(fx)

fft(fft(x), inverse = TRUE)/length(x)



X.k <- fft(c(4,0,0,0))                   # get amount of each frequency k

time     <- 4                            # measuring time interval (seconds)
acq.freq <- 1                          # data acquisition frequency (Hz)
ts  <- seq(0,time-1/acq.freq,1/acq.freq) # vector of sampling time-points (s) 
x.n <- get.trajectory(X.k,ts,acq.freq)   # create time wave

plot(ts,x.n,type="l",ylim=c(-2,4),lwd=2)
abline(v=0:time,h=-2:4,lty=3); abline(h=0)

plot.harmonic(X.k,1,ts,acq.freq,"red")
plot.harmonic(X.k,2,ts,acq.freq,"green")
plot.harmonic(X.k,3,ts,acq.freq,"blue")


y     <- gtemp[1:50]
tempdf <- data.frame(y[2:50],y[1:49])
names(tempdf) <- c("temp","templag")
tempdf$g <- -0.5 * tempdf$temp 
tempdf$l <- -0.5 * tempdf$temp- (-0.5)*tempdf$templag

X.k <- fft(tempdf$l)                   # get amount of each frequency k

time     <- 50                            # measuring time interval (seconds)
acq.freq <- 1                          # data acquisition frequency (Hz)
ts  <- seq(0,time-1/acq.freq,1/acq.freq) # vector of sampling time-points (s) 
x.n <- get.trajectory(X.k,ts,acq.freq)   # create time wave

plot(ts,x.n,type="l",ylim=c(-0.15,0.15),lwd=2)
abline(v=0:time,h=-2:4,lty=3); abline(h=0)

plot.harmonic(X.k,1,ts,acq.freq,"red")

X.k <- fft(tempdf$g)                   # get amount of each frequency k

time     <- 50                            # measuring time interval (seconds)
acq.freq <- 1                          # data acquisition frequency (Hz)
ts  <- seq(0,time-1/acq.freq,1/acq.freq) # vector of sampling time-points (s) 
x.n <- get.trajectory(X.k,ts,acq.freq)   # create time wave

plot.harmonic(X.k,1,ts,acq.freq,"blue")



plot.harmonic(X.k,2,ts,acq.freq,"green")
plot.harmonic(X.k,3,ts,acq.freq,"blue")
plot.harmonic(X.k,14,ts,acq.freq,"black")
50/15

plot(Mod(X.k),type="l")
plot(Arg(X.k), type="l")
i <- 1
j <- 2
Xk.h <- rep(0,length(X.k))
Xk.h[(i+1):(i+j)] <- X.k[(i+1):(i+j)] # i-th harmonic
harmonic.trajectory <- get.trajectory(Xk.h, ts, acq.freq=acq.freq)

points(ts, Re(harmonic.trajectory), type="l", col='blue')
t1 <- data.frame(cbind(Re(harmonic.trajectory)[2:length(harmonic.trajectory)],Re(harmonic.trajectory)[1:(length(harmonic.trajectory)-1)]))
names(t1) <- c("temp","templag")
t1$g <- -0.5 * t1$temp - (-0.5)*t1$templag
t1$l <- -0.5 * t1$temp
plot(t1$temp,type="l")

points(t1$l, type="l", col='blue')

points(t1$g, type="l", col='red')

y     <- gtemp[1:50]
t     <- 1:50
rg    <- diff(range(y))

nff = function(x = NULL, n = NULL, up = 10L, plot = TRUE, add = FALSE, main = NULL, ...){
  #The direct transformation
  #The first frequency is DC, the rest are duplicated
  dff = fft(x)
  #The time
  t = seq(from = 1, to = length(x))
  #Upsampled time
  nt = seq(from = 1, to = length(x)+1-1/up, by = 1/up)
  #New spectrum
  ndff = array(data = 0, dim = c(length(nt), 1L))
  ndff[1] = dff[1] #Always, it's the DC component
  if(n != 0){
    ndff[2:(n+1)] = dff[2:(n+1)] #The positive frequencies always come first
    #The negative ones are trickier
    ndff[length(ndff):(length(ndff) - n + 1)] = dff[length(x):(length(x) - n + 1)]
  }
  #The inverses
  indff = fft(ndff/73, inverse = TRUE)
  idff = fft(dff/73, inverse = TRUE)
  if(plot){
    if(!add){
      plot(x = t, y = x, pch = 16L, xlab = "Time", ylab = "Measurement",
        main = ifelse(is.null(main), paste(n, "harmonics"), main))
      lines(y = Mod(idff), x = t, col = adjustcolor(1L, alpha = 0.5))
    }
    lines(y = Mod(indff), x = nt, ...)
  }
  ret = data.frame(time = nt, y = Mod(indff))
  return(ret)
}

res = nff(x = y, n = 18L, up = 100L, col = 2L)

sum5to18 = nff(x = y, n = 18L, up = 10L, plot = FALSE)
sum5to18$y = sum5to18$y - nff(x = y, n = 4L, up = 10L, plot = FALSE)$y
png("sum5to18.png")
plot(sum5to18, pch = 16L, xlab = "Time", ylab = "Measurement", main = "5th to 18th harmonics sum", type = "l", col = 2)
dev.off()

colors = rainbow(36L, alpha = 0.3)
#nff(x = y, n = 36L, up = 100L, col = colors[1])
png("all_waves.png")
for(i in 1:18){
  ad = ifelse(i == 1, FALSE, TRUE)
  nff(x = y, n = i, up = 100L, col = colors[i], add = ad, main = "All waves up to 18th harmonic")
}


sep = array(data = NA_real_, dim = c(7300L, 2 + 18), dimnames = list(NULL, c("t", paste0("H", 0:18))))
sep[,1:2] = as.matrix(nff(x = y, n = 0, up = 100L, plot = FALSE))

for(i in 1:18L){
  sep[,i+2] = nff(x = y, n = i, up = 100L, plot = FALSE)$y - nff(x = y, n = i-1, up = 100L, plot = FALSE)$y
} 

dev.off()
class(dd)
glimpse(dd)

#save(randomts,gtemp,file="C:/Users/fmoore/Box/Davis Stuff/Social Dynamics Paper/data/naturalvariability.Rdat")

  time=350
  basegr=0.01 #1% per year baseline growth rate
  start=100

  baseline=rep(basegr,time)

  coef=-0.05 #effect size - change in growth per degree warming()
  #note currently this effect size is ~10 times bigger than estimated becasue temperature variation of global mean temp is ~10 times smaller than at national level

  agg=c(1,5,10,25,50)
  p_l <- c(2,6,10,14,18)
  p_u <- c(5,9,13,17,33)
  p_l <- c(0,2,2,2,2)
  p_u <- c(0,3,4,5,6)

  growthsd=0.005 #standard deviation of growth variability unexplained by temperature

  nsims=1000
  sims=array(dim=c(nsims,length(agg),2))
  simssma=array(dim=c(nsims,length(agg),2))
  simsbktrend=array(dim=c(nsims,length(agg),2))
  simsbkcycle=array(dim=c(nsims,length(agg),2))
  for(i in 1:nsims){
    randomtemp=Re(randomts(gtemp))[1:time]
    randomgrowth_g=basegr #growth impact model
    randomgrowth_l=basegr #levels impact model
    for(j in 2:time){
      randomgrowth_g=c(randomgrowth_g, basegr+(randomtemp[j]*coef)+rnorm(1,mean=0,sd=growthsd))
      randomgrowth_l=c(randomgrowth_l, basegr+(randomtemp[j]-randomtemp[j-1])*coef+rnorm(1,mean=0,sd=growthsd))
    }
    dataset=data.frame(years=1:time,temp=randomtemp,g=randomgrowth_g,l=randomgrowth_l)
    for(k in 1:length(agg)){
      period=rep(1:(time/agg[k]),each=agg[k])
      dataset=cbind(dataset,period)
      aggdata=dataset%>%group_by(period)%>%summarise(temp=mean(temp),g=mean(g),l=mean(l))
      aggdata$templag=c(NA,aggdata$temp[1:(dim(aggdata)[1]-1)])
      mod_g=lm(g~temp+templag,data=aggdata)$coef[2]
      mod_l=lm(l~temp+templag,data=aggdata)$coef[2]
      dataset=dataset[,-which(colnames(dataset)=="period")]
      sims[i,k,]=c(mod_g,mod_l)

      #SMA (moving avg)
      tempsma <- SMA(randomtemp, n=agg[k])
      gsma <- SMA(randomgrowth_g, n=agg[k])
      lsma <- SMA(randomgrowth_l, n=agg[k])
      datasetsma <- data.frame(cbind(tempsma,gsma,lsma))
      datasetsma$templag=c(NA,datasetsma$temp[1:(dim(datasetsma)[1]-1)])
      mod_gsma=lm(gsma~tempsma,data=datasetsma)$coef[2]
      mod_lsma=lm(lsma~tempsma,data=datasetsma )$coef[2]
      simssma[i,k,]=c(mod_gsma,mod_lsma)

      #Bxter-King Filter
      #trend
      if (k==1){simsbktrend[i,k,]=c(mod_g,mod_l)
        }else{
          tempts.bk <- bkfilter(randomtemp, pl=p_l[k],pu=p_u[k])
          growth.bk <- bkfilter(randomgrowth_g, pl=p_l[k],pu=p_u[k])
          level.bk <- bkfilter(randomgrowth_l, pl=p_l[k],pu=p_u[k])
          #bktrend <- data.frame(temp = unclass(tempts.bk$trend), growth = unclass(growth.bk$trend),
          #  level = unclass(level.bk$trend))
          bktrend <- data.frame(temp = randomtemp, growth = unclass(growth.bk$trend),
            level = unclass(level.bk$trend))
          bktrend$templag=c(NA,bktrend$temp[1:(dim(bktrend)[1]-1)])
          names(bktrend) <- c("temp","growth","levels","templag")
          mod_gbktrend=lm(growth~temp+templag,data=bktrend)$coef[2]
          mod_lbktrend=lm(levels~temp+templag,data=bktrend)$coef[2]
          simsbktrend[i,k,]=c(mod_gbktrend,mod_lbktrend)
      }
      #cycle
      if (k==1){simsbkcycle[i,k,]=c(0,0)
        }else{
          tempts.bk <- bkfilter(randomtemp, pl=p_l[k],pu=p_u[k])
          growth.bk <- bkfilter(randomgrowth_g, pl=p_l[k],pu=p_u[k])
          level.bk <- bkfilter(randomgrowth_l, pl=p_l[k],pu=p_u[k])
          bkcycle <- data.frame(temp = unclass(tempts.bk$cycle), growth = unclass(growth.bk$cycle),
            level = unclass(level.bk$cycle))
          bkcycle$templag=c(NA,bkcycle$temp[1:(dim(bkcycle)[1]-1)])
          names(bkcycle) <- c("temp","growth","levels","templag")
          mod_gbkcycle=lm(growth~temp,data=bkcycle)$coef[2]
          mod_lbkcycle=lm(levels~temp,data=bkcycle)$coef[2]
          simsbkcycle[i,k,]=c(mod_gbkcycle,mod_lbkcycle)
        }
      

    }
    if(i%%50==0) print(i)
  }

  library('reshape2')
  simsmean=apply(sims,MARGIN=c(2,3),FUN=mean)
  simssd=apply(sims,MARGIN=c(2,3),FUN=sd)
  colnames(simsmean)=c("Growth","Level");rownames(simsmean)=agg
  simsdat=cbind(melt(simsmean),melt(simssd)[,3])
  colnames(simsdat)=c("AggregationPeriod","ImpactType","MeanEffect","SDEffect")
  theme_set(theme_bw(base_size = 20))
  a=ggplot(simsdat,aes(x=factor(AggregationPeriod),y=MeanEffect*100,col=ImpactType,group=ImpactType))+geom_line(lwd=1.25)
  a=a+geom_point(size=4)+geom_errorbar(aes(ymin=(MeanEffect-1.96*SDEffect)*100,ymax=(MeanEffect+1.96*SDEffect)*100),width=0.1,lwd=1.25)
  a=a+geom_hline(yintercept=0,lwd=1.5,lty=3)
  a=a+labs(x="Aggregation Interval (years)",y="Estimated Growth Impact\n(% per Degree)",col="Impact Model")
  a=a+scale_color_manual(values=c("#d3818c","#7375a4")) + ggtitle("Aggregation")
  A = a

  
   #ggarrange(A,b, common.legend = T)
  simsbktrendmean=apply(simsbktrend,MARGIN=c(2,3),FUN=mean)
  simsbktrendsd=apply(simsbktrend,MARGIN=c(2,3),FUN=sd)
  colnames(simsbktrendmean)=c("Growth","Level");rownames(simsbktrendmean)=paste(p_l,p_u, sep = "-")
  simsbktrenddat=cbind(melt(simsbktrendmean),melt(simsbktrendsd)[,3])
  colnames(simsbktrenddat)=c("AggregationPeriod","ImpactType","MeanEffect","SDEffect")
  theme_set(theme_bw(base_size = 20))
  a=ggplot(simsbktrenddat,aes(x=factor(AggregationPeriod),y=MeanEffect*100,col=ImpactType,group=ImpactType))+geom_line(lwd=1.25)
  a=a+geom_point(size=4)+geom_errorbar(aes(ymin=(MeanEffect-1.96*SDEffect)*100,ymax=(MeanEffect+1.96*SDEffect)*100),width=0.1,lwd=1.25)
  a=a+geom_hline(yintercept=0,lwd=1.5,lty=3)
  a=a+labs(x="Filter frequencies (years)",y="Estimated Growth Impact\n(% per Degree)",col="Impact Model")
  a=a+scale_color_manual(values=c("#d3818c","#7375a4")) + ggtitle("BK filtered trend") + ylim(-10,5)
  c = a

  ggarrange(A,c, common.legend = T, labels=c("a)","b)"))
  
  simssmamean=apply(simssma,MARGIN=c(2,3),FUN=mean)
  simssmasd=apply(simssma,MARGIN=c(2,3),FUN=sd)
  colnames(simssmamean)=c("Growth","Level");rownames(simssmamean)=agg
  simssmadat=cbind(melt(simssmamean),melt(simssmasd)[,3])
  colnames(simssmadat)=c("AggregationPeriod","ImpactType","MeanEffect","SDEffect")
  theme_set(theme_bw(base_size = 20))
  a=ggplot(simssmadat,aes(x=factor(AggregationPeriod),y=MeanEffect*100,col=ImpactType,group=ImpactType))+geom_line(lwd=1.25)
  a=a+geom_point(size=4)+geom_errorbar(aes(ymin=(MeanEffect-1.96*SDEffect)*100,ymax=(MeanEffect+1.96*SDEffect)*100),width=0.1,lwd=1.25)
  a=a+geom_hline(yintercept=0,lwd=1.5,lty=3)
  a=a+labs(x="Average Window (years)",y="Estimated Growth Impact\n(% per Degree)",col="Impact Model")
  a=a+scale_color_manual(values=c("#d3818c","#7375a4"))  + ggtitle("Moving average (with lag)")
  b = a

  
 

  
  simsbkcyclemean=apply(simsbkcycle,MARGIN=c(2,3),FUN=mean)
  simsbkcyclesd=apply(simsbkcycle,MARGIN=c(2,3),FUN=sd)
  colnames(simsbkcyclemean)=c("Growth","Level");rownames(simsbkcyclemean)=paste(p_l,p_u, sep = "-")
  simsbkcycledat=cbind(melt(simsbkcyclemean),melt(simsbkcyclesd)[,3])
  colnames(simsbkcycledat)=c("AggregationPeriod","ImpactType","MeanEffect","SDEffect")
  theme_set(theme_bw(base_size = 20))
  a=ggplot(simsbkcycledat,aes(x=factor(AggregationPeriod),y=MeanEffect*100,col=ImpactType,group=ImpactType))+geom_line(lwd=1.25)
  a=a+geom_point(size=4)+geom_errorbar(aes(ymin=(MeanEffect-1.96*SDEffect)*100,ymax=(MeanEffect+1.96*SDEffect)*100),width=0.1,lwd=1.25)
  a=a+geom_hline(yintercept=0,lwd=1.5,lty=3)
  a=a+labs(x="Filter frequencies (years)",y="Estimated Growth Impact\n(% per Degree)",col="Impact Model")
  a=a+scale_color_manual(values=c("#d3818c","#7375a4")) + ggtitle("BK filtered cycles")
  d = a


  library('ggpubr')
  ggarrange(A,b,c,d, common.legend = T)
