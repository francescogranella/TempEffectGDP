   library(plyr)  #sql
    library(raster)
    library(rworldmap) 
    library(reshape2)
    library(ncdf4)
    library(ncdf)
    library(lfe)
    library(dplyr)
    #barro GDP Data
    madd_gdp_lmr <- read.csv("C:/Users/bastien/Documents/GitHub/
        longrungdpgrowth/Long Run GDP Growth/Data/GDP_LMR.csv")

    #UDel Climate Data
    #UDel temperature data
    #list[[1]]
    filenames <- list.files("C:/Users/bastien/Box/Long Run GDP Growth/Data/UDel Data new",pattern="air_temp.*",full.names=TRUE)
    

    pop_r=raster("C:/Users/bastien/Box/Long Run GDP Growth/Data/Population Data/glp00ag30.asc")
    fullextent=extent(-180,180,-90,90)
    pop_r=extend(pop_r,fullextent) 
    pop_extract=extract(pop_r,countriesLow,na.rm=TRUE)
    
    for (j in 1:length(filenames)){
        temp <-file(filenames[[j]],open="r")
        linn <- readLines(temp) # 85794 lines(strings), each line contains lon, lat, and 12 observations
        close(temp)
        linn1 <-trimws(linn) # Remove leading and/or trailing whitespace from each strings
        tempsplit <- strsplit(linn1, "\\s+") # split by space in each line
        df <- data.frame(matrix(unlist(tempsplit), nrow=length(linn), byrow=T),stringsAsFactors=FALSE) # list to data frame
        df <- as.data.frame(sapply(df, as.numeric)) # character to numeric
        df <- df[,c(1,2,15)]
        ras <- rasterFromXYZ(df) # dataframe to raster
        fullextent=extent(-180,180,-90,90)
        temp_r=extend(ras,fullextent) 
        temp_extract=extract(temp_r,countriesLow,na.rm=TRUE)
        popaverage <- rep(0,244)*NA
        for(i in 1:244){    
            popaverage[i] <- weighted.mean(x=temp_extract[[i]], w=pop_extract[[i]], na.rm=TRUE)
        }
        if (j==1){
            UDel_t_pop <-data.frame(countrycode=countriesLow$ISO3,UDel_t_pop=popaverage,year=j+1899)
        }else{
            UDel_t_pop <- rbind(UDel_t_pop,data.frame(countrycode=countriesLow$ISO3,UDel_t_pop=popaverage,year=j+1899))
        }
        print(paste("year=",j+1))
    }
    write.csv(spatialaverage,file=("C:/Users/bastien/Box/Long Run GDP Growth/Data/UDel Data new/spatialaverage_corrected.csv"))
    write.csv(UDel_t_pop,file=("C:/Users/bastien/Box/Long Run GDP Growth/Data/UDel Data new/popaverage_corrected.csv"))
    



    ## Precipitation
    filenames <- list.files("C:/Users/bastien/Box/Long Run GDP Growth/Data/UDel Data new",pattern="precip.*",full.names=TRUE)
    for (j in 1:length(filenames)){
        temp <-file(filenames[[j]],open="r")
        linn <- readLines(temp) # 85794 lines(strings), each line contains lon, lat, and 12 observations
        close(temp)
        linn1 <-trimws(linn) # Remove leading and/or trailing whitespace from each strings
        tempsplit <- strsplit(linn1, "\\s+") # split by space in each line
        df <- data.frame(matrix(unlist(tempsplit), nrow=length(linn), byrow=T),stringsAsFactors=FALSE) # list to data frame
        df <- as.data.frame(sapply(df, as.numeric)) # character to numeric
        df <- df[,c(1,2,15)]
        ras <- rasterFromXYZ(df) # dataframe to raster
        fullextent=extent(-180,180,-90,90)
        temp_r=extend(ras,fullextent) 
        temp_extract=raster::extract(temp_r,countriesLow,na.rm=TRUE)
        popaverage <- rep(0,244)*NA
        for(i in 1:244){    
            popaverage[i] <- weighted.mean(x=temp_extract[[i]], w=pop_extract[[i]], na.rm=TRUE)
        }
        if (j==1){
            UDel_p_pop <-data.frame(countrycode=countriesLow$ISO3,UDel_p_pop=popaverage,year=j+1899)
        }else{
            UDel_p_pop <- rbind(UDel_t_pop,data.frame(countrycode=countriesLow$ISO3,UDel_p_pop=popaverage,year=j+1899))
        }
        print(paste("year=",j+1))
    }
    write.csv(UDel_p_pop,file=("C:/Users/bastien/Box/Long Run GDP Growth/Data/UDel Data new/p_popaverage_corrected.csv"))
    
