library(RNetCDF)

uP <- NULL
uPtime <- NULL
cloud <- NULL

flnm <- dir(di<-'D:/ERA5/SurfPres/')

for (i in flnm) {

f<-open.nc(paste(di,i,sep=''))
d<-read.nc(f)
close.nc(f)

u<-d$sp
cc <- d$tcc
lat<-d$latitude
lon<-d$longitude
uP <- c(uP,u[lon==132.25,lat==43.5,]) #around Primorskaya
cloud <- c(cloud,cc[lon==132.25,lat==43.5,])
uPtime <- c(uPtime,d$time)
}

zust <- array(dim=c(length(uP),2))
zust[,1] <- uPtime
zust[,2] <- uP

f <- file('D:/ERA5/SP.txt','w')
write.table(zust,f,row.names=FALSE,col.names=FALSE)
close(f)

zust <- array(dim=c(length(uP),2))
zust[,1] <- uPtime
zust[,2] <- cloud

f <- file('D:/ERA5/CLOUD.txt','w')
write.table(zust,f,row.names=FALSE,col.names=FALSE)
close(f)




#__________________

library(RNetCDF)
fl <- open.nc('D:/DryDepPrim/ERA20162017_Prim.nc') #FarEast
print.nc(fl)
data <- read.nc(fl)
close.nc(fl)

zust_grid <- data$zust*1.13479407883108e-05+0.37856334226269 #friction velocity, [m/s]

Ts_grid <- data$stl1*0.000744095606531251+277.0618227764162 #soil temperature,[K]

##What about heat fluxes? J-W per second hour... how to covert one from ERA to one from the model?
## H <- data$sshf*24.4170417957365-267773.208520898  #surface sensible heat flux J/m^2

Tp_grid <- data$d2m*0.000888134263060333+269.188513147224  #dewpoint , [K]

SW_grid <- data$swvl1*4.92883289214858e-06+0.278993317720471 #soil water, [m^3/m^3]

T2_grid <- data$t2m*0.000888391087016751+277.857008197035 #temperature at 2m, [K]

u10_grid <- data$u10*0.000218942859423744-1.00328425322903 #wind, [m/s]
v10_grid <- data$v10*0.00032756866773605-1.87641738922278
U_grid <- sqrt((u10_grid)^2+(v10_grid)^2) #wind speed, [m/s]

p_grid <- data$sp*0.118753147269315+99779.8781234264 #surface pressure, [Pa]

f_cloud_grid <- data$tcc*1.5259536701369e-05+0.499993979557058 #cloud fraction, [0-1]
f_cloud_grid[f_cloud_grid>1] <- 1

LAI_grid <- data$lai_hv*5.09464336737133e-05+3.54980400981539  #Leaf Area Index, [m^2/m^2]

library(solarPos)
library(chron)

lon <- data$longitude
lat <- data$latitude
time <- data$time
dn <- trunc(div<-((time-1016832)/24))+1 
date <- as.Date('2016-01-01')+dn-1
hour <- (div-trunc(div))*24
#date_full <- date+hour



SZA_grid <- array(0, dim=c(length(lon),length(lat),length(time)))

#for Primorskaya, elevation=85m

for (i in 1:length(time)) {
  for (j in 1:length(lon)) {
    for (k in 1:length(lat)) {
      jd <- julianDay(unlist(month.day.year(date[i]))[3], unlist(month.day.year(date[i]))[1], unlist(month.day.year(date[i]))[2], hour[i], 0, 0, 0, 0) 
      s <- solarPosition(jd,lon[j],lat[k],delta_t=32.184, elev=85,temp=T2_grid[j,k,i]-273,pres=p_grid[j,k,i]/100)[1] 
      SZA_grid[j,k,i] <- s 
    } #k
  } #j
} #i

 
