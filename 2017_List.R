###2017
###Dry DEP LISTVYANKA

#__________________

library(RNetCDF)
fl <- open.nc('D:/DryDepPrim/ERA20162017_List.nc') #FarEast
print.nc(fl)
data <- read.nc(fl)
close.nc(fl)

RH_grid <- da$r*0.00145247894548357+68.5666107966601

zust_grid <- data$zust*9.71179389179416e-06+0.327385591362003 #friction velocity, [m/s]

Ts_grid <- data$stl1*0.000740411824287835+275.594150362936 #soil temperature,[K]

##What about heat fluxes? J-W per second hour... how to covert one from ERA to one from the model?
## H <- data$sshf*24.4170417957365-267773.208520898  #surface sensible heat flux J/m^2

Tp_grid <- data$d2m*0.00080037451349952+265.659702656982  #dewpoint , [K]

SW_grid <- data$swvl1*4.92883289214858e-06+0.278993317720471 #soil water, [m^3/m^3]

T2_grid <- data$t2m*0.000888186652353347+271.974668821713 #temperature at 2m, [K]

u10_grid <- data$u10*0.000305494540992825+2.4124677971631 #wind, [m/s]
v10_grid <- data$v10*0.000243188536914966-2.01869916995205
U_grid <- sqrt((u10_grid)^2+(v10_grid)^2) #wind speed, [m/s]

p_grid <- data$sp*0.114372720614042+95436.5990636397 #surface pressure, [Pa]

f_cloud_grid <- data$tcc*1.52594875864068e-05+0.499992370256207 #cloud fraction, [0-1]
f_cloud_grid[f_cloud_grid>1] <- 1

LAI_grid <- data$lai_hv*7.48277436320251e-05+2.45180584784693  #Leaf Area Index, [m^2/m^2]

sol_rad_grid <- ((data$ssrd)*49.8557978423085+1633575.07210108)/3600 #solar radiation for stability class

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

#for LIstvyanka, elevation=646m

for (i in 1:length(time)) {
  for (j in 1:length(lon)) {
    for (k in 1:length(lat)) {
      jd <- julianDay(unlist(month.day.year(date[i]))[3], unlist(month.day.year(date[i]))[1], unlist(month.day.year(date[i]))[2], hour[i], 0, 0, 0, 0) 
      s <- solarPosition(jd,lon[j],lat[k],delta_t=32.184, elev=646,temp=T2_grid[j,k,i]-273,pres=p_grid[j,k,i]/100)[1] 
      SZA_grid[j,k,i] <- s 
    } #k
  } #j
} #i

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#DO FROM THIS MOMENT AFTER extract_from_nc

omega <- time #omega - hour angle
omega[hour==0] <- pi
omega[(hour==6)|(hour==18)] <- pi/2
omega[hour==12] <- 0 
gamma <- 2+pi*(dn-1)/365




Pr<-0.72 #Prandtl number
k<-0.41 # Karman const

#nu<- 0.15 # cm^2 s^-1 at 20C kinetic viscosity of air 
#// nu=0.0093*TCels+1.32


#--------------------------------------------------------------------------------------
#END OF INPUTS
#START OF CALCULATIONS
#--------------------------------------------------------------------------------------


#initialize deposition velocities

duration <- 1465:2924 #from 1465 tp 2924 = 1460 points=2017 year

V_LIST_NH3_2017 <- array(0,dim=c(length(lon),length(lat),length(duration)))
V_LIST_O3_2017 <- array(0,dim=c(length(lon),length(lat),length(duration)))
V_LIST_SO2_2017 <- array(0,dim=c(length(lon),length(lat),length(duration)))


dnum <- 365 #day numbers ro the year of calculations - currently 2016 year, 366 days
#Change RH,change parameters from Iqbal in f_phen


for (lon_num in 1:length(lon)) { #huge cycle on spatial grid
for (lat_num in 1:length(lat)) {

#reassignment of grid variables to point variables


zust <- zust_grid[lon_num,lat_num,duration]
Ts <- Ts_grid[lon_num,lat_num,duration]
Tp <- Tp_grid[lon_num,lat_num,duration]
SW <- SW_grid[lon_num,lat_num,duration]
T2 <- T2_grid[lon_num,lat_num,duration] #take T2 from point EANET measurements
u10 <- u10_grid[lon_num,lat_num,duration]
v10 <- v10_grid[lon_num,lat_num,duration]
U <- U_grid [lon_num,lat_num,duration]
p <- p_grid[lon_num,lat_num,duration]
f_cloud <- f_cloud_grid[lon_num,lat_num,duration]
LAI <- LAI_grid[lon_num,lat_num,duration]
SZA <- SZA_grid[lon_num,lat_num,duration]
sol_rad <- sol_rad_grid[lon_num,lat_num,duration]
RH <- RH_grid[lon_num,lat_num,duration]


##---------------------------T2 from EANET 
##suppose min at 00,06 max at 12, mean at 18

T2_EANET <- rep(NULL,times=365*4)

f <- file('D:/ERA5/T2List2017.txt','r')
d <- read.table(f)
close(f)

for (i in 1:365) {
l <- d[i,3]
T2_EANET[((i-1)*4+1):(i*4-2)] <- rep(l,times=2)
l <- d[i,1]
T2_EANET[(i*4-1)] <- l 
l <- (d[i,1]+d[i,2])/2
T2_EANET[(i*4)] <- l
}

T2 <- T2_EANET+273
#T2 <- T2+0.001
##-----------------------------------------------
#nu <- 0.0093*(T2-273)+1.32

#D_H2O<-0.242 # s cm**-2 molecualr diffusivity of H2O

#D_r_SO2<-1.9 #molecualr diffusivity of gas=D_H2O/D_r_GAS
#D_r_NH3<-1.0 #coef from EMEP model 

#D_SO2<-D_H2O/D_r_SO2
#D_NH3<- 0.259 #D_H2O/D_r_NH3

Sc_SO2 <- 1.22 #Schmidt number for SO2
Sc_NH3 <- 0.57 #for NH3

Pr <- 0.72 #Prandtl number

R_b_NH3<-(2)/(k*zust)*(Sc_NH3/Pr)^(2/3)
R_b_SO2<-(2)/(k*zust)*(Sc_SO2/Pr)^(2/3)

#-------------------------
#Ra calculations

h <- 40 #[m], height of the trees
d.h <- 0.78*h # displacement height
z0 <- 0.07*h # roughness length
Z <- 50 #z0+1 #file:///C:/Users/Alisa/AppData/Local/Temp/aer_scid.pdf 10 #[m], reference height
rho <- 1.2 #[kg/m^3], air density
C_p <- 1005 #[J kg^-1 K^-1]

#R_ah <- rho*C_p*(Ts-T2)/H # same for heat and momentum transfer in neutral conditions
#R_ah <- U/(zust^2)#

##Pasquil classes

class <- array(0,dim=length(duration))
a <- array(0,dim=length(duration))
b <- array(0,dim=length(duration))

class[(f_cloud>=0.4375)&(U<2)] <- 6
class[(f_cloud>=0.4375)&(U<3)&(U>=2)] <- 5
class[(f_cloud>=0.4375)&(U>=3)] <- 4
class[(f_cloud<0.4375)&(U<2)] <- 6
class[(f_cloud<0.4375)&(U>=2)&(U<3)] <- 6
class[(f_cloud<0.4375)&(U>=3)&(U<5)] <- 5
class[(f_cloud<0.4375)&(U>=5)] <- 4
class[(U<2)&(sol_rad>700)&(SZA<90)] <- 1
class[(U>=2)&(U<3)&(sol_rad>700)&(SZA<90)] <- 1
class[(U>=3)&(U<5)&(sol_rad>700)&(SZA<90)] <- 2
class[(U>=5)&(sol_rad>700)&(SZA<90)] <- 3
class[(U<2)&(sol_rad<=700)&(sol_rad>300)&(SZA<90)] <- 2#1
class[(U>=2)&(U<3)&(sol_rad<=700)&(sol_rad>300)&(SZA<90)] <- 2
class[(U>=3)&(U<5)&(sol_rad<=700)&(sol_rad>300)&(SZA<90)] <- 2
class[(U>=3)&(U<5)&(sol_rad<=700)&(sol_rad>300)&(SZA<90)] <- 3
class[(U>=5)&(sol_rad<=700)&(sol_rad>300)&(SZA<90)] <- 4
class[(U<2)&(sol_rad<=300)&(SZA<90)] <- 2
class[(U<5)&(U>=2)&(sol_rad<=300)&(SZA<90)] <- 3
class[(U>=5)&(sol_rad<=300)&(SZA<90)] <- 4

#Monin-Obuhov Length

L_1 <- array(0,dim=length(duration))

a[class==1] <- -0.096
a[class==2] <- -0.037
a[class==3] <- -0.002
a[class==4] <- 0
a[class==5] <- 0.004
a[class==6] <- 0.035
b[class==1] <- 0.029
b[class==2] <- 0.029
b[class==3] <- 0.018
b[class==4] <- 0
b[class==5] <- -0.018
b[class==6] <- -0.040


L_1 <- (a+b*log(z0))

#y_h <- array(0,dim=length(duration))
#y_h[((Z*L_1>-1)&(Z*L_1<0))] <- 2*log(0.5*(1+(1-16*Z*L_1[((Z*L_1>-1)&(Z*L_1<0))])^0.5)) #unstable conditions
#y_h[((Z*L_1<1)&(Z*L_1>0))] <- -5*Z*L_1[((Z*L_1<1)&(Z*L_1>0))] #stable conditions

#R_ah <- log(d/z0)/k/zust
#R_ah <- (log(Z/z0)-y_h)/(k*zust)

R_ah <- array(0,dim=length(duration))
psi_max <- array(0,dim=length(duration))

#R_ah[L_1>=0] <- (log(Z/z0)+5*Z*L_1[L_1>=0])/(k*zust[L_1>=0]) #stable
#R_ah[L_1<0] <- (log((sqrt(1-16*Z*L_1[L_1<0])-1)*(sqrt(1-16*z0*L_1[L_1<0])+1)/((sqrt(1-16*Z*L_1[L_1<0])+1)*(sqrt(1-16*z0*L_1[L_1<0])-1))))/(k*zust[L_1<0])
#R_ah[L_1<0] <- (log(Z/z0)-2*log(0.5*(1+(1-16*Z*L_1[L_1<0])^0.5))+2*log(0.5*(1+(1-16*z0*L_1[L_1<0])^0.5)))/(k*zust[L_1<0])

zeta <- (Z-d.h)*L_1
x <- (1-16*zeta)^0.25
psi_m <- log(((1+x)^2)*(1+x^2)/8)-2*atan(x)+pi/2
psi_h <- 2*log((1+x^2)/2)
Y_h <- log((Z-d.h)/z0)
Y_m <- Y_h
R_ah[L_1<0] <- (Y_h-psi_h[L_1<0])*(Y_m-psi_m[L_1<0])/((k^2)*zust[L_1<0])
psi_max[L_1>=0] <- max(-5, -5*zeta[L_1>=0])
R_ah[L_1>=0] <- ((Y_h-psi_max[L_1>=0])^2)/(zust[L_1>=0]*k^2)



#----------------------------
# Rc calculations

gm_max <- 150
g_max <- gm_max/41000

  #f_phen
phi_a <- 0
phi_b <- 0
phi_c <- 1
phi_d <- 0
phi_e <- 20
phi_f <- 30
phi_As <- 0
phi_Ae <- 0

dn_SGS <- 31+28+31+20+20 # day number start of growing season 20 may
dn_EGS <- 31+28+31+30+31+30+31+31+10 # day number end of growing season 10 sept

Astart <- dn_SGS+phi_As
Aend <- dn_EGS+phi_Ae

#dn1 <- current day number
dn1 <- dn[1465:2924]-366
f_phen <- rep(NULL,times=length(dn1))
f_phen[(dn1<=dn_SGS)|(dn1>dn_EGS)] <- 0
f_phen[(dn1>dn_SGS)&(dn1<=Astart)] <- phi_a
f_phen[(dn1>Astart)&(dn1<=Astart+phi_e)] <- phi_b+(phi_c-phi_b)*(dn1[(dn1>Astart)&(dn1<=Astart+phi_e)]-Astart)/phi_e
f_phen[(dn1>Astart+phi_e)&(dn1<=Aend-phi_f)] <- phi_c
f_phen[(dn1>Aend-phi_f)&(dn1<=Aend)] <- phi_d+(phi_c-phi_d)*(Aend-dn1[(dn1>Aend-phi_f)&(dn1<=Aend)])/phi_f
f_phen[(dn1>Aend)&(dn1<=dn_EGS)] <- phi_d



  #f_T
T_min <- 0
T_max <- 35
T_opt <- 20
beta <- (T_max-T_opt)/(T_opt-T_min)

f_T <- (T2-273-T_min)/(T_opt-T_min)*(((T_max-T2+273)/(T_max-T_opt))^beta) # T2 air temperature at 2m

  #relative humidity

#RH <- 1-0.05*(T2-Tp) #Tp dewpoint temperature
T <- T2-273
Tp1 <- Tp-273
Ps <- T
Ps_d <- T
#RH <- T

Ps_d[T>=0] <- 6.1121*exp((18.678-(Tp1[T>=0]/234.5))*Tp1[T>=0]/(257.14+Tp1[T>=0]))
Ps[T>=0] <- 6.1121*exp((18.678-(T[T>=0]/234.5))*T[T>=0]/(257.14+T[T>=0]))

Ps_d[T<0] <- 6.1115*exp((23.036-(Tp1[T<0]/333.7))*Tp1[T<0]/(279.82+Tp1[T<0]))
Ps[T<0] <- 6.1115*exp((23.036-(T[T<0]/333.7))*T[T<0]/(279.82+T[T<0]))

#RH <- Ps_d*100/Ps

###------------------------------------RH from EANET
f <- file('D:/ERA5/RHList2017.txt','r')
d <- read.table(f)
close(f)

for (i in 1:365) {
l <- d[i,1]
RH[((i-1)*4+1):(i*4)] <- rep(l,times=4) 
} #i
#RH <- RH+0.001
###------------------RH from EANET

  #f_D
f_min <- 0.1
D_max <- 1000 #[Pa]
D_min <-3250 #[Pa]

aa <- -1.044e4
bb <- -11.29
cc <- -2.7e-2
dd <- 1.289e-5
ee <- -2.478e-9
ff <- 6.456

#Trankin=Tcelcius*1.8+491.67

Tr <- (T2-273)*1.8+491.67
vp_sat <- exp(aa/Tr+bb+cc*Tr+dd*Tr^2+ee*Tr^3+ff*log(Tr))  #T=Trankin
D <- vp_sat*(1-RH/100)*6894.75 #6894.75 to convert from PSI to Pascal

f_D <- f_min+(1-f_min)*(D_min-D)/(D_min-D_max)

  #f_SW

#soil type fine

PWP <-0.279 #[m**3/m**3]
FC <- 0.448 #[m**3/m**3]

SMI <- (SW-PWP)/(FC-PWP) # soil moisture index
f_SW <- SW

f_SW[SMI>=0.5] <- 1
f_SW[SMI<0.5] <- 2*SMI[SMI<0.5]


  #f_light
A_iq <- c(1202,1187,1164,1130,1106,1092,1093,1107,1136,1136,1190,1204) # [W m**-2], parameters from Iqbal
B_iq <- c(0.141,0.142,0.149,0.164,0.177,0.185,0.186,0.182,0.165,0.152,0.144,0.141)
C_iq <- c(0.103,0.104,0.109,0.120,0.130,0.137,0.138,0.134,0.121,0.111,0.106,0.103)

Aiq <- rep(NULL, times=1460)
Aiq[1:(w<-31*4)] <- A_iq[1]
Aiq[(w+1):(w<-(w+28*4))] <- A_iq[2]
Aiq[(w+1):(w<-(w+31*4))] <- A_iq[3]
Aiq[(w+1):(w<-(w+30*4))] <- A_iq[4]
Aiq[(w+1):(w<-(w+31*4))] <- A_iq[5]
Aiq[(w+1):(w<-(w+30*4))] <- A_iq[6]
Aiq[(w+1):(w<-(w+31*4))] <- A_iq[7]
Aiq[(w+1):(w<-(w+31*4))] <- A_iq[8]
Aiq[(w+1):(w<-(w+30*4))] <- A_iq[9]
Aiq[(w+1):(w<-(w+31*4))] <- A_iq[10]
Aiq[(w+1):(w<-(w+30*4))] <- A_iq[11]
Aiq[(w+1):(w<-(w+31*4))] <- A_iq[12]
Biq <- rep(NULL, times=1460)
Biq[1:(w<-31*4)] <- B_iq[1]
Biq[(w+1):(w<-(w+28*4))] <- B_iq[2]
Biq[(w+1):(w<-(w+31*4))] <- B_iq[3]
Biq[(w+1):(w<-(w+30*4))] <- B_iq[4]
Biq[(w+1):(w<-(w+31*4))] <- B_iq[5]
Biq[(w+1):(w<-(w+30*4))] <- B_iq[6]
Biq[(w+1):(w<-(w+31*4))] <- B_iq[7]
Biq[(w+1):(w<-(w+31*4))] <- B_iq[8]
Biq[(w+1):(w<-(w+30*4))] <- B_iq[9]
Biq[(w+1):(w<-(w+31*4))] <- B_iq[10]
Biq[(w+1):(w<-(w+30*4))] <- B_iq[11]
Biq[(w+1):(w<-(w+31*4))] <- B_iq[12]
Ciq <- rep(NULL, times=1460)
Ciq[1:(w<-31*4)] <- C_iq[1]
Ciq[(w+1):(w<-(w+28*4))] <- C_iq[2]
Ciq[(w+1):(w<-(w+31*4))] <- C_iq[3]
Ciq[(w+1):(w<-(w+30*4))] <- C_iq[4]
Ciq[(w+1):(w<-(w+31*4))] <- C_iq[5]
Ciq[(w+1):(w<-(w+30*4))] <- C_iq[6]
Ciq[(w+1):(w<-(w+31*4))] <- C_iq[7]
Ciq[(w+1):(w<-(w+31*4))] <- C_iq[8]
Ciq[(w+1):(w<-(w+30*4))] <- C_iq[9]
Ciq[(w+1):(w<-(w+31*4))] <- C_iq[10]
Ciq[(w+1):(w<-(w+30*4))] <- C_iq[11]
Ciq[(w+1):(w<-(w+31*4))] <- C_iq[12]

p0 <- 101.3*1e3 #[Pa], standard sea-level pressure
al <- 0.06

#theta calculations

theta <- SZA*pi/180
cos_theta <- cos(theta)
cos_theta[cos_theta<0] <- 0 

Tk <- 1-0.7*(f_cloud)^3.4 # f_cloud - cloudiness

IN_dir <- Aiq*Tk*exp(-Biq*p/cos_theta/p0)
I_dir <- IN_dir*cos_theta
I_diff <- Ciq*I_dir

alpha <- pi/3 #inclination of leaves
LAI_sun <- (1-exp(-0.5*LAI/cos_theta))*cos_theta/cos(alpha)
LAI_shade <- LAI-LAI_sun
Ishade_PAR <- I_diff*exp(-0.5*(LAI)^0.7)+0.07*I_dir*(1.1-0.1*LAI)*exp(-cos_theta)
Isun_PAR <- I_dir*cos(alpha)/cos_theta+Ishade_PAR
Isun_PAR[cos_theta==0] <- 0
Ishade_PAR[cos_theta==0] <- 0
f_sun <- 1-exp(-al*Isun_PAR)
f_shade <- 1-exp(-al*Ishade_PAR)
fsun_LAI <- LAI_sun/LAI
f_light <- fsun_LAI*f_sun+(1-fsun_LAI)*f_shade


#------------------------------------------------------------------------------
#g_sto

maxx <- rep(NULL, times=length(f_T))
for (i in 1:length(f_T)) {
  maxx[i] <- max(f_min,f_T[i]*f_D[i]*f_SW[i])
}
g_sto <- g_max*f_phen*f_light*maxx

#-------------------------------------------------------------------------------
#non-stomatal NH3 SO2,NH3 [nmol/m**3]
beta_ns <- 1/22
SO2 <- rep(NULL,times=1460)

SO2[1:(w<-(8*4))] <- 200.49
SO2[(w+1):(w<-(w+10*4))] <- 144.4
SO2[(w+1):(w<-(w+7*4))] <- 289.91
SO2[(w+1):(w<-(w+7*4))] <- 129.42
SO2[(w+1):(w<-(w+7*4))] <- 265.14
SO2[(w+1):(w<-(w+7*4))] <- 457.46
SO2[(w+1):(w<-(w+11*4))] <- 114.05
SO2[(w+1):(w<-(w+3*4))] <- 254.17
SO2[(w+1):(w<-(w+7*4))] <- 141.27
SO2[(w+1):(w<-(w+7*4))] <- 131.68
SO2[(w+1):(w<-(w+7*4))] <- 161.71
SO2[(w+1):(w<-(w+7*4))] <- 93.61
SO2[(w+1):(w<-(w+7*4))] <- 113.96
SO2[(w+1):(w<-(w+7*4))] <- 42.83
SO2[(w+1):(w<-(w+7*4))] <- 54.94
SO2[(w+1):(w<-(w+7*4))] <- 74.87
SO2[(w+1):(w<-(w+7*4))] <- 84.46
SO2[(w+1):(w<-(w+7*4))] <- 51.14
SO2[(w+1):(w<-(w+7*4))] <- 20.53
SO2[(w+1):(w<-(w+7*4))] <- 80.68
SO2[(w+1):(w<-(w+7*4))] <- 19.17
SO2[(w+1):(w<-(w+7*4))] <- 17.38
SO2[(w+1):(w<-(w+7*4))] <- 7.27
SO2[(w+1):(w<-(w+7*4))] <- 74.34
SO2[(w+1):(w<-(w+7*4))] <- 23.27
SO2[(w+1):(w<-(w+7*4))] <- 51.56
SO2[(w+1):(w<-(w+7*4))] <- 25.08
SO2[(w+1):(w<-(w+7*4))] <- 51.94
SO2[(w+1):(w<-(w+5*4))] <- 35.69
SO2[(w+1):(w<-(w+8*4))] <- 74.96
SO2[(w+1):(w<-(w+8*4))] <- 207.86
SO2[(w+1):(w<-(w+7*4))] <- 52.77
SO2[(w+1):(w<-(w+7*4))] <- 64.18
SO2[(w+1):(w<-(w+7*4))] <- 30.63
SO2[(w+1):(w<-(w+7*4))] <- 43.81
SO2[(w+1):(w<-(w+6*4))] <- 33.75
SO2[(w+1):(w<-(w+7*4))] <- 11.42
SO2[(w+1):(w<-(w+7*4))] <- 56.79
SO2[(w+1):(w<-(w+7*4))] <- 100.41
SO2[(w+1):(w<-(w+7*4))] <- 77.17
SO2[(w+1):(w<-(w+7*4))] <- 99.45
SO2[(w+1):(w<-(w+7*4))] <- 89.37
SO2[(w+1):(w<-(w+14*4))] <- 170.11
SO2[(w+1):(w<-(w+7*4))] <- 113.77
SO2[(w+1):(w<-(w+7*4))] <- 340.63
SO2[(w+1):(w<-(w+7*4))] <- 230.95
SO2[(w+1):(w<-(w+7*4))] <- 294.6
SO2[(w+1):(w<-(w+7*4))] <- 201.53
SO2[(w+1):(w<-(w+7*4))] <- 225.57
SO2[(w+1):(w<-(w+7*4))] <- 226.79
SO2[(w+1):(w<-(w+5*4))] <- 461.45



NH3 <- rep(NULL,times=1460)
NH3[1:(w<-(8*4))] <- 118.31
NH3[(w+1):(w<-(w+10*4))] <- 117.19 
NH3[(w+1):(w<-(w+7*4))] <- 138.51
NH3[(w+1):(w<-(w+7*4))] <- 31.16
NH3[(w+1):(w<-(w+7*4))] <- 21.57
NH3[(w+1):(w<-(w+7*4))] <- 30.41
NH3[(w+1):(w<-(w+11*4))] <- 33.53
NH3[(w+1):(w<-(w+3*4))] <- 320.48
NH3[(w+1):(w<-(w+7*4))] <- 82.83
NH3[(w+1):(w<-(w+7*4))] <- 70.59
NH3[(w+1):(w<-(w+7*4))] <- 23.45
NH3[(w+1):(w<-(w+7*4))] <- 11.43
NH3[(w+1):(w<-(w+7*4))] <- 11.33
NH3[(w+1):(w<-(w+7*4))] <- 15.24
NH3[(w+1):(w<-(w+7*4))] <- 97.63
NH3[(w+1):(w<-(w+7*4))] <- 70.18
NH3[(w+1):(w<-(w+7*4))] <- 76.75
NH3[(w+1):(w<-(w+7*4))] <- 68.73
NH3[(w+1):(w<-(w+7*4))] <- 91.32
NH3[(w+1):(w<-(w+7*4))] <- 4.74
NH3[(w+1):(w<-(w+7*4))] <- 56.05
NH3[(w+1):(w<-(w+7*4))] <- 4.9
NH3[(w+1):(w<-(w+7*4))] <- 4.92
NH3[(w+1):(w<-(w+7*4))] <- 257.74
NH3[(w+1):(w<-(w+7*4))] <- 377.96
NH3[(w+1):(w<-(w+7*4))] <- 331.15
NH3[(w+1):(w<-(w+7*4))] <- 205.36
NH3[(w+1):(w<-(w+7*4))] <- 260.07
NH3[(w+1):(w<-(w+5*4))] <- 465.85
NH3[(w+1):(w<-(w+8*4))] <- 219.11
NH3[(w+1):(w<-(w+8*4))] <- 307.22
NH3[(w+1):(w<-(w+7*4))] <- 148.44
NH3[(w+1):(w<-(w+7*4))] <- 127.17
NH3[(w+1):(w<-(w+7*4))] <- 145.37
NH3[(w+1):(w<-(w+7*4))] <- 130.5
NH3[(w+1):(w<-(w+6*4))] <- 239.39
NH3[(w+1):(w<-(w+7*4))] <- 286.47
NH3[(w+1):(w<-(w+7*4))] <- 23.06
NH3[(w+1):(w<-(w+7*4))] <- 50.14
NH3[(w+1):(w<-(w+7*4))] <- 13.14
NH3[(w+1):(w<-(w+7*4))] <- 32.7
NH3[(w+1):(w<-(w+7*4))] <- 8.22
NH3[(w+1):(w<-(w+14*4))] <- 34.09
NH3[(w+1):(w<-(w+7*4))] <- 253.82
NH3[(w+1):(w<-(w+7*4))] <- 528.46
NH3[(w+1):(w<-(w+7*4))] <- 56.51
NH3[(w+1):(w<-(w+7*4))] <- 7.34
NH3[(w+1):(w<-(w+7*4))] <- 7.34
NH3[(w+1):(w<-(w+7*4))] <- 15.54
NH3[(w+1):(w<-(w+7*4))] <- 23.6
NH3[(w+1):(w<-(w+5*4))] <- 16.68



alfa_SN <- SO2/NH3#nh3.pr#NH3 #molar "acidity ratio"

R_NH3 <- rep(NULL,times=length(g_sto))
R_NH3[T>0] <- 1/22*10*log10((T[T>0]+2)*exp((100-RH[T>0])/7))*10^(-1.1099*alfa_SN[T>0]+1.6769)
R_NH3[(T>-5)&(T<=0)] <- 100
R_NH3[T<=-5] <- 500 

G_NH3_ns <- 1/R_NH3


#Rc

R_c_NH3 <- 1/(g_sto*LAI+G_NH3_ns)

#------------------------------------------------------------------
#Ozone G_O3_ns

F_t <- exp(-0.2*(1+T)) #temperature correction factor
F_t[F_t<1] <- 1
F_t[F_t>2] <- 2

r_ext <- 2500*F_t

SAI <- LAI+1 #for forests - I suppose there is a forest around Primorskaya station

b <- 14 #[s^-1] empirical constant

R_inc <- b*SAI*h/zust

#in the abscence of snow f_snow <-0

f_snow <- 0 #should be (snow depth)/(expected maximum snow depth)

R_O3_gs_base <- 200   #[s m^-1] from the table EMEP model

R_O3_snow <- 2000 #[s m^-1]

R_O3_gs <- 1/((1-2*f_snow)/(F_t*R_O3_gs_base)+2*f_snow/(R_O3_snow))

G_O3_ns <- SAI/r_ext+1/(R_inc+R_O3_gs)

R_c_O3 <- 1/(g_sto*LAI+G_O3_ns)

#------------------------------------------------------------------------------
#SO2
R_SO2_ns <- rep(NULL,times=length(g_sto))
R_SO2_ns[T>0] <- 11.84*exp(1.1*alfa_SN[T>0])*(RH[T>0]/100)^(-1.67)
R_SO2_ns[(T>-5)&(T<=0)] <- 100
R_SO2_ns[T<=-5] <- 500

R_c_SO2 <- 1/(g_sto*LAI+1/R_SO2_ns)

#---------------------------------
#DRY DEP VELOCITY

V_LIST_NH3_2017[lon_num,lat_num,] <- 1/(R_ah+R_b_NH3+R_c_NH3) # deposition velocity [m/s]
#V_LIST_O3_2017[lon_num,lat_num,] <- 1/(R_ah+R_b+R_c_O3)
V_LIST_SO2_2017[lon_num,lat_num,] <- 1/(R_ah+R_b_SO2+R_c_SO2)


} #lat_num
} #lon_num

#END OF HUGE CYCLE FOR GRID

#-----------------------------------
#TOTAL AMOUNT CALCULATION

#flux_NH3 <- NH3*V_NH3*6*60*60 #nmol/m**2

