library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(dlnm)
library(zoo)
library(splines)
library(Epi)
library(tsModel)
library(rnoaa)


############ Read in GHCN NOAA data###############################
stations <- ghcnd_stations()  ### THIS WILL TAKE A LONG TIME- AVOID RE-CALLING AFTER CREATING

View(stations %>% filter(grepl("ALBANY", name) & state=="NY" & last_year=="2020"))  ## input keywords and information to locate weather station

####loc:USW00014735####
loc=ghcnd_search(stationid="USW00014735", date_min= "2020-03-07", date_max = "2020-04-29",
                var="all")

## temperature reported in tenths of degrees Celsius (be sure to divide by 10)
## to convert C to F; use (T*1.8) + 32
loc_tmin= loc$tmin
loc_tmin= loc_tmin %>% mutate(tmin=tmin/10) %>% select(id, tmin, date)

## rainfall
#Precipitation, in tenths of mm

loc_rain=loc$prcp
loc_rain= loc_rain %>% mutate(prcp=prcp/10) %>% select(date, prcp)


## WIND
loc_wsf2=loc$wsf2 %>% select(date, wsf2)
loc_wsf5=loc$wsf5 %>% select(date, wsf5)

## Weather Type

loc_wt01= loc$wt01 %>% select(id, wt01, date)
loc_wt02= loc$wt02 %>% select(wt02, date)
loc_wt03= loc$wt03 %>% select(wt03, date)
#loc_wt06= loc$wt06 %>% select(wt06, date)
loc_wt08= loc$wt08 %>% select(wt08, date)

loc_wt=left_join(loc_wt01, loc_wt02, by="date")
loc_wt=left_join(loc_wt, loc_wt03, by="date")
loc_wt= left_join(loc_wt, loc_wt08, by="date")

loc_wt=loc_wt %>% select(date, wt01, wt02, wt03, wt08)

## combine all meterological data
loc_weather=left_join(loc_tmin, loc_rain, by="date")
loc_weather=left_join(loc_weather, loc_wsf2, by="date")
loc_weather=left_join(loc_weather, loc_wsf5, by="date")
loc_weather=left_join(loc_weather, loc_wt, by="date")


############################### read in COVID-Data######################################
# access here: https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data/csse_covid_19_time_series


library(httr)
library(tidyr)
time_series_covid19_confirmed_US <-read.csv(("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv"),header=T)


head(time_series_covid19_confirmed_US)
time_series_covid19_confirmed_US=time_series_covid19_confirmed_US %>% gather(date, cases, `X1.22.20`:`X4.30.20`)

time_series_covid19_confirmed_US$date<-substring(time_series_covid19_confirmed_US$date, 2)
time_series_covid19_confirmed_US$date= as.Date(time_series_covid19_confirmed_US$date, format= "%m.%d.%y")
time_series_covid19_confirmed_US$DOW<-as.factor(weekdays(time_series_covid19_confirmed_US$date))


#time_series_covid19_deaths_US <- read_csv("C:/Users/Richard Remigio/Desktop/covid19/time_series_covid19_deaths_US.csv")


time_series_covid19_deaths_US <-read.csv(("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_US.csv"),header=T)


head(time_series_covid19_deaths_US )
time_series_covid19_deaths_US =time_series_covid19_deaths_US  %>% gather(date, cases, `X1.22.20`:`X4.30.20`)

time_series_covid19_deaths_US $date<-substring(time_series_covid19_deaths_US $date, 2)
time_series_covid19_deaths_US $date= as.Date(time_series_covid19_deaths_US $date, format= "%m.%d.%y")
time_series_covid19_deaths_US $DOW<-as.factor(weekdays(time_series_covid19_deaths_US $date))

############################loc##################################
#loc
## use FIPS code to specify a location
## FIPS look up: https://www.nrcs.usda.gov/wps/portal/nrcs/detail/national/home/?cid=nrcs143_013697
## Albany= 36001

## clean up time series data
loc_cases=time_series_covid19_confirmed_US %>% filter(FIPS=="36001") %>% select(date, UID, Admin2, Province_State, FIPS, cases, iso2, iso3, code3, Country_Region, Lat, Long_, Combined_Key, DOW)
loc_deaths=time_series_covid19_deaths_US %>% filter(FIPS=="36001") %>% select(date, death_cases=cases)
#loc_2121903 <- read_csv("C:/Users/Richard Remigio/Desktop/covid19/weather/loc_2121903.csv", 
#   col_types = cols(DATE = col_date(format = "%Y-%m-%d")))

##decompose cumulative counts
loc_cases$cases_cum<-loc_cases$cases
loc_cases= loc_cases %>% mutate(cases=cases_cum- lag(cases_cum,1 ))

loc_deaths$death_cases_cum<-loc_deaths$death_cases
loc_deaths= loc_deaths %>% mutate(death_cases=death_cases_cum- lag(death_cases_cum,1 ))


## join cases and death data frames
loc_df= left_join(loc_cases, loc_deaths, by="date")
loc_weather=loc_weather %>% transmute(id=id, DATE=date, TMIN=tmin, PRCP=prcp, WSF2=wsf2, WSF5=wsf5, WT01=wt01, WT02=wt02, WT03=wt03, WT08=wt08)

## join weather and cases
loc_df=left_join(loc_weather,loc_df, by=c("DATE"="date")) %>% 
  mutate(cases_5dayroll=rollapply(cases,5,mean,align='center',fill=NA), death_cases_5dayroll=rollapply(death_cases,5,mean,align='center',fill=NA),
         bad_weather=ifelse(!is.na(WT01),"1", ifelse(!is.na(WT02), "1", ifelse(!is.na(WT03), "1",  ifelse(!is.na(WT08), "1", "0")))))

loc_df$bad_weather.c= as.numeric(loc_df$bad_weather)




## loc analysis
######################################################
setwd("~/Folders")
#############################################################3

summary(glm(cases_5dayroll~DATE + TMIN + DOW, family=quasipoisson, data=loc_df))

varknots <- equalknots(loc_df$TMIN, nk=3)
lagknots <- logknots(fun="ns", x= 30, nk=4)

cb_tmin_lin <- crossbasis(loc_df$TMIN, lag=c(0,30), argvar=list(fun="lin"),  arglag=list(knots=lagknots))

mod=glm(cases_5dayroll~ cb_tmin_lin + DATE + DOW, family=quasipoisson, data=loc_df)
cp_tmin<- crosspred(cb_tmin_lin, mod, by=0.1)
quantile(loc_df$TMIN, probs = 0.1, na.rm = T)


png("loc/loc_temp_3D.png")
#plot(cp_tmin,  main="Min Temp lag-response functions")
mar <- par()$mar 
par(mar=c(2.5,1,0.5,1))
d3 <- plot(cp_tmin,xlab="Temperature (C)",zlab="RR", phi=35, ltheta=170,  main="Location")
lines(trans3d(x=1.7, y=0:30, z=cp_tmin$matRRfit[as.character(1.7),], pmat=d3), col=2, lwd=2)
lines(trans3d(x=cp_tmin$predvar, y=15, z=cp_tmin$matRRfit[,"lag15"], pmat=d3), col=3, lwd=2)
par(mar=mar)
dev.off()

## specify TMIN in "var" cariable
png("loc/loc_temp_slices.png")
plot(cp_tmin, "slices",  var=(-3), col=3, ylab="RR", main="Location")
dev.off()

## specify lags in "lag" cariable

png("loc/loc_temp_multi_slices.png")
plot(cp_tmin, "slices",  lag=c(0, 1,2), col=3, ylab="RR", main="Temperature, C\n 10%ile TMIN\nLocation")
dev.off()


#### loc death#####
summary(glm(death_cases~DATE + TMIN + DOW, family=quasipoisson, data=loc_df))

mod=glm(death_cases~ cb_tmin_lin + DATE + DOW, family=quasipoisson, data=loc_df)
cp_tmin<- crosspred(cb_tmin_lin, mod,  cen=median(loc_df$TMIN, na.rm=T), by=0.5, from = c(min(loc_df$TMIN, na.rm=T), max(loc_df$TMIN, na.rm=T)))
quantile(loc_df$TMIN, probs = 0.5, na.rm = T)

png("loc/loc_temp_3D_death.png")
#plot(cp_tmin,  main="Min Temp lag-response functions")
mar <- par()$mar 
par(mar=c(2.5,1,0.5,1))
d3 <- plot(cp_tmin,xlab="Temperature (C)",zlab="RR", main="loc", shade=0.9)
#lines(trans3d(x=0, y=0:30, z=cp_tmin$matRRfit[as.character(5),], pmat=d3), col=2, lwd=2)
#lines(trans3d(x=cp_tmin$predvar, y=15, z=cp_tmin$matRRfit[,"lag15"], pmat=d3), col=3, lwd=2)
par(mar=mar)
dev.off()

png("loc/loc_temp_slices_death.png")
## specify TMIN in "var" variable
plot(cp_tmin, "slices",  var=c(0), col=3, ylab="RR", main= "loc_death")
dev.off()

## specify lags in "lag" variable

png("loc/loc_temp_multi_slices_lag_death.png")
plot(cp_tmin, "slices",  lag=c(0, 1, 2), col=3, ylab="RR",  main="Temperature, C\n 10%ile TMIN\nloc_death")
dev.off()
