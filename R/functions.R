# functions


#' get data
#'
#' @param file The data file (it needs to be updated annually)
#' @importFrom rio import_list
#' @return
#' @export
#'

get_data<-function(file){
  data_list <-withr::with_locale(c(LC_CTYPE="Norwegian Bokmål.1252"), rio::import_list(file))
}


#' make_correction 
#'
#' @param lokalitetsdata 
#' @param transektdata 
#' @param rutedata 
#' @param artsdata 
#' @param korreksjonsdata 
#'
#' @return
#' @export
#'
make_correction<-function(lokalitetsdata,
                          transektdata,
                          rutedata,
                          artsdata,
                          korreksjonsdata){
  rutedata <- rutedata[!is.na(rutedata$Lokalitet), ]
  
  ### Datajusteringer
  
  lokalitetsdata <- lokalitetsdata[!apply(lokalitetsdata,1,function(x){all(is.na(x))}),]
  transektdata <- transektdata[!apply(transektdata,1,function(x){all(is.na(x))}),]
  rutedata <- rutedata[!apply(rutedata,1,function(x){all(is.na(x))}),]
  
  
  
  lokalitetsdata$Lokalitet <- factor(lokalitetsdata$Lokalitet)
  
  lok.utelates <- c("Nedre Røykenvik","Øvre Røykenvik")
  lokalitetsdata <- lokalitetsdata[!(lokalitetsdata$Lokalitet%in%lok.utelates),]; lokalitetsdata$Lokalitet <- factor(lokalitetsdata$Lokalitet)
  transektdata <- transektdata[!(transektdata$Lokalitet%in%lok.utelates),]; transektdata$Lokalitet <- factor(transektdata$Lokalitet)
  rutedata <- rutedata[!(rutedata$Lokalitet%in%lok.utelates),]; rutedata$Lokalitet <- factor(rutedata$Lokalitet)
  
  lokalitetsdata <- lokalitetsdata[!(lokalitetsdata$Lokalitet%in%c("Falang") & lokalitetsdata$År<2020),]
  transektdata <- transektdata[!(transektdata$Lokalitet%in%c("Falang") & transektdata$År<2020),]
  rutedata <- rutedata[!(rutedata$Lokalitet%in%c("Falang") & rutedata$År<2020),]
  
  lokalitetsdata$Hovednaturtype <- sapply(strsplit(lokalitetsdata$Naturtype,"-"),function(x){x[[1]]})
  
  lokalitetsdata$Hovedskjøtsel <- as.character(lokalitetsdata$Krattrydding_binær)
  lokalitetsdata$Hovedskjøtsel[lokalitetsdata$Krattrydding_binær=="Ja"] <- "Krattrydding"
  lokalitetsdata$Hovedskjøtsel[lokalitetsdata$Slått_binær=="Ja"] <- "Slått"
  lokalitetsdata$Hovedskjøtsel[lokalitetsdata$Hovedskjøtsel=="Nei"] <- "Ingen/annet"
  lokalitetsdata$Hovedskjøtsel <- factor(lokalitetsdata$Hovedskjøtsel)
  
  transektdata$`Forekomst dragehode (m)`[transektdata$`Forekomst dragehode (m)`%in%c("ingen"," ")] <- NA
  
  rutedata$Vedpl.dekn <- beregn.dekning(rutedata,artsdata,"Vedplante")
  
  lokalitetsnavn <- levels(lokalitetsdata$Lokalitet)
  regionnavn <- levels(factor(lokalitetsdata$Region))
  naturtypenavn <- levels(factor(lokalitetsdata$Hovednaturtype))
  
  
  # # Horgen og Møllerenga: tellinger av skudd (vegetative og fertile) korrigeres til estimert antall individer
  names(korreksjonsdata)[6:11] <- c("x1","x2","x3","y1","y2","y3")
  m1 <- lm(y1~-1+x1,data=korreksjonsdata)
  m2 <- lm(y2~-1+x2,data=korreksjonsdata)
  m3 <- lm(y3~-1+x3,data=korreksjonsdata)

  i <- rutedata$Lokalitet%in%c("Horgen","Møllerenga") & rutedata$År%in%c(2019,2020)
  newdat <- rutedata[i,c("Veg.planter","Fert.planter")]
  names(newdat) <- c("x2","x3")
  newdat$x2<-as.numeric(newdat$x2)
  newdat$x3<-as.numeric(newdat$x3)
  y2 <- round(predict(m2,newdat))
  y3 <- round(predict(m3,newdat))
  y4 <- as.numeric(rutedata[i,"Småplanter"])+y2+y3
  
  rutedata[i,"Veg.planter"] <- y2
  rutedata[i,"Fert.planter"] <- y3
  rutedata[i,"Ant.DR"] <- y4
  
  rutedata$Småplanter<-as.numeric(rutedata$Småplanter)
  rutedata$Veg.planter<-as.numeric(rutedata$Veg.planter)
  rutedata$Fert.planter<-as.numeric(rutedata$Fert.planter)
  rutedata$Ant.DR<-as.numeric(rutedata$Ant.DR)
  rutedata
}





#' Calculate coverage
#'
#' @param rutedata 
#' @param artsdata 
#' @param artsgruppe 
#'
#' @return
#' @export
#' 

beregn.dekning <- function(rutedata,artsdata,artsgruppe)
{
  #withr::with_locale(c(LC_CTYPE="Norwegian Bokmål.1252"), {
  dekning <- rep(0,nrow(rutedata))
  for(i in 1:nrow(rutedata))
  {
    utvalgte.data <- artsdata[,artsgruppe]==1 & artsdata$Rute==rutedata$RuteID[i] & artsdata$Ãr==rutedata$Ãr[i]
    if(length(utvalgte.data)>0) dekning[i] <- sum(artsdata$Dekning[utvalgte.data],na.rm=T)
  }
  return(dekning)
#})
}

beregn.arealvekt <- function(x,design,lokalitetsbredde=NA,rutebredde=1)
{
  if(design=="Linje")
  {
    return(rep(lokalitetsbredde/rutebredde,length(x)))
  }
  if(design=="Rose")
  {
    r1 <- x             # radius for indre sirkel (rutas avstand fra startpunktet)
    r2 <- x+rutebredde  # radius for ytre sirkel
    return( (pi*r2^2 - pi*r1^2)/8 )
  }
}

veid.sum <- function(x,v) sum(x*v)
 
boot.pop <- function(tetthet,forekomstareal,nboot=2000)
 {
   #print(tetthet)
   tetthet <- tetthet[!is.na(tetthet)]
   #print(tetthet)
   #print(forekomstareal)
   if(length(tetthet)==0) {print("Ingen individer"); return(rep(NA,nboot))}
   tetthet.boot <- boot(tetthet,statistic=function(x,k){mean(x[k])},R=nboot)
   forekomstareal.boot <- boot(forekomstareal,statistic=function(x,k){sum(x[k])},R=nboot)
   tetthet.boot$t * forekomstareal.boot$t
 }
 
beregn.forekomstareal <- function(transektdata,utvalgt.lokalitet,year,design,lokalitetsbredde,rutebredde=1,
                                   forekomst_transekt="Dragehode",forekomst_avstand="Forekomst dragehode (m)")
 {
   x <- list()
   for(i in 1:length(year))
   {
     transekter <- transektdata[transektdata$Lokalitet==utvalgt.lokalitet & transektdata$År==year[i],]
     #print(year[i])
     #print(transekter)
     if(all(is.na(transekter[,forekomst_transekt])))
     {
       #print("Forekomstruter ikke registrert")
       next
     }
     if(all(transekter[,forekomst_transekt]==0))
     {
       #print("Forekomstruter ikke registrert")
       next
     }
     forekomstruter <- as.numeric(unlist(sapply(transekter[,forekomst_avstand],strsplit,split=',')))
     #print(forekomstruter)
     forekomstruter <- forekomstruter[!is.na(forekomstruter)]
     #print(forekomstruter)
     forekomstareal <- beregn.arealvekt(forekomstruter,design[i],lokalitetsbredde[i],rutebredde=rutebredde)
     #print(forekomstareal)
     x[[i]] <- forekomstareal
   }
   mangler <- unlist(lapply(x,is.null))
   if(any(mangler))
   {
     erstattes <- which(mangler)
     registrert <- which(!mangler)
     #print(erstattes)
     #print(registrert)
     for(i in erstattes)
     {
       delta <- abs((registrert+0.5) - i)
       #print(delta)
       j <- registrert[delta==min(delta)]
       #print(j)
       x[[i]] <- x[[j]]
     }
   }
   return(x)
 }
# 
# 
 beregn.popstruktur <- function(utvalgt.lokalitet,lokalitetsdata,transektdata,rutedata,rutebredde=1,
                                fert="Fert.planter",veg="Veg.planter",sma="Småplanter",tot="Ant.DR",
                                forekomst_transekt="Dragehode",forekomst_avstand="Forekomst dragehode (m)",
                                quantiles=c(0.025,0.975))
 {
   year <- lokalitetsdata[lokalitetsdata$Lokalitet==utvalgt.lokalitet,]$År
   design <- lokalitetsdata[lokalitetsdata$Lokalitet==utvalgt.lokalitet,]$Design
   lokalitetsbredde <- lokalitetsdata[lokalitetsdata$Lokalitet==utvalgt.lokalitet,]$Lokalitetsbredde
#   print(paste(utvalgt.lokalitet, year, design,"bredde:",lokalitetsbredde))
   nboot <- 2000
   nq <- length(quantiles)
   nyear <- length(year)
   nFert <- numeric(nyear)
   nVeg <- numeric(nyear)
   nSma <- numeric(nyear)
   nTot <- numeric(nyear)
   Fert.CI <- matrix(NA,nyear,nq)
   Veg.CI <- matrix(NA,nyear,nq)
   Sma.CI <- matrix(NA,nyear,nq)
   Tot.CI <- matrix(NA,nyear,nq)
 
forekomstareal.liste <- beregn.forekomstareal(transektdata,utvalgt.lokalitet,year,design,lokalitetsbredde,rutebredde,forekomst_transekt,forekomst_avstand)
   print(forekomstareal.liste)
 
   for(i in 1:length(year))
   {
     forekomstareal <- forekomstareal.liste[[i]]
     saveRDS(forekomstareal, paste0("data/derived_data/", utvalgt.lokalitet,"_forekomstareal","_",i,".RDS"))
     #print(year[i])
    ruter <- rutedata[rutedata$Lokalitet==utvalgt.lokalitet & rutedata$År==year[i] & rutedata[,tot]>0,]
    saveRDS(ruter, paste0("data/derived_data/", utvalgt.lokalitet,"_ruter","_",i ,".RDS"))

    #print("Fertile")
    if(!is.na(fert)) x <- boot.pop(ruter[,fert],forekomstareal,nboot)
    else x <- NA
    #print(summary(x))
    nFert[i] <- mean(x); Fert.CI[i,] <- quantile(x,quantiles,na.rm=T)

    #print("Vegetative")
    if(!is.na(veg)) x <- boot.pop(ruter[,veg],forekomstareal,nboot)
    else x <- NA
    #print(summary(x))
    nVeg[i] <- mean(x); Veg.CI[i,] <- quantile(x,quantiles,na.rm=T)

    #print("SmÃ¥planter")
    if(!is.na(sma)) x <- boot.pop(ruter[,sma],forekomstareal,nboot)
    else x <- NA
    #print(summary(x))
    nSma[i] <- mean(x); Sma.CI[i,] <- quantile(x,quantiles,na.rm=T)

    #print("Totalt")
    if(!is.na(tot)) x <- boot.pop(ruter[,tot],forekomstareal,nboot)
    else x <- NA
    #print(summary(x))
    nTot[i] <- mean(x); Tot.CI[i,] <- quantile(x,quantiles,na.rm=T)
  }

  popstr <- list(lokalitet=utvalgt.lokalitet,year=year,nFert=nFert,nVeg=nVeg,nSma=nSma,nTot=nTot,Fert.CI=Fert.CI,Veg.CI=Veg.CI,Sma.CI=Sma.CI,Tot.CI=Tot.CI)
  return(popstr)
}

plot.popstruktur <- function(popstr,tidsrom=c(2016.5,2021.5))
{
  plot(popstr$year,popstr$nTot,ylim=c(0,max(popstr$Tot.CI,na.rm=T)),type="o",main=popstr$lokalitet,xlab="Ãr",ylab="Antall individer",xlim=tidsrom)
  polygon(c(popstr$year,rev(popstr$year)),c(popstr$Tot.CI[,1],rev(popstr$Tot.CI[,2])),col=rgb(0.5,0.5,0.5,alpha=0.2),border=NA)
  lines(popstr$year,popstr$nFert,col="red",type="o"); polygon(c(popstr$year,rev(popstr$year)),c(popstr$Fert.CI[,1],rev(popstr$Fert.CI[,2])),col=rgb(1,0,0,alpha=0.2),border=NA)
  lines(popstr$year,popstr$nVeg,col="green",type="o"); polygon(c(popstr$year,rev(popstr$year)),c(popstr$Veg.CI[,1],rev(popstr$Veg.CI[,2])),col=rgb(0,1,0,alpha=0.2),border=NA)
  lines(popstr$year,popstr$nSma,col="blue",,type="o"); polygon(c(popstr$year,rev(popstr$year)),c(popstr$Sma.CI[,1],rev(popstr$Sma.CI[,2])),col=rgb(0,0,1,alpha=0.2),border=NA)
  #lines(popstr$year,popstr$nSma+popstr$nVeg+popstr$nFert,lty=2,type="o")
  legend("topleft",c("Totalt","Fertile","Vegetative","Småplanter"),lty=c(1,1,1,1),pch=c(1,1,1,1),col=c("black","red","green","blue"),bty="n")
}

plot.gruppetrender <- function(lokalitetsdata,gruppevariabel,lokalitetsestimater,tidsrom=c(2016.5,2021.5),reverser=F,lokalitetsnavn,lokfarge,regsymbol,natlinje)
{
  gruppe <- unique(lokalitetsdata[,gruppevariabel])
  if(reverser) gruppe <- rev(gruppe)
  ngrp <- length(gruppe)
  par(mfcol=c(4,ngrp),mar=c(2,4,2,0.2))
  for(i in 1:ngrp)
  {
    grp <- as.character(gruppe[i])
    print(grp)
    lok <- as.character(unique(lokalitetsdata[lokalitetsdata[,gruppevariabel]==grp,"Lokalitet"]))
    farge <- lokfarge[lok]
    #print(lok)
    #print(farge

    variabler <- c("nTot","nFert","nVeg","nSma")
    if(i==1) ylab <- c("Antall totalt","Antall fertile","Antall vegetative","Antall småplanter")
    else ylab <- rep("",length(variabler))
    xlab <- c(rep("",length(variabler)-1),"Ã…r")
    main <- c(grp,rep("",length(variabler)-1))
    for(i in 1:length(variabler))
    {
      noplot <- TRUE
      for(j in 1:length(lok))
      {
        #print(lok[j])
        reg <- as.character(unique(lokalitetsdata[lokalitetsdata$Lokalitet==lok,"Region"]))
        symbol <- regsymbol[reg]
        nat <- as.character(unique(lokalitetsdata[lokalitetsdata$Lokalitet==lok,"Hovednaturtype"]))
        linje <- natlinje[nat]
        #print(linje)

        popstr <- lokalitetsestimater[lok[j]][[1]]
        if(length(popstr)==1) {print("Ingen estimater"); next}
        if(noplot)
        {
          plot(popstr$year,popstr[variabler[i]][[1]],xlim=tidsrom,ylim=c(1,70000),type="o",col=farge[j],     #pch=symbol,lty=linje,
               main=main[i],xlab=xlab[i],ylab=ylab[i],log="y")
          noplot <- FALSE
        }
        else lines(popstr$year,popstr[variabler[i]][[1]],type="o",col=farge[j]) #,pch=symbol,lty=linje)
      }
    }
  }
}


beregn.vekstrate <- function(rutedata,popvar,idvar)
{
  vekstrate <- rep(NA,nrow(rutedata))
  year <- sort(as.numeric(unique(rutedata$Ãr)))
  nyear <- length(year)
  if(nyear<2) return(NA)
  for(i in 1:(nyear-1))
  {
    t1 <- rutedata$Ãr==year[i]
    t2 <- rutedata$Ãr==year[i+1]
    idmatch <- match(rutedata[t1,idvar],rutedata[t2,idvar])
    vekstrate[t1] <- rutedata[t2,popvar][idmatch] / rutedata[t1,popvar]
  }
  rutedata$vekstrate <- vekstrate
  return(rutedata)
}


#######################
# Populasjonsstørrelser
#######################

Populasjon<-function(rutedata,lokalitetsdata, transektdata){
  lokalitetsnavn <- levels(lokalitetsdata$Lokalitet)
  
  ### Beregninger
  
  # # Sjekk om forekomst varierer med avstand fra midtpunkt
  # par(mfrow=c(1,1))
  # j <- rutedata$Lokalitet!="Ekebergskråningen"
  # plot(rutedata$Avst[j],rutedata$Ant.DR[j]>0,main="Alle lokaliteter")
  # par(mfrow=c(5,5))
  # for(i in 1:length(lokalitetsnavn))
  # {
  #   j <- rutedata$Lokalitet==lokalitetsnavn[i]
  #   plot(rutedata$Avst[j],rutedata$Ant.DR[j]>0,main=lokalitetsnavn[i])
  # }
  # 
  # # Sjekk om tetthet varierer med avstand fra midtpunkt
  # par(mfrow=c(1,1))
  # j <- rutedata$Ant.DR>0 & rutedata$Lokalitet!="Ekebergskråningen"
  # plot(rutedata$Avst[j],rutedata$Ant.DR[j],main="Alle lokaliteter")
  # par(mfrow=c(5,5))
  # for(i in 1:length(lokalitetsnavn))
  # {
  #   j <- rutedata$Lokalitet==lokalitetsnavn[i] & rutedata$Ant.DR>0
  #   plot(rutedata$Avst[j],rutedata$Ant.DR[j],main=lokalitetsnavn[i])
  # }
  # 
  
  # For hver populasjon (fordi de har særegenheter i design, og det er ønskelig å presentere resultater for enkelt-populasjoner på nett)
  # Også for hvert år (Hensiktsmessig hvis f.eks. design endrer seg...
  
  lokalitetsestimater <- list()
  for(i in 1:length(lokalitetsnavn))
  {
    tryCatch({
      print(lokalitetsnavn[i])
      popstr <- beregn.popstruktur(lokalitetsnavn[i],lokalitetsdata,transektdata,rutedata,quantiles=c(0.25,0.75))
      lokalitetsestimater[[i]] <- popstr
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  names(lokalitetsestimater) <- lokalitetsnavn
  
  
  ### Figurer
  
  if(!dir.exists("Figurer")) dir.create("Figurer")
  
  # "Rå" figurer av alle lokaliteter med populasjonsestimater m/CI over tid for fertile, vegetative, små og totalt
  #par(mfrow=c(5,5))
  pdf(paste("Figurer/","Alle lokaliteter",".pdf",sep=""))
  for(i in 1:length(lokalitetsestimater))
  {
    popstr <- lokalitetsestimater[[i]]
    plot.popstruktur(popstr)
  }
  dev.off()
  
  for(i in 1:length(lokalitetsestimater))
  {
    #pdf(paste("Figurer/",lokalitetsnavn[i],".pdf",sep=""))
    #png(paste("Figurer/",lokalitetsnavn[i],".png",sep=""))
    jpeg(paste("Figurer/",lokalitetsnavn[i],".jpg",sep=""))
    popstr <- lokalitetsestimater[[i]]
    plot.popstruktur(popstr)
    dev.off()
  }
  
  
  # Populasjonstrender per gruppe, totalt, fertile, vegetative og små
  
  lokfarge <- rainbow(length(lokalitetsnavn)); names(lokfarge) <- lokalitetsnavn
  regsymbol <- c(1,2,3); names(regsymbol) <- regionnavn
  natlinje <- c(1,2); names(natlinje) <- naturtypenavn
  
  # Regioner
  pdf(paste("Figurer/Regioner/","Alle regioner",".pdf",sep=""))
  plot.gruppetrender(lokalitetsdata,gruppevariabel="Region",lokalitetsestimater,tidsrom=c(2016.5,2021.5),reverser=T,lokalitetsnavn,lokfarge,regsymbol,natlinje)
  par(mfrow=c(1,3))
  plot(0,0,type="n",axes=F,xlab="",ylab=""); legend("topleft",lty=1,col=lokfarge,legend=lokalitetsnavn,bty="n")
  dev.off()
  
  # Naturtyper
  pdf(paste("Figurer/Naturtyper/","Alle hovedtyper",".pdf",sep=""))
  plot.gruppetrender(lokalitetsdata,gruppevariabel="Hovednaturtype",lokalitetsestimater,tidsrom=c(2016.5,2021.5),reverser=T,lokalitetsnavn,lokfarge,regsymbol,natlinje)
  par(mfrow=c(1,3))
  plot(0,0,type="n",axes=F,xlab="",ylab=""); legend("topleft",lty=1,col=lokfarge,legend=lokalitetsnavn,bty="n")
  dev.off()
  
  
  pdf(paste("Figurer/Skjøtsel/","Hovedskjøtsel",".pdf",sep=""))
  plot.gruppetrender(lokalitetsdata,gruppevariabel="Hovedskjøtsel",lokalitetsestimater,tidsrom=c(2016.5,2021.5),reverser=T,lokalitetsnavn,lokfarge,regsymbol,natlinje)
  par(mfrow=c(1,3))
  plot(0,0,type="n",axes=F,xlab="",ylab=""); legend("topleft",lty=1,col=lokfarge,legend=lokalitetsnavn,bty="n")
  dev.off()
  
  
  # Populasjonsestimater 2017-2020
  utvalgte.lok <- lokalitetsdata$Lokalitet[lokalitetsdata$År==2017]
  for(i in 1:length(utvalgte.lok))
  {
    lok <- utvalgte.lok[i]
    print(lokalitetsestimater[lok])
  }
  
}

