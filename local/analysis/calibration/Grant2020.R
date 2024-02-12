######### Load R workspace ######### 
setwd("~/Code/Calibration")
load("Grant2020_PAGES_workspace.RData") ### Has core datasets listed 

##PlioSeaNZ<-read.csv("PlioSeaNZ.csv")
##Miller2020=read.csv("Miller2020.csv")
##Dwyer=read.csv('DwyerChandler2009.csv')
##Sosdian=read.csv('SosdianRosenthal2009.csv')
##Rohling=read.csv("Rohlingetal2014.csv")
## LR04=read.csv("LR04.csv")
## Naish<-read.csv("Miller2012_Naishcal.csv",header=TRUE,sep=",")


library(astrochron)
library(ggplot2)
library(dplyr)

rsl <- PlioSeaNZ$norm.mean
prior <- density(rsl)
plot(prior)

t20 <- na.omit(PlioSeaNZ_amp$`t=20`)
t40 <- na.omit(PlioSeaNZ_amp$`t=40`)
t100 <-na.omit(PlioSeaNZ_amp$`t=100`)

plot(density(t20))
plot(density(t40))
plot(density(t100))

quantile(t20)

###### PlioSeaNZ Grant et al 2019 RSL amplitude assessment ######

iso=iso(PlioSeaNZ,3032,3310) # isolate record to mPWP
strats(iso) # a strat assessment to determiene sampling variability
lin_PlioSeaNZ=linterp(iso,2) # linear interpolation at mean or median sampling

  site<-lin_PlioSeaNZ
  outmain<-lin_PlioSeaNZ[1]
  n=2 #dt sampling
  
  #### This calculates the minima and maxima at various window durations (Milankovitch periodicities)

for(i in 1:nrow(site)) {
  #period = 20 
  celln <- (20/n) 
  
  low <- min(site[i:(i+celln),2])  # lowest and highest levels for target section/group combination
  high <- max(site[i:(i+celln),2])
  outmain[i,2]<-high-low 
  
}

  for(i in 1:nrow(site)) {
  #period = 41 
  celln <- (41/n) 
  
  low <- min(site[i:(i+celln),2])  # lowest and highest levels for target section/group combination
  high <- max(site[i:(i+celln),2])
  outmain[i,3]<-high-low 
  
  }

  for(i in 1:nrow(site)) {
  #period = 100 
  celln <- (100/n) 
  
  low <- min(site[i:(i+celln),2])  # lowest and highest levels for target section/group combination
  high <- max(site[i:(i+celln),2])
  outmain[i,4]<-high-low 
  
  }
  
colnames(outmain)<-c("time.ka","t=20","t=40","t=100")
PlioSeaNZ_amp<-outmain



########### Miller et al 2020 amplitude assesment ##################

Milleriso=iso(Miller2020,3032,3330)
strats(Milleriso)
lin=linterp(Milleriso,2.5)

site<-lin
outmain<-lin[1]
n=2.5
for(i in 1:nrow(site)) {
  #period = 20 
  celln <- (20/n) 
  
  low <- min(site[i:(i+celln),2])  # lowest and highest levels for target section/group combination
  high <- max(site[i:(i+celln),2])
  outmain[i,2]<-high-low 
  
}


for(i in 1:nrow(site)) {
  #period = 41 
  celln <- (41/n)
  
  low <- min(site[i:(i+celln),2])  # lowest and highest levels for target section/group combination
  high <- max(site[i:(i+celln),2])
  outmain[i,3]<-high-low 
  
}

for(i in 1:nrow(site)) {
  #period = 100 
  celln <- (100/n) 
  
  low <- min(site[i:(i+celln),2])  # lowest and highest levels for target section/group combination
  high <- max(site[i:(i+celln),2])
  outmain[i,4]<-high-low 
  
}
colnames(outmain)<-c("time.ka","t=20","t=40","t=100")
Miller_amp<-outmain


########### Dwyer and Chandler 2009 amplitude assesment ##################

Dwyeriso=iso(Dwyer,3032,3330)
strats(Dwyeriso)
lin=linterp(Dwyeriso,4)

site<-lin
outmain<-lin[1]
n=4

for(i in 1:nrow(site)) {
  #period = 20 
  celln <- (20/n) 
  
  low <- min(site[i:(i+celln),2])  # lowest and highest levels for target section/group combination
  high <- max(site[i:(i+celln),2])
  outmain[i,2]<-high-low 
  
}


for(i in 1:nrow(site)) {
  #period = 41 
  celln <- (41/n) 
  
  low <- min(site[i:(i+celln),2])  # lowest and highest levels for target section/group combination
  high <- max(site[i:(i+celln),2])
  outmain[i,3]<-high-low 
  
}

for(i in 1:nrow(site)) {
  #period = 100 
  celln <- (100/n) 
  
  low <- min(site[i:(i+celln),2])  # lowest and highest levels for target section/group combination
  high <- max(site[i:(i+celln),2])
  outmain[i,4]<-high-low 
  
}
colnames(outmain)<-c("time.ka","t=20","t=40","t=100")
Dwyer_amp<-outmain




########### Sosdian and Rosenthal 2009 amplitude assesment ##################

Sosdianiso=iso(Sosdian,3032,3330)
strats(Sosdianiso)
lin=linterp(Sosdianiso,5.5)

site<-lin
outmain<-lin[1]
n=5.5

for(i in 1:nrow(site)) {
  #period = 20 
  celln <- (20/n)  
  
  low <- min(site[i:(i+celln),2])  # lowest and highest levels for target section/group combination
  high <- max(site[i:(i+celln),2])
  outmain[i,2]<-high-low 
  
}


for(i in 1:nrow(site)) {
  #period = 41 
  celln <- (41/n) 
  
  low <- min(site[i:(i+celln),2])  # lowest and highest levels for target section/group combination
  high <- max(site[i:(i+celln),2])
  outmain[i,3]<-high-low 
  
}

for(i in 1:nrow(site)) {
  #period = 100 
  celln <- (100/n) 
  
  low <- min(site[i:(i+celln),2])  # lowest and highest levels for target section/group combination
  high <- max(site[i:(i+celln),2])
  outmain[i,4]<-high-low 
  
}
colnames(outmain)<-c("time.ka","t=20","t=40","t=100")
Sosdian_amp<-outmain


########### Rohling et al., 2014 amplitude assesment ##################

Rohlingiso=iso(Rohling,3032,3330)
strats(Rohlingiso)

lin=linterp(Rohlingiso,1)

site<-lin
outmain<-lin[1]
n=1

for(i in 1:nrow(site)) {
  #period = 20 
  celln <- (20/n) 
  
  low <- min(site[i:(i+celln),2])  # lowest and highest levels for target section/group combination
  high <- max(site[i:(i+celln),2])
  outmain[i,2]<-high-low 
  
}


for(i in 1:nrow(site)) {
  #period = 41 
  celln <- (41/n)
  
  low <- min(site[i:(i+celln),2])  # lowest and highest levels for target section/group combination
  high <- max(site[i:(i+celln),2])
  outmain[i,3]<-high-low 
  
}

for(i in 1:nrow(site)) {
  #period = 100 
  celln <- (100/n)
  
  low <- min(site[i:(i+celln),2])  # lowest and highest levels for target section/group combination
  high <- max(site[i:(i+celln),2])
  outmain[i,4]<-high-low 
  
}
colnames(outmain)<-c("time.ka","t=20","t=40","t=100")
Rohling_amp<-outmain


################# Miller et al 2012 LR04 calibration #####


LR04iso=iso(LR04,3032,3330)
Miller2012=LR04iso
Miller2012[,2]=(((LR04iso[,2]-3.23)*0.66)/0.01)+50


site<-Miller2012
outmain<-Miller2012[1]
n=5

for(i in 1:nrow(site)) {
  #period = 20 
  celln <- (20/n) 
  
  low <- min(site[i:(i+celln),2])  # lowest and highest levels for target section/group combination
  high <- max(site[i:(i+celln),2])
  outmain[i,2]<-high-low 
  
}


for(i in 1:nrow(site)) {
  #period = 41 
  celln <- (41/n) 
  
  low <- min(site[i:(i+celln),2])  # lowest and highest levels for target section/group combination
  high <- max(site[i:(i+celln),2])
  outmain[i,3]<-high-low 
  
}

for(i in 1:nrow(site)) {
  #period = 100 
  celln <- (100/n) 
  
  low <- min(site[i:(i+celln),2])  # lowest and highest levels for target section/group combination
  high <- max(site[i:(i+celln),2])
  outmain[i,4]<-high-low 
  
}
colnames(outmain)<-c("time.ka","t=20","t=40","t=100")
Miller2012_amp<-outmain

########## Miller et al., 2012 : Naish NZ backstrip LR04 scale ############

iso=iso(Naish,3032,3330)

site<-iso
outmain<-iso[1]
n=5

for(i in 1:nrow(site)) {
  #period = 20 
  celln <- (20/n) 
  
  low <- min(site[i:(i+celln),2])  # lowest and highest levels for target section/group combination
  high <- max(site[i:(i+celln),2])
  outmain[i,2]<-high-low 
  
}


for(i in 1:nrow(site)) {
  #period = 41 
  celln <- (41/n) 
  
  low <- min(site[i:(i+celln),2])  # lowest and highest levels for target section/group combination
  high <- max(site[i:(i+celln),2])
  outmain[i,3]<-high-low 
  
}

for(i in 1:nrow(site)) {
  #period = 100 
  celln <- (100/n) 
  
  low <- min(site[i:(i+celln),2])  # lowest and highest levels for target section/group combination
  high <- max(site[i:(i+celln),2])
  outmain[i,4]<-high-low 
  
}
colnames(outmain)<-c("time.ka","t=20","t=40","t=100")
Naish2012_amp<-outmain

############## Figure 1a ############

Amplitude<-list(Dwyer_amp[,2],Sosdian_amp[,2],Miller_amp[,2],Miller2012_amp[,2],Rohling_amp[,2],Naish2012_amp[,2],PlioSeaNZ_amp[,2])
labs<-c("Dwyr.Chndlr2009","Sosdn.Rsnthl2009","MillerEtAl2020","MillerEtAl2012","Rohling2014","Naish2012","PlioSeaNZ")

cols3<-c("hotpink4","hotpink3","firebrick","firebrick2","lightsteelblue4","dodgerblue3","deepskyblue1")

Amplitude.df<-data.frame(RSL.m = unlist(Amplitude), Study = rep(labs, lengths(Amplitude)))
Amplitude.df[,3]<-rep(cols3,lengths(Amplitude))
Ampl.df<-na.omit(Amplitude.df)


f <- function(x) {
  r <- quantile(x, probs = c(0.1, 0.159, 0.5, 0.841, 0.9))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

m <- function(x) {
  return(quantile(x,probs=0.5))
}

n=lengths(Amplitude)

Ampl.df$Study<-factor(Ampl.df$Study,levels=c("Dwyr.Chndlr2009","Sosdn.Rsnthl2009","MillerEtAl2020",
                                             "MillerEtAl2012","Rohling2014","Naish2012","PlioSeaNZ")) #%>%

#pdf("PAGES_Fig2a.pdf", width=11.7,height=8.25)
ggplot(Ampl.df,aes(x=Study,y=RSL.m,fill=Study)) +
  geom_violin(trim=TRUE,scale="width", adjust=0.5, colour=NA) + #scale="count", draw_quantiles = c(0.01,0.159,0.50,0.841,0.999)
  scale_fill_manual(values=cols3) + #c("blue4","blue1","dodgerblue3","deepskyblue2","lightgoldenrod4","deeppink4","darkred","red"))+
  stat_summary(fun.data = f,geom="boxplot",width=0.05,fill="white",lwd=1) +
  stat_summary(fun=m,geom="line",aes(group=1),linewidth=2,col="white")+
  stat_summary(fun=m,geom="point",col="black",linewidth=2)+
  theme(axis.title.x=element_blank(), legend.position="none", axis.text.x=element_text(face="bold"),axis.text.y = element_text(face="bold"),panel.grid.major.x  = element_blank())+
  scale_y_continuous(name="Relative sea-level amplitudes (m)",limits=c(0,70),expand=expansion(mult=c(0,0)),breaks = seq(from=0,to=70,by=5))+
  geom_hline(yintercept=22.7, linetype="dashed", color = "black", size=1)+
  annotate("text", x = 7, y = 24, label = "MBIS (AIS)")+
  geom_hline(yintercept=30, linetype="dotted", color = "black", size=1)+ #GIS 7.3 m SLE Bamber et al., 2001
  annotate("text", x = 6.5, y = 31.5, label = "GIS")+
  geom_hline(yintercept=65.6, linetype="solid", color = "black", size=1)+
  annotate("text", x = 6, y = 64.5, label = "Terrestrial AIS")+
  annotate("text", x = 5.5, y = 67, label = "Present-day global ice sheet budget") +
 annotate("text", x = 1, y = 6, label = "paste(italic(n),\"= 82\")", parse = TRUE)+
annotate("text", x = 2, y = 17.5, label = "paste(italic(n),\"=28\")", parse = TRUE)+
annotate("text", x = 3, y = 5, label = "paste(italic(n),\"=132\")", parse = TRUE)+
annotate("text", x = 4, y =6, label = "paste(italic(n),\"=67\")", parse = TRUE)+
annotate("text", x = 5, y = 6, label = "paste(italic(n),\"=331\")", parse = TRUE)+
annotate("text", x = 6, y = 3, label = "paste(italic(n),\"=67\")", parse = TRUE)+
annotate("text", x = 7, y = 3.5, label = "paste(italic(n),\"=151\")", parse = TRUE)
#dev.off()

##### Figure 2b #################
f <- function(x) {
  r <- quantile(x, probs = c(0, 0.1, 0.159, 0.5, 0.841, 0.9,1),na.rm=TRUE)
  names(r) <- c("ymin", "very likely","likely", "median", "unlikely", "very unlikely","ymax")
  r
}
quant<-f(PlioSeaNZ_amp[,2])


CA<-rgb(red=100,green=150,blue=255,alpha=55,maxColorValue = 255)
CB<-rgb(red=0,green=0,blue=255,alpha=255,maxColorValue = 255)

Dimutru=cbind(x=c(1,1),y=c(2.9,20.4)) #Dimitru et al., 2019 range of 10th and 90th percentiles
Hearty<-cbind(x=c(2,2),y=c(16.2,26.3)) #South Africa paleo shorelines Hearty et al., 2020
Enewetak<-cbind(x=c(3,3),y=c(10,32.5)) #Enewetak Atoll Wardlaw and Quinn, 1990 
Miller<- cbind(x=c(4,4),y=c(12,32))#Miller et al., 2012 22 +- 10m
Rovere<-cbind(x=5,y=36.74) #minimum height Atlantic coastal Scarp Rovere et al., 2015, Dowsett and Cronin, 1990
Shakun<- cbind(6,22.7) #Terrestrial AIS limiting Shakun et al., 2018
Naish <- cbind(7,3.4) #Evidence for West Antarctic Ice Sheet limiting Naish et al., 2009
Bertram<- cbind(8,3) #Evidence for Aurora Basin retreat Bertram et al., 2018
Bierman<-cbind(9,6.4) # Greenland ice sheet maximum retreat Bierman et al., 2016
Paleoshorelines<-list(Dimutru,Hearty,Enewetak,Miller,Rovere,Shakun,Naish,Bertram, Bierman)
names(Paleoshorelines)<-c("Dimutru","Hearty","Enewetak", "Miller","Rovere","Shakun", "Naish", "Bertram","Bierman")


#pdf("PAGES_Fig1b_V2.pdf", width=6,height=8.25)
par(bg="lightgrey")
plot.new()
plot.window(xlim=c(0.75,8.5),ylim=c(0,70),xlab="Study",ylab="RSL (m)",yaxs="i",xaxs="i")
axis(1,at=c(1:8),labels=c(names(Paleoshorelines))) #"Dimutru","Hearty","Enewetak", "Dolan", "Koenig","Rovere","Shakun", "Naish", "Bertram"))
axis(2)
panel.min = rect(0,quant["ymin"], 9,quant["very likely"], col=CB, border=NA,lwd=2)
panel.2sd = rect(0,quant["very likely"], 9,quant["likely"], col=CA, border=NA,lwd=2) 
panel.med = rect(0,quant["likely"], 9,quant["unlikely"], col="white", lwd=2,border=NA) 
panel.2sd = rect(0,quant["unlikely"], 9,quant["very unlikely"], col=CA, border=NA,lwd=2)
panel.max = rect(0,quant["very unlikely"], 9,quant["ymax"], col=CB, border=NA,lwd=2)
abline(h=quant["median"],col="red")

lines(Paleoshorelines[[1]],lwd=3,lend=2)
lines(Paleoshorelines[[2]],lwd=3,lend=2)
lines(Paleoshorelines[[3]],lwd=3,lend=2)
lines(Paleoshorelines[[4]],lwd=3,lend=2)
points(Paleoshorelines[[5]],lwd=3,lend=2)
points(Paleoshorelines[[6]],lwd=3,lend=2)
points(Paleoshorelines[[7]],lwd=3,lend=2)
points(Paleoshorelines[[8]],lwd=3,lend=2)
points(Paleoshorelines[[9]],lwd=3,lend=2)


text(6,12.5,round(quant["median"],digits=1), col="black") #50th
text(6,3.5,round(quant["likely"],digits=1), col="black") #33rd
text(6,7,round(quant["very likely"],digits=1), col="black") #10th
text(6,16.5,round(quant["unlikely"],digits=1), col="black") #66th
text(6,21.5,round(quant["very unlikely"],digits=1), col="black") #90th
legend(x=1,y=55,fill=c(col="white",CA,CB),legend=c("PlioSeaNZ Likely (68% PD)","PlioSeaNZ Unlikely (33% PD)","PlioSeaNZ Very-unlikely (1% PD"))

#dev.off()


