#install.packages("linpk")
library(linpk)
library(ggplot2)
#Artemether-lumefantrine model
t.obs <- seq(0, 42, 1)

y <- pkprofile(t.obs, cl = 0.1155, vc = 11.2, ka = 1.6, 
               dose = list(t.dose=c(0,0.5,1,1.5,2,2.5), amt=360))

#1-compartment model
y <- pkprofile(t.obs, cl = 32.4, vc = 11.2, ka = 0.9264, 
               dose = list(t.dose=c(0,0.5,1,1.5,2,2.5), amt=240))
plot(y, xlab="Time(day)", ylab="Concentration(mg/l)", main="Artemether-lumefantrine PK plot")
as.data.frame(y)

y_df <- as.data.frame(y)
sp <- ggplot(data=y_df, aes(x=time, y=conc))+
  geom_line() +
  scale_y_continuous(trans='log10') +
  xlab("Time(day)") + ylab("Concentration(mg/l)")
sp <- sp + geom_segment(aes(x=7,xend=28,y=3.950833e-01,yend=3.950833e-01),linetype="dotted") +
  geom_segment(aes(x=7,xend=42,y=1.405343e-09,yend=1.405343e-09),linetype="dotted") +
  geom_segment(aes(x=28,xend=42,y=3.274588e-15,yend=3.274588e-15),linetype="dotted") +
  geom_segment(aes(x=7,xend=28,y=0.1975416507,yend=0.1975416507),linetype="dashed") +
  geom_segment(aes(x=28,xend=42,y=7.0267314e-10,yend=7.0267314e-10),linetype="dashed")
sp <- sp + geom_vline(xintercept = c(7,28,42),  linetype="dotted") 
sp <- sp + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"))


#2-compartment model
y2 <- pkprofile(t.obs, cl = 32.4, vc = 11.2, q = 8.256, vp = 59, ka = 0.9264, 
               dose = list(t.dose=c(0,0.5,1,1.5,2,2.5), amt=360))
plot(y2, xlab="Time(day)", ylab="Concentration(mg/l)", main="Artemether-lumefantrine PK plot")


df <- data.frame(as.data.frame(y),as.data.frame(y2))
names(df) <- c("time","y1","time.1","y2")

library(ggplot2)

g <- ggplot(df, aes(time))
g <- g + geom_line(aes(y=y1), colour="black")
g <- g + geom_line(aes(y=y2), colour="blue")
g <- g + ylab("Concentration(mg/l)") + xlab("Time(day)")
g + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))



#Artesunate-mefloquine model
t.obs <- seq(0, 42, 1)

#mefloquine
#1-compartment model
y <- pkprofile(t.obs, cl = 27.6, vc = 279, ka = 3.984, 
               dose = list(t.dose=c(0,1), amt=c(510,340)))
plot(y, xlab="Time(day)", ylab="Concentration(mg/l)", main="Artesunate-mefloquine PK plot")

#2-compartment model
y2 <- pkprofile(t.obs, cl = 27.6, vc = 279, q=31.92, vp = 341, ka = 3.984,
               dose = list(t.dose=c(0,1), amt=c(510,340)))
plot(y2, xlab="Time(day)", ylab="Concentration(mg/l)", main="Artesunate-mefloquine PK plot")

df <- data.frame(as.data.frame(y),as.data.frame(y2))
names(df) <- c("time","y1","time.1","y2")

#library(ggplot2)

g <- ggplot(df, aes(time))
g <- g + geom_line(aes(y=y1), colour="black")
g <- g + geom_line(aes(y=y2), colour="blue")
g <- g + ylab("Concentration(mg/l)") + xlab("Time(day)")
g + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))

#Artesunate-mefloquine model 2
t.obs <- seq(0, 42, 1)


#Artesunate
y <- pkprofile(t.obs, cl = 3504, vc = 139, ka = 76.8, 
               dose = list(t.dose=c(0,1,2), amt=c(136,136, 136), f = 100))
plot(y,  xlab="Time(day)", ylab="Concentration(mg/l)", main="Artesunate PK plot")
abline(v = 7, col="red", lwd=3, lty=2)
abline(v = 28, col="red", lwd=3, lty=2)
abline(v = 42, col="red", lwd=3, lty=2)

#mefloquine
y2 <- pkprofile(t.obs, cl = 10.8, vc = 95,  ka = 9.6,
                dose = list(t.dose=c(0,1), amt=c(510,340), f = 1))
plot(y2,  xlab="Time(day)", ylab="Concentration(mg/l)", main="Mefloquine PK plot")
abline(v = 7, col="red", lwd=3, lty=2)
abline(v = 28, col="red", lwd=3, lty=2)
abline(v = 42, col="red", lwd=3, lty=2)

y2_df <- as.data.frame(y2)
sp <- ggplot(data=y2_df, aes(x=time, y=conc))+
  geom_line() +
  scale_y_continuous(trans='log10') +
  xlab("Time(day)") + ylab("Concentration(mg/l)")
sp <- sp + geom_segment(aes(x=7,xend=28,y=4.28242033,yend=4.28242033),linetype="dotted") +
  geom_segment(aes(x=7,xend=42,y=0.39343080,yend=0.39343080),linetype="dotted") +
  geom_segment(aes(x=28,xend=42,y=0.08010404,yend=0.08010404),linetype="dotted") +
  geom_segment(aes(x=7,xend=28,y=2.337925565,yend=2.337925565),linetype="dashed") +
  geom_segment(aes(x=28,xend=42,y=0.23676742,yend=0.23676742),linetype="dashed")
sp <- sp + geom_vline(xintercept = c(7,28,42),  linetype="dotted") 
sp <- sp + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"))

df <- data.frame(as.data.frame(y),as.data.frame(y2))
names(df) <- c("time","y1","time.1","y2")

#library(ggplot2)

g <- ggplot(df, aes(time))
g <- g + geom_line(aes(y=y1), colour="black")
g <- g + geom_line(aes(y=y2), colour="blue")
g <- g + coord_trans(y="log10")
g <- g + ylab("Concentration(mg/l)") + xlab("Time(day)")
g + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))

#Artesunate vs. Mefloquine
t.obs <- seq(0, 42, 1)
#Artesunate
y <- pkprofile(t.obs, cl = 3504, vc = 139, ka = 76.8, 
               dose = list(t.dose=c(0,1,2), amt=c(136,136, 136), f = 100))
y_df <- as.data.frame(y)
y_df$conc_relative <- y_df$conc/max(y_df$conc)
y_df$conc_relative[y_df$conc_relative < 0.001] <- 0
#Mefloquine
y2 <- pkprofile(t.obs, cl = 10.8, vc = 95,  ka = 9.6,
                dose = list(t.dose=c(0,1), amt=c(510,340), f = 1))
y2_df <- as.data.frame(y2)
y2_df$conc_relative <- y2_df$conc/max(y2_df$conc)
y2_df$conc_relative[y2_df$conc_relative < 0.001] <- 0

df_combined <- data.frame(y_df,y2_df)

sp <- ggplot(df_combined, aes(time))+
  geom_line(aes(y=conc_relative), color='blue') +
  geom_line(aes(y=conc_relative.1), color='green') +
  scale_y_continuous(trans='log10', limits = c(0.001,1)) +
  #scale_y_continuous(trans='log10') +
  annotation_logticks(sides = "l") +
  xlab("Time(day)") + ylab("Relative Concentration") +
  geom_vline(xintercept = c(7,28,42),  linetype="dotted")

sp + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black"))

#Lumefantrine vs. Mefloquine
t.obs <- seq(0, 42, 1)
#Lumefantrine
y <- pkprofile(t.obs, cl = 32.4, vc = 11.2, ka = 0.9264, 
               dose = list(t.dose=c(0,0.5,1,1.5,2,2.5), amt=240))
y_df <- as.data.frame(y)
y_df$conc_relative <- y_df$conc/max(y_df$conc)
y_df$conc_relative[y_df$conc_relative < 0.001] <- 0
#Mefloquine
y2 <- pkprofile(t.obs, cl = 10.8, vc = 95,  ka = 9.6,
                dose = list(t.dose=c(0,1), amt=c(510,340), f = 1))
y2_df <- as.data.frame(y2)
y2_df$conc_relative <- y2_df$conc/max(y2_df$conc)
y2_df$conc_relative[y2_df$conc_relative < 0.001] <- 0

df_combined <- data.frame(y_df,y2_df)

sp <- ggplot(df_combined, aes(time))+
  geom_line(aes(y=conc_relative), color='cyan') +
  geom_line(aes(y=conc_relative.1), color='green') +
  scale_y_continuous(trans='log10', limits = c(0.001,1)) +
  #scale_y_continuous(trans='log10') +
  annotation_logticks(sides = "l") +
  xlab("Time(day)") + ylab("Relative Concentration") +
  geom_vline(xintercept = c(7,28,42),  linetype="dotted")

sp + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black"))

#Artemether vs. Lumefantrine
t.obs <- seq(0, 42, 1)
#Artemether
y <- pkprofile(t.obs, cl = 21000, vc = 2160, ka = 34.934, 
               dose = list(t.dose=c(0,0.5,1,1.5,2,2.5), amt=40))
y_df <- as.data.frame(y)
y_df$conc_relative <- y_df$conc/max(y_df$conc)
y_df$conc_relative[y_df$conc_relative < 0.001] <- 0

#Lumefantrine
y2 <- pkprofile(t.obs, cl = 32.4, vc = 11.2, ka = 0.9264, 
               dose = list(t.dose=c(0,0.5,1,1.5,2,2.5), amt=240))
y2_df <- as.data.frame(y2)
y2_df$conc_relative <- y2_df$conc/max(y2_df$conc)
y2_df$conc_relative[y2_df$conc_relative < 0.001] <- 0

df_combined <- data.frame(y_df,y2_df)

sp <- ggplot(df_combined, aes(time))+
  geom_line(aes(y=conc_relative), color='magenta') +
  geom_line(aes(y=conc_relative.1), color='cyan') +
  scale_y_continuous(trans='log10', limits = c(0.001,1)) + 
  #scale_y_continuous(trans='log10') + 
  annotation_logticks(sides = "l") +
  xlab("Time(day)") + ylab("Relative Concentration") +
  geom_vline(xintercept = c(7,28,42),  linetype="dotted")

sp + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black"))
