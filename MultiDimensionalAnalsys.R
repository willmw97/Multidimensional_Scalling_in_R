#Author: Marshal Will
#Here I demonstrate some important uses of Multidimensional Scaling and measuring 
#distances using methods like the euclidan method

#Analysis on DACA program
#The data here give the number of applications, 
#number of accepted applications, number of approved applications and some additional 
#related statistics on a by state basis (including the District of Columbia) in the United States.


library(corrplot)
library(vegan)
library(ade4)
library(MASS)
library(dplyr)
daca = read.csv("DACA.csv")
#Created 2 more metrics
#Finds the proportion of accepted renewal applicants that were approved
daca.data = mutate(daca,PercOfRenewApprove = ApprovedRenew/AcceptRenew)
head(daca.data$PercOfRenewApprove)
#Finds the proportion of applicants that were accpepted for renewal that weere ultimately approved
daca.data = mutate(daca.data, PercOfReviewAccp = ApproveTotal/AcceptTotal)
head(daca.data$PercOfReviewAccp)
#Creates a metrix MDS with k=2 dimensions

daca.mat = daca.data[,c(8,11,17,18)]
daca.scale = scale(daca.mat)
daca.dist = dist(daca.scale)
daca.dist = dist(daca.dist, method = "euclidean")
daca.iso = isoMDS(daca.dist,k=2)

daca.sammon = sammon(daca.dist,k=2)
plot(daca.sammon$points, type = "n", xlim = c(-20,20), ylim = c(-10,10))
text(daca.sammon$points,as.character(daca$State), cex = .5,col=as.numeric(daca$Lawsuit)+2)
abline(h=0,v=0,lty=2,col="red",lwd=2)


#Creates the same type of graph with k=3 dimensions
daca.scale = scale(daca.mat)
daca.dist = dist(daca.scale)
library(rgl)
daca.sammon3D = sammon(daca.dist,k=3)
plot3d(daca.sammon3D$points,type = "n",xlim = c(-18,15), ylim = c(-3,3))
text3d(daca.sammon3D$points,text = as.character(daca$State),col=as.numeric(daca$Lawsuit)+3,cex=.5)


#Creates a plot based on MDS where bubble size shows size of the percentage of applicants
#who were reviewed and ultametly approved.

#plot(daca.sammon$points,type = "n", main = "MDS plot for percentage of DACA recipients whose application was reviewed was accepeted",
  #   xlim = c(-15,10), ylim = c(-2,2))
symbols(daca.sammon$points,
circles=daca.data$State,inches = .10)
main = "Circle Size of Daca Recipients Percentage Approved"
text(daca.iso$points, labels = as.character(daca.data$State),cex = .5,col="black")
abline(h=0,v=0,lty=2,col="red",lwd=2)


library(ggplot2)
library(tidyverse)
df = as.data.frame(daca.iso)
ggplot(df$points, aes(x=df$points.1, y=df$points.2)) +
  geom_point(aes(size = daca.data$PercOfRenewApprove, alpha=.5), show.legend = FALSE ) + 
  geom_text(aes(label=daca.data$State),hjust=0, vjust=0,size=3)+
  theme_bw()+
  scale_y_continuous(limits=c(-14,14))+
  scale_x_continuous(limits=c(-14,14))
#Size of circle shows the percentage of applicants approved
#The larger the circle the larger percentage of applicants who were approved

#Gene Expression analysis for levels in Colon Tissue Sampels
#These data come from a microarray experiment where gene expression 
#levels were measured for colon tissues samples from 40 individuals 
#with colon cancer and 22 individuals without colon cancer.  
#Below is a blurb about microarray experiments taken from Modern Multivariate 
#Statistical Techniques by Alan Izenman.  

library(corrplot)
library(vegan)
library(ade4)
library(MASS)
#read data set
Alanton = read.csv("Alontop.csv")
Alanton.mat = Alanton[,-93]
str(Alonton.mat)
summary(Alanton.mat)

#Creates a distance matrix for 62 tissue samples using euclidean distance
Alanton.dist = dist(Alanton.mat,method = "euclidean")

summary(as.numeric(Alanton.dist))

#Uses a classical multidimensional scaling using k=2 dimensions. 
#Creates a MDS plot based on the scaling
Alanton.MDS = cmdscale(Alanton.dist,k=2,eig=T)
plot(Alanton.MDS$points,type="n")
text(Alanton.MDS$points,as.character(as.numeric(Alanton$Tissue.Type)),col=as.numeric(Alanton$Tissue.Type)+2)

#Repeats previous part with k = 3
Alanton.MDS = cmdscale(Alanton.dist,k=3,eig=T)
plot(Alanton.MDS$points,type="n")
text(Alanton.MDS$points,as.character(as.numeric(Alanton$Tissue.Type)),col=as.numeric(Alanton$Tissue.Type)+2)


#Creates a correlation matrix using a non-metric MDS with k = 2 dimensions
Alanton.corr = cor(Alanton.mat)

corrplot(Alanton.corr)


#In their paper “Data Analysis in Community and Landscape Ecology” 
#Jongman, R.H.G, ter Braak, C.J.F & van Tongeren, O.F.R., present 
#the results from a study examining 20 dune meadow sites and counting 
#the presence of 30 different species found within these sites.
library(corrplot)
library(vegan)
library(ade4)
library(MASS)
library(proxy)

Dune = read.csv("Dune%20with%20Environmentals.csv")


site.dist = dist(dune.mat,method="Bray")
site.iso = isoMDS(site.dist,k=2)


#Plot for Site

Tdune.mat = t(dune.mat)

Dune.dist = dist(Tdune.mat,method = "Bray")
Dune.iso = isoMDS(Dune.dist,k=2)
Radius = Dune$A1
plot(site.iso$points,type = "n", main = "MDS plot for species")
#symbols(site.iso$points,
      #  circles=Dune$A1,add=T,bg=as.numeric(Dune$Management)+2,inches=.20)
text(Dune.iso$points, as.character(row.names(Tdune.mat)),cex = .75)
abline(h=0,v=0,lty=2,col="red",lwd=2)


plot(site.iso$points,type = "n", main = "MDS plot for species")
#symbols(site.iso$points,
#circles=Dune$A1,add=T,bg=as.numeric(Dune$Management)+2,inches=.20)
        text(site.iso$points, as.character(row.names(dune.mat)),cex = .75)
        abline(h=0,v=0,lty=2,col="red",lwd=2)

#Plot for Moisture
plot(site.iso$points,type = "n", main = "MDS plot for Moisture")
symbols(site.iso$points,
        circles=Dune$Site,add=T,bg=as.numeric(Dune$Moisture[order(Dune$Moisture)])+2,inches=.20)
        text(site.iso$points, as.character(row.names(dune.mat)),cex = .75)
        abline(h=0,v=0,lty=2,col="red",lwd=2)           
    


#Plot for Management
plot(site.iso$points,type = "n", main = "MDS plot for Management")
    symbols(site.iso$points,
    circles=Dune$Site,add=T,bg=as.numeric(Dune$Management)+2,inches=.20)
    text(site.iso$points, as.character(row.names(dune.mat)),cex = .75)
    #abline(h=0,v=0,lty=2,col="red",lwd=2)   
    

#Plot for Use
plot(site.iso$points,type = "n", main = "MDS plot for Use")
    symbols(site.iso$points,
    circles=Dune$Use,add=T,bg=as.numeric(Dune$Use)+2,inches=.10)
    text(site.iso$points, as.character(row.names(dune.mat)),cex = .75)   
    

#Plot for Manure
plot(site.iso$points,type = "n", main = "MDS plot for Manure")
    symbols(site.iso$points,
            circles=Dune$Manure,add=T,bg=as.numeric(Dune$Manure)+2,inches=.10)
    text(site.iso$points, as.character(row.names(dune.mat)),cex = .75)     
    
    
#Plot for Achimill and Callcusp
    
    names(Dune)
    dune.mat = Dune[,7:35]
    rowSums(dune.mat)
    colSums(dune.mat)
    spec.counts = colSums(dune.mat)
    quantile(spec.counts,probs=seq(0,1,.01))
    spec100.mat = spec.mat[,colSums(spec.mat)>100]
    dim(spec100.mat)     
    tspec100.mat = t(spec100.mat)
    spec.dist = dist(tspec100.mat,method="Bray")
    spec.iso = isoMDS(spec.dist,k=2)
    
    site.meta = metaMDS(spec.counts,k=2)
    plot(site.meta,type="n",xlab="Dim 1",ylab="Dim 2",main="Site and Species Ordination")
    text(site.meta,display=c("species"),as.character(names(spec100.mat)),cex=.7,col="green")
    text(site.meta,display=c("site"),as.character(Barro$Site),cex=.5,col="blue")
    abline(h=0,v=0,lty=2,col="red",lwd=2)    
        
        
        






