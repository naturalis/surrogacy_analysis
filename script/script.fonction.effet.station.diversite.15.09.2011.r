targeted.station.results<-matrix(0,97,28)  #999

 #LE tableau de résultats à la forme : 
 
#  ALL stations       - station TS249        -station TS250   -station TS251 ........
#    SAIg130m
#    SAI g160m
#    SAI g1100m
#              etc....
#
#


 #ALL STATIONS
 
 #Load du tableau : "site" "poissons"

repertoire<-"C:/Users/Simon/Documents/Travail/SIM2011a/R"
setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
recensement=read.csv2("REQ_site_fonction.IRD.csv")

#J'enlève les doublons
recensement_unique=unique(recensement)

#création d'un tableau de présence absence
tablo_pres_abs<-table(recensement_unique$Code_Site,recensement_unique$Code_CT_SC_M)

#calcul du nombre d'espèce par station
vec_poisson<-apply(tablo_pres_abs,1,sum)

#selection aléatoire d'un pool de site
setwd(repertoire)
setwd(input)
source('script.aleatoire.r')
selection.aleatoire<-matrix(0,dim(tablo_pres_abs)[1],999) #999
for (i in 1:999){ #999
selection.aleatoire[,i]<-run.alea(tablo_pres_abs)}

#Selection des sites selon la diversité de poissons

setwd(repertoire)
 setwd(input)
requette.abondance.poissons<-read.csv2('abondance.midif_sim_output.csv')

input.shannon.abondance<-xtabs(formula=requette.abondance.poissons$Abondance~requette.abondance.poissons$Code_Site+requette.abondance.poissons$Code_CT_SC_M,requette.abondance.poissons)

source('script.shannon.poissons.r')
selection.shannon.poisson<-shannon.poisson(input.shannon.abondance)

#Diversité de poisson incluse dans les sites selectionnés (aléatoirement et selon diversité de poisson)
 source('script.diversite.shannon.selectionnee.fonction.r')
 diversite.shannon.selection.aleatoire<-matrix(0,dim(tablo_pres_abs)[1],999)   #999
 diversite.shannon.selection.aleatoire<-apply(as.matrix(selection.aleatoire),2,function(x) species.shannon.evol(x,requette.abondance.poissons))
 diversite.shannon.moyenne.selection.aleatoire<-apply(diversite.shannon.selection.aleatoire,1,mean)
 
 diversite.selection.shannon.poisson<-species.shannon.evol(selection.shannon.poisson,requette.abondance.poissons)
 
 #Calcul des courbes aléatore moyenne, uper, et low
 diversite.shannon.moyenne.selection.aleatoire
diversite.shannon.borne.sup.selection.aleatoire=apply(diversite.shannon.selection.aleatoire,1,function(x) sort(x)[ceiling(0.975*length(x))])
diversite.shannon.borne.inf.selection.aleatoire=apply(diversite.shannon.selection.aleatoire,1,function(x) sort(x)[ceiling(0.025*length(x))])

#nombre total d'espèces de poissons observés
nombre.poissons.total <- nlevels(factor(recensement_unique$Code_CT_SC_M))



  #############Habitat Geo1###################
  
  
  
   ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances

setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
relation.habitats.sites.geo1<-read.csv2("fc3_intersect_geo1.csv")

sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')

relation.habitats.sites.geo1<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$Code_Site%in%sites.etudies,]

#création d'une nouvelle table site-habitat par distance
relation.habitats.sites.geo1.30m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==30,]
relation.habitats.sites.geo1.60m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==60,]
relation.habitats.sites.geo1.100m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==100,]
relation.habitats.sites.geo1.250m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==250,]
relation.habitats.sites.geo1.500m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==500,]
relation.habitats.sites.geo1.750m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==750,]
relation.habitats.sites.geo1.1000m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==1000,]
relation.habitats.sites.geo1.3000m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==3000,]
relation.habitats.sites.geo1.5000m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==5000,]

#cration d'un tableau adapté à l'entrée du script shannon
relation.habitats.sites.geo1.30m.unique<-aggregate(relation.habitats.sites.geo1.30m$surface,list(Code_Site=relation.habitats.sites.geo1.30m$Code_Site,habitat=relation.habitats.sites.geo1.30m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1.30m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_30m<-xtabs(formula=relation.habitats.sites.geo1.30m.unique$surface~relation.habitats.sites.geo1.30m.unique$Code_Site+relation.habitats.sites.geo1.30m.unique$habitat,relation.habitats.sites.geo1.30m.unique)

relation.habitats.sites.geo1.60m.unique<-aggregate(relation.habitats.sites.geo1.60m$surface,list(Code_Site=relation.habitats.sites.geo1.60m$Code_Site,habitat=relation.habitats.sites.geo1.60m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1.60m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_60m<-xtabs(formula=relation.habitats.sites.geo1.60m.unique$surface~relation.habitats.sites.geo1.60m.unique$Code_Site+relation.habitats.sites.geo1.60m.unique$habitat,relation.habitats.sites.geo1.60m.unique)

relation.habitats.sites.geo1.100m.unique<-aggregate(relation.habitats.sites.geo1.100m$surface,list(Code_Site=relation.habitats.sites.geo1.100m$Code_Site,habitat=relation.habitats.sites.geo1.100m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1.100m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_100m<-xtabs(formula=relation.habitats.sites.geo1.100m.unique$surface~relation.habitats.sites.geo1.100m.unique$Code_Site+relation.habitats.sites.geo1.100m.unique$habitat,relation.habitats.sites.geo1.100m.unique)

relation.habitats.sites.geo1.250m.unique<-aggregate(relation.habitats.sites.geo1.250m$surface,list(Code_Site=relation.habitats.sites.geo1.250m$Code_Site,habitat=relation.habitats.sites.geo1.250m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1.250m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_250m<-xtabs(formula=relation.habitats.sites.geo1.250m.unique$surface~relation.habitats.sites.geo1.250m.unique$Code_Site+relation.habitats.sites.geo1.250m.unique$habitat,relation.habitats.sites.geo1.250m.unique)

relation.habitats.sites.geo1.500m.unique<-aggregate(relation.habitats.sites.geo1.500m$surface,list(Code_Site=relation.habitats.sites.geo1.500m$Code_Site,habitat=relation.habitats.sites.geo1.500m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1.500m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_500m<-xtabs(formula=relation.habitats.sites.geo1.500m.unique$surface~relation.habitats.sites.geo1.500m.unique$Code_Site+relation.habitats.sites.geo1.500m.unique$habitat,relation.habitats.sites.geo1.500m.unique)

relation.habitats.sites.geo1.750m.unique<-aggregate(relation.habitats.sites.geo1.750m$surface,list(Code_Site=relation.habitats.sites.geo1.750m$Code_Site,habitat=relation.habitats.sites.geo1.750m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1.750m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_750m<-xtabs(formula=relation.habitats.sites.geo1.750m.unique$surface~relation.habitats.sites.geo1.750m.unique$Code_Site+relation.habitats.sites.geo1.750m.unique$habitat,relation.habitats.sites.geo1.750m.unique)

relation.habitats.sites.geo1.1000m.unique<-aggregate(relation.habitats.sites.geo1.1000m$surface,list(Code_Site=relation.habitats.sites.geo1.1000m$Code_Site,habitat=relation.habitats.sites.geo1.1000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1.1000m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_1000m<-xtabs(formula=relation.habitats.sites.geo1.1000m.unique$surface~relation.habitats.sites.geo1.1000m.unique$Code_Site+relation.habitats.sites.geo1.1000m.unique$habitat,relation.habitats.sites.geo1.1000m.unique)

relation.habitats.sites.geo1.3000m.unique<-aggregate(relation.habitats.sites.geo1.3000m$surface,list(Code_Site=relation.habitats.sites.geo1.3000m$Code_Site,habitat=relation.habitats.sites.geo1.3000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1.3000m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_3000m<-xtabs(formula=relation.habitats.sites.geo1.3000m.unique$surface~relation.habitats.sites.geo1.3000m.unique$Code_Site+relation.habitats.sites.geo1.3000m.unique$habitat,relation.habitats.sites.geo1.3000m.unique)

relation.habitats.sites.geo1.5000m.unique<-aggregate(relation.habitats.sites.geo1.5000m$surface,list(Code_Site=relation.habitats.sites.geo1.5000m$Code_Site,habitat=relation.habitats.sites.geo1.5000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1.5000m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_5000m<-xtabs(formula=relation.habitats.sites.geo1.5000m.unique$surface~relation.habitats.sites.geo1.5000m.unique$Code_Site+relation.habitats.sites.geo1.5000m.unique$habitat,relation.habitats.sites.geo1.5000m.unique)

tablo_geo1_30m<-tablo_geo1_30m[25:51,]
tablo_geo1_60m<-tablo_geo1_60m[25:51,]
tablo_geo1_100m<-tablo_geo1_100m[25:51,]
tablo_geo1_250m<-tablo_geo1_250m[25:51,]
tablo_geo1_500m<-tablo_geo1_500m[25:51,]
tablo_geo1_750m<-tablo_geo1_750m[25:51,]
tablo_geo1_1000m<-tablo_geo1_1000m[25:51,]
tablo_geo1_3000m<-tablo_geo1_3000m[25:51,]
tablo_geo1_5000m<-tablo_geo1_5000m[25:51,]

#Load de l'algorithme shannon
source("script.shannon.r")

#tirage des sites poissons en optimisant la diversité shannon d'habitat
 selection.shannon.hab.geo1.30m<-shannon(tablo_geo1_30m)
 selection.shannon.hab.geo1.60m<-shannon(tablo_geo1_60m)
 selection.shannon.hab.geo1.100m<-shannon(tablo_geo1_100m)
 selection.shannon.hab.geo1.250m<-shannon(tablo_geo1_250m)
 selection.shannon.hab.geo1.500m<-shannon(tablo_geo1_500m)
 selection.shannon.hab.geo1.750m<-shannon(tablo_geo1_750m)
 selection.shannon.hab.geo1.1000m<-shannon(tablo_geo1_1000m)
 selection.shannon.hab.geo1.3000m<-shannon(tablo_geo1_3000m)
 selection.shannon.hab.geo1.5000m<-shannon(tablo_geo1_5000m)

  #DIVERSITE SELECTIONNE PAR LE SCRIPT SHANNON
  source('script.diversite.shannon.selectionnee.fonction.r')
diversite.selection.shannon.hab.geo1.30m<-species.shannon.evol(selection.shannon.hab.geo1.30m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1.60m<-species.shannon.evol(selection.shannon.hab.geo1.60m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1.100m<-species.shannon.evol(selection.shannon.hab.geo1.100m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1.250m<-species.shannon.evol(selection.shannon.hab.geo1.250m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1.500m<-species.shannon.evol(selection.shannon.hab.geo1.500m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1.750m<-species.shannon.evol(selection.shannon.hab.geo1.750m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1.1000m<-species.shannon.evol(selection.shannon.hab.geo1.1000m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1.3000m<-species.shannon.evol(selection.shannon.hab.geo1.3000m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1.5000m<-species.shannon.evol(selection.shannon.hab.geo1.5000m,requette.abondance.poissons)
  
 #CALCUL DES SAI#
source('script.sai.r')
SAI_shannon.geo1.30m<-SAI(diversite.selection.shannon.hab.geo1.30m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1.60m<-SAI(diversite.selection.shannon.hab.geo1.60m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1.100m<-SAI(diversite.selection.shannon.hab.geo1.100m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1.250m<-SAI(diversite.selection.shannon.hab.geo1.250m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1.500m<-SAI(diversite.selection.shannon.hab.geo1.500m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1.750m<-SAI(diversite.selection.shannon.hab.geo1.750m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1.1000m<-SAI(diversite.selection.shannon.hab.geo1.1000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1.3000m<-SAI(diversite.selection.shannon.hab.geo1.3000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1.5000m<-SAI(diversite.selection.shannon.hab.geo1.5000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)


  #############Habitat geo2###################



   ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances

setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
relation.habitats.sites.geo2<-read.csv2("fc3_intersect_geo2.csv")

sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')

relation.habitats.sites.geo2<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$Code_Site%in%sites.etudies,]

#création d'une nouvelle table site-habitat par distance
relation.habitats.sites.geo2.30m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==30,]
relation.habitats.sites.geo2.60m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==60,]
relation.habitats.sites.geo2.100m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==100,]
relation.habitats.sites.geo2.250m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==250,]
relation.habitats.sites.geo2.500m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==500,]
relation.habitats.sites.geo2.750m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==750,]
relation.habitats.sites.geo2.1000m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==1000,]
relation.habitats.sites.geo2.3000m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==3000,]
relation.habitats.sites.geo2.5000m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==5000,]

#cration d'un tableau adapté à l'entrée du script shannon
relation.habitats.sites.geo2.30m.unique<-aggregate(relation.habitats.sites.geo2.30m$surface,list(Code_Site=relation.habitats.sites.geo2.30m$Code_Site,habitat=relation.habitats.sites.geo2.30m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo2.30m.unique)=c("Code_Site","habitat","surface")
tablo_geo2_30m<-xtabs(formula=relation.habitats.sites.geo2.30m.unique$surface~relation.habitats.sites.geo2.30m.unique$Code_Site+relation.habitats.sites.geo2.30m.unique$habitat,relation.habitats.sites.geo2.30m.unique)

relation.habitats.sites.geo2.60m.unique<-aggregate(relation.habitats.sites.geo2.60m$surface,list(Code_Site=relation.habitats.sites.geo2.60m$Code_Site,habitat=relation.habitats.sites.geo2.60m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo2.60m.unique)=c("Code_Site","habitat","surface")
tablo_geo2_60m<-xtabs(formula=relation.habitats.sites.geo2.60m.unique$surface~relation.habitats.sites.geo2.60m.unique$Code_Site+relation.habitats.sites.geo2.60m.unique$habitat,relation.habitats.sites.geo2.60m.unique)

relation.habitats.sites.geo2.100m.unique<-aggregate(relation.habitats.sites.geo2.100m$surface,list(Code_Site=relation.habitats.sites.geo2.100m$Code_Site,habitat=relation.habitats.sites.geo2.100m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo2.100m.unique)=c("Code_Site","habitat","surface")
tablo_geo2_100m<-xtabs(formula=relation.habitats.sites.geo2.100m.unique$surface~relation.habitats.sites.geo2.100m.unique$Code_Site+relation.habitats.sites.geo2.100m.unique$habitat,relation.habitats.sites.geo2.100m.unique)

relation.habitats.sites.geo2.250m.unique<-aggregate(relation.habitats.sites.geo2.250m$surface,list(Code_Site=relation.habitats.sites.geo2.250m$Code_Site,habitat=relation.habitats.sites.geo2.250m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo2.250m.unique)=c("Code_Site","habitat","surface")
tablo_geo2_250m<-xtabs(formula=relation.habitats.sites.geo2.250m.unique$surface~relation.habitats.sites.geo2.250m.unique$Code_Site+relation.habitats.sites.geo2.250m.unique$habitat,relation.habitats.sites.geo2.250m.unique)

relation.habitats.sites.geo2.500m.unique<-aggregate(relation.habitats.sites.geo2.500m$surface,list(Code_Site=relation.habitats.sites.geo2.500m$Code_Site,habitat=relation.habitats.sites.geo2.500m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo2.500m.unique)=c("Code_Site","habitat","surface")
tablo_geo2_500m<-xtabs(formula=relation.habitats.sites.geo2.500m.unique$surface~relation.habitats.sites.geo2.500m.unique$Code_Site+relation.habitats.sites.geo2.500m.unique$habitat,relation.habitats.sites.geo2.500m.unique)

relation.habitats.sites.geo2.750m.unique<-aggregate(relation.habitats.sites.geo2.750m$surface,list(Code_Site=relation.habitats.sites.geo2.750m$Code_Site,habitat=relation.habitats.sites.geo2.750m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo2.750m.unique)=c("Code_Site","habitat","surface")
tablo_geo2_750m<-xtabs(formula=relation.habitats.sites.geo2.750m.unique$surface~relation.habitats.sites.geo2.750m.unique$Code_Site+relation.habitats.sites.geo2.750m.unique$habitat,relation.habitats.sites.geo2.750m.unique)

relation.habitats.sites.geo2.1000m.unique<-aggregate(relation.habitats.sites.geo2.1000m$surface,list(Code_Site=relation.habitats.sites.geo2.1000m$Code_Site,habitat=relation.habitats.sites.geo2.1000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo2.1000m.unique)=c("Code_Site","habitat","surface")
tablo_geo2_1000m<-xtabs(formula=relation.habitats.sites.geo2.1000m.unique$surface~relation.habitats.sites.geo2.1000m.unique$Code_Site+relation.habitats.sites.geo2.1000m.unique$habitat,relation.habitats.sites.geo2.1000m.unique)

relation.habitats.sites.geo2.3000m.unique<-aggregate(relation.habitats.sites.geo2.3000m$surface,list(Code_Site=relation.habitats.sites.geo2.3000m$Code_Site,habitat=relation.habitats.sites.geo2.3000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo2.3000m.unique)=c("Code_Site","habitat","surface")
tablo_geo2_3000m<-xtabs(formula=relation.habitats.sites.geo2.3000m.unique$surface~relation.habitats.sites.geo2.3000m.unique$Code_Site+relation.habitats.sites.geo2.3000m.unique$habitat,relation.habitats.sites.geo2.3000m.unique)

relation.habitats.sites.geo2.5000m.unique<-aggregate(relation.habitats.sites.geo2.5000m$surface,list(Code_Site=relation.habitats.sites.geo2.5000m$Code_Site,habitat=relation.habitats.sites.geo2.5000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo2.5000m.unique)=c("Code_Site","habitat","surface")
tablo_geo2_5000m<-xtabs(formula=relation.habitats.sites.geo2.5000m.unique$surface~relation.habitats.sites.geo2.5000m.unique$Code_Site+relation.habitats.sites.geo2.5000m.unique$habitat,relation.habitats.sites.geo2.5000m.unique)

tablo_geo2_30m<-tablo_geo2_30m[25:51,]
tablo_geo2_60m<-tablo_geo2_60m[25:51,]
tablo_geo2_100m<-tablo_geo2_100m[25:51,]
tablo_geo2_250m<-tablo_geo2_250m[25:51,]
tablo_geo2_500m<-tablo_geo2_500m[25:51,]
tablo_geo2_750m<-tablo_geo2_750m[25:51,]
tablo_geo2_1000m<-tablo_geo2_1000m[25:51,]
tablo_geo2_3000m<-tablo_geo2_3000m[25:51,]
tablo_geo2_5000m<-tablo_geo2_5000m[25:51,]

#Load de l'algorithme shannon
source("script.shannon.r")

#tirage des sites poissons en optimisant la diversité shannon d'habitat
 selection.shannon.hab.geo2.30m<-shannon(tablo_geo2_30m)
 selection.shannon.hab.geo2.60m<-shannon(tablo_geo2_60m)
 selection.shannon.hab.geo2.100m<-shannon(tablo_geo2_100m)
 selection.shannon.hab.geo2.250m<-shannon(tablo_geo2_250m)
 selection.shannon.hab.geo2.500m<-shannon(tablo_geo2_500m)
 selection.shannon.hab.geo2.750m<-shannon(tablo_geo2_750m)
 selection.shannon.hab.geo2.1000m<-shannon(tablo_geo2_1000m)
 selection.shannon.hab.geo2.3000m<-shannon(tablo_geo2_3000m)
 selection.shannon.hab.geo2.5000m<-shannon(tablo_geo2_5000m)

  #DIVERSITE SELECTIONNE PAR LE SCRIPT SHANNON
  source('script.diversite.shannon.selectionnee.fonction.r')
diversite.selection.shannon.hab.geo2.30m<-species.shannon.evol(selection.shannon.hab.geo2.30m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo2.60m<-species.shannon.evol(selection.shannon.hab.geo2.60m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo2.100m<-species.shannon.evol(selection.shannon.hab.geo2.100m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo2.250m<-species.shannon.evol(selection.shannon.hab.geo2.250m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo2.500m<-species.shannon.evol(selection.shannon.hab.geo2.500m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo2.750m<-species.shannon.evol(selection.shannon.hab.geo2.750m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo2.1000m<-species.shannon.evol(selection.shannon.hab.geo2.1000m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo2.3000m<-species.shannon.evol(selection.shannon.hab.geo2.3000m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo2.5000m<-species.shannon.evol(selection.shannon.hab.geo2.5000m,requette.abondance.poissons)

 #CALCUL DES SAI#
source('script.sai.r')
SAI_shannon.geo2.30m<-SAI(diversite.selection.shannon.hab.geo2.30m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo2.60m<-SAI(diversite.selection.shannon.hab.geo2.60m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo2.100m<-SAI(diversite.selection.shannon.hab.geo2.100m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo2.250m<-SAI(diversite.selection.shannon.hab.geo2.250m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo2.500m<-SAI(diversite.selection.shannon.hab.geo2.500m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo2.750m<-SAI(diversite.selection.shannon.hab.geo2.750m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo2.1000m<-SAI(diversite.selection.shannon.hab.geo2.1000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo2.3000m<-SAI(diversite.selection.shannon.hab.geo2.3000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo2.5000m<-SAI(diversite.selection.shannon.hab.geo2.5000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)


  #############Habitat geo3###################



   ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances

setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
relation.habitats.sites.geo3<-read.csv2("fc3_intersect_geo3.csv")

sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')

relation.habitats.sites.geo3<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$Code_Site%in%sites.etudies,]

#création d'une nouvelle table site-habitat par distance

relation.habitats.sites.geo3.100m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==100,]
relation.habitats.sites.geo3.250m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==250,]
relation.habitats.sites.geo3.500m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==500,]
relation.habitats.sites.geo3.750m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==750,]
relation.habitats.sites.geo3.1000m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==1000,]
relation.habitats.sites.geo3.3000m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==3000,]
relation.habitats.sites.geo3.5000m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==5000,]

#cration d'un tableau adapté à l'entrée du script shannon

relation.habitats.sites.geo3.100m.unique<-aggregate(relation.habitats.sites.geo3.100m$surface,list(Code_Site=relation.habitats.sites.geo3.100m$Code_Site,habitat=relation.habitats.sites.geo3.100m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo3.100m.unique)=c("Code_Site","habitat","surface")
tablo_geo3_100m<-xtabs(formula=relation.habitats.sites.geo3.100m.unique$surface~relation.habitats.sites.geo3.100m.unique$Code_Site+relation.habitats.sites.geo3.100m.unique$habitat,relation.habitats.sites.geo3.100m.unique)

relation.habitats.sites.geo3.250m.unique<-aggregate(relation.habitats.sites.geo3.250m$surface,list(Code_Site=relation.habitats.sites.geo3.250m$Code_Site,habitat=relation.habitats.sites.geo3.250m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo3.250m.unique)=c("Code_Site","habitat","surface")
tablo_geo3_250m<-xtabs(formula=relation.habitats.sites.geo3.250m.unique$surface~relation.habitats.sites.geo3.250m.unique$Code_Site+relation.habitats.sites.geo3.250m.unique$habitat,relation.habitats.sites.geo3.250m.unique)

relation.habitats.sites.geo3.500m.unique<-aggregate(relation.habitats.sites.geo3.500m$surface,list(Code_Site=relation.habitats.sites.geo3.500m$Code_Site,habitat=relation.habitats.sites.geo3.500m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo3.500m.unique)=c("Code_Site","habitat","surface")
tablo_geo3_500m<-xtabs(formula=relation.habitats.sites.geo3.500m.unique$surface~relation.habitats.sites.geo3.500m.unique$Code_Site+relation.habitats.sites.geo3.500m.unique$habitat,relation.habitats.sites.geo3.500m.unique)

relation.habitats.sites.geo3.750m.unique<-aggregate(relation.habitats.sites.geo3.750m$surface,list(Code_Site=relation.habitats.sites.geo3.750m$Code_Site,habitat=relation.habitats.sites.geo3.750m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo3.750m.unique)=c("Code_Site","habitat","surface")
tablo_geo3_750m<-xtabs(formula=relation.habitats.sites.geo3.750m.unique$surface~relation.habitats.sites.geo3.750m.unique$Code_Site+relation.habitats.sites.geo3.750m.unique$habitat,relation.habitats.sites.geo3.750m.unique)

relation.habitats.sites.geo3.1000m.unique<-aggregate(relation.habitats.sites.geo3.1000m$surface,list(Code_Site=relation.habitats.sites.geo3.1000m$Code_Site,habitat=relation.habitats.sites.geo3.1000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo3.1000m.unique)=c("Code_Site","habitat","surface")
tablo_geo3_1000m<-xtabs(formula=relation.habitats.sites.geo3.1000m.unique$surface~relation.habitats.sites.geo3.1000m.unique$Code_Site+relation.habitats.sites.geo3.1000m.unique$habitat,relation.habitats.sites.geo3.1000m.unique)

relation.habitats.sites.geo3.3000m.unique<-aggregate(relation.habitats.sites.geo3.3000m$surface,list(Code_Site=relation.habitats.sites.geo3.3000m$Code_Site,habitat=relation.habitats.sites.geo3.3000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo3.3000m.unique)=c("Code_Site","habitat","surface")
tablo_geo3_3000m<-xtabs(formula=relation.habitats.sites.geo3.3000m.unique$surface~relation.habitats.sites.geo3.3000m.unique$Code_Site+relation.habitats.sites.geo3.3000m.unique$habitat,relation.habitats.sites.geo3.3000m.unique)

relation.habitats.sites.geo3.5000m.unique<-aggregate(relation.habitats.sites.geo3.5000m$surface,list(Code_Site=relation.habitats.sites.geo3.5000m$Code_Site,habitat=relation.habitats.sites.geo3.5000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo3.5000m.unique)=c("Code_Site","habitat","surface")
tablo_geo3_5000m<-xtabs(formula=relation.habitats.sites.geo3.5000m.unique$surface~relation.habitats.sites.geo3.5000m.unique$Code_Site+relation.habitats.sites.geo3.5000m.unique$habitat,relation.habitats.sites.geo3.5000m.unique)


tablo_geo3_100m<-tablo_geo3_100m[25:51,]
tablo_geo3_250m<-tablo_geo3_250m[25:51,]
tablo_geo3_500m<-tablo_geo3_500m[25:51,]
tablo_geo3_750m<-tablo_geo3_750m[25:51,]
tablo_geo3_1000m<-tablo_geo3_1000m[25:51,]
tablo_geo3_3000m<-tablo_geo3_3000m[25:51,]
tablo_geo3_5000m<-tablo_geo3_5000m[25:51,]

#Load de l'algorithme shannon
source("script.shannon.r")

#tirage des sites poissons en optimisant la diversité shannon d'habitat

 selection.shannon.hab.geo3.100m<-shannon(tablo_geo3_100m)
 selection.shannon.hab.geo3.250m<-shannon(tablo_geo3_250m)
 selection.shannon.hab.geo3.500m<-shannon(tablo_geo3_500m)
 selection.shannon.hab.geo3.750m<-shannon(tablo_geo3_750m)
 selection.shannon.hab.geo3.1000m<-shannon(tablo_geo3_1000m)
 selection.shannon.hab.geo3.3000m<-shannon(tablo_geo3_3000m)
 selection.shannon.hab.geo3.5000m<-shannon(tablo_geo3_5000m)

  #DIVERSITE SELECTIONNE PAR LE SCRIPT SHANNON
  source('script.diversite.shannon.selectionnee.fonction.r')

diversite.selection.shannon.hab.geo3.100m<-species.shannon.evol(selection.shannon.hab.geo3.100m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo3.250m<-species.shannon.evol(selection.shannon.hab.geo3.250m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo3.500m<-species.shannon.evol(selection.shannon.hab.geo3.500m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo3.750m<-species.shannon.evol(selection.shannon.hab.geo3.750m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo3.1000m<-species.shannon.evol(selection.shannon.hab.geo3.1000m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo3.3000m<-species.shannon.evol(selection.shannon.hab.geo3.3000m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo3.5000m<-species.shannon.evol(selection.shannon.hab.geo3.5000m,requette.abondance.poissons)

 #CALCUL DES SAI#
source('script.sai.r')
SAI_shannon.geo3.100m<-SAI(diversite.selection.shannon.hab.geo3.100m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo3.250m<-SAI(diversite.selection.shannon.hab.geo3.250m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo3.500m<-SAI(diversite.selection.shannon.hab.geo3.500m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo3.750m<-SAI(diversite.selection.shannon.hab.geo3.750m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo3.1000m<-SAI(diversite.selection.shannon.hab.geo3.1000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo3.3000m<-SAI(diversite.selection.shannon.hab.geo3.3000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo3.5000m<-SAI(diversite.selection.shannon.hab.geo3.5000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)


  #############Habitat bent1###################



   ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances

setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
relation.habitats.sites.bent1<-read.csv2("fc3_intersect_bent1.csv")

sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')

relation.habitats.sites.bent1<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$Code_Site%in%sites.etudies,]

#création d'une nouvelle table site-habitat par distance
relation.habitats.sites.bent1.30m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==30,]
relation.habitats.sites.bent1.60m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==60,]
relation.habitats.sites.bent1.100m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==100,]
relation.habitats.sites.bent1.250m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==250,]
relation.habitats.sites.bent1.500m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==500,]
relation.habitats.sites.bent1.750m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==750,]
relation.habitats.sites.bent1.1000m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==1000,]
relation.habitats.sites.bent1.3000m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==3000,]
relation.habitats.sites.bent1.5000m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==5000,]

#cration d'un tableau adapté à l'entrée du script shannon
relation.habitats.sites.bent1.30m.unique<-aggregate(relation.habitats.sites.bent1.30m$surface,list(Code_Site=relation.habitats.sites.bent1.30m$Code_Site,habitat=relation.habitats.sites.bent1.30m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent1.30m.unique)=c("Code_Site","habitat","surface")
tablo_bent1_30m<-xtabs(formula=relation.habitats.sites.bent1.30m.unique$surface~relation.habitats.sites.bent1.30m.unique$Code_Site+relation.habitats.sites.bent1.30m.unique$habitat,relation.habitats.sites.bent1.30m.unique)

relation.habitats.sites.bent1.60m.unique<-aggregate(relation.habitats.sites.bent1.60m$surface,list(Code_Site=relation.habitats.sites.bent1.60m$Code_Site,habitat=relation.habitats.sites.bent1.60m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent1.60m.unique)=c("Code_Site","habitat","surface")
tablo_bent1_60m<-xtabs(formula=relation.habitats.sites.bent1.60m.unique$surface~relation.habitats.sites.bent1.60m.unique$Code_Site+relation.habitats.sites.bent1.60m.unique$habitat,relation.habitats.sites.bent1.60m.unique)

relation.habitats.sites.bent1.100m.unique<-aggregate(relation.habitats.sites.bent1.100m$surface,list(Code_Site=relation.habitats.sites.bent1.100m$Code_Site,habitat=relation.habitats.sites.bent1.100m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent1.100m.unique)=c("Code_Site","habitat","surface")
tablo_bent1_100m<-xtabs(formula=relation.habitats.sites.bent1.100m.unique$surface~relation.habitats.sites.bent1.100m.unique$Code_Site+relation.habitats.sites.bent1.100m.unique$habitat,relation.habitats.sites.bent1.100m.unique)

relation.habitats.sites.bent1.250m.unique<-aggregate(relation.habitats.sites.bent1.250m$surface,list(Code_Site=relation.habitats.sites.bent1.250m$Code_Site,habitat=relation.habitats.sites.bent1.250m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent1.250m.unique)=c("Code_Site","habitat","surface")
tablo_bent1_250m<-xtabs(formula=relation.habitats.sites.bent1.250m.unique$surface~relation.habitats.sites.bent1.250m.unique$Code_Site+relation.habitats.sites.bent1.250m.unique$habitat,relation.habitats.sites.bent1.250m.unique)

relation.habitats.sites.bent1.500m.unique<-aggregate(relation.habitats.sites.bent1.500m$surface,list(Code_Site=relation.habitats.sites.bent1.500m$Code_Site,habitat=relation.habitats.sites.bent1.500m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent1.500m.unique)=c("Code_Site","habitat","surface")
tablo_bent1_500m<-xtabs(formula=relation.habitats.sites.bent1.500m.unique$surface~relation.habitats.sites.bent1.500m.unique$Code_Site+relation.habitats.sites.bent1.500m.unique$habitat,relation.habitats.sites.bent1.500m.unique)

relation.habitats.sites.bent1.750m.unique<-aggregate(relation.habitats.sites.bent1.750m$surface,list(Code_Site=relation.habitats.sites.bent1.750m$Code_Site,habitat=relation.habitats.sites.bent1.750m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent1.750m.unique)=c("Code_Site","habitat","surface")
tablo_bent1_750m<-xtabs(formula=relation.habitats.sites.bent1.750m.unique$surface~relation.habitats.sites.bent1.750m.unique$Code_Site+relation.habitats.sites.bent1.750m.unique$habitat,relation.habitats.sites.bent1.750m.unique)

relation.habitats.sites.bent1.1000m.unique<-aggregate(relation.habitats.sites.bent1.1000m$surface,list(Code_Site=relation.habitats.sites.bent1.1000m$Code_Site,habitat=relation.habitats.sites.bent1.1000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent1.1000m.unique)=c("Code_Site","habitat","surface")
tablo_bent1_1000m<-xtabs(formula=relation.habitats.sites.bent1.1000m.unique$surface~relation.habitats.sites.bent1.1000m.unique$Code_Site+relation.habitats.sites.bent1.1000m.unique$habitat,relation.habitats.sites.bent1.1000m.unique)

relation.habitats.sites.bent1.3000m.unique<-aggregate(relation.habitats.sites.bent1.3000m$surface,list(Code_Site=relation.habitats.sites.bent1.3000m$Code_Site,habitat=relation.habitats.sites.bent1.3000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent1.3000m.unique)=c("Code_Site","habitat","surface")
tablo_bent1_3000m<-xtabs(formula=relation.habitats.sites.bent1.3000m.unique$surface~relation.habitats.sites.bent1.3000m.unique$Code_Site+relation.habitats.sites.bent1.3000m.unique$habitat,relation.habitats.sites.bent1.3000m.unique)

relation.habitats.sites.bent1.5000m.unique<-aggregate(relation.habitats.sites.bent1.5000m$surface,list(Code_Site=relation.habitats.sites.bent1.5000m$Code_Site,habitat=relation.habitats.sites.bent1.5000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent1.5000m.unique)=c("Code_Site","habitat","surface")
tablo_bent1_5000m<-xtabs(formula=relation.habitats.sites.bent1.5000m.unique$surface~relation.habitats.sites.bent1.5000m.unique$Code_Site+relation.habitats.sites.bent1.5000m.unique$habitat,relation.habitats.sites.bent1.5000m.unique)

tablo_bent1_30m<-tablo_bent1_30m[25:51,]
tablo_bent1_60m<-tablo_bent1_60m[25:51,]
tablo_bent1_100m<-tablo_bent1_100m[25:51,]
tablo_bent1_250m<-tablo_bent1_250m[25:51,]
tablo_bent1_500m<-tablo_bent1_500m[25:51,]
tablo_bent1_750m<-tablo_bent1_750m[25:51,]
tablo_bent1_1000m<-tablo_bent1_1000m[25:51,]
tablo_bent1_3000m<-tablo_bent1_3000m[25:51,]
tablo_bent1_5000m<-tablo_bent1_5000m[25:51,]

#Load de l'algorithme shannon
source("script.shannon.r")

#tirage des sites poissons en optimisant la diversité shannon d'habitat
 selection.shannon.hab.bent1.30m<-shannon(tablo_bent1_30m)
 selection.shannon.hab.bent1.60m<-shannon(tablo_bent1_60m)
 selection.shannon.hab.bent1.100m<-shannon(tablo_bent1_100m)
 selection.shannon.hab.bent1.250m<-shannon(tablo_bent1_250m)
 selection.shannon.hab.bent1.500m<-shannon(tablo_bent1_500m)
 selection.shannon.hab.bent1.750m<-shannon(tablo_bent1_750m)
 selection.shannon.hab.bent1.1000m<-shannon(tablo_bent1_1000m)
 selection.shannon.hab.bent1.3000m<-shannon(tablo_bent1_3000m)
 selection.shannon.hab.bent1.5000m<-shannon(tablo_bent1_5000m)

  #DIVERSITE SELECTIONNE PAR LE SCRIPT SHANNON
  source('script.diversite.shannon.selectionnee.fonction.r')
diversite.selection.shannon.hab.bent1.30m<-species.shannon.evol(selection.shannon.hab.bent1.30m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent1.60m<-species.shannon.evol(selection.shannon.hab.bent1.60m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent1.100m<-species.shannon.evol(selection.shannon.hab.bent1.100m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent1.250m<-species.shannon.evol(selection.shannon.hab.bent1.250m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent1.500m<-species.shannon.evol(selection.shannon.hab.bent1.500m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent1.750m<-species.shannon.evol(selection.shannon.hab.bent1.750m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent1.1000m<-species.shannon.evol(selection.shannon.hab.bent1.1000m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent1.3000m<-species.shannon.evol(selection.shannon.hab.bent1.3000m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent1.5000m<-species.shannon.evol(selection.shannon.hab.bent1.5000m,requette.abondance.poissons)

 #CALCUL DES SAI#
source('script.sai.r')
SAI_shannon.bent1.30m<-SAI(diversite.selection.shannon.hab.bent1.30m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent1.60m<-SAI(diversite.selection.shannon.hab.bent1.60m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent1.100m<-SAI(diversite.selection.shannon.hab.bent1.100m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent1.250m<-SAI(diversite.selection.shannon.hab.bent1.250m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent1.500m<-SAI(diversite.selection.shannon.hab.bent1.500m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent1.750m<-SAI(diversite.selection.shannon.hab.bent1.750m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent1.1000m<-SAI(diversite.selection.shannon.hab.bent1.1000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent1.3000m<-SAI(diversite.selection.shannon.hab.bent1.3000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent1.5000m<-SAI(diversite.selection.shannon.hab.bent1.5000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)


  #############Habitat bent2###################



   ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances

setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
relation.habitats.sites.bent2<-read.csv2("fc3_intersect_bent2.csv")

sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')

relation.habitats.sites.bent2<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$Code_Site%in%sites.etudies,]

#création d'une nouvelle table site-habitat par distance
relation.habitats.sites.bent2.30m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==30,]
relation.habitats.sites.bent2.60m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==60,]
relation.habitats.sites.bent2.100m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==100,]
relation.habitats.sites.bent2.250m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==250,]
relation.habitats.sites.bent2.500m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==500,]
relation.habitats.sites.bent2.750m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==750,]
relation.habitats.sites.bent2.1000m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==1000,]
relation.habitats.sites.bent2.3000m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==3000,]
relation.habitats.sites.bent2.5000m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==5000,]

#cration d'un tableau adapté à l'entrée du script shannon
relation.habitats.sites.bent2.30m.unique<-aggregate(relation.habitats.sites.bent2.30m$surface,list(Code_Site=relation.habitats.sites.bent2.30m$Code_Site,habitat=relation.habitats.sites.bent2.30m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent2.30m.unique)=c("Code_Site","habitat","surface")
tablo_bent2_30m<-xtabs(formula=relation.habitats.sites.bent2.30m.unique$surface~relation.habitats.sites.bent2.30m.unique$Code_Site+relation.habitats.sites.bent2.30m.unique$habitat,relation.habitats.sites.bent2.30m.unique)

relation.habitats.sites.bent2.60m.unique<-aggregate(relation.habitats.sites.bent2.60m$surface,list(Code_Site=relation.habitats.sites.bent2.60m$Code_Site,habitat=relation.habitats.sites.bent2.60m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent2.60m.unique)=c("Code_Site","habitat","surface")
tablo_bent2_60m<-xtabs(formula=relation.habitats.sites.bent2.60m.unique$surface~relation.habitats.sites.bent2.60m.unique$Code_Site+relation.habitats.sites.bent2.60m.unique$habitat,relation.habitats.sites.bent2.60m.unique)

relation.habitats.sites.bent2.100m.unique<-aggregate(relation.habitats.sites.bent2.100m$surface,list(Code_Site=relation.habitats.sites.bent2.100m$Code_Site,habitat=relation.habitats.sites.bent2.100m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent2.100m.unique)=c("Code_Site","habitat","surface")
tablo_bent2_100m<-xtabs(formula=relation.habitats.sites.bent2.100m.unique$surface~relation.habitats.sites.bent2.100m.unique$Code_Site+relation.habitats.sites.bent2.100m.unique$habitat,relation.habitats.sites.bent2.100m.unique)

relation.habitats.sites.bent2.250m.unique<-aggregate(relation.habitats.sites.bent2.250m$surface,list(Code_Site=relation.habitats.sites.bent2.250m$Code_Site,habitat=relation.habitats.sites.bent2.250m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent2.250m.unique)=c("Code_Site","habitat","surface")
tablo_bent2_250m<-xtabs(formula=relation.habitats.sites.bent2.250m.unique$surface~relation.habitats.sites.bent2.250m.unique$Code_Site+relation.habitats.sites.bent2.250m.unique$habitat,relation.habitats.sites.bent2.250m.unique)

relation.habitats.sites.bent2.500m.unique<-aggregate(relation.habitats.sites.bent2.500m$surface,list(Code_Site=relation.habitats.sites.bent2.500m$Code_Site,habitat=relation.habitats.sites.bent2.500m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent2.500m.unique)=c("Code_Site","habitat","surface")
tablo_bent2_500m<-xtabs(formula=relation.habitats.sites.bent2.500m.unique$surface~relation.habitats.sites.bent2.500m.unique$Code_Site+relation.habitats.sites.bent2.500m.unique$habitat,relation.habitats.sites.bent2.500m.unique)

relation.habitats.sites.bent2.750m.unique<-aggregate(relation.habitats.sites.bent2.750m$surface,list(Code_Site=relation.habitats.sites.bent2.750m$Code_Site,habitat=relation.habitats.sites.bent2.750m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent2.750m.unique)=c("Code_Site","habitat","surface")
tablo_bent2_750m<-xtabs(formula=relation.habitats.sites.bent2.750m.unique$surface~relation.habitats.sites.bent2.750m.unique$Code_Site+relation.habitats.sites.bent2.750m.unique$habitat,relation.habitats.sites.bent2.750m.unique)

relation.habitats.sites.bent2.1000m.unique<-aggregate(relation.habitats.sites.bent2.1000m$surface,list(Code_Site=relation.habitats.sites.bent2.1000m$Code_Site,habitat=relation.habitats.sites.bent2.1000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent2.1000m.unique)=c("Code_Site","habitat","surface")
tablo_bent2_1000m<-xtabs(formula=relation.habitats.sites.bent2.1000m.unique$surface~relation.habitats.sites.bent2.1000m.unique$Code_Site+relation.habitats.sites.bent2.1000m.unique$habitat,relation.habitats.sites.bent2.1000m.unique)

relation.habitats.sites.bent2.3000m.unique<-aggregate(relation.habitats.sites.bent2.3000m$surface,list(Code_Site=relation.habitats.sites.bent2.3000m$Code_Site,habitat=relation.habitats.sites.bent2.3000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent2.3000m.unique)=c("Code_Site","habitat","surface")
tablo_bent2_3000m<-xtabs(formula=relation.habitats.sites.bent2.3000m.unique$surface~relation.habitats.sites.bent2.3000m.unique$Code_Site+relation.habitats.sites.bent2.3000m.unique$habitat,relation.habitats.sites.bent2.3000m.unique)

relation.habitats.sites.bent2.5000m.unique<-aggregate(relation.habitats.sites.bent2.5000m$surface,list(Code_Site=relation.habitats.sites.bent2.5000m$Code_Site,habitat=relation.habitats.sites.bent2.5000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent2.5000m.unique)=c("Code_Site","habitat","surface")
tablo_bent2_5000m<-xtabs(formula=relation.habitats.sites.bent2.5000m.unique$surface~relation.habitats.sites.bent2.5000m.unique$Code_Site+relation.habitats.sites.bent2.5000m.unique$habitat,relation.habitats.sites.bent2.5000m.unique)

tablo_bent2_30m<-tablo_bent2_30m[25:51,]
tablo_bent2_60m<-tablo_bent2_60m[25:51,]
tablo_bent2_100m<-tablo_bent2_100m[25:51,]
tablo_bent2_250m<-tablo_bent2_250m[25:51,]
tablo_bent2_500m<-tablo_bent2_500m[25:51,]
tablo_bent2_750m<-tablo_bent2_750m[25:51,]
tablo_bent2_1000m<-tablo_bent2_1000m[25:51,]
tablo_bent2_3000m<-tablo_bent2_3000m[25:51,]
tablo_bent2_5000m<-tablo_bent2_5000m[25:51,]

#Load de l'algorithme shannon
source("script.shannon.r")

#tirage des sites poissons en optimisant la diversité shannon d'habitat
 selection.shannon.hab.bent2.30m<-shannon(tablo_bent2_30m)
 selection.shannon.hab.bent2.60m<-shannon(tablo_bent2_60m)
 selection.shannon.hab.bent2.100m<-shannon(tablo_bent2_100m)
 selection.shannon.hab.bent2.250m<-shannon(tablo_bent2_250m)
 selection.shannon.hab.bent2.500m<-shannon(tablo_bent2_500m)
 selection.shannon.hab.bent2.750m<-shannon(tablo_bent2_750m)
 selection.shannon.hab.bent2.1000m<-shannon(tablo_bent2_1000m)
 selection.shannon.hab.bent2.3000m<-shannon(tablo_bent2_3000m)
 selection.shannon.hab.bent2.5000m<-shannon(tablo_bent2_5000m)

  #DIVERSITE SELECTIONNE PAR LE SCRIPT SHANNON
 source('script.diversite.shannon.selectionnee.fonction.r')
diversite.selection.shannon.hab.bent2.30m<-species.shannon.evol(selection.shannon.hab.bent2.30m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent2.60m<-species.shannon.evol(selection.shannon.hab.bent2.60m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent2.100m<-species.shannon.evol(selection.shannon.hab.bent2.100m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent2.250m<-species.shannon.evol(selection.shannon.hab.bent2.250m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent2.500m<-species.shannon.evol(selection.shannon.hab.bent2.500m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent2.750m<-species.shannon.evol(selection.shannon.hab.bent2.750m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent2.1000m<-species.shannon.evol(selection.shannon.hab.bent2.1000m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent2.3000m<-species.shannon.evol(selection.shannon.hab.bent2.3000m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent2.5000m<-species.shannon.evol(selection.shannon.hab.bent2.5000m,requette.abondance.poissons)

 #CALCUL DES SAI#
source('script.sai.r')
SAI_shannon.bent2.30m<-SAI(diversite.selection.shannon.hab.bent2.30m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent2.60m<-SAI(diversite.selection.shannon.hab.bent2.60m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent2.100m<-SAI(diversite.selection.shannon.hab.bent2.100m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent2.250m<-SAI(diversite.selection.shannon.hab.bent2.250m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent2.500m<-SAI(diversite.selection.shannon.hab.bent2.500m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent2.750m<-SAI(diversite.selection.shannon.hab.bent2.750m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent2.1000m<-SAI(diversite.selection.shannon.hab.bent2.1000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent2.3000m<-SAI(diversite.selection.shannon.hab.bent2.3000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent2.5000m<-SAI(diversite.selection.shannon.hab.bent2.5000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)


  #############Habitat bent3###################



   ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances

setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
relation.habitats.sites.bent3<-read.csv2("fc3_intersect_bent3.csv")

sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')

relation.habitats.sites.bent3<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$Code_Site%in%sites.etudies,]

#création d'une nouvelle table site-habitat par distance
relation.habitats.sites.bent3.30m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==30,]
relation.habitats.sites.bent3.60m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==60,]
relation.habitats.sites.bent3.100m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==100,]
relation.habitats.sites.bent3.250m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==250,]
relation.habitats.sites.bent3.500m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==500,]
relation.habitats.sites.bent3.750m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==750,]
relation.habitats.sites.bent3.1000m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==1000,]
relation.habitats.sites.bent3.3000m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==3000,]
relation.habitats.sites.bent3.5000m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==5000,]

#cration d'un tableau adapté à l'entrée du script shannon
relation.habitats.sites.bent3.30m.unique<-aggregate(relation.habitats.sites.bent3.30m$surface,list(Code_Site=relation.habitats.sites.bent3.30m$Code_Site,habitat=relation.habitats.sites.bent3.30m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent3.30m.unique)=c("Code_Site","habitat","surface")
tablo_bent3_30m<-xtabs(formula=relation.habitats.sites.bent3.30m.unique$surface~relation.habitats.sites.bent3.30m.unique$Code_Site+relation.habitats.sites.bent3.30m.unique$habitat,relation.habitats.sites.bent3.30m.unique)

relation.habitats.sites.bent3.60m.unique<-aggregate(relation.habitats.sites.bent3.60m$surface,list(Code_Site=relation.habitats.sites.bent3.60m$Code_Site,habitat=relation.habitats.sites.bent3.60m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent3.60m.unique)=c("Code_Site","habitat","surface")
tablo_bent3_60m<-xtabs(formula=relation.habitats.sites.bent3.60m.unique$surface~relation.habitats.sites.bent3.60m.unique$Code_Site+relation.habitats.sites.bent3.60m.unique$habitat,relation.habitats.sites.bent3.60m.unique)

relation.habitats.sites.bent3.100m.unique<-aggregate(relation.habitats.sites.bent3.100m$surface,list(Code_Site=relation.habitats.sites.bent3.100m$Code_Site,habitat=relation.habitats.sites.bent3.100m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent3.100m.unique)=c("Code_Site","habitat","surface")
tablo_bent3_100m<-xtabs(formula=relation.habitats.sites.bent3.100m.unique$surface~relation.habitats.sites.bent3.100m.unique$Code_Site+relation.habitats.sites.bent3.100m.unique$habitat,relation.habitats.sites.bent3.100m.unique)

relation.habitats.sites.bent3.250m.unique<-aggregate(relation.habitats.sites.bent3.250m$surface,list(Code_Site=relation.habitats.sites.bent3.250m$Code_Site,habitat=relation.habitats.sites.bent3.250m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent3.250m.unique)=c("Code_Site","habitat","surface")
tablo_bent3_250m<-xtabs(formula=relation.habitats.sites.bent3.250m.unique$surface~relation.habitats.sites.bent3.250m.unique$Code_Site+relation.habitats.sites.bent3.250m.unique$habitat,relation.habitats.sites.bent3.250m.unique)

relation.habitats.sites.bent3.500m.unique<-aggregate(relation.habitats.sites.bent3.500m$surface,list(Code_Site=relation.habitats.sites.bent3.500m$Code_Site,habitat=relation.habitats.sites.bent3.500m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent3.500m.unique)=c("Code_Site","habitat","surface")
tablo_bent3_500m<-xtabs(formula=relation.habitats.sites.bent3.500m.unique$surface~relation.habitats.sites.bent3.500m.unique$Code_Site+relation.habitats.sites.bent3.500m.unique$habitat,relation.habitats.sites.bent3.500m.unique)

relation.habitats.sites.bent3.750m.unique<-aggregate(relation.habitats.sites.bent3.750m$surface,list(Code_Site=relation.habitats.sites.bent3.750m$Code_Site,habitat=relation.habitats.sites.bent3.750m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent3.750m.unique)=c("Code_Site","habitat","surface")
tablo_bent3_750m<-xtabs(formula=relation.habitats.sites.bent3.750m.unique$surface~relation.habitats.sites.bent3.750m.unique$Code_Site+relation.habitats.sites.bent3.750m.unique$habitat,relation.habitats.sites.bent3.750m.unique)

relation.habitats.sites.bent3.1000m.unique<-aggregate(relation.habitats.sites.bent3.1000m$surface,list(Code_Site=relation.habitats.sites.bent3.1000m$Code_Site,habitat=relation.habitats.sites.bent3.1000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent3.1000m.unique)=c("Code_Site","habitat","surface")
tablo_bent3_1000m<-xtabs(formula=relation.habitats.sites.bent3.1000m.unique$surface~relation.habitats.sites.bent3.1000m.unique$Code_Site+relation.habitats.sites.bent3.1000m.unique$habitat,relation.habitats.sites.bent3.1000m.unique)

relation.habitats.sites.bent3.3000m.unique<-aggregate(relation.habitats.sites.bent3.3000m$surface,list(Code_Site=relation.habitats.sites.bent3.3000m$Code_Site,habitat=relation.habitats.sites.bent3.3000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent3.3000m.unique)=c("Code_Site","habitat","surface")
tablo_bent3_3000m<-xtabs(formula=relation.habitats.sites.bent3.3000m.unique$surface~relation.habitats.sites.bent3.3000m.unique$Code_Site+relation.habitats.sites.bent3.3000m.unique$habitat,relation.habitats.sites.bent3.3000m.unique)

relation.habitats.sites.bent3.5000m.unique<-aggregate(relation.habitats.sites.bent3.5000m$surface,list(Code_Site=relation.habitats.sites.bent3.5000m$Code_Site,habitat=relation.habitats.sites.bent3.5000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent3.5000m.unique)=c("Code_Site","habitat","surface")
tablo_bent3_5000m<-xtabs(formula=relation.habitats.sites.bent3.5000m.unique$surface~relation.habitats.sites.bent3.5000m.unique$Code_Site+relation.habitats.sites.bent3.5000m.unique$habitat,relation.habitats.sites.bent3.5000m.unique)

tablo_bent3_30m<-tablo_bent3_30m[25:51,]
tablo_bent3_60m<-tablo_bent3_60m[25:51,]
tablo_bent3_100m<-tablo_bent3_100m[25:51,]
tablo_bent3_250m<-tablo_bent3_250m[25:51,]
tablo_bent3_500m<-tablo_bent3_500m[25:51,]
tablo_bent3_750m<-tablo_bent3_750m[25:51,]
tablo_bent3_1000m<-tablo_bent3_1000m[25:51,]
tablo_bent3_3000m<-tablo_bent3_3000m[25:51,]
tablo_bent3_5000m<-tablo_bent3_5000m[25:51,]

#Load de l'algorithme shannon
source("script.shannon.r")

#tirage des sites poissons en optimisant la diversité shannon d'habitat
 selection.shannon.hab.bent3.30m<-shannon(tablo_bent3_30m)
 selection.shannon.hab.bent3.60m<-shannon(tablo_bent3_60m)
 selection.shannon.hab.bent3.100m<-shannon(tablo_bent3_100m)
 selection.shannon.hab.bent3.250m<-shannon(tablo_bent3_250m)
 selection.shannon.hab.bent3.500m<-shannon(tablo_bent3_500m)
 selection.shannon.hab.bent3.750m<-shannon(tablo_bent3_750m)
 selection.shannon.hab.bent3.1000m<-shannon(tablo_bent3_1000m)
 selection.shannon.hab.bent3.3000m<-shannon(tablo_bent3_3000m)
 selection.shannon.hab.bent3.5000m<-shannon(tablo_bent3_5000m)

  #DIVERSITE SELECTIONNE PAR LE SCRIPT SHANNON
  source('script.diversite.shannon.selectionnee.fonction.r')
diversite.selection.shannon.hab.bent3.30m<-species.shannon.evol(selection.shannon.hab.bent3.30m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent3.60m<-species.shannon.evol(selection.shannon.hab.bent3.60m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent3.100m<-species.shannon.evol(selection.shannon.hab.bent3.100m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent3.250m<-species.shannon.evol(selection.shannon.hab.bent3.250m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent3.500m<-species.shannon.evol(selection.shannon.hab.bent3.500m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent3.750m<-species.shannon.evol(selection.shannon.hab.bent3.750m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent3.1000m<-species.shannon.evol(selection.shannon.hab.bent3.1000m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent3.3000m<-species.shannon.evol(selection.shannon.hab.bent3.3000m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent3.5000m<-species.shannon.evol(selection.shannon.hab.bent3.5000m,requette.abondance.poissons)

 #CALCUL DES SAI#
source('script.sai.r')
SAI_shannon.bent3.30m<-SAI(diversite.selection.shannon.hab.bent3.30m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent3.60m<-SAI(diversite.selection.shannon.hab.bent3.60m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent3.100m<-SAI(diversite.selection.shannon.hab.bent3.100m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent3.250m<-SAI(diversite.selection.shannon.hab.bent3.250m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent3.500m<-SAI(diversite.selection.shannon.hab.bent3.500m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent3.750m<-SAI(diversite.selection.shannon.hab.bent3.750m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent3.1000m<-SAI(diversite.selection.shannon.hab.bent3.1000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent3.3000m<-SAI(diversite.selection.shannon.hab.bent3.3000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent3.5000m<-SAI(diversite.selection.shannon.hab.bent3.5000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)


  #############Habitat geo1_geo2###################



   ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances

setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
relation.habitats.sites.geo1_geo2<-read.csv2("fc3_intersect_geo1_geo2.csv")

sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')

relation.habitats.sites.geo1_geo2<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$Code_Site%in%sites.etudies,]

#création d'une nouvelle table site-habitat par distance
relation.habitats.sites.geo1_geo2.30m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==30,]
relation.habitats.sites.geo1_geo2.60m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==60,]
relation.habitats.sites.geo1_geo2.100m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==100,]
relation.habitats.sites.geo1_geo2.250m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==250,]
relation.habitats.sites.geo1_geo2.500m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==500,]
relation.habitats.sites.geo1_geo2.750m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==750,]
relation.habitats.sites.geo1_geo2.1000m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==1000,]
relation.habitats.sites.geo1_geo2.3000m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==3000,]
relation.habitats.sites.geo1_geo2.5000m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==5000,]

#cration d'un tableau adapté à l'entrée du script shannon
relation.habitats.sites.geo1_geo2.30m.unique<-aggregate(relation.habitats.sites.geo1_geo2.30m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2.30m$Code_Site,habitat=relation.habitats.sites.geo1_geo2.30m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2.30m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_30m<-xtabs(formula=relation.habitats.sites.geo1_geo2.30m.unique$surface~relation.habitats.sites.geo1_geo2.30m.unique$Code_Site+relation.habitats.sites.geo1_geo2.30m.unique$habitat,relation.habitats.sites.geo1_geo2.30m.unique)

relation.habitats.sites.geo1_geo2.60m.unique<-aggregate(relation.habitats.sites.geo1_geo2.60m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2.60m$Code_Site,habitat=relation.habitats.sites.geo1_geo2.60m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2.60m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_60m<-xtabs(formula=relation.habitats.sites.geo1_geo2.60m.unique$surface~relation.habitats.sites.geo1_geo2.60m.unique$Code_Site+relation.habitats.sites.geo1_geo2.60m.unique$habitat,relation.habitats.sites.geo1_geo2.60m.unique)

relation.habitats.sites.geo1_geo2.100m.unique<-aggregate(relation.habitats.sites.geo1_geo2.100m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2.100m$Code_Site,habitat=relation.habitats.sites.geo1_geo2.100m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2.100m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_100m<-xtabs(formula=relation.habitats.sites.geo1_geo2.100m.unique$surface~relation.habitats.sites.geo1_geo2.100m.unique$Code_Site+relation.habitats.sites.geo1_geo2.100m.unique$habitat,relation.habitats.sites.geo1_geo2.100m.unique)

relation.habitats.sites.geo1_geo2.250m.unique<-aggregate(relation.habitats.sites.geo1_geo2.250m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2.250m$Code_Site,habitat=relation.habitats.sites.geo1_geo2.250m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2.250m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_250m<-xtabs(formula=relation.habitats.sites.geo1_geo2.250m.unique$surface~relation.habitats.sites.geo1_geo2.250m.unique$Code_Site+relation.habitats.sites.geo1_geo2.250m.unique$habitat,relation.habitats.sites.geo1_geo2.250m.unique)

relation.habitats.sites.geo1_geo2.500m.unique<-aggregate(relation.habitats.sites.geo1_geo2.500m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2.500m$Code_Site,habitat=relation.habitats.sites.geo1_geo2.500m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2.500m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_500m<-xtabs(formula=relation.habitats.sites.geo1_geo2.500m.unique$surface~relation.habitats.sites.geo1_geo2.500m.unique$Code_Site+relation.habitats.sites.geo1_geo2.500m.unique$habitat,relation.habitats.sites.geo1_geo2.500m.unique)

relation.habitats.sites.geo1_geo2.750m.unique<-aggregate(relation.habitats.sites.geo1_geo2.750m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2.750m$Code_Site,habitat=relation.habitats.sites.geo1_geo2.750m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2.750m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_750m<-xtabs(formula=relation.habitats.sites.geo1_geo2.750m.unique$surface~relation.habitats.sites.geo1_geo2.750m.unique$Code_Site+relation.habitats.sites.geo1_geo2.750m.unique$habitat,relation.habitats.sites.geo1_geo2.750m.unique)

relation.habitats.sites.geo1_geo2.1000m.unique<-aggregate(relation.habitats.sites.geo1_geo2.1000m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2.1000m$Code_Site,habitat=relation.habitats.sites.geo1_geo2.1000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2.1000m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_1000m<-xtabs(formula=relation.habitats.sites.geo1_geo2.1000m.unique$surface~relation.habitats.sites.geo1_geo2.1000m.unique$Code_Site+relation.habitats.sites.geo1_geo2.1000m.unique$habitat,relation.habitats.sites.geo1_geo2.1000m.unique)

relation.habitats.sites.geo1_geo2.3000m.unique<-aggregate(relation.habitats.sites.geo1_geo2.3000m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2.3000m$Code_Site,habitat=relation.habitats.sites.geo1_geo2.3000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2.3000m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_3000m<-xtabs(formula=relation.habitats.sites.geo1_geo2.3000m.unique$surface~relation.habitats.sites.geo1_geo2.3000m.unique$Code_Site+relation.habitats.sites.geo1_geo2.3000m.unique$habitat,relation.habitats.sites.geo1_geo2.3000m.unique)

relation.habitats.sites.geo1_geo2.5000m.unique<-aggregate(relation.habitats.sites.geo1_geo2.5000m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2.5000m$Code_Site,habitat=relation.habitats.sites.geo1_geo2.5000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2.5000m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_5000m<-xtabs(formula=relation.habitats.sites.geo1_geo2.5000m.unique$surface~relation.habitats.sites.geo1_geo2.5000m.unique$Code_Site+relation.habitats.sites.geo1_geo2.5000m.unique$habitat,relation.habitats.sites.geo1_geo2.5000m.unique)

tablo_geo1_geo2_30m<-tablo_geo1_geo2_30m[25:51,]
tablo_geo1_geo2_60m<-tablo_geo1_geo2_60m[25:51,]
tablo_geo1_geo2_100m<-tablo_geo1_geo2_100m[25:51,]
tablo_geo1_geo2_250m<-tablo_geo1_geo2_250m[25:51,]
tablo_geo1_geo2_500m<-tablo_geo1_geo2_500m[25:51,]
tablo_geo1_geo2_750m<-tablo_geo1_geo2_750m[25:51,]
tablo_geo1_geo2_1000m<-tablo_geo1_geo2_1000m[25:51,]
tablo_geo1_geo2_3000m<-tablo_geo1_geo2_3000m[25:51,]
tablo_geo1_geo2_5000m<-tablo_geo1_geo2_5000m[25:51,]

#Load de l'algorithme shannon
source("script.shannon.r")

#tirage des sites poissons en optimisant la diversité shannon d'habitat
 selection.shannon.hab.geo1_geo2.30m<-shannon(tablo_geo1_geo2_30m)
 selection.shannon.hab.geo1_geo2.60m<-shannon(tablo_geo1_geo2_60m)
 selection.shannon.hab.geo1_geo2.100m<-shannon(tablo_geo1_geo2_100m)
 selection.shannon.hab.geo1_geo2.250m<-shannon(tablo_geo1_geo2_250m)
 selection.shannon.hab.geo1_geo2.500m<-shannon(tablo_geo1_geo2_500m)
 selection.shannon.hab.geo1_geo2.750m<-shannon(tablo_geo1_geo2_750m)
 selection.shannon.hab.geo1_geo2.1000m<-shannon(tablo_geo1_geo2_1000m)
 selection.shannon.hab.geo1_geo2.3000m<-shannon(tablo_geo1_geo2_3000m)
 selection.shannon.hab.geo1_geo2.5000m<-shannon(tablo_geo1_geo2_5000m)

  #DIVERSITE SELECTIONNE PAR LE SCRIPT SHANNON
  source('script.diversite.shannon.selectionnee.fonction.r')
diversite.selection.shannon.hab.geo1_geo2.30m<-species.shannon.evol(selection.shannon.hab.geo1_geo2.30m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1_geo2.60m<-species.shannon.evol(selection.shannon.hab.geo1_geo2.60m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1_geo2.100m<-species.shannon.evol(selection.shannon.hab.geo1_geo2.100m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1_geo2.250m<-species.shannon.evol(selection.shannon.hab.geo1_geo2.250m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1_geo2.500m<-species.shannon.evol(selection.shannon.hab.geo1_geo2.500m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1_geo2.750m<-species.shannon.evol(selection.shannon.hab.geo1_geo2.750m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1_geo2.1000m<-species.shannon.evol(selection.shannon.hab.geo1_geo2.1000m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1_geo2.3000m<-species.shannon.evol(selection.shannon.hab.geo1_geo2.3000m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1_geo2.5000m<-species.shannon.evol(selection.shannon.hab.geo1_geo2.5000m,requette.abondance.poissons)

 #CALCUL DES SAI#
source('script.sai.r')
SAI_shannon.geo1_geo2.30m<-SAI(diversite.selection.shannon.hab.geo1_geo2.30m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1_geo2.60m<-SAI(diversite.selection.shannon.hab.geo1_geo2.60m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1_geo2.100m<-SAI(diversite.selection.shannon.hab.geo1_geo2.100m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1_geo2.250m<-SAI(diversite.selection.shannon.hab.geo1_geo2.250m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1_geo2.500m<-SAI(diversite.selection.shannon.hab.geo1_geo2.500m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1_geo2.750m<-SAI(diversite.selection.shannon.hab.geo1_geo2.750m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1_geo2.1000m<-SAI(diversite.selection.shannon.hab.geo1_geo2.1000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1_geo2.3000m<-SAI(diversite.selection.shannon.hab.geo1_geo2.3000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1_geo2.5000m<-SAI(diversite.selection.shannon.hab.geo1_geo2.5000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)


  #############Habitat geo1_geo2_geo3###################



   ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances

setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
relation.habitats.sites.geo1_geo2_geo3<-read.csv2("fc3_intersect_geo1_geo2_geo3.csv")

sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')

relation.habitats.sites.geo1_geo2_geo3<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$Code_Site%in%sites.etudies,]

#création d'une nouvelle table site-habitat par distance
relation.habitats.sites.geo1_geo2_geo3.30m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==30,]
relation.habitats.sites.geo1_geo2_geo3.60m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==60,]
relation.habitats.sites.geo1_geo2_geo3.100m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==100,]
relation.habitats.sites.geo1_geo2_geo3.250m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==250,]
relation.habitats.sites.geo1_geo2_geo3.500m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==500,]
relation.habitats.sites.geo1_geo2_geo3.750m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==750,]
relation.habitats.sites.geo1_geo2_geo3.1000m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==1000,]
relation.habitats.sites.geo1_geo2_geo3.3000m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==3000,]
relation.habitats.sites.geo1_geo2_geo3.5000m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==5000,]

#cration d'un tableau adapté à l'entrée du script shannon
relation.habitats.sites.geo1_geo2_geo3.30m.unique<-aggregate(relation.habitats.sites.geo1_geo2_geo3.30m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2_geo3.30m$Code_Site,habitat=relation.habitats.sites.geo1_geo2_geo3.30m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2_geo3.30m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_geo3_30m<-xtabs(formula=relation.habitats.sites.geo1_geo2_geo3.30m.unique$surface~relation.habitats.sites.geo1_geo2_geo3.30m.unique$Code_Site+relation.habitats.sites.geo1_geo2_geo3.30m.unique$habitat,relation.habitats.sites.geo1_geo2_geo3.30m.unique)

relation.habitats.sites.geo1_geo2_geo3.60m.unique<-aggregate(relation.habitats.sites.geo1_geo2_geo3.60m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2_geo3.60m$Code_Site,habitat=relation.habitats.sites.geo1_geo2_geo3.60m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2_geo3.60m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_geo3_60m<-xtabs(formula=relation.habitats.sites.geo1_geo2_geo3.60m.unique$surface~relation.habitats.sites.geo1_geo2_geo3.60m.unique$Code_Site+relation.habitats.sites.geo1_geo2_geo3.60m.unique$habitat,relation.habitats.sites.geo1_geo2_geo3.60m.unique)

relation.habitats.sites.geo1_geo2_geo3.100m.unique<-aggregate(relation.habitats.sites.geo1_geo2_geo3.100m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2_geo3.100m$Code_Site,habitat=relation.habitats.sites.geo1_geo2_geo3.100m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2_geo3.100m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_geo3_100m<-xtabs(formula=relation.habitats.sites.geo1_geo2_geo3.100m.unique$surface~relation.habitats.sites.geo1_geo2_geo3.100m.unique$Code_Site+relation.habitats.sites.geo1_geo2_geo3.100m.unique$habitat,relation.habitats.sites.geo1_geo2_geo3.100m.unique)

relation.habitats.sites.geo1_geo2_geo3.250m.unique<-aggregate(relation.habitats.sites.geo1_geo2_geo3.250m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2_geo3.250m$Code_Site,habitat=relation.habitats.sites.geo1_geo2_geo3.250m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2_geo3.250m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_geo3_250m<-xtabs(formula=relation.habitats.sites.geo1_geo2_geo3.250m.unique$surface~relation.habitats.sites.geo1_geo2_geo3.250m.unique$Code_Site+relation.habitats.sites.geo1_geo2_geo3.250m.unique$habitat,relation.habitats.sites.geo1_geo2_geo3.250m.unique)

relation.habitats.sites.geo1_geo2_geo3.500m.unique<-aggregate(relation.habitats.sites.geo1_geo2_geo3.500m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2_geo3.500m$Code_Site,habitat=relation.habitats.sites.geo1_geo2_geo3.500m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2_geo3.500m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_geo3_500m<-xtabs(formula=relation.habitats.sites.geo1_geo2_geo3.500m.unique$surface~relation.habitats.sites.geo1_geo2_geo3.500m.unique$Code_Site+relation.habitats.sites.geo1_geo2_geo3.500m.unique$habitat,relation.habitats.sites.geo1_geo2_geo3.500m.unique)

relation.habitats.sites.geo1_geo2_geo3.750m.unique<-aggregate(relation.habitats.sites.geo1_geo2_geo3.750m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2_geo3.750m$Code_Site,habitat=relation.habitats.sites.geo1_geo2_geo3.750m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2_geo3.750m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_geo3_750m<-xtabs(formula=relation.habitats.sites.geo1_geo2_geo3.750m.unique$surface~relation.habitats.sites.geo1_geo2_geo3.750m.unique$Code_Site+relation.habitats.sites.geo1_geo2_geo3.750m.unique$habitat,relation.habitats.sites.geo1_geo2_geo3.750m.unique)

relation.habitats.sites.geo1_geo2_geo3.1000m.unique<-aggregate(relation.habitats.sites.geo1_geo2_geo3.1000m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2_geo3.1000m$Code_Site,habitat=relation.habitats.sites.geo1_geo2_geo3.1000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2_geo3.1000m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_geo3_1000m<-xtabs(formula=relation.habitats.sites.geo1_geo2_geo3.1000m.unique$surface~relation.habitats.sites.geo1_geo2_geo3.1000m.unique$Code_Site+relation.habitats.sites.geo1_geo2_geo3.1000m.unique$habitat,relation.habitats.sites.geo1_geo2_geo3.1000m.unique)

relation.habitats.sites.geo1_geo2_geo3.3000m.unique<-aggregate(relation.habitats.sites.geo1_geo2_geo3.3000m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2_geo3.3000m$Code_Site,habitat=relation.habitats.sites.geo1_geo2_geo3.3000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2_geo3.3000m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_geo3_3000m<-xtabs(formula=relation.habitats.sites.geo1_geo2_geo3.3000m.unique$surface~relation.habitats.sites.geo1_geo2_geo3.3000m.unique$Code_Site+relation.habitats.sites.geo1_geo2_geo3.3000m.unique$habitat,relation.habitats.sites.geo1_geo2_geo3.3000m.unique)

relation.habitats.sites.geo1_geo2_geo3.5000m.unique<-aggregate(relation.habitats.sites.geo1_geo2_geo3.5000m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2_geo3.5000m$Code_Site,habitat=relation.habitats.sites.geo1_geo2_geo3.5000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2_geo3.5000m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_geo3_5000m<-xtabs(formula=relation.habitats.sites.geo1_geo2_geo3.5000m.unique$surface~relation.habitats.sites.geo1_geo2_geo3.5000m.unique$Code_Site+relation.habitats.sites.geo1_geo2_geo3.5000m.unique$habitat,relation.habitats.sites.geo1_geo2_geo3.5000m.unique)

tablo_geo1_geo2_geo3_30m<-tablo_geo1_geo2_geo3_30m[25:51,]
tablo_geo1_geo2_geo3_60m<-tablo_geo1_geo2_geo3_60m[25:51,]
tablo_geo1_geo2_geo3_100m<-tablo_geo1_geo2_geo3_100m[25:51,]
tablo_geo1_geo2_geo3_250m<-tablo_geo1_geo2_geo3_250m[25:51,]
tablo_geo1_geo2_geo3_500m<-tablo_geo1_geo2_geo3_500m[25:51,]
tablo_geo1_geo2_geo3_750m<-tablo_geo1_geo2_geo3_750m[25:51,]
tablo_geo1_geo2_geo3_1000m<-tablo_geo1_geo2_geo3_1000m[25:51,]
tablo_geo1_geo2_geo3_3000m<-tablo_geo1_geo2_geo3_3000m[25:51,]
tablo_geo1_geo2_geo3_5000m<-tablo_geo1_geo2_geo3_5000m[25:51,]

#Load de l'algorithme shannon
source("script.shannon.r")

#tirage des sites poissons en optimisant la diversité shannon d'habitat
 selection.shannon.hab.geo1_geo2_geo3.30m<-shannon(tablo_geo1_geo2_geo3_30m)
 selection.shannon.hab.geo1_geo2_geo3.60m<-shannon(tablo_geo1_geo2_geo3_60m)
 selection.shannon.hab.geo1_geo2_geo3.100m<-shannon(tablo_geo1_geo2_geo3_100m)
 selection.shannon.hab.geo1_geo2_geo3.250m<-shannon(tablo_geo1_geo2_geo3_250m)
 selection.shannon.hab.geo1_geo2_geo3.500m<-shannon(tablo_geo1_geo2_geo3_500m)
 selection.shannon.hab.geo1_geo2_geo3.750m<-shannon(tablo_geo1_geo2_geo3_750m)
 selection.shannon.hab.geo1_geo2_geo3.1000m<-shannon(tablo_geo1_geo2_geo3_1000m)
 selection.shannon.hab.geo1_geo2_geo3.3000m<-shannon(tablo_geo1_geo2_geo3_3000m)
 selection.shannon.hab.geo1_geo2_geo3.5000m<-shannon(tablo_geo1_geo2_geo3_5000m)

  #DIVERSITE SELECTIONNE PAR LE SCRIPT SHANNON
  source('script.diversite.shannon.selectionnee.fonction.r')
diversite.selection.shannon.hab.geo1_geo2_geo3.30m<-species.shannon.evol(selection.shannon.hab.geo1_geo2_geo3.30m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1_geo2_geo3.60m<-species.shannon.evol(selection.shannon.hab.geo1_geo2_geo3.60m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1_geo2_geo3.100m<-species.shannon.evol(selection.shannon.hab.geo1_geo2_geo3.100m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1_geo2_geo3.250m<-species.shannon.evol(selection.shannon.hab.geo1_geo2_geo3.250m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1_geo2_geo3.500m<-species.shannon.evol(selection.shannon.hab.geo1_geo2_geo3.500m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1_geo2_geo3.750m<-species.shannon.evol(selection.shannon.hab.geo1_geo2_geo3.750m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1_geo2_geo3.1000m<-species.shannon.evol(selection.shannon.hab.geo1_geo2_geo3.1000m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1_geo2_geo3.3000m<-species.shannon.evol(selection.shannon.hab.geo1_geo2_geo3.3000m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1_geo2_geo3.5000m<-species.shannon.evol(selection.shannon.hab.geo1_geo2_geo3.5000m,requette.abondance.poissons)

 #CALCUL DES SAI#
source('script.sai.r')
SAI_shannon.geo1_geo2_geo3.30m<-SAI(diversite.selection.shannon.hab.geo1_geo2_geo3.30m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1_geo2_geo3.60m<-SAI(diversite.selection.shannon.hab.geo1_geo2_geo3.60m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1_geo2_geo3.100m<-SAI(diversite.selection.shannon.hab.geo1_geo2_geo3.100m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1_geo2_geo3.250m<-SAI(diversite.selection.shannon.hab.geo1_geo2_geo3.250m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1_geo2_geo3.500m<-SAI(diversite.selection.shannon.hab.geo1_geo2_geo3.500m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1_geo2_geo3.750m<-SAI(diversite.selection.shannon.hab.geo1_geo2_geo3.750m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1_geo2_geo3.1000m<-SAI(diversite.selection.shannon.hab.geo1_geo2_geo3.1000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1_geo2_geo3.3000m<-SAI(diversite.selection.shannon.hab.geo1_geo2_geo3.3000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1_geo2_geo3.5000m<-SAI(diversite.selection.shannon.hab.geo1_geo2_geo3.5000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)


  #############Habitat g1_g2_g3_b1###################



   ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances

setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
relation.habitats.sites.g1_g2_g3_b1<-read.csv2("fc3_intersect_g1_g2_g3_b1.csv")

sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')

relation.habitats.sites.g1_g2_g3_b1<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$Code_Site%in%sites.etudies,]

#création d'une nouvelle table site-habitat par distance
relation.habitats.sites.g1_g2_g3_b1.30m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==30,]
relation.habitats.sites.g1_g2_g3_b1.60m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==60,]
relation.habitats.sites.g1_g2_g3_b1.100m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==100,]
relation.habitats.sites.g1_g2_g3_b1.250m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==250,]
relation.habitats.sites.g1_g2_g3_b1.500m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==500,]
relation.habitats.sites.g1_g2_g3_b1.750m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==750,]
relation.habitats.sites.g1_g2_g3_b1.1000m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==1000,]
relation.habitats.sites.g1_g2_g3_b1.3000m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==3000,]
relation.habitats.sites.g1_g2_g3_b1.5000m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==5000,]

#cration d'un tableau adapté à l'entrée du script shannon
relation.habitats.sites.g1_g2_g3_b1.30m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1.30m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1.30m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1.30m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1.30m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_30m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1.30m.unique$surface~relation.habitats.sites.g1_g2_g3_b1.30m.unique$Code_Site+relation.habitats.sites.g1_g2_g3_b1.30m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1.30m.unique)

relation.habitats.sites.g1_g2_g3_b1.60m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1.60m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1.60m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1.60m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1.60m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_60m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1.60m.unique$surface~relation.habitats.sites.g1_g2_g3_b1.60m.unique$Code_Site+relation.habitats.sites.g1_g2_g3_b1.60m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1.60m.unique)

relation.habitats.sites.g1_g2_g3_b1.100m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1.100m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1.100m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1.100m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1.100m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_100m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1.100m.unique$surface~relation.habitats.sites.g1_g2_g3_b1.100m.unique$Code_Site+relation.habitats.sites.g1_g2_g3_b1.100m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1.100m.unique)

relation.habitats.sites.g1_g2_g3_b1.250m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1.250m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1.250m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1.250m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1.250m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_250m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1.250m.unique$surface~relation.habitats.sites.g1_g2_g3_b1.250m.unique$Code_Site+relation.habitats.sites.g1_g2_g3_b1.250m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1.250m.unique)

relation.habitats.sites.g1_g2_g3_b1.500m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1.500m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1.500m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1.500m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1.500m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_500m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1.500m.unique$surface~relation.habitats.sites.g1_g2_g3_b1.500m.unique$Code_Site+relation.habitats.sites.g1_g2_g3_b1.500m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1.500m.unique)

relation.habitats.sites.g1_g2_g3_b1.750m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1.750m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1.750m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1.750m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1.750m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_750m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1.750m.unique$surface~relation.habitats.sites.g1_g2_g3_b1.750m.unique$Code_Site+relation.habitats.sites.g1_g2_g3_b1.750m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1.750m.unique)

relation.habitats.sites.g1_g2_g3_b1.1000m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1.1000m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1.1000m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1.1000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1.1000m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_1000m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1.1000m.unique$surface~relation.habitats.sites.g1_g2_g3_b1.1000m.unique$Code_Site+relation.habitats.sites.g1_g2_g3_b1.1000m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1.1000m.unique)

relation.habitats.sites.g1_g2_g3_b1.3000m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1.3000m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1.3000m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1.3000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1.3000m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_3000m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1.3000m.unique$surface~relation.habitats.sites.g1_g2_g3_b1.3000m.unique$Code_Site+relation.habitats.sites.g1_g2_g3_b1.3000m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1.3000m.unique)

relation.habitats.sites.g1_g2_g3_b1.5000m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1.5000m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1.5000m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1.5000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1.5000m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_5000m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1.5000m.unique$surface~relation.habitats.sites.g1_g2_g3_b1.5000m.unique$Code_Site+relation.habitats.sites.g1_g2_g3_b1.5000m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1.5000m.unique)

tablo_g1_g2_g3_b1_30m<-tablo_g1_g2_g3_b1_30m[25:51,]
tablo_g1_g2_g3_b1_60m<-tablo_g1_g2_g3_b1_60m[25:51,]
tablo_g1_g2_g3_b1_100m<-tablo_g1_g2_g3_b1_100m[25:51,]
tablo_g1_g2_g3_b1_250m<-tablo_g1_g2_g3_b1_250m[25:51,]
tablo_g1_g2_g3_b1_500m<-tablo_g1_g2_g3_b1_500m[25:51,]
tablo_g1_g2_g3_b1_750m<-tablo_g1_g2_g3_b1_750m[25:51,]
tablo_g1_g2_g3_b1_1000m<-tablo_g1_g2_g3_b1_1000m[25:51,]
tablo_g1_g2_g3_b1_3000m<-tablo_g1_g2_g3_b1_3000m[25:51,]
tablo_g1_g2_g3_b1_5000m<-tablo_g1_g2_g3_b1_5000m[25:51,]

#Load de l'algorithme shannon
source("script.shannon.r")

#tirage des sites poissons en optimisant la diversité shannon d'habitat
 selection.shannon.hab.g1_g2_g3_b1.30m<-shannon(tablo_g1_g2_g3_b1_30m)
 selection.shannon.hab.g1_g2_g3_b1.60m<-shannon(tablo_g1_g2_g3_b1_60m)
 selection.shannon.hab.g1_g2_g3_b1.100m<-shannon(tablo_g1_g2_g3_b1_100m)
 selection.shannon.hab.g1_g2_g3_b1.250m<-shannon(tablo_g1_g2_g3_b1_250m)
 selection.shannon.hab.g1_g2_g3_b1.500m<-shannon(tablo_g1_g2_g3_b1_500m)
 selection.shannon.hab.g1_g2_g3_b1.750m<-shannon(tablo_g1_g2_g3_b1_750m)
 selection.shannon.hab.g1_g2_g3_b1.1000m<-shannon(tablo_g1_g2_g3_b1_1000m)
 selection.shannon.hab.g1_g2_g3_b1.3000m<-shannon(tablo_g1_g2_g3_b1_3000m)
 selection.shannon.hab.g1_g2_g3_b1.5000m<-shannon(tablo_g1_g2_g3_b1_5000m)

  #DIVERSITE SELECTIONNE PAR LE SCRIPT SHANNON
  source('script.diversite.shannon.selectionnee.fonction.r')
diversite.selection.shannon.hab.g1_g2_g3_b1.30m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1.30m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1.60m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1.60m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1.100m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1.100m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1.250m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1.250m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1.500m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1.500m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1.750m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1.750m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1.1000m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1.1000m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1.3000m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1.3000m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1.5000m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1.5000m,requette.abondance.poissons)

 #CALCUL DES SAI#
source('script.sai.r')
SAI_shannon.g1_g2_g3_b1.30m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1.30m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1.60m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1.60m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1.100m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1.100m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1.250m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1.250m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1.500m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1.500m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1.750m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1.750m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1.1000m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1.1000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1.3000m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1.3000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1.5000m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1.5000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)


  #############Habitat g1_g2_g3_b1_b2###################



   ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances

setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
relation.habitats.sites.g1_g2_g3_b1_b2<-read.csv2("fc3_intersect_g1_g2_g3_b1_b2.csv")

sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')

relation.habitats.sites.g1_g2_g3_b1_b2<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$Code_Site%in%sites.etudies,]

#création d'une nouvelle table site-habitat par distance
relation.habitats.sites.g1_g2_g3_b1_b2.30m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==30,]
relation.habitats.sites.g1_g2_g3_b1_b2.60m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==60,]
relation.habitats.sites.g1_g2_g3_b1_b2.100m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==100,]
relation.habitats.sites.g1_g2_g3_b1_b2.250m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==250,]
relation.habitats.sites.g1_g2_g3_b1_b2.500m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==500,]
relation.habitats.sites.g1_g2_g3_b1_b2.750m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==750,]
relation.habitats.sites.g1_g2_g3_b1_b2.1000m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==1000,]
relation.habitats.sites.g1_g2_g3_b1_b2.3000m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==3000,]
relation.habitats.sites.g1_g2_g3_b1_b2.5000m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==5000,]

#cration d'un tableau adapté à l'entrée du script shannon
relation.habitats.sites.g1_g2_g3_b1_b2.30m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2.30m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2.30m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2.30m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2.30m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_30m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2.30m.unique$surface~relation.habitats.sites.g1_g2_g3_b1_b2.30m.unique$Code_Site+relation.habitats.sites.g1_g2_g3_b1_b2.30m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2.30m.unique)

relation.habitats.sites.g1_g2_g3_b1_b2.60m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2.60m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2.60m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2.60m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2.60m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_60m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2.60m.unique$surface~relation.habitats.sites.g1_g2_g3_b1_b2.60m.unique$Code_Site+relation.habitats.sites.g1_g2_g3_b1_b2.60m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2.60m.unique)

relation.habitats.sites.g1_g2_g3_b1_b2.100m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2.100m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2.100m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2.100m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2.100m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_100m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2.100m.unique$surface~relation.habitats.sites.g1_g2_g3_b1_b2.100m.unique$Code_Site+relation.habitats.sites.g1_g2_g3_b1_b2.100m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2.100m.unique)

relation.habitats.sites.g1_g2_g3_b1_b2.250m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2.250m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2.250m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2.250m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2.250m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_250m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2.250m.unique$surface~relation.habitats.sites.g1_g2_g3_b1_b2.250m.unique$Code_Site+relation.habitats.sites.g1_g2_g3_b1_b2.250m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2.250m.unique)

relation.habitats.sites.g1_g2_g3_b1_b2.500m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2.500m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2.500m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2.500m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2.500m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_500m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2.500m.unique$surface~relation.habitats.sites.g1_g2_g3_b1_b2.500m.unique$Code_Site+relation.habitats.sites.g1_g2_g3_b1_b2.500m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2.500m.unique)

relation.habitats.sites.g1_g2_g3_b1_b2.750m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2.750m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2.750m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2.750m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2.750m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_750m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2.750m.unique$surface~relation.habitats.sites.g1_g2_g3_b1_b2.750m.unique$Code_Site+relation.habitats.sites.g1_g2_g3_b1_b2.750m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2.750m.unique)

relation.habitats.sites.g1_g2_g3_b1_b2.1000m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2.1000m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2.1000m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2.1000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2.1000m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_1000m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2.1000m.unique$surface~relation.habitats.sites.g1_g2_g3_b1_b2.1000m.unique$Code_Site+relation.habitats.sites.g1_g2_g3_b1_b2.1000m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2.1000m.unique)

relation.habitats.sites.g1_g2_g3_b1_b2.3000m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2.3000m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2.3000m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2.3000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2.3000m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_3000m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2.3000m.unique$surface~relation.habitats.sites.g1_g2_g3_b1_b2.3000m.unique$Code_Site+relation.habitats.sites.g1_g2_g3_b1_b2.3000m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2.3000m.unique)

relation.habitats.sites.g1_g2_g3_b1_b2.5000m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2.5000m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2.5000m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2.5000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2.5000m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_5000m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2.5000m.unique$surface~relation.habitats.sites.g1_g2_g3_b1_b2.5000m.unique$Code_Site+relation.habitats.sites.g1_g2_g3_b1_b2.5000m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2.5000m.unique)

tablo_g1_g2_g3_b1_b2_30m<-tablo_g1_g2_g3_b1_b2_30m[25:51,]
tablo_g1_g2_g3_b1_b2_60m<-tablo_g1_g2_g3_b1_b2_60m[25:51,]
tablo_g1_g2_g3_b1_b2_100m<-tablo_g1_g2_g3_b1_b2_100m[25:51,]
tablo_g1_g2_g3_b1_b2_250m<-tablo_g1_g2_g3_b1_b2_250m[25:51,]
tablo_g1_g2_g3_b1_b2_500m<-tablo_g1_g2_g3_b1_b2_500m[25:51,]
tablo_g1_g2_g3_b1_b2_750m<-tablo_g1_g2_g3_b1_b2_750m[25:51,]
tablo_g1_g2_g3_b1_b2_1000m<-tablo_g1_g2_g3_b1_b2_1000m[25:51,]
tablo_g1_g2_g3_b1_b2_3000m<-tablo_g1_g2_g3_b1_b2_3000m[25:51,]
tablo_g1_g2_g3_b1_b2_5000m<-tablo_g1_g2_g3_b1_b2_5000m[25:51,]

#Load de l'algorithme shannon
source("script.shannon.r")

#tirage des sites poissons en optimisant la diversité shannon d'habitat
 selection.shannon.hab.g1_g2_g3_b1_b2.30m<-shannon(tablo_g1_g2_g3_b1_b2_30m)
 selection.shannon.hab.g1_g2_g3_b1_b2.60m<-shannon(tablo_g1_g2_g3_b1_b2_60m)
 selection.shannon.hab.g1_g2_g3_b1_b2.100m<-shannon(tablo_g1_g2_g3_b1_b2_100m)
 selection.shannon.hab.g1_g2_g3_b1_b2.250m<-shannon(tablo_g1_g2_g3_b1_b2_250m)
 selection.shannon.hab.g1_g2_g3_b1_b2.500m<-shannon(tablo_g1_g2_g3_b1_b2_500m)
 selection.shannon.hab.g1_g2_g3_b1_b2.750m<-shannon(tablo_g1_g2_g3_b1_b2_750m)
 selection.shannon.hab.g1_g2_g3_b1_b2.1000m<-shannon(tablo_g1_g2_g3_b1_b2_1000m)
 selection.shannon.hab.g1_g2_g3_b1_b2.3000m<-shannon(tablo_g1_g2_g3_b1_b2_3000m)
 selection.shannon.hab.g1_g2_g3_b1_b2.5000m<-shannon(tablo_g1_g2_g3_b1_b2_5000m)

  #DIVERSITE SELECTIONNE PAR LE SCRIPT SHANNON
  source('script.diversite.shannon.selectionnee.fonction.r')
diversite.selection.shannon.hab.g1_g2_g3_b1_b2.30m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2.30m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1_b2.60m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2.60m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1_b2.100m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2.100m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1_b2.250m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2.250m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1_b2.500m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2.500m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1_b2.750m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2.750m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1_b2.1000m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2.1000m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1_b2.3000m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2.3000m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1_b2.5000m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2.5000m,requette.abondance.poissons)

 #CALCUL DES SAI#
source('script.sai.r')
SAI_shannon.g1_g2_g3_b1_b2.30m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2.30m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1_b2.60m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2.60m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1_b2.100m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2.100m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1_b2.250m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2.250m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1_b2.500m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2.500m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1_b2.750m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2.750m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1_b2.1000m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2.1000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1_b2.3000m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2.3000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1_b2.5000m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2.5000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)


  #############Habitat g1_g2_g3_b1_b2_b3###################



   ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances

setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
relation.habitats.sites.g1_g2_g3_b1_b2_b3<-read.csv2("fc3_intersect_g1_g2_g3_b1_b2_b3.csv")

sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')

relation.habitats.sites.g1_g2_g3_b1_b2_b3<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$Code_Site%in%sites.etudies,]

#création d'une nouvelle table site-habitat par distance
relation.habitats.sites.g1_g2_g3_b1_b2_b3.30m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==30,]
relation.habitats.sites.g1_g2_g3_b1_b2_b3.60m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==60,]
relation.habitats.sites.g1_g2_g3_b1_b2_b3.100m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==100,]
relation.habitats.sites.g1_g2_g3_b1_b2_b3.250m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==250,]
relation.habitats.sites.g1_g2_g3_b1_b2_b3.500m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==500,]
relation.habitats.sites.g1_g2_g3_b1_b2_b3.750m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==750,]
relation.habitats.sites.g1_g2_g3_b1_b2_b3.1000m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==1000,]
relation.habitats.sites.g1_g2_g3_b1_b2_b3.3000m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==3000,]
relation.habitats.sites.g1_g2_g3_b1_b2_b3.5000m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==5000,]

#cration d'un tableau adapté à l'entrée du script shannon
relation.habitats.sites.g1_g2_g3_b1_b2_b3.30m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2_b3.30m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2_b3.30m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2_b3.30m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2_b3.30m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_b3_30m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2_b3.30m.unique$surface~relation.habitats.sites.g1_g2_g3_b1_b2_b3.30m.unique$Code_Site+relation.habitats.sites.g1_g2_g3_b1_b2_b3.30m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2_b3.30m.unique)

relation.habitats.sites.g1_g2_g3_b1_b2_b3.60m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2_b3.60m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2_b3.60m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2_b3.60m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2_b3.60m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_b3_60m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2_b3.60m.unique$surface~relation.habitats.sites.g1_g2_g3_b1_b2_b3.60m.unique$Code_Site+relation.habitats.sites.g1_g2_g3_b1_b2_b3.60m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2_b3.60m.unique)

relation.habitats.sites.g1_g2_g3_b1_b2_b3.100m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2_b3.100m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2_b3.100m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2_b3.100m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2_b3.100m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_b3_100m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2_b3.100m.unique$surface~relation.habitats.sites.g1_g2_g3_b1_b2_b3.100m.unique$Code_Site+relation.habitats.sites.g1_g2_g3_b1_b2_b3.100m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2_b3.100m.unique)

relation.habitats.sites.g1_g2_g3_b1_b2_b3.250m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2_b3.250m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2_b3.250m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2_b3.250m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2_b3.250m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_b3_250m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2_b3.250m.unique$surface~relation.habitats.sites.g1_g2_g3_b1_b2_b3.250m.unique$Code_Site+relation.habitats.sites.g1_g2_g3_b1_b2_b3.250m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2_b3.250m.unique)

relation.habitats.sites.g1_g2_g3_b1_b2_b3.500m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2_b3.500m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2_b3.500m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2_b3.500m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2_b3.500m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_b3_500m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2_b3.500m.unique$surface~relation.habitats.sites.g1_g2_g3_b1_b2_b3.500m.unique$Code_Site+relation.habitats.sites.g1_g2_g3_b1_b2_b3.500m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2_b3.500m.unique)

relation.habitats.sites.g1_g2_g3_b1_b2_b3.750m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2_b3.750m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2_b3.750m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2_b3.750m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2_b3.750m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_b3_750m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2_b3.750m.unique$surface~relation.habitats.sites.g1_g2_g3_b1_b2_b3.750m.unique$Code_Site+relation.habitats.sites.g1_g2_g3_b1_b2_b3.750m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2_b3.750m.unique)

relation.habitats.sites.g1_g2_g3_b1_b2_b3.1000m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2_b3.1000m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2_b3.1000m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2_b3.1000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2_b3.1000m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_b3_1000m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2_b3.1000m.unique$surface~relation.habitats.sites.g1_g2_g3_b1_b2_b3.1000m.unique$Code_Site+relation.habitats.sites.g1_g2_g3_b1_b2_b3.1000m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2_b3.1000m.unique)

relation.habitats.sites.g1_g2_g3_b1_b2_b3.3000m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2_b3.3000m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2_b3.3000m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2_b3.3000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2_b3.3000m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_b3_3000m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2_b3.3000m.unique$surface~relation.habitats.sites.g1_g2_g3_b1_b2_b3.3000m.unique$Code_Site+relation.habitats.sites.g1_g2_g3_b1_b2_b3.3000m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2_b3.3000m.unique)

relation.habitats.sites.g1_g2_g3_b1_b2_b3.5000m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2_b3.5000m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2_b3.5000m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2_b3.5000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2_b3.5000m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_b3_5000m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2_b3.5000m.unique$surface~relation.habitats.sites.g1_g2_g3_b1_b2_b3.5000m.unique$Code_Site+relation.habitats.sites.g1_g2_g3_b1_b2_b3.5000m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2_b3.5000m.unique)

tablo_g1_g2_g3_b1_b2_b3_30m<-tablo_g1_g2_g3_b1_b2_b3_30m[25:51,]
tablo_g1_g2_g3_b1_b2_b3_60m<-tablo_g1_g2_g3_b1_b2_b3_60m[25:51,]
tablo_g1_g2_g3_b1_b2_b3_100m<-tablo_g1_g2_g3_b1_b2_b3_100m[25:51,]
tablo_g1_g2_g3_b1_b2_b3_250m<-tablo_g1_g2_g3_b1_b2_b3_250m[25:51,]
tablo_g1_g2_g3_b1_b2_b3_500m<-tablo_g1_g2_g3_b1_b2_b3_500m[25:51,]
tablo_g1_g2_g3_b1_b2_b3_750m<-tablo_g1_g2_g3_b1_b2_b3_750m[25:51,]
tablo_g1_g2_g3_b1_b2_b3_1000m<-tablo_g1_g2_g3_b1_b2_b3_1000m[25:51,]
tablo_g1_g2_g3_b1_b2_b3_3000m<-tablo_g1_g2_g3_b1_b2_b3_3000m[25:51,]
tablo_g1_g2_g3_b1_b2_b3_5000m<-tablo_g1_g2_g3_b1_b2_b3_5000m[25:51,]

#Load de l'algorithme shannon
source("script.shannon.r")

#tirage des sites poissons en optimisant la diversité shannon d'habitat
 selection.shannon.hab.g1_g2_g3_b1_b2_b3.30m<-shannon(tablo_g1_g2_g3_b1_b2_b3_30m)
 selection.shannon.hab.g1_g2_g3_b1_b2_b3.60m<-shannon(tablo_g1_g2_g3_b1_b2_b3_60m)
 selection.shannon.hab.g1_g2_g3_b1_b2_b3.100m<-shannon(tablo_g1_g2_g3_b1_b2_b3_100m)
 selection.shannon.hab.g1_g2_g3_b1_b2_b3.250m<-shannon(tablo_g1_g2_g3_b1_b2_b3_250m)
 selection.shannon.hab.g1_g2_g3_b1_b2_b3.500m<-shannon(tablo_g1_g2_g3_b1_b2_b3_500m)
 selection.shannon.hab.g1_g2_g3_b1_b2_b3.750m<-shannon(tablo_g1_g2_g3_b1_b2_b3_750m)
 selection.shannon.hab.g1_g2_g3_b1_b2_b3.1000m<-shannon(tablo_g1_g2_g3_b1_b2_b3_1000m)
 selection.shannon.hab.g1_g2_g3_b1_b2_b3.3000m<-shannon(tablo_g1_g2_g3_b1_b2_b3_3000m)
 selection.shannon.hab.g1_g2_g3_b1_b2_b3.5000m<-shannon(tablo_g1_g2_g3_b1_b2_b3_5000m)

  #DIVERSITE SELECTIONNE PAR LE SCRIPT SHANNON
  source('script.diversite.shannon.selectionnee.fonction.r')
diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.30m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2_b3.30m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.60m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2_b3.60m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.100m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2_b3.100m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.250m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2_b3.250m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.500m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2_b3.500m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.750m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2_b3.750m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.1000m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2_b3.1000m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.3000m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2_b3.3000m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.5000m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2_b3.5000m,requette.abondance.poissons)

 #CALCUL DES SAI#
source('script.sai.r')
SAI_shannon.g1_g2_g3_b1_b2_b3.30m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.30m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1_b2_b3.60m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.60m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1_b2_b3.100m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.100m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1_b2_b3.250m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.250m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1_b2_b3.500m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.500m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1_b2_b3.750m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.750m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1_b2_b3.1000m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.1000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1_b2_b3.3000m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.3000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1_b2_b3.5000m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.5000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)


  #STOCKAGE DES RESULTATS

targeted.station.results[,1]=c(SAI_shannon.geo1.30m,SAI_shannon.geo1.60m,SAI_shannon.geo1.100m,SAI_shannon.geo1.250m,
SAI_shannon.geo1.500m,SAI_shannon.geo1.750m,SAI_shannon.geo1.1000m,SAI_shannon.geo1.3000m,SAI_shannon.geo1.5000m,
SAI_shannon.geo2.30m,SAI_shannon.geo2.60m,SAI_shannon.geo2.100m,SAI_shannon.geo2.250m,
SAI_shannon.geo2.500m,SAI_shannon.geo2.750m,SAI_shannon.geo2.1000m,SAI_shannon.geo2.3000m,SAI_shannon.geo2.5000m,
SAI_shannon.geo3.100m,SAI_shannon.geo3.250m,
SAI_shannon.geo3.500m,SAI_shannon.geo3.750m,SAI_shannon.geo3.1000m,SAI_shannon.geo3.3000m,SAI_shannon.geo3.5000m,
SAI_shannon.bent1.30m,SAI_shannon.bent1.60m,SAI_shannon.bent1.100m,SAI_shannon.bent1.250m,
SAI_shannon.bent1.500m,SAI_shannon.bent1.750m,SAI_shannon.bent1.1000m,SAI_shannon.bent1.3000m,SAI_shannon.bent1.5000m,
SAI_shannon.bent2.30m,SAI_shannon.bent2.60m,SAI_shannon.bent2.100m,SAI_shannon.bent2.250m,
SAI_shannon.bent2.500m,SAI_shannon.bent2.750m,SAI_shannon.bent2.1000m,SAI_shannon.bent2.3000m,SAI_shannon.bent2.5000m,
SAI_shannon.bent3.30m,SAI_shannon.bent3.60m,SAI_shannon.bent3.100m,SAI_shannon.bent3.250m,
SAI_shannon.bent3.500m,SAI_shannon.bent3.750m,SAI_shannon.bent3.1000m,SAI_shannon.bent3.3000m,SAI_shannon.bent3.5000m,
SAI_shannon.geo1_geo2.30m,SAI_shannon.geo1_geo2.60m,SAI_shannon.geo1_geo2.100m,SAI_shannon.geo1_geo2.250m,
SAI_shannon.geo1_geo2.500m,SAI_shannon.geo1_geo2.750m,SAI_shannon.geo1_geo2.1000m,SAI_shannon.geo1_geo2.3000m,SAI_shannon.geo1_geo2.5000m,
SAI_shannon.geo1_geo2_geo3.30m,SAI_shannon.geo1_geo2_geo3.60m,SAI_shannon.geo1_geo2_geo3.100m,SAI_shannon.geo1_geo2_geo3.250m,
SAI_shannon.geo1_geo2_geo3.500m,SAI_shannon.geo1_geo2_geo3.750m,SAI_shannon.geo1_geo2_geo3.1000m,SAI_shannon.geo1_geo2_geo3.3000m,SAI_shannon.geo1_geo2_geo3.5000m,
SAI_shannon.g1_g2_g3_b1.30m,SAI_shannon.g1_g2_g3_b1.60m,SAI_shannon.g1_g2_g3_b1.100m,SAI_shannon.g1_g2_g3_b1.250m,
SAI_shannon.g1_g2_g3_b1.500m,SAI_shannon.g1_g2_g3_b1.750m,SAI_shannon.g1_g2_g3_b1.1000m,SAI_shannon.g1_g2_g3_b1.3000m,SAI_shannon.g1_g2_g3_b1.5000m,
SAI_shannon.g1_g2_g3_b1_b2.30m,SAI_shannon.g1_g2_g3_b1_b2.60m,SAI_shannon.g1_g2_g3_b1_b2.100m,SAI_shannon.g1_g2_g3_b1_b2.250m,
SAI_shannon.g1_g2_g3_b1_b2.500m,SAI_shannon.g1_g2_g3_b1_b2.750m,SAI_shannon.g1_g2_g3_b1_b2.1000m,SAI_shannon.g1_g2_g3_b1_b2.3000m,SAI_shannon.g1_g2_g3_b1_b2.5000m,
SAI_shannon.g1_g2_g3_b1_b2_b3.30m,SAI_shannon.g1_g2_g3_b1_b2_b3.60m,SAI_shannon.g1_g2_g3_b1_b2_b3.100m,SAI_shannon.g1_g2_g3_b1_b2_b3.250m,
SAI_shannon.g1_g2_g3_b1_b2_b3.500m,SAI_shannon.g1_g2_g3_b1_b2_b3.750m,SAI_shannon.g1_g2_g3_b1_b2_b3.1000m,SAI_shannon.g1_g2_g3_b1_b2_b3.3000m,SAI_shannon.g1_g2_g3_b1_b2_b3.5000m)

 
 
 
 
 
 
 
 
 
 #-1 STATION TESTS



targeted.stations<-c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')

for (p in 1:length(targeted.stations)){
  
suppressed.station<-targeted.stations[p]

#Load du tableau : "site" "poissons"


setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
recensement=read.csv2("REQ_site_fonction.IRD.csv")

#j'enlève la station supprimée
recensement<-recensement[!recensement$Code_Site%in%suppressed.station,]

#J'enlève les doublons
recensement_unique=unique(recensement)

#création d'un tableau de présence absence
tablo_pres_abs<-table(as.character(recensement_unique$Code_Site),recensement_unique$Code_CT_SC_M)

#calcul du nombre d'espèce par station
vec_poisson<-apply(tablo_pres_abs,1,sum)

#selection aléatoire d'un pool de site
setwd(repertoire)
setwd(input)
source('script.aleatoire.r')
selection.aleatoire<-matrix(0,dim(tablo_pres_abs)[1],999) #999
for (i in 1:999){ #999
selection.aleatoire[,i]<-run.alea(tablo_pres_abs)}

#Selection des sites selon la diversité de poissons

setwd(repertoire)
 setwd(input)
 requette.abondance.poissons<-read.csv2('abondance.midif_sim_output.csv')
  requette.abondance.poissons<-requette.abondance.poissons[!requette.abondance.poissons$Code_Site%in%suppressed.station,]
  requette.abondance.poissons$Code_Site<-factor(requette.abondance.poissons$Code_Site)
input.shannon.abondance<-xtabs(formula=requette.abondance.poissons$Abondance~requette.abondance.poissons$Code_Site+requette.abondance.poissons$Code_CT_SC_M,requette.abondance.poissons)

source('script.shannon.poissons.r')
selection.shannon.poisson<-shannon.poisson(input.shannon.abondance)

#Diversité de poisson incluse dans les sites selectionnés (aléatoirement et selon diversité de poisson)
  source('script.diversite.shannon.selectionnee.fonction.r')
 diversite.shannon.selection.aleatoire<-matrix(0,dim(tablo_pres_abs)[1],999)   #999
 diversite.shannon.selection.aleatoire<-apply(as.matrix(selection.aleatoire),2,function(x) species.shannon.evol(x,requette.abondance.poissons))
 diversite.shannon.moyenne.selection.aleatoire<-apply(diversite.shannon.selection.aleatoire,1,mean)

 diversite.selection.shannon.poisson<-species.shannon.evol(selection.shannon.poisson,requette.abondance.poissons)

 #Calcul des courbes aléatore moyenne, uper, et low
 diversite.shannon.moyenne.selection.aleatoire
diversite.shannon.borne.sup.selection.aleatoire=apply(diversite.shannon.selection.aleatoire,1,function(x) sort(x)[ceiling(0.975*length(x))])
diversite.shannon.borne.inf.selection.aleatoire=apply(diversite.shannon.selection.aleatoire,1,function(x) sort(x)[ceiling(0.025*length(x))])

#nombre total d'espèces de poissons observés
nombre.poissons.total <- nlevels(factor(recensement_unique$Code_CT_SC_M))



  #############Habitat Geo1###################



   ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances

setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
relation.habitats.sites.geo1<-read.csv2("fc3_intersect_geo1.csv")

sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')
sites.etudies<-sites.etudies[!sites.etudies%in%suppressed.station]

relation.habitats.sites.geo1<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$Code_Site%in%sites.etudies,]

#création d'une nouvelle table site-habitat par distance
relation.habitats.sites.geo1.30m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==30,]
relation.habitats.sites.geo1.60m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==60,]
relation.habitats.sites.geo1.100m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==100,]
relation.habitats.sites.geo1.250m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==250,]
relation.habitats.sites.geo1.500m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==500,]
relation.habitats.sites.geo1.750m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==750,]
relation.habitats.sites.geo1.1000m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==1000,]
relation.habitats.sites.geo1.3000m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==3000,]
relation.habitats.sites.geo1.5000m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==5000,]

#cration d'un tableau adapté à l'entrée du script shannon
relation.habitats.sites.geo1.30m.unique<-aggregate(relation.habitats.sites.geo1.30m$surface,list(Code_Site=relation.habitats.sites.geo1.30m$Code_Site,habitat=relation.habitats.sites.geo1.30m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1.30m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_30m<-xtabs(formula=relation.habitats.sites.geo1.30m.unique$surface~as.character(relation.habitats.sites.geo1.30m.unique$Code_Site)+relation.habitats.sites.geo1.30m.unique$habitat,relation.habitats.sites.geo1.30m.unique)

relation.habitats.sites.geo1.60m.unique<-aggregate(relation.habitats.sites.geo1.60m$surface,list(Code_Site=relation.habitats.sites.geo1.60m$Code_Site,habitat=relation.habitats.sites.geo1.60m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1.60m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_60m<-xtabs(formula=relation.habitats.sites.geo1.60m.unique$surface~as.character(relation.habitats.sites.geo1.60m.unique$Code_Site)+relation.habitats.sites.geo1.60m.unique$habitat,relation.habitats.sites.geo1.60m.unique)

relation.habitats.sites.geo1.100m.unique<-aggregate(relation.habitats.sites.geo1.100m$surface,list(Code_Site=relation.habitats.sites.geo1.100m$Code_Site,habitat=relation.habitats.sites.geo1.100m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1.100m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_100m<-xtabs(formula=relation.habitats.sites.geo1.100m.unique$surface~as.character(relation.habitats.sites.geo1.100m.unique$Code_Site)+relation.habitats.sites.geo1.100m.unique$habitat,relation.habitats.sites.geo1.100m.unique)

relation.habitats.sites.geo1.250m.unique<-aggregate(relation.habitats.sites.geo1.250m$surface,list(Code_Site=relation.habitats.sites.geo1.250m$Code_Site,habitat=relation.habitats.sites.geo1.250m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1.250m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_250m<-xtabs(formula=relation.habitats.sites.geo1.250m.unique$surface~as.character(relation.habitats.sites.geo1.250m.unique$Code_Site)+relation.habitats.sites.geo1.250m.unique$habitat,relation.habitats.sites.geo1.250m.unique)

relation.habitats.sites.geo1.500m.unique<-aggregate(relation.habitats.sites.geo1.500m$surface,list(Code_Site=relation.habitats.sites.geo1.500m$Code_Site,habitat=relation.habitats.sites.geo1.500m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1.500m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_500m<-xtabs(formula=relation.habitats.sites.geo1.500m.unique$surface~as.character(relation.habitats.sites.geo1.500m.unique$Code_Site)+relation.habitats.sites.geo1.500m.unique$habitat,relation.habitats.sites.geo1.500m.unique)

relation.habitats.sites.geo1.750m.unique<-aggregate(relation.habitats.sites.geo1.750m$surface,list(Code_Site=relation.habitats.sites.geo1.750m$Code_Site,habitat=relation.habitats.sites.geo1.750m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1.750m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_750m<-xtabs(formula=relation.habitats.sites.geo1.750m.unique$surface~as.character(relation.habitats.sites.geo1.750m.unique$Code_Site)+relation.habitats.sites.geo1.750m.unique$habitat,relation.habitats.sites.geo1.750m.unique)

relation.habitats.sites.geo1.1000m.unique<-aggregate(relation.habitats.sites.geo1.1000m$surface,list(Code_Site=relation.habitats.sites.geo1.1000m$Code_Site,habitat=relation.habitats.sites.geo1.1000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1.1000m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_1000m<-xtabs(formula=relation.habitats.sites.geo1.1000m.unique$surface~as.character(relation.habitats.sites.geo1.1000m.unique$Code_Site)+relation.habitats.sites.geo1.1000m.unique$habitat,relation.habitats.sites.geo1.1000m.unique)

relation.habitats.sites.geo1.3000m.unique<-aggregate(relation.habitats.sites.geo1.3000m$surface,list(Code_Site=relation.habitats.sites.geo1.3000m$Code_Site,habitat=relation.habitats.sites.geo1.3000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1.3000m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_3000m<-xtabs(formula=relation.habitats.sites.geo1.3000m.unique$surface~as.character(relation.habitats.sites.geo1.3000m.unique$Code_Site)+relation.habitats.sites.geo1.3000m.unique$habitat,relation.habitats.sites.geo1.3000m.unique)

relation.habitats.sites.geo1.5000m.unique<-aggregate(relation.habitats.sites.geo1.5000m$surface,list(Code_Site=relation.habitats.sites.geo1.5000m$Code_Site,habitat=relation.habitats.sites.geo1.5000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1.5000m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_5000m<-xtabs(formula=relation.habitats.sites.geo1.5000m.unique$surface~as.character(relation.habitats.sites.geo1.5000m.unique$Code_Site)+relation.habitats.sites.geo1.5000m.unique$habitat,relation.habitats.sites.geo1.5000m.unique)

#Load de l'algorithme shannon
source("script.shannon.r")

#tirage des sites poissons en optimisant la diversité shannon d'habitat
 selection.shannon.hab.geo1.30m<-shannon(tablo_geo1_30m)
 selection.shannon.hab.geo1.60m<-shannon(tablo_geo1_60m)
 selection.shannon.hab.geo1.100m<-shannon(tablo_geo1_100m)
 selection.shannon.hab.geo1.250m<-shannon(tablo_geo1_250m)
 selection.shannon.hab.geo1.500m<-shannon(tablo_geo1_500m)
 selection.shannon.hab.geo1.750m<-shannon(tablo_geo1_750m)
 selection.shannon.hab.geo1.1000m<-shannon(tablo_geo1_1000m)
 selection.shannon.hab.geo1.3000m<-shannon(tablo_geo1_3000m)
 selection.shannon.hab.geo1.5000m<-shannon(tablo_geo1_5000m)

  #DIVERSITE SELECTIONNE PAR LE SCRIPT SHANNON
  source('script.diversite.shannon.selectionnee.fonction.r')
diversite.selection.shannon.hab.geo1.30m<-species.shannon.evol(selection.shannon.hab.geo1.30m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1.60m<-species.shannon.evol(selection.shannon.hab.geo1.60m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1.100m<-species.shannon.evol(selection.shannon.hab.geo1.100m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1.250m<-species.shannon.evol(selection.shannon.hab.geo1.250m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1.500m<-species.shannon.evol(selection.shannon.hab.geo1.500m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1.750m<-species.shannon.evol(selection.shannon.hab.geo1.750m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1.1000m<-species.shannon.evol(selection.shannon.hab.geo1.1000m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1.3000m<-species.shannon.evol(selection.shannon.hab.geo1.3000m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1.5000m<-species.shannon.evol(selection.shannon.hab.geo1.5000m,requette.abondance.poissons)

 #CALCUL DES SAI#
source('script.sai.r')
SAI_shannon.geo1.30m<-SAI(diversite.selection.shannon.hab.geo1.30m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1.60m<-SAI(diversite.selection.shannon.hab.geo1.60m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1.100m<-SAI(diversite.selection.shannon.hab.geo1.100m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1.250m<-SAI(diversite.selection.shannon.hab.geo1.250m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1.500m<-SAI(diversite.selection.shannon.hab.geo1.500m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1.750m<-SAI(diversite.selection.shannon.hab.geo1.750m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1.1000m<-SAI(diversite.selection.shannon.hab.geo1.1000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1.3000m<-SAI(diversite.selection.shannon.hab.geo1.3000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1.5000m<-SAI(diversite.selection.shannon.hab.geo1.5000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)

  
         

  #############Habitat geo2###################



   ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances

setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
relation.habitats.sites.geo2<-read.csv2("fc3_intersect_geo2.csv")

sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')
sites.etudies<-sites.etudies[!sites.etudies%in%suppressed.station]

relation.habitats.sites.geo2<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$Code_Site%in%sites.etudies,]

#création d'une nouvelle table site-habitat par distance
relation.habitats.sites.geo2.30m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==30,]
relation.habitats.sites.geo2.60m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==60,]
relation.habitats.sites.geo2.100m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==100,]
relation.habitats.sites.geo2.250m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==250,]
relation.habitats.sites.geo2.500m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==500,]
relation.habitats.sites.geo2.750m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==750,]
relation.habitats.sites.geo2.1000m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==1000,]
relation.habitats.sites.geo2.3000m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==3000,]
relation.habitats.sites.geo2.5000m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==5000,]

#cration d'un tableau adapté à l'entrée du script shannon
relation.habitats.sites.geo2.30m.unique<-aggregate(relation.habitats.sites.geo2.30m$surface,list(Code_Site=relation.habitats.sites.geo2.30m$Code_Site,habitat=relation.habitats.sites.geo2.30m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo2.30m.unique)=c("Code_Site","habitat","surface")
tablo_geo2_30m<-xtabs(formula=relation.habitats.sites.geo2.30m.unique$surface~as.character(relation.habitats.sites.geo2.30m.unique$Code_Site)+relation.habitats.sites.geo2.30m.unique$habitat,relation.habitats.sites.geo2.30m.unique)

relation.habitats.sites.geo2.60m.unique<-aggregate(relation.habitats.sites.geo2.60m$surface,list(Code_Site=relation.habitats.sites.geo2.60m$Code_Site,habitat=relation.habitats.sites.geo2.60m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo2.60m.unique)=c("Code_Site","habitat","surface")
tablo_geo2_60m<-xtabs(formula=relation.habitats.sites.geo2.60m.unique$surface~as.character(relation.habitats.sites.geo2.60m.unique$Code_Site)+relation.habitats.sites.geo2.60m.unique$habitat,relation.habitats.sites.geo2.60m.unique)

relation.habitats.sites.geo2.100m.unique<-aggregate(relation.habitats.sites.geo2.100m$surface,list(Code_Site=relation.habitats.sites.geo2.100m$Code_Site,habitat=relation.habitats.sites.geo2.100m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo2.100m.unique)=c("Code_Site","habitat","surface")
tablo_geo2_100m<-xtabs(formula=relation.habitats.sites.geo2.100m.unique$surface~as.character(relation.habitats.sites.geo2.100m.unique$Code_Site)+relation.habitats.sites.geo2.100m.unique$habitat,relation.habitats.sites.geo2.100m.unique)

relation.habitats.sites.geo2.250m.unique<-aggregate(relation.habitats.sites.geo2.250m$surface,list(Code_Site=relation.habitats.sites.geo2.250m$Code_Site,habitat=relation.habitats.sites.geo2.250m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo2.250m.unique)=c("Code_Site","habitat","surface")
tablo_geo2_250m<-xtabs(formula=relation.habitats.sites.geo2.250m.unique$surface~as.character(relation.habitats.sites.geo2.250m.unique$Code_Site)+relation.habitats.sites.geo2.250m.unique$habitat,relation.habitats.sites.geo2.250m.unique)

relation.habitats.sites.geo2.500m.unique<-aggregate(relation.habitats.sites.geo2.500m$surface,list(Code_Site=relation.habitats.sites.geo2.500m$Code_Site,habitat=relation.habitats.sites.geo2.500m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo2.500m.unique)=c("Code_Site","habitat","surface")
tablo_geo2_500m<-xtabs(formula=relation.habitats.sites.geo2.500m.unique$surface~as.character(relation.habitats.sites.geo2.500m.unique$Code_Site)+relation.habitats.sites.geo2.500m.unique$habitat,relation.habitats.sites.geo2.500m.unique)

relation.habitats.sites.geo2.750m.unique<-aggregate(relation.habitats.sites.geo2.750m$surface,list(Code_Site=relation.habitats.sites.geo2.750m$Code_Site,habitat=relation.habitats.sites.geo2.750m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo2.750m.unique)=c("Code_Site","habitat","surface")
tablo_geo2_750m<-xtabs(formula=relation.habitats.sites.geo2.750m.unique$surface~as.character(relation.habitats.sites.geo2.750m.unique$Code_Site)+relation.habitats.sites.geo2.750m.unique$habitat,relation.habitats.sites.geo2.750m.unique)

relation.habitats.sites.geo2.1000m.unique<-aggregate(relation.habitats.sites.geo2.1000m$surface,list(Code_Site=relation.habitats.sites.geo2.1000m$Code_Site,habitat=relation.habitats.sites.geo2.1000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo2.1000m.unique)=c("Code_Site","habitat","surface")
tablo_geo2_1000m<-xtabs(formula=relation.habitats.sites.geo2.1000m.unique$surface~as.character(relation.habitats.sites.geo2.1000m.unique$Code_Site)+relation.habitats.sites.geo2.1000m.unique$habitat,relation.habitats.sites.geo2.1000m.unique)

relation.habitats.sites.geo2.3000m.unique<-aggregate(relation.habitats.sites.geo2.3000m$surface,list(Code_Site=relation.habitats.sites.geo2.3000m$Code_Site,habitat=relation.habitats.sites.geo2.3000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo2.3000m.unique)=c("Code_Site","habitat","surface")
tablo_geo2_3000m<-xtabs(formula=relation.habitats.sites.geo2.3000m.unique$surface~as.character(relation.habitats.sites.geo2.3000m.unique$Code_Site)+relation.habitats.sites.geo2.3000m.unique$habitat,relation.habitats.sites.geo2.3000m.unique)

relation.habitats.sites.geo2.5000m.unique<-aggregate(relation.habitats.sites.geo2.5000m$surface,list(Code_Site=relation.habitats.sites.geo2.5000m$Code_Site,habitat=relation.habitats.sites.geo2.5000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo2.5000m.unique)=c("Code_Site","habitat","surface")
tablo_geo2_5000m<-xtabs(formula=relation.habitats.sites.geo2.5000m.unique$surface~as.character(relation.habitats.sites.geo2.5000m.unique$Code_Site)+relation.habitats.sites.geo2.5000m.unique$habitat,relation.habitats.sites.geo2.5000m.unique)

#Load de l'algorithme shannon
source("script.shannon.r")

#tirage des sites poissons en optimisant la diversité shannon d'habitat
 selection.shannon.hab.geo2.30m<-shannon(tablo_geo2_30m)
 selection.shannon.hab.geo2.60m<-shannon(tablo_geo2_60m)
 selection.shannon.hab.geo2.100m<-shannon(tablo_geo2_100m)
 selection.shannon.hab.geo2.250m<-shannon(tablo_geo2_250m)
 selection.shannon.hab.geo2.500m<-shannon(tablo_geo2_500m)
 selection.shannon.hab.geo2.750m<-shannon(tablo_geo2_750m)
 selection.shannon.hab.geo2.1000m<-shannon(tablo_geo2_1000m)
 selection.shannon.hab.geo2.3000m<-shannon(tablo_geo2_3000m)
 selection.shannon.hab.geo2.5000m<-shannon(tablo_geo2_5000m)

  #DIVERSITE SELECTIONNE PAR LE SCRIPT SHANNON
  source('script.diversite.shannon.selectionnee.fonction.r')
diversite.selection.shannon.hab.geo2.30m<-species.shannon.evol(selection.shannon.hab.geo2.30m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo2.60m<-species.shannon.evol(selection.shannon.hab.geo2.60m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo2.100m<-species.shannon.evol(selection.shannon.hab.geo2.100m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo2.250m<-species.shannon.evol(selection.shannon.hab.geo2.250m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo2.500m<-species.shannon.evol(selection.shannon.hab.geo2.500m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo2.750m<-species.shannon.evol(selection.shannon.hab.geo2.750m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo2.1000m<-species.shannon.evol(selection.shannon.hab.geo2.1000m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo2.3000m<-species.shannon.evol(selection.shannon.hab.geo2.3000m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo2.5000m<-species.shannon.evol(selection.shannon.hab.geo2.5000m,requette.abondance.poissons)

 #CALCUL DES SAI#
source('script.sai.r')
SAI_shannon.geo2.30m<-SAI(diversite.selection.shannon.hab.geo2.30m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo2.60m<-SAI(diversite.selection.shannon.hab.geo2.60m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo2.100m<-SAI(diversite.selection.shannon.hab.geo2.100m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo2.250m<-SAI(diversite.selection.shannon.hab.geo2.250m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo2.500m<-SAI(diversite.selection.shannon.hab.geo2.500m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo2.750m<-SAI(diversite.selection.shannon.hab.geo2.750m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo2.1000m<-SAI(diversite.selection.shannon.hab.geo2.1000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo2.3000m<-SAI(diversite.selection.shannon.hab.geo2.3000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo2.5000m<-SAI(diversite.selection.shannon.hab.geo2.5000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)



  #############Habitat geo3###################



   ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances

setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
relation.habitats.sites.geo3<-read.csv2("fc3_intersect_geo3.csv")

sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')
sites.etudies<-sites.etudies[!sites.etudies%in%suppressed.station]

relation.habitats.sites.geo3<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$Code_Site%in%sites.etudies,]

#création d'une nouvelle table site-habitat par distance
relation.habitats.sites.geo3.100m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==100,]
relation.habitats.sites.geo3.250m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==250,]
relation.habitats.sites.geo3.500m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==500,]
relation.habitats.sites.geo3.750m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==750,]
relation.habitats.sites.geo3.1000m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==1000,]
relation.habitats.sites.geo3.3000m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==3000,]
relation.habitats.sites.geo3.5000m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==5000,]

#cration d'un tableau adapté à l'entrée du script shannon

relation.habitats.sites.geo3.100m.unique<-aggregate(relation.habitats.sites.geo3.100m$surface,list(Code_Site=relation.habitats.sites.geo3.100m$Code_Site,habitat=relation.habitats.sites.geo3.100m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo3.100m.unique)=c("Code_Site","habitat","surface")
tablo_geo3_100m<-xtabs(formula=relation.habitats.sites.geo3.100m.unique$surface~as.character(relation.habitats.sites.geo3.100m.unique$Code_Site)+relation.habitats.sites.geo3.100m.unique$habitat,relation.habitats.sites.geo3.100m.unique)

relation.habitats.sites.geo3.250m.unique<-aggregate(relation.habitats.sites.geo3.250m$surface,list(Code_Site=relation.habitats.sites.geo3.250m$Code_Site,habitat=relation.habitats.sites.geo3.250m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo3.250m.unique)=c("Code_Site","habitat","surface")
tablo_geo3_250m<-xtabs(formula=relation.habitats.sites.geo3.250m.unique$surface~as.character(relation.habitats.sites.geo3.250m.unique$Code_Site)+relation.habitats.sites.geo3.250m.unique$habitat,relation.habitats.sites.geo3.250m.unique)

relation.habitats.sites.geo3.500m.unique<-aggregate(relation.habitats.sites.geo3.500m$surface,list(Code_Site=relation.habitats.sites.geo3.500m$Code_Site,habitat=relation.habitats.sites.geo3.500m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo3.500m.unique)=c("Code_Site","habitat","surface")
tablo_geo3_500m<-xtabs(formula=relation.habitats.sites.geo3.500m.unique$surface~as.character(relation.habitats.sites.geo3.500m.unique$Code_Site)+relation.habitats.sites.geo3.500m.unique$habitat,relation.habitats.sites.geo3.500m.unique)

relation.habitats.sites.geo3.750m.unique<-aggregate(relation.habitats.sites.geo3.750m$surface,list(Code_Site=relation.habitats.sites.geo3.750m$Code_Site,habitat=relation.habitats.sites.geo3.750m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo3.750m.unique)=c("Code_Site","habitat","surface")
tablo_geo3_750m<-xtabs(formula=relation.habitats.sites.geo3.750m.unique$surface~as.character(relation.habitats.sites.geo3.750m.unique$Code_Site)+relation.habitats.sites.geo3.750m.unique$habitat,relation.habitats.sites.geo3.750m.unique)

relation.habitats.sites.geo3.1000m.unique<-aggregate(relation.habitats.sites.geo3.1000m$surface,list(Code_Site=relation.habitats.sites.geo3.1000m$Code_Site,habitat=relation.habitats.sites.geo3.1000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo3.1000m.unique)=c("Code_Site","habitat","surface")
tablo_geo3_1000m<-xtabs(formula=relation.habitats.sites.geo3.1000m.unique$surface~as.character(relation.habitats.sites.geo3.1000m.unique$Code_Site)+relation.habitats.sites.geo3.1000m.unique$habitat,relation.habitats.sites.geo3.1000m.unique)

relation.habitats.sites.geo3.3000m.unique<-aggregate(relation.habitats.sites.geo3.3000m$surface,list(Code_Site=relation.habitats.sites.geo3.3000m$Code_Site,habitat=relation.habitats.sites.geo3.3000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo3.3000m.unique)=c("Code_Site","habitat","surface")
tablo_geo3_3000m<-xtabs(formula=relation.habitats.sites.geo3.3000m.unique$surface~as.character(relation.habitats.sites.geo3.3000m.unique$Code_Site)+relation.habitats.sites.geo3.3000m.unique$habitat,relation.habitats.sites.geo3.3000m.unique)

relation.habitats.sites.geo3.5000m.unique<-aggregate(relation.habitats.sites.geo3.5000m$surface,list(Code_Site=relation.habitats.sites.geo3.5000m$Code_Site,habitat=relation.habitats.sites.geo3.5000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo3.5000m.unique)=c("Code_Site","habitat","surface")
tablo_geo3_5000m<-xtabs(formula=relation.habitats.sites.geo3.5000m.unique$surface~as.character(relation.habitats.sites.geo3.5000m.unique$Code_Site)+relation.habitats.sites.geo3.5000m.unique$habitat,relation.habitats.sites.geo3.5000m.unique)

#Load de l'algorithme shannon
source("script.shannon.r")

#tirage des sites poissons en optimisant la diversité shannon d'habitat
 selection.shannon.hab.geo3.100m<-shannon(tablo_geo3_100m)
 selection.shannon.hab.geo3.250m<-shannon(tablo_geo3_250m)
 selection.shannon.hab.geo3.500m<-shannon(tablo_geo3_500m)
 selection.shannon.hab.geo3.750m<-shannon(tablo_geo3_750m)
 selection.shannon.hab.geo3.1000m<-shannon(tablo_geo3_1000m)
 selection.shannon.hab.geo3.3000m<-shannon(tablo_geo3_3000m)
 selection.shannon.hab.geo3.5000m<-shannon(tablo_geo3_5000m)

  #DIVERSITE SELECTIONNE PAR LE SCRIPT SHANNON
  source('script.diversite.shannon.selectionnee.fonction.r')
diversite.selection.shannon.hab.geo3.100m<-species.shannon.evol(selection.shannon.hab.geo3.100m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo3.250m<-species.shannon.evol(selection.shannon.hab.geo3.250m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo3.500m<-species.shannon.evol(selection.shannon.hab.geo3.500m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo3.750m<-species.shannon.evol(selection.shannon.hab.geo3.750m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo3.1000m<-species.shannon.evol(selection.shannon.hab.geo3.1000m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo3.3000m<-species.shannon.evol(selection.shannon.hab.geo3.3000m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo3.5000m<-species.shannon.evol(selection.shannon.hab.geo3.5000m,requette.abondance.poissons)

 #CALCUL DES SAI#
source('script.sai.r')
SAI_shannon.geo3.100m<-SAI(diversite.selection.shannon.hab.geo3.100m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo3.250m<-SAI(diversite.selection.shannon.hab.geo3.250m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo3.500m<-SAI(diversite.selection.shannon.hab.geo3.500m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo3.750m<-SAI(diversite.selection.shannon.hab.geo3.750m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo3.1000m<-SAI(diversite.selection.shannon.hab.geo3.1000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo3.3000m<-SAI(diversite.selection.shannon.hab.geo3.3000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo3.5000m<-SAI(diversite.selection.shannon.hab.geo3.5000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)




  #############Habitat bent1###################



   ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances

setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
relation.habitats.sites.bent1<-read.csv2("fc3_intersect_bent1.csv")

sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')
sites.etudies<-sites.etudies[!sites.etudies%in%suppressed.station]

relation.habitats.sites.bent1<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$Code_Site%in%sites.etudies,]

#création d'une nouvelle table site-habitat par distance
relation.habitats.sites.bent1.30m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==30,]
relation.habitats.sites.bent1.60m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==60,]
relation.habitats.sites.bent1.100m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==100,]
relation.habitats.sites.bent1.250m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==250,]
relation.habitats.sites.bent1.500m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==500,]
relation.habitats.sites.bent1.750m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==750,]
relation.habitats.sites.bent1.1000m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==1000,]
relation.habitats.sites.bent1.3000m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==3000,]
relation.habitats.sites.bent1.5000m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==5000,]

#cration d'un tableau adapté à l'entrée du script shannon
relation.habitats.sites.bent1.30m.unique<-aggregate(relation.habitats.sites.bent1.30m$surface,list(Code_Site=relation.habitats.sites.bent1.30m$Code_Site,habitat=relation.habitats.sites.bent1.30m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent1.30m.unique)=c("Code_Site","habitat","surface")
tablo_bent1_30m<-xtabs(formula=relation.habitats.sites.bent1.30m.unique$surface~as.character(relation.habitats.sites.bent1.30m.unique$Code_Site)+relation.habitats.sites.bent1.30m.unique$habitat,relation.habitats.sites.bent1.30m.unique)

relation.habitats.sites.bent1.60m.unique<-aggregate(relation.habitats.sites.bent1.60m$surface,list(Code_Site=relation.habitats.sites.bent1.60m$Code_Site,habitat=relation.habitats.sites.bent1.60m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent1.60m.unique)=c("Code_Site","habitat","surface")
tablo_bent1_60m<-xtabs(formula=relation.habitats.sites.bent1.60m.unique$surface~as.character(relation.habitats.sites.bent1.60m.unique$Code_Site)+relation.habitats.sites.bent1.60m.unique$habitat,relation.habitats.sites.bent1.60m.unique)

relation.habitats.sites.bent1.100m.unique<-aggregate(relation.habitats.sites.bent1.100m$surface,list(Code_Site=relation.habitats.sites.bent1.100m$Code_Site,habitat=relation.habitats.sites.bent1.100m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent1.100m.unique)=c("Code_Site","habitat","surface")
tablo_bent1_100m<-xtabs(formula=relation.habitats.sites.bent1.100m.unique$surface~as.character(relation.habitats.sites.bent1.100m.unique$Code_Site)+relation.habitats.sites.bent1.100m.unique$habitat,relation.habitats.sites.bent1.100m.unique)

relation.habitats.sites.bent1.250m.unique<-aggregate(relation.habitats.sites.bent1.250m$surface,list(Code_Site=relation.habitats.sites.bent1.250m$Code_Site,habitat=relation.habitats.sites.bent1.250m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent1.250m.unique)=c("Code_Site","habitat","surface")
tablo_bent1_250m<-xtabs(formula=relation.habitats.sites.bent1.250m.unique$surface~as.character(relation.habitats.sites.bent1.250m.unique$Code_Site)+relation.habitats.sites.bent1.250m.unique$habitat,relation.habitats.sites.bent1.250m.unique)

relation.habitats.sites.bent1.500m.unique<-aggregate(relation.habitats.sites.bent1.500m$surface,list(Code_Site=relation.habitats.sites.bent1.500m$Code_Site,habitat=relation.habitats.sites.bent1.500m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent1.500m.unique)=c("Code_Site","habitat","surface")
tablo_bent1_500m<-xtabs(formula=relation.habitats.sites.bent1.500m.unique$surface~as.character(relation.habitats.sites.bent1.500m.unique$Code_Site)+relation.habitats.sites.bent1.500m.unique$habitat,relation.habitats.sites.bent1.500m.unique)

relation.habitats.sites.bent1.750m.unique<-aggregate(relation.habitats.sites.bent1.750m$surface,list(Code_Site=relation.habitats.sites.bent1.750m$Code_Site,habitat=relation.habitats.sites.bent1.750m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent1.750m.unique)=c("Code_Site","habitat","surface")
tablo_bent1_750m<-xtabs(formula=relation.habitats.sites.bent1.750m.unique$surface~as.character(relation.habitats.sites.bent1.750m.unique$Code_Site)+relation.habitats.sites.bent1.750m.unique$habitat,relation.habitats.sites.bent1.750m.unique)

relation.habitats.sites.bent1.1000m.unique<-aggregate(relation.habitats.sites.bent1.1000m$surface,list(Code_Site=relation.habitats.sites.bent1.1000m$Code_Site,habitat=relation.habitats.sites.bent1.1000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent1.1000m.unique)=c("Code_Site","habitat","surface")
tablo_bent1_1000m<-xtabs(formula=relation.habitats.sites.bent1.1000m.unique$surface~as.character(relation.habitats.sites.bent1.1000m.unique$Code_Site)+relation.habitats.sites.bent1.1000m.unique$habitat,relation.habitats.sites.bent1.1000m.unique)

relation.habitats.sites.bent1.3000m.unique<-aggregate(relation.habitats.sites.bent1.3000m$surface,list(Code_Site=relation.habitats.sites.bent1.3000m$Code_Site,habitat=relation.habitats.sites.bent1.3000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent1.3000m.unique)=c("Code_Site","habitat","surface")
tablo_bent1_3000m<-xtabs(formula=relation.habitats.sites.bent1.3000m.unique$surface~as.character(relation.habitats.sites.bent1.3000m.unique$Code_Site)+relation.habitats.sites.bent1.3000m.unique$habitat,relation.habitats.sites.bent1.3000m.unique)

relation.habitats.sites.bent1.5000m.unique<-aggregate(relation.habitats.sites.bent1.5000m$surface,list(Code_Site=relation.habitats.sites.bent1.5000m$Code_Site,habitat=relation.habitats.sites.bent1.5000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent1.5000m.unique)=c("Code_Site","habitat","surface")
tablo_bent1_5000m<-xtabs(formula=relation.habitats.sites.bent1.5000m.unique$surface~as.character(relation.habitats.sites.bent1.5000m.unique$Code_Site)+relation.habitats.sites.bent1.5000m.unique$habitat,relation.habitats.sites.bent1.5000m.unique)

#Load de l'algorithme shannon
source("script.shannon.r")

#tirage des sites poissons en optimisant la diversité shannon d'habitat
 selection.shannon.hab.bent1.30m<-shannon(tablo_bent1_30m)
 selection.shannon.hab.bent1.60m<-shannon(tablo_bent1_60m)
 selection.shannon.hab.bent1.100m<-shannon(tablo_bent1_100m)
 selection.shannon.hab.bent1.250m<-shannon(tablo_bent1_250m)
 selection.shannon.hab.bent1.500m<-shannon(tablo_bent1_500m)
 selection.shannon.hab.bent1.750m<-shannon(tablo_bent1_750m)
 selection.shannon.hab.bent1.1000m<-shannon(tablo_bent1_1000m)
 selection.shannon.hab.bent1.3000m<-shannon(tablo_bent1_3000m)
 selection.shannon.hab.bent1.5000m<-shannon(tablo_bent1_5000m)

  #DIVERSITE SELECTIONNE PAR LE SCRIPT SHANNON
  source('script.diversite.shannon.selectionnee.fonction.r')
diversite.selection.shannon.hab.bent1.30m<-species.shannon.evol(selection.shannon.hab.bent1.30m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent1.60m<-species.shannon.evol(selection.shannon.hab.bent1.60m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent1.100m<-species.shannon.evol(selection.shannon.hab.bent1.100m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent1.250m<-species.shannon.evol(selection.shannon.hab.bent1.250m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent1.500m<-species.shannon.evol(selection.shannon.hab.bent1.500m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent1.750m<-species.shannon.evol(selection.shannon.hab.bent1.750m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent1.1000m<-species.shannon.evol(selection.shannon.hab.bent1.1000m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent1.3000m<-species.shannon.evol(selection.shannon.hab.bent1.3000m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent1.5000m<-species.shannon.evol(selection.shannon.hab.bent1.5000m,requette.abondance.poissons)

 #CALCUL DES SAI#
source('script.sai.r')
SAI_shannon.bent1.30m<-SAI(diversite.selection.shannon.hab.bent1.30m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent1.60m<-SAI(diversite.selection.shannon.hab.bent1.60m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent1.100m<-SAI(diversite.selection.shannon.hab.bent1.100m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent1.250m<-SAI(diversite.selection.shannon.hab.bent1.250m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent1.500m<-SAI(diversite.selection.shannon.hab.bent1.500m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent1.750m<-SAI(diversite.selection.shannon.hab.bent1.750m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent1.1000m<-SAI(diversite.selection.shannon.hab.bent1.1000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent1.3000m<-SAI(diversite.selection.shannon.hab.bent1.3000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent1.5000m<-SAI(diversite.selection.shannon.hab.bent1.5000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)




  #############Habitat bent2###################



   ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances

setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
relation.habitats.sites.bent2<-read.csv2("fc3_intersect_bent2.csv")

sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')
sites.etudies<-sites.etudies[!sites.etudies%in%suppressed.station]

relation.habitats.sites.bent2<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$Code_Site%in%sites.etudies,]

#création d'une nouvelle table site-habitat par distance
relation.habitats.sites.bent2.30m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==30,]
relation.habitats.sites.bent2.60m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==60,]
relation.habitats.sites.bent2.100m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==100,]
relation.habitats.sites.bent2.250m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==250,]
relation.habitats.sites.bent2.500m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==500,]
relation.habitats.sites.bent2.750m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==750,]
relation.habitats.sites.bent2.1000m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==1000,]
relation.habitats.sites.bent2.3000m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==3000,]
relation.habitats.sites.bent2.5000m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==5000,]

#cration d'un tableau adapté à l'entrée du script shannon
relation.habitats.sites.bent2.30m.unique<-aggregate(relation.habitats.sites.bent2.30m$surface,list(Code_Site=relation.habitats.sites.bent2.30m$Code_Site,habitat=relation.habitats.sites.bent2.30m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent2.30m.unique)=c("Code_Site","habitat","surface")
tablo_bent2_30m<-xtabs(formula=relation.habitats.sites.bent2.30m.unique$surface~as.character(relation.habitats.sites.bent2.30m.unique$Code_Site)+relation.habitats.sites.bent2.30m.unique$habitat,relation.habitats.sites.bent2.30m.unique)

relation.habitats.sites.bent2.60m.unique<-aggregate(relation.habitats.sites.bent2.60m$surface,list(Code_Site=relation.habitats.sites.bent2.60m$Code_Site,habitat=relation.habitats.sites.bent2.60m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent2.60m.unique)=c("Code_Site","habitat","surface")
tablo_bent2_60m<-xtabs(formula=relation.habitats.sites.bent2.60m.unique$surface~as.character(relation.habitats.sites.bent2.60m.unique$Code_Site)+relation.habitats.sites.bent2.60m.unique$habitat,relation.habitats.sites.bent2.60m.unique)

relation.habitats.sites.bent2.100m.unique<-aggregate(relation.habitats.sites.bent2.100m$surface,list(Code_Site=relation.habitats.sites.bent2.100m$Code_Site,habitat=relation.habitats.sites.bent2.100m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent2.100m.unique)=c("Code_Site","habitat","surface")
tablo_bent2_100m<-xtabs(formula=relation.habitats.sites.bent2.100m.unique$surface~as.character(relation.habitats.sites.bent2.100m.unique$Code_Site)+relation.habitats.sites.bent2.100m.unique$habitat,relation.habitats.sites.bent2.100m.unique)

relation.habitats.sites.bent2.250m.unique<-aggregate(relation.habitats.sites.bent2.250m$surface,list(Code_Site=relation.habitats.sites.bent2.250m$Code_Site,habitat=relation.habitats.sites.bent2.250m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent2.250m.unique)=c("Code_Site","habitat","surface")
tablo_bent2_250m<-xtabs(formula=relation.habitats.sites.bent2.250m.unique$surface~as.character(relation.habitats.sites.bent2.250m.unique$Code_Site)+relation.habitats.sites.bent2.250m.unique$habitat,relation.habitats.sites.bent2.250m.unique)

relation.habitats.sites.bent2.500m.unique<-aggregate(relation.habitats.sites.bent2.500m$surface,list(Code_Site=relation.habitats.sites.bent2.500m$Code_Site,habitat=relation.habitats.sites.bent2.500m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent2.500m.unique)=c("Code_Site","habitat","surface")
tablo_bent2_500m<-xtabs(formula=relation.habitats.sites.bent2.500m.unique$surface~as.character(relation.habitats.sites.bent2.500m.unique$Code_Site)+relation.habitats.sites.bent2.500m.unique$habitat,relation.habitats.sites.bent2.500m.unique)

relation.habitats.sites.bent2.750m.unique<-aggregate(relation.habitats.sites.bent2.750m$surface,list(Code_Site=relation.habitats.sites.bent2.750m$Code_Site,habitat=relation.habitats.sites.bent2.750m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent2.750m.unique)=c("Code_Site","habitat","surface")
tablo_bent2_750m<-xtabs(formula=relation.habitats.sites.bent2.750m.unique$surface~as.character(relation.habitats.sites.bent2.750m.unique$Code_Site)+relation.habitats.sites.bent2.750m.unique$habitat,relation.habitats.sites.bent2.750m.unique)

relation.habitats.sites.bent2.1000m.unique<-aggregate(relation.habitats.sites.bent2.1000m$surface,list(Code_Site=relation.habitats.sites.bent2.1000m$Code_Site,habitat=relation.habitats.sites.bent2.1000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent2.1000m.unique)=c("Code_Site","habitat","surface")
tablo_bent2_1000m<-xtabs(formula=relation.habitats.sites.bent2.1000m.unique$surface~as.character(relation.habitats.sites.bent2.1000m.unique$Code_Site)+relation.habitats.sites.bent2.1000m.unique$habitat,relation.habitats.sites.bent2.1000m.unique)

relation.habitats.sites.bent2.3000m.unique<-aggregate(relation.habitats.sites.bent2.3000m$surface,list(Code_Site=relation.habitats.sites.bent2.3000m$Code_Site,habitat=relation.habitats.sites.bent2.3000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent2.3000m.unique)=c("Code_Site","habitat","surface")
tablo_bent2_3000m<-xtabs(formula=relation.habitats.sites.bent2.3000m.unique$surface~as.character(relation.habitats.sites.bent2.3000m.unique$Code_Site)+relation.habitats.sites.bent2.3000m.unique$habitat,relation.habitats.sites.bent2.3000m.unique)

relation.habitats.sites.bent2.5000m.unique<-aggregate(relation.habitats.sites.bent2.5000m$surface,list(Code_Site=relation.habitats.sites.bent2.5000m$Code_Site,habitat=relation.habitats.sites.bent2.5000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent2.5000m.unique)=c("Code_Site","habitat","surface")
tablo_bent2_5000m<-xtabs(formula=relation.habitats.sites.bent2.5000m.unique$surface~as.character(relation.habitats.sites.bent2.5000m.unique$Code_Site)+relation.habitats.sites.bent2.5000m.unique$habitat,relation.habitats.sites.bent2.5000m.unique)

#Load de l'algorithme shannon
source("script.shannon.r")

#tirage des sites poissons en optimisant la diversité shannon d'habitat
 selection.shannon.hab.bent2.30m<-shannon(tablo_bent2_30m)
 selection.shannon.hab.bent2.60m<-shannon(tablo_bent2_60m)
 selection.shannon.hab.bent2.100m<-shannon(tablo_bent2_100m)
 selection.shannon.hab.bent2.250m<-shannon(tablo_bent2_250m)
 selection.shannon.hab.bent2.500m<-shannon(tablo_bent2_500m)
 selection.shannon.hab.bent2.750m<-shannon(tablo_bent2_750m)
 selection.shannon.hab.bent2.1000m<-shannon(tablo_bent2_1000m)
 selection.shannon.hab.bent2.3000m<-shannon(tablo_bent2_3000m)
 selection.shannon.hab.bent2.5000m<-shannon(tablo_bent2_5000m)

  #DIVERSITE SELECTIONNE PAR LE SCRIPT SHANNON
  source('script.diversite.shannon.selectionnee.fonction.r')
diversite.selection.shannon.hab.bent2.30m<-species.shannon.evol(selection.shannon.hab.bent2.30m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent2.60m<-species.shannon.evol(selection.shannon.hab.bent2.60m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent2.100m<-species.shannon.evol(selection.shannon.hab.bent2.100m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent2.250m<-species.shannon.evol(selection.shannon.hab.bent2.250m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent2.500m<-species.shannon.evol(selection.shannon.hab.bent2.500m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent2.750m<-species.shannon.evol(selection.shannon.hab.bent2.750m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent2.1000m<-species.shannon.evol(selection.shannon.hab.bent2.1000m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent2.3000m<-species.shannon.evol(selection.shannon.hab.bent2.3000m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent2.5000m<-species.shannon.evol(selection.shannon.hab.bent2.5000m,requette.abondance.poissons)

 #CALCUL DES SAI#
source('script.sai.r')
SAI_shannon.bent2.30m<-SAI(diversite.selection.shannon.hab.bent2.30m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent2.60m<-SAI(diversite.selection.shannon.hab.bent2.60m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent2.100m<-SAI(diversite.selection.shannon.hab.bent2.100m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent2.250m<-SAI(diversite.selection.shannon.hab.bent2.250m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent2.500m<-SAI(diversite.selection.shannon.hab.bent2.500m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent2.750m<-SAI(diversite.selection.shannon.hab.bent2.750m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent2.1000m<-SAI(diversite.selection.shannon.hab.bent2.1000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent2.3000m<-SAI(diversite.selection.shannon.hab.bent2.3000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent2.5000m<-SAI(diversite.selection.shannon.hab.bent2.5000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)




  #############Habitat bent3###################



   ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances

setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
relation.habitats.sites.bent3<-read.csv2("fc3_intersect_bent3.csv")

sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')
sites.etudies<-sites.etudies[!sites.etudies%in%suppressed.station]

relation.habitats.sites.bent3<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$Code_Site%in%sites.etudies,]

#création d'une nouvelle table site-habitat par distance
relation.habitats.sites.bent3.30m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==30,]
relation.habitats.sites.bent3.60m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==60,]
relation.habitats.sites.bent3.100m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==100,]
relation.habitats.sites.bent3.250m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==250,]
relation.habitats.sites.bent3.500m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==500,]
relation.habitats.sites.bent3.750m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==750,]
relation.habitats.sites.bent3.1000m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==1000,]
relation.habitats.sites.bent3.3000m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==3000,]
relation.habitats.sites.bent3.5000m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==5000,]

#cration d'un tableau adapté à l'entrée du script shannon
relation.habitats.sites.bent3.30m.unique<-aggregate(relation.habitats.sites.bent3.30m$surface,list(Code_Site=relation.habitats.sites.bent3.30m$Code_Site,habitat=relation.habitats.sites.bent3.30m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent3.30m.unique)=c("Code_Site","habitat","surface")
tablo_bent3_30m<-xtabs(formula=relation.habitats.sites.bent3.30m.unique$surface~as.character(relation.habitats.sites.bent3.30m.unique$Code_Site)+relation.habitats.sites.bent3.30m.unique$habitat,relation.habitats.sites.bent3.30m.unique)

relation.habitats.sites.bent3.60m.unique<-aggregate(relation.habitats.sites.bent3.60m$surface,list(Code_Site=relation.habitats.sites.bent3.60m$Code_Site,habitat=relation.habitats.sites.bent3.60m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent3.60m.unique)=c("Code_Site","habitat","surface")
tablo_bent3_60m<-xtabs(formula=relation.habitats.sites.bent3.60m.unique$surface~as.character(relation.habitats.sites.bent3.60m.unique$Code_Site)+relation.habitats.sites.bent3.60m.unique$habitat,relation.habitats.sites.bent3.60m.unique)

relation.habitats.sites.bent3.100m.unique<-aggregate(relation.habitats.sites.bent3.100m$surface,list(Code_Site=relation.habitats.sites.bent3.100m$Code_Site,habitat=relation.habitats.sites.bent3.100m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent3.100m.unique)=c("Code_Site","habitat","surface")
tablo_bent3_100m<-xtabs(formula=relation.habitats.sites.bent3.100m.unique$surface~as.character(relation.habitats.sites.bent3.100m.unique$Code_Site)+relation.habitats.sites.bent3.100m.unique$habitat,relation.habitats.sites.bent3.100m.unique)

relation.habitats.sites.bent3.250m.unique<-aggregate(relation.habitats.sites.bent3.250m$surface,list(Code_Site=relation.habitats.sites.bent3.250m$Code_Site,habitat=relation.habitats.sites.bent3.250m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent3.250m.unique)=c("Code_Site","habitat","surface")
tablo_bent3_250m<-xtabs(formula=relation.habitats.sites.bent3.250m.unique$surface~as.character(relation.habitats.sites.bent3.250m.unique$Code_Site)+relation.habitats.sites.bent3.250m.unique$habitat,relation.habitats.sites.bent3.250m.unique)

relation.habitats.sites.bent3.500m.unique<-aggregate(relation.habitats.sites.bent3.500m$surface,list(Code_Site=relation.habitats.sites.bent3.500m$Code_Site,habitat=relation.habitats.sites.bent3.500m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent3.500m.unique)=c("Code_Site","habitat","surface")
tablo_bent3_500m<-xtabs(formula=relation.habitats.sites.bent3.500m.unique$surface~as.character(relation.habitats.sites.bent3.500m.unique$Code_Site)+relation.habitats.sites.bent3.500m.unique$habitat,relation.habitats.sites.bent3.500m.unique)

relation.habitats.sites.bent3.750m.unique<-aggregate(relation.habitats.sites.bent3.750m$surface,list(Code_Site=relation.habitats.sites.bent3.750m$Code_Site,habitat=relation.habitats.sites.bent3.750m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent3.750m.unique)=c("Code_Site","habitat","surface")
tablo_bent3_750m<-xtabs(formula=relation.habitats.sites.bent3.750m.unique$surface~as.character(relation.habitats.sites.bent3.750m.unique$Code_Site)+relation.habitats.sites.bent3.750m.unique$habitat,relation.habitats.sites.bent3.750m.unique)

relation.habitats.sites.bent3.1000m.unique<-aggregate(relation.habitats.sites.bent3.1000m$surface,list(Code_Site=relation.habitats.sites.bent3.1000m$Code_Site,habitat=relation.habitats.sites.bent3.1000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent3.1000m.unique)=c("Code_Site","habitat","surface")
tablo_bent3_1000m<-xtabs(formula=relation.habitats.sites.bent3.1000m.unique$surface~as.character(relation.habitats.sites.bent3.1000m.unique$Code_Site)+relation.habitats.sites.bent3.1000m.unique$habitat,relation.habitats.sites.bent3.1000m.unique)

relation.habitats.sites.bent3.3000m.unique<-aggregate(relation.habitats.sites.bent3.3000m$surface,list(Code_Site=relation.habitats.sites.bent3.3000m$Code_Site,habitat=relation.habitats.sites.bent3.3000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent3.3000m.unique)=c("Code_Site","habitat","surface")
tablo_bent3_3000m<-xtabs(formula=relation.habitats.sites.bent3.3000m.unique$surface~as.character(relation.habitats.sites.bent3.3000m.unique$Code_Site)+relation.habitats.sites.bent3.3000m.unique$habitat,relation.habitats.sites.bent3.3000m.unique)

relation.habitats.sites.bent3.5000m.unique<-aggregate(relation.habitats.sites.bent3.5000m$surface,list(Code_Site=relation.habitats.sites.bent3.5000m$Code_Site,habitat=relation.habitats.sites.bent3.5000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.bent3.5000m.unique)=c("Code_Site","habitat","surface")
tablo_bent3_5000m<-xtabs(formula=relation.habitats.sites.bent3.5000m.unique$surface~as.character(relation.habitats.sites.bent3.5000m.unique$Code_Site)+relation.habitats.sites.bent3.5000m.unique$habitat,relation.habitats.sites.bent3.5000m.unique)

#Load de l'algorithme shannon
source("script.shannon.r")

#tirage des sites poissons en optimisant la diversité shannon d'habitat
 selection.shannon.hab.bent3.30m<-shannon(tablo_bent3_30m)
 selection.shannon.hab.bent3.60m<-shannon(tablo_bent3_60m)
 selection.shannon.hab.bent3.100m<-shannon(tablo_bent3_100m)
 selection.shannon.hab.bent3.250m<-shannon(tablo_bent3_250m)
 selection.shannon.hab.bent3.500m<-shannon(tablo_bent3_500m)
 selection.shannon.hab.bent3.750m<-shannon(tablo_bent3_750m)
 selection.shannon.hab.bent3.1000m<-shannon(tablo_bent3_1000m)
 selection.shannon.hab.bent3.3000m<-shannon(tablo_bent3_3000m)
 selection.shannon.hab.bent3.5000m<-shannon(tablo_bent3_5000m)

  #DIVERSITE SELECTIONNE PAR LE SCRIPT SHANNON
  source('script.diversite.shannon.selectionnee.fonction.r')
diversite.selection.shannon.hab.bent3.30m<-species.shannon.evol(selection.shannon.hab.bent3.30m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent3.60m<-species.shannon.evol(selection.shannon.hab.bent3.60m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent3.100m<-species.shannon.evol(selection.shannon.hab.bent3.100m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent3.250m<-species.shannon.evol(selection.shannon.hab.bent3.250m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent3.500m<-species.shannon.evol(selection.shannon.hab.bent3.500m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent3.750m<-species.shannon.evol(selection.shannon.hab.bent3.750m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent3.1000m<-species.shannon.evol(selection.shannon.hab.bent3.1000m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent3.3000m<-species.shannon.evol(selection.shannon.hab.bent3.3000m,requette.abondance.poissons)
diversite.selection.shannon.hab.bent3.5000m<-species.shannon.evol(selection.shannon.hab.bent3.5000m,requette.abondance.poissons)

 #CALCUL DES SAI#
source('script.sai.r')
SAI_shannon.bent3.30m<-SAI(diversite.selection.shannon.hab.bent3.30m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent3.60m<-SAI(diversite.selection.shannon.hab.bent3.60m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent3.100m<-SAI(diversite.selection.shannon.hab.bent3.100m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent3.250m<-SAI(diversite.selection.shannon.hab.bent3.250m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent3.500m<-SAI(diversite.selection.shannon.hab.bent3.500m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent3.750m<-SAI(diversite.selection.shannon.hab.bent3.750m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent3.1000m<-SAI(diversite.selection.shannon.hab.bent3.1000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent3.3000m<-SAI(diversite.selection.shannon.hab.bent3.3000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.bent3.5000m<-SAI(diversite.selection.shannon.hab.bent3.5000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)




  #############Habitat geo1_geo2###################



   ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances

setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
relation.habitats.sites.geo1_geo2<-read.csv2("fc3_intersect_geo1_geo2.csv")

sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')
sites.etudies<-sites.etudies[!sites.etudies%in%suppressed.station]

relation.habitats.sites.geo1_geo2<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$Code_Site%in%sites.etudies,]

#création d'une nouvelle table site-habitat par distance
relation.habitats.sites.geo1_geo2.30m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==30,]
relation.habitats.sites.geo1_geo2.60m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==60,]
relation.habitats.sites.geo1_geo2.100m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==100,]
relation.habitats.sites.geo1_geo2.250m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==250,]
relation.habitats.sites.geo1_geo2.500m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==500,]
relation.habitats.sites.geo1_geo2.750m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==750,]
relation.habitats.sites.geo1_geo2.1000m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==1000,]
relation.habitats.sites.geo1_geo2.3000m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==3000,]
relation.habitats.sites.geo1_geo2.5000m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==5000,]

#cration d'un tableau adapté à l'entrée du script shannon
relation.habitats.sites.geo1_geo2.30m.unique<-aggregate(relation.habitats.sites.geo1_geo2.30m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2.30m$Code_Site,habitat=relation.habitats.sites.geo1_geo2.30m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2.30m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_30m<-xtabs(formula=relation.habitats.sites.geo1_geo2.30m.unique$surface~as.character(relation.habitats.sites.geo1_geo2.30m.unique$Code_Site)+relation.habitats.sites.geo1_geo2.30m.unique$habitat,relation.habitats.sites.geo1_geo2.30m.unique)

relation.habitats.sites.geo1_geo2.60m.unique<-aggregate(relation.habitats.sites.geo1_geo2.60m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2.60m$Code_Site,habitat=relation.habitats.sites.geo1_geo2.60m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2.60m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_60m<-xtabs(formula=relation.habitats.sites.geo1_geo2.60m.unique$surface~as.character(relation.habitats.sites.geo1_geo2.60m.unique$Code_Site)+relation.habitats.sites.geo1_geo2.60m.unique$habitat,relation.habitats.sites.geo1_geo2.60m.unique)

relation.habitats.sites.geo1_geo2.100m.unique<-aggregate(relation.habitats.sites.geo1_geo2.100m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2.100m$Code_Site,habitat=relation.habitats.sites.geo1_geo2.100m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2.100m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_100m<-xtabs(formula=relation.habitats.sites.geo1_geo2.100m.unique$surface~as.character(relation.habitats.sites.geo1_geo2.100m.unique$Code_Site)+relation.habitats.sites.geo1_geo2.100m.unique$habitat,relation.habitats.sites.geo1_geo2.100m.unique)

relation.habitats.sites.geo1_geo2.250m.unique<-aggregate(relation.habitats.sites.geo1_geo2.250m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2.250m$Code_Site,habitat=relation.habitats.sites.geo1_geo2.250m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2.250m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_250m<-xtabs(formula=relation.habitats.sites.geo1_geo2.250m.unique$surface~as.character(relation.habitats.sites.geo1_geo2.250m.unique$Code_Site)+relation.habitats.sites.geo1_geo2.250m.unique$habitat,relation.habitats.sites.geo1_geo2.250m.unique)

relation.habitats.sites.geo1_geo2.500m.unique<-aggregate(relation.habitats.sites.geo1_geo2.500m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2.500m$Code_Site,habitat=relation.habitats.sites.geo1_geo2.500m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2.500m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_500m<-xtabs(formula=relation.habitats.sites.geo1_geo2.500m.unique$surface~as.character(relation.habitats.sites.geo1_geo2.500m.unique$Code_Site)+relation.habitats.sites.geo1_geo2.500m.unique$habitat,relation.habitats.sites.geo1_geo2.500m.unique)

relation.habitats.sites.geo1_geo2.750m.unique<-aggregate(relation.habitats.sites.geo1_geo2.750m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2.750m$Code_Site,habitat=relation.habitats.sites.geo1_geo2.750m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2.750m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_750m<-xtabs(formula=relation.habitats.sites.geo1_geo2.750m.unique$surface~as.character(relation.habitats.sites.geo1_geo2.750m.unique$Code_Site)+relation.habitats.sites.geo1_geo2.750m.unique$habitat,relation.habitats.sites.geo1_geo2.750m.unique)

relation.habitats.sites.geo1_geo2.1000m.unique<-aggregate(relation.habitats.sites.geo1_geo2.1000m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2.1000m$Code_Site,habitat=relation.habitats.sites.geo1_geo2.1000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2.1000m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_1000m<-xtabs(formula=relation.habitats.sites.geo1_geo2.1000m.unique$surface~as.character(relation.habitats.sites.geo1_geo2.1000m.unique$Code_Site)+relation.habitats.sites.geo1_geo2.1000m.unique$habitat,relation.habitats.sites.geo1_geo2.1000m.unique)

relation.habitats.sites.geo1_geo2.3000m.unique<-aggregate(relation.habitats.sites.geo1_geo2.3000m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2.3000m$Code_Site,habitat=relation.habitats.sites.geo1_geo2.3000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2.3000m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_3000m<-xtabs(formula=relation.habitats.sites.geo1_geo2.3000m.unique$surface~as.character(relation.habitats.sites.geo1_geo2.3000m.unique$Code_Site)+relation.habitats.sites.geo1_geo2.3000m.unique$habitat,relation.habitats.sites.geo1_geo2.3000m.unique)

relation.habitats.sites.geo1_geo2.5000m.unique<-aggregate(relation.habitats.sites.geo1_geo2.5000m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2.5000m$Code_Site,habitat=relation.habitats.sites.geo1_geo2.5000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2.5000m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_5000m<-xtabs(formula=relation.habitats.sites.geo1_geo2.5000m.unique$surface~as.character(relation.habitats.sites.geo1_geo2.5000m.unique$Code_Site)+relation.habitats.sites.geo1_geo2.5000m.unique$habitat,relation.habitats.sites.geo1_geo2.5000m.unique)

#Load de l'algorithme shannon
source("script.shannon.r")

#tirage des sites poissons en optimisant la diversité shannon d'habitat
 selection.shannon.hab.geo1_geo2.30m<-shannon(tablo_geo1_geo2_30m)
 selection.shannon.hab.geo1_geo2.60m<-shannon(tablo_geo1_geo2_60m)
 selection.shannon.hab.geo1_geo2.100m<-shannon(tablo_geo1_geo2_100m)
 selection.shannon.hab.geo1_geo2.250m<-shannon(tablo_geo1_geo2_250m)
 selection.shannon.hab.geo1_geo2.500m<-shannon(tablo_geo1_geo2_500m)
 selection.shannon.hab.geo1_geo2.750m<-shannon(tablo_geo1_geo2_750m)
 selection.shannon.hab.geo1_geo2.1000m<-shannon(tablo_geo1_geo2_1000m)
 selection.shannon.hab.geo1_geo2.3000m<-shannon(tablo_geo1_geo2_3000m)
 selection.shannon.hab.geo1_geo2.5000m<-shannon(tablo_geo1_geo2_5000m)

  #DIVERSITE SELECTIONNE PAR LE SCRIPT SHANNON
  source('script.diversite.shannon.selectionnee.fonction.r')
diversite.selection.shannon.hab.geo1_geo2.30m<-species.shannon.evol(selection.shannon.hab.geo1_geo2.30m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1_geo2.60m<-species.shannon.evol(selection.shannon.hab.geo1_geo2.60m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1_geo2.100m<-species.shannon.evol(selection.shannon.hab.geo1_geo2.100m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1_geo2.250m<-species.shannon.evol(selection.shannon.hab.geo1_geo2.250m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1_geo2.500m<-species.shannon.evol(selection.shannon.hab.geo1_geo2.500m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1_geo2.750m<-species.shannon.evol(selection.shannon.hab.geo1_geo2.750m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1_geo2.1000m<-species.shannon.evol(selection.shannon.hab.geo1_geo2.1000m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1_geo2.3000m<-species.shannon.evol(selection.shannon.hab.geo1_geo2.3000m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1_geo2.5000m<-species.shannon.evol(selection.shannon.hab.geo1_geo2.5000m,requette.abondance.poissons)

 #CALCUL DES SAI#
source('script.sai.r')
SAI_shannon.geo1_geo2.30m<-SAI(diversite.selection.shannon.hab.geo1_geo2.30m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1_geo2.60m<-SAI(diversite.selection.shannon.hab.geo1_geo2.60m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1_geo2.100m<-SAI(diversite.selection.shannon.hab.geo1_geo2.100m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1_geo2.250m<-SAI(diversite.selection.shannon.hab.geo1_geo2.250m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1_geo2.500m<-SAI(diversite.selection.shannon.hab.geo1_geo2.500m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1_geo2.750m<-SAI(diversite.selection.shannon.hab.geo1_geo2.750m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1_geo2.1000m<-SAI(diversite.selection.shannon.hab.geo1_geo2.1000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1_geo2.3000m<-SAI(diversite.selection.shannon.hab.geo1_geo2.3000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1_geo2.5000m<-SAI(diversite.selection.shannon.hab.geo1_geo2.5000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)




  #############Habitat geo1_geo2_geo3###################



   ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances

setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
relation.habitats.sites.geo1_geo2_geo3<-read.csv2("fc3_intersect_geo1_geo2_geo3.csv")

sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')
sites.etudies<-sites.etudies[!sites.etudies%in%suppressed.station]

relation.habitats.sites.geo1_geo2_geo3<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$Code_Site%in%sites.etudies,]

#création d'une nouvelle table site-habitat par distance
relation.habitats.sites.geo1_geo2_geo3.30m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==30,]
relation.habitats.sites.geo1_geo2_geo3.60m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==60,]
relation.habitats.sites.geo1_geo2_geo3.100m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==100,]
relation.habitats.sites.geo1_geo2_geo3.250m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==250,]
relation.habitats.sites.geo1_geo2_geo3.500m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==500,]
relation.habitats.sites.geo1_geo2_geo3.750m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==750,]
relation.habitats.sites.geo1_geo2_geo3.1000m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==1000,]
relation.habitats.sites.geo1_geo2_geo3.3000m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==3000,]
relation.habitats.sites.geo1_geo2_geo3.5000m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==5000,]

#cration d'un tableau adapté à l'entrée du script shannon
relation.habitats.sites.geo1_geo2_geo3.30m.unique<-aggregate(relation.habitats.sites.geo1_geo2_geo3.30m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2_geo3.30m$Code_Site,habitat=relation.habitats.sites.geo1_geo2_geo3.30m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2_geo3.30m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_geo3_30m<-xtabs(formula=relation.habitats.sites.geo1_geo2_geo3.30m.unique$surface~as.character(relation.habitats.sites.geo1_geo2_geo3.30m.unique$Code_Site)+relation.habitats.sites.geo1_geo2_geo3.30m.unique$habitat,relation.habitats.sites.geo1_geo2_geo3.30m.unique)

relation.habitats.sites.geo1_geo2_geo3.60m.unique<-aggregate(relation.habitats.sites.geo1_geo2_geo3.60m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2_geo3.60m$Code_Site,habitat=relation.habitats.sites.geo1_geo2_geo3.60m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2_geo3.60m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_geo3_60m<-xtabs(formula=relation.habitats.sites.geo1_geo2_geo3.60m.unique$surface~as.character(relation.habitats.sites.geo1_geo2_geo3.60m.unique$Code_Site)+relation.habitats.sites.geo1_geo2_geo3.60m.unique$habitat,relation.habitats.sites.geo1_geo2_geo3.60m.unique)

relation.habitats.sites.geo1_geo2_geo3.100m.unique<-aggregate(relation.habitats.sites.geo1_geo2_geo3.100m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2_geo3.100m$Code_Site,habitat=relation.habitats.sites.geo1_geo2_geo3.100m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2_geo3.100m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_geo3_100m<-xtabs(formula=relation.habitats.sites.geo1_geo2_geo3.100m.unique$surface~as.character(relation.habitats.sites.geo1_geo2_geo3.100m.unique$Code_Site)+relation.habitats.sites.geo1_geo2_geo3.100m.unique$habitat,relation.habitats.sites.geo1_geo2_geo3.100m.unique)

relation.habitats.sites.geo1_geo2_geo3.250m.unique<-aggregate(relation.habitats.sites.geo1_geo2_geo3.250m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2_geo3.250m$Code_Site,habitat=relation.habitats.sites.geo1_geo2_geo3.250m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2_geo3.250m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_geo3_250m<-xtabs(formula=relation.habitats.sites.geo1_geo2_geo3.250m.unique$surface~as.character(relation.habitats.sites.geo1_geo2_geo3.250m.unique$Code_Site)+relation.habitats.sites.geo1_geo2_geo3.250m.unique$habitat,relation.habitats.sites.geo1_geo2_geo3.250m.unique)

relation.habitats.sites.geo1_geo2_geo3.500m.unique<-aggregate(relation.habitats.sites.geo1_geo2_geo3.500m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2_geo3.500m$Code_Site,habitat=relation.habitats.sites.geo1_geo2_geo3.500m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2_geo3.500m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_geo3_500m<-xtabs(formula=relation.habitats.sites.geo1_geo2_geo3.500m.unique$surface~as.character(relation.habitats.sites.geo1_geo2_geo3.500m.unique$Code_Site)+relation.habitats.sites.geo1_geo2_geo3.500m.unique$habitat,relation.habitats.sites.geo1_geo2_geo3.500m.unique)

relation.habitats.sites.geo1_geo2_geo3.750m.unique<-aggregate(relation.habitats.sites.geo1_geo2_geo3.750m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2_geo3.750m$Code_Site,habitat=relation.habitats.sites.geo1_geo2_geo3.750m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2_geo3.750m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_geo3_750m<-xtabs(formula=relation.habitats.sites.geo1_geo2_geo3.750m.unique$surface~as.character(relation.habitats.sites.geo1_geo2_geo3.750m.unique$Code_Site)+relation.habitats.sites.geo1_geo2_geo3.750m.unique$habitat,relation.habitats.sites.geo1_geo2_geo3.750m.unique)

relation.habitats.sites.geo1_geo2_geo3.1000m.unique<-aggregate(relation.habitats.sites.geo1_geo2_geo3.1000m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2_geo3.1000m$Code_Site,habitat=relation.habitats.sites.geo1_geo2_geo3.1000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2_geo3.1000m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_geo3_1000m<-xtabs(formula=relation.habitats.sites.geo1_geo2_geo3.1000m.unique$surface~as.character(relation.habitats.sites.geo1_geo2_geo3.1000m.unique$Code_Site)+relation.habitats.sites.geo1_geo2_geo3.1000m.unique$habitat,relation.habitats.sites.geo1_geo2_geo3.1000m.unique)

relation.habitats.sites.geo1_geo2_geo3.3000m.unique<-aggregate(relation.habitats.sites.geo1_geo2_geo3.3000m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2_geo3.3000m$Code_Site,habitat=relation.habitats.sites.geo1_geo2_geo3.3000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2_geo3.3000m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_geo3_3000m<-xtabs(formula=relation.habitats.sites.geo1_geo2_geo3.3000m.unique$surface~as.character(relation.habitats.sites.geo1_geo2_geo3.3000m.unique$Code_Site)+relation.habitats.sites.geo1_geo2_geo3.3000m.unique$habitat,relation.habitats.sites.geo1_geo2_geo3.3000m.unique)

relation.habitats.sites.geo1_geo2_geo3.5000m.unique<-aggregate(relation.habitats.sites.geo1_geo2_geo3.5000m$surface,list(Code_Site=relation.habitats.sites.geo1_geo2_geo3.5000m$Code_Site,habitat=relation.habitats.sites.geo1_geo2_geo3.5000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.geo1_geo2_geo3.5000m.unique)=c("Code_Site","habitat","surface")
tablo_geo1_geo2_geo3_5000m<-xtabs(formula=relation.habitats.sites.geo1_geo2_geo3.5000m.unique$surface~as.character(relation.habitats.sites.geo1_geo2_geo3.5000m.unique$Code_Site)+relation.habitats.sites.geo1_geo2_geo3.5000m.unique$habitat,relation.habitats.sites.geo1_geo2_geo3.5000m.unique)

#Load de l'algorithme shannon
source("script.shannon.r")

#tirage des sites poissons en optimisant la diversité shannon d'habitat
 selection.shannon.hab.geo1_geo2_geo3.30m<-shannon(tablo_geo1_geo2_geo3_30m)
 selection.shannon.hab.geo1_geo2_geo3.60m<-shannon(tablo_geo1_geo2_geo3_60m)
 selection.shannon.hab.geo1_geo2_geo3.100m<-shannon(tablo_geo1_geo2_geo3_100m)
 selection.shannon.hab.geo1_geo2_geo3.250m<-shannon(tablo_geo1_geo2_geo3_250m)
 selection.shannon.hab.geo1_geo2_geo3.500m<-shannon(tablo_geo1_geo2_geo3_500m)
 selection.shannon.hab.geo1_geo2_geo3.750m<-shannon(tablo_geo1_geo2_geo3_750m)
 selection.shannon.hab.geo1_geo2_geo3.1000m<-shannon(tablo_geo1_geo2_geo3_1000m)
 selection.shannon.hab.geo1_geo2_geo3.3000m<-shannon(tablo_geo1_geo2_geo3_3000m)
 selection.shannon.hab.geo1_geo2_geo3.5000m<-shannon(tablo_geo1_geo2_geo3_5000m)

  #DIVERSITE SELECTIONNE PAR LE SCRIPT SHANNON
  source('script.diversite.shannon.selectionnee.fonction.r')
diversite.selection.shannon.hab.geo1_geo2_geo3.30m<-species.shannon.evol(selection.shannon.hab.geo1_geo2_geo3.30m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1_geo2_geo3.60m<-species.shannon.evol(selection.shannon.hab.geo1_geo2_geo3.60m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1_geo2_geo3.100m<-species.shannon.evol(selection.shannon.hab.geo1_geo2_geo3.100m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1_geo2_geo3.250m<-species.shannon.evol(selection.shannon.hab.geo1_geo2_geo3.250m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1_geo2_geo3.500m<-species.shannon.evol(selection.shannon.hab.geo1_geo2_geo3.500m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1_geo2_geo3.750m<-species.shannon.evol(selection.shannon.hab.geo1_geo2_geo3.750m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1_geo2_geo3.1000m<-species.shannon.evol(selection.shannon.hab.geo1_geo2_geo3.1000m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1_geo2_geo3.3000m<-species.shannon.evol(selection.shannon.hab.geo1_geo2_geo3.3000m,requette.abondance.poissons)
diversite.selection.shannon.hab.geo1_geo2_geo3.5000m<-species.shannon.evol(selection.shannon.hab.geo1_geo2_geo3.5000m,requette.abondance.poissons)

 #CALCUL DES SAI#
source('script.sai.r')
SAI_shannon.geo1_geo2_geo3.30m<-SAI(diversite.selection.shannon.hab.geo1_geo2_geo3.30m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1_geo2_geo3.60m<-SAI(diversite.selection.shannon.hab.geo1_geo2_geo3.60m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1_geo2_geo3.100m<-SAI(diversite.selection.shannon.hab.geo1_geo2_geo3.100m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1_geo2_geo3.250m<-SAI(diversite.selection.shannon.hab.geo1_geo2_geo3.250m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1_geo2_geo3.500m<-SAI(diversite.selection.shannon.hab.geo1_geo2_geo3.500m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1_geo2_geo3.750m<-SAI(diversite.selection.shannon.hab.geo1_geo2_geo3.750m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1_geo2_geo3.1000m<-SAI(diversite.selection.shannon.hab.geo1_geo2_geo3.1000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1_geo2_geo3.3000m<-SAI(diversite.selection.shannon.hab.geo1_geo2_geo3.3000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.geo1_geo2_geo3.5000m<-SAI(diversite.selection.shannon.hab.geo1_geo2_geo3.5000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)



  #############Habitat g1_g2_g3_b1###################



   ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances

setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
relation.habitats.sites.g1_g2_g3_b1<-read.csv2("fc3_intersect_g1_g2_g3_b1.csv")

sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')
sites.etudies<-sites.etudies[!sites.etudies%in%suppressed.station]

relation.habitats.sites.g1_g2_g3_b1<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$Code_Site%in%sites.etudies,]

#création d'une nouvelle table site-habitat par distance
relation.habitats.sites.g1_g2_g3_b1.30m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==30,]
relation.habitats.sites.g1_g2_g3_b1.60m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==60,]
relation.habitats.sites.g1_g2_g3_b1.100m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==100,]
relation.habitats.sites.g1_g2_g3_b1.250m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==250,]
relation.habitats.sites.g1_g2_g3_b1.500m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==500,]
relation.habitats.sites.g1_g2_g3_b1.750m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==750,]
relation.habitats.sites.g1_g2_g3_b1.1000m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==1000,]
relation.habitats.sites.g1_g2_g3_b1.3000m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==3000,]
relation.habitats.sites.g1_g2_g3_b1.5000m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==5000,]

#cration d'un tableau adapté à l'entrée du script shannon
relation.habitats.sites.g1_g2_g3_b1.30m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1.30m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1.30m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1.30m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1.30m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_30m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1.30m.unique$surface~as.character(relation.habitats.sites.g1_g2_g3_b1.30m.unique$Code_Site)+relation.habitats.sites.g1_g2_g3_b1.30m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1.30m.unique)

relation.habitats.sites.g1_g2_g3_b1.60m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1.60m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1.60m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1.60m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1.60m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_60m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1.60m.unique$surface~as.character(relation.habitats.sites.g1_g2_g3_b1.60m.unique$Code_Site)+relation.habitats.sites.g1_g2_g3_b1.60m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1.60m.unique)

relation.habitats.sites.g1_g2_g3_b1.100m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1.100m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1.100m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1.100m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1.100m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_100m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1.100m.unique$surface~as.character(relation.habitats.sites.g1_g2_g3_b1.100m.unique$Code_Site)+relation.habitats.sites.g1_g2_g3_b1.100m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1.100m.unique)

relation.habitats.sites.g1_g2_g3_b1.250m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1.250m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1.250m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1.250m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1.250m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_250m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1.250m.unique$surface~as.character(relation.habitats.sites.g1_g2_g3_b1.250m.unique$Code_Site)+relation.habitats.sites.g1_g2_g3_b1.250m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1.250m.unique)

relation.habitats.sites.g1_g2_g3_b1.500m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1.500m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1.500m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1.500m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1.500m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_500m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1.500m.unique$surface~as.character(relation.habitats.sites.g1_g2_g3_b1.500m.unique$Code_Site)+relation.habitats.sites.g1_g2_g3_b1.500m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1.500m.unique)

relation.habitats.sites.g1_g2_g3_b1.750m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1.750m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1.750m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1.750m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1.750m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_750m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1.750m.unique$surface~as.character(relation.habitats.sites.g1_g2_g3_b1.750m.unique$Code_Site)+relation.habitats.sites.g1_g2_g3_b1.750m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1.750m.unique)

relation.habitats.sites.g1_g2_g3_b1.1000m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1.1000m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1.1000m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1.1000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1.1000m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_1000m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1.1000m.unique$surface~as.character(relation.habitats.sites.g1_g2_g3_b1.1000m.unique$Code_Site)+relation.habitats.sites.g1_g2_g3_b1.1000m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1.1000m.unique)

relation.habitats.sites.g1_g2_g3_b1.3000m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1.3000m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1.3000m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1.3000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1.3000m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_3000m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1.3000m.unique$surface~as.character(relation.habitats.sites.g1_g2_g3_b1.3000m.unique$Code_Site)+relation.habitats.sites.g1_g2_g3_b1.3000m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1.3000m.unique)

relation.habitats.sites.g1_g2_g3_b1.5000m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1.5000m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1.5000m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1.5000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1.5000m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_5000m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1.5000m.unique$surface~as.character(relation.habitats.sites.g1_g2_g3_b1.5000m.unique$Code_Site)+relation.habitats.sites.g1_g2_g3_b1.5000m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1.5000m.unique)

#Load de l'algorithme shannon
source("script.shannon.r")

#tirage des sites poissons en optimisant la diversité shannon d'habitat
 selection.shannon.hab.g1_g2_g3_b1.30m<-shannon(tablo_g1_g2_g3_b1_30m)
 selection.shannon.hab.g1_g2_g3_b1.60m<-shannon(tablo_g1_g2_g3_b1_60m)
 selection.shannon.hab.g1_g2_g3_b1.100m<-shannon(tablo_g1_g2_g3_b1_100m)
 selection.shannon.hab.g1_g2_g3_b1.250m<-shannon(tablo_g1_g2_g3_b1_250m)
 selection.shannon.hab.g1_g2_g3_b1.500m<-shannon(tablo_g1_g2_g3_b1_500m)
 selection.shannon.hab.g1_g2_g3_b1.750m<-shannon(tablo_g1_g2_g3_b1_750m)
 selection.shannon.hab.g1_g2_g3_b1.1000m<-shannon(tablo_g1_g2_g3_b1_1000m)
 selection.shannon.hab.g1_g2_g3_b1.3000m<-shannon(tablo_g1_g2_g3_b1_3000m)
 selection.shannon.hab.g1_g2_g3_b1.5000m<-shannon(tablo_g1_g2_g3_b1_5000m)

  #DIVERSITE SELECTIONNE PAR LE SCRIPT SHANNON
  source('script.diversite.shannon.selectionnee.fonction.r')
diversite.selection.shannon.hab.g1_g2_g3_b1.30m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1.30m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1.60m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1.60m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1.100m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1.100m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1.250m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1.250m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1.500m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1.500m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1.750m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1.750m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1.1000m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1.1000m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1.3000m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1.3000m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1.5000m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1.5000m,requette.abondance.poissons)

 #CALCUL DES SAI#
source('script.sai.r')
SAI_shannon.g1_g2_g3_b1.30m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1.30m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1.60m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1.60m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1.100m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1.100m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1.250m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1.250m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1.500m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1.500m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1.750m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1.750m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1.1000m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1.1000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1.3000m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1.3000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1.5000m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1.5000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)




  #############Habitat g1_g2_g3_b1_b2###################



   ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances

setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
relation.habitats.sites.g1_g2_g3_b1_b2<-read.csv2("fc3_intersect_g1_g2_g3_b1_b2.csv")

sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')
sites.etudies<-sites.etudies[!sites.etudies%in%suppressed.station]

relation.habitats.sites.g1_g2_g3_b1_b2<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$Code_Site%in%sites.etudies,]

#création d'une nouvelle table site-habitat par distance
relation.habitats.sites.g1_g2_g3_b1_b2.30m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==30,]
relation.habitats.sites.g1_g2_g3_b1_b2.60m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==60,]
relation.habitats.sites.g1_g2_g3_b1_b2.100m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==100,]
relation.habitats.sites.g1_g2_g3_b1_b2.250m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==250,]
relation.habitats.sites.g1_g2_g3_b1_b2.500m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==500,]
relation.habitats.sites.g1_g2_g3_b1_b2.750m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==750,]
relation.habitats.sites.g1_g2_g3_b1_b2.1000m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==1000,]
relation.habitats.sites.g1_g2_g3_b1_b2.3000m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==3000,]
relation.habitats.sites.g1_g2_g3_b1_b2.5000m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==5000,]

#cration d'un tableau adapté à l'entrée du script shannon
relation.habitats.sites.g1_g2_g3_b1_b2.30m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2.30m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2.30m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2.30m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2.30m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_30m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2.30m.unique$surface~as.character(relation.habitats.sites.g1_g2_g3_b1_b2.30m.unique$Code_Site)+relation.habitats.sites.g1_g2_g3_b1_b2.30m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2.30m.unique)

relation.habitats.sites.g1_g2_g3_b1_b2.60m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2.60m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2.60m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2.60m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2.60m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_60m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2.60m.unique$surface~as.character(relation.habitats.sites.g1_g2_g3_b1_b2.60m.unique$Code_Site)+relation.habitats.sites.g1_g2_g3_b1_b2.60m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2.60m.unique)

relation.habitats.sites.g1_g2_g3_b1_b2.100m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2.100m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2.100m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2.100m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2.100m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_100m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2.100m.unique$surface~as.character(relation.habitats.sites.g1_g2_g3_b1_b2.100m.unique$Code_Site)+relation.habitats.sites.g1_g2_g3_b1_b2.100m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2.100m.unique)

relation.habitats.sites.g1_g2_g3_b1_b2.250m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2.250m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2.250m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2.250m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2.250m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_250m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2.250m.unique$surface~as.character(relation.habitats.sites.g1_g2_g3_b1_b2.250m.unique$Code_Site)+relation.habitats.sites.g1_g2_g3_b1_b2.250m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2.250m.unique)

relation.habitats.sites.g1_g2_g3_b1_b2.500m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2.500m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2.500m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2.500m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2.500m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_500m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2.500m.unique$surface~as.character(relation.habitats.sites.g1_g2_g3_b1_b2.500m.unique$Code_Site)+relation.habitats.sites.g1_g2_g3_b1_b2.500m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2.500m.unique)

relation.habitats.sites.g1_g2_g3_b1_b2.750m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2.750m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2.750m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2.750m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2.750m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_750m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2.750m.unique$surface~as.character(relation.habitats.sites.g1_g2_g3_b1_b2.750m.unique$Code_Site)+relation.habitats.sites.g1_g2_g3_b1_b2.750m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2.750m.unique)

relation.habitats.sites.g1_g2_g3_b1_b2.1000m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2.1000m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2.1000m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2.1000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2.1000m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_1000m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2.1000m.unique$surface~as.character(relation.habitats.sites.g1_g2_g3_b1_b2.1000m.unique$Code_Site)+relation.habitats.sites.g1_g2_g3_b1_b2.1000m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2.1000m.unique)

relation.habitats.sites.g1_g2_g3_b1_b2.3000m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2.3000m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2.3000m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2.3000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2.3000m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_3000m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2.3000m.unique$surface~as.character(relation.habitats.sites.g1_g2_g3_b1_b2.3000m.unique$Code_Site)+relation.habitats.sites.g1_g2_g3_b1_b2.3000m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2.3000m.unique)

relation.habitats.sites.g1_g2_g3_b1_b2.5000m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2.5000m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2.5000m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2.5000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2.5000m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_5000m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2.5000m.unique$surface~as.character(relation.habitats.sites.g1_g2_g3_b1_b2.5000m.unique$Code_Site)+relation.habitats.sites.g1_g2_g3_b1_b2.5000m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2.5000m.unique)

#Load de l'algorithme shannon
source("script.shannon.r")

#tirage des sites poissons en optimisant la diversité shannon d'habitat
 selection.shannon.hab.g1_g2_g3_b1_b2.30m<-shannon(tablo_g1_g2_g3_b1_b2_30m)
 selection.shannon.hab.g1_g2_g3_b1_b2.60m<-shannon(tablo_g1_g2_g3_b1_b2_60m)
 selection.shannon.hab.g1_g2_g3_b1_b2.100m<-shannon(tablo_g1_g2_g3_b1_b2_100m)
 selection.shannon.hab.g1_g2_g3_b1_b2.250m<-shannon(tablo_g1_g2_g3_b1_b2_250m)
 selection.shannon.hab.g1_g2_g3_b1_b2.500m<-shannon(tablo_g1_g2_g3_b1_b2_500m)
 selection.shannon.hab.g1_g2_g3_b1_b2.750m<-shannon(tablo_g1_g2_g3_b1_b2_750m)
 selection.shannon.hab.g1_g2_g3_b1_b2.1000m<-shannon(tablo_g1_g2_g3_b1_b2_1000m)
 selection.shannon.hab.g1_g2_g3_b1_b2.3000m<-shannon(tablo_g1_g2_g3_b1_b2_3000m)
 selection.shannon.hab.g1_g2_g3_b1_b2.5000m<-shannon(tablo_g1_g2_g3_b1_b2_5000m)

  #DIVERSITE SELECTIONNE PAR LE SCRIPT SHANNON
 source('script.diversite.shannon.selectionnee.fonction.r')
diversite.selection.shannon.hab.g1_g2_g3_b1_b2.30m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2.30m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1_b2.60m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2.60m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1_b2.100m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2.100m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1_b2.250m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2.250m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1_b2.500m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2.500m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1_b2.750m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2.750m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1_b2.1000m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2.1000m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1_b2.3000m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2.3000m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1_b2.5000m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2.5000m,requette.abondance.poissons)

 #CALCUL DES SAI#
source('script.sai.r')
SAI_shannon.g1_g2_g3_b1_b2.30m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2.30m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1_b2.60m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2.60m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1_b2.100m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2.100m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1_b2.250m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2.250m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1_b2.500m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2.500m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1_b2.750m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2.750m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1_b2.1000m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2.1000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1_b2.3000m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2.3000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1_b2.5000m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2.5000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)



  #############Habitat g1_g2_g3_b1_b2_b3###################



   ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances

setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
relation.habitats.sites.g1_g2_g3_b1_b2_b3<-read.csv2("fc3_intersect_g1_g2_g3_b1_b2_b3.csv")

sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')
sites.etudies<-sites.etudies[!sites.etudies%in%suppressed.station]

relation.habitats.sites.g1_g2_g3_b1_b2_b3<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$Code_Site%in%sites.etudies,]

#création d'une nouvelle table site-habitat par distance
relation.habitats.sites.g1_g2_g3_b1_b2_b3.30m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==30,]
relation.habitats.sites.g1_g2_g3_b1_b2_b3.60m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==60,]
relation.habitats.sites.g1_g2_g3_b1_b2_b3.100m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==100,]
relation.habitats.sites.g1_g2_g3_b1_b2_b3.250m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==250,]
relation.habitats.sites.g1_g2_g3_b1_b2_b3.500m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==500,]
relation.habitats.sites.g1_g2_g3_b1_b2_b3.750m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==750,]
relation.habitats.sites.g1_g2_g3_b1_b2_b3.1000m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==1000,]
relation.habitats.sites.g1_g2_g3_b1_b2_b3.3000m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==3000,]
relation.habitats.sites.g1_g2_g3_b1_b2_b3.5000m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==5000,]

#cration d'un tableau adapté à l'entrée du script shannon
relation.habitats.sites.g1_g2_g3_b1_b2_b3.30m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2_b3.30m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2_b3.30m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2_b3.30m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2_b3.30m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_b3_30m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2_b3.30m.unique$surface~as.character(relation.habitats.sites.g1_g2_g3_b1_b2_b3.30m.unique$Code_Site)+relation.habitats.sites.g1_g2_g3_b1_b2_b3.30m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2_b3.30m.unique)

relation.habitats.sites.g1_g2_g3_b1_b2_b3.60m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2_b3.60m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2_b3.60m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2_b3.60m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2_b3.60m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_b3_60m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2_b3.60m.unique$surface~as.character(relation.habitats.sites.g1_g2_g3_b1_b2_b3.60m.unique$Code_Site)+relation.habitats.sites.g1_g2_g3_b1_b2_b3.60m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2_b3.60m.unique)

relation.habitats.sites.g1_g2_g3_b1_b2_b3.100m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2_b3.100m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2_b3.100m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2_b3.100m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2_b3.100m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_b3_100m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2_b3.100m.unique$surface~as.character(relation.habitats.sites.g1_g2_g3_b1_b2_b3.100m.unique$Code_Site)+relation.habitats.sites.g1_g2_g3_b1_b2_b3.100m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2_b3.100m.unique)

relation.habitats.sites.g1_g2_g3_b1_b2_b3.250m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2_b3.250m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2_b3.250m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2_b3.250m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2_b3.250m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_b3_250m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2_b3.250m.unique$surface~as.character(relation.habitats.sites.g1_g2_g3_b1_b2_b3.250m.unique$Code_Site)+relation.habitats.sites.g1_g2_g3_b1_b2_b3.250m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2_b3.250m.unique)

relation.habitats.sites.g1_g2_g3_b1_b2_b3.500m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2_b3.500m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2_b3.500m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2_b3.500m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2_b3.500m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_b3_500m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2_b3.500m.unique$surface~as.character(relation.habitats.sites.g1_g2_g3_b1_b2_b3.500m.unique$Code_Site)+relation.habitats.sites.g1_g2_g3_b1_b2_b3.500m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2_b3.500m.unique)

relation.habitats.sites.g1_g2_g3_b1_b2_b3.750m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2_b3.750m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2_b3.750m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2_b3.750m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2_b3.750m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_b3_750m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2_b3.750m.unique$surface~as.character(relation.habitats.sites.g1_g2_g3_b1_b2_b3.750m.unique$Code_Site)+relation.habitats.sites.g1_g2_g3_b1_b2_b3.750m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2_b3.750m.unique)

relation.habitats.sites.g1_g2_g3_b1_b2_b3.1000m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2_b3.1000m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2_b3.1000m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2_b3.1000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2_b3.1000m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_b3_1000m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2_b3.1000m.unique$surface~as.character(relation.habitats.sites.g1_g2_g3_b1_b2_b3.1000m.unique$Code_Site)+relation.habitats.sites.g1_g2_g3_b1_b2_b3.1000m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2_b3.1000m.unique)

relation.habitats.sites.g1_g2_g3_b1_b2_b3.3000m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2_b3.3000m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2_b3.3000m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2_b3.3000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2_b3.3000m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_b3_3000m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2_b3.3000m.unique$surface~as.character(relation.habitats.sites.g1_g2_g3_b1_b2_b3.3000m.unique$Code_Site)+relation.habitats.sites.g1_g2_g3_b1_b2_b3.3000m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2_b3.3000m.unique)

relation.habitats.sites.g1_g2_g3_b1_b2_b3.5000m.unique<-aggregate(relation.habitats.sites.g1_g2_g3_b1_b2_b3.5000m$surface,list(Code_Site=relation.habitats.sites.g1_g2_g3_b1_b2_b3.5000m$Code_Site,habitat=relation.habitats.sites.g1_g2_g3_b1_b2_b3.5000m$Code_hab),sum,na.rm=TRUE)
colnames(relation.habitats.sites.g1_g2_g3_b1_b2_b3.5000m.unique)=c("Code_Site","habitat","surface")
tablo_g1_g2_g3_b1_b2_b3_5000m<-xtabs(formula=relation.habitats.sites.g1_g2_g3_b1_b2_b3.5000m.unique$surface~as.character(relation.habitats.sites.g1_g2_g3_b1_b2_b3.5000m.unique$Code_Site)+relation.habitats.sites.g1_g2_g3_b1_b2_b3.5000m.unique$habitat,relation.habitats.sites.g1_g2_g3_b1_b2_b3.5000m.unique)

#Load de l'algorithme shannon
source("script.shannon.r")

#tirage des sites poissons en optimisant la diversité shannon d'habitat
 selection.shannon.hab.g1_g2_g3_b1_b2_b3.30m<-shannon(tablo_g1_g2_g3_b1_b2_b3_30m)
 selection.shannon.hab.g1_g2_g3_b1_b2_b3.60m<-shannon(tablo_g1_g2_g3_b1_b2_b3_60m)
 selection.shannon.hab.g1_g2_g3_b1_b2_b3.100m<-shannon(tablo_g1_g2_g3_b1_b2_b3_100m)
 selection.shannon.hab.g1_g2_g3_b1_b2_b3.250m<-shannon(tablo_g1_g2_g3_b1_b2_b3_250m)
 selection.shannon.hab.g1_g2_g3_b1_b2_b3.500m<-shannon(tablo_g1_g2_g3_b1_b2_b3_500m)
 selection.shannon.hab.g1_g2_g3_b1_b2_b3.750m<-shannon(tablo_g1_g2_g3_b1_b2_b3_750m)
 selection.shannon.hab.g1_g2_g3_b1_b2_b3.1000m<-shannon(tablo_g1_g2_g3_b1_b2_b3_1000m)
 selection.shannon.hab.g1_g2_g3_b1_b2_b3.3000m<-shannon(tablo_g1_g2_g3_b1_b2_b3_3000m)
 selection.shannon.hab.g1_g2_g3_b1_b2_b3.5000m<-shannon(tablo_g1_g2_g3_b1_b2_b3_5000m)

  #DIVERSITE SELECTIONNE PAR LE SCRIPT SHANNON
  source('script.diversite.shannon.selectionnee.fonction.r')
diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.30m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2_b3.30m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.60m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2_b3.60m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.100m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2_b3.100m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.250m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2_b3.250m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.500m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2_b3.500m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.750m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2_b3.750m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.1000m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2_b3.1000m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.3000m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2_b3.3000m,requette.abondance.poissons)
diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.5000m<-species.shannon.evol(selection.shannon.hab.g1_g2_g3_b1_b2_b3.5000m,requette.abondance.poissons)

 #CALCUL DES SAI#
source('script.sai.r')
SAI_shannon.g1_g2_g3_b1_b2_b3.30m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.30m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1_b2_b3.60m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.60m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1_b2_b3.100m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.100m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1_b2_b3.250m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.250m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1_b2_b3.500m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.500m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1_b2_b3.750m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.750m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1_b2_b3.1000m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.1000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1_b2_b3.3000m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.3000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)
SAI_shannon.g1_g2_g3_b1_b2_b3.5000m<-SAI(diversite.selection.shannon.hab.g1_g2_g3_b1_b2_b3.5000m,diversite.shannon.moyenne.selection.aleatoire,diversite.selection.shannon.poisson)


#STOCKAGE DES RESULTATS

targeted.station.results[,p+1]=c(SAI_shannon.geo1.30m,SAI_shannon.geo1.60m,SAI_shannon.geo1.100m,SAI_shannon.geo1.250m,
SAI_shannon.geo1.500m,SAI_shannon.geo1.750m,SAI_shannon.geo1.1000m,SAI_shannon.geo1.3000m,SAI_shannon.geo1.5000m,
SAI_shannon.geo2.30m,SAI_shannon.geo2.60m,SAI_shannon.geo2.100m,SAI_shannon.geo2.250m,
SAI_shannon.geo2.500m,SAI_shannon.geo2.750m,SAI_shannon.geo2.1000m,SAI_shannon.geo2.3000m,SAI_shannon.geo2.5000m,
SAI_shannon.geo3.100m,SAI_shannon.geo3.250m,
SAI_shannon.geo3.500m,SAI_shannon.geo3.750m,SAI_shannon.geo3.1000m,SAI_shannon.geo3.3000m,SAI_shannon.geo3.5000m,
SAI_shannon.bent1.30m,SAI_shannon.bent1.60m,SAI_shannon.bent1.100m,SAI_shannon.bent1.250m,
SAI_shannon.bent1.500m,SAI_shannon.bent1.750m,SAI_shannon.bent1.1000m,SAI_shannon.bent1.3000m,SAI_shannon.bent1.5000m,
SAI_shannon.bent2.30m,SAI_shannon.bent2.60m,SAI_shannon.bent2.100m,SAI_shannon.bent2.250m,
SAI_shannon.bent2.500m,SAI_shannon.bent2.750m,SAI_shannon.bent2.1000m,SAI_shannon.bent2.3000m,SAI_shannon.bent2.5000m,
SAI_shannon.bent3.30m,SAI_shannon.bent3.60m,SAI_shannon.bent3.100m,SAI_shannon.bent3.250m,
SAI_shannon.bent3.500m,SAI_shannon.bent3.750m,SAI_shannon.bent3.1000m,SAI_shannon.bent3.3000m,SAI_shannon.bent3.5000m,
SAI_shannon.geo1_geo2.30m,SAI_shannon.geo1_geo2.60m,SAI_shannon.geo1_geo2.100m,SAI_shannon.geo1_geo2.250m,
SAI_shannon.geo1_geo2.500m,SAI_shannon.geo1_geo2.750m,SAI_shannon.geo1_geo2.1000m,SAI_shannon.geo1_geo2.3000m,SAI_shannon.geo1_geo2.5000m,
SAI_shannon.geo1_geo2_geo3.30m,SAI_shannon.geo1_geo2_geo3.60m,SAI_shannon.geo1_geo2_geo3.100m,SAI_shannon.geo1_geo2_geo3.250m,
SAI_shannon.geo1_geo2_geo3.500m,SAI_shannon.geo1_geo2_geo3.750m,SAI_shannon.geo1_geo2_geo3.1000m,SAI_shannon.geo1_geo2_geo3.3000m,SAI_shannon.geo1_geo2_geo3.5000m,
SAI_shannon.g1_g2_g3_b1.30m,SAI_shannon.g1_g2_g3_b1.60m,SAI_shannon.g1_g2_g3_b1.100m,SAI_shannon.g1_g2_g3_b1.250m,
SAI_shannon.g1_g2_g3_b1.500m,SAI_shannon.g1_g2_g3_b1.750m,SAI_shannon.g1_g2_g3_b1.1000m,SAI_shannon.g1_g2_g3_b1.3000m,SAI_shannon.g1_g2_g3_b1.5000m,
SAI_shannon.g1_g2_g3_b1_b2.30m,SAI_shannon.g1_g2_g3_b1_b2.60m,SAI_shannon.g1_g2_g3_b1_b2.100m,SAI_shannon.g1_g2_g3_b1_b2.250m,
SAI_shannon.g1_g2_g3_b1_b2.500m,SAI_shannon.g1_g2_g3_b1_b2.750m,SAI_shannon.g1_g2_g3_b1_b2.1000m,SAI_shannon.g1_g2_g3_b1_b2.3000m,SAI_shannon.g1_g2_g3_b1_b2.5000m,
SAI_shannon.g1_g2_g3_b1_b2_b3.30m,SAI_shannon.g1_g2_g3_b1_b2_b3.60m,SAI_shannon.g1_g2_g3_b1_b2_b3.100m,SAI_shannon.g1_g2_g3_b1_b2_b3.250m,
SAI_shannon.g1_g2_g3_b1_b2_b3.500m,SAI_shannon.g1_g2_g3_b1_b2_b3.750m,SAI_shannon.g1_g2_g3_b1_b2_b3.1000m,SAI_shannon.g1_g2_g3_b1_b2_b3.3000m,SAI_shannon.g1_g2_g3_b1_b2_b3.5000m)
}

 targeted.station.results
