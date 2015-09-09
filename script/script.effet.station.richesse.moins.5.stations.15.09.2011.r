#Script effet station richesse
targeted.station.results<-matrix(0,97,100)  #999

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

repertoire<-"C:/Documents and Settings/vanw/Bureau/R"
setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
recensement=read.csv2("REQ_site_poisson.IRD.csv")

#J'enlève les doublons
recensement_unique=unique(recensement)

#création d'un tableau de présence absence
tablo_pres_abs<-table(recensement_unique$Code_Site,recensement_unique$Code_poisson)

#calcul du nombre d'espèce par station
vec_poisson<-apply(tablo_pres_abs,1,sum)

#selection aléatoire d'un pool de site
setwd(repertoire)
setwd(input)
source('script.aleatoire.r')
selection.aleatoire<-matrix(0,dim(tablo_pres_abs)[1],999) #999
for (i in 1:999){ #999
selection.aleatoire[,i]<-run.alea(tablo_pres_abs)}


#Selection des dites selon un scenario de richesse complementarité basé sur les poissons
setwd(repertoire)
setwd(input)
source('script.richesse.complementarite.r')

selection.richesse.complementarite.poisson<-run.rich(tablo_pres_abs)

#Richesse de poisson incluse dans les sites selectionnés (aléatoirement et selon diversité de poisson)
source('script.diversite.selectionnee.r')
diversite.selection.aleatoire<-matrix(0,dim(tablo_pres_abs)[1],999)    #999
 diversite.selection.aleatoire<-apply(as.matrix(selection.aleatoire),2,function(x) species.richness.evol(x,recensement_unique))
 diversite.moyenne.selection.aleatoire<-apply(diversite.selection.aleatoire,1,mean)
 

diversite.selection.richesse.complementarite.poisson<-species.richness.evol(selection.richesse.complementarite.poisson,recensement_unique)

#Calcul des courbes aléatore moyenne, uper, et low
diversite.moyenne.selection.aleatoire
diversite.borne.sup.selection.aleatoire=apply(diversite.selection.aleatoire,1,function(x) sort(x)[ceiling(0.975*length(x))])
diversite.borne.inf.selection.aleatoire=apply(diversite.selection.aleatoire,1,function(x) sort(x)[ceiling(0.025*length(x))])

#nombre total d'espèces de poissons observés
nombre.poissons.total <- nlevels(factor(recensement_unique$Code_poisson))

# fonction pourcentage
percentage=function(total,truc){
  truc * 100 / total
  }
  
# Pourcentages, RANDOM, upper, low
diversite.moyenne.selection.aleatoire.pourcentage <- c(0,percentage(nombre.poissons.total,diversite.moyenne.selection.aleatoire))
diversite.borne.sup.selection.aleatoire.pourcentage <- c(0,percentage(nombre.poissons.total,diversite.borne.sup.selection.aleatoire))
diversite.borne.inf.selection.aleatoire.pourcentage <- c(0,percentage(nombre.poissons.total,diversite.borne.inf.selection.aleatoire))

diversite.selection.richesse.complementarite.poisson.pourcentage<- c(0,percentage(nombre.poissons.total,diversite.selection.richesse.complementarite.poisson))




############################HABITAT geo1###########################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances


setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
relation.habitats.sites.geo1<-read.csv2("fc3_intersect_geo1.csv")

sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')

relation.habitats.sites.geo1<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$Code_Site%in%sites.etudies,]

#création d'une nouvelle table site-habitat par distance
#habitat géo1
relation.habitats.sites.geo1.30m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==30,]
relation.habitats.sites.geo1.60m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==60,]
relation.habitats.sites.geo1.100m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==100,]
relation.habitats.sites.geo1.250m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==250,]
relation.habitats.sites.geo1.500m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==500,]
relation.habitats.sites.geo1.750m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==750,]
relation.habitats.sites.geo1.1000m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==1000,]
relation.habitats.sites.geo1.3000m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==3000,]
relation.habitats.sites.geo1.5000m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==5000,]

#preparation des donnée d'entrée

 m=matrix(0,length(relation.habitats.sites.geo1.30m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1.30m$Code_Site)
  m[,2]=relation.habitats.sites.geo1.30m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_30m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1.60m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1.60m$Code_Site)
  m[,2]=relation.habitats.sites.geo1.60m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_60m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1.100m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1.100m$Code_Site)
  m[,2]=relation.habitats.sites.geo1.100m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_100m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1.250m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1.250m$Code_Site)
  m[,2]=relation.habitats.sites.geo1.250m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_250m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1.500m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1.500m$Code_Site)
  m[,2]=relation.habitats.sites.geo1.500m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_500m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1.750m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1.750m$Code_Site)
  m[,2]=relation.habitats.sites.geo1.750m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_750m=as.data.frame(unique(m))

 m=matrix(0,length(relation.habitats.sites.geo1.1000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1.1000m$Code_Site)
  m[,2]=relation.habitats.sites.geo1.1000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_1000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1.3000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1.3000m$Code_Site)
  m[,2]=relation.habitats.sites.geo1.3000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_3000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1.5000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1.5000m$Code_Site)
  m[,2]=relation.habitats.sites.geo1.5000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_5000m=as.data.frame(unique(m))


habitat_par_site_geo1_30m<-habitat_par_site_geo1_30m[habitat_par_site_geo1_30m$Code_Site%in%sites.etudies,]
habitat_par_site_geo1_60m<-habitat_par_site_geo1_60m[habitat_par_site_geo1_60m$Code_Site%in%sites.etudies,]
habitat_par_site_geo1_100m<-habitat_par_site_geo1_100m[habitat_par_site_geo1_100m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_250m<-habitat_par_site_geo1_250m[habitat_par_site_geo1_250m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_500m<-habitat_par_site_geo1_500m[habitat_par_site_geo1_500m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_750m<-habitat_par_site_geo1_750m[habitat_par_site_geo1_750m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_1000m<-habitat_par_site_geo1_1000m[habitat_par_site_geo1_1000m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_3000m<-habitat_par_site_geo1_3000m[habitat_par_site_geo1_3000m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_5000m<-habitat_par_site_geo1_5000m[habitat_par_site_geo1_5000m$Code_Site%in%sites.etudies,]


#création de table de présence absence
#habitat géo1

habitat_par_site_geo1_30m_pres_abs<-table(habitat_par_site_geo1_30m$Code_Site,habitat_par_site_geo1_30m$Code_hab)
habitat_par_site_geo1_60m_pres_abs<-table(habitat_par_site_geo1_60m$Code_Site,habitat_par_site_geo1_60m$Code_hab)
habitat_par_site_geo1_100m_pres_abs<-table(habitat_par_site_geo1_100m$Code_Site,habitat_par_site_geo1_100m$Code_hab)
habitat_par_site_geo1_250m_pres_abs<-table(habitat_par_site_geo1_250m$Code_Site,habitat_par_site_geo1_250m$Code_hab)
habitat_par_site_geo1_500m_pres_abs<-table(habitat_par_site_geo1_500m$Code_Site,habitat_par_site_geo1_500m$Code_hab)
habitat_par_site_geo1_750m_pres_abs<-table(habitat_par_site_geo1_750m$Code_Site,habitat_par_site_geo1_750m$Code_hab)
habitat_par_site_geo1_1000m_pres_abs<-table(habitat_par_site_geo1_1000m$Code_Site,habitat_par_site_geo1_1000m$Code_hab)
habitat_par_site_geo1_3000m_pres_abs<-table(habitat_par_site_geo1_3000m$Code_Site,habitat_par_site_geo1_3000m$Code_hab)
habitat_par_site_geo1_5000m_pres_abs<-table(habitat_par_site_geo1_5000m$Code_Site,habitat_par_site_geo1_5000m$Code_hab)


#selection des sites
setwd(repertoire)
setwd(input)
source('script.richesse.complementarite.r')

selection.richesse.complementarite.geo1.30m<-run.rich(habitat_par_site_geo1_30m_pres_abs)
selection.richesse.complementarite.geo1.60m<-run.rich(habitat_par_site_geo1_60m_pres_abs)
selection.richesse.complementarite.geo1.100m<-run.rich(habitat_par_site_geo1_100m_pres_abs)
selection.richesse.complementarite.geo1.250m<-run.rich(habitat_par_site_geo1_250m_pres_abs)
selection.richesse.complementarite.geo1.500m<-run.rich(habitat_par_site_geo1_500m_pres_abs)
selection.richesse.complementarite.geo1.750m<-run.rich(habitat_par_site_geo1_750m_pres_abs)
selection.richesse.complementarite.geo1.1000m<-run.rich(habitat_par_site_geo1_1000m_pres_abs)
selection.richesse.complementarite.geo1.3000m<-run.rich(habitat_par_site_geo1_3000m_pres_abs)
selection.richesse.complementarite.geo1.5000m<-run.rich(habitat_par_site_geo1_5000m_pres_abs)

#RICHESSE ET DIVERSITE DE POISSON CONTENUE DANS LES SELECTIONS DE SITES#
source('script.diversite.selectionnee.r')

diversite.selection.richesse.complementarite.geo1.30m<-species.richness.evol(selection.richesse.complementarite.geo1.30m,recensement_unique)
diversite.selection.richesse.complementarite.geo1.60m<-species.richness.evol(selection.richesse.complementarite.geo1.60m,recensement_unique)
diversite.selection.richesse.complementarite.geo1.100m<-species.richness.evol(selection.richesse.complementarite.geo1.100m,recensement_unique)
diversite.selection.richesse.complementarite.geo1.250m<-species.richness.evol(selection.richesse.complementarite.geo1.250m,recensement_unique)
diversite.selection.richesse.complementarite.geo1.500m<-species.richness.evol(selection.richesse.complementarite.geo1.500m,recensement_unique)
diversite.selection.richesse.complementarite.geo1.750m<-species.richness.evol(selection.richesse.complementarite.geo1.750m,recensement_unique)
diversite.selection.richesse.complementarite.geo1.1000m<-species.richness.evol(selection.richesse.complementarite.geo1.1000m,recensement_unique)
diversite.selection.richesse.complementarite.geo1.3000m<-species.richness.evol(selection.richesse.complementarite.geo1.3000m,recensement_unique)
diversite.selection.richesse.complementarite.geo1.5000m<-species.richness.evol(selection.richesse.complementarite.geo1.5000m,recensement_unique)

#Calcul des SAI
source('script.sai.r')
SAI_richesse.complementarite.geo1.30m<-SAI(diversite.selection.richesse.complementarite.geo1.30m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1.60m<-SAI(diversite.selection.richesse.complementarite.geo1.60m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1.100m<-SAI(diversite.selection.richesse.complementarite.geo1.100m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1.250m<-SAI(diversite.selection.richesse.complementarite.geo1.250m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1.500m<-SAI(diversite.selection.richesse.complementarite.geo1.500m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1.750m<-SAI(diversite.selection.richesse.complementarite.geo1.750m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1.1000m<-SAI(diversite.selection.richesse.complementarite.geo1.1000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1.3000m<-SAI(diversite.selection.richesse.complementarite.geo1.3000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1.5000m<-SAI(diversite.selection.richesse.complementarite.geo1.5000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)



############################HABITAT geo2###########################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances


setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
relation.habitats.sites.geo2<-read.csv2("fc3_intersect_geo2.csv")

sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')

relation.habitats.sites.geo2<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$Code_Site%in%sites.etudies,]

#création d'une nouvelle table site-habitat par distance
#habitat géo1
relation.habitats.sites.geo2.30m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==30,]
relation.habitats.sites.geo2.60m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==60,]
relation.habitats.sites.geo2.100m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==100,]
relation.habitats.sites.geo2.250m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==250,]
relation.habitats.sites.geo2.500m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==500,]
relation.habitats.sites.geo2.750m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==750,]
relation.habitats.sites.geo2.1000m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==1000,]
relation.habitats.sites.geo2.3000m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==3000,]
relation.habitats.sites.geo2.5000m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==5000,]

#preparation des donnée d'entrée

 m=matrix(0,length(relation.habitats.sites.geo2.30m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo2.30m$Code_Site)
  m[,2]=relation.habitats.sites.geo2.30m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo2_30m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo2.60m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo2.60m$Code_Site)
  m[,2]=relation.habitats.sites.geo2.60m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo2_60m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo2.100m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo2.100m$Code_Site)
  m[,2]=relation.habitats.sites.geo2.100m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo2_100m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo2.250m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo2.250m$Code_Site)
  m[,2]=relation.habitats.sites.geo2.250m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo2_250m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo2.500m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo2.500m$Code_Site)
  m[,2]=relation.habitats.sites.geo2.500m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo2_500m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo2.750m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo2.750m$Code_Site)
  m[,2]=relation.habitats.sites.geo2.750m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo2_750m=as.data.frame(unique(m))

 m=matrix(0,length(relation.habitats.sites.geo2.1000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo2.1000m$Code_Site)
  m[,2]=relation.habitats.sites.geo2.1000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo2_1000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo2.3000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo2.3000m$Code_Site)
  m[,2]=relation.habitats.sites.geo2.3000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo2_3000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo2.5000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo2.5000m$Code_Site)
  m[,2]=relation.habitats.sites.geo2.5000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo2_5000m=as.data.frame(unique(m))


habitat_par_site_geo2_30m<-habitat_par_site_geo2_30m[habitat_par_site_geo2_30m$Code_Site%in%sites.etudies,]
habitat_par_site_geo2_60m<-habitat_par_site_geo2_60m[habitat_par_site_geo2_60m$Code_Site%in%sites.etudies,]
habitat_par_site_geo2_100m<-habitat_par_site_geo2_100m[habitat_par_site_geo2_100m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo2_250m<-habitat_par_site_geo2_250m[habitat_par_site_geo2_250m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo2_500m<-habitat_par_site_geo2_500m[habitat_par_site_geo2_500m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo2_750m<-habitat_par_site_geo2_750m[habitat_par_site_geo2_750m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo2_1000m<-habitat_par_site_geo2_1000m[habitat_par_site_geo2_1000m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo2_3000m<-habitat_par_site_geo2_3000m[habitat_par_site_geo2_3000m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo2_5000m<-habitat_par_site_geo2_5000m[habitat_par_site_geo2_5000m$Code_Site%in%sites.etudies,]


#création de table de présence absence
#habitat géo1

habitat_par_site_geo2_30m_pres_abs<-table(habitat_par_site_geo2_30m$Code_Site,habitat_par_site_geo2_30m$Code_hab)
habitat_par_site_geo2_60m_pres_abs<-table(habitat_par_site_geo2_60m$Code_Site,habitat_par_site_geo2_60m$Code_hab)
habitat_par_site_geo2_100m_pres_abs<-table(habitat_par_site_geo2_100m$Code_Site,habitat_par_site_geo2_100m$Code_hab)
habitat_par_site_geo2_250m_pres_abs<-table(habitat_par_site_geo2_250m$Code_Site,habitat_par_site_geo2_250m$Code_hab)
habitat_par_site_geo2_500m_pres_abs<-table(habitat_par_site_geo2_500m$Code_Site,habitat_par_site_geo2_500m$Code_hab)
habitat_par_site_geo2_750m_pres_abs<-table(habitat_par_site_geo2_750m$Code_Site,habitat_par_site_geo2_750m$Code_hab)
habitat_par_site_geo2_1000m_pres_abs<-table(habitat_par_site_geo2_1000m$Code_Site,habitat_par_site_geo2_1000m$Code_hab)
habitat_par_site_geo2_3000m_pres_abs<-table(habitat_par_site_geo2_3000m$Code_Site,habitat_par_site_geo2_3000m$Code_hab)
habitat_par_site_geo2_5000m_pres_abs<-table(habitat_par_site_geo2_5000m$Code_Site,habitat_par_site_geo2_5000m$Code_hab)


#selection des sites
setwd(repertoire)
setwd(input)
source('script.richesse.complementarite.r')

selection.richesse.complementarite.geo2.30m<-run.rich(habitat_par_site_geo2_30m_pres_abs)
selection.richesse.complementarite.geo2.60m<-run.rich(habitat_par_site_geo2_60m_pres_abs)
selection.richesse.complementarite.geo2.100m<-run.rich(habitat_par_site_geo2_100m_pres_abs)
selection.richesse.complementarite.geo2.250m<-run.rich(habitat_par_site_geo2_250m_pres_abs)
selection.richesse.complementarite.geo2.500m<-run.rich(habitat_par_site_geo2_500m_pres_abs)
selection.richesse.complementarite.geo2.750m<-run.rich(habitat_par_site_geo2_750m_pres_abs)
selection.richesse.complementarite.geo2.1000m<-run.rich(habitat_par_site_geo2_1000m_pres_abs)
selection.richesse.complementarite.geo2.3000m<-run.rich(habitat_par_site_geo2_3000m_pres_abs)
selection.richesse.complementarite.geo2.5000m<-run.rich(habitat_par_site_geo2_5000m_pres_abs)

#RICHESSE ET DIVERSITE DE POISSON CONTENUE DANS LES SELECTIONS DE SITES#
source('script.diversite.selectionnee.r')

diversite.selection.richesse.complementarite.geo2.30m<-species.richness.evol(selection.richesse.complementarite.geo2.30m,recensement_unique)
diversite.selection.richesse.complementarite.geo2.60m<-species.richness.evol(selection.richesse.complementarite.geo2.60m,recensement_unique)
diversite.selection.richesse.complementarite.geo2.100m<-species.richness.evol(selection.richesse.complementarite.geo2.100m,recensement_unique)
diversite.selection.richesse.complementarite.geo2.250m<-species.richness.evol(selection.richesse.complementarite.geo2.250m,recensement_unique)
diversite.selection.richesse.complementarite.geo2.500m<-species.richness.evol(selection.richesse.complementarite.geo2.500m,recensement_unique)
diversite.selection.richesse.complementarite.geo2.750m<-species.richness.evol(selection.richesse.complementarite.geo2.750m,recensement_unique)
diversite.selection.richesse.complementarite.geo2.1000m<-species.richness.evol(selection.richesse.complementarite.geo2.1000m,recensement_unique)
diversite.selection.richesse.complementarite.geo2.3000m<-species.richness.evol(selection.richesse.complementarite.geo2.3000m,recensement_unique)
diversite.selection.richesse.complementarite.geo2.5000m<-species.richness.evol(selection.richesse.complementarite.geo2.5000m,recensement_unique)

#Calcul des SAI
source('script.sai.r')
SAI_richesse.complementarite.geo2.30m<-SAI(diversite.selection.richesse.complementarite.geo2.30m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo2.60m<-SAI(diversite.selection.richesse.complementarite.geo2.60m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo2.100m<-SAI(diversite.selection.richesse.complementarite.geo2.100m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo2.250m<-SAI(diversite.selection.richesse.complementarite.geo2.250m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo2.500m<-SAI(diversite.selection.richesse.complementarite.geo2.500m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo2.750m<-SAI(diversite.selection.richesse.complementarite.geo2.750m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo2.1000m<-SAI(diversite.selection.richesse.complementarite.geo2.1000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo2.3000m<-SAI(diversite.selection.richesse.complementarite.geo2.3000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo2.5000m<-SAI(diversite.selection.richesse.complementarite.geo2.5000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)



############################HABITAT geo3###########################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances


setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
relation.habitats.sites.geo3<-read.csv2("fc3_intersect_geo3.csv")

sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')

relation.habitats.sites.geo3<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$Code_Site%in%sites.etudies,]

#création d'une nouvelle table site-habitat par distance
#habitat géo1

relation.habitats.sites.geo3.100m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==100,]
relation.habitats.sites.geo3.250m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==250,]
relation.habitats.sites.geo3.500m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==500,]
relation.habitats.sites.geo3.750m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==750,]
relation.habitats.sites.geo3.1000m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==1000,]
relation.habitats.sites.geo3.3000m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==3000,]
relation.habitats.sites.geo3.5000m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==5000,]

#preparation des donnée d'entrée



 m=matrix(0,length(relation.habitats.sites.geo3.100m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo3.100m$Code_Site)
  m[,2]=relation.habitats.sites.geo3.100m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo3_100m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo3.250m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo3.250m$Code_Site)
  m[,2]=relation.habitats.sites.geo3.250m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo3_250m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo3.500m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo3.500m$Code_Site)
  m[,2]=relation.habitats.sites.geo3.500m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo3_500m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo3.750m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo3.750m$Code_Site)
  m[,2]=relation.habitats.sites.geo3.750m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo3_750m=as.data.frame(unique(m))

 m=matrix(0,length(relation.habitats.sites.geo3.1000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo3.1000m$Code_Site)
  m[,2]=relation.habitats.sites.geo3.1000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo3_1000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo3.3000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo3.3000m$Code_Site)
  m[,2]=relation.habitats.sites.geo3.3000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo3_3000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo3.5000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo3.5000m$Code_Site)
  m[,2]=relation.habitats.sites.geo3.5000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo3_5000m=as.data.frame(unique(m))


habitat_par_site_geo3_100m<-habitat_par_site_geo3_100m[habitat_par_site_geo3_100m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo3_250m<-habitat_par_site_geo3_250m[habitat_par_site_geo3_250m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo3_500m<-habitat_par_site_geo3_500m[habitat_par_site_geo3_500m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo3_750m<-habitat_par_site_geo3_750m[habitat_par_site_geo3_750m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo3_1000m<-habitat_par_site_geo3_1000m[habitat_par_site_geo3_1000m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo3_3000m<-habitat_par_site_geo3_3000m[habitat_par_site_geo3_3000m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo3_5000m<-habitat_par_site_geo3_5000m[habitat_par_site_geo3_5000m$Code_Site%in%sites.etudies,]


#création de table de présence absence
#habitat géo1

habitat_par_site_geo3_100m_pres_abs<-table(habitat_par_site_geo3_100m$Code_Site,habitat_par_site_geo3_100m$Code_hab)
habitat_par_site_geo3_250m_pres_abs<-table(habitat_par_site_geo3_250m$Code_Site,habitat_par_site_geo3_250m$Code_hab)
habitat_par_site_geo3_500m_pres_abs<-table(habitat_par_site_geo3_500m$Code_Site,habitat_par_site_geo3_500m$Code_hab)
habitat_par_site_geo3_750m_pres_abs<-table(habitat_par_site_geo3_750m$Code_Site,habitat_par_site_geo3_750m$Code_hab)
habitat_par_site_geo3_1000m_pres_abs<-table(habitat_par_site_geo3_1000m$Code_Site,habitat_par_site_geo3_1000m$Code_hab)
habitat_par_site_geo3_3000m_pres_abs<-table(habitat_par_site_geo3_3000m$Code_Site,habitat_par_site_geo3_3000m$Code_hab)
habitat_par_site_geo3_5000m_pres_abs<-table(habitat_par_site_geo3_5000m$Code_Site,habitat_par_site_geo3_5000m$Code_hab)


#selection des sites
setwd(repertoire)
setwd(input)
source('script.richesse.complementarite.r')

selection.richesse.complementarite.geo3.100m<-run.rich(habitat_par_site_geo3_100m_pres_abs)
selection.richesse.complementarite.geo3.250m<-run.rich(habitat_par_site_geo3_250m_pres_abs)
selection.richesse.complementarite.geo3.500m<-run.rich(habitat_par_site_geo3_500m_pres_abs)
selection.richesse.complementarite.geo3.750m<-run.rich(habitat_par_site_geo3_750m_pres_abs)
selection.richesse.complementarite.geo3.1000m<-run.rich(habitat_par_site_geo3_1000m_pres_abs)
selection.richesse.complementarite.geo3.3000m<-run.rich(habitat_par_site_geo3_3000m_pres_abs)
selection.richesse.complementarite.geo3.5000m<-run.rich(habitat_par_site_geo3_5000m_pres_abs)

#RICHESSE ET DIVERSITE DE POISSON CONTENUE DANS LES SELECTIONS DE SITES#
source('script.diversite.selectionnee.r')

diversite.selection.richesse.complementarite.geo3.100m<-species.richness.evol(selection.richesse.complementarite.geo3.100m,recensement_unique)
diversite.selection.richesse.complementarite.geo3.250m<-species.richness.evol(selection.richesse.complementarite.geo3.250m,recensement_unique)
diversite.selection.richesse.complementarite.geo3.500m<-species.richness.evol(selection.richesse.complementarite.geo3.500m,recensement_unique)
diversite.selection.richesse.complementarite.geo3.750m<-species.richness.evol(selection.richesse.complementarite.geo3.750m,recensement_unique)
diversite.selection.richesse.complementarite.geo3.1000m<-species.richness.evol(selection.richesse.complementarite.geo3.1000m,recensement_unique)
diversite.selection.richesse.complementarite.geo3.3000m<-species.richness.evol(selection.richesse.complementarite.geo3.3000m,recensement_unique)
diversite.selection.richesse.complementarite.geo3.5000m<-species.richness.evol(selection.richesse.complementarite.geo3.5000m,recensement_unique)

#Calcul des SAI
source('script.sai.r')
SAI_richesse.complementarite.geo3.100m<-SAI(diversite.selection.richesse.complementarite.geo3.100m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo3.250m<-SAI(diversite.selection.richesse.complementarite.geo3.250m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo3.500m<-SAI(diversite.selection.richesse.complementarite.geo3.500m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo3.750m<-SAI(diversite.selection.richesse.complementarite.geo3.750m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo3.1000m<-SAI(diversite.selection.richesse.complementarite.geo3.1000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo3.3000m<-SAI(diversite.selection.richesse.complementarite.geo3.3000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo3.5000m<-SAI(diversite.selection.richesse.complementarite.geo3.5000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)



############################HABITAT bent1###########################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances


setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
relation.habitats.sites.bent1<-read.csv2("fc3_intersect_bent1.csv")

sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')

relation.habitats.sites.bent1<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$Code_Site%in%sites.etudies,]

#création d'une nouvelle table site-habitat par distance
#habitat géo1
relation.habitats.sites.bent1.30m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==30,]
relation.habitats.sites.bent1.60m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==60,]
relation.habitats.sites.bent1.100m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==100,]
relation.habitats.sites.bent1.250m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==250,]
relation.habitats.sites.bent1.500m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==500,]
relation.habitats.sites.bent1.750m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==750,]
relation.habitats.sites.bent1.1000m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==1000,]
relation.habitats.sites.bent1.3000m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==3000,]
relation.habitats.sites.bent1.5000m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==5000,]

#preparation des donnée d'entrée

 m=matrix(0,length(relation.habitats.sites.bent1.30m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent1.30m$Code_Site)
  m[,2]=relation.habitats.sites.bent1.30m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent1_30m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent1.60m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent1.60m$Code_Site)
  m[,2]=relation.habitats.sites.bent1.60m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent1_60m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent1.100m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent1.100m$Code_Site)
  m[,2]=relation.habitats.sites.bent1.100m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent1_100m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent1.250m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent1.250m$Code_Site)
  m[,2]=relation.habitats.sites.bent1.250m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent1_250m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent1.500m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent1.500m$Code_Site)
  m[,2]=relation.habitats.sites.bent1.500m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent1_500m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent1.750m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent1.750m$Code_Site)
  m[,2]=relation.habitats.sites.bent1.750m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent1_750m=as.data.frame(unique(m))

 m=matrix(0,length(relation.habitats.sites.bent1.1000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent1.1000m$Code_Site)
  m[,2]=relation.habitats.sites.bent1.1000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent1_1000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent1.3000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent1.3000m$Code_Site)
  m[,2]=relation.habitats.sites.bent1.3000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent1_3000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent1.5000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent1.5000m$Code_Site)
  m[,2]=relation.habitats.sites.bent1.5000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent1_5000m=as.data.frame(unique(m))


habitat_par_site_bent1_30m<-habitat_par_site_bent1_30m[habitat_par_site_bent1_30m$Code_Site%in%sites.etudies,]
habitat_par_site_bent1_60m<-habitat_par_site_bent1_60m[habitat_par_site_bent1_60m$Code_Site%in%sites.etudies,]
habitat_par_site_bent1_100m<-habitat_par_site_bent1_100m[habitat_par_site_bent1_100m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent1_250m<-habitat_par_site_bent1_250m[habitat_par_site_bent1_250m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent1_500m<-habitat_par_site_bent1_500m[habitat_par_site_bent1_500m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent1_750m<-habitat_par_site_bent1_750m[habitat_par_site_bent1_750m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent1_1000m<-habitat_par_site_bent1_1000m[habitat_par_site_bent1_1000m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent1_3000m<-habitat_par_site_bent1_3000m[habitat_par_site_bent1_3000m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent1_5000m<-habitat_par_site_bent1_5000m[habitat_par_site_bent1_5000m$Code_Site%in%sites.etudies,]


#création de table de présence absence
#habitat géo1

habitat_par_site_bent1_30m_pres_abs<-table(habitat_par_site_bent1_30m$Code_Site,habitat_par_site_bent1_30m$Code_hab)
habitat_par_site_bent1_60m_pres_abs<-table(habitat_par_site_bent1_60m$Code_Site,habitat_par_site_bent1_60m$Code_hab)
habitat_par_site_bent1_100m_pres_abs<-table(habitat_par_site_bent1_100m$Code_Site,habitat_par_site_bent1_100m$Code_hab)
habitat_par_site_bent1_250m_pres_abs<-table(habitat_par_site_bent1_250m$Code_Site,habitat_par_site_bent1_250m$Code_hab)
habitat_par_site_bent1_500m_pres_abs<-table(habitat_par_site_bent1_500m$Code_Site,habitat_par_site_bent1_500m$Code_hab)
habitat_par_site_bent1_750m_pres_abs<-table(habitat_par_site_bent1_750m$Code_Site,habitat_par_site_bent1_750m$Code_hab)
habitat_par_site_bent1_1000m_pres_abs<-table(habitat_par_site_bent1_1000m$Code_Site,habitat_par_site_bent1_1000m$Code_hab)
habitat_par_site_bent1_3000m_pres_abs<-table(habitat_par_site_bent1_3000m$Code_Site,habitat_par_site_bent1_3000m$Code_hab)
habitat_par_site_bent1_5000m_pres_abs<-table(habitat_par_site_bent1_5000m$Code_Site,habitat_par_site_bent1_5000m$Code_hab)


#selection des sites
setwd(repertoire)
setwd(input)
source('script.richesse.complementarite.r')

selection.richesse.complementarite.bent1.30m<-run.rich(habitat_par_site_bent1_30m_pres_abs)
selection.richesse.complementarite.bent1.60m<-run.rich(habitat_par_site_bent1_60m_pres_abs)
selection.richesse.complementarite.bent1.100m<-run.rich(habitat_par_site_bent1_100m_pres_abs)
selection.richesse.complementarite.bent1.250m<-run.rich(habitat_par_site_bent1_250m_pres_abs)
selection.richesse.complementarite.bent1.500m<-run.rich(habitat_par_site_bent1_500m_pres_abs)
selection.richesse.complementarite.bent1.750m<-run.rich(habitat_par_site_bent1_750m_pres_abs)
selection.richesse.complementarite.bent1.1000m<-run.rich(habitat_par_site_bent1_1000m_pres_abs)
selection.richesse.complementarite.bent1.3000m<-run.rich(habitat_par_site_bent1_3000m_pres_abs)
selection.richesse.complementarite.bent1.5000m<-run.rich(habitat_par_site_bent1_5000m_pres_abs)

#RICHESSE ET DIVERSITE DE POISSON CONTENUE DANS LES SELECTIONS DE SITES#
source('script.diversite.selectionnee.r')

diversite.selection.richesse.complementarite.bent1.30m<-species.richness.evol(selection.richesse.complementarite.bent1.30m,recensement_unique)
diversite.selection.richesse.complementarite.bent1.60m<-species.richness.evol(selection.richesse.complementarite.bent1.60m,recensement_unique)
diversite.selection.richesse.complementarite.bent1.100m<-species.richness.evol(selection.richesse.complementarite.bent1.100m,recensement_unique)
diversite.selection.richesse.complementarite.bent1.250m<-species.richness.evol(selection.richesse.complementarite.bent1.250m,recensement_unique)
diversite.selection.richesse.complementarite.bent1.500m<-species.richness.evol(selection.richesse.complementarite.bent1.500m,recensement_unique)
diversite.selection.richesse.complementarite.bent1.750m<-species.richness.evol(selection.richesse.complementarite.bent1.750m,recensement_unique)
diversite.selection.richesse.complementarite.bent1.1000m<-species.richness.evol(selection.richesse.complementarite.bent1.1000m,recensement_unique)
diversite.selection.richesse.complementarite.bent1.3000m<-species.richness.evol(selection.richesse.complementarite.bent1.3000m,recensement_unique)
diversite.selection.richesse.complementarite.bent1.5000m<-species.richness.evol(selection.richesse.complementarite.bent1.5000m,recensement_unique)

#Calcul des SAI
source('script.sai.r')
SAI_richesse.complementarite.bent1.30m<-SAI(diversite.selection.richesse.complementarite.bent1.30m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent1.60m<-SAI(diversite.selection.richesse.complementarite.bent1.60m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent1.100m<-SAI(diversite.selection.richesse.complementarite.bent1.100m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent1.250m<-SAI(diversite.selection.richesse.complementarite.bent1.250m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent1.500m<-SAI(diversite.selection.richesse.complementarite.bent1.500m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent1.750m<-SAI(diversite.selection.richesse.complementarite.bent1.750m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent1.1000m<-SAI(diversite.selection.richesse.complementarite.bent1.1000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent1.3000m<-SAI(diversite.selection.richesse.complementarite.bent1.3000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent1.5000m<-SAI(diversite.selection.richesse.complementarite.bent1.5000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)



############################HABITAT bent2###########################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances

setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
relation.habitats.sites.bent2<-read.csv2("fc3_intersect_bent2.csv")

sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')

relation.habitats.sites.bent2<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$Code_Site%in%sites.etudies,]

#création d'une nouvelle table site-habitat par distance
#habitat géo1
relation.habitats.sites.bent2.30m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==30,]
relation.habitats.sites.bent2.60m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==60,]
relation.habitats.sites.bent2.100m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==100,]
relation.habitats.sites.bent2.250m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==250,]
relation.habitats.sites.bent2.500m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==500,]
relation.habitats.sites.bent2.750m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==750,]
relation.habitats.sites.bent2.1000m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==1000,]
relation.habitats.sites.bent2.3000m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==3000,]
relation.habitats.sites.bent2.5000m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==5000,]

#preparation des donnée d'entrée

 m=matrix(0,length(relation.habitats.sites.bent2.30m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent2.30m$Code_Site)
  m[,2]=relation.habitats.sites.bent2.30m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent2_30m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent2.60m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent2.60m$Code_Site)
  m[,2]=relation.habitats.sites.bent2.60m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent2_60m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent2.100m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent2.100m$Code_Site)
  m[,2]=relation.habitats.sites.bent2.100m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent2_100m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent2.250m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent2.250m$Code_Site)
  m[,2]=relation.habitats.sites.bent2.250m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent2_250m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent2.500m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent2.500m$Code_Site)
  m[,2]=relation.habitats.sites.bent2.500m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent2_500m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent2.750m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent2.750m$Code_Site)
  m[,2]=relation.habitats.sites.bent2.750m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent2_750m=as.data.frame(unique(m))

 m=matrix(0,length(relation.habitats.sites.bent2.1000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent2.1000m$Code_Site)
  m[,2]=relation.habitats.sites.bent2.1000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent2_1000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent2.3000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent2.3000m$Code_Site)
  m[,2]=relation.habitats.sites.bent2.3000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent2_3000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent2.5000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent2.5000m$Code_Site)
  m[,2]=relation.habitats.sites.bent2.5000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent2_5000m=as.data.frame(unique(m))


habitat_par_site_bent2_30m<-habitat_par_site_bent2_30m[habitat_par_site_bent2_30m$Code_Site%in%sites.etudies,]
habitat_par_site_bent2_60m<-habitat_par_site_bent2_60m[habitat_par_site_bent2_60m$Code_Site%in%sites.etudies,]
habitat_par_site_bent2_100m<-habitat_par_site_bent2_100m[habitat_par_site_bent2_100m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent2_250m<-habitat_par_site_bent2_250m[habitat_par_site_bent2_250m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent2_500m<-habitat_par_site_bent2_500m[habitat_par_site_bent2_500m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent2_750m<-habitat_par_site_bent2_750m[habitat_par_site_bent2_750m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent2_1000m<-habitat_par_site_bent2_1000m[habitat_par_site_bent2_1000m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent2_3000m<-habitat_par_site_bent2_3000m[habitat_par_site_bent2_3000m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent2_5000m<-habitat_par_site_bent2_5000m[habitat_par_site_bent2_5000m$Code_Site%in%sites.etudies,]


#création de table de présence absence
#habitat géo1

habitat_par_site_bent2_30m_pres_abs<-table(habitat_par_site_bent2_30m$Code_Site,habitat_par_site_bent2_30m$Code_hab)
habitat_par_site_bent2_60m_pres_abs<-table(habitat_par_site_bent2_60m$Code_Site,habitat_par_site_bent2_60m$Code_hab)
habitat_par_site_bent2_100m_pres_abs<-table(habitat_par_site_bent2_100m$Code_Site,habitat_par_site_bent2_100m$Code_hab)
habitat_par_site_bent2_250m_pres_abs<-table(habitat_par_site_bent2_250m$Code_Site,habitat_par_site_bent2_250m$Code_hab)
habitat_par_site_bent2_500m_pres_abs<-table(habitat_par_site_bent2_500m$Code_Site,habitat_par_site_bent2_500m$Code_hab)
habitat_par_site_bent2_750m_pres_abs<-table(habitat_par_site_bent2_750m$Code_Site,habitat_par_site_bent2_750m$Code_hab)
habitat_par_site_bent2_1000m_pres_abs<-table(habitat_par_site_bent2_1000m$Code_Site,habitat_par_site_bent2_1000m$Code_hab)
habitat_par_site_bent2_3000m_pres_abs<-table(habitat_par_site_bent2_3000m$Code_Site,habitat_par_site_bent2_3000m$Code_hab)
habitat_par_site_bent2_5000m_pres_abs<-table(habitat_par_site_bent2_5000m$Code_Site,habitat_par_site_bent2_5000m$Code_hab)


#selection des sites
setwd(repertoire)
setwd(input)
source('script.richesse.complementarite.r')

selection.richesse.complementarite.bent2.30m<-run.rich(habitat_par_site_bent2_30m_pres_abs)
selection.richesse.complementarite.bent2.60m<-run.rich(habitat_par_site_bent2_60m_pres_abs)
selection.richesse.complementarite.bent2.100m<-run.rich(habitat_par_site_bent2_100m_pres_abs)
selection.richesse.complementarite.bent2.250m<-run.rich(habitat_par_site_bent2_250m_pres_abs)
selection.richesse.complementarite.bent2.500m<-run.rich(habitat_par_site_bent2_500m_pres_abs)
selection.richesse.complementarite.bent2.750m<-run.rich(habitat_par_site_bent2_750m_pres_abs)
selection.richesse.complementarite.bent2.1000m<-run.rich(habitat_par_site_bent2_1000m_pres_abs)
selection.richesse.complementarite.bent2.3000m<-run.rich(habitat_par_site_bent2_3000m_pres_abs)
selection.richesse.complementarite.bent2.5000m<-run.rich(habitat_par_site_bent2_5000m_pres_abs)

#RICHESSE ET DIVERSITE DE POISSON CONTENUE DANS LES SELECTIONS DE SITES#
source('script.diversite.selectionnee.r')

diversite.selection.richesse.complementarite.bent2.30m<-species.richness.evol(selection.richesse.complementarite.bent2.30m,recensement_unique)
diversite.selection.richesse.complementarite.bent2.60m<-species.richness.evol(selection.richesse.complementarite.bent2.60m,recensement_unique)
diversite.selection.richesse.complementarite.bent2.100m<-species.richness.evol(selection.richesse.complementarite.bent2.100m,recensement_unique)
diversite.selection.richesse.complementarite.bent2.250m<-species.richness.evol(selection.richesse.complementarite.bent2.250m,recensement_unique)
diversite.selection.richesse.complementarite.bent2.500m<-species.richness.evol(selection.richesse.complementarite.bent2.500m,recensement_unique)
diversite.selection.richesse.complementarite.bent2.750m<-species.richness.evol(selection.richesse.complementarite.bent2.750m,recensement_unique)
diversite.selection.richesse.complementarite.bent2.1000m<-species.richness.evol(selection.richesse.complementarite.bent2.1000m,recensement_unique)
diversite.selection.richesse.complementarite.bent2.3000m<-species.richness.evol(selection.richesse.complementarite.bent2.3000m,recensement_unique)
diversite.selection.richesse.complementarite.bent2.5000m<-species.richness.evol(selection.richesse.complementarite.bent2.5000m,recensement_unique)

#Calcul des SAI
source('script.sai.r')
SAI_richesse.complementarite.bent2.30m<-SAI(diversite.selection.richesse.complementarite.bent2.30m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent2.60m<-SAI(diversite.selection.richesse.complementarite.bent2.60m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent2.100m<-SAI(diversite.selection.richesse.complementarite.bent2.100m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent2.250m<-SAI(diversite.selection.richesse.complementarite.bent2.250m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent2.500m<-SAI(diversite.selection.richesse.complementarite.bent2.500m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent2.750m<-SAI(diversite.selection.richesse.complementarite.bent2.750m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent2.1000m<-SAI(diversite.selection.richesse.complementarite.bent2.1000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent2.3000m<-SAI(diversite.selection.richesse.complementarite.bent2.3000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent2.5000m<-SAI(diversite.selection.richesse.complementarite.bent2.5000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)



############################HABITAT bent3###########################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances

setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
relation.habitats.sites.bent3<-read.csv2("fc3_intersect_bent3.csv")

sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')

relation.habitats.sites.bent3<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$Code_Site%in%sites.etudies,]

#création d'une nouvelle table site-habitat par distance
#habitat géo1
relation.habitats.sites.bent3.30m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==30,]
relation.habitats.sites.bent3.60m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==60,]
relation.habitats.sites.bent3.100m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==100,]
relation.habitats.sites.bent3.250m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==250,]
relation.habitats.sites.bent3.500m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==500,]
relation.habitats.sites.bent3.750m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==750,]
relation.habitats.sites.bent3.1000m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==1000,]
relation.habitats.sites.bent3.3000m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==3000,]
relation.habitats.sites.bent3.5000m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==5000,]

#preparation des donnée d'entrée

 m=matrix(0,length(relation.habitats.sites.bent3.30m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent3.30m$Code_Site)
  m[,2]=relation.habitats.sites.bent3.30m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent3_30m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent3.60m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent3.60m$Code_Site)
  m[,2]=relation.habitats.sites.bent3.60m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent3_60m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent3.100m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent3.100m$Code_Site)
  m[,2]=relation.habitats.sites.bent3.100m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent3_100m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent3.250m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent3.250m$Code_Site)
  m[,2]=relation.habitats.sites.bent3.250m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent3_250m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent3.500m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent3.500m$Code_Site)
  m[,2]=relation.habitats.sites.bent3.500m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent3_500m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent3.750m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent3.750m$Code_Site)
  m[,2]=relation.habitats.sites.bent3.750m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent3_750m=as.data.frame(unique(m))

 m=matrix(0,length(relation.habitats.sites.bent3.1000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent3.1000m$Code_Site)
  m[,2]=relation.habitats.sites.bent3.1000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent3_1000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent3.3000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent3.3000m$Code_Site)
  m[,2]=relation.habitats.sites.bent3.3000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent3_3000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent3.5000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent3.5000m$Code_Site)
  m[,2]=relation.habitats.sites.bent3.5000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent3_5000m=as.data.frame(unique(m))


habitat_par_site_bent3_30m<-habitat_par_site_bent3_30m[habitat_par_site_bent3_30m$Code_Site%in%sites.etudies,]
habitat_par_site_bent3_60m<-habitat_par_site_bent3_60m[habitat_par_site_bent3_60m$Code_Site%in%sites.etudies,]
habitat_par_site_bent3_100m<-habitat_par_site_bent3_100m[habitat_par_site_bent3_100m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent3_250m<-habitat_par_site_bent3_250m[habitat_par_site_bent3_250m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent3_500m<-habitat_par_site_bent3_500m[habitat_par_site_bent3_500m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent3_750m<-habitat_par_site_bent3_750m[habitat_par_site_bent3_750m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent3_1000m<-habitat_par_site_bent3_1000m[habitat_par_site_bent3_1000m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent3_3000m<-habitat_par_site_bent3_3000m[habitat_par_site_bent3_3000m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent3_5000m<-habitat_par_site_bent3_5000m[habitat_par_site_bent3_5000m$Code_Site%in%sites.etudies,]


#création de table de présence absence
#habitat géo1

habitat_par_site_bent3_30m_pres_abs<-table(habitat_par_site_bent3_30m$Code_Site,habitat_par_site_bent3_30m$Code_hab)
habitat_par_site_bent3_60m_pres_abs<-table(habitat_par_site_bent3_60m$Code_Site,habitat_par_site_bent3_60m$Code_hab)
habitat_par_site_bent3_100m_pres_abs<-table(habitat_par_site_bent3_100m$Code_Site,habitat_par_site_bent3_100m$Code_hab)
habitat_par_site_bent3_250m_pres_abs<-table(habitat_par_site_bent3_250m$Code_Site,habitat_par_site_bent3_250m$Code_hab)
habitat_par_site_bent3_500m_pres_abs<-table(habitat_par_site_bent3_500m$Code_Site,habitat_par_site_bent3_500m$Code_hab)
habitat_par_site_bent3_750m_pres_abs<-table(habitat_par_site_bent3_750m$Code_Site,habitat_par_site_bent3_750m$Code_hab)
habitat_par_site_bent3_1000m_pres_abs<-table(habitat_par_site_bent3_1000m$Code_Site,habitat_par_site_bent3_1000m$Code_hab)
habitat_par_site_bent3_3000m_pres_abs<-table(habitat_par_site_bent3_3000m$Code_Site,habitat_par_site_bent3_3000m$Code_hab)
habitat_par_site_bent3_5000m_pres_abs<-table(habitat_par_site_bent3_5000m$Code_Site,habitat_par_site_bent3_5000m$Code_hab)


#selection des sites
setwd(repertoire)
setwd(input)
source('script.richesse.complementarite.r')

selection.richesse.complementarite.bent3.30m<-run.rich(habitat_par_site_bent3_30m_pres_abs)
selection.richesse.complementarite.bent3.60m<-run.rich(habitat_par_site_bent3_60m_pres_abs)
selection.richesse.complementarite.bent3.100m<-run.rich(habitat_par_site_bent3_100m_pres_abs)
selection.richesse.complementarite.bent3.250m<-run.rich(habitat_par_site_bent3_250m_pres_abs)
selection.richesse.complementarite.bent3.500m<-run.rich(habitat_par_site_bent3_500m_pres_abs)
selection.richesse.complementarite.bent3.750m<-run.rich(habitat_par_site_bent3_750m_pres_abs)
selection.richesse.complementarite.bent3.1000m<-run.rich(habitat_par_site_bent3_1000m_pres_abs)
selection.richesse.complementarite.bent3.3000m<-run.rich(habitat_par_site_bent3_3000m_pres_abs)
selection.richesse.complementarite.bent3.5000m<-run.rich(habitat_par_site_bent3_5000m_pres_abs)

#RICHESSE ET DIVERSITE DE POISSON CONTENUE DANS LES SELECTIONS DE SITES#
source('script.diversite.selectionnee.r')

diversite.selection.richesse.complementarite.bent3.30m<-species.richness.evol(selection.richesse.complementarite.bent3.30m,recensement_unique)
diversite.selection.richesse.complementarite.bent3.60m<-species.richness.evol(selection.richesse.complementarite.bent3.60m,recensement_unique)
diversite.selection.richesse.complementarite.bent3.100m<-species.richness.evol(selection.richesse.complementarite.bent3.100m,recensement_unique)
diversite.selection.richesse.complementarite.bent3.250m<-species.richness.evol(selection.richesse.complementarite.bent3.250m,recensement_unique)
diversite.selection.richesse.complementarite.bent3.500m<-species.richness.evol(selection.richesse.complementarite.bent3.500m,recensement_unique)
diversite.selection.richesse.complementarite.bent3.750m<-species.richness.evol(selection.richesse.complementarite.bent3.750m,recensement_unique)
diversite.selection.richesse.complementarite.bent3.1000m<-species.richness.evol(selection.richesse.complementarite.bent3.1000m,recensement_unique)
diversite.selection.richesse.complementarite.bent3.3000m<-species.richness.evol(selection.richesse.complementarite.bent3.3000m,recensement_unique)
diversite.selection.richesse.complementarite.bent3.5000m<-species.richness.evol(selection.richesse.complementarite.bent3.5000m,recensement_unique)

#Calcul des SAI
source('script.sai.r')
SAI_richesse.complementarite.bent3.30m<-SAI(diversite.selection.richesse.complementarite.bent3.30m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent3.60m<-SAI(diversite.selection.richesse.complementarite.bent3.60m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent3.100m<-SAI(diversite.selection.richesse.complementarite.bent3.100m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent3.250m<-SAI(diversite.selection.richesse.complementarite.bent3.250m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent3.500m<-SAI(diversite.selection.richesse.complementarite.bent3.500m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent3.750m<-SAI(diversite.selection.richesse.complementarite.bent3.750m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent3.1000m<-SAI(diversite.selection.richesse.complementarite.bent3.1000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent3.3000m<-SAI(diversite.selection.richesse.complementarite.bent3.3000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent3.5000m<-SAI(diversite.selection.richesse.complementarite.bent3.5000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)



############################HABITAT geo1_geo2###########################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances

setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
relation.habitats.sites.geo1_geo2<-read.csv2("fc3_intersect_geo1_geo2.csv")

sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')

relation.habitats.sites.geo1_geo2<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$Code_Site%in%sites.etudies,]

#création d'une nouvelle table site-habitat par distance
#habitat géo1
relation.habitats.sites.geo1_geo2.30m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==30,]
relation.habitats.sites.geo1_geo2.60m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==60,]
relation.habitats.sites.geo1_geo2.100m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==100,]
relation.habitats.sites.geo1_geo2.250m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==250,]
relation.habitats.sites.geo1_geo2.500m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==500,]
relation.habitats.sites.geo1_geo2.750m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==750,]
relation.habitats.sites.geo1_geo2.1000m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==1000,]
relation.habitats.sites.geo1_geo2.3000m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==3000,]
relation.habitats.sites.geo1_geo2.5000m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==5000,]

#preparation des donnée d'entrée

 m=matrix(0,length(relation.habitats.sites.geo1_geo2.30m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2.30m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2.30m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_30m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1_geo2.60m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2.60m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2.60m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_60m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1_geo2.100m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2.100m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2.100m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_100m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1_geo2.250m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2.250m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2.250m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_250m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1_geo2.500m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2.500m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2.500m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_500m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1_geo2.750m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2.750m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2.750m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_750m=as.data.frame(unique(m))

 m=matrix(0,length(relation.habitats.sites.geo1_geo2.1000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2.1000m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2.1000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_1000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1_geo2.3000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2.3000m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2.3000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_3000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1_geo2.5000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2.5000m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2.5000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_5000m=as.data.frame(unique(m))


habitat_par_site_geo1_geo2_30m<-habitat_par_site_geo1_geo2_30m[habitat_par_site_geo1_geo2_30m$Code_Site%in%sites.etudies,]
habitat_par_site_geo1_geo2_60m<-habitat_par_site_geo1_geo2_60m[habitat_par_site_geo1_geo2_60m$Code_Site%in%sites.etudies,]
habitat_par_site_geo1_geo2_100m<-habitat_par_site_geo1_geo2_100m[habitat_par_site_geo1_geo2_100m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_geo2_250m<-habitat_par_site_geo1_geo2_250m[habitat_par_site_geo1_geo2_250m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_geo2_500m<-habitat_par_site_geo1_geo2_500m[habitat_par_site_geo1_geo2_500m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_geo2_750m<-habitat_par_site_geo1_geo2_750m[habitat_par_site_geo1_geo2_750m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_geo2_1000m<-habitat_par_site_geo1_geo2_1000m[habitat_par_site_geo1_geo2_1000m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_geo2_3000m<-habitat_par_site_geo1_geo2_3000m[habitat_par_site_geo1_geo2_3000m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_geo2_5000m<-habitat_par_site_geo1_geo2_5000m[habitat_par_site_geo1_geo2_5000m$Code_Site%in%sites.etudies,]


#création de table de présence absence
#habitat géo1

habitat_par_site_geo1_geo2_30m_pres_abs<-table(habitat_par_site_geo1_geo2_30m$Code_Site,habitat_par_site_geo1_geo2_30m$Code_hab)
habitat_par_site_geo1_geo2_60m_pres_abs<-table(habitat_par_site_geo1_geo2_60m$Code_Site,habitat_par_site_geo1_geo2_60m$Code_hab)
habitat_par_site_geo1_geo2_100m_pres_abs<-table(habitat_par_site_geo1_geo2_100m$Code_Site,habitat_par_site_geo1_geo2_100m$Code_hab)
habitat_par_site_geo1_geo2_250m_pres_abs<-table(habitat_par_site_geo1_geo2_250m$Code_Site,habitat_par_site_geo1_geo2_250m$Code_hab)
habitat_par_site_geo1_geo2_500m_pres_abs<-table(habitat_par_site_geo1_geo2_500m$Code_Site,habitat_par_site_geo1_geo2_500m$Code_hab)
habitat_par_site_geo1_geo2_750m_pres_abs<-table(habitat_par_site_geo1_geo2_750m$Code_Site,habitat_par_site_geo1_geo2_750m$Code_hab)
habitat_par_site_geo1_geo2_1000m_pres_abs<-table(habitat_par_site_geo1_geo2_1000m$Code_Site,habitat_par_site_geo1_geo2_1000m$Code_hab)
habitat_par_site_geo1_geo2_3000m_pres_abs<-table(habitat_par_site_geo1_geo2_3000m$Code_Site,habitat_par_site_geo1_geo2_3000m$Code_hab)
habitat_par_site_geo1_geo2_5000m_pres_abs<-table(habitat_par_site_geo1_geo2_5000m$Code_Site,habitat_par_site_geo1_geo2_5000m$Code_hab)


#selection des sites
setwd(repertoire)
setwd(input)
source('script.richesse.complementarite.r')

selection.richesse.complementarite.geo1_geo2.30m<-run.rich(habitat_par_site_geo1_geo2_30m_pres_abs)
selection.richesse.complementarite.geo1_geo2.60m<-run.rich(habitat_par_site_geo1_geo2_60m_pres_abs)
selection.richesse.complementarite.geo1_geo2.100m<-run.rich(habitat_par_site_geo1_geo2_100m_pres_abs)
selection.richesse.complementarite.geo1_geo2.250m<-run.rich(habitat_par_site_geo1_geo2_250m_pres_abs)
selection.richesse.complementarite.geo1_geo2.500m<-run.rich(habitat_par_site_geo1_geo2_500m_pres_abs)
selection.richesse.complementarite.geo1_geo2.750m<-run.rich(habitat_par_site_geo1_geo2_750m_pres_abs)
selection.richesse.complementarite.geo1_geo2.1000m<-run.rich(habitat_par_site_geo1_geo2_1000m_pres_abs)
selection.richesse.complementarite.geo1_geo2.3000m<-run.rich(habitat_par_site_geo1_geo2_3000m_pres_abs)
selection.richesse.complementarite.geo1_geo2.5000m<-run.rich(habitat_par_site_geo1_geo2_5000m_pres_abs)

#RICHESSE ET DIVERSITE DE POISSON CONTENUE DANS LES SELECTIONS DE SITES#
source('script.diversite.selectionnee.r')

diversite.selection.richesse.complementarite.geo1_geo2.30m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2.30m,recensement_unique)
diversite.selection.richesse.complementarite.geo1_geo2.60m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2.60m,recensement_unique)
diversite.selection.richesse.complementarite.geo1_geo2.100m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2.100m,recensement_unique)
diversite.selection.richesse.complementarite.geo1_geo2.250m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2.250m,recensement_unique)
diversite.selection.richesse.complementarite.geo1_geo2.500m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2.500m,recensement_unique)
diversite.selection.richesse.complementarite.geo1_geo2.750m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2.750m,recensement_unique)
diversite.selection.richesse.complementarite.geo1_geo2.1000m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2.1000m,recensement_unique)
diversite.selection.richesse.complementarite.geo1_geo2.3000m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2.3000m,recensement_unique)
diversite.selection.richesse.complementarite.geo1_geo2.5000m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2.5000m,recensement_unique)

#Calcul des SAI
source('script.sai.r')
SAI_richesse.complementarite.geo1_geo2.30m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2.30m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1_geo2.60m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2.60m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1_geo2.100m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2.100m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1_geo2.250m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2.250m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1_geo2.500m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2.500m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1_geo2.750m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2.750m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1_geo2.1000m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2.1000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1_geo2.3000m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2.3000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1_geo2.5000m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2.5000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)



############################HABITAT geo1_geo2_geo3###########################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances

setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
relation.habitats.sites.geo1_geo2_geo3<-read.csv2("fc3_intersect_geo1_geo2_geo3.csv")

sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')

relation.habitats.sites.geo1_geo2_geo3<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$Code_Site%in%sites.etudies,]

#création d'une nouvelle table site-habitat par distance
#habitat géo1
relation.habitats.sites.geo1_geo2_geo3.30m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==30,]
relation.habitats.sites.geo1_geo2_geo3.60m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==60,]
relation.habitats.sites.geo1_geo2_geo3.100m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==100,]
relation.habitats.sites.geo1_geo2_geo3.250m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==250,]
relation.habitats.sites.geo1_geo2_geo3.500m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==500,]
relation.habitats.sites.geo1_geo2_geo3.750m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==750,]
relation.habitats.sites.geo1_geo2_geo3.1000m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==1000,]
relation.habitats.sites.geo1_geo2_geo3.3000m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==3000,]
relation.habitats.sites.geo1_geo2_geo3.5000m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==5000,]

#preparation des donnée d'entrée

 m=matrix(0,length(relation.habitats.sites.geo1_geo2_geo3.30m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2_geo3.30m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2_geo3.30m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_geo3_30m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1_geo2_geo3.60m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2_geo3.60m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2_geo3.60m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_geo3_60m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1_geo2_geo3.100m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2_geo3.100m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2_geo3.100m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_geo3_100m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1_geo2_geo3.250m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2_geo3.250m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2_geo3.250m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_geo3_250m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1_geo2_geo3.500m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2_geo3.500m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2_geo3.500m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_geo3_500m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1_geo2_geo3.750m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2_geo3.750m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2_geo3.750m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_geo3_750m=as.data.frame(unique(m))

 m=matrix(0,length(relation.habitats.sites.geo1_geo2_geo3.1000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2_geo3.1000m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2_geo3.1000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_geo3_1000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1_geo2_geo3.3000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2_geo3.3000m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2_geo3.3000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_geo3_3000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1_geo2_geo3.5000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2_geo3.5000m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2_geo3.5000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_geo3_5000m=as.data.frame(unique(m))


habitat_par_site_geo1_geo2_geo3_30m<-habitat_par_site_geo1_geo2_geo3_30m[habitat_par_site_geo1_geo2_geo3_30m$Code_Site%in%sites.etudies,]
habitat_par_site_geo1_geo2_geo3_60m<-habitat_par_site_geo1_geo2_geo3_60m[habitat_par_site_geo1_geo2_geo3_60m$Code_Site%in%sites.etudies,]
habitat_par_site_geo1_geo2_geo3_100m<-habitat_par_site_geo1_geo2_geo3_100m[habitat_par_site_geo1_geo2_geo3_100m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_geo2_geo3_250m<-habitat_par_site_geo1_geo2_geo3_250m[habitat_par_site_geo1_geo2_geo3_250m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_geo2_geo3_500m<-habitat_par_site_geo1_geo2_geo3_500m[habitat_par_site_geo1_geo2_geo3_500m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_geo2_geo3_750m<-habitat_par_site_geo1_geo2_geo3_750m[habitat_par_site_geo1_geo2_geo3_750m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_geo2_geo3_1000m<-habitat_par_site_geo1_geo2_geo3_1000m[habitat_par_site_geo1_geo2_geo3_1000m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_geo2_geo3_3000m<-habitat_par_site_geo1_geo2_geo3_3000m[habitat_par_site_geo1_geo2_geo3_3000m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_geo2_geo3_5000m<-habitat_par_site_geo1_geo2_geo3_5000m[habitat_par_site_geo1_geo2_geo3_5000m$Code_Site%in%sites.etudies,]


#création de table de présence absence
#habitat géo1

habitat_par_site_geo1_geo2_geo3_30m_pres_abs<-table(habitat_par_site_geo1_geo2_geo3_30m$Code_Site,habitat_par_site_geo1_geo2_geo3_30m$Code_hab)
habitat_par_site_geo1_geo2_geo3_60m_pres_abs<-table(habitat_par_site_geo1_geo2_geo3_60m$Code_Site,habitat_par_site_geo1_geo2_geo3_60m$Code_hab)
habitat_par_site_geo1_geo2_geo3_100m_pres_abs<-table(habitat_par_site_geo1_geo2_geo3_100m$Code_Site,habitat_par_site_geo1_geo2_geo3_100m$Code_hab)
habitat_par_site_geo1_geo2_geo3_250m_pres_abs<-table(habitat_par_site_geo1_geo2_geo3_250m$Code_Site,habitat_par_site_geo1_geo2_geo3_250m$Code_hab)
habitat_par_site_geo1_geo2_geo3_500m_pres_abs<-table(habitat_par_site_geo1_geo2_geo3_500m$Code_Site,habitat_par_site_geo1_geo2_geo3_500m$Code_hab)
habitat_par_site_geo1_geo2_geo3_750m_pres_abs<-table(habitat_par_site_geo1_geo2_geo3_750m$Code_Site,habitat_par_site_geo1_geo2_geo3_750m$Code_hab)
habitat_par_site_geo1_geo2_geo3_1000m_pres_abs<-table(habitat_par_site_geo1_geo2_geo3_1000m$Code_Site,habitat_par_site_geo1_geo2_geo3_1000m$Code_hab)
habitat_par_site_geo1_geo2_geo3_3000m_pres_abs<-table(habitat_par_site_geo1_geo2_geo3_3000m$Code_Site,habitat_par_site_geo1_geo2_geo3_3000m$Code_hab)
habitat_par_site_geo1_geo2_geo3_5000m_pres_abs<-table(habitat_par_site_geo1_geo2_geo3_5000m$Code_Site,habitat_par_site_geo1_geo2_geo3_5000m$Code_hab)


#selection des sites
setwd(repertoire)
setwd(input)
source('script.richesse.complementarite.r')

selection.richesse.complementarite.geo1_geo2_geo3.30m<-run.rich(habitat_par_site_geo1_geo2_geo3_30m_pres_abs)
selection.richesse.complementarite.geo1_geo2_geo3.60m<-run.rich(habitat_par_site_geo1_geo2_geo3_60m_pres_abs)
selection.richesse.complementarite.geo1_geo2_geo3.100m<-run.rich(habitat_par_site_geo1_geo2_geo3_100m_pres_abs)
selection.richesse.complementarite.geo1_geo2_geo3.250m<-run.rich(habitat_par_site_geo1_geo2_geo3_250m_pres_abs)
selection.richesse.complementarite.geo1_geo2_geo3.500m<-run.rich(habitat_par_site_geo1_geo2_geo3_500m_pres_abs)
selection.richesse.complementarite.geo1_geo2_geo3.750m<-run.rich(habitat_par_site_geo1_geo2_geo3_750m_pres_abs)
selection.richesse.complementarite.geo1_geo2_geo3.1000m<-run.rich(habitat_par_site_geo1_geo2_geo3_1000m_pres_abs)
selection.richesse.complementarite.geo1_geo2_geo3.3000m<-run.rich(habitat_par_site_geo1_geo2_geo3_3000m_pres_abs)
selection.richesse.complementarite.geo1_geo2_geo3.5000m<-run.rich(habitat_par_site_geo1_geo2_geo3_5000m_pres_abs)

#RICHESSE ET DIVERSITE DE POISSON CONTENUE DANS LES SELECTIONS DE SITES#
source('script.diversite.selectionnee.r')

diversite.selection.richesse.complementarite.geo1_geo2_geo3.30m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2_geo3.30m,recensement_unique)
diversite.selection.richesse.complementarite.geo1_geo2_geo3.60m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2_geo3.60m,recensement_unique)
diversite.selection.richesse.complementarite.geo1_geo2_geo3.100m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2_geo3.100m,recensement_unique)
diversite.selection.richesse.complementarite.geo1_geo2_geo3.250m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2_geo3.250m,recensement_unique)
diversite.selection.richesse.complementarite.geo1_geo2_geo3.500m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2_geo3.500m,recensement_unique)
diversite.selection.richesse.complementarite.geo1_geo2_geo3.750m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2_geo3.750m,recensement_unique)
diversite.selection.richesse.complementarite.geo1_geo2_geo3.1000m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2_geo3.1000m,recensement_unique)
diversite.selection.richesse.complementarite.geo1_geo2_geo3.3000m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2_geo3.3000m,recensement_unique)
diversite.selection.richesse.complementarite.geo1_geo2_geo3.5000m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2_geo3.5000m,recensement_unique)

#Calcul des SAI
source('script.sai.r')
SAI_richesse.complementarite.geo1_geo2_geo3.30m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2_geo3.30m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1_geo2_geo3.60m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2_geo3.60m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1_geo2_geo3.100m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2_geo3.100m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1_geo2_geo3.250m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2_geo3.250m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1_geo2_geo3.500m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2_geo3.500m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1_geo2_geo3.750m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2_geo3.750m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1_geo2_geo3.1000m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2_geo3.1000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1_geo2_geo3.3000m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2_geo3.3000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1_geo2_geo3.5000m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2_geo3.5000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)



############################HABITAT g1_g2_g3_b1###########################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances

setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
relation.habitats.sites.g1_g2_g3_b1<-read.csv2("fc3_intersect_g1_g2_g3_b1.csv")

sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')

relation.habitats.sites.g1_g2_g3_b1<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$Code_Site%in%sites.etudies,]

#création d'une nouvelle table site-habitat par distance
#habitat géo1
relation.habitats.sites.g1_g2_g3_b1.30m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==30,]
relation.habitats.sites.g1_g2_g3_b1.60m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==60,]
relation.habitats.sites.g1_g2_g3_b1.100m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==100,]
relation.habitats.sites.g1_g2_g3_b1.250m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==250,]
relation.habitats.sites.g1_g2_g3_b1.500m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==500,]
relation.habitats.sites.g1_g2_g3_b1.750m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==750,]
relation.habitats.sites.g1_g2_g3_b1.1000m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==1000,]
relation.habitats.sites.g1_g2_g3_b1.3000m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==3000,]
relation.habitats.sites.g1_g2_g3_b1.5000m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==5000,]

#preparation des donnée d'entrée

 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1.30m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1.30m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1.30m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_30m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1.60m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1.60m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1.60m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_60m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1.100m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1.100m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1.100m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_100m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1.250m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1.250m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1.250m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_250m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1.500m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1.500m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1.500m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_500m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1.750m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1.750m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1.750m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_750m=as.data.frame(unique(m))

 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1.1000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1.1000m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1.1000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_1000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1.3000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1.3000m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1.3000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_3000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1.5000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1.5000m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1.5000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_5000m=as.data.frame(unique(m))


habitat_par_site_g1_g2_g3_b1_30m<-habitat_par_site_g1_g2_g3_b1_30m[habitat_par_site_g1_g2_g3_b1_30m$Code_Site%in%sites.etudies,]
habitat_par_site_g1_g2_g3_b1_60m<-habitat_par_site_g1_g2_g3_b1_60m[habitat_par_site_g1_g2_g3_b1_60m$Code_Site%in%sites.etudies,]
habitat_par_site_g1_g2_g3_b1_100m<-habitat_par_site_g1_g2_g3_b1_100m[habitat_par_site_g1_g2_g3_b1_100m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_250m<-habitat_par_site_g1_g2_g3_b1_250m[habitat_par_site_g1_g2_g3_b1_250m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_500m<-habitat_par_site_g1_g2_g3_b1_500m[habitat_par_site_g1_g2_g3_b1_500m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_750m<-habitat_par_site_g1_g2_g3_b1_750m[habitat_par_site_g1_g2_g3_b1_750m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_1000m<-habitat_par_site_g1_g2_g3_b1_1000m[habitat_par_site_g1_g2_g3_b1_1000m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_3000m<-habitat_par_site_g1_g2_g3_b1_3000m[habitat_par_site_g1_g2_g3_b1_3000m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_5000m<-habitat_par_site_g1_g2_g3_b1_5000m[habitat_par_site_g1_g2_g3_b1_5000m$Code_Site%in%sites.etudies,]


#création de table de présence absence
#habitat géo1

habitat_par_site_g1_g2_g3_b1_30m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_30m$Code_Site,habitat_par_site_g1_g2_g3_b1_30m$Code_hab)
habitat_par_site_g1_g2_g3_b1_60m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_60m$Code_Site,habitat_par_site_g1_g2_g3_b1_60m$Code_hab)
habitat_par_site_g1_g2_g3_b1_100m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_100m$Code_Site,habitat_par_site_g1_g2_g3_b1_100m$Code_hab)
habitat_par_site_g1_g2_g3_b1_250m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_250m$Code_Site,habitat_par_site_g1_g2_g3_b1_250m$Code_hab)
habitat_par_site_g1_g2_g3_b1_500m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_500m$Code_Site,habitat_par_site_g1_g2_g3_b1_500m$Code_hab)
habitat_par_site_g1_g2_g3_b1_750m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_750m$Code_Site,habitat_par_site_g1_g2_g3_b1_750m$Code_hab)
habitat_par_site_g1_g2_g3_b1_1000m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_1000m$Code_Site,habitat_par_site_g1_g2_g3_b1_1000m$Code_hab)
habitat_par_site_g1_g2_g3_b1_3000m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_3000m$Code_Site,habitat_par_site_g1_g2_g3_b1_3000m$Code_hab)
habitat_par_site_g1_g2_g3_b1_5000m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_5000m$Code_Site,habitat_par_site_g1_g2_g3_b1_5000m$Code_hab)


#selection des sites
setwd(repertoire)
setwd(input)
source('script.richesse.complementarite.r')

selection.richesse.complementarite.g1_g2_g3_b1.30m<-run.rich(habitat_par_site_g1_g2_g3_b1_30m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1.60m<-run.rich(habitat_par_site_g1_g2_g3_b1_60m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1.100m<-run.rich(habitat_par_site_g1_g2_g3_b1_100m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1.250m<-run.rich(habitat_par_site_g1_g2_g3_b1_250m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1.500m<-run.rich(habitat_par_site_g1_g2_g3_b1_500m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1.750m<-run.rich(habitat_par_site_g1_g2_g3_b1_750m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1.1000m<-run.rich(habitat_par_site_g1_g2_g3_b1_1000m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1.3000m<-run.rich(habitat_par_site_g1_g2_g3_b1_3000m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1.5000m<-run.rich(habitat_par_site_g1_g2_g3_b1_5000m_pres_abs)

#RICHESSE ET DIVERSITE DE POISSON CONTENUE DANS LES SELECTIONS DE SITES#
source('script.diversite.selectionnee.r')

diversite.selection.richesse.complementarite.g1_g2_g3_b1.30m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1.30m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1.60m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1.60m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1.100m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1.100m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1.250m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1.250m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1.500m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1.500m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1.750m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1.750m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1.1000m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1.1000m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1.3000m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1.3000m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1.5000m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1.5000m,recensement_unique)

#Calcul des SAI
source('script.sai.r')
SAI_richesse.complementarite.g1_g2_g3_b1.30m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1.30m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1.60m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1.60m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1.100m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1.100m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1.250m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1.250m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1.500m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1.500m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1.750m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1.750m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1.1000m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1.1000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1.3000m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1.3000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1.5000m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1.5000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)




############################HABITAT g1_g2_g3_b1_b2###########################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances

setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
relation.habitats.sites.g1_g2_g3_b1_b2<-read.csv2("fc3_intersect_g1_g2_g3_b1_b2.csv")

sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')

relation.habitats.sites.g1_g2_g3_b1_b2<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$Code_Site%in%sites.etudies,]

#création d'une nouvelle table site-habitat par distance
#habitat géo1
relation.habitats.sites.g1_g2_g3_b1_b2.30m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==30,]
relation.habitats.sites.g1_g2_g3_b1_b2.60m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==60,]
relation.habitats.sites.g1_g2_g3_b1_b2.100m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==100,]
relation.habitats.sites.g1_g2_g3_b1_b2.250m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==250,]
relation.habitats.sites.g1_g2_g3_b1_b2.500m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==500,]
relation.habitats.sites.g1_g2_g3_b1_b2.750m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==750,]
relation.habitats.sites.g1_g2_g3_b1_b2.1000m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==1000,]
relation.habitats.sites.g1_g2_g3_b1_b2.3000m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==3000,]
relation.habitats.sites.g1_g2_g3_b1_b2.5000m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==5000,]

#preparation des donnée d'entrée

 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2.30m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2.30m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2.30m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_30m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2.60m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2.60m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2.60m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_60m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2.100m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2.100m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2.100m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_100m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2.250m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2.250m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2.250m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_250m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2.500m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2.500m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2.500m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_500m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2.750m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2.750m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2.750m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_750m=as.data.frame(unique(m))

 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2.1000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2.1000m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2.1000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_1000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2.3000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2.3000m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2.3000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_3000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2.5000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2.5000m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2.5000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_5000m=as.data.frame(unique(m))


habitat_par_site_g1_g2_g3_b1_b2_30m<-habitat_par_site_g1_g2_g3_b1_b2_30m[habitat_par_site_g1_g2_g3_b1_b2_30m$Code_Site%in%sites.etudies,]
habitat_par_site_g1_g2_g3_b1_b2_60m<-habitat_par_site_g1_g2_g3_b1_b2_60m[habitat_par_site_g1_g2_g3_b1_b2_60m$Code_Site%in%sites.etudies,]
habitat_par_site_g1_g2_g3_b1_b2_100m<-habitat_par_site_g1_g2_g3_b1_b2_100m[habitat_par_site_g1_g2_g3_b1_b2_100m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_b2_250m<-habitat_par_site_g1_g2_g3_b1_b2_250m[habitat_par_site_g1_g2_g3_b1_b2_250m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_b2_500m<-habitat_par_site_g1_g2_g3_b1_b2_500m[habitat_par_site_g1_g2_g3_b1_b2_500m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_b2_750m<-habitat_par_site_g1_g2_g3_b1_b2_750m[habitat_par_site_g1_g2_g3_b1_b2_750m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_b2_1000m<-habitat_par_site_g1_g2_g3_b1_b2_1000m[habitat_par_site_g1_g2_g3_b1_b2_1000m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_b2_3000m<-habitat_par_site_g1_g2_g3_b1_b2_3000m[habitat_par_site_g1_g2_g3_b1_b2_3000m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_b2_5000m<-habitat_par_site_g1_g2_g3_b1_b2_5000m[habitat_par_site_g1_g2_g3_b1_b2_5000m$Code_Site%in%sites.etudies,]


#création de table de présence absence
#habitat géo1

habitat_par_site_g1_g2_g3_b1_b2_30m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_30m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_30m$Code_hab)
habitat_par_site_g1_g2_g3_b1_b2_60m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_60m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_60m$Code_hab)
habitat_par_site_g1_g2_g3_b1_b2_100m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_100m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_100m$Code_hab)
habitat_par_site_g1_g2_g3_b1_b2_250m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_250m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_250m$Code_hab)
habitat_par_site_g1_g2_g3_b1_b2_500m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_500m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_500m$Code_hab)
habitat_par_site_g1_g2_g3_b1_b2_750m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_750m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_750m$Code_hab)
habitat_par_site_g1_g2_g3_b1_b2_1000m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_1000m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_1000m$Code_hab)
habitat_par_site_g1_g2_g3_b1_b2_3000m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_3000m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_3000m$Code_hab)
habitat_par_site_g1_g2_g3_b1_b2_5000m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_5000m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_5000m$Code_hab)


#selection des sites
setwd(repertoire)
setwd(input)
source('script.richesse.complementarite.r')

selection.richesse.complementarite.g1_g2_g3_b1_b2.30m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_30m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1_b2.60m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_60m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1_b2.100m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_100m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1_b2.250m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_250m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1_b2.500m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_500m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1_b2.750m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_750m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1_b2.1000m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_1000m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1_b2.3000m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_3000m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1_b2.5000m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_5000m_pres_abs)

#RICHESSE ET DIVERSITE DE POISSON CONTENUE DANS LES SELECTIONS DE SITES#
source('script.diversite.selectionnee.r')

diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.30m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2.30m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.60m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2.60m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.100m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2.100m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.250m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2.250m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.500m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2.500m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.750m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2.750m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.1000m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2.1000m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.3000m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2.3000m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.5000m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2.5000m,recensement_unique)

#Calcul des SAI
source('script.sai.r')
SAI_richesse.complementarite.g1_g2_g3_b1_b2.30m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.30m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1_b2.60m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.60m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1_b2.100m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.100m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1_b2.250m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.250m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1_b2.500m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.500m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1_b2.750m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.750m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1_b2.1000m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.1000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1_b2.3000m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.3000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1_b2.5000m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.5000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)



############################HABITAT g1_g2_g3_b1_b2_b3###########################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances

setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
relation.habitats.sites.g1_g2_g3_b1_b2_b3<-read.csv2("fc3_intersect_g1_g2_g3_b1_b2_b3.csv")

sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')

relation.habitats.sites.g1_g2_g3_b1_b2_b3<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$Code_Site%in%sites.etudies,]

#création d'une nouvelle table site-habitat par distance
#habitat géo1
relation.habitats.sites.g1_g2_g3_b1_b2_b3.30m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==30,]
relation.habitats.sites.g1_g2_g3_b1_b2_b3.60m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==60,]
relation.habitats.sites.g1_g2_g3_b1_b2_b3.100m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==100,]
relation.habitats.sites.g1_g2_g3_b1_b2_b3.250m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==250,]
relation.habitats.sites.g1_g2_g3_b1_b2_b3.500m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==500,]
relation.habitats.sites.g1_g2_g3_b1_b2_b3.750m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==750,]
relation.habitats.sites.g1_g2_g3_b1_b2_b3.1000m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==1000,]
relation.habitats.sites.g1_g2_g3_b1_b2_b3.3000m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==3000,]
relation.habitats.sites.g1_g2_g3_b1_b2_b3.5000m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==5000,]

#preparation des donnée d'entrée

 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2_b3.30m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2_b3.30m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2_b3.30m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_b3_30m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2_b3.60m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2_b3.60m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2_b3.60m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_b3_60m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2_b3.100m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2_b3.100m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2_b3.100m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_b3_100m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2_b3.250m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2_b3.250m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2_b3.250m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_b3_250m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2_b3.500m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2_b3.500m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2_b3.500m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_b3_500m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2_b3.750m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2_b3.750m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2_b3.750m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_b3_750m=as.data.frame(unique(m))

 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2_b3.1000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2_b3.1000m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2_b3.1000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_b3_1000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2_b3.3000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2_b3.3000m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2_b3.3000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_b3_3000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2_b3.5000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2_b3.5000m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2_b3.5000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_b3_5000m=as.data.frame(unique(m))


habitat_par_site_g1_g2_g3_b1_b2_b3_30m<-habitat_par_site_g1_g2_g3_b1_b2_b3_30m[habitat_par_site_g1_g2_g3_b1_b2_b3_30m$Code_Site%in%sites.etudies,]
habitat_par_site_g1_g2_g3_b1_b2_b3_60m<-habitat_par_site_g1_g2_g3_b1_b2_b3_60m[habitat_par_site_g1_g2_g3_b1_b2_b3_60m$Code_Site%in%sites.etudies,]
habitat_par_site_g1_g2_g3_b1_b2_b3_100m<-habitat_par_site_g1_g2_g3_b1_b2_b3_100m[habitat_par_site_g1_g2_g3_b1_b2_b3_100m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_b2_b3_250m<-habitat_par_site_g1_g2_g3_b1_b2_b3_250m[habitat_par_site_g1_g2_g3_b1_b2_b3_250m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_b2_b3_500m<-habitat_par_site_g1_g2_g3_b1_b2_b3_500m[habitat_par_site_g1_g2_g3_b1_b2_b3_500m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_b2_b3_750m<-habitat_par_site_g1_g2_g3_b1_b2_b3_750m[habitat_par_site_g1_g2_g3_b1_b2_b3_750m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_b2_b3_1000m<-habitat_par_site_g1_g2_g3_b1_b2_b3_1000m[habitat_par_site_g1_g2_g3_b1_b2_b3_1000m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_b2_b3_3000m<-habitat_par_site_g1_g2_g3_b1_b2_b3_3000m[habitat_par_site_g1_g2_g3_b1_b2_b3_3000m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_b2_b3_5000m<-habitat_par_site_g1_g2_g3_b1_b2_b3_5000m[habitat_par_site_g1_g2_g3_b1_b2_b3_5000m$Code_Site%in%sites.etudies,]


#création de table de présence absence
#habitat géo1

habitat_par_site_g1_g2_g3_b1_b2_b3_30m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_b3_30m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_b3_30m$Code_hab)
habitat_par_site_g1_g2_g3_b1_b2_b3_60m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_b3_60m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_b3_60m$Code_hab)
habitat_par_site_g1_g2_g3_b1_b2_b3_100m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_b3_100m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_b3_100m$Code_hab)
habitat_par_site_g1_g2_g3_b1_b2_b3_250m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_b3_250m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_b3_250m$Code_hab)
habitat_par_site_g1_g2_g3_b1_b2_b3_500m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_b3_500m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_b3_500m$Code_hab)
habitat_par_site_g1_g2_g3_b1_b2_b3_750m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_b3_750m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_b3_750m$Code_hab)
habitat_par_site_g1_g2_g3_b1_b2_b3_1000m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_b3_1000m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_b3_1000m$Code_hab)
habitat_par_site_g1_g2_g3_b1_b2_b3_3000m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_b3_3000m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_b3_3000m$Code_hab)
habitat_par_site_g1_g2_g3_b1_b2_b3_5000m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_b3_5000m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_b3_5000m$Code_hab)


#selection des sites
setwd(repertoire)
setwd(input)
source('script.richesse.complementarite.r')

selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.30m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_b3_30m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.60m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_b3_60m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.100m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_b3_100m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.250m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_b3_250m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.500m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_b3_500m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.750m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_b3_750m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.1000m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_b3_1000m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.3000m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_b3_3000m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.5000m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_b3_5000m_pres_abs)

#RICHESSE ET DIVERSITE DE POISSON CONTENUE DANS LES SELECTIONS DE SITES#
source('script.diversite.selectionnee.r')

diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.30m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.30m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.60m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.60m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.100m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.100m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.250m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.250m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.500m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.500m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.750m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.750m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.1000m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.1000m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.3000m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.3000m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.5000m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.5000m,recensement_unique)

#Calcul des SAI
source('script.sai.r')
SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.30m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.30m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.60m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.60m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.100m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.100m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.250m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.250m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.500m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.500m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.750m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.750m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.1000m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.1000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.3000m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.3000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.5000m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.5000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)


#STOCKAGE DES RESULTATS

targeted.station.results[,1]=c(SAI_richesse.complementarite.geo1.30m,SAI_richesse.complementarite.geo1.60m,SAI_richesse.complementarite.geo1.100m,SAI_richesse.complementarite.geo1.250m,
SAI_richesse.complementarite.geo1.500m,SAI_richesse.complementarite.geo1.750m,SAI_richesse.complementarite.geo1.1000m,SAI_richesse.complementarite.geo1.3000m,SAI_richesse.complementarite.geo1.5000m,
SAI_richesse.complementarite.geo2.30m,SAI_richesse.complementarite.geo2.60m,SAI_richesse.complementarite.geo2.100m,SAI_richesse.complementarite.geo2.250m,
SAI_richesse.complementarite.geo2.500m,SAI_richesse.complementarite.geo2.750m,SAI_richesse.complementarite.geo2.1000m,SAI_richesse.complementarite.geo2.3000m,SAI_richesse.complementarite.geo2.5000m,
SAI_richesse.complementarite.geo3.100m,SAI_richesse.complementarite.geo3.250m,
SAI_richesse.complementarite.geo3.500m,SAI_richesse.complementarite.geo3.750m,SAI_richesse.complementarite.geo3.1000m,SAI_richesse.complementarite.geo3.3000m,SAI_richesse.complementarite.geo3.5000m,
SAI_richesse.complementarite.bent1.30m,SAI_richesse.complementarite.bent1.60m,SAI_richesse.complementarite.bent1.100m,SAI_richesse.complementarite.bent1.250m,
SAI_richesse.complementarite.bent1.500m,SAI_richesse.complementarite.bent1.750m,SAI_richesse.complementarite.bent1.1000m,SAI_richesse.complementarite.bent1.3000m,SAI_richesse.complementarite.bent1.5000m,
SAI_richesse.complementarite.bent2.30m,SAI_richesse.complementarite.bent2.60m,SAI_richesse.complementarite.bent2.100m,SAI_richesse.complementarite.bent2.250m,
SAI_richesse.complementarite.bent2.500m,SAI_richesse.complementarite.bent2.750m,SAI_richesse.complementarite.bent2.1000m,SAI_richesse.complementarite.bent2.3000m,SAI_richesse.complementarite.bent2.5000m,
SAI_richesse.complementarite.bent3.30m,SAI_richesse.complementarite.bent3.60m,SAI_richesse.complementarite.bent3.100m,SAI_richesse.complementarite.bent3.250m,
SAI_richesse.complementarite.bent3.500m,SAI_richesse.complementarite.bent3.750m,SAI_richesse.complementarite.bent3.1000m,SAI_richesse.complementarite.bent3.3000m,SAI_richesse.complementarite.bent3.5000m,
SAI_richesse.complementarite.geo1_geo2.30m,SAI_richesse.complementarite.geo1_geo2.60m,SAI_richesse.complementarite.geo1_geo2.100m,SAI_richesse.complementarite.geo1_geo2.250m,
SAI_richesse.complementarite.geo1_geo2.500m,SAI_richesse.complementarite.geo1_geo2.750m,SAI_richesse.complementarite.geo1_geo2.1000m,SAI_richesse.complementarite.geo1_geo2.3000m,SAI_richesse.complementarite.geo1_geo2.5000m,
SAI_richesse.complementarite.geo1_geo2_geo3.30m,SAI_richesse.complementarite.geo1_geo2_geo3.60m,SAI_richesse.complementarite.geo1_geo2_geo3.100m,SAI_richesse.complementarite.geo1_geo2_geo3.250m,
SAI_richesse.complementarite.geo1_geo2_geo3.500m,SAI_richesse.complementarite.geo1_geo2_geo3.750m,SAI_richesse.complementarite.geo1_geo2_geo3.1000m,SAI_richesse.complementarite.geo1_geo2_geo3.3000m,SAI_richesse.complementarite.geo1_geo2_geo3.5000m,
SAI_richesse.complementarite.g1_g2_g3_b1.30m,SAI_richesse.complementarite.g1_g2_g3_b1.60m,SAI_richesse.complementarite.g1_g2_g3_b1.100m,SAI_richesse.complementarite.g1_g2_g3_b1.250m,
SAI_richesse.complementarite.g1_g2_g3_b1.500m,SAI_richesse.complementarite.g1_g2_g3_b1.750m,SAI_richesse.complementarite.g1_g2_g3_b1.1000m,SAI_richesse.complementarite.g1_g2_g3_b1.3000m,SAI_richesse.complementarite.g1_g2_g3_b1.5000m,
SAI_richesse.complementarite.g1_g2_g3_b1_b2.30m,SAI_richesse.complementarite.g1_g2_g3_b1_b2.60m,SAI_richesse.complementarite.g1_g2_g3_b1_b2.100m,SAI_richesse.complementarite.g1_g2_g3_b1_b2.250m,
SAI_richesse.complementarite.g1_g2_g3_b1_b2.500m,SAI_richesse.complementarite.g1_g2_g3_b1_b2.750m,SAI_richesse.complementarite.g1_g2_g3_b1_b2.1000m,SAI_richesse.complementarite.g1_g2_g3_b1_b2.3000m,SAI_richesse.complementarite.g1_g2_g3_b1_b2.5000m,
SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.30m,SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.60m,SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.100m,SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.250m,
SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.500m,SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.750m,SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.1000m,SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.3000m,SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.5000m)









#-1 STATION TESTS


 sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')


combinaisons=matrix(0,99,5)
combinaisons[1,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[2,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[3,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[4,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[5,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[6,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[7,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[8,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[9,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[10,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[11,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[12,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[13,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[14,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[15,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[16,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[17,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[18,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[19,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[20,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[21,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[22,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[23,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[24,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[25,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[26,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[27,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[28,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[29,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[30,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[31,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[32,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[33,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[34,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[35,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[36,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[37,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[38,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[39,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[40,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[41,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[42,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[43,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[44,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[45,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[46,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[47,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[48,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[49,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[50,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[51,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[52,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[53,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[54,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[55,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[56,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[57,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[58,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[59,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[60,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[61,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[62,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[63,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[64,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[65,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[66,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[67,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[68,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[69,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[70,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[71,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[72,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[73,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[74,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[75,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[76,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[77,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[78,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[79,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[80,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[81,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[82,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[83,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[84,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[85,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[86,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[87,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[88,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[89,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[90,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[91,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[92,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[93,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[94,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[95,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[96,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[97,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[98,]<-sample(sites.etudies,5,replace=FALSE)
combinaisons[99,]<-sample(sites.etudies,5,replace=FALSE)

 for (p in 1:dim(combinaisons)[1]){

suppressed.station<-combinaisons[p,]
                                        #Load du tableau : "site" "poissons"

setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
recensement=read.csv2("REQ_site_poisson.IRD.csv")

#j'enlève la station supprimée
                             
recensement<-recensement[!recensement$Code_Site%in%suppressed.station,]
recensement$Code_Site<-factor(recensement$Code_Site)

#J'enlève les doublons
recensement_unique=(unique(recensement))
 recensement_unique$Code_Site<-factor(recensement_unique$Code_Site)
 
#création d'un tableau de présence absence
tablo_pres_abs<-table(recensement_unique$Code_Site,recensement_unique$Code_poisson)
 

#calcul du nombre d'espèce par station
vec_poisson<-apply(tablo_pres_abs,1,sum)

#selection aléatoire d'un pool de site
setwd(repertoire)
setwd(input)
source('script.aleatoire.r')
selection.aleatoire<-matrix(0,dim(tablo_pres_abs)[1],999) #999
for (i in 1:999){ #999
selection.aleatoire[,i]<-run.alea(tablo_pres_abs)}


#Selection des dites selon un scenario de richesse complementarité basé sur les poissons
setwd(repertoire)
setwd(input)
source('script.richesse.complementarite.r')

selection.richesse.complementarite.poisson<-run.rich(tablo_pres_abs)

#Richesse de poisson incluse dans les sites selectionnés (aléatoirement et selon diversité de poisson)
source('script.diversite.selectionnee.r')
diversite.selection.aleatoire<-matrix(0,dim(tablo_pres_abs)[1],999)    #999
 diversite.selection.aleatoire<-apply(as.matrix(selection.aleatoire),2,function(x) species.richness.evol(x,recensement_unique))
 diversite.moyenne.selection.aleatoire<-apply(diversite.selection.aleatoire,1,mean)


diversite.selection.richesse.complementarite.poisson<-species.richness.evol(selection.richesse.complementarite.poisson,recensement_unique)

#Calcul des courbes aléatore moyenne, uper, et low
diversite.moyenne.selection.aleatoire
diversite.borne.sup.selection.aleatoire=apply(diversite.selection.aleatoire,1,function(x) sort(x)[ceiling(0.975*length(x))])
diversite.borne.inf.selection.aleatoire=apply(diversite.selection.aleatoire,1,function(x) sort(x)[ceiling(0.025*length(x))])

#nombre total d'espèces de poissons observés
nombre.poissons.total <- nlevels(factor(recensement_unique$Code_poisson))

# fonction pourcentage
percentage=function(total,truc){
  truc * 100 / total
  }

# Pourcentages, RANDOM, upper, low
diversite.moyenne.selection.aleatoire.pourcentage <- c(0,percentage(nombre.poissons.total,diversite.moyenne.selection.aleatoire))
diversite.borne.sup.selection.aleatoire.pourcentage <- c(0,percentage(nombre.poissons.total,diversite.borne.sup.selection.aleatoire))
diversite.borne.inf.selection.aleatoire.pourcentage <- c(0,percentage(nombre.poissons.total,diversite.borne.inf.selection.aleatoire))

diversite.selection.richesse.complementarite.poisson.pourcentage<- c(0,percentage(nombre.poissons.total,diversite.selection.richesse.complementarite.poisson))




############################HABITAT geo1###########################################

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
#habitat géo1
relation.habitats.sites.geo1.30m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==30,]
relation.habitats.sites.geo1.60m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==60,]
relation.habitats.sites.geo1.100m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==100,]
relation.habitats.sites.geo1.250m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==250,]
relation.habitats.sites.geo1.500m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==500,]
relation.habitats.sites.geo1.750m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==750,]
relation.habitats.sites.geo1.1000m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==1000,]
relation.habitats.sites.geo1.3000m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==3000,]
relation.habitats.sites.geo1.5000m<-relation.habitats.sites.geo1[relation.habitats.sites.geo1$distance==5000,]

#preparation des donnée d'entrée

 m=matrix(0,length(relation.habitats.sites.geo1.30m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1.30m$Code_Site)
  m[,2]=relation.habitats.sites.geo1.30m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_30m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1.60m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1.60m$Code_Site)
  m[,2]=relation.habitats.sites.geo1.60m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_60m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1.100m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1.100m$Code_Site)
  m[,2]=relation.habitats.sites.geo1.100m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_100m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1.250m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1.250m$Code_Site)
  m[,2]=relation.habitats.sites.geo1.250m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_250m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1.500m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1.500m$Code_Site)
  m[,2]=relation.habitats.sites.geo1.500m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_500m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1.750m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1.750m$Code_Site)
  m[,2]=relation.habitats.sites.geo1.750m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_750m=as.data.frame(unique(m))

 m=matrix(0,length(relation.habitats.sites.geo1.1000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1.1000m$Code_Site)
  m[,2]=relation.habitats.sites.geo1.1000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_1000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1.3000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1.3000m$Code_Site)
  m[,2]=relation.habitats.sites.geo1.3000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_3000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1.5000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1.5000m$Code_Site)
  m[,2]=relation.habitats.sites.geo1.5000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_5000m=as.data.frame(unique(m))


habitat_par_site_geo1_30m<-habitat_par_site_geo1_30m[habitat_par_site_geo1_30m$Code_Site%in%sites.etudies,]
habitat_par_site_geo1_60m<-habitat_par_site_geo1_60m[habitat_par_site_geo1_60m$Code_Site%in%sites.etudies,]
habitat_par_site_geo1_100m<-habitat_par_site_geo1_100m[habitat_par_site_geo1_100m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_250m<-habitat_par_site_geo1_250m[habitat_par_site_geo1_250m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_500m<-habitat_par_site_geo1_500m[habitat_par_site_geo1_500m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_750m<-habitat_par_site_geo1_750m[habitat_par_site_geo1_750m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_1000m<-habitat_par_site_geo1_1000m[habitat_par_site_geo1_1000m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_3000m<-habitat_par_site_geo1_3000m[habitat_par_site_geo1_3000m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_5000m<-habitat_par_site_geo1_5000m[habitat_par_site_geo1_5000m$Code_Site%in%sites.etudies,]


#création de table de présence absence
#habitat géo1

habitat_par_site_geo1_30m_pres_abs<-table(habitat_par_site_geo1_30m$Code_Site,habitat_par_site_geo1_30m$Code_hab)
habitat_par_site_geo1_60m_pres_abs<-table(habitat_par_site_geo1_60m$Code_Site,habitat_par_site_geo1_60m$Code_hab)
habitat_par_site_geo1_100m_pres_abs<-table(habitat_par_site_geo1_100m$Code_Site,habitat_par_site_geo1_100m$Code_hab)
habitat_par_site_geo1_250m_pres_abs<-table(habitat_par_site_geo1_250m$Code_Site,habitat_par_site_geo1_250m$Code_hab)
habitat_par_site_geo1_500m_pres_abs<-table(habitat_par_site_geo1_500m$Code_Site,habitat_par_site_geo1_500m$Code_hab)
habitat_par_site_geo1_750m_pres_abs<-table(habitat_par_site_geo1_750m$Code_Site,habitat_par_site_geo1_750m$Code_hab)
habitat_par_site_geo1_1000m_pres_abs<-table(habitat_par_site_geo1_1000m$Code_Site,habitat_par_site_geo1_1000m$Code_hab)
habitat_par_site_geo1_3000m_pres_abs<-table(habitat_par_site_geo1_3000m$Code_Site,habitat_par_site_geo1_3000m$Code_hab)
habitat_par_site_geo1_5000m_pres_abs<-table(habitat_par_site_geo1_5000m$Code_Site,habitat_par_site_geo1_5000m$Code_hab)


#selection des sites
setwd(repertoire)
setwd(input)
source('script.richesse.complementarite.r')

selection.richesse.complementarite.geo1.30m<-run.rich(habitat_par_site_geo1_30m_pres_abs)
selection.richesse.complementarite.geo1.60m<-run.rich(habitat_par_site_geo1_60m_pres_abs)
selection.richesse.complementarite.geo1.100m<-run.rich(habitat_par_site_geo1_100m_pres_abs)
selection.richesse.complementarite.geo1.250m<-run.rich(habitat_par_site_geo1_250m_pres_abs)
selection.richesse.complementarite.geo1.500m<-run.rich(habitat_par_site_geo1_500m_pres_abs)
selection.richesse.complementarite.geo1.750m<-run.rich(habitat_par_site_geo1_750m_pres_abs)
selection.richesse.complementarite.geo1.1000m<-run.rich(habitat_par_site_geo1_1000m_pres_abs)
selection.richesse.complementarite.geo1.3000m<-run.rich(habitat_par_site_geo1_3000m_pres_abs)
selection.richesse.complementarite.geo1.5000m<-run.rich(habitat_par_site_geo1_5000m_pres_abs)

#RICHESSE ET DIVERSITE DE POISSON CONTENUE DANS LES SELECTIONS DE SITES#
source('script.diversite.selectionnee.r')

diversite.selection.richesse.complementarite.geo1.30m<-species.richness.evol(selection.richesse.complementarite.geo1.30m,recensement_unique)
diversite.selection.richesse.complementarite.geo1.60m<-species.richness.evol(selection.richesse.complementarite.geo1.60m,recensement_unique)
diversite.selection.richesse.complementarite.geo1.100m<-species.richness.evol(selection.richesse.complementarite.geo1.100m,recensement_unique)
diversite.selection.richesse.complementarite.geo1.250m<-species.richness.evol(selection.richesse.complementarite.geo1.250m,recensement_unique)
diversite.selection.richesse.complementarite.geo1.500m<-species.richness.evol(selection.richesse.complementarite.geo1.500m,recensement_unique)
diversite.selection.richesse.complementarite.geo1.750m<-species.richness.evol(selection.richesse.complementarite.geo1.750m,recensement_unique)
diversite.selection.richesse.complementarite.geo1.1000m<-species.richness.evol(selection.richesse.complementarite.geo1.1000m,recensement_unique)
diversite.selection.richesse.complementarite.geo1.3000m<-species.richness.evol(selection.richesse.complementarite.geo1.3000m,recensement_unique)
diversite.selection.richesse.complementarite.geo1.5000m<-species.richness.evol(selection.richesse.complementarite.geo1.5000m,recensement_unique)

#Calcul des SAI
source('script.sai.r')
SAI_richesse.complementarite.geo1.30m<-SAI(diversite.selection.richesse.complementarite.geo1.30m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1.60m<-SAI(diversite.selection.richesse.complementarite.geo1.60m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1.100m<-SAI(diversite.selection.richesse.complementarite.geo1.100m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1.250m<-SAI(diversite.selection.richesse.complementarite.geo1.250m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1.500m<-SAI(diversite.selection.richesse.complementarite.geo1.500m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1.750m<-SAI(diversite.selection.richesse.complementarite.geo1.750m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1.1000m<-SAI(diversite.selection.richesse.complementarite.geo1.1000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1.3000m<-SAI(diversite.selection.richesse.complementarite.geo1.3000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1.5000m<-SAI(diversite.selection.richesse.complementarite.geo1.5000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)



############################HABITAT geo2###########################################

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
#habitat géo1
relation.habitats.sites.geo2.30m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==30,]
relation.habitats.sites.geo2.60m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==60,]
relation.habitats.sites.geo2.100m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==100,]
relation.habitats.sites.geo2.250m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==250,]
relation.habitats.sites.geo2.500m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==500,]
relation.habitats.sites.geo2.750m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==750,]
relation.habitats.sites.geo2.1000m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==1000,]
relation.habitats.sites.geo2.3000m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==3000,]
relation.habitats.sites.geo2.5000m<-relation.habitats.sites.geo2[relation.habitats.sites.geo2$distance==5000,]

#preparation des donnée d'entrée

 m=matrix(0,length(relation.habitats.sites.geo2.30m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo2.30m$Code_Site)
  m[,2]=relation.habitats.sites.geo2.30m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo2_30m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo2.60m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo2.60m$Code_Site)
  m[,2]=relation.habitats.sites.geo2.60m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo2_60m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo2.100m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo2.100m$Code_Site)
  m[,2]=relation.habitats.sites.geo2.100m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo2_100m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo2.250m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo2.250m$Code_Site)
  m[,2]=relation.habitats.sites.geo2.250m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo2_250m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo2.500m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo2.500m$Code_Site)
  m[,2]=relation.habitats.sites.geo2.500m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo2_500m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo2.750m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo2.750m$Code_Site)
  m[,2]=relation.habitats.sites.geo2.750m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo2_750m=as.data.frame(unique(m))

 m=matrix(0,length(relation.habitats.sites.geo2.1000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo2.1000m$Code_Site)
  m[,2]=relation.habitats.sites.geo2.1000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo2_1000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo2.3000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo2.3000m$Code_Site)
  m[,2]=relation.habitats.sites.geo2.3000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo2_3000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo2.5000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo2.5000m$Code_Site)
  m[,2]=relation.habitats.sites.geo2.5000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo2_5000m=as.data.frame(unique(m))


habitat_par_site_geo2_30m<-habitat_par_site_geo2_30m[habitat_par_site_geo2_30m$Code_Site%in%sites.etudies,]
habitat_par_site_geo2_60m<-habitat_par_site_geo2_60m[habitat_par_site_geo2_60m$Code_Site%in%sites.etudies,]
habitat_par_site_geo2_100m<-habitat_par_site_geo2_100m[habitat_par_site_geo2_100m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo2_250m<-habitat_par_site_geo2_250m[habitat_par_site_geo2_250m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo2_500m<-habitat_par_site_geo2_500m[habitat_par_site_geo2_500m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo2_750m<-habitat_par_site_geo2_750m[habitat_par_site_geo2_750m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo2_1000m<-habitat_par_site_geo2_1000m[habitat_par_site_geo2_1000m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo2_3000m<-habitat_par_site_geo2_3000m[habitat_par_site_geo2_3000m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo2_5000m<-habitat_par_site_geo2_5000m[habitat_par_site_geo2_5000m$Code_Site%in%sites.etudies,]


#création de table de présence absence
#habitat géo1

habitat_par_site_geo2_30m_pres_abs<-table(habitat_par_site_geo2_30m$Code_Site,habitat_par_site_geo2_30m$Code_hab)
habitat_par_site_geo2_60m_pres_abs<-table(habitat_par_site_geo2_60m$Code_Site,habitat_par_site_geo2_60m$Code_hab)
habitat_par_site_geo2_100m_pres_abs<-table(habitat_par_site_geo2_100m$Code_Site,habitat_par_site_geo2_100m$Code_hab)
habitat_par_site_geo2_250m_pres_abs<-table(habitat_par_site_geo2_250m$Code_Site,habitat_par_site_geo2_250m$Code_hab)
habitat_par_site_geo2_500m_pres_abs<-table(habitat_par_site_geo2_500m$Code_Site,habitat_par_site_geo2_500m$Code_hab)
habitat_par_site_geo2_750m_pres_abs<-table(habitat_par_site_geo2_750m$Code_Site,habitat_par_site_geo2_750m$Code_hab)
habitat_par_site_geo2_1000m_pres_abs<-table(habitat_par_site_geo2_1000m$Code_Site,habitat_par_site_geo2_1000m$Code_hab)
habitat_par_site_geo2_3000m_pres_abs<-table(habitat_par_site_geo2_3000m$Code_Site,habitat_par_site_geo2_3000m$Code_hab)
habitat_par_site_geo2_5000m_pres_abs<-table(habitat_par_site_geo2_5000m$Code_Site,habitat_par_site_geo2_5000m$Code_hab)


#selection des sites
setwd(repertoire)
setwd(input)
source('script.richesse.complementarite.r')

selection.richesse.complementarite.geo2.30m<-run.rich(habitat_par_site_geo2_30m_pres_abs)
selection.richesse.complementarite.geo2.60m<-run.rich(habitat_par_site_geo2_60m_pres_abs)
selection.richesse.complementarite.geo2.100m<-run.rich(habitat_par_site_geo2_100m_pres_abs)
selection.richesse.complementarite.geo2.250m<-run.rich(habitat_par_site_geo2_250m_pres_abs)
selection.richesse.complementarite.geo2.500m<-run.rich(habitat_par_site_geo2_500m_pres_abs)
selection.richesse.complementarite.geo2.750m<-run.rich(habitat_par_site_geo2_750m_pres_abs)
selection.richesse.complementarite.geo2.1000m<-run.rich(habitat_par_site_geo2_1000m_pres_abs)
selection.richesse.complementarite.geo2.3000m<-run.rich(habitat_par_site_geo2_3000m_pres_abs)
selection.richesse.complementarite.geo2.5000m<-run.rich(habitat_par_site_geo2_5000m_pres_abs)

#RICHESSE ET DIVERSITE DE POISSON CONTENUE DANS LES SELECTIONS DE SITES#
source('script.diversite.selectionnee.r')

diversite.selection.richesse.complementarite.geo2.30m<-species.richness.evol(selection.richesse.complementarite.geo2.30m,recensement_unique)
diversite.selection.richesse.complementarite.geo2.60m<-species.richness.evol(selection.richesse.complementarite.geo2.60m,recensement_unique)
diversite.selection.richesse.complementarite.geo2.100m<-species.richness.evol(selection.richesse.complementarite.geo2.100m,recensement_unique)
diversite.selection.richesse.complementarite.geo2.250m<-species.richness.evol(selection.richesse.complementarite.geo2.250m,recensement_unique)
diversite.selection.richesse.complementarite.geo2.500m<-species.richness.evol(selection.richesse.complementarite.geo2.500m,recensement_unique)
diversite.selection.richesse.complementarite.geo2.750m<-species.richness.evol(selection.richesse.complementarite.geo2.750m,recensement_unique)
diversite.selection.richesse.complementarite.geo2.1000m<-species.richness.evol(selection.richesse.complementarite.geo2.1000m,recensement_unique)
diversite.selection.richesse.complementarite.geo2.3000m<-species.richness.evol(selection.richesse.complementarite.geo2.3000m,recensement_unique)
diversite.selection.richesse.complementarite.geo2.5000m<-species.richness.evol(selection.richesse.complementarite.geo2.5000m,recensement_unique)

#Calcul des SAI
source('script.sai.r')
SAI_richesse.complementarite.geo2.30m<-SAI(diversite.selection.richesse.complementarite.geo2.30m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo2.60m<-SAI(diversite.selection.richesse.complementarite.geo2.60m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo2.100m<-SAI(diversite.selection.richesse.complementarite.geo2.100m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo2.250m<-SAI(diversite.selection.richesse.complementarite.geo2.250m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo2.500m<-SAI(diversite.selection.richesse.complementarite.geo2.500m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo2.750m<-SAI(diversite.selection.richesse.complementarite.geo2.750m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo2.1000m<-SAI(diversite.selection.richesse.complementarite.geo2.1000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo2.3000m<-SAI(diversite.selection.richesse.complementarite.geo2.3000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo2.5000m<-SAI(diversite.selection.richesse.complementarite.geo2.5000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)



############################HABITAT geo3###########################################

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
#habitat géo1

relation.habitats.sites.geo3.100m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==100,]
relation.habitats.sites.geo3.250m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==250,]
relation.habitats.sites.geo3.500m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==500,]
relation.habitats.sites.geo3.750m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==750,]
relation.habitats.sites.geo3.1000m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==1000,]
relation.habitats.sites.geo3.3000m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==3000,]
relation.habitats.sites.geo3.5000m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==5000,]

#preparation des donnée d'entrée



 m=matrix(0,length(relation.habitats.sites.geo3.100m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo3.100m$Code_Site)
  m[,2]=relation.habitats.sites.geo3.100m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo3_100m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo3.250m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo3.250m$Code_Site)
  m[,2]=relation.habitats.sites.geo3.250m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo3_250m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo3.500m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo3.500m$Code_Site)
  m[,2]=relation.habitats.sites.geo3.500m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo3_500m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo3.750m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo3.750m$Code_Site)
  m[,2]=relation.habitats.sites.geo3.750m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo3_750m=as.data.frame(unique(m))

 m=matrix(0,length(relation.habitats.sites.geo3.1000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo3.1000m$Code_Site)
  m[,2]=relation.habitats.sites.geo3.1000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo3_1000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo3.3000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo3.3000m$Code_Site)
  m[,2]=relation.habitats.sites.geo3.3000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo3_3000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo3.5000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo3.5000m$Code_Site)
  m[,2]=relation.habitats.sites.geo3.5000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo3_5000m=as.data.frame(unique(m))


habitat_par_site_geo3_100m<-habitat_par_site_geo3_100m[habitat_par_site_geo3_100m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo3_250m<-habitat_par_site_geo3_250m[habitat_par_site_geo3_250m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo3_500m<-habitat_par_site_geo3_500m[habitat_par_site_geo3_500m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo3_750m<-habitat_par_site_geo3_750m[habitat_par_site_geo3_750m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo3_1000m<-habitat_par_site_geo3_1000m[habitat_par_site_geo3_1000m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo3_3000m<-habitat_par_site_geo3_3000m[habitat_par_site_geo3_3000m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo3_5000m<-habitat_par_site_geo3_5000m[habitat_par_site_geo3_5000m$Code_Site%in%sites.etudies,]


#création de table de présence absence
#habitat géo1

habitat_par_site_geo3_100m_pres_abs<-table(habitat_par_site_geo3_100m$Code_Site,habitat_par_site_geo3_100m$Code_hab)
habitat_par_site_geo3_250m_pres_abs<-table(habitat_par_site_geo3_250m$Code_Site,habitat_par_site_geo3_250m$Code_hab)
habitat_par_site_geo3_500m_pres_abs<-table(habitat_par_site_geo3_500m$Code_Site,habitat_par_site_geo3_500m$Code_hab)
habitat_par_site_geo3_750m_pres_abs<-table(habitat_par_site_geo3_750m$Code_Site,habitat_par_site_geo3_750m$Code_hab)
habitat_par_site_geo3_1000m_pres_abs<-table(habitat_par_site_geo3_1000m$Code_Site,habitat_par_site_geo3_1000m$Code_hab)
habitat_par_site_geo3_3000m_pres_abs<-table(habitat_par_site_geo3_3000m$Code_Site,habitat_par_site_geo3_3000m$Code_hab)
habitat_par_site_geo3_5000m_pres_abs<-table(habitat_par_site_geo3_5000m$Code_Site,habitat_par_site_geo3_5000m$Code_hab)


#selection des sites
setwd(repertoire)
setwd(input)
source('script.richesse.complementarite.r')

selection.richesse.complementarite.geo3.100m<-run.rich(habitat_par_site_geo3_100m_pres_abs)
selection.richesse.complementarite.geo3.250m<-run.rich(habitat_par_site_geo3_250m_pres_abs)
selection.richesse.complementarite.geo3.500m<-run.rich(habitat_par_site_geo3_500m_pres_abs)
selection.richesse.complementarite.geo3.750m<-run.rich(habitat_par_site_geo3_750m_pres_abs)
selection.richesse.complementarite.geo3.1000m<-run.rich(habitat_par_site_geo3_1000m_pres_abs)
selection.richesse.complementarite.geo3.3000m<-run.rich(habitat_par_site_geo3_3000m_pres_abs)
selection.richesse.complementarite.geo3.5000m<-run.rich(habitat_par_site_geo3_5000m_pres_abs)

#RICHESSE ET DIVERSITE DE POISSON CONTENUE DANS LES SELECTIONS DE SITES#
source('script.diversite.selectionnee.r')

diversite.selection.richesse.complementarite.geo3.100m<-species.richness.evol(selection.richesse.complementarite.geo3.100m,recensement_unique)
diversite.selection.richesse.complementarite.geo3.250m<-species.richness.evol(selection.richesse.complementarite.geo3.250m,recensement_unique)
diversite.selection.richesse.complementarite.geo3.500m<-species.richness.evol(selection.richesse.complementarite.geo3.500m,recensement_unique)
diversite.selection.richesse.complementarite.geo3.750m<-species.richness.evol(selection.richesse.complementarite.geo3.750m,recensement_unique)
diversite.selection.richesse.complementarite.geo3.1000m<-species.richness.evol(selection.richesse.complementarite.geo3.1000m,recensement_unique)
diversite.selection.richesse.complementarite.geo3.3000m<-species.richness.evol(selection.richesse.complementarite.geo3.3000m,recensement_unique)
diversite.selection.richesse.complementarite.geo3.5000m<-species.richness.evol(selection.richesse.complementarite.geo3.5000m,recensement_unique)

#Calcul des SAI
source('script.sai.r')
SAI_richesse.complementarite.geo3.100m<-SAI(diversite.selection.richesse.complementarite.geo3.100m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo3.250m<-SAI(diversite.selection.richesse.complementarite.geo3.250m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo3.500m<-SAI(diversite.selection.richesse.complementarite.geo3.500m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo3.750m<-SAI(diversite.selection.richesse.complementarite.geo3.750m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo3.1000m<-SAI(diversite.selection.richesse.complementarite.geo3.1000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo3.3000m<-SAI(diversite.selection.richesse.complementarite.geo3.3000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo3.5000m<-SAI(diversite.selection.richesse.complementarite.geo3.5000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)



############################HABITAT bent1###########################################

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
#habitat géo1
relation.habitats.sites.bent1.30m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==30,]
relation.habitats.sites.bent1.60m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==60,]
relation.habitats.sites.bent1.100m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==100,]
relation.habitats.sites.bent1.250m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==250,]
relation.habitats.sites.bent1.500m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==500,]
relation.habitats.sites.bent1.750m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==750,]
relation.habitats.sites.bent1.1000m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==1000,]
relation.habitats.sites.bent1.3000m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==3000,]
relation.habitats.sites.bent1.5000m<-relation.habitats.sites.bent1[relation.habitats.sites.bent1$distance==5000,]

#preparation des donnée d'entrée

 m=matrix(0,length(relation.habitats.sites.bent1.30m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent1.30m$Code_Site)
  m[,2]=relation.habitats.sites.bent1.30m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent1_30m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent1.60m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent1.60m$Code_Site)
  m[,2]=relation.habitats.sites.bent1.60m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent1_60m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent1.100m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent1.100m$Code_Site)
  m[,2]=relation.habitats.sites.bent1.100m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent1_100m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent1.250m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent1.250m$Code_Site)
  m[,2]=relation.habitats.sites.bent1.250m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent1_250m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent1.500m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent1.500m$Code_Site)
  m[,2]=relation.habitats.sites.bent1.500m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent1_500m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent1.750m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent1.750m$Code_Site)
  m[,2]=relation.habitats.sites.bent1.750m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent1_750m=as.data.frame(unique(m))

 m=matrix(0,length(relation.habitats.sites.bent1.1000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent1.1000m$Code_Site)
  m[,2]=relation.habitats.sites.bent1.1000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent1_1000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent1.3000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent1.3000m$Code_Site)
  m[,2]=relation.habitats.sites.bent1.3000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent1_3000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent1.5000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent1.5000m$Code_Site)
  m[,2]=relation.habitats.sites.bent1.5000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent1_5000m=as.data.frame(unique(m))


habitat_par_site_bent1_30m<-habitat_par_site_bent1_30m[habitat_par_site_bent1_30m$Code_Site%in%sites.etudies,]
habitat_par_site_bent1_60m<-habitat_par_site_bent1_60m[habitat_par_site_bent1_60m$Code_Site%in%sites.etudies,]
habitat_par_site_bent1_100m<-habitat_par_site_bent1_100m[habitat_par_site_bent1_100m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent1_250m<-habitat_par_site_bent1_250m[habitat_par_site_bent1_250m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent1_500m<-habitat_par_site_bent1_500m[habitat_par_site_bent1_500m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent1_750m<-habitat_par_site_bent1_750m[habitat_par_site_bent1_750m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent1_1000m<-habitat_par_site_bent1_1000m[habitat_par_site_bent1_1000m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent1_3000m<-habitat_par_site_bent1_3000m[habitat_par_site_bent1_3000m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent1_5000m<-habitat_par_site_bent1_5000m[habitat_par_site_bent1_5000m$Code_Site%in%sites.etudies,]


#création de table de présence absence
#habitat géo1

habitat_par_site_bent1_30m_pres_abs<-table(habitat_par_site_bent1_30m$Code_Site,habitat_par_site_bent1_30m$Code_hab)
habitat_par_site_bent1_60m_pres_abs<-table(habitat_par_site_bent1_60m$Code_Site,habitat_par_site_bent1_60m$Code_hab)
habitat_par_site_bent1_100m_pres_abs<-table(habitat_par_site_bent1_100m$Code_Site,habitat_par_site_bent1_100m$Code_hab)
habitat_par_site_bent1_250m_pres_abs<-table(habitat_par_site_bent1_250m$Code_Site,habitat_par_site_bent1_250m$Code_hab)
habitat_par_site_bent1_500m_pres_abs<-table(habitat_par_site_bent1_500m$Code_Site,habitat_par_site_bent1_500m$Code_hab)
habitat_par_site_bent1_750m_pres_abs<-table(habitat_par_site_bent1_750m$Code_Site,habitat_par_site_bent1_750m$Code_hab)
habitat_par_site_bent1_1000m_pres_abs<-table(habitat_par_site_bent1_1000m$Code_Site,habitat_par_site_bent1_1000m$Code_hab)
habitat_par_site_bent1_3000m_pres_abs<-table(habitat_par_site_bent1_3000m$Code_Site,habitat_par_site_bent1_3000m$Code_hab)
habitat_par_site_bent1_5000m_pres_abs<-table(habitat_par_site_bent1_5000m$Code_Site,habitat_par_site_bent1_5000m$Code_hab)


#selection des sites
setwd(repertoire)
setwd(input)
source('script.richesse.complementarite.r')

selection.richesse.complementarite.bent1.30m<-run.rich(habitat_par_site_bent1_30m_pres_abs)
selection.richesse.complementarite.bent1.60m<-run.rich(habitat_par_site_bent1_60m_pres_abs)
selection.richesse.complementarite.bent1.100m<-run.rich(habitat_par_site_bent1_100m_pres_abs)
selection.richesse.complementarite.bent1.250m<-run.rich(habitat_par_site_bent1_250m_pres_abs)
selection.richesse.complementarite.bent1.500m<-run.rich(habitat_par_site_bent1_500m_pres_abs)
selection.richesse.complementarite.bent1.750m<-run.rich(habitat_par_site_bent1_750m_pres_abs)
selection.richesse.complementarite.bent1.1000m<-run.rich(habitat_par_site_bent1_1000m_pres_abs)
selection.richesse.complementarite.bent1.3000m<-run.rich(habitat_par_site_bent1_3000m_pres_abs)
selection.richesse.complementarite.bent1.5000m<-run.rich(habitat_par_site_bent1_5000m_pres_abs)

#RICHESSE ET DIVERSITE DE POISSON CONTENUE DANS LES SELECTIONS DE SITES#
source('script.diversite.selectionnee.r')

diversite.selection.richesse.complementarite.bent1.30m<-species.richness.evol(selection.richesse.complementarite.bent1.30m,recensement_unique)
diversite.selection.richesse.complementarite.bent1.60m<-species.richness.evol(selection.richesse.complementarite.bent1.60m,recensement_unique)
diversite.selection.richesse.complementarite.bent1.100m<-species.richness.evol(selection.richesse.complementarite.bent1.100m,recensement_unique)
diversite.selection.richesse.complementarite.bent1.250m<-species.richness.evol(selection.richesse.complementarite.bent1.250m,recensement_unique)
diversite.selection.richesse.complementarite.bent1.500m<-species.richness.evol(selection.richesse.complementarite.bent1.500m,recensement_unique)
diversite.selection.richesse.complementarite.bent1.750m<-species.richness.evol(selection.richesse.complementarite.bent1.750m,recensement_unique)
diversite.selection.richesse.complementarite.bent1.1000m<-species.richness.evol(selection.richesse.complementarite.bent1.1000m,recensement_unique)
diversite.selection.richesse.complementarite.bent1.3000m<-species.richness.evol(selection.richesse.complementarite.bent1.3000m,recensement_unique)
diversite.selection.richesse.complementarite.bent1.5000m<-species.richness.evol(selection.richesse.complementarite.bent1.5000m,recensement_unique)

#Calcul des SAI
source('script.sai.r')
SAI_richesse.complementarite.bent1.30m<-SAI(diversite.selection.richesse.complementarite.bent1.30m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent1.60m<-SAI(diversite.selection.richesse.complementarite.bent1.60m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent1.100m<-SAI(diversite.selection.richesse.complementarite.bent1.100m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent1.250m<-SAI(diversite.selection.richesse.complementarite.bent1.250m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent1.500m<-SAI(diversite.selection.richesse.complementarite.bent1.500m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent1.750m<-SAI(diversite.selection.richesse.complementarite.bent1.750m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent1.1000m<-SAI(diversite.selection.richesse.complementarite.bent1.1000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent1.3000m<-SAI(diversite.selection.richesse.complementarite.bent1.3000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent1.5000m<-SAI(diversite.selection.richesse.complementarite.bent1.5000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)



############################HABITAT bent2###########################################

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
#habitat géo1
relation.habitats.sites.bent2.30m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==30,]
relation.habitats.sites.bent2.60m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==60,]
relation.habitats.sites.bent2.100m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==100,]
relation.habitats.sites.bent2.250m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==250,]
relation.habitats.sites.bent2.500m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==500,]
relation.habitats.sites.bent2.750m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==750,]
relation.habitats.sites.bent2.1000m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==1000,]
relation.habitats.sites.bent2.3000m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==3000,]
relation.habitats.sites.bent2.5000m<-relation.habitats.sites.bent2[relation.habitats.sites.bent2$distance==5000,]

#preparation des donnée d'entrée

 m=matrix(0,length(relation.habitats.sites.bent2.30m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent2.30m$Code_Site)
  m[,2]=relation.habitats.sites.bent2.30m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent2_30m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent2.60m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent2.60m$Code_Site)
  m[,2]=relation.habitats.sites.bent2.60m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent2_60m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent2.100m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent2.100m$Code_Site)
  m[,2]=relation.habitats.sites.bent2.100m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent2_100m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent2.250m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent2.250m$Code_Site)
  m[,2]=relation.habitats.sites.bent2.250m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent2_250m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent2.500m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent2.500m$Code_Site)
  m[,2]=relation.habitats.sites.bent2.500m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent2_500m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent2.750m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent2.750m$Code_Site)
  m[,2]=relation.habitats.sites.bent2.750m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent2_750m=as.data.frame(unique(m))

 m=matrix(0,length(relation.habitats.sites.bent2.1000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent2.1000m$Code_Site)
  m[,2]=relation.habitats.sites.bent2.1000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent2_1000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent2.3000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent2.3000m$Code_Site)
  m[,2]=relation.habitats.sites.bent2.3000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent2_3000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent2.5000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent2.5000m$Code_Site)
  m[,2]=relation.habitats.sites.bent2.5000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent2_5000m=as.data.frame(unique(m))


habitat_par_site_bent2_30m<-habitat_par_site_bent2_30m[habitat_par_site_bent2_30m$Code_Site%in%sites.etudies,]
habitat_par_site_bent2_60m<-habitat_par_site_bent2_60m[habitat_par_site_bent2_60m$Code_Site%in%sites.etudies,]
habitat_par_site_bent2_100m<-habitat_par_site_bent2_100m[habitat_par_site_bent2_100m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent2_250m<-habitat_par_site_bent2_250m[habitat_par_site_bent2_250m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent2_500m<-habitat_par_site_bent2_500m[habitat_par_site_bent2_500m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent2_750m<-habitat_par_site_bent2_750m[habitat_par_site_bent2_750m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent2_1000m<-habitat_par_site_bent2_1000m[habitat_par_site_bent2_1000m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent2_3000m<-habitat_par_site_bent2_3000m[habitat_par_site_bent2_3000m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent2_5000m<-habitat_par_site_bent2_5000m[habitat_par_site_bent2_5000m$Code_Site%in%sites.etudies,]


#création de table de présence absence
#habitat géo1

habitat_par_site_bent2_30m_pres_abs<-table(habitat_par_site_bent2_30m$Code_Site,habitat_par_site_bent2_30m$Code_hab)
habitat_par_site_bent2_60m_pres_abs<-table(habitat_par_site_bent2_60m$Code_Site,habitat_par_site_bent2_60m$Code_hab)
habitat_par_site_bent2_100m_pres_abs<-table(habitat_par_site_bent2_100m$Code_Site,habitat_par_site_bent2_100m$Code_hab)
habitat_par_site_bent2_250m_pres_abs<-table(habitat_par_site_bent2_250m$Code_Site,habitat_par_site_bent2_250m$Code_hab)
habitat_par_site_bent2_500m_pres_abs<-table(habitat_par_site_bent2_500m$Code_Site,habitat_par_site_bent2_500m$Code_hab)
habitat_par_site_bent2_750m_pres_abs<-table(habitat_par_site_bent2_750m$Code_Site,habitat_par_site_bent2_750m$Code_hab)
habitat_par_site_bent2_1000m_pres_abs<-table(habitat_par_site_bent2_1000m$Code_Site,habitat_par_site_bent2_1000m$Code_hab)
habitat_par_site_bent2_3000m_pres_abs<-table(habitat_par_site_bent2_3000m$Code_Site,habitat_par_site_bent2_3000m$Code_hab)
habitat_par_site_bent2_5000m_pres_abs<-table(habitat_par_site_bent2_5000m$Code_Site,habitat_par_site_bent2_5000m$Code_hab)


#selection des sites
setwd(repertoire)
setwd(input)
source('script.richesse.complementarite.r')

selection.richesse.complementarite.bent2.30m<-run.rich(habitat_par_site_bent2_30m_pres_abs)
selection.richesse.complementarite.bent2.60m<-run.rich(habitat_par_site_bent2_60m_pres_abs)
selection.richesse.complementarite.bent2.100m<-run.rich(habitat_par_site_bent2_100m_pres_abs)
selection.richesse.complementarite.bent2.250m<-run.rich(habitat_par_site_bent2_250m_pres_abs)
selection.richesse.complementarite.bent2.500m<-run.rich(habitat_par_site_bent2_500m_pres_abs)
selection.richesse.complementarite.bent2.750m<-run.rich(habitat_par_site_bent2_750m_pres_abs)
selection.richesse.complementarite.bent2.1000m<-run.rich(habitat_par_site_bent2_1000m_pres_abs)
selection.richesse.complementarite.bent2.3000m<-run.rich(habitat_par_site_bent2_3000m_pres_abs)
selection.richesse.complementarite.bent2.5000m<-run.rich(habitat_par_site_bent2_5000m_pres_abs)

#RICHESSE ET DIVERSITE DE POISSON CONTENUE DANS LES SELECTIONS DE SITES#
source('script.diversite.selectionnee.r')

diversite.selection.richesse.complementarite.bent2.30m<-species.richness.evol(selection.richesse.complementarite.bent2.30m,recensement_unique)
diversite.selection.richesse.complementarite.bent2.60m<-species.richness.evol(selection.richesse.complementarite.bent2.60m,recensement_unique)
diversite.selection.richesse.complementarite.bent2.100m<-species.richness.evol(selection.richesse.complementarite.bent2.100m,recensement_unique)
diversite.selection.richesse.complementarite.bent2.250m<-species.richness.evol(selection.richesse.complementarite.bent2.250m,recensement_unique)
diversite.selection.richesse.complementarite.bent2.500m<-species.richness.evol(selection.richesse.complementarite.bent2.500m,recensement_unique)
diversite.selection.richesse.complementarite.bent2.750m<-species.richness.evol(selection.richesse.complementarite.bent2.750m,recensement_unique)
diversite.selection.richesse.complementarite.bent2.1000m<-species.richness.evol(selection.richesse.complementarite.bent2.1000m,recensement_unique)
diversite.selection.richesse.complementarite.bent2.3000m<-species.richness.evol(selection.richesse.complementarite.bent2.3000m,recensement_unique)
diversite.selection.richesse.complementarite.bent2.5000m<-species.richness.evol(selection.richesse.complementarite.bent2.5000m,recensement_unique)

#Calcul des SAI
source('script.sai.r')
SAI_richesse.complementarite.bent2.30m<-SAI(diversite.selection.richesse.complementarite.bent2.30m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent2.60m<-SAI(diversite.selection.richesse.complementarite.bent2.60m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent2.100m<-SAI(diversite.selection.richesse.complementarite.bent2.100m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent2.250m<-SAI(diversite.selection.richesse.complementarite.bent2.250m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent2.500m<-SAI(diversite.selection.richesse.complementarite.bent2.500m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent2.750m<-SAI(diversite.selection.richesse.complementarite.bent2.750m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent2.1000m<-SAI(diversite.selection.richesse.complementarite.bent2.1000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent2.3000m<-SAI(diversite.selection.richesse.complementarite.bent2.3000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent2.5000m<-SAI(diversite.selection.richesse.complementarite.bent2.5000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)



############################HABITAT bent3###########################################

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
#habitat géo1
relation.habitats.sites.bent3.30m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==30,]
relation.habitats.sites.bent3.60m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==60,]
relation.habitats.sites.bent3.100m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==100,]
relation.habitats.sites.bent3.250m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==250,]
relation.habitats.sites.bent3.500m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==500,]
relation.habitats.sites.bent3.750m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==750,]
relation.habitats.sites.bent3.1000m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==1000,]
relation.habitats.sites.bent3.3000m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==3000,]
relation.habitats.sites.bent3.5000m<-relation.habitats.sites.bent3[relation.habitats.sites.bent3$distance==5000,]

#preparation des donnée d'entrée

 m=matrix(0,length(relation.habitats.sites.bent3.30m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent3.30m$Code_Site)
  m[,2]=relation.habitats.sites.bent3.30m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent3_30m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent3.60m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent3.60m$Code_Site)
  m[,2]=relation.habitats.sites.bent3.60m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent3_60m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent3.100m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent3.100m$Code_Site)
  m[,2]=relation.habitats.sites.bent3.100m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent3_100m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent3.250m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent3.250m$Code_Site)
  m[,2]=relation.habitats.sites.bent3.250m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent3_250m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent3.500m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent3.500m$Code_Site)
  m[,2]=relation.habitats.sites.bent3.500m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent3_500m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent3.750m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent3.750m$Code_Site)
  m[,2]=relation.habitats.sites.bent3.750m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent3_750m=as.data.frame(unique(m))

 m=matrix(0,length(relation.habitats.sites.bent3.1000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent3.1000m$Code_Site)
  m[,2]=relation.habitats.sites.bent3.1000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent3_1000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent3.3000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent3.3000m$Code_Site)
  m[,2]=relation.habitats.sites.bent3.3000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent3_3000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.bent3.5000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.bent3.5000m$Code_Site)
  m[,2]=relation.habitats.sites.bent3.5000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_bent3_5000m=as.data.frame(unique(m))


habitat_par_site_bent3_30m<-habitat_par_site_bent3_30m[habitat_par_site_bent3_30m$Code_Site%in%sites.etudies,]
habitat_par_site_bent3_60m<-habitat_par_site_bent3_60m[habitat_par_site_bent3_60m$Code_Site%in%sites.etudies,]
habitat_par_site_bent3_100m<-habitat_par_site_bent3_100m[habitat_par_site_bent3_100m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent3_250m<-habitat_par_site_bent3_250m[habitat_par_site_bent3_250m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent3_500m<-habitat_par_site_bent3_500m[habitat_par_site_bent3_500m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent3_750m<-habitat_par_site_bent3_750m[habitat_par_site_bent3_750m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent3_1000m<-habitat_par_site_bent3_1000m[habitat_par_site_bent3_1000m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent3_3000m<-habitat_par_site_bent3_3000m[habitat_par_site_bent3_3000m$Code_Site%in%sites.etudies,]
 habitat_par_site_bent3_5000m<-habitat_par_site_bent3_5000m[habitat_par_site_bent3_5000m$Code_Site%in%sites.etudies,]


#création de table de présence absence
#habitat géo1

habitat_par_site_bent3_30m_pres_abs<-table(habitat_par_site_bent3_30m$Code_Site,habitat_par_site_bent3_30m$Code_hab)
habitat_par_site_bent3_60m_pres_abs<-table(habitat_par_site_bent3_60m$Code_Site,habitat_par_site_bent3_60m$Code_hab)
habitat_par_site_bent3_100m_pres_abs<-table(habitat_par_site_bent3_100m$Code_Site,habitat_par_site_bent3_100m$Code_hab)
habitat_par_site_bent3_250m_pres_abs<-table(habitat_par_site_bent3_250m$Code_Site,habitat_par_site_bent3_250m$Code_hab)
habitat_par_site_bent3_500m_pres_abs<-table(habitat_par_site_bent3_500m$Code_Site,habitat_par_site_bent3_500m$Code_hab)
habitat_par_site_bent3_750m_pres_abs<-table(habitat_par_site_bent3_750m$Code_Site,habitat_par_site_bent3_750m$Code_hab)
habitat_par_site_bent3_1000m_pres_abs<-table(habitat_par_site_bent3_1000m$Code_Site,habitat_par_site_bent3_1000m$Code_hab)
habitat_par_site_bent3_3000m_pres_abs<-table(habitat_par_site_bent3_3000m$Code_Site,habitat_par_site_bent3_3000m$Code_hab)
habitat_par_site_bent3_5000m_pres_abs<-table(habitat_par_site_bent3_5000m$Code_Site,habitat_par_site_bent3_5000m$Code_hab)


#selection des sites
setwd(repertoire)
setwd(input)
source('script.richesse.complementarite.r')

selection.richesse.complementarite.bent3.30m<-run.rich(habitat_par_site_bent3_30m_pres_abs)
selection.richesse.complementarite.bent3.60m<-run.rich(habitat_par_site_bent3_60m_pres_abs)
selection.richesse.complementarite.bent3.100m<-run.rich(habitat_par_site_bent3_100m_pres_abs)
selection.richesse.complementarite.bent3.250m<-run.rich(habitat_par_site_bent3_250m_pres_abs)
selection.richesse.complementarite.bent3.500m<-run.rich(habitat_par_site_bent3_500m_pres_abs)
selection.richesse.complementarite.bent3.750m<-run.rich(habitat_par_site_bent3_750m_pres_abs)
selection.richesse.complementarite.bent3.1000m<-run.rich(habitat_par_site_bent3_1000m_pres_abs)
selection.richesse.complementarite.bent3.3000m<-run.rich(habitat_par_site_bent3_3000m_pres_abs)
selection.richesse.complementarite.bent3.5000m<-run.rich(habitat_par_site_bent3_5000m_pres_abs)

#RICHESSE ET DIVERSITE DE POISSON CONTENUE DANS LES SELECTIONS DE SITES#
source('script.diversite.selectionnee.r')

diversite.selection.richesse.complementarite.bent3.30m<-species.richness.evol(selection.richesse.complementarite.bent3.30m,recensement_unique)
diversite.selection.richesse.complementarite.bent3.60m<-species.richness.evol(selection.richesse.complementarite.bent3.60m,recensement_unique)
diversite.selection.richesse.complementarite.bent3.100m<-species.richness.evol(selection.richesse.complementarite.bent3.100m,recensement_unique)
diversite.selection.richesse.complementarite.bent3.250m<-species.richness.evol(selection.richesse.complementarite.bent3.250m,recensement_unique)
diversite.selection.richesse.complementarite.bent3.500m<-species.richness.evol(selection.richesse.complementarite.bent3.500m,recensement_unique)
diversite.selection.richesse.complementarite.bent3.750m<-species.richness.evol(selection.richesse.complementarite.bent3.750m,recensement_unique)
diversite.selection.richesse.complementarite.bent3.1000m<-species.richness.evol(selection.richesse.complementarite.bent3.1000m,recensement_unique)
diversite.selection.richesse.complementarite.bent3.3000m<-species.richness.evol(selection.richesse.complementarite.bent3.3000m,recensement_unique)
diversite.selection.richesse.complementarite.bent3.5000m<-species.richness.evol(selection.richesse.complementarite.bent3.5000m,recensement_unique)

#Calcul des SAI
source('script.sai.r')
SAI_richesse.complementarite.bent3.30m<-SAI(diversite.selection.richesse.complementarite.bent3.30m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent3.60m<-SAI(diversite.selection.richesse.complementarite.bent3.60m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent3.100m<-SAI(diversite.selection.richesse.complementarite.bent3.100m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent3.250m<-SAI(diversite.selection.richesse.complementarite.bent3.250m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent3.500m<-SAI(diversite.selection.richesse.complementarite.bent3.500m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent3.750m<-SAI(diversite.selection.richesse.complementarite.bent3.750m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent3.1000m<-SAI(diversite.selection.richesse.complementarite.bent3.1000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent3.3000m<-SAI(diversite.selection.richesse.complementarite.bent3.3000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.bent3.5000m<-SAI(diversite.selection.richesse.complementarite.bent3.5000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)



############################HABITAT geo1_geo2###########################################

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
#habitat géo1
relation.habitats.sites.geo1_geo2.30m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==30,]
relation.habitats.sites.geo1_geo2.60m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==60,]
relation.habitats.sites.geo1_geo2.100m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==100,]
relation.habitats.sites.geo1_geo2.250m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==250,]
relation.habitats.sites.geo1_geo2.500m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==500,]
relation.habitats.sites.geo1_geo2.750m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==750,]
relation.habitats.sites.geo1_geo2.1000m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==1000,]
relation.habitats.sites.geo1_geo2.3000m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==3000,]
relation.habitats.sites.geo1_geo2.5000m<-relation.habitats.sites.geo1_geo2[relation.habitats.sites.geo1_geo2$distance==5000,]

#preparation des donnée d'entrée

 m=matrix(0,length(relation.habitats.sites.geo1_geo2.30m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2.30m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2.30m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_30m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1_geo2.60m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2.60m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2.60m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_60m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1_geo2.100m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2.100m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2.100m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_100m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1_geo2.250m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2.250m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2.250m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_250m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1_geo2.500m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2.500m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2.500m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_500m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1_geo2.750m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2.750m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2.750m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_750m=as.data.frame(unique(m))

 m=matrix(0,length(relation.habitats.sites.geo1_geo2.1000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2.1000m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2.1000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_1000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1_geo2.3000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2.3000m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2.3000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_3000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1_geo2.5000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2.5000m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2.5000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_5000m=as.data.frame(unique(m))


habitat_par_site_geo1_geo2_30m<-habitat_par_site_geo1_geo2_30m[habitat_par_site_geo1_geo2_30m$Code_Site%in%sites.etudies,]
habitat_par_site_geo1_geo2_60m<-habitat_par_site_geo1_geo2_60m[habitat_par_site_geo1_geo2_60m$Code_Site%in%sites.etudies,]
habitat_par_site_geo1_geo2_100m<-habitat_par_site_geo1_geo2_100m[habitat_par_site_geo1_geo2_100m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_geo2_250m<-habitat_par_site_geo1_geo2_250m[habitat_par_site_geo1_geo2_250m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_geo2_500m<-habitat_par_site_geo1_geo2_500m[habitat_par_site_geo1_geo2_500m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_geo2_750m<-habitat_par_site_geo1_geo2_750m[habitat_par_site_geo1_geo2_750m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_geo2_1000m<-habitat_par_site_geo1_geo2_1000m[habitat_par_site_geo1_geo2_1000m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_geo2_3000m<-habitat_par_site_geo1_geo2_3000m[habitat_par_site_geo1_geo2_3000m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_geo2_5000m<-habitat_par_site_geo1_geo2_5000m[habitat_par_site_geo1_geo2_5000m$Code_Site%in%sites.etudies,]


#création de table de présence absence
#habitat géo1

habitat_par_site_geo1_geo2_30m_pres_abs<-table(habitat_par_site_geo1_geo2_30m$Code_Site,habitat_par_site_geo1_geo2_30m$Code_hab)
habitat_par_site_geo1_geo2_60m_pres_abs<-table(habitat_par_site_geo1_geo2_60m$Code_Site,habitat_par_site_geo1_geo2_60m$Code_hab)
habitat_par_site_geo1_geo2_100m_pres_abs<-table(habitat_par_site_geo1_geo2_100m$Code_Site,habitat_par_site_geo1_geo2_100m$Code_hab)
habitat_par_site_geo1_geo2_250m_pres_abs<-table(habitat_par_site_geo1_geo2_250m$Code_Site,habitat_par_site_geo1_geo2_250m$Code_hab)
habitat_par_site_geo1_geo2_500m_pres_abs<-table(habitat_par_site_geo1_geo2_500m$Code_Site,habitat_par_site_geo1_geo2_500m$Code_hab)
habitat_par_site_geo1_geo2_750m_pres_abs<-table(habitat_par_site_geo1_geo2_750m$Code_Site,habitat_par_site_geo1_geo2_750m$Code_hab)
habitat_par_site_geo1_geo2_1000m_pres_abs<-table(habitat_par_site_geo1_geo2_1000m$Code_Site,habitat_par_site_geo1_geo2_1000m$Code_hab)
habitat_par_site_geo1_geo2_3000m_pres_abs<-table(habitat_par_site_geo1_geo2_3000m$Code_Site,habitat_par_site_geo1_geo2_3000m$Code_hab)
habitat_par_site_geo1_geo2_5000m_pres_abs<-table(habitat_par_site_geo1_geo2_5000m$Code_Site,habitat_par_site_geo1_geo2_5000m$Code_hab)


#selection des sites
setwd(repertoire)
setwd(input)
source('script.richesse.complementarite.r')

selection.richesse.complementarite.geo1_geo2.30m<-run.rich(habitat_par_site_geo1_geo2_30m_pres_abs)
selection.richesse.complementarite.geo1_geo2.60m<-run.rich(habitat_par_site_geo1_geo2_60m_pres_abs)
selection.richesse.complementarite.geo1_geo2.100m<-run.rich(habitat_par_site_geo1_geo2_100m_pres_abs)
selection.richesse.complementarite.geo1_geo2.250m<-run.rich(habitat_par_site_geo1_geo2_250m_pres_abs)
selection.richesse.complementarite.geo1_geo2.500m<-run.rich(habitat_par_site_geo1_geo2_500m_pres_abs)
selection.richesse.complementarite.geo1_geo2.750m<-run.rich(habitat_par_site_geo1_geo2_750m_pres_abs)
selection.richesse.complementarite.geo1_geo2.1000m<-run.rich(habitat_par_site_geo1_geo2_1000m_pres_abs)
selection.richesse.complementarite.geo1_geo2.3000m<-run.rich(habitat_par_site_geo1_geo2_3000m_pres_abs)
selection.richesse.complementarite.geo1_geo2.5000m<-run.rich(habitat_par_site_geo1_geo2_5000m_pres_abs)

#RICHESSE ET DIVERSITE DE POISSON CONTENUE DANS LES SELECTIONS DE SITES#
source('script.diversite.selectionnee.r')

diversite.selection.richesse.complementarite.geo1_geo2.30m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2.30m,recensement_unique)
diversite.selection.richesse.complementarite.geo1_geo2.60m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2.60m,recensement_unique)
diversite.selection.richesse.complementarite.geo1_geo2.100m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2.100m,recensement_unique)
diversite.selection.richesse.complementarite.geo1_geo2.250m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2.250m,recensement_unique)
diversite.selection.richesse.complementarite.geo1_geo2.500m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2.500m,recensement_unique)
diversite.selection.richesse.complementarite.geo1_geo2.750m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2.750m,recensement_unique)
diversite.selection.richesse.complementarite.geo1_geo2.1000m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2.1000m,recensement_unique)
diversite.selection.richesse.complementarite.geo1_geo2.3000m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2.3000m,recensement_unique)
diversite.selection.richesse.complementarite.geo1_geo2.5000m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2.5000m,recensement_unique)

#Calcul des SAI
source('script.sai.r')
SAI_richesse.complementarite.geo1_geo2.30m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2.30m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1_geo2.60m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2.60m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1_geo2.100m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2.100m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1_geo2.250m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2.250m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1_geo2.500m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2.500m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1_geo2.750m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2.750m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1_geo2.1000m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2.1000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1_geo2.3000m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2.3000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1_geo2.5000m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2.5000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)



############################HABITAT geo1_geo2_geo3###########################################

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
#habitat géo1
relation.habitats.sites.geo1_geo2_geo3.30m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==30,]
relation.habitats.sites.geo1_geo2_geo3.60m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==60,]
relation.habitats.sites.geo1_geo2_geo3.100m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==100,]
relation.habitats.sites.geo1_geo2_geo3.250m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==250,]
relation.habitats.sites.geo1_geo2_geo3.500m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==500,]
relation.habitats.sites.geo1_geo2_geo3.750m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==750,]
relation.habitats.sites.geo1_geo2_geo3.1000m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==1000,]
relation.habitats.sites.geo1_geo2_geo3.3000m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==3000,]
relation.habitats.sites.geo1_geo2_geo3.5000m<-relation.habitats.sites.geo1_geo2_geo3[relation.habitats.sites.geo1_geo2_geo3$distance==5000,]

#preparation des donnée d'entrée

 m=matrix(0,length(relation.habitats.sites.geo1_geo2_geo3.30m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2_geo3.30m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2_geo3.30m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_geo3_30m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1_geo2_geo3.60m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2_geo3.60m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2_geo3.60m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_geo3_60m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1_geo2_geo3.100m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2_geo3.100m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2_geo3.100m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_geo3_100m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1_geo2_geo3.250m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2_geo3.250m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2_geo3.250m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_geo3_250m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1_geo2_geo3.500m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2_geo3.500m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2_geo3.500m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_geo3_500m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1_geo2_geo3.750m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2_geo3.750m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2_geo3.750m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_geo3_750m=as.data.frame(unique(m))

 m=matrix(0,length(relation.habitats.sites.geo1_geo2_geo3.1000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2_geo3.1000m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2_geo3.1000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_geo3_1000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1_geo2_geo3.3000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2_geo3.3000m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2_geo3.3000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_geo3_3000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo1_geo2_geo3.5000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo1_geo2_geo3.5000m$Code_Site)
  m[,2]=relation.habitats.sites.geo1_geo2_geo3.5000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo1_geo2_geo3_5000m=as.data.frame(unique(m))


habitat_par_site_geo1_geo2_geo3_30m<-habitat_par_site_geo1_geo2_geo3_30m[habitat_par_site_geo1_geo2_geo3_30m$Code_Site%in%sites.etudies,]
habitat_par_site_geo1_geo2_geo3_60m<-habitat_par_site_geo1_geo2_geo3_60m[habitat_par_site_geo1_geo2_geo3_60m$Code_Site%in%sites.etudies,]
habitat_par_site_geo1_geo2_geo3_100m<-habitat_par_site_geo1_geo2_geo3_100m[habitat_par_site_geo1_geo2_geo3_100m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_geo2_geo3_250m<-habitat_par_site_geo1_geo2_geo3_250m[habitat_par_site_geo1_geo2_geo3_250m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_geo2_geo3_500m<-habitat_par_site_geo1_geo2_geo3_500m[habitat_par_site_geo1_geo2_geo3_500m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_geo2_geo3_750m<-habitat_par_site_geo1_geo2_geo3_750m[habitat_par_site_geo1_geo2_geo3_750m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_geo2_geo3_1000m<-habitat_par_site_geo1_geo2_geo3_1000m[habitat_par_site_geo1_geo2_geo3_1000m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_geo2_geo3_3000m<-habitat_par_site_geo1_geo2_geo3_3000m[habitat_par_site_geo1_geo2_geo3_3000m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo1_geo2_geo3_5000m<-habitat_par_site_geo1_geo2_geo3_5000m[habitat_par_site_geo1_geo2_geo3_5000m$Code_Site%in%sites.etudies,]


#création de table de présence absence
#habitat géo1

habitat_par_site_geo1_geo2_geo3_30m_pres_abs<-table(habitat_par_site_geo1_geo2_geo3_30m$Code_Site,habitat_par_site_geo1_geo2_geo3_30m$Code_hab)
habitat_par_site_geo1_geo2_geo3_60m_pres_abs<-table(habitat_par_site_geo1_geo2_geo3_60m$Code_Site,habitat_par_site_geo1_geo2_geo3_60m$Code_hab)
habitat_par_site_geo1_geo2_geo3_100m_pres_abs<-table(habitat_par_site_geo1_geo2_geo3_100m$Code_Site,habitat_par_site_geo1_geo2_geo3_100m$Code_hab)
habitat_par_site_geo1_geo2_geo3_250m_pres_abs<-table(habitat_par_site_geo1_geo2_geo3_250m$Code_Site,habitat_par_site_geo1_geo2_geo3_250m$Code_hab)
habitat_par_site_geo1_geo2_geo3_500m_pres_abs<-table(habitat_par_site_geo1_geo2_geo3_500m$Code_Site,habitat_par_site_geo1_geo2_geo3_500m$Code_hab)
habitat_par_site_geo1_geo2_geo3_750m_pres_abs<-table(habitat_par_site_geo1_geo2_geo3_750m$Code_Site,habitat_par_site_geo1_geo2_geo3_750m$Code_hab)
habitat_par_site_geo1_geo2_geo3_1000m_pres_abs<-table(habitat_par_site_geo1_geo2_geo3_1000m$Code_Site,habitat_par_site_geo1_geo2_geo3_1000m$Code_hab)
habitat_par_site_geo1_geo2_geo3_3000m_pres_abs<-table(habitat_par_site_geo1_geo2_geo3_3000m$Code_Site,habitat_par_site_geo1_geo2_geo3_3000m$Code_hab)
habitat_par_site_geo1_geo2_geo3_5000m_pres_abs<-table(habitat_par_site_geo1_geo2_geo3_5000m$Code_Site,habitat_par_site_geo1_geo2_geo3_5000m$Code_hab)


#selection des sites
setwd(repertoire)
setwd(input)
source('script.richesse.complementarite.r')

selection.richesse.complementarite.geo1_geo2_geo3.30m<-run.rich(habitat_par_site_geo1_geo2_geo3_30m_pres_abs)
selection.richesse.complementarite.geo1_geo2_geo3.60m<-run.rich(habitat_par_site_geo1_geo2_geo3_60m_pres_abs)
selection.richesse.complementarite.geo1_geo2_geo3.100m<-run.rich(habitat_par_site_geo1_geo2_geo3_100m_pres_abs)
selection.richesse.complementarite.geo1_geo2_geo3.250m<-run.rich(habitat_par_site_geo1_geo2_geo3_250m_pres_abs)
selection.richesse.complementarite.geo1_geo2_geo3.500m<-run.rich(habitat_par_site_geo1_geo2_geo3_500m_pres_abs)
selection.richesse.complementarite.geo1_geo2_geo3.750m<-run.rich(habitat_par_site_geo1_geo2_geo3_750m_pres_abs)
selection.richesse.complementarite.geo1_geo2_geo3.1000m<-run.rich(habitat_par_site_geo1_geo2_geo3_1000m_pres_abs)
selection.richesse.complementarite.geo1_geo2_geo3.3000m<-run.rich(habitat_par_site_geo1_geo2_geo3_3000m_pres_abs)
selection.richesse.complementarite.geo1_geo2_geo3.5000m<-run.rich(habitat_par_site_geo1_geo2_geo3_5000m_pres_abs)

#RICHESSE ET DIVERSITE DE POISSON CONTENUE DANS LES SELECTIONS DE SITES#
source('script.diversite.selectionnee.r')

diversite.selection.richesse.complementarite.geo1_geo2_geo3.30m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2_geo3.30m,recensement_unique)
diversite.selection.richesse.complementarite.geo1_geo2_geo3.60m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2_geo3.60m,recensement_unique)
diversite.selection.richesse.complementarite.geo1_geo2_geo3.100m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2_geo3.100m,recensement_unique)
diversite.selection.richesse.complementarite.geo1_geo2_geo3.250m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2_geo3.250m,recensement_unique)
diversite.selection.richesse.complementarite.geo1_geo2_geo3.500m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2_geo3.500m,recensement_unique)
diversite.selection.richesse.complementarite.geo1_geo2_geo3.750m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2_geo3.750m,recensement_unique)
diversite.selection.richesse.complementarite.geo1_geo2_geo3.1000m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2_geo3.1000m,recensement_unique)
diversite.selection.richesse.complementarite.geo1_geo2_geo3.3000m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2_geo3.3000m,recensement_unique)
diversite.selection.richesse.complementarite.geo1_geo2_geo3.5000m<-species.richness.evol(selection.richesse.complementarite.geo1_geo2_geo3.5000m,recensement_unique)

#Calcul des SAI
source('script.sai.r')
SAI_richesse.complementarite.geo1_geo2_geo3.30m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2_geo3.30m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1_geo2_geo3.60m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2_geo3.60m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1_geo2_geo3.100m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2_geo3.100m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1_geo2_geo3.250m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2_geo3.250m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1_geo2_geo3.500m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2_geo3.500m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1_geo2_geo3.750m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2_geo3.750m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1_geo2_geo3.1000m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2_geo3.1000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1_geo2_geo3.3000m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2_geo3.3000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.geo1_geo2_geo3.5000m<-SAI(diversite.selection.richesse.complementarite.geo1_geo2_geo3.5000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)



############################HABITAT g1_g2_g3_b1###########################################

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
#habitat géo1
relation.habitats.sites.g1_g2_g3_b1.30m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==30,]
relation.habitats.sites.g1_g2_g3_b1.60m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==60,]
relation.habitats.sites.g1_g2_g3_b1.100m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==100,]
relation.habitats.sites.g1_g2_g3_b1.250m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==250,]
relation.habitats.sites.g1_g2_g3_b1.500m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==500,]
relation.habitats.sites.g1_g2_g3_b1.750m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==750,]
relation.habitats.sites.g1_g2_g3_b1.1000m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==1000,]
relation.habitats.sites.g1_g2_g3_b1.3000m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==3000,]
relation.habitats.sites.g1_g2_g3_b1.5000m<-relation.habitats.sites.g1_g2_g3_b1[relation.habitats.sites.g1_g2_g3_b1$distance==5000,]

#preparation des donnée d'entrée

 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1.30m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1.30m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1.30m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_30m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1.60m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1.60m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1.60m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_60m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1.100m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1.100m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1.100m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_100m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1.250m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1.250m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1.250m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_250m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1.500m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1.500m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1.500m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_500m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1.750m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1.750m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1.750m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_750m=as.data.frame(unique(m))

 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1.1000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1.1000m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1.1000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_1000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1.3000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1.3000m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1.3000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_3000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1.5000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1.5000m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1.5000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_5000m=as.data.frame(unique(m))


habitat_par_site_g1_g2_g3_b1_30m<-habitat_par_site_g1_g2_g3_b1_30m[habitat_par_site_g1_g2_g3_b1_30m$Code_Site%in%sites.etudies,]
habitat_par_site_g1_g2_g3_b1_60m<-habitat_par_site_g1_g2_g3_b1_60m[habitat_par_site_g1_g2_g3_b1_60m$Code_Site%in%sites.etudies,]
habitat_par_site_g1_g2_g3_b1_100m<-habitat_par_site_g1_g2_g3_b1_100m[habitat_par_site_g1_g2_g3_b1_100m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_250m<-habitat_par_site_g1_g2_g3_b1_250m[habitat_par_site_g1_g2_g3_b1_250m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_500m<-habitat_par_site_g1_g2_g3_b1_500m[habitat_par_site_g1_g2_g3_b1_500m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_750m<-habitat_par_site_g1_g2_g3_b1_750m[habitat_par_site_g1_g2_g3_b1_750m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_1000m<-habitat_par_site_g1_g2_g3_b1_1000m[habitat_par_site_g1_g2_g3_b1_1000m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_3000m<-habitat_par_site_g1_g2_g3_b1_3000m[habitat_par_site_g1_g2_g3_b1_3000m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_5000m<-habitat_par_site_g1_g2_g3_b1_5000m[habitat_par_site_g1_g2_g3_b1_5000m$Code_Site%in%sites.etudies,]


#création de table de présence absence
#habitat géo1

habitat_par_site_g1_g2_g3_b1_30m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_30m$Code_Site,habitat_par_site_g1_g2_g3_b1_30m$Code_hab)
habitat_par_site_g1_g2_g3_b1_60m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_60m$Code_Site,habitat_par_site_g1_g2_g3_b1_60m$Code_hab)
habitat_par_site_g1_g2_g3_b1_100m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_100m$Code_Site,habitat_par_site_g1_g2_g3_b1_100m$Code_hab)
habitat_par_site_g1_g2_g3_b1_250m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_250m$Code_Site,habitat_par_site_g1_g2_g3_b1_250m$Code_hab)
habitat_par_site_g1_g2_g3_b1_500m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_500m$Code_Site,habitat_par_site_g1_g2_g3_b1_500m$Code_hab)
habitat_par_site_g1_g2_g3_b1_750m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_750m$Code_Site,habitat_par_site_g1_g2_g3_b1_750m$Code_hab)
habitat_par_site_g1_g2_g3_b1_1000m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_1000m$Code_Site,habitat_par_site_g1_g2_g3_b1_1000m$Code_hab)
habitat_par_site_g1_g2_g3_b1_3000m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_3000m$Code_Site,habitat_par_site_g1_g2_g3_b1_3000m$Code_hab)
habitat_par_site_g1_g2_g3_b1_5000m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_5000m$Code_Site,habitat_par_site_g1_g2_g3_b1_5000m$Code_hab)


#selection des sites
setwd(repertoire)
setwd(input)
source('script.richesse.complementarite.r')

selection.richesse.complementarite.g1_g2_g3_b1.30m<-run.rich(habitat_par_site_g1_g2_g3_b1_30m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1.60m<-run.rich(habitat_par_site_g1_g2_g3_b1_60m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1.100m<-run.rich(habitat_par_site_g1_g2_g3_b1_100m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1.250m<-run.rich(habitat_par_site_g1_g2_g3_b1_250m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1.500m<-run.rich(habitat_par_site_g1_g2_g3_b1_500m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1.750m<-run.rich(habitat_par_site_g1_g2_g3_b1_750m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1.1000m<-run.rich(habitat_par_site_g1_g2_g3_b1_1000m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1.3000m<-run.rich(habitat_par_site_g1_g2_g3_b1_3000m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1.5000m<-run.rich(habitat_par_site_g1_g2_g3_b1_5000m_pres_abs)

#RICHESSE ET DIVERSITE DE POISSON CONTENUE DANS LES SELECTIONS DE SITES#
source('script.diversite.selectionnee.r')

diversite.selection.richesse.complementarite.g1_g2_g3_b1.30m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1.30m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1.60m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1.60m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1.100m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1.100m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1.250m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1.250m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1.500m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1.500m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1.750m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1.750m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1.1000m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1.1000m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1.3000m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1.3000m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1.5000m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1.5000m,recensement_unique)

#Calcul des SAI
source('script.sai.r')
SAI_richesse.complementarite.g1_g2_g3_b1.30m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1.30m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1.60m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1.60m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1.100m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1.100m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1.250m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1.250m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1.500m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1.500m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1.750m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1.750m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1.1000m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1.1000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1.3000m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1.3000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1.5000m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1.5000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)




############################HABITAT g1_g2_g3_b1_b2###########################################

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
#habitat géo1
relation.habitats.sites.g1_g2_g3_b1_b2.30m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==30,]
relation.habitats.sites.g1_g2_g3_b1_b2.60m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==60,]
relation.habitats.sites.g1_g2_g3_b1_b2.100m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==100,]
relation.habitats.sites.g1_g2_g3_b1_b2.250m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==250,]
relation.habitats.sites.g1_g2_g3_b1_b2.500m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==500,]
relation.habitats.sites.g1_g2_g3_b1_b2.750m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==750,]
relation.habitats.sites.g1_g2_g3_b1_b2.1000m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==1000,]
relation.habitats.sites.g1_g2_g3_b1_b2.3000m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==3000,]
relation.habitats.sites.g1_g2_g3_b1_b2.5000m<-relation.habitats.sites.g1_g2_g3_b1_b2[relation.habitats.sites.g1_g2_g3_b1_b2$distance==5000,]

#preparation des donnée d'entrée

 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2.30m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2.30m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2.30m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_30m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2.60m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2.60m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2.60m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_60m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2.100m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2.100m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2.100m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_100m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2.250m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2.250m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2.250m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_250m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2.500m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2.500m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2.500m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_500m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2.750m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2.750m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2.750m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_750m=as.data.frame(unique(m))

 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2.1000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2.1000m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2.1000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_1000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2.3000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2.3000m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2.3000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_3000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2.5000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2.5000m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2.5000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_5000m=as.data.frame(unique(m))


habitat_par_site_g1_g2_g3_b1_b2_30m<-habitat_par_site_g1_g2_g3_b1_b2_30m[habitat_par_site_g1_g2_g3_b1_b2_30m$Code_Site%in%sites.etudies,]
habitat_par_site_g1_g2_g3_b1_b2_60m<-habitat_par_site_g1_g2_g3_b1_b2_60m[habitat_par_site_g1_g2_g3_b1_b2_60m$Code_Site%in%sites.etudies,]
habitat_par_site_g1_g2_g3_b1_b2_100m<-habitat_par_site_g1_g2_g3_b1_b2_100m[habitat_par_site_g1_g2_g3_b1_b2_100m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_b2_250m<-habitat_par_site_g1_g2_g3_b1_b2_250m[habitat_par_site_g1_g2_g3_b1_b2_250m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_b2_500m<-habitat_par_site_g1_g2_g3_b1_b2_500m[habitat_par_site_g1_g2_g3_b1_b2_500m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_b2_750m<-habitat_par_site_g1_g2_g3_b1_b2_750m[habitat_par_site_g1_g2_g3_b1_b2_750m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_b2_1000m<-habitat_par_site_g1_g2_g3_b1_b2_1000m[habitat_par_site_g1_g2_g3_b1_b2_1000m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_b2_3000m<-habitat_par_site_g1_g2_g3_b1_b2_3000m[habitat_par_site_g1_g2_g3_b1_b2_3000m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_b2_5000m<-habitat_par_site_g1_g2_g3_b1_b2_5000m[habitat_par_site_g1_g2_g3_b1_b2_5000m$Code_Site%in%sites.etudies,]


#création de table de présence absence
#habitat géo1

habitat_par_site_g1_g2_g3_b1_b2_30m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_30m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_30m$Code_hab)
habitat_par_site_g1_g2_g3_b1_b2_60m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_60m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_60m$Code_hab)
habitat_par_site_g1_g2_g3_b1_b2_100m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_100m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_100m$Code_hab)
habitat_par_site_g1_g2_g3_b1_b2_250m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_250m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_250m$Code_hab)
habitat_par_site_g1_g2_g3_b1_b2_500m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_500m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_500m$Code_hab)
habitat_par_site_g1_g2_g3_b1_b2_750m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_750m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_750m$Code_hab)
habitat_par_site_g1_g2_g3_b1_b2_1000m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_1000m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_1000m$Code_hab)
habitat_par_site_g1_g2_g3_b1_b2_3000m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_3000m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_3000m$Code_hab)
habitat_par_site_g1_g2_g3_b1_b2_5000m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_5000m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_5000m$Code_hab)


#selection des sites
setwd(repertoire)
setwd(input)
source('script.richesse.complementarite.r')

selection.richesse.complementarite.g1_g2_g3_b1_b2.30m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_30m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1_b2.60m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_60m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1_b2.100m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_100m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1_b2.250m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_250m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1_b2.500m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_500m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1_b2.750m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_750m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1_b2.1000m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_1000m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1_b2.3000m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_3000m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1_b2.5000m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_5000m_pres_abs)

#RICHESSE ET DIVERSITE DE POISSON CONTENUE DANS LES SELECTIONS DE SITES#
source('script.diversite.selectionnee.r')

diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.30m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2.30m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.60m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2.60m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.100m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2.100m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.250m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2.250m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.500m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2.500m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.750m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2.750m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.1000m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2.1000m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.3000m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2.3000m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.5000m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2.5000m,recensement_unique)

#Calcul des SAI
source('script.sai.r')
SAI_richesse.complementarite.g1_g2_g3_b1_b2.30m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.30m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1_b2.60m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.60m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1_b2.100m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.100m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1_b2.250m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.250m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1_b2.500m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.500m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1_b2.750m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.750m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1_b2.1000m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.1000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1_b2.3000m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.3000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1_b2.5000m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2.5000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)



############################HABITAT g1_g2_g3_b1_b2_b3###########################################

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
#habitat géo1
relation.habitats.sites.g1_g2_g3_b1_b2_b3.30m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==30,]
relation.habitats.sites.g1_g2_g3_b1_b2_b3.60m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==60,]
relation.habitats.sites.g1_g2_g3_b1_b2_b3.100m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==100,]
relation.habitats.sites.g1_g2_g3_b1_b2_b3.250m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==250,]
relation.habitats.sites.g1_g2_g3_b1_b2_b3.500m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==500,]
relation.habitats.sites.g1_g2_g3_b1_b2_b3.750m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==750,]
relation.habitats.sites.g1_g2_g3_b1_b2_b3.1000m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==1000,]
relation.habitats.sites.g1_g2_g3_b1_b2_b3.3000m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==3000,]
relation.habitats.sites.g1_g2_g3_b1_b2_b3.5000m<-relation.habitats.sites.g1_g2_g3_b1_b2_b3[relation.habitats.sites.g1_g2_g3_b1_b2_b3$distance==5000,]

#preparation des donnée d'entrée

 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2_b3.30m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2_b3.30m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2_b3.30m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_b3_30m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2_b3.60m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2_b3.60m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2_b3.60m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_b3_60m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2_b3.100m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2_b3.100m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2_b3.100m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_b3_100m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2_b3.250m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2_b3.250m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2_b3.250m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_b3_250m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2_b3.500m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2_b3.500m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2_b3.500m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_b3_500m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2_b3.750m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2_b3.750m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2_b3.750m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_b3_750m=as.data.frame(unique(m))

 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2_b3.1000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2_b3.1000m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2_b3.1000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_b3_1000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2_b3.3000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2_b3.3000m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2_b3.3000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_b3_3000m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.g1_g2_g3_b1_b2_b3.5000m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.g1_g2_g3_b1_b2_b3.5000m$Code_Site)
  m[,2]=relation.habitats.sites.g1_g2_g3_b1_b2_b3.5000m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_g1_g2_g3_b1_b2_b3_5000m=as.data.frame(unique(m))


habitat_par_site_g1_g2_g3_b1_b2_b3_30m<-habitat_par_site_g1_g2_g3_b1_b2_b3_30m[habitat_par_site_g1_g2_g3_b1_b2_b3_30m$Code_Site%in%sites.etudies,]
habitat_par_site_g1_g2_g3_b1_b2_b3_60m<-habitat_par_site_g1_g2_g3_b1_b2_b3_60m[habitat_par_site_g1_g2_g3_b1_b2_b3_60m$Code_Site%in%sites.etudies,]
habitat_par_site_g1_g2_g3_b1_b2_b3_100m<-habitat_par_site_g1_g2_g3_b1_b2_b3_100m[habitat_par_site_g1_g2_g3_b1_b2_b3_100m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_b2_b3_250m<-habitat_par_site_g1_g2_g3_b1_b2_b3_250m[habitat_par_site_g1_g2_g3_b1_b2_b3_250m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_b2_b3_500m<-habitat_par_site_g1_g2_g3_b1_b2_b3_500m[habitat_par_site_g1_g2_g3_b1_b2_b3_500m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_b2_b3_750m<-habitat_par_site_g1_g2_g3_b1_b2_b3_750m[habitat_par_site_g1_g2_g3_b1_b2_b3_750m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_b2_b3_1000m<-habitat_par_site_g1_g2_g3_b1_b2_b3_1000m[habitat_par_site_g1_g2_g3_b1_b2_b3_1000m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_b2_b3_3000m<-habitat_par_site_g1_g2_g3_b1_b2_b3_3000m[habitat_par_site_g1_g2_g3_b1_b2_b3_3000m$Code_Site%in%sites.etudies,]
 habitat_par_site_g1_g2_g3_b1_b2_b3_5000m<-habitat_par_site_g1_g2_g3_b1_b2_b3_5000m[habitat_par_site_g1_g2_g3_b1_b2_b3_5000m$Code_Site%in%sites.etudies,]


#création de table de présence absence
#habitat géo1

habitat_par_site_g1_g2_g3_b1_b2_b3_30m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_b3_30m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_b3_30m$Code_hab)
habitat_par_site_g1_g2_g3_b1_b2_b3_60m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_b3_60m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_b3_60m$Code_hab)
habitat_par_site_g1_g2_g3_b1_b2_b3_100m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_b3_100m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_b3_100m$Code_hab)
habitat_par_site_g1_g2_g3_b1_b2_b3_250m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_b3_250m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_b3_250m$Code_hab)
habitat_par_site_g1_g2_g3_b1_b2_b3_500m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_b3_500m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_b3_500m$Code_hab)
habitat_par_site_g1_g2_g3_b1_b2_b3_750m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_b3_750m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_b3_750m$Code_hab)
habitat_par_site_g1_g2_g3_b1_b2_b3_1000m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_b3_1000m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_b3_1000m$Code_hab)
habitat_par_site_g1_g2_g3_b1_b2_b3_3000m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_b3_3000m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_b3_3000m$Code_hab)
habitat_par_site_g1_g2_g3_b1_b2_b3_5000m_pres_abs<-table(habitat_par_site_g1_g2_g3_b1_b2_b3_5000m$Code_Site,habitat_par_site_g1_g2_g3_b1_b2_b3_5000m$Code_hab)


#selection des sites
setwd(repertoire)
setwd(input)
source('script.richesse.complementarite.r')

selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.30m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_b3_30m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.60m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_b3_60m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.100m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_b3_100m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.250m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_b3_250m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.500m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_b3_500m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.750m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_b3_750m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.1000m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_b3_1000m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.3000m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_b3_3000m_pres_abs)
selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.5000m<-run.rich(habitat_par_site_g1_g2_g3_b1_b2_b3_5000m_pres_abs)

#RICHESSE ET DIVERSITE DE POISSON CONTENUE DANS LES SELECTIONS DE SITES#
source('script.diversite.selectionnee.r')

diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.30m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.30m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.60m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.60m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.100m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.100m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.250m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.250m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.500m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.500m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.750m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.750m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.1000m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.1000m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.3000m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.3000m,recensement_unique)
diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.5000m<-species.richness.evol(selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.5000m,recensement_unique)

#Calcul des SAI
source('script.sai.r')
SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.30m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.30m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.60m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.60m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.100m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.100m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.250m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.250m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.500m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.500m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.750m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.750m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.1000m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.1000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.3000m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.3000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.5000m<-SAI(diversite.selection.richesse.complementarite.g1_g2_g3_b1_b2_b3.5000m,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)


#STOCKAGE DES RESULTATS

targeted.station.results[,p+1]=c(SAI_richesse.complementarite.geo1.30m,SAI_richesse.complementarite.geo1.60m,SAI_richesse.complementarite.geo1.100m,SAI_richesse.complementarite.geo1.250m,
SAI_richesse.complementarite.geo1.500m,SAI_richesse.complementarite.geo1.750m,SAI_richesse.complementarite.geo1.1000m,SAI_richesse.complementarite.geo1.3000m,SAI_richesse.complementarite.geo1.5000m,
SAI_richesse.complementarite.geo2.30m,SAI_richesse.complementarite.geo2.60m,SAI_richesse.complementarite.geo2.100m,SAI_richesse.complementarite.geo2.250m,
SAI_richesse.complementarite.geo2.500m,SAI_richesse.complementarite.geo2.750m,SAI_richesse.complementarite.geo2.1000m,SAI_richesse.complementarite.geo2.3000m,SAI_richesse.complementarite.geo2.5000m,
SAI_richesse.complementarite.geo3.100m,SAI_richesse.complementarite.geo3.250m,
SAI_richesse.complementarite.geo3.500m,SAI_richesse.complementarite.geo3.750m,SAI_richesse.complementarite.geo3.1000m,SAI_richesse.complementarite.geo3.3000m,SAI_richesse.complementarite.geo3.5000m,
SAI_richesse.complementarite.bent1.30m,SAI_richesse.complementarite.bent1.60m,SAI_richesse.complementarite.bent1.100m,SAI_richesse.complementarite.bent1.250m,
SAI_richesse.complementarite.bent1.500m,SAI_richesse.complementarite.bent1.750m,SAI_richesse.complementarite.bent1.1000m,SAI_richesse.complementarite.bent1.3000m,SAI_richesse.complementarite.bent1.5000m,
SAI_richesse.complementarite.bent2.30m,SAI_richesse.complementarite.bent2.60m,SAI_richesse.complementarite.bent2.100m,SAI_richesse.complementarite.bent2.250m,
SAI_richesse.complementarite.bent2.500m,SAI_richesse.complementarite.bent2.750m,SAI_richesse.complementarite.bent2.1000m,SAI_richesse.complementarite.bent2.3000m,SAI_richesse.complementarite.bent2.5000m,
SAI_richesse.complementarite.bent3.30m,SAI_richesse.complementarite.bent3.60m,SAI_richesse.complementarite.bent3.100m,SAI_richesse.complementarite.bent3.250m,
SAI_richesse.complementarite.bent3.500m,SAI_richesse.complementarite.bent3.750m,SAI_richesse.complementarite.bent3.1000m,SAI_richesse.complementarite.bent3.3000m,SAI_richesse.complementarite.bent3.5000m,
SAI_richesse.complementarite.geo1_geo2.30m,SAI_richesse.complementarite.geo1_geo2.60m,SAI_richesse.complementarite.geo1_geo2.100m,SAI_richesse.complementarite.geo1_geo2.250m,
SAI_richesse.complementarite.geo1_geo2.500m,SAI_richesse.complementarite.geo1_geo2.750m,SAI_richesse.complementarite.geo1_geo2.1000m,SAI_richesse.complementarite.geo1_geo2.3000m,SAI_richesse.complementarite.geo1_geo2.5000m,
SAI_richesse.complementarite.geo1_geo2_geo3.30m,SAI_richesse.complementarite.geo1_geo2_geo3.60m,SAI_richesse.complementarite.geo1_geo2_geo3.100m,SAI_richesse.complementarite.geo1_geo2_geo3.250m,
SAI_richesse.complementarite.geo1_geo2_geo3.500m,SAI_richesse.complementarite.geo1_geo2_geo3.750m,SAI_richesse.complementarite.geo1_geo2_geo3.1000m,SAI_richesse.complementarite.geo1_geo2_geo3.3000m,SAI_richesse.complementarite.geo1_geo2_geo3.5000m,
SAI_richesse.complementarite.g1_g2_g3_b1.30m,SAI_richesse.complementarite.g1_g2_g3_b1.60m,SAI_richesse.complementarite.g1_g2_g3_b1.100m,SAI_richesse.complementarite.g1_g2_g3_b1.250m,
SAI_richesse.complementarite.g1_g2_g3_b1.500m,SAI_richesse.complementarite.g1_g2_g3_b1.750m,SAI_richesse.complementarite.g1_g2_g3_b1.1000m,SAI_richesse.complementarite.g1_g2_g3_b1.3000m,SAI_richesse.complementarite.g1_g2_g3_b1.5000m,
SAI_richesse.complementarite.g1_g2_g3_b1_b2.30m,SAI_richesse.complementarite.g1_g2_g3_b1_b2.60m,SAI_richesse.complementarite.g1_g2_g3_b1_b2.100m,SAI_richesse.complementarite.g1_g2_g3_b1_b2.250m,
SAI_richesse.complementarite.g1_g2_g3_b1_b2.500m,SAI_richesse.complementarite.g1_g2_g3_b1_b2.750m,SAI_richesse.complementarite.g1_g2_g3_b1_b2.1000m,SAI_richesse.complementarite.g1_g2_g3_b1_b2.3000m,SAI_richesse.complementarite.g1_g2_g3_b1_b2.5000m,
SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.30m,SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.60m,SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.100m,SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.250m,
SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.500m,SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.750m,SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.1000m,SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.3000m,SAI_richesse.complementarite.g1_g2_g3_b1_b2_b3.5000m)
}

 targeted.station.results
