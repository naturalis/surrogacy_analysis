targeted.station.results<-matrix(0,97,28)  #999

 #LE tableau de résultats à la forme :

#  ALL stations       - station TS249        -station TS250   -station TS251 ........
#    SAIg130m
#    SAI g160m
#    SAI g1100m
#              etc....
#
#


 #######################        ALL STATIONS     #############################

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


#calcul indice de rareté des espèces

  b=unique(recensement_unique$Code_poisson)
 
  Indice_rarete_sp=NULL
  site_poisson<-matrix(0,27,length(recensement_unique$Code_poisson))

  for (n in 1:length(b)){
  site_poisson[,n]<-c(as.vector(recensement_unique$Code_Site[recensement_unique$Code_poisson==b[n]]),rep(0,27-length(recensement_unique$Code_Site[recensement_unique$Code_poisson==b[n]])))
  names(site_poisson[,n])=b[n]
  Indice_rarete_sp[n]<-2^(27-length(site_poisson[,n][site_poisson[,n]!=0]))
  names(Indice_rarete_sp)[n]=as.vector(b[n])     #Indice_rarete_sp est un vecteur qui contien l'indice de rareté pour chaque espèce
  } 
    
  #calcul indice de rareté des stations
  
  sites.etudies=c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')

  
    poissons_dans_station<-matrix(0,335,27)
    Indice_rarete_station<-NULL
    recensement_unique$Code_Site
   
   for (m in 1:27){             
    poissons_dans_station[,m]<-c(as.vector(recensement_unique$Code_poisson[recensement_unique$Code_Site==sites.etudies[m]]),rep(0,335-length(recensement_unique$Code_poisson[recensement_unique$Code_Site==sites.etudies[m]])))
  Indice_rarete_station[m]=sum(Indice_rarete_sp[names(Indice_rarete_sp)%in%poissons_dans_station[,m]])            
  }              
                           
 
  names(Indice_rarete_station)<-sites.etudies

   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot=sum(Indice_rarete_station)

#Selection des sites selon la rareté de poissons
source('script.rarete.complementarite.r')
selection.rarete.complementarite.poisson<-run.bin(tablo_pres_abs,Indice_rarete_station)


#rareté de poisson incluse dans les sites selectionnés (aléatoirement et selon diversité de poisson)
source('script.diversite.rarete.selectionnee.r')
 diversite.rarete.selection.aleatoire<-matrix(0,dim(tablo_pres_abs)[1],999)  #999
 diversite.rarete.selection.aleatoire<-apply(as.matrix(selection.aleatoire),2,function(x) species.rarete.evol(x,Indice_rarete_station))
  diversite.rarete.moyenne.selection.aleatoire<-apply(diversite.rarete.selection.aleatoire,1,mean)


diversite.selection.rarete.complementarite.poisson<-species.rarete.evol(selection.rarete.complementarite.poisson,Indice_rarete_station)


# Calcul des courbes aléatoires (moyenne, quantile supérieur à 97.5, quantile inférieur à 2.5% sur les 1000 tirages)

diversite.rarete.moyenne.selection.aleatoire
diversite.rarete.borne.sup.selection.aleatoire=apply(diversite.rarete.selection.aleatoire,1,function(x) sort(x)[ceiling(0.975*length(x))])
diversite.rarete.borne.inf.selection.aleatoire=apply(diversite.rarete.selection.aleatoire,1,function(x) sort(x)[ceiling(0.025*length(x))])


#nombre total d'espèces de poissons observés
nombre.poissons.total <- nlevels(factor(recensement_unique$Code_poisson))

# fonction pourcentage
percentage=function(total,truc){
  truc * 100 / total
  }

 # Pourcentages, RANDOM, upper, low
diversite.rarete.moyenne.selection.aleatoire.pourcentage <- c(0,percentage(tot,diversite.rarete.moyenne.selection.aleatoire))
diversite.rarete.borne.sup.selection.aleatoire.pourcentage <- c(0,percentage(tot,diversite.rarete.borne.sup.selection.aleatoire))
diversite.rarete.borne.inf.selection.aleatoire.pourcentage <- c(0,percentage(tot,diversite.rarete.borne.inf.selection.aleatoire))

 diversite.selection.rarete.complementarite.poisson.pourcentage<- c(0,percentage(tot,diversite.selection.rarete.complementarite.poisson))




 
############################HABITAT geo1###########################################


 ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances


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



#Load de la requette "REQ_site_geo1_30m.csv"

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

#calcul indice de rareté des habitats

#30m

  b=unique(habitat_par_site_geo1_30m$Code_hab)

    Indice_rarete_geo1_30m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo1_30m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_30m$Code_Site[habitat_par_site_geo1_30m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo1_30m$Code_Site[habitat_par_site_geo1_30m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_30m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_30m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_30m_station<-matrix(0,500,27)
Indice_rarete_geo1_30m_station<-NULL

    for (m in 1:27){
  geo1_30m_station[,m]<-c(as.vector(habitat_par_site_geo1_30m$Code_hab[habitat_par_site_geo1_30m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_30m$Code_hab[habitat_par_site_geo1_30m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_30m_station[m]=sum(Indice_rarete_geo1_30m[names(Indice_rarete_geo1_30m)%in%geo1_30m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_30m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_30m=sum(Indice_rarete_geo1_30m_station)                #somme des indices de rareté  de ttes les sp


   
#60m

  b=unique(habitat_par_site_geo1_60m$Code_hab)

    Indice_rarete_geo1_60m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo1_60m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_60m$Code_Site[habitat_par_site_geo1_60m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo1_60m$Code_Site[habitat_par_site_geo1_60m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_60m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_60m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_60m_station<-matrix(0,500,27)
Indice_rarete_geo1_60m_station<-NULL

    for (m in 1:27){
  geo1_60m_station[,m]<-c(as.vector(habitat_par_site_geo1_60m$Code_hab[habitat_par_site_geo1_60m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_60m$Code_hab[habitat_par_site_geo1_60m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_60m_station[m]=sum(Indice_rarete_geo1_60m[names(Indice_rarete_geo1_60m)%in%geo1_60m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_60m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_60m=sum(Indice_rarete_geo1_60m_station)                #somme des indices de rareté  de ttes les sp


    
#100m

  b=unique(habitat_par_site_geo1_100m$Code_hab)

    Indice_rarete_geo1_100m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo1_100m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_100m$Code_Site[habitat_par_site_geo1_100m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo1_100m$Code_Site[habitat_par_site_geo1_100m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_100m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_100m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_100m_station<-matrix(0,500,27)
Indice_rarete_geo1_100m_station<-NULL

    for (m in 1:27){
  geo1_100m_station[,m]<-c(as.vector(habitat_par_site_geo1_100m$Code_hab[habitat_par_site_geo1_100m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_100m$Code_hab[habitat_par_site_geo1_100m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_100m_station[m]=sum(Indice_rarete_geo1_100m[names(Indice_rarete_geo1_100m)%in%geo1_100m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_100m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_100m=sum(Indice_rarete_geo1_100m_station)                #somme des indices de rareté  de ttes les sp

              
#250m

  b=unique(habitat_par_site_geo1_250m$Code_hab)

    Indice_rarete_geo1_250m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo1_250m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_250m$Code_Site[habitat_par_site_geo1_250m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo1_250m$Code_Site[habitat_par_site_geo1_250m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_250m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_250m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_250m_station<-matrix(0,500,27)
Indice_rarete_geo1_250m_station<-NULL

    for (m in 1:27){
  geo1_250m_station[,m]<-c(as.vector(habitat_par_site_geo1_250m$Code_hab[habitat_par_site_geo1_250m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_250m$Code_hab[habitat_par_site_geo1_250m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_250m_station[m]=sum(Indice_rarete_geo1_250m[names(Indice_rarete_geo1_250m)%in%geo1_250m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_250m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_250m=sum(Indice_rarete_geo1_250m_station)                #somme des indices de rareté  de ttes les sp

     
#500m

  b=unique(habitat_par_site_geo1_500m$Code_hab)

    Indice_rarete_geo1_500m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo1_500m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_500m$Code_Site[habitat_par_site_geo1_500m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo1_500m$Code_Site[habitat_par_site_geo1_500m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_500m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_500m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_500m_station<-matrix(0,500,27)
Indice_rarete_geo1_500m_station<-NULL

    for (m in 1:27){
  geo1_500m_station[,m]<-c(as.vector(habitat_par_site_geo1_500m$Code_hab[habitat_par_site_geo1_500m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_500m$Code_hab[habitat_par_site_geo1_500m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_500m_station[m]=sum(Indice_rarete_geo1_500m[names(Indice_rarete_geo1_500m)%in%geo1_500m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_500m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_500m=sum(Indice_rarete_geo1_500m_station)                #somme des indices de rareté  de ttes les sp

     
#750m

  b=unique(habitat_par_site_geo1_750m$Code_hab)

    Indice_rarete_geo1_750m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo1_750m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_750m$Code_Site[habitat_par_site_geo1_750m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo1_750m$Code_Site[habitat_par_site_geo1_750m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_750m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_750m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_750m_station<-matrix(0,500,27)
Indice_rarete_geo1_750m_station<-NULL

    for (m in 1:27){
  geo1_750m_station[,m]<-c(as.vector(habitat_par_site_geo1_750m$Code_hab[habitat_par_site_geo1_750m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_750m$Code_hab[habitat_par_site_geo1_750m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_750m_station[m]=sum(Indice_rarete_geo1_750m[names(Indice_rarete_geo1_750m)%in%geo1_750m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_750m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_750m=sum(Indice_rarete_geo1_750m_station)                #somme des indices de rareté  de ttes les sp

      
#1000m

  b=unique(habitat_par_site_geo1_1000m$Code_hab)

    Indice_rarete_geo1_1000m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo1_1000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_1000m$Code_Site[habitat_par_site_geo1_1000m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo1_1000m$Code_Site[habitat_par_site_geo1_1000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_1000m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_1000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_1000m_station<-matrix(0,500,27)
Indice_rarete_geo1_1000m_station<-NULL

    for (m in 1:27){
  geo1_1000m_station[,m]<-c(as.vector(habitat_par_site_geo1_1000m$Code_hab[habitat_par_site_geo1_1000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_1000m$Code_hab[habitat_par_site_geo1_1000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_1000m_station[m]=sum(Indice_rarete_geo1_1000m[names(Indice_rarete_geo1_1000m)%in%geo1_1000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_1000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_1000m=sum(Indice_rarete_geo1_1000m_station)                #somme des indices de rareté  de ttes les sp

                  
#3000m

  b=unique(habitat_par_site_geo1_3000m$Code_hab)

    Indice_rarete_geo1_3000m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo1_3000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_3000m$Code_Site[habitat_par_site_geo1_3000m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo1_3000m$Code_Site[habitat_par_site_geo1_3000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_3000m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_3000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_3000m_station<-matrix(0,500,27)
Indice_rarete_geo1_3000m_station<-NULL

    for (m in 1:27){
  geo1_3000m_station[,m]<-c(as.vector(habitat_par_site_geo1_3000m$Code_hab[habitat_par_site_geo1_3000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_3000m$Code_hab[habitat_par_site_geo1_3000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_3000m_station[m]=sum(Indice_rarete_geo1_3000m[names(Indice_rarete_geo1_3000m)%in%geo1_3000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_3000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_3000m=sum(Indice_rarete_geo1_3000m_station)                #somme des indices de rareté  de ttes les sp

         
#5000m

  b=unique(habitat_par_site_geo1_5000m$Code_hab)

    Indice_rarete_geo1_5000m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo1_5000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_5000m$Code_Site[habitat_par_site_geo1_5000m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo1_5000m$Code_Site[habitat_par_site_geo1_5000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_5000m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_5000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_5000m_station<-matrix(0,500,27)
Indice_rarete_geo1_5000m_station<-NULL

    for (m in 1:27){
  geo1_5000m_station[,m]<-c(as.vector(habitat_par_site_geo1_5000m$Code_hab[habitat_par_site_geo1_5000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_5000m$Code_hab[habitat_par_site_geo1_5000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_5000m_station[m]=sum(Indice_rarete_geo1_5000m[names(Indice_rarete_geo1_5000m)%in%geo1_5000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_5000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_5000m=sum(Indice_rarete_geo1_5000m_station)                #somme des indices de rareté  de ttes les sp


source('script.rarete.complementarite.r')

#Selection des dites selon un scenario de rareté complementarité basé sur les habitats

#habitat géo1
selection.rarete.complementarite.geo1.30m<-run.bin(habitat_par_site_geo1_30m_pres_abs,Indice_rarete_geo1_30m_station)
selection.rarete.complementarite.geo1.60m<-run.bin(habitat_par_site_geo1_60m_pres_abs,Indice_rarete_geo1_60m_station)
selection.rarete.complementarite.geo1.100m<-run.bin(habitat_par_site_geo1_100m_pres_abs,Indice_rarete_geo1_100m_station)
selection.rarete.complementarite.geo1.250m<-run.bin(habitat_par_site_geo1_250m_pres_abs,Indice_rarete_geo1_250m_station)
selection.rarete.complementarite.geo1.500m<-run.bin(habitat_par_site_geo1_500m_pres_abs,Indice_rarete_geo1_500m_station)
selection.rarete.complementarite.geo1.750m<-run.bin(habitat_par_site_geo1_750m_pres_abs,Indice_rarete_geo1_750m_station)
selection.rarete.complementarite.geo1.1000m<-run.bin(habitat_par_site_geo1_1000m_pres_abs,Indice_rarete_geo1_1000m_station)
selection.rarete.complementarite.geo1.3000m<-run.bin(habitat_par_site_geo1_3000m_pres_abs,Indice_rarete_geo1_3000m_station)
selection.rarete.complementarite.geo1.5000m<-run.bin(habitat_par_site_geo1_5000m_pres_abs,Indice_rarete_geo1_5000m_station)

#DIVERSITE SELECTIONNE PAR LE SCRIPT RARETE COMPLEMENTARITE#
 source('script.diversite.rarete.selectionnee.r')
 
 diversite.selection.rarete.complementarite.geo1.30m<-species.rarete.evol(selection.rarete.complementarite.geo1.30m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1.60m<-species.rarete.evol(selection.rarete.complementarite.geo1.60m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1.100m<-species.rarete.evol(selection.rarete.complementarite.geo1.100m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1.250m<-species.rarete.evol(selection.rarete.complementarite.geo1.250m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1.500m<-species.rarete.evol(selection.rarete.complementarite.geo1.500m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1.750m<-species.rarete.evol(selection.rarete.complementarite.geo1.750m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1.1000m<-species.rarete.evol(selection.rarete.complementarite.geo1.1000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1.3000m<-species.rarete.evol(selection.rarete.complementarite.geo1.3000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1.5000m<-species.rarete.evol(selection.rarete.complementarite.geo1.5000m,Indice_rarete_station)

#CALCUL DES SAI


source('script.sai.r')
SAI_rarete.complementarite.geo1.30m<-SAI(diversite.selection.rarete.complementarite.geo1.30m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1.60m<-SAI(diversite.selection.rarete.complementarite.geo1.60m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1.100m<-SAI(diversite.selection.rarete.complementarite.geo1.100m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1.250m<-SAI(diversite.selection.rarete.complementarite.geo1.250m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1.500m<-SAI(diversite.selection.rarete.complementarite.geo1.500m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1.750m<-SAI(diversite.selection.rarete.complementarite.geo1.750m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1.1000m<-SAI(diversite.selection.rarete.complementarite.geo1.1000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1.3000m<-SAI(diversite.selection.rarete.complementarite.geo1.3000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1.5000m<-SAI(diversite.selection.rarete.complementarite.geo1.5000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)


  
 
############################HABITAT geo2###########################################


 ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances


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



#Load de la requette "REQ_site_geo2_30m.csv"

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

#calcul indice de rareté des habitats

#30m

  b=unique(habitat_par_site_geo2_30m$Code_hab)

    Indice_rarete_geo2_30m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo2_30m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo2_30m$Code_Site[habitat_par_site_geo2_30m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo2_30m$Code_Site[habitat_par_site_geo2_30m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo2_30m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo2_30m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo2_30m_station<-matrix(0,500,27)
Indice_rarete_geo2_30m_station<-NULL

    for (m in 1:27){
  geo2_30m_station[,m]<-c(as.vector(habitat_par_site_geo2_30m$Code_hab[habitat_par_site_geo2_30m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo2_30m$Code_hab[habitat_par_site_geo2_30m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo2_30m_station[m]=sum(Indice_rarete_geo2_30m[names(Indice_rarete_geo2_30m)%in%geo2_30m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo2_30m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo2_30m=sum(Indice_rarete_geo2_30m_station)                #somme des indices de rareté  de ttes les sp


   
#60m

  b=unique(habitat_par_site_geo2_60m$Code_hab)

    Indice_rarete_geo2_60m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo2_60m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo2_60m$Code_Site[habitat_par_site_geo2_60m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo2_60m$Code_Site[habitat_par_site_geo2_60m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo2_60m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo2_60m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo2_60m_station<-matrix(0,500,27)
Indice_rarete_geo2_60m_station<-NULL

    for (m in 1:27){
  geo2_60m_station[,m]<-c(as.vector(habitat_par_site_geo2_60m$Code_hab[habitat_par_site_geo2_60m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo2_60m$Code_hab[habitat_par_site_geo2_60m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo2_60m_station[m]=sum(Indice_rarete_geo2_60m[names(Indice_rarete_geo2_60m)%in%geo2_60m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo2_60m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo2_60m=sum(Indice_rarete_geo2_60m_station)                #somme des indices de rareté  de ttes les sp


    
#100m

  b=unique(habitat_par_site_geo2_100m$Code_hab)

    Indice_rarete_geo2_100m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo2_100m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo2_100m$Code_Site[habitat_par_site_geo2_100m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo2_100m$Code_Site[habitat_par_site_geo2_100m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo2_100m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo2_100m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo2_100m_station<-matrix(0,500,27)
Indice_rarete_geo2_100m_station<-NULL

    for (m in 1:27){
  geo2_100m_station[,m]<-c(as.vector(habitat_par_site_geo2_100m$Code_hab[habitat_par_site_geo2_100m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo2_100m$Code_hab[habitat_par_site_geo2_100m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo2_100m_station[m]=sum(Indice_rarete_geo2_100m[names(Indice_rarete_geo2_100m)%in%geo2_100m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo2_100m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo2_100m=sum(Indice_rarete_geo2_100m_station)                #somme des indices de rareté  de ttes les sp

              
#250m

  b=unique(habitat_par_site_geo2_250m$Code_hab)

    Indice_rarete_geo2_250m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo2_250m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo2_250m$Code_Site[habitat_par_site_geo2_250m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo2_250m$Code_Site[habitat_par_site_geo2_250m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo2_250m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo2_250m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo2_250m_station<-matrix(0,500,27)
Indice_rarete_geo2_250m_station<-NULL

    for (m in 1:27){
  geo2_250m_station[,m]<-c(as.vector(habitat_par_site_geo2_250m$Code_hab[habitat_par_site_geo2_250m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo2_250m$Code_hab[habitat_par_site_geo2_250m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo2_250m_station[m]=sum(Indice_rarete_geo2_250m[names(Indice_rarete_geo2_250m)%in%geo2_250m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo2_250m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo2_250m=sum(Indice_rarete_geo2_250m_station)                #somme des indices de rareté  de ttes les sp

     
#500m

  b=unique(habitat_par_site_geo2_500m$Code_hab)

    Indice_rarete_geo2_500m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo2_500m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo2_500m$Code_Site[habitat_par_site_geo2_500m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo2_500m$Code_Site[habitat_par_site_geo2_500m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo2_500m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo2_500m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo2_500m_station<-matrix(0,500,27)
Indice_rarete_geo2_500m_station<-NULL

    for (m in 1:27){
  geo2_500m_station[,m]<-c(as.vector(habitat_par_site_geo2_500m$Code_hab[habitat_par_site_geo2_500m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo2_500m$Code_hab[habitat_par_site_geo2_500m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo2_500m_station[m]=sum(Indice_rarete_geo2_500m[names(Indice_rarete_geo2_500m)%in%geo2_500m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo2_500m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo2_500m=sum(Indice_rarete_geo2_500m_station)                #somme des indices de rareté  de ttes les sp

     
#750m

  b=unique(habitat_par_site_geo2_750m$Code_hab)

    Indice_rarete_geo2_750m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo2_750m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo2_750m$Code_Site[habitat_par_site_geo2_750m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo2_750m$Code_Site[habitat_par_site_geo2_750m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo2_750m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo2_750m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo2_750m_station<-matrix(0,500,27)
Indice_rarete_geo2_750m_station<-NULL

    for (m in 1:27){
  geo2_750m_station[,m]<-c(as.vector(habitat_par_site_geo2_750m$Code_hab[habitat_par_site_geo2_750m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo2_750m$Code_hab[habitat_par_site_geo2_750m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo2_750m_station[m]=sum(Indice_rarete_geo2_750m[names(Indice_rarete_geo2_750m)%in%geo2_750m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo2_750m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo2_750m=sum(Indice_rarete_geo2_750m_station)                #somme des indices de rareté  de ttes les sp

      
#1000m

  b=unique(habitat_par_site_geo2_1000m$Code_hab)

    Indice_rarete_geo2_1000m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo2_1000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo2_1000m$Code_Site[habitat_par_site_geo2_1000m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo2_1000m$Code_Site[habitat_par_site_geo2_1000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo2_1000m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo2_1000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo2_1000m_station<-matrix(0,500,27)
Indice_rarete_geo2_1000m_station<-NULL

    for (m in 1:27){
  geo2_1000m_station[,m]<-c(as.vector(habitat_par_site_geo2_1000m$Code_hab[habitat_par_site_geo2_1000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo2_1000m$Code_hab[habitat_par_site_geo2_1000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo2_1000m_station[m]=sum(Indice_rarete_geo2_1000m[names(Indice_rarete_geo2_1000m)%in%geo2_1000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo2_1000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo2_1000m=sum(Indice_rarete_geo2_1000m_station)                #somme des indices de rareté  de ttes les sp

                  
#3000m

  b=unique(habitat_par_site_geo2_3000m$Code_hab)

    Indice_rarete_geo2_3000m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo2_3000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo2_3000m$Code_Site[habitat_par_site_geo2_3000m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo2_3000m$Code_Site[habitat_par_site_geo2_3000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo2_3000m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo2_3000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo2_3000m_station<-matrix(0,500,27)
Indice_rarete_geo2_3000m_station<-NULL

    for (m in 1:27){
  geo2_3000m_station[,m]<-c(as.vector(habitat_par_site_geo2_3000m$Code_hab[habitat_par_site_geo2_3000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo2_3000m$Code_hab[habitat_par_site_geo2_3000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo2_3000m_station[m]=sum(Indice_rarete_geo2_3000m[names(Indice_rarete_geo2_3000m)%in%geo2_3000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo2_3000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo2_3000m=sum(Indice_rarete_geo2_3000m_station)                #somme des indices de rareté  de ttes les sp

         
#5000m

  b=unique(habitat_par_site_geo2_5000m$Code_hab)

    Indice_rarete_geo2_5000m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo2_5000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo2_5000m$Code_Site[habitat_par_site_geo2_5000m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo2_5000m$Code_Site[habitat_par_site_geo2_5000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo2_5000m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo2_5000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo2_5000m_station<-matrix(0,500,27)
Indice_rarete_geo2_5000m_station<-NULL

    for (m in 1:27){
  geo2_5000m_station[,m]<-c(as.vector(habitat_par_site_geo2_5000m$Code_hab[habitat_par_site_geo2_5000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo2_5000m$Code_hab[habitat_par_site_geo2_5000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo2_5000m_station[m]=sum(Indice_rarete_geo2_5000m[names(Indice_rarete_geo2_5000m)%in%geo2_5000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo2_5000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo2_5000m=sum(Indice_rarete_geo2_5000m_station)                #somme des indices de rareté  de ttes les sp


source('script.rarete.complementarite.r')

#Selection des dites selon un scenario de rareté complementarité basé sur les habitats

#habitat géo1
selection.rarete.complementarite.geo2.30m<-run.bin(habitat_par_site_geo2_30m_pres_abs,Indice_rarete_geo2_30m_station)
selection.rarete.complementarite.geo2.60m<-run.bin(habitat_par_site_geo2_60m_pres_abs,Indice_rarete_geo2_60m_station)
selection.rarete.complementarite.geo2.100m<-run.bin(habitat_par_site_geo2_100m_pres_abs,Indice_rarete_geo2_100m_station)
selection.rarete.complementarite.geo2.250m<-run.bin(habitat_par_site_geo2_250m_pres_abs,Indice_rarete_geo2_250m_station)
selection.rarete.complementarite.geo2.500m<-run.bin(habitat_par_site_geo2_500m_pres_abs,Indice_rarete_geo2_500m_station)
selection.rarete.complementarite.geo2.750m<-run.bin(habitat_par_site_geo2_750m_pres_abs,Indice_rarete_geo2_750m_station)
selection.rarete.complementarite.geo2.1000m<-run.bin(habitat_par_site_geo2_1000m_pres_abs,Indice_rarete_geo2_1000m_station)
selection.rarete.complementarite.geo2.3000m<-run.bin(habitat_par_site_geo2_3000m_pres_abs,Indice_rarete_geo2_3000m_station)
selection.rarete.complementarite.geo2.5000m<-run.bin(habitat_par_site_geo2_5000m_pres_abs,Indice_rarete_geo2_5000m_station)

#DIVERSITE SELECTIONNE PAR LE SCRIPT RARETE COMPLEMENTARITE#
 source('script.diversite.rarete.selectionnee.r')
 
 diversite.selection.rarete.complementarite.geo2.30m<-species.rarete.evol(selection.rarete.complementarite.geo2.30m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo2.60m<-species.rarete.evol(selection.rarete.complementarite.geo2.60m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo2.100m<-species.rarete.evol(selection.rarete.complementarite.geo2.100m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo2.250m<-species.rarete.evol(selection.rarete.complementarite.geo2.250m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo2.500m<-species.rarete.evol(selection.rarete.complementarite.geo2.500m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo2.750m<-species.rarete.evol(selection.rarete.complementarite.geo2.750m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo2.1000m<-species.rarete.evol(selection.rarete.complementarite.geo2.1000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo2.3000m<-species.rarete.evol(selection.rarete.complementarite.geo2.3000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo2.5000m<-species.rarete.evol(selection.rarete.complementarite.geo2.5000m,Indice_rarete_station)

#CALCUL DES SAI


source('script.sai.r')
SAI_rarete.complementarite.geo2.30m<-SAI(diversite.selection.rarete.complementarite.geo2.30m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo2.60m<-SAI(diversite.selection.rarete.complementarite.geo2.60m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo2.100m<-SAI(diversite.selection.rarete.complementarite.geo2.100m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo2.250m<-SAI(diversite.selection.rarete.complementarite.geo2.250m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo2.500m<-SAI(diversite.selection.rarete.complementarite.geo2.500m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo2.750m<-SAI(diversite.selection.rarete.complementarite.geo2.750m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo2.1000m<-SAI(diversite.selection.rarete.complementarite.geo2.1000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo2.3000m<-SAI(diversite.selection.rarete.complementarite.geo2.3000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo2.5000m<-SAI(diversite.selection.rarete.complementarite.geo2.5000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)





 
############################HABITAT geo3###########################################


 ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances


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



#Load de la requette "REQ_site_geo3_30m.csv"

 


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

#calcul indice de rareté des habitats
  
#100m

  b=unique(habitat_par_site_geo3_100m$Code_hab)

    Indice_rarete_geo3_100m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo3_100m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo3_100m$Code_Site[habitat_par_site_geo3_100m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo3_100m$Code_Site[habitat_par_site_geo3_100m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo3_100m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo3_100m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo3_100m_station<-matrix(0,500,27)
Indice_rarete_geo3_100m_station<-NULL

    for (m in 1:27){
  geo3_100m_station[,m]<-c(as.vector(habitat_par_site_geo3_100m$Code_hab[habitat_par_site_geo3_100m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo3_100m$Code_hab[habitat_par_site_geo3_100m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo3_100m_station[m]=sum(Indice_rarete_geo3_100m[names(Indice_rarete_geo3_100m)%in%geo3_100m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo3_100m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo3_100m=sum(Indice_rarete_geo3_100m_station)                #somme des indices de rareté  de ttes les sp

              
#250m

  b=unique(habitat_par_site_geo3_250m$Code_hab)

    Indice_rarete_geo3_250m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo3_250m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo3_250m$Code_Site[habitat_par_site_geo3_250m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo3_250m$Code_Site[habitat_par_site_geo3_250m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo3_250m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo3_250m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo3_250m_station<-matrix(0,500,27)
Indice_rarete_geo3_250m_station<-NULL

    for (m in 1:27){
  geo3_250m_station[,m]<-c(as.vector(habitat_par_site_geo3_250m$Code_hab[habitat_par_site_geo3_250m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo3_250m$Code_hab[habitat_par_site_geo3_250m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo3_250m_station[m]=sum(Indice_rarete_geo3_250m[names(Indice_rarete_geo3_250m)%in%geo3_250m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo3_250m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo3_250m=sum(Indice_rarete_geo3_250m_station)                #somme des indices de rareté  de ttes les sp

     
#500m

  b=unique(habitat_par_site_geo3_500m$Code_hab)

    Indice_rarete_geo3_500m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo3_500m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo3_500m$Code_Site[habitat_par_site_geo3_500m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo3_500m$Code_Site[habitat_par_site_geo3_500m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo3_500m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo3_500m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo3_500m_station<-matrix(0,500,27)
Indice_rarete_geo3_500m_station<-NULL

    for (m in 1:27){
  geo3_500m_station[,m]<-c(as.vector(habitat_par_site_geo3_500m$Code_hab[habitat_par_site_geo3_500m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo3_500m$Code_hab[habitat_par_site_geo3_500m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo3_500m_station[m]=sum(Indice_rarete_geo3_500m[names(Indice_rarete_geo3_500m)%in%geo3_500m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo3_500m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo3_500m=sum(Indice_rarete_geo3_500m_station)                #somme des indices de rareté  de ttes les sp

     
#750m

  b=unique(habitat_par_site_geo3_750m$Code_hab)

    Indice_rarete_geo3_750m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo3_750m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo3_750m$Code_Site[habitat_par_site_geo3_750m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo3_750m$Code_Site[habitat_par_site_geo3_750m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo3_750m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo3_750m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo3_750m_station<-matrix(0,500,27)
Indice_rarete_geo3_750m_station<-NULL

    for (m in 1:27){
  geo3_750m_station[,m]<-c(as.vector(habitat_par_site_geo3_750m$Code_hab[habitat_par_site_geo3_750m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo3_750m$Code_hab[habitat_par_site_geo3_750m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo3_750m_station[m]=sum(Indice_rarete_geo3_750m[names(Indice_rarete_geo3_750m)%in%geo3_750m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo3_750m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo3_750m=sum(Indice_rarete_geo3_750m_station)                #somme des indices de rareté  de ttes les sp

      
#1000m

  b=unique(habitat_par_site_geo3_1000m$Code_hab)

    Indice_rarete_geo3_1000m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo3_1000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo3_1000m$Code_Site[habitat_par_site_geo3_1000m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo3_1000m$Code_Site[habitat_par_site_geo3_1000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo3_1000m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo3_1000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo3_1000m_station<-matrix(0,500,27)
Indice_rarete_geo3_1000m_station<-NULL

    for (m in 1:27){
  geo3_1000m_station[,m]<-c(as.vector(habitat_par_site_geo3_1000m$Code_hab[habitat_par_site_geo3_1000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo3_1000m$Code_hab[habitat_par_site_geo3_1000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo3_1000m_station[m]=sum(Indice_rarete_geo3_1000m[names(Indice_rarete_geo3_1000m)%in%geo3_1000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo3_1000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo3_1000m=sum(Indice_rarete_geo3_1000m_station)                #somme des indices de rareté  de ttes les sp

                  
#3000m

  b=unique(habitat_par_site_geo3_3000m$Code_hab)

    Indice_rarete_geo3_3000m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo3_3000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo3_3000m$Code_Site[habitat_par_site_geo3_3000m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo3_3000m$Code_Site[habitat_par_site_geo3_3000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo3_3000m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo3_3000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo3_3000m_station<-matrix(0,500,27)
Indice_rarete_geo3_3000m_station<-NULL

    for (m in 1:27){
  geo3_3000m_station[,m]<-c(as.vector(habitat_par_site_geo3_3000m$Code_hab[habitat_par_site_geo3_3000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo3_3000m$Code_hab[habitat_par_site_geo3_3000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo3_3000m_station[m]=sum(Indice_rarete_geo3_3000m[names(Indice_rarete_geo3_3000m)%in%geo3_3000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo3_3000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo3_3000m=sum(Indice_rarete_geo3_3000m_station)                #somme des indices de rareté  de ttes les sp

         
#5000m

  b=unique(habitat_par_site_geo3_5000m$Code_hab)

    Indice_rarete_geo3_5000m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo3_5000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo3_5000m$Code_Site[habitat_par_site_geo3_5000m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo3_5000m$Code_Site[habitat_par_site_geo3_5000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo3_5000m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo3_5000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo3_5000m_station<-matrix(0,500,27)
Indice_rarete_geo3_5000m_station<-NULL

    for (m in 1:27){
  geo3_5000m_station[,m]<-c(as.vector(habitat_par_site_geo3_5000m$Code_hab[habitat_par_site_geo3_5000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo3_5000m$Code_hab[habitat_par_site_geo3_5000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo3_5000m_station[m]=sum(Indice_rarete_geo3_5000m[names(Indice_rarete_geo3_5000m)%in%geo3_5000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo3_5000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo3_5000m=sum(Indice_rarete_geo3_5000m_station)                #somme des indices de rareté  de ttes les sp


source('script.rarete.complementarite.r')

#Selection des dites selon un scenario de rareté complementarité basé sur les habitats

#habitat géo1
selection.rarete.complementarite.geo3.100m<-run.bin(habitat_par_site_geo3_100m_pres_abs,Indice_rarete_geo3_100m_station)
selection.rarete.complementarite.geo3.250m<-run.bin(habitat_par_site_geo3_250m_pres_abs,Indice_rarete_geo3_250m_station)
selection.rarete.complementarite.geo3.500m<-run.bin(habitat_par_site_geo3_500m_pres_abs,Indice_rarete_geo3_500m_station)
selection.rarete.complementarite.geo3.750m<-run.bin(habitat_par_site_geo3_750m_pres_abs,Indice_rarete_geo3_750m_station)
selection.rarete.complementarite.geo3.1000m<-run.bin(habitat_par_site_geo3_1000m_pres_abs,Indice_rarete_geo3_1000m_station)
selection.rarete.complementarite.geo3.3000m<-run.bin(habitat_par_site_geo3_3000m_pres_abs,Indice_rarete_geo3_3000m_station)
selection.rarete.complementarite.geo3.5000m<-run.bin(habitat_par_site_geo3_5000m_pres_abs,Indice_rarete_geo3_5000m_station)

#DIVERSITE SELECTIONNE PAR LE SCRIPT RARETE COMPLEMENTARITE#
 source('script.diversite.rarete.selectionnee.r')
 
diversite.selection.rarete.complementarite.geo3.100m<-species.rarete.evol(selection.rarete.complementarite.geo3.100m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo3.250m<-species.rarete.evol(selection.rarete.complementarite.geo3.250m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo3.500m<-species.rarete.evol(selection.rarete.complementarite.geo3.500m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo3.750m<-species.rarete.evol(selection.rarete.complementarite.geo3.750m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo3.1000m<-species.rarete.evol(selection.rarete.complementarite.geo3.1000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo3.3000m<-species.rarete.evol(selection.rarete.complementarite.geo3.3000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo3.5000m<-species.rarete.evol(selection.rarete.complementarite.geo3.5000m,Indice_rarete_station)

#CALCUL DES SAI


source('script.sai.r')
SAI_rarete.complementarite.geo3.100m<-SAI(diversite.selection.rarete.complementarite.geo3.100m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo3.250m<-SAI(diversite.selection.rarete.complementarite.geo3.250m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo3.500m<-SAI(diversite.selection.rarete.complementarite.geo3.500m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo3.750m<-SAI(diversite.selection.rarete.complementarite.geo3.750m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo3.1000m<-SAI(diversite.selection.rarete.complementarite.geo3.1000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo3.3000m<-SAI(diversite.selection.rarete.complementarite.geo3.3000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo3.5000m<-SAI(diversite.selection.rarete.complementarite.geo3.5000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)



 
############################HABITAT bent1###########################################


 ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances


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



#Load de la requette "REQ_site_bent1_30m.csv"

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

#calcul indice de rareté des habitats

#30m

  b=unique(habitat_par_site_bent1_30m$Code_hab)

    Indice_rarete_bent1_30m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_bent1_30m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent1_30m$Code_Site[habitat_par_site_bent1_30m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_bent1_30m$Code_Site[habitat_par_site_bent1_30m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent1_30m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent1_30m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent1_30m_station<-matrix(0,500,27)
Indice_rarete_bent1_30m_station<-NULL

    for (m in 1:27){
  bent1_30m_station[,m]<-c(as.vector(habitat_par_site_bent1_30m$Code_hab[habitat_par_site_bent1_30m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent1_30m$Code_hab[habitat_par_site_bent1_30m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent1_30m_station[m]=sum(Indice_rarete_bent1_30m[names(Indice_rarete_bent1_30m)%in%bent1_30m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent1_30m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent1_30m=sum(Indice_rarete_bent1_30m_station)                #somme des indices de rareté  de ttes les sp


   
#60m

  b=unique(habitat_par_site_bent1_60m$Code_hab)

    Indice_rarete_bent1_60m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_bent1_60m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent1_60m$Code_Site[habitat_par_site_bent1_60m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_bent1_60m$Code_Site[habitat_par_site_bent1_60m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent1_60m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent1_60m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent1_60m_station<-matrix(0,500,27)
Indice_rarete_bent1_60m_station<-NULL

    for (m in 1:27){
  bent1_60m_station[,m]<-c(as.vector(habitat_par_site_bent1_60m$Code_hab[habitat_par_site_bent1_60m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent1_60m$Code_hab[habitat_par_site_bent1_60m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent1_60m_station[m]=sum(Indice_rarete_bent1_60m[names(Indice_rarete_bent1_60m)%in%bent1_60m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent1_60m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent1_60m=sum(Indice_rarete_bent1_60m_station)                #somme des indices de rareté  de ttes les sp


    
#100m

  b=unique(habitat_par_site_bent1_100m$Code_hab)

    Indice_rarete_bent1_100m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_bent1_100m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent1_100m$Code_Site[habitat_par_site_bent1_100m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_bent1_100m$Code_Site[habitat_par_site_bent1_100m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent1_100m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent1_100m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent1_100m_station<-matrix(0,500,27)
Indice_rarete_bent1_100m_station<-NULL

    for (m in 1:27){
  bent1_100m_station[,m]<-c(as.vector(habitat_par_site_bent1_100m$Code_hab[habitat_par_site_bent1_100m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent1_100m$Code_hab[habitat_par_site_bent1_100m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent1_100m_station[m]=sum(Indice_rarete_bent1_100m[names(Indice_rarete_bent1_100m)%in%bent1_100m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent1_100m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent1_100m=sum(Indice_rarete_bent1_100m_station)                #somme des indices de rareté  de ttes les sp

              
#250m

  b=unique(habitat_par_site_bent1_250m$Code_hab)

    Indice_rarete_bent1_250m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_bent1_250m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent1_250m$Code_Site[habitat_par_site_bent1_250m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_bent1_250m$Code_Site[habitat_par_site_bent1_250m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent1_250m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent1_250m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent1_250m_station<-matrix(0,500,27)
Indice_rarete_bent1_250m_station<-NULL

    for (m in 1:27){
  bent1_250m_station[,m]<-c(as.vector(habitat_par_site_bent1_250m$Code_hab[habitat_par_site_bent1_250m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent1_250m$Code_hab[habitat_par_site_bent1_250m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent1_250m_station[m]=sum(Indice_rarete_bent1_250m[names(Indice_rarete_bent1_250m)%in%bent1_250m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent1_250m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent1_250m=sum(Indice_rarete_bent1_250m_station)                #somme des indices de rareté  de ttes les sp

     
#500m

  b=unique(habitat_par_site_bent1_500m$Code_hab)

    Indice_rarete_bent1_500m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_bent1_500m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent1_500m$Code_Site[habitat_par_site_bent1_500m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_bent1_500m$Code_Site[habitat_par_site_bent1_500m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent1_500m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent1_500m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent1_500m_station<-matrix(0,500,27)
Indice_rarete_bent1_500m_station<-NULL

    for (m in 1:27){
  bent1_500m_station[,m]<-c(as.vector(habitat_par_site_bent1_500m$Code_hab[habitat_par_site_bent1_500m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent1_500m$Code_hab[habitat_par_site_bent1_500m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent1_500m_station[m]=sum(Indice_rarete_bent1_500m[names(Indice_rarete_bent1_500m)%in%bent1_500m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent1_500m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent1_500m=sum(Indice_rarete_bent1_500m_station)                #somme des indices de rareté  de ttes les sp

     
#750m

  b=unique(habitat_par_site_bent1_750m$Code_hab)

    Indice_rarete_bent1_750m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_bent1_750m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent1_750m$Code_Site[habitat_par_site_bent1_750m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_bent1_750m$Code_Site[habitat_par_site_bent1_750m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent1_750m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent1_750m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent1_750m_station<-matrix(0,500,27)
Indice_rarete_bent1_750m_station<-NULL

    for (m in 1:27){
  bent1_750m_station[,m]<-c(as.vector(habitat_par_site_bent1_750m$Code_hab[habitat_par_site_bent1_750m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent1_750m$Code_hab[habitat_par_site_bent1_750m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent1_750m_station[m]=sum(Indice_rarete_bent1_750m[names(Indice_rarete_bent1_750m)%in%bent1_750m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent1_750m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent1_750m=sum(Indice_rarete_bent1_750m_station)                #somme des indices de rareté  de ttes les sp

      
#1000m

  b=unique(habitat_par_site_bent1_1000m$Code_hab)

    Indice_rarete_bent1_1000m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_bent1_1000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent1_1000m$Code_Site[habitat_par_site_bent1_1000m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_bent1_1000m$Code_Site[habitat_par_site_bent1_1000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent1_1000m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent1_1000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent1_1000m_station<-matrix(0,500,27)
Indice_rarete_bent1_1000m_station<-NULL

    for (m in 1:27){
  bent1_1000m_station[,m]<-c(as.vector(habitat_par_site_bent1_1000m$Code_hab[habitat_par_site_bent1_1000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent1_1000m$Code_hab[habitat_par_site_bent1_1000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent1_1000m_station[m]=sum(Indice_rarete_bent1_1000m[names(Indice_rarete_bent1_1000m)%in%bent1_1000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent1_1000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent1_1000m=sum(Indice_rarete_bent1_1000m_station)                #somme des indices de rareté  de ttes les sp

                  
#3000m

  b=unique(habitat_par_site_bent1_3000m$Code_hab)

    Indice_rarete_bent1_3000m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_bent1_3000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent1_3000m$Code_Site[habitat_par_site_bent1_3000m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_bent1_3000m$Code_Site[habitat_par_site_bent1_3000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent1_3000m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent1_3000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent1_3000m_station<-matrix(0,500,27)
Indice_rarete_bent1_3000m_station<-NULL

    for (m in 1:27){
  bent1_3000m_station[,m]<-c(as.vector(habitat_par_site_bent1_3000m$Code_hab[habitat_par_site_bent1_3000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent1_3000m$Code_hab[habitat_par_site_bent1_3000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent1_3000m_station[m]=sum(Indice_rarete_bent1_3000m[names(Indice_rarete_bent1_3000m)%in%bent1_3000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent1_3000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent1_3000m=sum(Indice_rarete_bent1_3000m_station)                #somme des indices de rareté  de ttes les sp

         
#5000m

  b=unique(habitat_par_site_bent1_5000m$Code_hab)

    Indice_rarete_bent1_5000m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_bent1_5000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent1_5000m$Code_Site[habitat_par_site_bent1_5000m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_bent1_5000m$Code_Site[habitat_par_site_bent1_5000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent1_5000m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent1_5000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent1_5000m_station<-matrix(0,500,27)
Indice_rarete_bent1_5000m_station<-NULL

    for (m in 1:27){
  bent1_5000m_station[,m]<-c(as.vector(habitat_par_site_bent1_5000m$Code_hab[habitat_par_site_bent1_5000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent1_5000m$Code_hab[habitat_par_site_bent1_5000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent1_5000m_station[m]=sum(Indice_rarete_bent1_5000m[names(Indice_rarete_bent1_5000m)%in%bent1_5000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent1_5000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent1_5000m=sum(Indice_rarete_bent1_5000m_station)                #somme des indices de rareté  de ttes les sp


source('script.rarete.complementarite.r')

#Selection des dites selon un scenario de rareté complementarité basé sur les habitats

#habitat géo1
selection.rarete.complementarite.bent1.30m<-run.bin(habitat_par_site_bent1_30m_pres_abs,Indice_rarete_bent1_30m_station)
selection.rarete.complementarite.bent1.60m<-run.bin(habitat_par_site_bent1_60m_pres_abs,Indice_rarete_bent1_60m_station)
selection.rarete.complementarite.bent1.100m<-run.bin(habitat_par_site_bent1_100m_pres_abs,Indice_rarete_bent1_100m_station)
selection.rarete.complementarite.bent1.250m<-run.bin(habitat_par_site_bent1_250m_pres_abs,Indice_rarete_bent1_250m_station)
selection.rarete.complementarite.bent1.500m<-run.bin(habitat_par_site_bent1_500m_pres_abs,Indice_rarete_bent1_500m_station)
selection.rarete.complementarite.bent1.750m<-run.bin(habitat_par_site_bent1_750m_pres_abs,Indice_rarete_bent1_750m_station)
selection.rarete.complementarite.bent1.1000m<-run.bin(habitat_par_site_bent1_1000m_pres_abs,Indice_rarete_bent1_1000m_station)
selection.rarete.complementarite.bent1.3000m<-run.bin(habitat_par_site_bent1_3000m_pres_abs,Indice_rarete_bent1_3000m_station)
selection.rarete.complementarite.bent1.5000m<-run.bin(habitat_par_site_bent1_5000m_pres_abs,Indice_rarete_bent1_5000m_station)

#DIVERSITE SELECTIONNE PAR LE SCRIPT RARETE COMPLEMENTARITE#
 source('script.diversite.rarete.selectionnee.r')
 
 diversite.selection.rarete.complementarite.bent1.30m<-species.rarete.evol(selection.rarete.complementarite.bent1.30m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent1.60m<-species.rarete.evol(selection.rarete.complementarite.bent1.60m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent1.100m<-species.rarete.evol(selection.rarete.complementarite.bent1.100m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent1.250m<-species.rarete.evol(selection.rarete.complementarite.bent1.250m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent1.500m<-species.rarete.evol(selection.rarete.complementarite.bent1.500m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent1.750m<-species.rarete.evol(selection.rarete.complementarite.bent1.750m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent1.1000m<-species.rarete.evol(selection.rarete.complementarite.bent1.1000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent1.3000m<-species.rarete.evol(selection.rarete.complementarite.bent1.3000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent1.5000m<-species.rarete.evol(selection.rarete.complementarite.bent1.5000m,Indice_rarete_station)

#CALCUL DES SAI


source('script.sai.r')
SAI_rarete.complementarite.bent1.30m<-SAI(diversite.selection.rarete.complementarite.bent1.30m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent1.60m<-SAI(diversite.selection.rarete.complementarite.bent1.60m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent1.100m<-SAI(diversite.selection.rarete.complementarite.bent1.100m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent1.250m<-SAI(diversite.selection.rarete.complementarite.bent1.250m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent1.500m<-SAI(diversite.selection.rarete.complementarite.bent1.500m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent1.750m<-SAI(diversite.selection.rarete.complementarite.bent1.750m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent1.1000m<-SAI(diversite.selection.rarete.complementarite.bent1.1000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent1.3000m<-SAI(diversite.selection.rarete.complementarite.bent1.3000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent1.5000m<-SAI(diversite.selection.rarete.complementarite.bent1.5000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)



 
############################HABITAT bent2###########################################


 ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances


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



#Load de la requette "REQ_site_bent2_30m.csv"

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

#calcul indice de rareté des habitats

#30m

  b=unique(habitat_par_site_bent2_30m$Code_hab)

    Indice_rarete_bent2_30m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_bent2_30m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent2_30m$Code_Site[habitat_par_site_bent2_30m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_bent2_30m$Code_Site[habitat_par_site_bent2_30m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent2_30m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent2_30m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent2_30m_station<-matrix(0,500,27)
Indice_rarete_bent2_30m_station<-NULL

    for (m in 1:27){
  bent2_30m_station[,m]<-c(as.vector(habitat_par_site_bent2_30m$Code_hab[habitat_par_site_bent2_30m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent2_30m$Code_hab[habitat_par_site_bent2_30m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent2_30m_station[m]=sum(Indice_rarete_bent2_30m[names(Indice_rarete_bent2_30m)%in%bent2_30m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent2_30m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent2_30m=sum(Indice_rarete_bent2_30m_station)                #somme des indices de rareté  de ttes les sp


   
#60m

  b=unique(habitat_par_site_bent2_60m$Code_hab)

    Indice_rarete_bent2_60m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_bent2_60m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent2_60m$Code_Site[habitat_par_site_bent2_60m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_bent2_60m$Code_Site[habitat_par_site_bent2_60m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent2_60m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent2_60m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent2_60m_station<-matrix(0,500,27)
Indice_rarete_bent2_60m_station<-NULL

    for (m in 1:27){
  bent2_60m_station[,m]<-c(as.vector(habitat_par_site_bent2_60m$Code_hab[habitat_par_site_bent2_60m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent2_60m$Code_hab[habitat_par_site_bent2_60m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent2_60m_station[m]=sum(Indice_rarete_bent2_60m[names(Indice_rarete_bent2_60m)%in%bent2_60m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent2_60m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent2_60m=sum(Indice_rarete_bent2_60m_station)                #somme des indices de rareté  de ttes les sp


    
#100m

  b=unique(habitat_par_site_bent2_100m$Code_hab)

    Indice_rarete_bent2_100m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_bent2_100m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent2_100m$Code_Site[habitat_par_site_bent2_100m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_bent2_100m$Code_Site[habitat_par_site_bent2_100m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent2_100m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent2_100m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent2_100m_station<-matrix(0,500,27)
Indice_rarete_bent2_100m_station<-NULL

    for (m in 1:27){
  bent2_100m_station[,m]<-c(as.vector(habitat_par_site_bent2_100m$Code_hab[habitat_par_site_bent2_100m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent2_100m$Code_hab[habitat_par_site_bent2_100m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent2_100m_station[m]=sum(Indice_rarete_bent2_100m[names(Indice_rarete_bent2_100m)%in%bent2_100m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent2_100m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent2_100m=sum(Indice_rarete_bent2_100m_station)                #somme des indices de rareté  de ttes les sp

              
#250m

  b=unique(habitat_par_site_bent2_250m$Code_hab)

    Indice_rarete_bent2_250m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_bent2_250m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent2_250m$Code_Site[habitat_par_site_bent2_250m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_bent2_250m$Code_Site[habitat_par_site_bent2_250m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent2_250m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent2_250m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent2_250m_station<-matrix(0,500,27)
Indice_rarete_bent2_250m_station<-NULL

    for (m in 1:27){
  bent2_250m_station[,m]<-c(as.vector(habitat_par_site_bent2_250m$Code_hab[habitat_par_site_bent2_250m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent2_250m$Code_hab[habitat_par_site_bent2_250m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent2_250m_station[m]=sum(Indice_rarete_bent2_250m[names(Indice_rarete_bent2_250m)%in%bent2_250m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent2_250m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent2_250m=sum(Indice_rarete_bent2_250m_station)                #somme des indices de rareté  de ttes les sp

     
#500m

  b=unique(habitat_par_site_bent2_500m$Code_hab)

    Indice_rarete_bent2_500m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_bent2_500m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent2_500m$Code_Site[habitat_par_site_bent2_500m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_bent2_500m$Code_Site[habitat_par_site_bent2_500m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent2_500m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent2_500m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent2_500m_station<-matrix(0,500,27)
Indice_rarete_bent2_500m_station<-NULL

    for (m in 1:27){
  bent2_500m_station[,m]<-c(as.vector(habitat_par_site_bent2_500m$Code_hab[habitat_par_site_bent2_500m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent2_500m$Code_hab[habitat_par_site_bent2_500m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent2_500m_station[m]=sum(Indice_rarete_bent2_500m[names(Indice_rarete_bent2_500m)%in%bent2_500m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent2_500m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent2_500m=sum(Indice_rarete_bent2_500m_station)                #somme des indices de rareté  de ttes les sp

     
#750m

  b=unique(habitat_par_site_bent2_750m$Code_hab)

    Indice_rarete_bent2_750m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_bent2_750m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent2_750m$Code_Site[habitat_par_site_bent2_750m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_bent2_750m$Code_Site[habitat_par_site_bent2_750m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent2_750m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent2_750m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent2_750m_station<-matrix(0,500,27)
Indice_rarete_bent2_750m_station<-NULL

    for (m in 1:27){
  bent2_750m_station[,m]<-c(as.vector(habitat_par_site_bent2_750m$Code_hab[habitat_par_site_bent2_750m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent2_750m$Code_hab[habitat_par_site_bent2_750m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent2_750m_station[m]=sum(Indice_rarete_bent2_750m[names(Indice_rarete_bent2_750m)%in%bent2_750m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent2_750m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent2_750m=sum(Indice_rarete_bent2_750m_station)                #somme des indices de rareté  de ttes les sp

      
#1000m

  b=unique(habitat_par_site_bent2_1000m$Code_hab)

    Indice_rarete_bent2_1000m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_bent2_1000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent2_1000m$Code_Site[habitat_par_site_bent2_1000m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_bent2_1000m$Code_Site[habitat_par_site_bent2_1000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent2_1000m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent2_1000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent2_1000m_station<-matrix(0,500,27)
Indice_rarete_bent2_1000m_station<-NULL

    for (m in 1:27){
  bent2_1000m_station[,m]<-c(as.vector(habitat_par_site_bent2_1000m$Code_hab[habitat_par_site_bent2_1000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent2_1000m$Code_hab[habitat_par_site_bent2_1000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent2_1000m_station[m]=sum(Indice_rarete_bent2_1000m[names(Indice_rarete_bent2_1000m)%in%bent2_1000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent2_1000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent2_1000m=sum(Indice_rarete_bent2_1000m_station)                #somme des indices de rareté  de ttes les sp

                  
#3000m

  b=unique(habitat_par_site_bent2_3000m$Code_hab)

    Indice_rarete_bent2_3000m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_bent2_3000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent2_3000m$Code_Site[habitat_par_site_bent2_3000m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_bent2_3000m$Code_Site[habitat_par_site_bent2_3000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent2_3000m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent2_3000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent2_3000m_station<-matrix(0,500,27)
Indice_rarete_bent2_3000m_station<-NULL

    for (m in 1:27){
  bent2_3000m_station[,m]<-c(as.vector(habitat_par_site_bent2_3000m$Code_hab[habitat_par_site_bent2_3000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent2_3000m$Code_hab[habitat_par_site_bent2_3000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent2_3000m_station[m]=sum(Indice_rarete_bent2_3000m[names(Indice_rarete_bent2_3000m)%in%bent2_3000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent2_3000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent2_3000m=sum(Indice_rarete_bent2_3000m_station)                #somme des indices de rareté  de ttes les sp

         
#5000m

  b=unique(habitat_par_site_bent2_5000m$Code_hab)

    Indice_rarete_bent2_5000m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_bent2_5000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent2_5000m$Code_Site[habitat_par_site_bent2_5000m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_bent2_5000m$Code_Site[habitat_par_site_bent2_5000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent2_5000m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent2_5000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent2_5000m_station<-matrix(0,500,27)
Indice_rarete_bent2_5000m_station<-NULL

    for (m in 1:27){
  bent2_5000m_station[,m]<-c(as.vector(habitat_par_site_bent2_5000m$Code_hab[habitat_par_site_bent2_5000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent2_5000m$Code_hab[habitat_par_site_bent2_5000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent2_5000m_station[m]=sum(Indice_rarete_bent2_5000m[names(Indice_rarete_bent2_5000m)%in%bent2_5000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent2_5000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent2_5000m=sum(Indice_rarete_bent2_5000m_station)                #somme des indices de rareté  de ttes les sp


source('script.rarete.complementarite.r')

#Selection des dites selon un scenario de rareté complementarité basé sur les habitats

#habitat géo1
selection.rarete.complementarite.bent2.30m<-run.bin(habitat_par_site_bent2_30m_pres_abs,Indice_rarete_bent2_30m_station)
selection.rarete.complementarite.bent2.60m<-run.bin(habitat_par_site_bent2_60m_pres_abs,Indice_rarete_bent2_60m_station)
selection.rarete.complementarite.bent2.100m<-run.bin(habitat_par_site_bent2_100m_pres_abs,Indice_rarete_bent2_100m_station)
selection.rarete.complementarite.bent2.250m<-run.bin(habitat_par_site_bent2_250m_pres_abs,Indice_rarete_bent2_250m_station)
selection.rarete.complementarite.bent2.500m<-run.bin(habitat_par_site_bent2_500m_pres_abs,Indice_rarete_bent2_500m_station)
selection.rarete.complementarite.bent2.750m<-run.bin(habitat_par_site_bent2_750m_pres_abs,Indice_rarete_bent2_750m_station)
selection.rarete.complementarite.bent2.1000m<-run.bin(habitat_par_site_bent2_1000m_pres_abs,Indice_rarete_bent2_1000m_station)
selection.rarete.complementarite.bent2.3000m<-run.bin(habitat_par_site_bent2_3000m_pres_abs,Indice_rarete_bent2_3000m_station)
selection.rarete.complementarite.bent2.5000m<-run.bin(habitat_par_site_bent2_5000m_pres_abs,Indice_rarete_bent2_5000m_station)

#DIVERSITE SELECTIONNE PAR LE SCRIPT RARETE COMPLEMENTARITE#
 source('script.diversite.rarete.selectionnee.r')
 
 diversite.selection.rarete.complementarite.bent2.30m<-species.rarete.evol(selection.rarete.complementarite.bent2.30m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent2.60m<-species.rarete.evol(selection.rarete.complementarite.bent2.60m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent2.100m<-species.rarete.evol(selection.rarete.complementarite.bent2.100m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent2.250m<-species.rarete.evol(selection.rarete.complementarite.bent2.250m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent2.500m<-species.rarete.evol(selection.rarete.complementarite.bent2.500m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent2.750m<-species.rarete.evol(selection.rarete.complementarite.bent2.750m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent2.1000m<-species.rarete.evol(selection.rarete.complementarite.bent2.1000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent2.3000m<-species.rarete.evol(selection.rarete.complementarite.bent2.3000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent2.5000m<-species.rarete.evol(selection.rarete.complementarite.bent2.5000m,Indice_rarete_station)

#CALCUL DES SAI


source('script.sai.r')
SAI_rarete.complementarite.bent2.30m<-SAI(diversite.selection.rarete.complementarite.bent2.30m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent2.60m<-SAI(diversite.selection.rarete.complementarite.bent2.60m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent2.100m<-SAI(diversite.selection.rarete.complementarite.bent2.100m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent2.250m<-SAI(diversite.selection.rarete.complementarite.bent2.250m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent2.500m<-SAI(diversite.selection.rarete.complementarite.bent2.500m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent2.750m<-SAI(diversite.selection.rarete.complementarite.bent2.750m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent2.1000m<-SAI(diversite.selection.rarete.complementarite.bent2.1000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent2.3000m<-SAI(diversite.selection.rarete.complementarite.bent2.3000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent2.5000m<-SAI(diversite.selection.rarete.complementarite.bent2.5000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)



 
############################HABITAT bent3###########################################


 ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances


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



#Load de la requette "REQ_site_bent3_30m.csv"

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

#calcul indice de rareté des habitats

#30m

  b=unique(habitat_par_site_bent3_30m$Code_hab)

    Indice_rarete_bent3_30m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_bent3_30m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent3_30m$Code_Site[habitat_par_site_bent3_30m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_bent3_30m$Code_Site[habitat_par_site_bent3_30m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent3_30m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent3_30m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent3_30m_station<-matrix(0,500,27)
Indice_rarete_bent3_30m_station<-NULL

    for (m in 1:27){
  bent3_30m_station[,m]<-c(as.vector(habitat_par_site_bent3_30m$Code_hab[habitat_par_site_bent3_30m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent3_30m$Code_hab[habitat_par_site_bent3_30m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent3_30m_station[m]=sum(Indice_rarete_bent3_30m[names(Indice_rarete_bent3_30m)%in%bent3_30m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent3_30m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent3_30m=sum(Indice_rarete_bent3_30m_station)                #somme des indices de rareté  de ttes les sp


   
#60m

  b=unique(habitat_par_site_bent3_60m$Code_hab)

    Indice_rarete_bent3_60m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_bent3_60m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent3_60m$Code_Site[habitat_par_site_bent3_60m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_bent3_60m$Code_Site[habitat_par_site_bent3_60m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent3_60m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent3_60m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent3_60m_station<-matrix(0,500,27)
Indice_rarete_bent3_60m_station<-NULL

    for (m in 1:27){
  bent3_60m_station[,m]<-c(as.vector(habitat_par_site_bent3_60m$Code_hab[habitat_par_site_bent3_60m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent3_60m$Code_hab[habitat_par_site_bent3_60m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent3_60m_station[m]=sum(Indice_rarete_bent3_60m[names(Indice_rarete_bent3_60m)%in%bent3_60m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent3_60m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent3_60m=sum(Indice_rarete_bent3_60m_station)                #somme des indices de rareté  de ttes les sp


    
#100m

  b=unique(habitat_par_site_bent3_100m$Code_hab)

    Indice_rarete_bent3_100m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_bent3_100m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent3_100m$Code_Site[habitat_par_site_bent3_100m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_bent3_100m$Code_Site[habitat_par_site_bent3_100m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent3_100m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent3_100m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent3_100m_station<-matrix(0,500,27)
Indice_rarete_bent3_100m_station<-NULL

    for (m in 1:27){
  bent3_100m_station[,m]<-c(as.vector(habitat_par_site_bent3_100m$Code_hab[habitat_par_site_bent3_100m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent3_100m$Code_hab[habitat_par_site_bent3_100m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent3_100m_station[m]=sum(Indice_rarete_bent3_100m[names(Indice_rarete_bent3_100m)%in%bent3_100m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent3_100m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent3_100m=sum(Indice_rarete_bent3_100m_station)                #somme des indices de rareté  de ttes les sp

              
#250m

  b=unique(habitat_par_site_bent3_250m$Code_hab)

    Indice_rarete_bent3_250m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_bent3_250m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent3_250m$Code_Site[habitat_par_site_bent3_250m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_bent3_250m$Code_Site[habitat_par_site_bent3_250m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent3_250m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent3_250m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent3_250m_station<-matrix(0,500,27)
Indice_rarete_bent3_250m_station<-NULL

    for (m in 1:27){
  bent3_250m_station[,m]<-c(as.vector(habitat_par_site_bent3_250m$Code_hab[habitat_par_site_bent3_250m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent3_250m$Code_hab[habitat_par_site_bent3_250m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent3_250m_station[m]=sum(Indice_rarete_bent3_250m[names(Indice_rarete_bent3_250m)%in%bent3_250m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent3_250m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent3_250m=sum(Indice_rarete_bent3_250m_station)                #somme des indices de rareté  de ttes les sp

     
#500m

  b=unique(habitat_par_site_bent3_500m$Code_hab)

    Indice_rarete_bent3_500m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_bent3_500m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent3_500m$Code_Site[habitat_par_site_bent3_500m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_bent3_500m$Code_Site[habitat_par_site_bent3_500m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent3_500m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent3_500m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent3_500m_station<-matrix(0,500,27)
Indice_rarete_bent3_500m_station<-NULL

    for (m in 1:27){
  bent3_500m_station[,m]<-c(as.vector(habitat_par_site_bent3_500m$Code_hab[habitat_par_site_bent3_500m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent3_500m$Code_hab[habitat_par_site_bent3_500m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent3_500m_station[m]=sum(Indice_rarete_bent3_500m[names(Indice_rarete_bent3_500m)%in%bent3_500m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent3_500m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent3_500m=sum(Indice_rarete_bent3_500m_station)                #somme des indices de rareté  de ttes les sp

     
#750m

  b=unique(habitat_par_site_bent3_750m$Code_hab)

    Indice_rarete_bent3_750m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_bent3_750m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent3_750m$Code_Site[habitat_par_site_bent3_750m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_bent3_750m$Code_Site[habitat_par_site_bent3_750m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent3_750m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent3_750m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent3_750m_station<-matrix(0,500,27)
Indice_rarete_bent3_750m_station<-NULL

    for (m in 1:27){
  bent3_750m_station[,m]<-c(as.vector(habitat_par_site_bent3_750m$Code_hab[habitat_par_site_bent3_750m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent3_750m$Code_hab[habitat_par_site_bent3_750m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent3_750m_station[m]=sum(Indice_rarete_bent3_750m[names(Indice_rarete_bent3_750m)%in%bent3_750m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent3_750m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent3_750m=sum(Indice_rarete_bent3_750m_station)                #somme des indices de rareté  de ttes les sp

      
#1000m

  b=unique(habitat_par_site_bent3_1000m$Code_hab)

    Indice_rarete_bent3_1000m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_bent3_1000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent3_1000m$Code_Site[habitat_par_site_bent3_1000m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_bent3_1000m$Code_Site[habitat_par_site_bent3_1000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent3_1000m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent3_1000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent3_1000m_station<-matrix(0,500,27)
Indice_rarete_bent3_1000m_station<-NULL

    for (m in 1:27){
  bent3_1000m_station[,m]<-c(as.vector(habitat_par_site_bent3_1000m$Code_hab[habitat_par_site_bent3_1000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent3_1000m$Code_hab[habitat_par_site_bent3_1000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent3_1000m_station[m]=sum(Indice_rarete_bent3_1000m[names(Indice_rarete_bent3_1000m)%in%bent3_1000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent3_1000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent3_1000m=sum(Indice_rarete_bent3_1000m_station)                #somme des indices de rareté  de ttes les sp

                  
#3000m

  b=unique(habitat_par_site_bent3_3000m$Code_hab)

    Indice_rarete_bent3_3000m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_bent3_3000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent3_3000m$Code_Site[habitat_par_site_bent3_3000m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_bent3_3000m$Code_Site[habitat_par_site_bent3_3000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent3_3000m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent3_3000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent3_3000m_station<-matrix(0,500,27)
Indice_rarete_bent3_3000m_station<-NULL

    for (m in 1:27){
  bent3_3000m_station[,m]<-c(as.vector(habitat_par_site_bent3_3000m$Code_hab[habitat_par_site_bent3_3000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent3_3000m$Code_hab[habitat_par_site_bent3_3000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent3_3000m_station[m]=sum(Indice_rarete_bent3_3000m[names(Indice_rarete_bent3_3000m)%in%bent3_3000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent3_3000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent3_3000m=sum(Indice_rarete_bent3_3000m_station)                #somme des indices de rareté  de ttes les sp

         
#5000m

  b=unique(habitat_par_site_bent3_5000m$Code_hab)

    Indice_rarete_bent3_5000m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_bent3_5000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent3_5000m$Code_Site[habitat_par_site_bent3_5000m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_bent3_5000m$Code_Site[habitat_par_site_bent3_5000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent3_5000m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent3_5000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent3_5000m_station<-matrix(0,500,27)
Indice_rarete_bent3_5000m_station<-NULL

    for (m in 1:27){
  bent3_5000m_station[,m]<-c(as.vector(habitat_par_site_bent3_5000m$Code_hab[habitat_par_site_bent3_5000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent3_5000m$Code_hab[habitat_par_site_bent3_5000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent3_5000m_station[m]=sum(Indice_rarete_bent3_5000m[names(Indice_rarete_bent3_5000m)%in%bent3_5000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent3_5000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent3_5000m=sum(Indice_rarete_bent3_5000m_station)                #somme des indices de rareté  de ttes les sp


source('script.rarete.complementarite.r')

#Selection des dites selon un scenario de rareté complementarité basé sur les habitats

#habitat géo1
selection.rarete.complementarite.bent3.30m<-run.bin(habitat_par_site_bent3_30m_pres_abs,Indice_rarete_bent3_30m_station)
selection.rarete.complementarite.bent3.60m<-run.bin(habitat_par_site_bent3_60m_pres_abs,Indice_rarete_bent3_60m_station)
selection.rarete.complementarite.bent3.100m<-run.bin(habitat_par_site_bent3_100m_pres_abs,Indice_rarete_bent3_100m_station)
selection.rarete.complementarite.bent3.250m<-run.bin(habitat_par_site_bent3_250m_pres_abs,Indice_rarete_bent3_250m_station)
selection.rarete.complementarite.bent3.500m<-run.bin(habitat_par_site_bent3_500m_pres_abs,Indice_rarete_bent3_500m_station)
selection.rarete.complementarite.bent3.750m<-run.bin(habitat_par_site_bent3_750m_pres_abs,Indice_rarete_bent3_750m_station)
selection.rarete.complementarite.bent3.1000m<-run.bin(habitat_par_site_bent3_1000m_pres_abs,Indice_rarete_bent3_1000m_station)
selection.rarete.complementarite.bent3.3000m<-run.bin(habitat_par_site_bent3_3000m_pres_abs,Indice_rarete_bent3_3000m_station)
selection.rarete.complementarite.bent3.5000m<-run.bin(habitat_par_site_bent3_5000m_pres_abs,Indice_rarete_bent3_5000m_station)

#DIVERSITE SELECTIONNE PAR LE SCRIPT RARETE COMPLEMENTARITE#
 source('script.diversite.rarete.selectionnee.r')
 
 diversite.selection.rarete.complementarite.bent3.30m<-species.rarete.evol(selection.rarete.complementarite.bent3.30m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent3.60m<-species.rarete.evol(selection.rarete.complementarite.bent3.60m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent3.100m<-species.rarete.evol(selection.rarete.complementarite.bent3.100m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent3.250m<-species.rarete.evol(selection.rarete.complementarite.bent3.250m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent3.500m<-species.rarete.evol(selection.rarete.complementarite.bent3.500m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent3.750m<-species.rarete.evol(selection.rarete.complementarite.bent3.750m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent3.1000m<-species.rarete.evol(selection.rarete.complementarite.bent3.1000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent3.3000m<-species.rarete.evol(selection.rarete.complementarite.bent3.3000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent3.5000m<-species.rarete.evol(selection.rarete.complementarite.bent3.5000m,Indice_rarete_station)

#CALCUL DES SAI


source('script.sai.r')
SAI_rarete.complementarite.bent3.30m<-SAI(diversite.selection.rarete.complementarite.bent3.30m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent3.60m<-SAI(diversite.selection.rarete.complementarite.bent3.60m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent3.100m<-SAI(diversite.selection.rarete.complementarite.bent3.100m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent3.250m<-SAI(diversite.selection.rarete.complementarite.bent3.250m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent3.500m<-SAI(diversite.selection.rarete.complementarite.bent3.500m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent3.750m<-SAI(diversite.selection.rarete.complementarite.bent3.750m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent3.1000m<-SAI(diversite.selection.rarete.complementarite.bent3.1000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent3.3000m<-SAI(diversite.selection.rarete.complementarite.bent3.3000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent3.5000m<-SAI(diversite.selection.rarete.complementarite.bent3.5000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)



 
############################HABITAT geo1_geo2###########################################


 ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances


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



#Load de la requette "REQ_site_geo1_geo2_30m.csv"

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

#calcul indice de rareté des habitats

#30m

  b=unique(habitat_par_site_geo1_geo2_30m$Code_hab)

    Indice_rarete_geo1_geo2_30m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo1_geo2_30m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_30m$Code_Site[habitat_par_site_geo1_geo2_30m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo1_geo2_30m$Code_Site[habitat_par_site_geo1_geo2_30m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_30m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_30m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_30m_station<-matrix(0,500,27)
Indice_rarete_geo1_geo2_30m_station<-NULL

    for (m in 1:27){
  geo1_geo2_30m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_30m$Code_hab[habitat_par_site_geo1_geo2_30m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_30m$Code_hab[habitat_par_site_geo1_geo2_30m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_30m_station[m]=sum(Indice_rarete_geo1_geo2_30m[names(Indice_rarete_geo1_geo2_30m)%in%geo1_geo2_30m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_30m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_30m=sum(Indice_rarete_geo1_geo2_30m_station)                #somme des indices de rareté  de ttes les sp


   
#60m

  b=unique(habitat_par_site_geo1_geo2_60m$Code_hab)

    Indice_rarete_geo1_geo2_60m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo1_geo2_60m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_60m$Code_Site[habitat_par_site_geo1_geo2_60m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo1_geo2_60m$Code_Site[habitat_par_site_geo1_geo2_60m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_60m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_60m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_60m_station<-matrix(0,500,27)
Indice_rarete_geo1_geo2_60m_station<-NULL

    for (m in 1:27){
  geo1_geo2_60m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_60m$Code_hab[habitat_par_site_geo1_geo2_60m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_60m$Code_hab[habitat_par_site_geo1_geo2_60m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_60m_station[m]=sum(Indice_rarete_geo1_geo2_60m[names(Indice_rarete_geo1_geo2_60m)%in%geo1_geo2_60m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_60m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_60m=sum(Indice_rarete_geo1_geo2_60m_station)                #somme des indices de rareté  de ttes les sp


    
#100m

  b=unique(habitat_par_site_geo1_geo2_100m$Code_hab)

    Indice_rarete_geo1_geo2_100m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo1_geo2_100m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_100m$Code_Site[habitat_par_site_geo1_geo2_100m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo1_geo2_100m$Code_Site[habitat_par_site_geo1_geo2_100m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_100m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_100m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_100m_station<-matrix(0,500,27)
Indice_rarete_geo1_geo2_100m_station<-NULL

    for (m in 1:27){
  geo1_geo2_100m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_100m$Code_hab[habitat_par_site_geo1_geo2_100m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_100m$Code_hab[habitat_par_site_geo1_geo2_100m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_100m_station[m]=sum(Indice_rarete_geo1_geo2_100m[names(Indice_rarete_geo1_geo2_100m)%in%geo1_geo2_100m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_100m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_100m=sum(Indice_rarete_geo1_geo2_100m_station)                #somme des indices de rareté  de ttes les sp

              
#250m

  b=unique(habitat_par_site_geo1_geo2_250m$Code_hab)

    Indice_rarete_geo1_geo2_250m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo1_geo2_250m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_250m$Code_Site[habitat_par_site_geo1_geo2_250m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo1_geo2_250m$Code_Site[habitat_par_site_geo1_geo2_250m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_250m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_250m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_250m_station<-matrix(0,500,27)
Indice_rarete_geo1_geo2_250m_station<-NULL

    for (m in 1:27){
  geo1_geo2_250m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_250m$Code_hab[habitat_par_site_geo1_geo2_250m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_250m$Code_hab[habitat_par_site_geo1_geo2_250m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_250m_station[m]=sum(Indice_rarete_geo1_geo2_250m[names(Indice_rarete_geo1_geo2_250m)%in%geo1_geo2_250m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_250m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_250m=sum(Indice_rarete_geo1_geo2_250m_station)                #somme des indices de rareté  de ttes les sp

     
#500m

  b=unique(habitat_par_site_geo1_geo2_500m$Code_hab)

    Indice_rarete_geo1_geo2_500m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo1_geo2_500m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_500m$Code_Site[habitat_par_site_geo1_geo2_500m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo1_geo2_500m$Code_Site[habitat_par_site_geo1_geo2_500m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_500m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_500m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_500m_station<-matrix(0,500,27)
Indice_rarete_geo1_geo2_500m_station<-NULL

    for (m in 1:27){
  geo1_geo2_500m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_500m$Code_hab[habitat_par_site_geo1_geo2_500m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_500m$Code_hab[habitat_par_site_geo1_geo2_500m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_500m_station[m]=sum(Indice_rarete_geo1_geo2_500m[names(Indice_rarete_geo1_geo2_500m)%in%geo1_geo2_500m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_500m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_500m=sum(Indice_rarete_geo1_geo2_500m_station)                #somme des indices de rareté  de ttes les sp

     
#750m

  b=unique(habitat_par_site_geo1_geo2_750m$Code_hab)

    Indice_rarete_geo1_geo2_750m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo1_geo2_750m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_750m$Code_Site[habitat_par_site_geo1_geo2_750m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo1_geo2_750m$Code_Site[habitat_par_site_geo1_geo2_750m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_750m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_750m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_750m_station<-matrix(0,500,27)
Indice_rarete_geo1_geo2_750m_station<-NULL

    for (m in 1:27){
  geo1_geo2_750m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_750m$Code_hab[habitat_par_site_geo1_geo2_750m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_750m$Code_hab[habitat_par_site_geo1_geo2_750m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_750m_station[m]=sum(Indice_rarete_geo1_geo2_750m[names(Indice_rarete_geo1_geo2_750m)%in%geo1_geo2_750m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_750m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_750m=sum(Indice_rarete_geo1_geo2_750m_station)                #somme des indices de rareté  de ttes les sp

      
#1000m

  b=unique(habitat_par_site_geo1_geo2_1000m$Code_hab)

    Indice_rarete_geo1_geo2_1000m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo1_geo2_1000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_1000m$Code_Site[habitat_par_site_geo1_geo2_1000m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo1_geo2_1000m$Code_Site[habitat_par_site_geo1_geo2_1000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_1000m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_1000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_1000m_station<-matrix(0,500,27)
Indice_rarete_geo1_geo2_1000m_station<-NULL

    for (m in 1:27){
  geo1_geo2_1000m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_1000m$Code_hab[habitat_par_site_geo1_geo2_1000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_1000m$Code_hab[habitat_par_site_geo1_geo2_1000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_1000m_station[m]=sum(Indice_rarete_geo1_geo2_1000m[names(Indice_rarete_geo1_geo2_1000m)%in%geo1_geo2_1000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_1000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_1000m=sum(Indice_rarete_geo1_geo2_1000m_station)                #somme des indices de rareté  de ttes les sp

                  
#3000m

  b=unique(habitat_par_site_geo1_geo2_3000m$Code_hab)

    Indice_rarete_geo1_geo2_3000m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo1_geo2_3000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_3000m$Code_Site[habitat_par_site_geo1_geo2_3000m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo1_geo2_3000m$Code_Site[habitat_par_site_geo1_geo2_3000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_3000m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_3000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_3000m_station<-matrix(0,500,27)
Indice_rarete_geo1_geo2_3000m_station<-NULL

    for (m in 1:27){
  geo1_geo2_3000m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_3000m$Code_hab[habitat_par_site_geo1_geo2_3000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_3000m$Code_hab[habitat_par_site_geo1_geo2_3000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_3000m_station[m]=sum(Indice_rarete_geo1_geo2_3000m[names(Indice_rarete_geo1_geo2_3000m)%in%geo1_geo2_3000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_3000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_3000m=sum(Indice_rarete_geo1_geo2_3000m_station)                #somme des indices de rareté  de ttes les sp

         
#5000m

  b=unique(habitat_par_site_geo1_geo2_5000m$Code_hab)

    Indice_rarete_geo1_geo2_5000m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo1_geo2_5000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_5000m$Code_Site[habitat_par_site_geo1_geo2_5000m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo1_geo2_5000m$Code_Site[habitat_par_site_geo1_geo2_5000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_5000m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_5000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_5000m_station<-matrix(0,500,27)
Indice_rarete_geo1_geo2_5000m_station<-NULL

    for (m in 1:27){
  geo1_geo2_5000m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_5000m$Code_hab[habitat_par_site_geo1_geo2_5000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_5000m$Code_hab[habitat_par_site_geo1_geo2_5000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_5000m_station[m]=sum(Indice_rarete_geo1_geo2_5000m[names(Indice_rarete_geo1_geo2_5000m)%in%geo1_geo2_5000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_5000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_5000m=sum(Indice_rarete_geo1_geo2_5000m_station)                #somme des indices de rareté  de ttes les sp


source('script.rarete.complementarite.r')

#Selection des dites selon un scenario de rareté complementarité basé sur les habitats

#habitat géo1
selection.rarete.complementarite.geo1_geo2.30m<-run.bin(habitat_par_site_geo1_geo2_30m_pres_abs,Indice_rarete_geo1_geo2_30m_station)
selection.rarete.complementarite.geo1_geo2.60m<-run.bin(habitat_par_site_geo1_geo2_60m_pres_abs,Indice_rarete_geo1_geo2_60m_station)
selection.rarete.complementarite.geo1_geo2.100m<-run.bin(habitat_par_site_geo1_geo2_100m_pres_abs,Indice_rarete_geo1_geo2_100m_station)
selection.rarete.complementarite.geo1_geo2.250m<-run.bin(habitat_par_site_geo1_geo2_250m_pres_abs,Indice_rarete_geo1_geo2_250m_station)
selection.rarete.complementarite.geo1_geo2.500m<-run.bin(habitat_par_site_geo1_geo2_500m_pres_abs,Indice_rarete_geo1_geo2_500m_station)
selection.rarete.complementarite.geo1_geo2.750m<-run.bin(habitat_par_site_geo1_geo2_750m_pres_abs,Indice_rarete_geo1_geo2_750m_station)
selection.rarete.complementarite.geo1_geo2.1000m<-run.bin(habitat_par_site_geo1_geo2_1000m_pres_abs,Indice_rarete_geo1_geo2_1000m_station)
selection.rarete.complementarite.geo1_geo2.3000m<-run.bin(habitat_par_site_geo1_geo2_3000m_pres_abs,Indice_rarete_geo1_geo2_3000m_station)
selection.rarete.complementarite.geo1_geo2.5000m<-run.bin(habitat_par_site_geo1_geo2_5000m_pres_abs,Indice_rarete_geo1_geo2_5000m_station)

#DIVERSITE SELECTIONNE PAR LE SCRIPT RARETE COMPLEMENTARITE#
 source('script.diversite.rarete.selectionnee.r')
 
 diversite.selection.rarete.complementarite.geo1_geo2.30m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2.30m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1_geo2.60m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2.60m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1_geo2.100m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2.100m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1_geo2.250m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2.250m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1_geo2.500m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2.500m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1_geo2.750m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2.750m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1_geo2.1000m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2.1000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1_geo2.3000m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2.3000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1_geo2.5000m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2.5000m,Indice_rarete_station)

#CALCUL DES SAI


source('script.sai.r')
SAI_rarete.complementarite.geo1_geo2.30m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2.30m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1_geo2.60m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2.60m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1_geo2.100m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2.100m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1_geo2.250m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2.250m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1_geo2.500m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2.500m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1_geo2.750m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2.750m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1_geo2.1000m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2.1000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1_geo2.3000m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2.3000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1_geo2.5000m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2.5000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)



 
############################HABITAT geo1_geo2_geo3###########################################


 ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances


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



#Load de la requette "REQ_site_geo1_geo2_geo3_30m.csv"

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

#calcul indice de rareté des habitats

#30m

  b=unique(habitat_par_site_geo1_geo2_geo3_30m$Code_hab)

    Indice_rarete_geo1_geo2_geo3_30m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo1_geo2_geo3_30m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_30m$Code_Site[habitat_par_site_geo1_geo2_geo3_30m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo1_geo2_geo3_30m$Code_Site[habitat_par_site_geo1_geo2_geo3_30m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_geo3_30m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_geo3_30m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_geo3_30m_station<-matrix(0,500,27)
Indice_rarete_geo1_geo2_geo3_30m_station<-NULL

    for (m in 1:27){
  geo1_geo2_geo3_30m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_30m$Code_hab[habitat_par_site_geo1_geo2_geo3_30m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_geo3_30m$Code_hab[habitat_par_site_geo1_geo2_geo3_30m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_geo3_30m_station[m]=sum(Indice_rarete_geo1_geo2_geo3_30m[names(Indice_rarete_geo1_geo2_geo3_30m)%in%geo1_geo2_geo3_30m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_geo3_30m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_geo3_30m=sum(Indice_rarete_geo1_geo2_geo3_30m_station)                #somme des indices de rareté  de ttes les sp


   
#60m

  b=unique(habitat_par_site_geo1_geo2_geo3_60m$Code_hab)

    Indice_rarete_geo1_geo2_geo3_60m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo1_geo2_geo3_60m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_60m$Code_Site[habitat_par_site_geo1_geo2_geo3_60m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo1_geo2_geo3_60m$Code_Site[habitat_par_site_geo1_geo2_geo3_60m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_geo3_60m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_geo3_60m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_geo3_60m_station<-matrix(0,500,27)
Indice_rarete_geo1_geo2_geo3_60m_station<-NULL

    for (m in 1:27){
  geo1_geo2_geo3_60m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_60m$Code_hab[habitat_par_site_geo1_geo2_geo3_60m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_geo3_60m$Code_hab[habitat_par_site_geo1_geo2_geo3_60m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_geo3_60m_station[m]=sum(Indice_rarete_geo1_geo2_geo3_60m[names(Indice_rarete_geo1_geo2_geo3_60m)%in%geo1_geo2_geo3_60m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_geo3_60m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_geo3_60m=sum(Indice_rarete_geo1_geo2_geo3_60m_station)                #somme des indices de rareté  de ttes les sp


    
#100m

  b=unique(habitat_par_site_geo1_geo2_geo3_100m$Code_hab)

    Indice_rarete_geo1_geo2_geo3_100m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo1_geo2_geo3_100m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_100m$Code_Site[habitat_par_site_geo1_geo2_geo3_100m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo1_geo2_geo3_100m$Code_Site[habitat_par_site_geo1_geo2_geo3_100m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_geo3_100m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_geo3_100m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_geo3_100m_station<-matrix(0,500,27)
Indice_rarete_geo1_geo2_geo3_100m_station<-NULL

    for (m in 1:27){
  geo1_geo2_geo3_100m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_100m$Code_hab[habitat_par_site_geo1_geo2_geo3_100m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_geo3_100m$Code_hab[habitat_par_site_geo1_geo2_geo3_100m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_geo3_100m_station[m]=sum(Indice_rarete_geo1_geo2_geo3_100m[names(Indice_rarete_geo1_geo2_geo3_100m)%in%geo1_geo2_geo3_100m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_geo3_100m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_geo3_100m=sum(Indice_rarete_geo1_geo2_geo3_100m_station)                #somme des indices de rareté  de ttes les sp

              
#250m

  b=unique(habitat_par_site_geo1_geo2_geo3_250m$Code_hab)

    Indice_rarete_geo1_geo2_geo3_250m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo1_geo2_geo3_250m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_250m$Code_Site[habitat_par_site_geo1_geo2_geo3_250m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo1_geo2_geo3_250m$Code_Site[habitat_par_site_geo1_geo2_geo3_250m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_geo3_250m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_geo3_250m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_geo3_250m_station<-matrix(0,500,27)
Indice_rarete_geo1_geo2_geo3_250m_station<-NULL

    for (m in 1:27){
  geo1_geo2_geo3_250m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_250m$Code_hab[habitat_par_site_geo1_geo2_geo3_250m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_geo3_250m$Code_hab[habitat_par_site_geo1_geo2_geo3_250m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_geo3_250m_station[m]=sum(Indice_rarete_geo1_geo2_geo3_250m[names(Indice_rarete_geo1_geo2_geo3_250m)%in%geo1_geo2_geo3_250m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_geo3_250m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_geo3_250m=sum(Indice_rarete_geo1_geo2_geo3_250m_station)                #somme des indices de rareté  de ttes les sp

     
#500m

  b=unique(habitat_par_site_geo1_geo2_geo3_500m$Code_hab)

    Indice_rarete_geo1_geo2_geo3_500m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo1_geo2_geo3_500m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_500m$Code_Site[habitat_par_site_geo1_geo2_geo3_500m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo1_geo2_geo3_500m$Code_Site[habitat_par_site_geo1_geo2_geo3_500m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_geo3_500m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_geo3_500m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_geo3_500m_station<-matrix(0,500,27)
Indice_rarete_geo1_geo2_geo3_500m_station<-NULL

    for (m in 1:27){
  geo1_geo2_geo3_500m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_500m$Code_hab[habitat_par_site_geo1_geo2_geo3_500m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_geo3_500m$Code_hab[habitat_par_site_geo1_geo2_geo3_500m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_geo3_500m_station[m]=sum(Indice_rarete_geo1_geo2_geo3_500m[names(Indice_rarete_geo1_geo2_geo3_500m)%in%geo1_geo2_geo3_500m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_geo3_500m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_geo3_500m=sum(Indice_rarete_geo1_geo2_geo3_500m_station)                #somme des indices de rareté  de ttes les sp

     
#750m

  b=unique(habitat_par_site_geo1_geo2_geo3_750m$Code_hab)

    Indice_rarete_geo1_geo2_geo3_750m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo1_geo2_geo3_750m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_750m$Code_Site[habitat_par_site_geo1_geo2_geo3_750m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo1_geo2_geo3_750m$Code_Site[habitat_par_site_geo1_geo2_geo3_750m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_geo3_750m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_geo3_750m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_geo3_750m_station<-matrix(0,500,27)
Indice_rarete_geo1_geo2_geo3_750m_station<-NULL

    for (m in 1:27){
  geo1_geo2_geo3_750m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_750m$Code_hab[habitat_par_site_geo1_geo2_geo3_750m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_geo3_750m$Code_hab[habitat_par_site_geo1_geo2_geo3_750m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_geo3_750m_station[m]=sum(Indice_rarete_geo1_geo2_geo3_750m[names(Indice_rarete_geo1_geo2_geo3_750m)%in%geo1_geo2_geo3_750m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_geo3_750m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_geo3_750m=sum(Indice_rarete_geo1_geo2_geo3_750m_station)                #somme des indices de rareté  de ttes les sp

      
#1000m

  b=unique(habitat_par_site_geo1_geo2_geo3_1000m$Code_hab)

    Indice_rarete_geo1_geo2_geo3_1000m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo1_geo2_geo3_1000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_1000m$Code_Site[habitat_par_site_geo1_geo2_geo3_1000m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo1_geo2_geo3_1000m$Code_Site[habitat_par_site_geo1_geo2_geo3_1000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_geo3_1000m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_geo3_1000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_geo3_1000m_station<-matrix(0,500,27)
Indice_rarete_geo1_geo2_geo3_1000m_station<-NULL

    for (m in 1:27){
  geo1_geo2_geo3_1000m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_1000m$Code_hab[habitat_par_site_geo1_geo2_geo3_1000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_geo3_1000m$Code_hab[habitat_par_site_geo1_geo2_geo3_1000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_geo3_1000m_station[m]=sum(Indice_rarete_geo1_geo2_geo3_1000m[names(Indice_rarete_geo1_geo2_geo3_1000m)%in%geo1_geo2_geo3_1000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_geo3_1000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_geo3_1000m=sum(Indice_rarete_geo1_geo2_geo3_1000m_station)                #somme des indices de rareté  de ttes les sp

                  
#3000m

  b=unique(habitat_par_site_geo1_geo2_geo3_3000m$Code_hab)

    Indice_rarete_geo1_geo2_geo3_3000m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo1_geo2_geo3_3000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_3000m$Code_Site[habitat_par_site_geo1_geo2_geo3_3000m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo1_geo2_geo3_3000m$Code_Site[habitat_par_site_geo1_geo2_geo3_3000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_geo3_3000m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_geo3_3000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_geo3_3000m_station<-matrix(0,500,27)
Indice_rarete_geo1_geo2_geo3_3000m_station<-NULL

    for (m in 1:27){
  geo1_geo2_geo3_3000m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_3000m$Code_hab[habitat_par_site_geo1_geo2_geo3_3000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_geo3_3000m$Code_hab[habitat_par_site_geo1_geo2_geo3_3000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_geo3_3000m_station[m]=sum(Indice_rarete_geo1_geo2_geo3_3000m[names(Indice_rarete_geo1_geo2_geo3_3000m)%in%geo1_geo2_geo3_3000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_geo3_3000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_geo3_3000m=sum(Indice_rarete_geo1_geo2_geo3_3000m_station)                #somme des indices de rareté  de ttes les sp

         
#5000m

  b=unique(habitat_par_site_geo1_geo2_geo3_5000m$Code_hab)

    Indice_rarete_geo1_geo2_geo3_5000m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_geo1_geo2_geo3_5000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_5000m$Code_Site[habitat_par_site_geo1_geo2_geo3_5000m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_geo1_geo2_geo3_5000m$Code_Site[habitat_par_site_geo1_geo2_geo3_5000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_geo3_5000m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_geo3_5000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_geo3_5000m_station<-matrix(0,500,27)
Indice_rarete_geo1_geo2_geo3_5000m_station<-NULL

    for (m in 1:27){
  geo1_geo2_geo3_5000m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_5000m$Code_hab[habitat_par_site_geo1_geo2_geo3_5000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_geo3_5000m$Code_hab[habitat_par_site_geo1_geo2_geo3_5000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_geo3_5000m_station[m]=sum(Indice_rarete_geo1_geo2_geo3_5000m[names(Indice_rarete_geo1_geo2_geo3_5000m)%in%geo1_geo2_geo3_5000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_geo3_5000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_geo3_5000m=sum(Indice_rarete_geo1_geo2_geo3_5000m_station)                #somme des indices de rareté  de ttes les sp


source('script.rarete.complementarite.r')

#Selection des dites selon un scenario de rareté complementarité basé sur les habitats

#habitat géo1
selection.rarete.complementarite.geo1_geo2_geo3.30m<-run.bin(habitat_par_site_geo1_geo2_geo3_30m_pres_abs,Indice_rarete_geo1_geo2_geo3_30m_station)
selection.rarete.complementarite.geo1_geo2_geo3.60m<-run.bin(habitat_par_site_geo1_geo2_geo3_60m_pres_abs,Indice_rarete_geo1_geo2_geo3_60m_station)
selection.rarete.complementarite.geo1_geo2_geo3.100m<-run.bin(habitat_par_site_geo1_geo2_geo3_100m_pres_abs,Indice_rarete_geo1_geo2_geo3_100m_station)
selection.rarete.complementarite.geo1_geo2_geo3.250m<-run.bin(habitat_par_site_geo1_geo2_geo3_250m_pres_abs,Indice_rarete_geo1_geo2_geo3_250m_station)
selection.rarete.complementarite.geo1_geo2_geo3.500m<-run.bin(habitat_par_site_geo1_geo2_geo3_500m_pres_abs,Indice_rarete_geo1_geo2_geo3_500m_station)
selection.rarete.complementarite.geo1_geo2_geo3.750m<-run.bin(habitat_par_site_geo1_geo2_geo3_750m_pres_abs,Indice_rarete_geo1_geo2_geo3_750m_station)
selection.rarete.complementarite.geo1_geo2_geo3.1000m<-run.bin(habitat_par_site_geo1_geo2_geo3_1000m_pres_abs,Indice_rarete_geo1_geo2_geo3_1000m_station)
selection.rarete.complementarite.geo1_geo2_geo3.3000m<-run.bin(habitat_par_site_geo1_geo2_geo3_3000m_pres_abs,Indice_rarete_geo1_geo2_geo3_3000m_station)
selection.rarete.complementarite.geo1_geo2_geo3.5000m<-run.bin(habitat_par_site_geo1_geo2_geo3_5000m_pres_abs,Indice_rarete_geo1_geo2_geo3_5000m_station)

#DIVERSITE SELECTIONNE PAR LE SCRIPT RARETE COMPLEMENTARITE#
 source('script.diversite.rarete.selectionnee.r')
 
 diversite.selection.rarete.complementarite.geo1_geo2_geo3.30m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2_geo3.30m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1_geo2_geo3.60m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2_geo3.60m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1_geo2_geo3.100m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2_geo3.100m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1_geo2_geo3.250m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2_geo3.250m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1_geo2_geo3.500m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2_geo3.500m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1_geo2_geo3.750m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2_geo3.750m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1_geo2_geo3.1000m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2_geo3.1000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1_geo2_geo3.3000m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2_geo3.3000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1_geo2_geo3.5000m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2_geo3.5000m,Indice_rarete_station)

#CALCUL DES SAI


source('script.sai.r')
SAI_rarete.complementarite.geo1_geo2_geo3.30m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2_geo3.30m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1_geo2_geo3.60m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2_geo3.60m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1_geo2_geo3.100m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2_geo3.100m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1_geo2_geo3.250m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2_geo3.250m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1_geo2_geo3.500m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2_geo3.500m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1_geo2_geo3.750m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2_geo3.750m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1_geo2_geo3.1000m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2_geo3.1000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1_geo2_geo3.3000m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2_geo3.3000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1_geo2_geo3.5000m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2_geo3.5000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)



 
############################HABITAT g1_g2_g3_b1###########################################


 ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances


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



#Load de la requette "REQ_site_g1_g2_g3_b1_30m.csv"

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

#calcul indice de rareté des habitats

#30m

  b=unique(habitat_par_site_g1_g2_g3_b1_30m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_30m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_g1_g2_g3_b1_30m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_30m$Code_Site[habitat_par_site_g1_g2_g3_b1_30m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_g1_g2_g3_b1_30m$Code_Site[habitat_par_site_g1_g2_g3_b1_30m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_30m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_30m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_30m_station<-matrix(0,500,27)
Indice_rarete_g1_g2_g3_b1_30m_station<-NULL

    for (m in 1:27){
  g1_g2_g3_b1_30m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_30m$Code_hab[habitat_par_site_g1_g2_g3_b1_30m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_30m$Code_hab[habitat_par_site_g1_g2_g3_b1_30m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_30m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_30m[names(Indice_rarete_g1_g2_g3_b1_30m)%in%g1_g2_g3_b1_30m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_30m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_30m=sum(Indice_rarete_g1_g2_g3_b1_30m_station)                #somme des indices de rareté  de ttes les sp


   
#60m

  b=unique(habitat_par_site_g1_g2_g3_b1_60m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_60m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_g1_g2_g3_b1_60m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_60m$Code_Site[habitat_par_site_g1_g2_g3_b1_60m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_g1_g2_g3_b1_60m$Code_Site[habitat_par_site_g1_g2_g3_b1_60m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_60m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_60m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_60m_station<-matrix(0,500,27)
Indice_rarete_g1_g2_g3_b1_60m_station<-NULL

    for (m in 1:27){
  g1_g2_g3_b1_60m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_60m$Code_hab[habitat_par_site_g1_g2_g3_b1_60m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_60m$Code_hab[habitat_par_site_g1_g2_g3_b1_60m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_60m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_60m[names(Indice_rarete_g1_g2_g3_b1_60m)%in%g1_g2_g3_b1_60m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_60m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_60m=sum(Indice_rarete_g1_g2_g3_b1_60m_station)                #somme des indices de rareté  de ttes les sp


    
#100m

  b=unique(habitat_par_site_g1_g2_g3_b1_100m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_100m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_g1_g2_g3_b1_100m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_100m$Code_Site[habitat_par_site_g1_g2_g3_b1_100m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_g1_g2_g3_b1_100m$Code_Site[habitat_par_site_g1_g2_g3_b1_100m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_100m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_100m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_100m_station<-matrix(0,500,27)
Indice_rarete_g1_g2_g3_b1_100m_station<-NULL

    for (m in 1:27){
  g1_g2_g3_b1_100m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_100m$Code_hab[habitat_par_site_g1_g2_g3_b1_100m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_100m$Code_hab[habitat_par_site_g1_g2_g3_b1_100m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_100m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_100m[names(Indice_rarete_g1_g2_g3_b1_100m)%in%g1_g2_g3_b1_100m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_100m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_100m=sum(Indice_rarete_g1_g2_g3_b1_100m_station)                #somme des indices de rareté  de ttes les sp

              
#250m

  b=unique(habitat_par_site_g1_g2_g3_b1_250m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_250m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_g1_g2_g3_b1_250m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_250m$Code_Site[habitat_par_site_g1_g2_g3_b1_250m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_g1_g2_g3_b1_250m$Code_Site[habitat_par_site_g1_g2_g3_b1_250m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_250m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_250m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_250m_station<-matrix(0,500,27)
Indice_rarete_g1_g2_g3_b1_250m_station<-NULL

    for (m in 1:27){
  g1_g2_g3_b1_250m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_250m$Code_hab[habitat_par_site_g1_g2_g3_b1_250m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_250m$Code_hab[habitat_par_site_g1_g2_g3_b1_250m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_250m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_250m[names(Indice_rarete_g1_g2_g3_b1_250m)%in%g1_g2_g3_b1_250m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_250m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_250m=sum(Indice_rarete_g1_g2_g3_b1_250m_station)                #somme des indices de rareté  de ttes les sp

     
#500m

  b=unique(habitat_par_site_g1_g2_g3_b1_500m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_500m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_g1_g2_g3_b1_500m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_500m$Code_Site[habitat_par_site_g1_g2_g3_b1_500m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_g1_g2_g3_b1_500m$Code_Site[habitat_par_site_g1_g2_g3_b1_500m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_500m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_500m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_500m_station<-matrix(0,500,27)
Indice_rarete_g1_g2_g3_b1_500m_station<-NULL

    for (m in 1:27){
  g1_g2_g3_b1_500m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_500m$Code_hab[habitat_par_site_g1_g2_g3_b1_500m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_500m$Code_hab[habitat_par_site_g1_g2_g3_b1_500m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_500m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_500m[names(Indice_rarete_g1_g2_g3_b1_500m)%in%g1_g2_g3_b1_500m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_500m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_500m=sum(Indice_rarete_g1_g2_g3_b1_500m_station)                #somme des indices de rareté  de ttes les sp

     
#750m

  b=unique(habitat_par_site_g1_g2_g3_b1_750m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_750m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_g1_g2_g3_b1_750m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_750m$Code_Site[habitat_par_site_g1_g2_g3_b1_750m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_g1_g2_g3_b1_750m$Code_Site[habitat_par_site_g1_g2_g3_b1_750m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_750m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_750m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_750m_station<-matrix(0,500,27)
Indice_rarete_g1_g2_g3_b1_750m_station<-NULL

    for (m in 1:27){
  g1_g2_g3_b1_750m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_750m$Code_hab[habitat_par_site_g1_g2_g3_b1_750m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_750m$Code_hab[habitat_par_site_g1_g2_g3_b1_750m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_750m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_750m[names(Indice_rarete_g1_g2_g3_b1_750m)%in%g1_g2_g3_b1_750m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_750m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_750m=sum(Indice_rarete_g1_g2_g3_b1_750m_station)                #somme des indices de rareté  de ttes les sp

      
#1000m

  b=unique(habitat_par_site_g1_g2_g3_b1_1000m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_1000m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_g1_g2_g3_b1_1000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_1000m$Code_Site[habitat_par_site_g1_g2_g3_b1_1000m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_g1_g2_g3_b1_1000m$Code_Site[habitat_par_site_g1_g2_g3_b1_1000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_1000m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_1000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_1000m_station<-matrix(0,500,27)
Indice_rarete_g1_g2_g3_b1_1000m_station<-NULL

    for (m in 1:27){
  g1_g2_g3_b1_1000m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_1000m$Code_hab[habitat_par_site_g1_g2_g3_b1_1000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_1000m$Code_hab[habitat_par_site_g1_g2_g3_b1_1000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_1000m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_1000m[names(Indice_rarete_g1_g2_g3_b1_1000m)%in%g1_g2_g3_b1_1000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_1000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_1000m=sum(Indice_rarete_g1_g2_g3_b1_1000m_station)                #somme des indices de rareté  de ttes les sp

                  
#3000m

  b=unique(habitat_par_site_g1_g2_g3_b1_3000m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_3000m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_g1_g2_g3_b1_3000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_3000m$Code_Site[habitat_par_site_g1_g2_g3_b1_3000m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_g1_g2_g3_b1_3000m$Code_Site[habitat_par_site_g1_g2_g3_b1_3000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_3000m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_3000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_3000m_station<-matrix(0,500,27)
Indice_rarete_g1_g2_g3_b1_3000m_station<-NULL

    for (m in 1:27){
  g1_g2_g3_b1_3000m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_3000m$Code_hab[habitat_par_site_g1_g2_g3_b1_3000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_3000m$Code_hab[habitat_par_site_g1_g2_g3_b1_3000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_3000m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_3000m[names(Indice_rarete_g1_g2_g3_b1_3000m)%in%g1_g2_g3_b1_3000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_3000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_3000m=sum(Indice_rarete_g1_g2_g3_b1_3000m_station)                #somme des indices de rareté  de ttes les sp

         
#5000m

  b=unique(habitat_par_site_g1_g2_g3_b1_5000m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_5000m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_g1_g2_g3_b1_5000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_5000m$Code_Site[habitat_par_site_g1_g2_g3_b1_5000m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_g1_g2_g3_b1_5000m$Code_Site[habitat_par_site_g1_g2_g3_b1_5000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_5000m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_5000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_5000m_station<-matrix(0,500,27)
Indice_rarete_g1_g2_g3_b1_5000m_station<-NULL

    for (m in 1:27){
  g1_g2_g3_b1_5000m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_5000m$Code_hab[habitat_par_site_g1_g2_g3_b1_5000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_5000m$Code_hab[habitat_par_site_g1_g2_g3_b1_5000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_5000m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_5000m[names(Indice_rarete_g1_g2_g3_b1_5000m)%in%g1_g2_g3_b1_5000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_5000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_5000m=sum(Indice_rarete_g1_g2_g3_b1_5000m_station)                #somme des indices de rareté  de ttes les sp


source('script.rarete.complementarite.r')

#Selection des dites selon un scenario de rareté complementarité basé sur les habitats

#habitat géo1
selection.rarete.complementarite.g1_g2_g3_b1.30m<-run.bin(habitat_par_site_g1_g2_g3_b1_30m_pres_abs,Indice_rarete_g1_g2_g3_b1_30m_station)
selection.rarete.complementarite.g1_g2_g3_b1.60m<-run.bin(habitat_par_site_g1_g2_g3_b1_60m_pres_abs,Indice_rarete_g1_g2_g3_b1_60m_station)
selection.rarete.complementarite.g1_g2_g3_b1.100m<-run.bin(habitat_par_site_g1_g2_g3_b1_100m_pres_abs,Indice_rarete_g1_g2_g3_b1_100m_station)
selection.rarete.complementarite.g1_g2_g3_b1.250m<-run.bin(habitat_par_site_g1_g2_g3_b1_250m_pres_abs,Indice_rarete_g1_g2_g3_b1_250m_station)
selection.rarete.complementarite.g1_g2_g3_b1.500m<-run.bin(habitat_par_site_g1_g2_g3_b1_500m_pres_abs,Indice_rarete_g1_g2_g3_b1_500m_station)
selection.rarete.complementarite.g1_g2_g3_b1.750m<-run.bin(habitat_par_site_g1_g2_g3_b1_750m_pres_abs,Indice_rarete_g1_g2_g3_b1_750m_station)
selection.rarete.complementarite.g1_g2_g3_b1.1000m<-run.bin(habitat_par_site_g1_g2_g3_b1_1000m_pres_abs,Indice_rarete_g1_g2_g3_b1_1000m_station)
selection.rarete.complementarite.g1_g2_g3_b1.3000m<-run.bin(habitat_par_site_g1_g2_g3_b1_3000m_pres_abs,Indice_rarete_g1_g2_g3_b1_3000m_station)
selection.rarete.complementarite.g1_g2_g3_b1.5000m<-run.bin(habitat_par_site_g1_g2_g3_b1_5000m_pres_abs,Indice_rarete_g1_g2_g3_b1_5000m_station)

#DIVERSITE SELECTIONNE PAR LE SCRIPT RARETE COMPLEMENTARITE#
 source('script.diversite.rarete.selectionnee.r')
 
 diversite.selection.rarete.complementarite.g1_g2_g3_b1.30m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1.30m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1.60m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1.60m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1.100m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1.100m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1.250m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1.250m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1.500m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1.500m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1.750m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1.750m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1.1000m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1.1000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1.3000m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1.3000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1.5000m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1.5000m,Indice_rarete_station)

#CALCUL DES SAI


source('script.sai.r')
SAI_rarete.complementarite.g1_g2_g3_b1.30m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1.30m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1.60m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1.60m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1.100m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1.100m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1.250m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1.250m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1.500m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1.500m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1.750m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1.750m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1.1000m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1.1000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1.3000m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1.3000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1.5000m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1.5000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)



 
############################HABITAT g1_g2_g3_b1_b2###########################################


 ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances


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



#Load de la requette "REQ_site_g1_g2_g3_b1_b2_30m.csv"

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

#calcul indice de rareté des habitats

#30m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_30m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_30m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_g1_g2_g3_b1_b2_30m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_30m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_30m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_g1_g2_g3_b1_b2_30m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_30m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_30m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_30m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_30m_station<-matrix(0,500,27)
Indice_rarete_g1_g2_g3_b1_b2_30m_station<-NULL

    for (m in 1:27){
  g1_g2_g3_b1_b2_30m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_30m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_30m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_30m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_30m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_30m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_30m[names(Indice_rarete_g1_g2_g3_b1_b2_30m)%in%g1_g2_g3_b1_b2_30m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_30m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_30m=sum(Indice_rarete_g1_g2_g3_b1_b2_30m_station)                #somme des indices de rareté  de ttes les sp


   
#60m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_60m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_60m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_g1_g2_g3_b1_b2_60m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_60m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_60m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_g1_g2_g3_b1_b2_60m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_60m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_60m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_60m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_60m_station<-matrix(0,500,27)
Indice_rarete_g1_g2_g3_b1_b2_60m_station<-NULL

    for (m in 1:27){
  g1_g2_g3_b1_b2_60m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_60m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_60m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_60m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_60m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_60m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_60m[names(Indice_rarete_g1_g2_g3_b1_b2_60m)%in%g1_g2_g3_b1_b2_60m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_60m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_60m=sum(Indice_rarete_g1_g2_g3_b1_b2_60m_station)                #somme des indices de rareté  de ttes les sp


    
#100m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_100m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_100m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_g1_g2_g3_b1_b2_100m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_100m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_100m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_g1_g2_g3_b1_b2_100m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_100m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_100m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_100m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_100m_station<-matrix(0,500,27)
Indice_rarete_g1_g2_g3_b1_b2_100m_station<-NULL

    for (m in 1:27){
  g1_g2_g3_b1_b2_100m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_100m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_100m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_100m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_100m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_100m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_100m[names(Indice_rarete_g1_g2_g3_b1_b2_100m)%in%g1_g2_g3_b1_b2_100m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_100m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_100m=sum(Indice_rarete_g1_g2_g3_b1_b2_100m_station)                #somme des indices de rareté  de ttes les sp

              
#250m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_250m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_250m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_g1_g2_g3_b1_b2_250m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_250m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_250m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_g1_g2_g3_b1_b2_250m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_250m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_250m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_250m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_250m_station<-matrix(0,500,27)
Indice_rarete_g1_g2_g3_b1_b2_250m_station<-NULL

    for (m in 1:27){
  g1_g2_g3_b1_b2_250m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_250m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_250m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_250m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_250m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_250m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_250m[names(Indice_rarete_g1_g2_g3_b1_b2_250m)%in%g1_g2_g3_b1_b2_250m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_250m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_250m=sum(Indice_rarete_g1_g2_g3_b1_b2_250m_station)                #somme des indices de rareté  de ttes les sp

     
#500m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_500m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_500m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_g1_g2_g3_b1_b2_500m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_500m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_500m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_g1_g2_g3_b1_b2_500m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_500m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_500m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_500m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_500m_station<-matrix(0,500,27)
Indice_rarete_g1_g2_g3_b1_b2_500m_station<-NULL

    for (m in 1:27){
  g1_g2_g3_b1_b2_500m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_500m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_500m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_500m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_500m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_500m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_500m[names(Indice_rarete_g1_g2_g3_b1_b2_500m)%in%g1_g2_g3_b1_b2_500m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_500m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_500m=sum(Indice_rarete_g1_g2_g3_b1_b2_500m_station)                #somme des indices de rareté  de ttes les sp

     
#750m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_750m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_750m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_g1_g2_g3_b1_b2_750m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_750m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_750m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_g1_g2_g3_b1_b2_750m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_750m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_750m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_750m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_750m_station<-matrix(0,500,27)
Indice_rarete_g1_g2_g3_b1_b2_750m_station<-NULL

    for (m in 1:27){
  g1_g2_g3_b1_b2_750m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_750m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_750m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_750m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_750m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_750m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_750m[names(Indice_rarete_g1_g2_g3_b1_b2_750m)%in%g1_g2_g3_b1_b2_750m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_750m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_750m=sum(Indice_rarete_g1_g2_g3_b1_b2_750m_station)                #somme des indices de rareté  de ttes les sp

      
#1000m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_1000m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_1000m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_g1_g2_g3_b1_b2_1000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_1000m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_1000m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_g1_g2_g3_b1_b2_1000m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_1000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_1000m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_1000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_1000m_station<-matrix(0,500,27)
Indice_rarete_g1_g2_g3_b1_b2_1000m_station<-NULL

    for (m in 1:27){
  g1_g2_g3_b1_b2_1000m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_1000m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_1000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_1000m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_1000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_1000m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_1000m[names(Indice_rarete_g1_g2_g3_b1_b2_1000m)%in%g1_g2_g3_b1_b2_1000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_1000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_1000m=sum(Indice_rarete_g1_g2_g3_b1_b2_1000m_station)                #somme des indices de rareté  de ttes les sp

                  
#3000m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_3000m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_3000m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_g1_g2_g3_b1_b2_3000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_3000m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_3000m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_g1_g2_g3_b1_b2_3000m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_3000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_3000m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_3000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_3000m_station<-matrix(0,500,27)
Indice_rarete_g1_g2_g3_b1_b2_3000m_station<-NULL

    for (m in 1:27){
  g1_g2_g3_b1_b2_3000m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_3000m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_3000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_3000m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_3000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_3000m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_3000m[names(Indice_rarete_g1_g2_g3_b1_b2_3000m)%in%g1_g2_g3_b1_b2_3000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_3000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_3000m=sum(Indice_rarete_g1_g2_g3_b1_b2_3000m_station)                #somme des indices de rareté  de ttes les sp

         
#5000m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_5000m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_5000m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_g1_g2_g3_b1_b2_5000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_5000m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_5000m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_g1_g2_g3_b1_b2_5000m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_5000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_5000m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_5000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_5000m_station<-matrix(0,500,27)
Indice_rarete_g1_g2_g3_b1_b2_5000m_station<-NULL

    for (m in 1:27){
  g1_g2_g3_b1_b2_5000m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_5000m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_5000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_5000m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_5000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_5000m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_5000m[names(Indice_rarete_g1_g2_g3_b1_b2_5000m)%in%g1_g2_g3_b1_b2_5000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_5000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_5000m=sum(Indice_rarete_g1_g2_g3_b1_b2_5000m_station)                #somme des indices de rareté  de ttes les sp


source('script.rarete.complementarite.r')

#Selection des dites selon un scenario de rareté complementarité basé sur les habitats

#habitat géo1
selection.rarete.complementarite.g1_g2_g3_b1_b2.30m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_30m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_30m_station)
selection.rarete.complementarite.g1_g2_g3_b1_b2.60m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_60m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_60m_station)
selection.rarete.complementarite.g1_g2_g3_b1_b2.100m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_100m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_100m_station)
selection.rarete.complementarite.g1_g2_g3_b1_b2.250m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_250m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_250m_station)
selection.rarete.complementarite.g1_g2_g3_b1_b2.500m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_500m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_500m_station)
selection.rarete.complementarite.g1_g2_g3_b1_b2.750m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_750m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_750m_station)
selection.rarete.complementarite.g1_g2_g3_b1_b2.1000m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_1000m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_1000m_station)
selection.rarete.complementarite.g1_g2_g3_b1_b2.3000m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_3000m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_3000m_station)
selection.rarete.complementarite.g1_g2_g3_b1_b2.5000m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_5000m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_5000m_station)

#DIVERSITE SELECTIONNE PAR LE SCRIPT RARETE COMPLEMENTARITE#
 source('script.diversite.rarete.selectionnee.r')
 
 diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.30m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2.30m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.60m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2.60m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.100m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2.100m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.250m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2.250m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.500m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2.500m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.750m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2.750m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.1000m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2.1000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.3000m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2.3000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.5000m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2.5000m,Indice_rarete_station)

#CALCUL DES SAI


source('script.sai.r')
SAI_rarete.complementarite.g1_g2_g3_b1_b2.30m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.30m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1_b2.60m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.60m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1_b2.100m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.100m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1_b2.250m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.250m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1_b2.500m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.500m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1_b2.750m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.750m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1_b2.1000m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.1000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1_b2.3000m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.3000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1_b2.5000m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.5000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)



 
############################HABITAT g1_g2_g3_b1_b2_b3###########################################


 ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances


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



#Load de la requette "REQ_site_g1_g2_g3_b1_b2_b3_30m.csv"

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

#calcul indice de rareté des habitats

#30m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_b3_30m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_b3_30m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_g1_g2_g3_b1_b2_b3_30m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_30m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_30m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_g1_g2_g3_b1_b2_b3_30m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_30m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_b3_30m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_b3_30m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_b3_30m_station<-matrix(0,500,27)
Indice_rarete_g1_g2_g3_b1_b2_b3_30m_station<-NULL

    for (m in 1:27){
  g1_g2_g3_b1_b2_b3_30m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_30m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_30m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_b3_30m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_30m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_b3_30m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_30m[names(Indice_rarete_g1_g2_g3_b1_b2_b3_30m)%in%g1_g2_g3_b1_b2_b3_30m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_b3_30m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_b3_30m=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_30m_station)                #somme des indices de rareté  de ttes les sp


   
#60m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_b3_60m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_b3_60m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_g1_g2_g3_b1_b2_b3_60m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_60m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_60m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_g1_g2_g3_b1_b2_b3_60m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_60m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_b3_60m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_b3_60m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_b3_60m_station<-matrix(0,500,27)
Indice_rarete_g1_g2_g3_b1_b2_b3_60m_station<-NULL

    for (m in 1:27){
  g1_g2_g3_b1_b2_b3_60m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_60m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_60m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_b3_60m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_60m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_b3_60m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_60m[names(Indice_rarete_g1_g2_g3_b1_b2_b3_60m)%in%g1_g2_g3_b1_b2_b3_60m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_b3_60m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_b3_60m=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_60m_station)                #somme des indices de rareté  de ttes les sp


    
#100m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_b3_100m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_b3_100m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_g1_g2_g3_b1_b2_b3_100m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_100m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_100m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_g1_g2_g3_b1_b2_b3_100m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_100m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_b3_100m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_b3_100m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_b3_100m_station<-matrix(0,500,27)
Indice_rarete_g1_g2_g3_b1_b2_b3_100m_station<-NULL

    for (m in 1:27){
  g1_g2_g3_b1_b2_b3_100m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_100m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_100m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_b3_100m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_100m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_b3_100m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_100m[names(Indice_rarete_g1_g2_g3_b1_b2_b3_100m)%in%g1_g2_g3_b1_b2_b3_100m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_b3_100m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_b3_100m=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_100m_station)                #somme des indices de rareté  de ttes les sp

              
#250m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_b3_250m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_b3_250m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_g1_g2_g3_b1_b2_b3_250m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_250m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_250m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_g1_g2_g3_b1_b2_b3_250m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_250m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_b3_250m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_b3_250m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_b3_250m_station<-matrix(0,500,27)
Indice_rarete_g1_g2_g3_b1_b2_b3_250m_station<-NULL

    for (m in 1:27){
  g1_g2_g3_b1_b2_b3_250m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_250m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_250m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_b3_250m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_250m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_b3_250m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_250m[names(Indice_rarete_g1_g2_g3_b1_b2_b3_250m)%in%g1_g2_g3_b1_b2_b3_250m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_b3_250m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_b3_250m=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_250m_station)                #somme des indices de rareté  de ttes les sp

     
#500m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_b3_500m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_b3_500m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_g1_g2_g3_b1_b2_b3_500m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_500m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_500m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_g1_g2_g3_b1_b2_b3_500m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_500m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_b3_500m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_b3_500m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_b3_500m_station<-matrix(0,500,27)
Indice_rarete_g1_g2_g3_b1_b2_b3_500m_station<-NULL

    for (m in 1:27){
  g1_g2_g3_b1_b2_b3_500m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_500m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_500m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_b3_500m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_500m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_b3_500m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_500m[names(Indice_rarete_g1_g2_g3_b1_b2_b3_500m)%in%g1_g2_g3_b1_b2_b3_500m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_b3_500m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_b3_500m=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_500m_station)                #somme des indices de rareté  de ttes les sp

     
#750m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_b3_750m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_b3_750m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_g1_g2_g3_b1_b2_b3_750m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_750m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_750m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_g1_g2_g3_b1_b2_b3_750m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_750m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_b3_750m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_b3_750m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_b3_750m_station<-matrix(0,500,27)
Indice_rarete_g1_g2_g3_b1_b2_b3_750m_station<-NULL

    for (m in 1:27){
  g1_g2_g3_b1_b2_b3_750m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_750m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_750m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_b3_750m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_750m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_b3_750m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_750m[names(Indice_rarete_g1_g2_g3_b1_b2_b3_750m)%in%g1_g2_g3_b1_b2_b3_750m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_b3_750m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_b3_750m=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_750m_station)                #somme des indices de rareté  de ttes les sp

      
#1000m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_b3_1000m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_b3_1000m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_g1_g2_g3_b1_b2_b3_1000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_1000m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_1000m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_g1_g2_g3_b1_b2_b3_1000m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_1000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_b3_1000m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_b3_1000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_b3_1000m_station<-matrix(0,500,27)
Indice_rarete_g1_g2_g3_b1_b2_b3_1000m_station<-NULL

    for (m in 1:27){
  g1_g2_g3_b1_b2_b3_1000m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_1000m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_1000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_b3_1000m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_1000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_b3_1000m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_1000m[names(Indice_rarete_g1_g2_g3_b1_b2_b3_1000m)%in%g1_g2_g3_b1_b2_b3_1000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_b3_1000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_b3_1000m=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_1000m_station)                #somme des indices de rareté  de ttes les sp

                  
#3000m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_b3_3000m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_b3_3000m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_g1_g2_g3_b1_b2_b3_3000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_3000m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_3000m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_g1_g2_g3_b1_b2_b3_3000m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_3000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_b3_3000m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_b3_3000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_b3_3000m_station<-matrix(0,500,27)
Indice_rarete_g1_g2_g3_b1_b2_b3_3000m_station<-NULL

    for (m in 1:27){
  g1_g2_g3_b1_b2_b3_3000m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_3000m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_3000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_b3_3000m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_3000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_b3_3000m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_3000m[names(Indice_rarete_g1_g2_g3_b1_b2_b3_3000m)%in%g1_g2_g3_b1_b2_b3_3000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_b3_3000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_b3_3000m=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_3000m_station)                #somme des indices de rareté  de ttes les sp

         
#5000m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_b3_5000m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_b3_5000m=NULL
  site_habitat<-matrix(0,27,length(habitat_par_site_g1_g2_g3_b1_b2_b3_5000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_5000m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_5000m$Code_hab==b[n]]),rep(0,27-length(habitat_par_site_g1_g2_g3_b1_b2_b3_5000m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_5000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_b3_5000m[n]<-2^(27-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_b3_5000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_b3_5000m_station<-matrix(0,500,27)
Indice_rarete_g1_g2_g3_b1_b2_b3_5000m_station<-NULL

    for (m in 1:27){
  g1_g2_g3_b1_b2_b3_5000m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_5000m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_5000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_b3_5000m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_5000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_b3_5000m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_5000m[names(Indice_rarete_g1_g2_g3_b1_b2_b3_5000m)%in%g1_g2_g3_b1_b2_b3_5000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_b3_5000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_b3_5000m=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_5000m_station)                #somme des indices de rareté  de ttes les sp


source('script.rarete.complementarite.r')

#Selection des dites selon un scenario de rareté complementarité basé sur les habitats

#habitat géo1
selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.30m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_b3_30m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_b3_30m_station)
selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.60m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_b3_60m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_b3_60m_station)
selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.100m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_b3_100m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_b3_100m_station)
selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.250m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_b3_250m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_b3_250m_station)
selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.500m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_b3_500m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_b3_500m_station)
selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.750m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_b3_750m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_b3_750m_station)
selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.1000m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_b3_1000m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_b3_1000m_station)
selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.3000m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_b3_3000m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_b3_3000m_station)
selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.5000m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_b3_5000m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_b3_5000m_station)

#DIVERSITE SELECTIONNE PAR LE SCRIPT RARETE COMPLEMENTARITE#
 source('script.diversite.rarete.selectionnee.r')
 
 diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.30m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.30m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.60m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.60m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.100m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.100m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.250m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.250m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.500m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.500m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.750m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.750m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.1000m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.1000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.3000m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.3000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.5000m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.5000m,Indice_rarete_station)

#CALCUL DES SAI


source('script.sai.r')
SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.30m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.30m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.60m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.60m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.100m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.100m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.250m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.250m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.500m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.500m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.750m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.750m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.1000m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.1000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.3000m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.3000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.5000m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.5000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)


  #STOCKAGE DES RESULTATS

targeted.station.results[,1]=c(SAI_rarete.complementarite.geo1.30m,SAI_rarete.complementarite.geo1.60m,SAI_rarete.complementarite.geo1.100m,SAI_rarete.complementarite.geo1.250m,
SAI_rarete.complementarite.geo1.500m,SAI_rarete.complementarite.geo1.750m,SAI_rarete.complementarite.geo1.1000m,SAI_rarete.complementarite.geo1.3000m,SAI_rarete.complementarite.geo1.5000m,
SAI_rarete.complementarite.geo2.30m,SAI_rarete.complementarite.geo2.60m,SAI_rarete.complementarite.geo2.100m,SAI_rarete.complementarite.geo2.250m,
SAI_rarete.complementarite.geo2.500m,SAI_rarete.complementarite.geo2.750m,SAI_rarete.complementarite.geo2.1000m,SAI_rarete.complementarite.geo2.3000m,SAI_rarete.complementarite.geo2.5000m,
SAI_rarete.complementarite.geo3.100m,SAI_rarete.complementarite.geo3.250m,
SAI_rarete.complementarite.geo3.500m,SAI_rarete.complementarite.geo3.750m,SAI_rarete.complementarite.geo3.1000m,SAI_rarete.complementarite.geo3.3000m,SAI_rarete.complementarite.geo3.5000m,
SAI_rarete.complementarite.bent1.30m,SAI_rarete.complementarite.bent1.60m,SAI_rarete.complementarite.bent1.100m,SAI_rarete.complementarite.bent1.250m,
SAI_rarete.complementarite.bent1.500m,SAI_rarete.complementarite.bent1.750m,SAI_rarete.complementarite.bent1.1000m,SAI_rarete.complementarite.bent1.3000m,SAI_rarete.complementarite.bent1.5000m,
SAI_rarete.complementarite.bent2.30m,SAI_rarete.complementarite.bent2.60m,SAI_rarete.complementarite.bent2.100m,SAI_rarete.complementarite.bent2.250m,
SAI_rarete.complementarite.bent2.500m,SAI_rarete.complementarite.bent2.750m,SAI_rarete.complementarite.bent2.1000m,SAI_rarete.complementarite.bent2.3000m,SAI_rarete.complementarite.bent2.5000m,
SAI_rarete.complementarite.bent3.30m,SAI_rarete.complementarite.bent3.60m,SAI_rarete.complementarite.bent3.100m,SAI_rarete.complementarite.bent3.250m,
SAI_rarete.complementarite.bent3.500m,SAI_rarete.complementarite.bent3.750m,SAI_rarete.complementarite.bent3.1000m,SAI_rarete.complementarite.bent3.3000m,SAI_rarete.complementarite.bent3.5000m,
SAI_rarete.complementarite.geo1_geo2.30m,SAI_rarete.complementarite.geo1_geo2.60m,SAI_rarete.complementarite.geo1_geo2.100m,SAI_rarete.complementarite.geo1_geo2.250m,
SAI_rarete.complementarite.geo1_geo2.500m,SAI_rarete.complementarite.geo1_geo2.750m,SAI_rarete.complementarite.geo1_geo2.1000m,SAI_rarete.complementarite.geo1_geo2.3000m,SAI_rarete.complementarite.geo1_geo2.5000m,
SAI_rarete.complementarite.geo1_geo2_geo3.30m,SAI_rarete.complementarite.geo1_geo2_geo3.60m,SAI_rarete.complementarite.geo1_geo2_geo3.100m,SAI_rarete.complementarite.geo1_geo2_geo3.250m,
SAI_rarete.complementarite.geo1_geo2_geo3.500m,SAI_rarete.complementarite.geo1_geo2_geo3.750m,SAI_rarete.complementarite.geo1_geo2_geo3.1000m,SAI_rarete.complementarite.geo1_geo2_geo3.3000m,SAI_rarete.complementarite.geo1_geo2_geo3.5000m,
SAI_rarete.complementarite.g1_g2_g3_b1.30m,SAI_rarete.complementarite.g1_g2_g3_b1.60m,SAI_rarete.complementarite.g1_g2_g3_b1.100m,SAI_rarete.complementarite.g1_g2_g3_b1.250m,
SAI_rarete.complementarite.g1_g2_g3_b1.500m,SAI_rarete.complementarite.g1_g2_g3_b1.750m,SAI_rarete.complementarite.g1_g2_g3_b1.1000m,SAI_rarete.complementarite.g1_g2_g3_b1.3000m,SAI_rarete.complementarite.g1_g2_g3_b1.5000m,
SAI_rarete.complementarite.g1_g2_g3_b1_b2.30m,SAI_rarete.complementarite.g1_g2_g3_b1_b2.60m,SAI_rarete.complementarite.g1_g2_g3_b1_b2.100m,SAI_rarete.complementarite.g1_g2_g3_b1_b2.250m,
SAI_rarete.complementarite.g1_g2_g3_b1_b2.500m,SAI_rarete.complementarite.g1_g2_g3_b1_b2.750m,SAI_rarete.complementarite.g1_g2_g3_b1_b2.1000m,SAI_rarete.complementarite.g1_g2_g3_b1_b2.3000m,SAI_rarete.complementarite.g1_g2_g3_b1_b2.5000m,
SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.30m,SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.60m,SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.100m,SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.250m,
SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.500m,SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.750m,SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.1000m,SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.3000m,SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.5000m)




 #######################        -1 STATIONS     #############################

  targeted.stations<-c('T_S249','T_S250','T_S251','T_S252','T_S253','T_S254','T_S255','T_S256','T_S257','T_S258','T_S259','T_S260','T_S261','T_S262','T_S263','T_S264','T_S265','T_S266','T_S267','T_S268','T_S269','T_S270','T_S271','T_S272','T_S273','T_S274','T_S275')

for (p in 1:length(targeted.stations)){
  
suppressed.station<-targeted.stations[p]



 #Load du tableau : "site" "poissons"

repertoire<-"C:/Documents and Settings/vanw/Bureau/R"
setwd(repertoire)
input<-"INPUT/"
output<-"OUTPUT/"
setwd(input)
recensement=read.csv2("REQ_site_poisson.IRD.csv")

#j'enlève la station supprimée
recensement<-recensement[!recensement$Code_Site%in%suppressed.station,]
recensement$Code_Site<-factor(recensement$Code_Site)

#J'enlève les doublons
recensement_unique=unique(recensement)
  recensement_unique$Code_Site<-factor(recensement_unique$Code_Site)
  
#création d'un tableau de présence absence
tablo_pres_abs<-table(recensement_unique$Code_Site,recensement_unique$Code_poisson)
 recensement$Code_Site<-factor(recensement$Code_Site)

#calcul du nombre d'espèce par station
vec_poisson<-apply(tablo_pres_abs,1,sum)

#selection aléatoire d'un pool de site
setwd(repertoire)
setwd(input)
source('script.aleatoire.r')
selection.aleatoire<-matrix(0,dim(tablo_pres_abs)[1],999) #999
for (i in 1:999){ #999
selection.aleatoire[,i]<-run.alea(tablo_pres_abs)}


#calcul indice de rareté des espèces

  b=unique(recensement_unique$Code_poisson)
 
  Indice_rarete_sp=NULL
  site_poisson<-matrix(0,26,length(recensement_unique$Code_poisson))

  for (n in 1:length(b)){
  site_poisson[,n]<-c(as.vector(recensement_unique$Code_Site[recensement_unique$Code_poisson==b[n]]),rep(0,26-length(recensement_unique$Code_Site[recensement_unique$Code_poisson==b[n]])))
  names(site_poisson[,n])=b[n]
  Indice_rarete_sp[n]<-2^(26-length(site_poisson[,n][site_poisson[,n]!=0]))
  names(Indice_rarete_sp)[n]=as.vector(b[n])     #Indice_rarete_sp est un vecteur qui contien l'indice de rareté pour chaque espèce
  } 
    
  #calcul indice de rareté des stations
  
  sites.etudies=targeted.stations[targeted.stations!=suppressed.station]
  
  
    poissons_dans_station<-matrix(0,335,26)
    Indice_rarete_station<-NULL
    recensement_unique$Code_Site
   
   for (m in 1:26){             
    poissons_dans_station[,m]<-c(as.vector(recensement_unique$Code_poisson[recensement_unique$Code_Site==sites.etudies[m]]),rep(0,335-length(recensement_unique$Code_poisson[recensement_unique$Code_Site==sites.etudies[m]])))
  Indice_rarete_station[m]=sum(Indice_rarete_sp[names(Indice_rarete_sp)%in%poissons_dans_station[,m]])            
  }              
                           
 
  names(Indice_rarete_station)<-sites.etudies

   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot=sum(Indice_rarete_station)

#Selection des sites selon la rareté de poissons
source('script.rarete.complementarite.r')
selection.rarete.complementarite.poisson<-run.bin(tablo_pres_abs,Indice_rarete_station)


#rareté de poisson incluse dans les sites selectionnés (aléatoirement et selon diversité de poisson)
source('script.diversite.rarete.selectionnee.r')
 diversite.rarete.selection.aleatoire<-matrix(0,dim(tablo_pres_abs)[1],999)  #999
 diversite.rarete.selection.aleatoire<-apply(as.matrix(selection.aleatoire),2,function(x) species.rarete.evol(x,Indice_rarete_station))
  diversite.rarete.moyenne.selection.aleatoire<-apply(diversite.rarete.selection.aleatoire,1,mean)


diversite.selection.rarete.complementarite.poisson<-species.rarete.evol(selection.rarete.complementarite.poisson,Indice_rarete_station)


# Calcul des courbes aléatoires (moyenne, quantile supérieur à 97.5, quantile inférieur à 2.5% sur les 1000 tirages)

diversite.rarete.moyenne.selection.aleatoire
diversite.rarete.borne.sup.selection.aleatoire=apply(diversite.rarete.selection.aleatoire,1,function(x) sort(x)[ceiling(0.975*length(x))])
diversite.rarete.borne.inf.selection.aleatoire=apply(diversite.rarete.selection.aleatoire,1,function(x) sort(x)[ceiling(0.025*length(x))])


#nombre total d'espèces de poissons observés
nombre.poissons.total <- nlevels(factor(recensement_unique$Code_poisson))

# fonction pourcentage
percentage=function(total,truc){
  truc * 100 / total
  }

 # Pourcentages, RANDOM, upper, low
diversite.rarete.moyenne.selection.aleatoire.pourcentage <- c(0,percentage(tot,diversite.rarete.moyenne.selection.aleatoire))
diversite.rarete.borne.sup.selection.aleatoire.pourcentage <- c(0,percentage(tot,diversite.rarete.borne.sup.selection.aleatoire))
diversite.rarete.borne.inf.selection.aleatoire.pourcentage <- c(0,percentage(tot,diversite.rarete.borne.inf.selection.aleatoire))

 diversite.selection.rarete.complementarite.poisson.pourcentage<- c(0,percentage(tot,diversite.selection.rarete.complementarite.poisson))



 
 
############################HABITAT geo1###########################################


 ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances


relation.habitats.sites.geo1<-read.csv2("fc3_intersect_geo1.csv")


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



#Load de la requette "REQ_site_geo1_30m.csv"

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

#calcul indice de rareté des habitats

#30m

  b=unique(habitat_par_site_geo1_30m$Code_hab)

    Indice_rarete_geo1_30m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo1_30m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_30m$Code_Site[habitat_par_site_geo1_30m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo1_30m$Code_Site[habitat_par_site_geo1_30m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_30m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_30m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_30m_station<-matrix(0,500,26)
Indice_rarete_geo1_30m_station<-NULL

    for (m in 1:26){
  geo1_30m_station[,m]<-c(as.vector(habitat_par_site_geo1_30m$Code_hab[habitat_par_site_geo1_30m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_30m$Code_hab[habitat_par_site_geo1_30m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_30m_station[m]=sum(Indice_rarete_geo1_30m[names(Indice_rarete_geo1_30m)%in%geo1_30m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_30m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_30m=sum(Indice_rarete_geo1_30m_station)                #somme des indices de rareté  de ttes les sp


   
#60m

  b=unique(habitat_par_site_geo1_60m$Code_hab)

    Indice_rarete_geo1_60m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo1_60m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_60m$Code_Site[habitat_par_site_geo1_60m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo1_60m$Code_Site[habitat_par_site_geo1_60m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_60m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_60m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_60m_station<-matrix(0,500,26)
Indice_rarete_geo1_60m_station<-NULL

    for (m in 1:26){
  geo1_60m_station[,m]<-c(as.vector(habitat_par_site_geo1_60m$Code_hab[habitat_par_site_geo1_60m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_60m$Code_hab[habitat_par_site_geo1_60m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_60m_station[m]=sum(Indice_rarete_geo1_60m[names(Indice_rarete_geo1_60m)%in%geo1_60m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_60m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_60m=sum(Indice_rarete_geo1_60m_station)                #somme des indices de rareté  de ttes les sp


    
#100m

  b=unique(habitat_par_site_geo1_100m$Code_hab)

    Indice_rarete_geo1_100m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo1_100m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_100m$Code_Site[habitat_par_site_geo1_100m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo1_100m$Code_Site[habitat_par_site_geo1_100m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_100m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_100m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_100m_station<-matrix(0,500,26)
Indice_rarete_geo1_100m_station<-NULL

    for (m in 1:26){
  geo1_100m_station[,m]<-c(as.vector(habitat_par_site_geo1_100m$Code_hab[habitat_par_site_geo1_100m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_100m$Code_hab[habitat_par_site_geo1_100m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_100m_station[m]=sum(Indice_rarete_geo1_100m[names(Indice_rarete_geo1_100m)%in%geo1_100m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_100m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_100m=sum(Indice_rarete_geo1_100m_station)                #somme des indices de rareté  de ttes les sp

              
#250m

  b=unique(habitat_par_site_geo1_250m$Code_hab)

    Indice_rarete_geo1_250m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo1_250m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_250m$Code_Site[habitat_par_site_geo1_250m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo1_250m$Code_Site[habitat_par_site_geo1_250m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_250m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_250m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_250m_station<-matrix(0,500,26)
Indice_rarete_geo1_250m_station<-NULL

    for (m in 1:26){
  geo1_250m_station[,m]<-c(as.vector(habitat_par_site_geo1_250m$Code_hab[habitat_par_site_geo1_250m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_250m$Code_hab[habitat_par_site_geo1_250m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_250m_station[m]=sum(Indice_rarete_geo1_250m[names(Indice_rarete_geo1_250m)%in%geo1_250m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_250m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_250m=sum(Indice_rarete_geo1_250m_station)                #somme des indices de rareté  de ttes les sp

     
#500m

  b=unique(habitat_par_site_geo1_500m$Code_hab)

    Indice_rarete_geo1_500m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo1_500m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_500m$Code_Site[habitat_par_site_geo1_500m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo1_500m$Code_Site[habitat_par_site_geo1_500m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_500m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_500m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_500m_station<-matrix(0,500,26)
Indice_rarete_geo1_500m_station<-NULL

    for (m in 1:26){
  geo1_500m_station[,m]<-c(as.vector(habitat_par_site_geo1_500m$Code_hab[habitat_par_site_geo1_500m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_500m$Code_hab[habitat_par_site_geo1_500m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_500m_station[m]=sum(Indice_rarete_geo1_500m[names(Indice_rarete_geo1_500m)%in%geo1_500m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_500m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_500m=sum(Indice_rarete_geo1_500m_station)                #somme des indices de rareté  de ttes les sp

     
#750m

  b=unique(habitat_par_site_geo1_750m$Code_hab)

    Indice_rarete_geo1_750m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo1_750m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_750m$Code_Site[habitat_par_site_geo1_750m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo1_750m$Code_Site[habitat_par_site_geo1_750m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_750m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_750m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_750m_station<-matrix(0,500,26)
Indice_rarete_geo1_750m_station<-NULL

    for (m in 1:26){
  geo1_750m_station[,m]<-c(as.vector(habitat_par_site_geo1_750m$Code_hab[habitat_par_site_geo1_750m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_750m$Code_hab[habitat_par_site_geo1_750m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_750m_station[m]=sum(Indice_rarete_geo1_750m[names(Indice_rarete_geo1_750m)%in%geo1_750m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_750m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_750m=sum(Indice_rarete_geo1_750m_station)                #somme des indices de rareté  de ttes les sp

      
#1000m

  b=unique(habitat_par_site_geo1_1000m$Code_hab)

    Indice_rarete_geo1_1000m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo1_1000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_1000m$Code_Site[habitat_par_site_geo1_1000m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo1_1000m$Code_Site[habitat_par_site_geo1_1000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_1000m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_1000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_1000m_station<-matrix(0,500,26)
Indice_rarete_geo1_1000m_station<-NULL

    for (m in 1:26){
  geo1_1000m_station[,m]<-c(as.vector(habitat_par_site_geo1_1000m$Code_hab[habitat_par_site_geo1_1000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_1000m$Code_hab[habitat_par_site_geo1_1000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_1000m_station[m]=sum(Indice_rarete_geo1_1000m[names(Indice_rarete_geo1_1000m)%in%geo1_1000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_1000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_1000m=sum(Indice_rarete_geo1_1000m_station)                #somme des indices de rareté  de ttes les sp

                  
#3000m

  b=unique(habitat_par_site_geo1_3000m$Code_hab)

    Indice_rarete_geo1_3000m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo1_3000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_3000m$Code_Site[habitat_par_site_geo1_3000m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo1_3000m$Code_Site[habitat_par_site_geo1_3000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_3000m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_3000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_3000m_station<-matrix(0,500,26)
Indice_rarete_geo1_3000m_station<-NULL

    for (m in 1:26){
  geo1_3000m_station[,m]<-c(as.vector(habitat_par_site_geo1_3000m$Code_hab[habitat_par_site_geo1_3000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_3000m$Code_hab[habitat_par_site_geo1_3000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_3000m_station[m]=sum(Indice_rarete_geo1_3000m[names(Indice_rarete_geo1_3000m)%in%geo1_3000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_3000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_3000m=sum(Indice_rarete_geo1_3000m_station)                #somme des indices de rareté  de ttes les sp

         
#5000m

  b=unique(habitat_par_site_geo1_5000m$Code_hab)

    Indice_rarete_geo1_5000m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo1_5000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_5000m$Code_Site[habitat_par_site_geo1_5000m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo1_5000m$Code_Site[habitat_par_site_geo1_5000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_5000m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_5000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_5000m_station<-matrix(0,500,26)
Indice_rarete_geo1_5000m_station<-NULL

    for (m in 1:26){
  geo1_5000m_station[,m]<-c(as.vector(habitat_par_site_geo1_5000m$Code_hab[habitat_par_site_geo1_5000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_5000m$Code_hab[habitat_par_site_geo1_5000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_5000m_station[m]=sum(Indice_rarete_geo1_5000m[names(Indice_rarete_geo1_5000m)%in%geo1_5000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_5000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_5000m=sum(Indice_rarete_geo1_5000m_station)                #somme des indices de rareté  de ttes les sp


source('script.rarete.complementarite.r')

#Selection des dites selon un scenario de rareté complementarité basé sur les habitats

#habitat géo1
selection.rarete.complementarite.geo1.30m<-run.bin(habitat_par_site_geo1_30m_pres_abs,Indice_rarete_geo1_30m_station)
selection.rarete.complementarite.geo1.60m<-run.bin(habitat_par_site_geo1_60m_pres_abs,Indice_rarete_geo1_60m_station)
selection.rarete.complementarite.geo1.100m<-run.bin(habitat_par_site_geo1_100m_pres_abs,Indice_rarete_geo1_100m_station)
selection.rarete.complementarite.geo1.250m<-run.bin(habitat_par_site_geo1_250m_pres_abs,Indice_rarete_geo1_250m_station)
selection.rarete.complementarite.geo1.500m<-run.bin(habitat_par_site_geo1_500m_pres_abs,Indice_rarete_geo1_500m_station)
selection.rarete.complementarite.geo1.750m<-run.bin(habitat_par_site_geo1_750m_pres_abs,Indice_rarete_geo1_750m_station)
selection.rarete.complementarite.geo1.1000m<-run.bin(habitat_par_site_geo1_1000m_pres_abs,Indice_rarete_geo1_1000m_station)
selection.rarete.complementarite.geo1.3000m<-run.bin(habitat_par_site_geo1_3000m_pres_abs,Indice_rarete_geo1_3000m_station)
selection.rarete.complementarite.geo1.5000m<-run.bin(habitat_par_site_geo1_5000m_pres_abs,Indice_rarete_geo1_5000m_station)

#DIVERSITE SELECTIONNE PAR LE SCRIPT RARETE COMPLEMENTARITE#
 source('script.diversite.rarete.selectionnee.r')
 
 diversite.selection.rarete.complementarite.geo1.30m<-species.rarete.evol(selection.rarete.complementarite.geo1.30m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1.60m<-species.rarete.evol(selection.rarete.complementarite.geo1.60m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1.100m<-species.rarete.evol(selection.rarete.complementarite.geo1.100m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1.250m<-species.rarete.evol(selection.rarete.complementarite.geo1.250m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1.500m<-species.rarete.evol(selection.rarete.complementarite.geo1.500m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1.750m<-species.rarete.evol(selection.rarete.complementarite.geo1.750m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1.1000m<-species.rarete.evol(selection.rarete.complementarite.geo1.1000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1.3000m<-species.rarete.evol(selection.rarete.complementarite.geo1.3000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1.5000m<-species.rarete.evol(selection.rarete.complementarite.geo1.5000m,Indice_rarete_station)

#CALCUL DES SAI


source('script.sai.r')
SAI_rarete.complementarite.geo1.30m<-SAI(diversite.selection.rarete.complementarite.geo1.30m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1.60m<-SAI(diversite.selection.rarete.complementarite.geo1.60m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1.100m<-SAI(diversite.selection.rarete.complementarite.geo1.100m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1.250m<-SAI(diversite.selection.rarete.complementarite.geo1.250m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1.500m<-SAI(diversite.selection.rarete.complementarite.geo1.500m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1.750m<-SAI(diversite.selection.rarete.complementarite.geo1.750m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1.1000m<-SAI(diversite.selection.rarete.complementarite.geo1.1000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1.3000m<-SAI(diversite.selection.rarete.complementarite.geo1.3000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1.5000m<-SAI(diversite.selection.rarete.complementarite.geo1.5000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)



 
############################HABITAT geo2###########################################


 ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances


relation.habitats.sites.geo2<-read.csv2("fc3_intersect_geo2.csv")


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



#Load de la requette "REQ_site_geo2_30m.csv"

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

#calcul indice de rareté des habitats

#30m

  b=unique(habitat_par_site_geo2_30m$Code_hab)

    Indice_rarete_geo2_30m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo2_30m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo2_30m$Code_Site[habitat_par_site_geo2_30m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo2_30m$Code_Site[habitat_par_site_geo2_30m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo2_30m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo2_30m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo2_30m_station<-matrix(0,500,26)
Indice_rarete_geo2_30m_station<-NULL

    for (m in 1:26){
  geo2_30m_station[,m]<-c(as.vector(habitat_par_site_geo2_30m$Code_hab[habitat_par_site_geo2_30m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo2_30m$Code_hab[habitat_par_site_geo2_30m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo2_30m_station[m]=sum(Indice_rarete_geo2_30m[names(Indice_rarete_geo2_30m)%in%geo2_30m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo2_30m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo2_30m=sum(Indice_rarete_geo2_30m_station)                #somme des indices de rareté  de ttes les sp


   
#60m

  b=unique(habitat_par_site_geo2_60m$Code_hab)

    Indice_rarete_geo2_60m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo2_60m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo2_60m$Code_Site[habitat_par_site_geo2_60m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo2_60m$Code_Site[habitat_par_site_geo2_60m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo2_60m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo2_60m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo2_60m_station<-matrix(0,500,26)
Indice_rarete_geo2_60m_station<-NULL

    for (m in 1:26){
  geo2_60m_station[,m]<-c(as.vector(habitat_par_site_geo2_60m$Code_hab[habitat_par_site_geo2_60m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo2_60m$Code_hab[habitat_par_site_geo2_60m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo2_60m_station[m]=sum(Indice_rarete_geo2_60m[names(Indice_rarete_geo2_60m)%in%geo2_60m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo2_60m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo2_60m=sum(Indice_rarete_geo2_60m_station)                #somme des indices de rareté  de ttes les sp


    
#100m

  b=unique(habitat_par_site_geo2_100m$Code_hab)

    Indice_rarete_geo2_100m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo2_100m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo2_100m$Code_Site[habitat_par_site_geo2_100m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo2_100m$Code_Site[habitat_par_site_geo2_100m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo2_100m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo2_100m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo2_100m_station<-matrix(0,500,26)
Indice_rarete_geo2_100m_station<-NULL

    for (m in 1:26){
  geo2_100m_station[,m]<-c(as.vector(habitat_par_site_geo2_100m$Code_hab[habitat_par_site_geo2_100m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo2_100m$Code_hab[habitat_par_site_geo2_100m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo2_100m_station[m]=sum(Indice_rarete_geo2_100m[names(Indice_rarete_geo2_100m)%in%geo2_100m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo2_100m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo2_100m=sum(Indice_rarete_geo2_100m_station)                #somme des indices de rareté  de ttes les sp

              
#250m

  b=unique(habitat_par_site_geo2_250m$Code_hab)

    Indice_rarete_geo2_250m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo2_250m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo2_250m$Code_Site[habitat_par_site_geo2_250m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo2_250m$Code_Site[habitat_par_site_geo2_250m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo2_250m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo2_250m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo2_250m_station<-matrix(0,500,26)
Indice_rarete_geo2_250m_station<-NULL

    for (m in 1:26){
  geo2_250m_station[,m]<-c(as.vector(habitat_par_site_geo2_250m$Code_hab[habitat_par_site_geo2_250m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo2_250m$Code_hab[habitat_par_site_geo2_250m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo2_250m_station[m]=sum(Indice_rarete_geo2_250m[names(Indice_rarete_geo2_250m)%in%geo2_250m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo2_250m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo2_250m=sum(Indice_rarete_geo2_250m_station)                #somme des indices de rareté  de ttes les sp

     
#500m

  b=unique(habitat_par_site_geo2_500m$Code_hab)

    Indice_rarete_geo2_500m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo2_500m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo2_500m$Code_Site[habitat_par_site_geo2_500m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo2_500m$Code_Site[habitat_par_site_geo2_500m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo2_500m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo2_500m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo2_500m_station<-matrix(0,500,26)
Indice_rarete_geo2_500m_station<-NULL

    for (m in 1:26){
  geo2_500m_station[,m]<-c(as.vector(habitat_par_site_geo2_500m$Code_hab[habitat_par_site_geo2_500m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo2_500m$Code_hab[habitat_par_site_geo2_500m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo2_500m_station[m]=sum(Indice_rarete_geo2_500m[names(Indice_rarete_geo2_500m)%in%geo2_500m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo2_500m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo2_500m=sum(Indice_rarete_geo2_500m_station)                #somme des indices de rareté  de ttes les sp

     
#750m

  b=unique(habitat_par_site_geo2_750m$Code_hab)

    Indice_rarete_geo2_750m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo2_750m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo2_750m$Code_Site[habitat_par_site_geo2_750m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo2_750m$Code_Site[habitat_par_site_geo2_750m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo2_750m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo2_750m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo2_750m_station<-matrix(0,500,26)
Indice_rarete_geo2_750m_station<-NULL

    for (m in 1:26){
  geo2_750m_station[,m]<-c(as.vector(habitat_par_site_geo2_750m$Code_hab[habitat_par_site_geo2_750m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo2_750m$Code_hab[habitat_par_site_geo2_750m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo2_750m_station[m]=sum(Indice_rarete_geo2_750m[names(Indice_rarete_geo2_750m)%in%geo2_750m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo2_750m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo2_750m=sum(Indice_rarete_geo2_750m_station)                #somme des indices de rareté  de ttes les sp

      
#1000m

  b=unique(habitat_par_site_geo2_1000m$Code_hab)

    Indice_rarete_geo2_1000m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo2_1000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo2_1000m$Code_Site[habitat_par_site_geo2_1000m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo2_1000m$Code_Site[habitat_par_site_geo2_1000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo2_1000m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo2_1000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo2_1000m_station<-matrix(0,500,26)
Indice_rarete_geo2_1000m_station<-NULL

    for (m in 1:26){
  geo2_1000m_station[,m]<-c(as.vector(habitat_par_site_geo2_1000m$Code_hab[habitat_par_site_geo2_1000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo2_1000m$Code_hab[habitat_par_site_geo2_1000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo2_1000m_station[m]=sum(Indice_rarete_geo2_1000m[names(Indice_rarete_geo2_1000m)%in%geo2_1000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo2_1000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo2_1000m=sum(Indice_rarete_geo2_1000m_station)                #somme des indices de rareté  de ttes les sp

                  
#3000m

  b=unique(habitat_par_site_geo2_3000m$Code_hab)

    Indice_rarete_geo2_3000m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo2_3000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo2_3000m$Code_Site[habitat_par_site_geo2_3000m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo2_3000m$Code_Site[habitat_par_site_geo2_3000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo2_3000m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo2_3000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo2_3000m_station<-matrix(0,500,26)
Indice_rarete_geo2_3000m_station<-NULL

    for (m in 1:26){
  geo2_3000m_station[,m]<-c(as.vector(habitat_par_site_geo2_3000m$Code_hab[habitat_par_site_geo2_3000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo2_3000m$Code_hab[habitat_par_site_geo2_3000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo2_3000m_station[m]=sum(Indice_rarete_geo2_3000m[names(Indice_rarete_geo2_3000m)%in%geo2_3000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo2_3000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo2_3000m=sum(Indice_rarete_geo2_3000m_station)                #somme des indices de rareté  de ttes les sp

         
#5000m

  b=unique(habitat_par_site_geo2_5000m$Code_hab)

    Indice_rarete_geo2_5000m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo2_5000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo2_5000m$Code_Site[habitat_par_site_geo2_5000m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo2_5000m$Code_Site[habitat_par_site_geo2_5000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo2_5000m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo2_5000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo2_5000m_station<-matrix(0,500,26)
Indice_rarete_geo2_5000m_station<-NULL

    for (m in 1:26){
  geo2_5000m_station[,m]<-c(as.vector(habitat_par_site_geo2_5000m$Code_hab[habitat_par_site_geo2_5000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo2_5000m$Code_hab[habitat_par_site_geo2_5000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo2_5000m_station[m]=sum(Indice_rarete_geo2_5000m[names(Indice_rarete_geo2_5000m)%in%geo2_5000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo2_5000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo2_5000m=sum(Indice_rarete_geo2_5000m_station)                #somme des indices de rareté  de ttes les sp


source('script.rarete.complementarite.r')

#Selection des dites selon un scenario de rareté complementarité basé sur les habitats

#habitat géo1
selection.rarete.complementarite.geo2.30m<-run.bin(habitat_par_site_geo2_30m_pres_abs,Indice_rarete_geo2_30m_station)
selection.rarete.complementarite.geo2.60m<-run.bin(habitat_par_site_geo2_60m_pres_abs,Indice_rarete_geo2_60m_station)
selection.rarete.complementarite.geo2.100m<-run.bin(habitat_par_site_geo2_100m_pres_abs,Indice_rarete_geo2_100m_station)
selection.rarete.complementarite.geo2.250m<-run.bin(habitat_par_site_geo2_250m_pres_abs,Indice_rarete_geo2_250m_station)
selection.rarete.complementarite.geo2.500m<-run.bin(habitat_par_site_geo2_500m_pres_abs,Indice_rarete_geo2_500m_station)
selection.rarete.complementarite.geo2.750m<-run.bin(habitat_par_site_geo2_750m_pres_abs,Indice_rarete_geo2_750m_station)
selection.rarete.complementarite.geo2.1000m<-run.bin(habitat_par_site_geo2_1000m_pres_abs,Indice_rarete_geo2_1000m_station)
selection.rarete.complementarite.geo2.3000m<-run.bin(habitat_par_site_geo2_3000m_pres_abs,Indice_rarete_geo2_3000m_station)
selection.rarete.complementarite.geo2.5000m<-run.bin(habitat_par_site_geo2_5000m_pres_abs,Indice_rarete_geo2_5000m_station)

#DIVERSITE SELECTIONNE PAR LE SCRIPT RARETE COMPLEMENTARITE#
 source('script.diversite.rarete.selectionnee.r')
 
 diversite.selection.rarete.complementarite.geo2.30m<-species.rarete.evol(selection.rarete.complementarite.geo2.30m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo2.60m<-species.rarete.evol(selection.rarete.complementarite.geo2.60m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo2.100m<-species.rarete.evol(selection.rarete.complementarite.geo2.100m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo2.250m<-species.rarete.evol(selection.rarete.complementarite.geo2.250m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo2.500m<-species.rarete.evol(selection.rarete.complementarite.geo2.500m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo2.750m<-species.rarete.evol(selection.rarete.complementarite.geo2.750m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo2.1000m<-species.rarete.evol(selection.rarete.complementarite.geo2.1000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo2.3000m<-species.rarete.evol(selection.rarete.complementarite.geo2.3000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo2.5000m<-species.rarete.evol(selection.rarete.complementarite.geo2.5000m,Indice_rarete_station)

#CALCUL DES SAI


source('script.sai.r')
SAI_rarete.complementarite.geo2.30m<-SAI(diversite.selection.rarete.complementarite.geo2.30m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo2.60m<-SAI(diversite.selection.rarete.complementarite.geo2.60m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo2.100m<-SAI(diversite.selection.rarete.complementarite.geo2.100m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo2.250m<-SAI(diversite.selection.rarete.complementarite.geo2.250m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo2.500m<-SAI(diversite.selection.rarete.complementarite.geo2.500m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo2.750m<-SAI(diversite.selection.rarete.complementarite.geo2.750m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo2.1000m<-SAI(diversite.selection.rarete.complementarite.geo2.1000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo2.3000m<-SAI(diversite.selection.rarete.complementarite.geo2.3000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo2.5000m<-SAI(diversite.selection.rarete.complementarite.geo2.5000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)


 
############################HABITAT geo3###########################################


 ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances


relation.habitats.sites.geo3<-read.csv2("fc3_intersect_geo3.csv")


relation.habitats.sites.geo3<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$Code_Site%in%sites.etudies,]

#création d'une nouvelle table site-habitat par distance
#habitat géo1
relation.habitats.sites.geo3.30m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==30,]
relation.habitats.sites.geo3.60m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==60,]
relation.habitats.sites.geo3.100m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==100,]
relation.habitats.sites.geo3.250m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==250,]
relation.habitats.sites.geo3.500m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==500,]
relation.habitats.sites.geo3.750m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==750,]
relation.habitats.sites.geo3.1000m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==1000,]
relation.habitats.sites.geo3.3000m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==3000,]
relation.habitats.sites.geo3.5000m<-relation.habitats.sites.geo3[relation.habitats.sites.geo3$distance==5000,]



#Load de la requette "REQ_site_geo3_30m.csv"

 m=matrix(0,length(relation.habitats.sites.geo3.30m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo3.30m$Code_Site)
  m[,2]=relation.habitats.sites.geo3.30m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo3_30m=as.data.frame(unique(m))


 m=matrix(0,length(relation.habitats.sites.geo3.60m$Code_Site),2)
 m[,1]=as.vector(relation.habitats.sites.geo3.60m$Code_Site)
  m[,2]=relation.habitats.sites.geo3.60m$Code_hab
 colnames(m)=c('Code_Site','Code_hab')
 habitat_par_site_geo3_60m=as.data.frame(unique(m))


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


habitat_par_site_geo3_30m<-habitat_par_site_geo3_30m[habitat_par_site_geo3_30m$Code_Site%in%sites.etudies,]
habitat_par_site_geo3_60m<-habitat_par_site_geo3_60m[habitat_par_site_geo3_60m$Code_Site%in%sites.etudies,]
habitat_par_site_geo3_100m<-habitat_par_site_geo3_100m[habitat_par_site_geo3_100m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo3_250m<-habitat_par_site_geo3_250m[habitat_par_site_geo3_250m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo3_500m<-habitat_par_site_geo3_500m[habitat_par_site_geo3_500m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo3_750m<-habitat_par_site_geo3_750m[habitat_par_site_geo3_750m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo3_1000m<-habitat_par_site_geo3_1000m[habitat_par_site_geo3_1000m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo3_3000m<-habitat_par_site_geo3_3000m[habitat_par_site_geo3_3000m$Code_Site%in%sites.etudies,]
 habitat_par_site_geo3_5000m<-habitat_par_site_geo3_5000m[habitat_par_site_geo3_5000m$Code_Site%in%sites.etudies,]


#création de table de présence absence
#habitat géo1

habitat_par_site_geo3_30m_pres_abs<-table(habitat_par_site_geo3_30m$Code_Site,habitat_par_site_geo3_30m$Code_hab)
habitat_par_site_geo3_60m_pres_abs<-table(habitat_par_site_geo3_60m$Code_Site,habitat_par_site_geo3_60m$Code_hab)
habitat_par_site_geo3_100m_pres_abs<-table(habitat_par_site_geo3_100m$Code_Site,habitat_par_site_geo3_100m$Code_hab)
habitat_par_site_geo3_250m_pres_abs<-table(habitat_par_site_geo3_250m$Code_Site,habitat_par_site_geo3_250m$Code_hab)
habitat_par_site_geo3_500m_pres_abs<-table(habitat_par_site_geo3_500m$Code_Site,habitat_par_site_geo3_500m$Code_hab)
habitat_par_site_geo3_750m_pres_abs<-table(habitat_par_site_geo3_750m$Code_Site,habitat_par_site_geo3_750m$Code_hab)
habitat_par_site_geo3_1000m_pres_abs<-table(habitat_par_site_geo3_1000m$Code_Site,habitat_par_site_geo3_1000m$Code_hab)
habitat_par_site_geo3_3000m_pres_abs<-table(habitat_par_site_geo3_3000m$Code_Site,habitat_par_site_geo3_3000m$Code_hab)
habitat_par_site_geo3_5000m_pres_abs<-table(habitat_par_site_geo3_5000m$Code_Site,habitat_par_site_geo3_5000m$Code_hab)

#calcul indice de rareté des habitats

#30m

  b=unique(habitat_par_site_geo3_30m$Code_hab)

    Indice_rarete_geo3_30m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo3_30m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo3_30m$Code_Site[habitat_par_site_geo3_30m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo3_30m$Code_Site[habitat_par_site_geo3_30m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo3_30m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo3_30m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo3_30m_station<-matrix(0,500,26)
Indice_rarete_geo3_30m_station<-NULL

    for (m in 1:26){
  geo3_30m_station[,m]<-c(as.vector(habitat_par_site_geo3_30m$Code_hab[habitat_par_site_geo3_30m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo3_30m$Code_hab[habitat_par_site_geo3_30m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo3_30m_station[m]=sum(Indice_rarete_geo3_30m[names(Indice_rarete_geo3_30m)%in%geo3_30m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo3_30m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo3_30m=sum(Indice_rarete_geo3_30m_station)                #somme des indices de rareté  de ttes les sp


   
#60m

  b=unique(habitat_par_site_geo3_60m$Code_hab)

    Indice_rarete_geo3_60m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo3_60m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo3_60m$Code_Site[habitat_par_site_geo3_60m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo3_60m$Code_Site[habitat_par_site_geo3_60m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo3_60m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo3_60m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo3_60m_station<-matrix(0,500,26)
Indice_rarete_geo3_60m_station<-NULL

    for (m in 1:26){
  geo3_60m_station[,m]<-c(as.vector(habitat_par_site_geo3_60m$Code_hab[habitat_par_site_geo3_60m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo3_60m$Code_hab[habitat_par_site_geo3_60m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo3_60m_station[m]=sum(Indice_rarete_geo3_60m[names(Indice_rarete_geo3_60m)%in%geo3_60m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo3_60m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo3_60m=sum(Indice_rarete_geo3_60m_station)                #somme des indices de rareté  de ttes les sp


    
#100m

  b=unique(habitat_par_site_geo3_100m$Code_hab)

    Indice_rarete_geo3_100m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo3_100m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo3_100m$Code_Site[habitat_par_site_geo3_100m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo3_100m$Code_Site[habitat_par_site_geo3_100m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo3_100m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo3_100m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo3_100m_station<-matrix(0,500,26)
Indice_rarete_geo3_100m_station<-NULL

    for (m in 1:26){
  geo3_100m_station[,m]<-c(as.vector(habitat_par_site_geo3_100m$Code_hab[habitat_par_site_geo3_100m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo3_100m$Code_hab[habitat_par_site_geo3_100m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo3_100m_station[m]=sum(Indice_rarete_geo3_100m[names(Indice_rarete_geo3_100m)%in%geo3_100m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo3_100m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo3_100m=sum(Indice_rarete_geo3_100m_station)                #somme des indices de rareté  de ttes les sp

              
#250m

  b=unique(habitat_par_site_geo3_250m$Code_hab)

    Indice_rarete_geo3_250m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo3_250m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo3_250m$Code_Site[habitat_par_site_geo3_250m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo3_250m$Code_Site[habitat_par_site_geo3_250m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo3_250m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo3_250m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo3_250m_station<-matrix(0,500,26)
Indice_rarete_geo3_250m_station<-NULL

    for (m in 1:26){
  geo3_250m_station[,m]<-c(as.vector(habitat_par_site_geo3_250m$Code_hab[habitat_par_site_geo3_250m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo3_250m$Code_hab[habitat_par_site_geo3_250m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo3_250m_station[m]=sum(Indice_rarete_geo3_250m[names(Indice_rarete_geo3_250m)%in%geo3_250m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo3_250m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo3_250m=sum(Indice_rarete_geo3_250m_station)                #somme des indices de rareté  de ttes les sp

     
#500m

  b=unique(habitat_par_site_geo3_500m$Code_hab)

    Indice_rarete_geo3_500m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo3_500m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo3_500m$Code_Site[habitat_par_site_geo3_500m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo3_500m$Code_Site[habitat_par_site_geo3_500m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo3_500m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo3_500m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo3_500m_station<-matrix(0,500,26)
Indice_rarete_geo3_500m_station<-NULL

    for (m in 1:26){
  geo3_500m_station[,m]<-c(as.vector(habitat_par_site_geo3_500m$Code_hab[habitat_par_site_geo3_500m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo3_500m$Code_hab[habitat_par_site_geo3_500m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo3_500m_station[m]=sum(Indice_rarete_geo3_500m[names(Indice_rarete_geo3_500m)%in%geo3_500m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo3_500m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo3_500m=sum(Indice_rarete_geo3_500m_station)                #somme des indices de rareté  de ttes les sp

     
#750m

  b=unique(habitat_par_site_geo3_750m$Code_hab)

    Indice_rarete_geo3_750m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo3_750m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo3_750m$Code_Site[habitat_par_site_geo3_750m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo3_750m$Code_Site[habitat_par_site_geo3_750m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo3_750m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo3_750m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo3_750m_station<-matrix(0,500,26)
Indice_rarete_geo3_750m_station<-NULL

    for (m in 1:26){
  geo3_750m_station[,m]<-c(as.vector(habitat_par_site_geo3_750m$Code_hab[habitat_par_site_geo3_750m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo3_750m$Code_hab[habitat_par_site_geo3_750m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo3_750m_station[m]=sum(Indice_rarete_geo3_750m[names(Indice_rarete_geo3_750m)%in%geo3_750m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo3_750m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo3_750m=sum(Indice_rarete_geo3_750m_station)                #somme des indices de rareté  de ttes les sp

      
#1000m

  b=unique(habitat_par_site_geo3_1000m$Code_hab)

    Indice_rarete_geo3_1000m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo3_1000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo3_1000m$Code_Site[habitat_par_site_geo3_1000m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo3_1000m$Code_Site[habitat_par_site_geo3_1000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo3_1000m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo3_1000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo3_1000m_station<-matrix(0,500,26)
Indice_rarete_geo3_1000m_station<-NULL

    for (m in 1:26){
  geo3_1000m_station[,m]<-c(as.vector(habitat_par_site_geo3_1000m$Code_hab[habitat_par_site_geo3_1000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo3_1000m$Code_hab[habitat_par_site_geo3_1000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo3_1000m_station[m]=sum(Indice_rarete_geo3_1000m[names(Indice_rarete_geo3_1000m)%in%geo3_1000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo3_1000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo3_1000m=sum(Indice_rarete_geo3_1000m_station)                #somme des indices de rareté  de ttes les sp

                  
#3000m

  b=unique(habitat_par_site_geo3_3000m$Code_hab)

    Indice_rarete_geo3_3000m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo3_3000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo3_3000m$Code_Site[habitat_par_site_geo3_3000m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo3_3000m$Code_Site[habitat_par_site_geo3_3000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo3_3000m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo3_3000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo3_3000m_station<-matrix(0,500,26)
Indice_rarete_geo3_3000m_station<-NULL

    for (m in 1:26){
  geo3_3000m_station[,m]<-c(as.vector(habitat_par_site_geo3_3000m$Code_hab[habitat_par_site_geo3_3000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo3_3000m$Code_hab[habitat_par_site_geo3_3000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo3_3000m_station[m]=sum(Indice_rarete_geo3_3000m[names(Indice_rarete_geo3_3000m)%in%geo3_3000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo3_3000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo3_3000m=sum(Indice_rarete_geo3_3000m_station)                #somme des indices de rareté  de ttes les sp

         
#5000m

  b=unique(habitat_par_site_geo3_5000m$Code_hab)

    Indice_rarete_geo3_5000m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo3_5000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo3_5000m$Code_Site[habitat_par_site_geo3_5000m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo3_5000m$Code_Site[habitat_par_site_geo3_5000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo3_5000m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo3_5000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo3_5000m_station<-matrix(0,500,26)
Indice_rarete_geo3_5000m_station<-NULL

    for (m in 1:26){
  geo3_5000m_station[,m]<-c(as.vector(habitat_par_site_geo3_5000m$Code_hab[habitat_par_site_geo3_5000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo3_5000m$Code_hab[habitat_par_site_geo3_5000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo3_5000m_station[m]=sum(Indice_rarete_geo3_5000m[names(Indice_rarete_geo3_5000m)%in%geo3_5000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo3_5000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo3_5000m=sum(Indice_rarete_geo3_5000m_station)                #somme des indices de rareté  de ttes les sp


source('script.rarete.complementarite.r')

#Selection des dites selon un scenario de rareté complementarité basé sur les habitats

#habitat géo1
selection.rarete.complementarite.geo3.30m<-run.bin(habitat_par_site_geo3_30m_pres_abs,Indice_rarete_geo3_30m_station)
selection.rarete.complementarite.geo3.60m<-run.bin(habitat_par_site_geo3_60m_pres_abs,Indice_rarete_geo3_60m_station)
selection.rarete.complementarite.geo3.100m<-run.bin(habitat_par_site_geo3_100m_pres_abs,Indice_rarete_geo3_100m_station)
selection.rarete.complementarite.geo3.250m<-run.bin(habitat_par_site_geo3_250m_pres_abs,Indice_rarete_geo3_250m_station)
selection.rarete.complementarite.geo3.500m<-run.bin(habitat_par_site_geo3_500m_pres_abs,Indice_rarete_geo3_500m_station)
selection.rarete.complementarite.geo3.750m<-run.bin(habitat_par_site_geo3_750m_pres_abs,Indice_rarete_geo3_750m_station)
selection.rarete.complementarite.geo3.1000m<-run.bin(habitat_par_site_geo3_1000m_pres_abs,Indice_rarete_geo3_1000m_station)
selection.rarete.complementarite.geo3.3000m<-run.bin(habitat_par_site_geo3_3000m_pres_abs,Indice_rarete_geo3_3000m_station)
selection.rarete.complementarite.geo3.5000m<-run.bin(habitat_par_site_geo3_5000m_pres_abs,Indice_rarete_geo3_5000m_station)

#DIVERSITE SELECTIONNE PAR LE SCRIPT RARETE COMPLEMENTARITE#
 source('script.diversite.rarete.selectionnee.r')
 
 diversite.selection.rarete.complementarite.geo3.30m<-species.rarete.evol(selection.rarete.complementarite.geo3.30m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo3.60m<-species.rarete.evol(selection.rarete.complementarite.geo3.60m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo3.100m<-species.rarete.evol(selection.rarete.complementarite.geo3.100m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo3.250m<-species.rarete.evol(selection.rarete.complementarite.geo3.250m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo3.500m<-species.rarete.evol(selection.rarete.complementarite.geo3.500m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo3.750m<-species.rarete.evol(selection.rarete.complementarite.geo3.750m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo3.1000m<-species.rarete.evol(selection.rarete.complementarite.geo3.1000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo3.3000m<-species.rarete.evol(selection.rarete.complementarite.geo3.3000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo3.5000m<-species.rarete.evol(selection.rarete.complementarite.geo3.5000m,Indice_rarete_station)

#CALCUL DES SAI


source('script.sai.r')
SAI_rarete.complementarite.geo3.30m<-SAI(diversite.selection.rarete.complementarite.geo3.30m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo3.60m<-SAI(diversite.selection.rarete.complementarite.geo3.60m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo3.100m<-SAI(diversite.selection.rarete.complementarite.geo3.100m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo3.250m<-SAI(diversite.selection.rarete.complementarite.geo3.250m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo3.500m<-SAI(diversite.selection.rarete.complementarite.geo3.500m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo3.750m<-SAI(diversite.selection.rarete.complementarite.geo3.750m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo3.1000m<-SAI(diversite.selection.rarete.complementarite.geo3.1000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo3.3000m<-SAI(diversite.selection.rarete.complementarite.geo3.3000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo3.5000m<-SAI(diversite.selection.rarete.complementarite.geo3.5000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)



 
############################HABITAT bent1###########################################


 ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances


relation.habitats.sites.bent1<-read.csv2("fc3_intersect_bent1.csv")


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



#Load de la requette "REQ_site_bent1_30m.csv"

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

#calcul indice de rareté des habitats

#30m

  b=unique(habitat_par_site_bent1_30m$Code_hab)

    Indice_rarete_bent1_30m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_bent1_30m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent1_30m$Code_Site[habitat_par_site_bent1_30m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_bent1_30m$Code_Site[habitat_par_site_bent1_30m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent1_30m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent1_30m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent1_30m_station<-matrix(0,500,26)
Indice_rarete_bent1_30m_station<-NULL

    for (m in 1:26){
  bent1_30m_station[,m]<-c(as.vector(habitat_par_site_bent1_30m$Code_hab[habitat_par_site_bent1_30m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent1_30m$Code_hab[habitat_par_site_bent1_30m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent1_30m_station[m]=sum(Indice_rarete_bent1_30m[names(Indice_rarete_bent1_30m)%in%bent1_30m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent1_30m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent1_30m=sum(Indice_rarete_bent1_30m_station)                #somme des indices de rareté  de ttes les sp


   
#60m

  b=unique(habitat_par_site_bent1_60m$Code_hab)

    Indice_rarete_bent1_60m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_bent1_60m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent1_60m$Code_Site[habitat_par_site_bent1_60m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_bent1_60m$Code_Site[habitat_par_site_bent1_60m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent1_60m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent1_60m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent1_60m_station<-matrix(0,500,26)
Indice_rarete_bent1_60m_station<-NULL

    for (m in 1:26){
  bent1_60m_station[,m]<-c(as.vector(habitat_par_site_bent1_60m$Code_hab[habitat_par_site_bent1_60m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent1_60m$Code_hab[habitat_par_site_bent1_60m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent1_60m_station[m]=sum(Indice_rarete_bent1_60m[names(Indice_rarete_bent1_60m)%in%bent1_60m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent1_60m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent1_60m=sum(Indice_rarete_bent1_60m_station)                #somme des indices de rareté  de ttes les sp


    
#100m

  b=unique(habitat_par_site_bent1_100m$Code_hab)

    Indice_rarete_bent1_100m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_bent1_100m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent1_100m$Code_Site[habitat_par_site_bent1_100m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_bent1_100m$Code_Site[habitat_par_site_bent1_100m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent1_100m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent1_100m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent1_100m_station<-matrix(0,500,26)
Indice_rarete_bent1_100m_station<-NULL

    for (m in 1:26){
  bent1_100m_station[,m]<-c(as.vector(habitat_par_site_bent1_100m$Code_hab[habitat_par_site_bent1_100m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent1_100m$Code_hab[habitat_par_site_bent1_100m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent1_100m_station[m]=sum(Indice_rarete_bent1_100m[names(Indice_rarete_bent1_100m)%in%bent1_100m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent1_100m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent1_100m=sum(Indice_rarete_bent1_100m_station)                #somme des indices de rareté  de ttes les sp

              
#250m

  b=unique(habitat_par_site_bent1_250m$Code_hab)

    Indice_rarete_bent1_250m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_bent1_250m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent1_250m$Code_Site[habitat_par_site_bent1_250m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_bent1_250m$Code_Site[habitat_par_site_bent1_250m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent1_250m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent1_250m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent1_250m_station<-matrix(0,500,26)
Indice_rarete_bent1_250m_station<-NULL

    for (m in 1:26){
  bent1_250m_station[,m]<-c(as.vector(habitat_par_site_bent1_250m$Code_hab[habitat_par_site_bent1_250m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent1_250m$Code_hab[habitat_par_site_bent1_250m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent1_250m_station[m]=sum(Indice_rarete_bent1_250m[names(Indice_rarete_bent1_250m)%in%bent1_250m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent1_250m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent1_250m=sum(Indice_rarete_bent1_250m_station)                #somme des indices de rareté  de ttes les sp

     
#500m

  b=unique(habitat_par_site_bent1_500m$Code_hab)

    Indice_rarete_bent1_500m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_bent1_500m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent1_500m$Code_Site[habitat_par_site_bent1_500m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_bent1_500m$Code_Site[habitat_par_site_bent1_500m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent1_500m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent1_500m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent1_500m_station<-matrix(0,500,26)
Indice_rarete_bent1_500m_station<-NULL

    for (m in 1:26){
  bent1_500m_station[,m]<-c(as.vector(habitat_par_site_bent1_500m$Code_hab[habitat_par_site_bent1_500m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent1_500m$Code_hab[habitat_par_site_bent1_500m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent1_500m_station[m]=sum(Indice_rarete_bent1_500m[names(Indice_rarete_bent1_500m)%in%bent1_500m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent1_500m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent1_500m=sum(Indice_rarete_bent1_500m_station)                #somme des indices de rareté  de ttes les sp

     
#750m

  b=unique(habitat_par_site_bent1_750m$Code_hab)

    Indice_rarete_bent1_750m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_bent1_750m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent1_750m$Code_Site[habitat_par_site_bent1_750m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_bent1_750m$Code_Site[habitat_par_site_bent1_750m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent1_750m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent1_750m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent1_750m_station<-matrix(0,500,26)
Indice_rarete_bent1_750m_station<-NULL

    for (m in 1:26){
  bent1_750m_station[,m]<-c(as.vector(habitat_par_site_bent1_750m$Code_hab[habitat_par_site_bent1_750m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent1_750m$Code_hab[habitat_par_site_bent1_750m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent1_750m_station[m]=sum(Indice_rarete_bent1_750m[names(Indice_rarete_bent1_750m)%in%bent1_750m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent1_750m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent1_750m=sum(Indice_rarete_bent1_750m_station)                #somme des indices de rareté  de ttes les sp

      
#1000m

  b=unique(habitat_par_site_bent1_1000m$Code_hab)

    Indice_rarete_bent1_1000m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_bent1_1000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent1_1000m$Code_Site[habitat_par_site_bent1_1000m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_bent1_1000m$Code_Site[habitat_par_site_bent1_1000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent1_1000m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent1_1000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent1_1000m_station<-matrix(0,500,26)
Indice_rarete_bent1_1000m_station<-NULL

    for (m in 1:26){
  bent1_1000m_station[,m]<-c(as.vector(habitat_par_site_bent1_1000m$Code_hab[habitat_par_site_bent1_1000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent1_1000m$Code_hab[habitat_par_site_bent1_1000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent1_1000m_station[m]=sum(Indice_rarete_bent1_1000m[names(Indice_rarete_bent1_1000m)%in%bent1_1000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent1_1000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent1_1000m=sum(Indice_rarete_bent1_1000m_station)                #somme des indices de rareté  de ttes les sp

                  
#3000m

  b=unique(habitat_par_site_bent1_3000m$Code_hab)

    Indice_rarete_bent1_3000m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_bent1_3000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent1_3000m$Code_Site[habitat_par_site_bent1_3000m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_bent1_3000m$Code_Site[habitat_par_site_bent1_3000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent1_3000m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent1_3000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent1_3000m_station<-matrix(0,500,26)
Indice_rarete_bent1_3000m_station<-NULL

    for (m in 1:26){
  bent1_3000m_station[,m]<-c(as.vector(habitat_par_site_bent1_3000m$Code_hab[habitat_par_site_bent1_3000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent1_3000m$Code_hab[habitat_par_site_bent1_3000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent1_3000m_station[m]=sum(Indice_rarete_bent1_3000m[names(Indice_rarete_bent1_3000m)%in%bent1_3000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent1_3000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent1_3000m=sum(Indice_rarete_bent1_3000m_station)                #somme des indices de rareté  de ttes les sp

         
#5000m

  b=unique(habitat_par_site_bent1_5000m$Code_hab)

    Indice_rarete_bent1_5000m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_bent1_5000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent1_5000m$Code_Site[habitat_par_site_bent1_5000m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_bent1_5000m$Code_Site[habitat_par_site_bent1_5000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent1_5000m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent1_5000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent1_5000m_station<-matrix(0,500,26)
Indice_rarete_bent1_5000m_station<-NULL

    for (m in 1:26){
  bent1_5000m_station[,m]<-c(as.vector(habitat_par_site_bent1_5000m$Code_hab[habitat_par_site_bent1_5000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent1_5000m$Code_hab[habitat_par_site_bent1_5000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent1_5000m_station[m]=sum(Indice_rarete_bent1_5000m[names(Indice_rarete_bent1_5000m)%in%bent1_5000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent1_5000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent1_5000m=sum(Indice_rarete_bent1_5000m_station)                #somme des indices de rareté  de ttes les sp


source('script.rarete.complementarite.r')

#Selection des dites selon un scenario de rareté complementarité basé sur les habitats

#habitat géo1
selection.rarete.complementarite.bent1.30m<-run.bin(habitat_par_site_bent1_30m_pres_abs,Indice_rarete_bent1_30m_station)
selection.rarete.complementarite.bent1.60m<-run.bin(habitat_par_site_bent1_60m_pres_abs,Indice_rarete_bent1_60m_station)
selection.rarete.complementarite.bent1.100m<-run.bin(habitat_par_site_bent1_100m_pres_abs,Indice_rarete_bent1_100m_station)
selection.rarete.complementarite.bent1.250m<-run.bin(habitat_par_site_bent1_250m_pres_abs,Indice_rarete_bent1_250m_station)
selection.rarete.complementarite.bent1.500m<-run.bin(habitat_par_site_bent1_500m_pres_abs,Indice_rarete_bent1_500m_station)
selection.rarete.complementarite.bent1.750m<-run.bin(habitat_par_site_bent1_750m_pres_abs,Indice_rarete_bent1_750m_station)
selection.rarete.complementarite.bent1.1000m<-run.bin(habitat_par_site_bent1_1000m_pres_abs,Indice_rarete_bent1_1000m_station)
selection.rarete.complementarite.bent1.3000m<-run.bin(habitat_par_site_bent1_3000m_pres_abs,Indice_rarete_bent1_3000m_station)
selection.rarete.complementarite.bent1.5000m<-run.bin(habitat_par_site_bent1_5000m_pres_abs,Indice_rarete_bent1_5000m_station)

#DIVERSITE SELECTIONNE PAR LE SCRIPT RARETE COMPLEMENTARITE#
 source('script.diversite.rarete.selectionnee.r')
 
 diversite.selection.rarete.complementarite.bent1.30m<-species.rarete.evol(selection.rarete.complementarite.bent1.30m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent1.60m<-species.rarete.evol(selection.rarete.complementarite.bent1.60m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent1.100m<-species.rarete.evol(selection.rarete.complementarite.bent1.100m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent1.250m<-species.rarete.evol(selection.rarete.complementarite.bent1.250m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent1.500m<-species.rarete.evol(selection.rarete.complementarite.bent1.500m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent1.750m<-species.rarete.evol(selection.rarete.complementarite.bent1.750m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent1.1000m<-species.rarete.evol(selection.rarete.complementarite.bent1.1000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent1.3000m<-species.rarete.evol(selection.rarete.complementarite.bent1.3000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent1.5000m<-species.rarete.evol(selection.rarete.complementarite.bent1.5000m,Indice_rarete_station)

#CALCUL DES SAI


source('script.sai.r')
SAI_rarete.complementarite.bent1.30m<-SAI(diversite.selection.rarete.complementarite.bent1.30m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent1.60m<-SAI(diversite.selection.rarete.complementarite.bent1.60m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent1.100m<-SAI(diversite.selection.rarete.complementarite.bent1.100m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent1.250m<-SAI(diversite.selection.rarete.complementarite.bent1.250m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent1.500m<-SAI(diversite.selection.rarete.complementarite.bent1.500m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent1.750m<-SAI(diversite.selection.rarete.complementarite.bent1.750m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent1.1000m<-SAI(diversite.selection.rarete.complementarite.bent1.1000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent1.3000m<-SAI(diversite.selection.rarete.complementarite.bent1.3000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent1.5000m<-SAI(diversite.selection.rarete.complementarite.bent1.5000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)



 
############################HABITAT bent2###########################################


 ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances


relation.habitats.sites.bent2<-read.csv2("fc3_intersect_bent2.csv")


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



#Load de la requette "REQ_site_bent2_30m.csv"

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

#calcul indice de rareté des habitats

#30m

  b=unique(habitat_par_site_bent2_30m$Code_hab)

    Indice_rarete_bent2_30m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_bent2_30m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent2_30m$Code_Site[habitat_par_site_bent2_30m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_bent2_30m$Code_Site[habitat_par_site_bent2_30m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent2_30m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent2_30m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent2_30m_station<-matrix(0,500,26)
Indice_rarete_bent2_30m_station<-NULL

    for (m in 1:26){
  bent2_30m_station[,m]<-c(as.vector(habitat_par_site_bent2_30m$Code_hab[habitat_par_site_bent2_30m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent2_30m$Code_hab[habitat_par_site_bent2_30m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent2_30m_station[m]=sum(Indice_rarete_bent2_30m[names(Indice_rarete_bent2_30m)%in%bent2_30m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent2_30m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent2_30m=sum(Indice_rarete_bent2_30m_station)                #somme des indices de rareté  de ttes les sp


   
#60m

  b=unique(habitat_par_site_bent2_60m$Code_hab)

    Indice_rarete_bent2_60m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_bent2_60m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent2_60m$Code_Site[habitat_par_site_bent2_60m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_bent2_60m$Code_Site[habitat_par_site_bent2_60m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent2_60m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent2_60m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent2_60m_station<-matrix(0,500,26)
Indice_rarete_bent2_60m_station<-NULL

    for (m in 1:26){
  bent2_60m_station[,m]<-c(as.vector(habitat_par_site_bent2_60m$Code_hab[habitat_par_site_bent2_60m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent2_60m$Code_hab[habitat_par_site_bent2_60m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent2_60m_station[m]=sum(Indice_rarete_bent2_60m[names(Indice_rarete_bent2_60m)%in%bent2_60m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent2_60m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent2_60m=sum(Indice_rarete_bent2_60m_station)                #somme des indices de rareté  de ttes les sp


    
#100m

  b=unique(habitat_par_site_bent2_100m$Code_hab)

    Indice_rarete_bent2_100m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_bent2_100m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent2_100m$Code_Site[habitat_par_site_bent2_100m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_bent2_100m$Code_Site[habitat_par_site_bent2_100m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent2_100m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent2_100m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent2_100m_station<-matrix(0,500,26)
Indice_rarete_bent2_100m_station<-NULL

    for (m in 1:26){
  bent2_100m_station[,m]<-c(as.vector(habitat_par_site_bent2_100m$Code_hab[habitat_par_site_bent2_100m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent2_100m$Code_hab[habitat_par_site_bent2_100m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent2_100m_station[m]=sum(Indice_rarete_bent2_100m[names(Indice_rarete_bent2_100m)%in%bent2_100m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent2_100m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent2_100m=sum(Indice_rarete_bent2_100m_station)                #somme des indices de rareté  de ttes les sp

              
#250m

  b=unique(habitat_par_site_bent2_250m$Code_hab)

    Indice_rarete_bent2_250m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_bent2_250m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent2_250m$Code_Site[habitat_par_site_bent2_250m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_bent2_250m$Code_Site[habitat_par_site_bent2_250m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent2_250m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent2_250m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent2_250m_station<-matrix(0,500,26)
Indice_rarete_bent2_250m_station<-NULL

    for (m in 1:26){
  bent2_250m_station[,m]<-c(as.vector(habitat_par_site_bent2_250m$Code_hab[habitat_par_site_bent2_250m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent2_250m$Code_hab[habitat_par_site_bent2_250m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent2_250m_station[m]=sum(Indice_rarete_bent2_250m[names(Indice_rarete_bent2_250m)%in%bent2_250m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent2_250m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent2_250m=sum(Indice_rarete_bent2_250m_station)                #somme des indices de rareté  de ttes les sp

     
#500m

  b=unique(habitat_par_site_bent2_500m$Code_hab)

    Indice_rarete_bent2_500m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_bent2_500m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent2_500m$Code_Site[habitat_par_site_bent2_500m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_bent2_500m$Code_Site[habitat_par_site_bent2_500m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent2_500m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent2_500m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent2_500m_station<-matrix(0,500,26)
Indice_rarete_bent2_500m_station<-NULL

    for (m in 1:26){
  bent2_500m_station[,m]<-c(as.vector(habitat_par_site_bent2_500m$Code_hab[habitat_par_site_bent2_500m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent2_500m$Code_hab[habitat_par_site_bent2_500m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent2_500m_station[m]=sum(Indice_rarete_bent2_500m[names(Indice_rarete_bent2_500m)%in%bent2_500m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent2_500m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent2_500m=sum(Indice_rarete_bent2_500m_station)                #somme des indices de rareté  de ttes les sp

     
#750m

  b=unique(habitat_par_site_bent2_750m$Code_hab)

    Indice_rarete_bent2_750m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_bent2_750m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent2_750m$Code_Site[habitat_par_site_bent2_750m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_bent2_750m$Code_Site[habitat_par_site_bent2_750m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent2_750m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent2_750m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent2_750m_station<-matrix(0,500,26)
Indice_rarete_bent2_750m_station<-NULL

    for (m in 1:26){
  bent2_750m_station[,m]<-c(as.vector(habitat_par_site_bent2_750m$Code_hab[habitat_par_site_bent2_750m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent2_750m$Code_hab[habitat_par_site_bent2_750m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent2_750m_station[m]=sum(Indice_rarete_bent2_750m[names(Indice_rarete_bent2_750m)%in%bent2_750m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent2_750m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent2_750m=sum(Indice_rarete_bent2_750m_station)                #somme des indices de rareté  de ttes les sp

      
#1000m

  b=unique(habitat_par_site_bent2_1000m$Code_hab)

    Indice_rarete_bent2_1000m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_bent2_1000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent2_1000m$Code_Site[habitat_par_site_bent2_1000m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_bent2_1000m$Code_Site[habitat_par_site_bent2_1000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent2_1000m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent2_1000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent2_1000m_station<-matrix(0,500,26)
Indice_rarete_bent2_1000m_station<-NULL

    for (m in 1:26){
  bent2_1000m_station[,m]<-c(as.vector(habitat_par_site_bent2_1000m$Code_hab[habitat_par_site_bent2_1000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent2_1000m$Code_hab[habitat_par_site_bent2_1000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent2_1000m_station[m]=sum(Indice_rarete_bent2_1000m[names(Indice_rarete_bent2_1000m)%in%bent2_1000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent2_1000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent2_1000m=sum(Indice_rarete_bent2_1000m_station)                #somme des indices de rareté  de ttes les sp

                  
#3000m

  b=unique(habitat_par_site_bent2_3000m$Code_hab)

    Indice_rarete_bent2_3000m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_bent2_3000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent2_3000m$Code_Site[habitat_par_site_bent2_3000m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_bent2_3000m$Code_Site[habitat_par_site_bent2_3000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent2_3000m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent2_3000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent2_3000m_station<-matrix(0,500,26)
Indice_rarete_bent2_3000m_station<-NULL

    for (m in 1:26){
  bent2_3000m_station[,m]<-c(as.vector(habitat_par_site_bent2_3000m$Code_hab[habitat_par_site_bent2_3000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent2_3000m$Code_hab[habitat_par_site_bent2_3000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent2_3000m_station[m]=sum(Indice_rarete_bent2_3000m[names(Indice_rarete_bent2_3000m)%in%bent2_3000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent2_3000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent2_3000m=sum(Indice_rarete_bent2_3000m_station)                #somme des indices de rareté  de ttes les sp

         
#5000m

  b=unique(habitat_par_site_bent2_5000m$Code_hab)

    Indice_rarete_bent2_5000m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_bent2_5000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent2_5000m$Code_Site[habitat_par_site_bent2_5000m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_bent2_5000m$Code_Site[habitat_par_site_bent2_5000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent2_5000m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent2_5000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent2_5000m_station<-matrix(0,500,26)
Indice_rarete_bent2_5000m_station<-NULL

    for (m in 1:26){
  bent2_5000m_station[,m]<-c(as.vector(habitat_par_site_bent2_5000m$Code_hab[habitat_par_site_bent2_5000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent2_5000m$Code_hab[habitat_par_site_bent2_5000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent2_5000m_station[m]=sum(Indice_rarete_bent2_5000m[names(Indice_rarete_bent2_5000m)%in%bent2_5000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent2_5000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent2_5000m=sum(Indice_rarete_bent2_5000m_station)                #somme des indices de rareté  de ttes les sp


source('script.rarete.complementarite.r')

#Selection des dites selon un scenario de rareté complementarité basé sur les habitats

#habitat géo1
selection.rarete.complementarite.bent2.30m<-run.bin(habitat_par_site_bent2_30m_pres_abs,Indice_rarete_bent2_30m_station)
selection.rarete.complementarite.bent2.60m<-run.bin(habitat_par_site_bent2_60m_pres_abs,Indice_rarete_bent2_60m_station)
selection.rarete.complementarite.bent2.100m<-run.bin(habitat_par_site_bent2_100m_pres_abs,Indice_rarete_bent2_100m_station)
selection.rarete.complementarite.bent2.250m<-run.bin(habitat_par_site_bent2_250m_pres_abs,Indice_rarete_bent2_250m_station)
selection.rarete.complementarite.bent2.500m<-run.bin(habitat_par_site_bent2_500m_pres_abs,Indice_rarete_bent2_500m_station)
selection.rarete.complementarite.bent2.750m<-run.bin(habitat_par_site_bent2_750m_pres_abs,Indice_rarete_bent2_750m_station)
selection.rarete.complementarite.bent2.1000m<-run.bin(habitat_par_site_bent2_1000m_pres_abs,Indice_rarete_bent2_1000m_station)
selection.rarete.complementarite.bent2.3000m<-run.bin(habitat_par_site_bent2_3000m_pres_abs,Indice_rarete_bent2_3000m_station)
selection.rarete.complementarite.bent2.5000m<-run.bin(habitat_par_site_bent2_5000m_pres_abs,Indice_rarete_bent2_5000m_station)

#DIVERSITE SELECTIONNE PAR LE SCRIPT RARETE COMPLEMENTARITE#
 source('script.diversite.rarete.selectionnee.r')
 
 diversite.selection.rarete.complementarite.bent2.30m<-species.rarete.evol(selection.rarete.complementarite.bent2.30m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent2.60m<-species.rarete.evol(selection.rarete.complementarite.bent2.60m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent2.100m<-species.rarete.evol(selection.rarete.complementarite.bent2.100m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent2.250m<-species.rarete.evol(selection.rarete.complementarite.bent2.250m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent2.500m<-species.rarete.evol(selection.rarete.complementarite.bent2.500m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent2.750m<-species.rarete.evol(selection.rarete.complementarite.bent2.750m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent2.1000m<-species.rarete.evol(selection.rarete.complementarite.bent2.1000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent2.3000m<-species.rarete.evol(selection.rarete.complementarite.bent2.3000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent2.5000m<-species.rarete.evol(selection.rarete.complementarite.bent2.5000m,Indice_rarete_station)

#CALCUL DES SAI


source('script.sai.r')
SAI_rarete.complementarite.bent2.30m<-SAI(diversite.selection.rarete.complementarite.bent2.30m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent2.60m<-SAI(diversite.selection.rarete.complementarite.bent2.60m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent2.100m<-SAI(diversite.selection.rarete.complementarite.bent2.100m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent2.250m<-SAI(diversite.selection.rarete.complementarite.bent2.250m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent2.500m<-SAI(diversite.selection.rarete.complementarite.bent2.500m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent2.750m<-SAI(diversite.selection.rarete.complementarite.bent2.750m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent2.1000m<-SAI(diversite.selection.rarete.complementarite.bent2.1000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent2.3000m<-SAI(diversite.selection.rarete.complementarite.bent2.3000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent2.5000m<-SAI(diversite.selection.rarete.complementarite.bent2.5000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)



 
############################HABITAT bent3###########################################


 ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances


relation.habitats.sites.bent3<-read.csv2("fc3_intersect_bent3.csv")


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



#Load de la requette "REQ_site_bent3_30m.csv"

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

#calcul indice de rareté des habitats

#30m

  b=unique(habitat_par_site_bent3_30m$Code_hab)

    Indice_rarete_bent3_30m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_bent3_30m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent3_30m$Code_Site[habitat_par_site_bent3_30m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_bent3_30m$Code_Site[habitat_par_site_bent3_30m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent3_30m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent3_30m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent3_30m_station<-matrix(0,500,26)
Indice_rarete_bent3_30m_station<-NULL

    for (m in 1:26){
  bent3_30m_station[,m]<-c(as.vector(habitat_par_site_bent3_30m$Code_hab[habitat_par_site_bent3_30m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent3_30m$Code_hab[habitat_par_site_bent3_30m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent3_30m_station[m]=sum(Indice_rarete_bent3_30m[names(Indice_rarete_bent3_30m)%in%bent3_30m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent3_30m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent3_30m=sum(Indice_rarete_bent3_30m_station)                #somme des indices de rareté  de ttes les sp


   
#60m

  b=unique(habitat_par_site_bent3_60m$Code_hab)

    Indice_rarete_bent3_60m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_bent3_60m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent3_60m$Code_Site[habitat_par_site_bent3_60m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_bent3_60m$Code_Site[habitat_par_site_bent3_60m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent3_60m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent3_60m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent3_60m_station<-matrix(0,500,26)
Indice_rarete_bent3_60m_station<-NULL

    for (m in 1:26){
  bent3_60m_station[,m]<-c(as.vector(habitat_par_site_bent3_60m$Code_hab[habitat_par_site_bent3_60m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent3_60m$Code_hab[habitat_par_site_bent3_60m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent3_60m_station[m]=sum(Indice_rarete_bent3_60m[names(Indice_rarete_bent3_60m)%in%bent3_60m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent3_60m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent3_60m=sum(Indice_rarete_bent3_60m_station)                #somme des indices de rareté  de ttes les sp


    
#100m

  b=unique(habitat_par_site_bent3_100m$Code_hab)

    Indice_rarete_bent3_100m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_bent3_100m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent3_100m$Code_Site[habitat_par_site_bent3_100m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_bent3_100m$Code_Site[habitat_par_site_bent3_100m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent3_100m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent3_100m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent3_100m_station<-matrix(0,500,26)
Indice_rarete_bent3_100m_station<-NULL

    for (m in 1:26){
  bent3_100m_station[,m]<-c(as.vector(habitat_par_site_bent3_100m$Code_hab[habitat_par_site_bent3_100m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent3_100m$Code_hab[habitat_par_site_bent3_100m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent3_100m_station[m]=sum(Indice_rarete_bent3_100m[names(Indice_rarete_bent3_100m)%in%bent3_100m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent3_100m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent3_100m=sum(Indice_rarete_bent3_100m_station)                #somme des indices de rareté  de ttes les sp

              
#250m

  b=unique(habitat_par_site_bent3_250m$Code_hab)

    Indice_rarete_bent3_250m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_bent3_250m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent3_250m$Code_Site[habitat_par_site_bent3_250m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_bent3_250m$Code_Site[habitat_par_site_bent3_250m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent3_250m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent3_250m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent3_250m_station<-matrix(0,500,26)
Indice_rarete_bent3_250m_station<-NULL

    for (m in 1:26){
  bent3_250m_station[,m]<-c(as.vector(habitat_par_site_bent3_250m$Code_hab[habitat_par_site_bent3_250m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent3_250m$Code_hab[habitat_par_site_bent3_250m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent3_250m_station[m]=sum(Indice_rarete_bent3_250m[names(Indice_rarete_bent3_250m)%in%bent3_250m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent3_250m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent3_250m=sum(Indice_rarete_bent3_250m_station)                #somme des indices de rareté  de ttes les sp

     
#500m

  b=unique(habitat_par_site_bent3_500m$Code_hab)

    Indice_rarete_bent3_500m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_bent3_500m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent3_500m$Code_Site[habitat_par_site_bent3_500m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_bent3_500m$Code_Site[habitat_par_site_bent3_500m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent3_500m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent3_500m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent3_500m_station<-matrix(0,500,26)
Indice_rarete_bent3_500m_station<-NULL

    for (m in 1:26){
  bent3_500m_station[,m]<-c(as.vector(habitat_par_site_bent3_500m$Code_hab[habitat_par_site_bent3_500m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent3_500m$Code_hab[habitat_par_site_bent3_500m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent3_500m_station[m]=sum(Indice_rarete_bent3_500m[names(Indice_rarete_bent3_500m)%in%bent3_500m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent3_500m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent3_500m=sum(Indice_rarete_bent3_500m_station)                #somme des indices de rareté  de ttes les sp

     
#750m

  b=unique(habitat_par_site_bent3_750m$Code_hab)

    Indice_rarete_bent3_750m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_bent3_750m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent3_750m$Code_Site[habitat_par_site_bent3_750m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_bent3_750m$Code_Site[habitat_par_site_bent3_750m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent3_750m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent3_750m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent3_750m_station<-matrix(0,500,26)
Indice_rarete_bent3_750m_station<-NULL

    for (m in 1:26){
  bent3_750m_station[,m]<-c(as.vector(habitat_par_site_bent3_750m$Code_hab[habitat_par_site_bent3_750m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent3_750m$Code_hab[habitat_par_site_bent3_750m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent3_750m_station[m]=sum(Indice_rarete_bent3_750m[names(Indice_rarete_bent3_750m)%in%bent3_750m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent3_750m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent3_750m=sum(Indice_rarete_bent3_750m_station)                #somme des indices de rareté  de ttes les sp

      
#1000m

  b=unique(habitat_par_site_bent3_1000m$Code_hab)

    Indice_rarete_bent3_1000m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_bent3_1000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent3_1000m$Code_Site[habitat_par_site_bent3_1000m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_bent3_1000m$Code_Site[habitat_par_site_bent3_1000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent3_1000m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent3_1000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent3_1000m_station<-matrix(0,500,26)
Indice_rarete_bent3_1000m_station<-NULL

    for (m in 1:26){
  bent3_1000m_station[,m]<-c(as.vector(habitat_par_site_bent3_1000m$Code_hab[habitat_par_site_bent3_1000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent3_1000m$Code_hab[habitat_par_site_bent3_1000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent3_1000m_station[m]=sum(Indice_rarete_bent3_1000m[names(Indice_rarete_bent3_1000m)%in%bent3_1000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent3_1000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent3_1000m=sum(Indice_rarete_bent3_1000m_station)                #somme des indices de rareté  de ttes les sp

                  
#3000m

  b=unique(habitat_par_site_bent3_3000m$Code_hab)

    Indice_rarete_bent3_3000m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_bent3_3000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent3_3000m$Code_Site[habitat_par_site_bent3_3000m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_bent3_3000m$Code_Site[habitat_par_site_bent3_3000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent3_3000m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent3_3000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent3_3000m_station<-matrix(0,500,26)
Indice_rarete_bent3_3000m_station<-NULL

    for (m in 1:26){
  bent3_3000m_station[,m]<-c(as.vector(habitat_par_site_bent3_3000m$Code_hab[habitat_par_site_bent3_3000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent3_3000m$Code_hab[habitat_par_site_bent3_3000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent3_3000m_station[m]=sum(Indice_rarete_bent3_3000m[names(Indice_rarete_bent3_3000m)%in%bent3_3000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent3_3000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent3_3000m=sum(Indice_rarete_bent3_3000m_station)                #somme des indices de rareté  de ttes les sp

         
#5000m

  b=unique(habitat_par_site_bent3_5000m$Code_hab)

    Indice_rarete_bent3_5000m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_bent3_5000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_bent3_5000m$Code_Site[habitat_par_site_bent3_5000m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_bent3_5000m$Code_Site[habitat_par_site_bent3_5000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_bent3_5000m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_bent3_5000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

bent3_5000m_station<-matrix(0,500,26)
Indice_rarete_bent3_5000m_station<-NULL

    for (m in 1:26){
  bent3_5000m_station[,m]<-c(as.vector(habitat_par_site_bent3_5000m$Code_hab[habitat_par_site_bent3_5000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_bent3_5000m$Code_hab[habitat_par_site_bent3_5000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_bent3_5000m_station[m]=sum(Indice_rarete_bent3_5000m[names(Indice_rarete_bent3_5000m)%in%bent3_5000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_bent3_5000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_bent3_5000m=sum(Indice_rarete_bent3_5000m_station)                #somme des indices de rareté  de ttes les sp


source('script.rarete.complementarite.r')

#Selection des dites selon un scenario de rareté complementarité basé sur les habitats

#habitat géo1
selection.rarete.complementarite.bent3.30m<-run.bin(habitat_par_site_bent3_30m_pres_abs,Indice_rarete_bent3_30m_station)
selection.rarete.complementarite.bent3.60m<-run.bin(habitat_par_site_bent3_60m_pres_abs,Indice_rarete_bent3_60m_station)
selection.rarete.complementarite.bent3.100m<-run.bin(habitat_par_site_bent3_100m_pres_abs,Indice_rarete_bent3_100m_station)
selection.rarete.complementarite.bent3.250m<-run.bin(habitat_par_site_bent3_250m_pres_abs,Indice_rarete_bent3_250m_station)
selection.rarete.complementarite.bent3.500m<-run.bin(habitat_par_site_bent3_500m_pres_abs,Indice_rarete_bent3_500m_station)
selection.rarete.complementarite.bent3.750m<-run.bin(habitat_par_site_bent3_750m_pres_abs,Indice_rarete_bent3_750m_station)
selection.rarete.complementarite.bent3.1000m<-run.bin(habitat_par_site_bent3_1000m_pres_abs,Indice_rarete_bent3_1000m_station)
selection.rarete.complementarite.bent3.3000m<-run.bin(habitat_par_site_bent3_3000m_pres_abs,Indice_rarete_bent3_3000m_station)
selection.rarete.complementarite.bent3.5000m<-run.bin(habitat_par_site_bent3_5000m_pres_abs,Indice_rarete_bent3_5000m_station)

#DIVERSITE SELECTIONNE PAR LE SCRIPT RARETE COMPLEMENTARITE#
 source('script.diversite.rarete.selectionnee.r')
 
 diversite.selection.rarete.complementarite.bent3.30m<-species.rarete.evol(selection.rarete.complementarite.bent3.30m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent3.60m<-species.rarete.evol(selection.rarete.complementarite.bent3.60m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent3.100m<-species.rarete.evol(selection.rarete.complementarite.bent3.100m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent3.250m<-species.rarete.evol(selection.rarete.complementarite.bent3.250m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent3.500m<-species.rarete.evol(selection.rarete.complementarite.bent3.500m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent3.750m<-species.rarete.evol(selection.rarete.complementarite.bent3.750m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent3.1000m<-species.rarete.evol(selection.rarete.complementarite.bent3.1000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent3.3000m<-species.rarete.evol(selection.rarete.complementarite.bent3.3000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.bent3.5000m<-species.rarete.evol(selection.rarete.complementarite.bent3.5000m,Indice_rarete_station)

#CALCUL DES SAI


source('script.sai.r')
SAI_rarete.complementarite.bent3.30m<-SAI(diversite.selection.rarete.complementarite.bent3.30m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent3.60m<-SAI(diversite.selection.rarete.complementarite.bent3.60m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent3.100m<-SAI(diversite.selection.rarete.complementarite.bent3.100m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent3.250m<-SAI(diversite.selection.rarete.complementarite.bent3.250m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent3.500m<-SAI(diversite.selection.rarete.complementarite.bent3.500m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent3.750m<-SAI(diversite.selection.rarete.complementarite.bent3.750m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent3.1000m<-SAI(diversite.selection.rarete.complementarite.bent3.1000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent3.3000m<-SAI(diversite.selection.rarete.complementarite.bent3.3000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.bent3.5000m<-SAI(diversite.selection.rarete.complementarite.bent3.5000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)



 
############################HABITAT geo1_geo2###########################################


 ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances


relation.habitats.sites.geo1_geo2<-read.csv2("fc3_intersect_geo1_geo2.csv")


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



#Load de la requette "REQ_site_geo1_geo2_30m.csv"

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

#calcul indice de rareté des habitats

#30m

  b=unique(habitat_par_site_geo1_geo2_30m$Code_hab)

    Indice_rarete_geo1_geo2_30m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo1_geo2_30m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_30m$Code_Site[habitat_par_site_geo1_geo2_30m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo1_geo2_30m$Code_Site[habitat_par_site_geo1_geo2_30m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_30m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_30m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_30m_station<-matrix(0,500,26)
Indice_rarete_geo1_geo2_30m_station<-NULL

    for (m in 1:26){
  geo1_geo2_30m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_30m$Code_hab[habitat_par_site_geo1_geo2_30m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_30m$Code_hab[habitat_par_site_geo1_geo2_30m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_30m_station[m]=sum(Indice_rarete_geo1_geo2_30m[names(Indice_rarete_geo1_geo2_30m)%in%geo1_geo2_30m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_30m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_30m=sum(Indice_rarete_geo1_geo2_30m_station)                #somme des indices de rareté  de ttes les sp


   
#60m

  b=unique(habitat_par_site_geo1_geo2_60m$Code_hab)

    Indice_rarete_geo1_geo2_60m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo1_geo2_60m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_60m$Code_Site[habitat_par_site_geo1_geo2_60m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo1_geo2_60m$Code_Site[habitat_par_site_geo1_geo2_60m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_60m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_60m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_60m_station<-matrix(0,500,26)
Indice_rarete_geo1_geo2_60m_station<-NULL

    for (m in 1:26){
  geo1_geo2_60m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_60m$Code_hab[habitat_par_site_geo1_geo2_60m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_60m$Code_hab[habitat_par_site_geo1_geo2_60m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_60m_station[m]=sum(Indice_rarete_geo1_geo2_60m[names(Indice_rarete_geo1_geo2_60m)%in%geo1_geo2_60m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_60m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_60m=sum(Indice_rarete_geo1_geo2_60m_station)                #somme des indices de rareté  de ttes les sp


    
#100m

  b=unique(habitat_par_site_geo1_geo2_100m$Code_hab)

    Indice_rarete_geo1_geo2_100m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo1_geo2_100m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_100m$Code_Site[habitat_par_site_geo1_geo2_100m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo1_geo2_100m$Code_Site[habitat_par_site_geo1_geo2_100m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_100m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_100m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_100m_station<-matrix(0,500,26)
Indice_rarete_geo1_geo2_100m_station<-NULL

    for (m in 1:26){
  geo1_geo2_100m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_100m$Code_hab[habitat_par_site_geo1_geo2_100m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_100m$Code_hab[habitat_par_site_geo1_geo2_100m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_100m_station[m]=sum(Indice_rarete_geo1_geo2_100m[names(Indice_rarete_geo1_geo2_100m)%in%geo1_geo2_100m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_100m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_100m=sum(Indice_rarete_geo1_geo2_100m_station)                #somme des indices de rareté  de ttes les sp

              
#250m

  b=unique(habitat_par_site_geo1_geo2_250m$Code_hab)

    Indice_rarete_geo1_geo2_250m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo1_geo2_250m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_250m$Code_Site[habitat_par_site_geo1_geo2_250m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo1_geo2_250m$Code_Site[habitat_par_site_geo1_geo2_250m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_250m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_250m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_250m_station<-matrix(0,500,26)
Indice_rarete_geo1_geo2_250m_station<-NULL

    for (m in 1:26){
  geo1_geo2_250m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_250m$Code_hab[habitat_par_site_geo1_geo2_250m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_250m$Code_hab[habitat_par_site_geo1_geo2_250m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_250m_station[m]=sum(Indice_rarete_geo1_geo2_250m[names(Indice_rarete_geo1_geo2_250m)%in%geo1_geo2_250m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_250m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_250m=sum(Indice_rarete_geo1_geo2_250m_station)                #somme des indices de rareté  de ttes les sp

     
#500m

  b=unique(habitat_par_site_geo1_geo2_500m$Code_hab)

    Indice_rarete_geo1_geo2_500m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo1_geo2_500m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_500m$Code_Site[habitat_par_site_geo1_geo2_500m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo1_geo2_500m$Code_Site[habitat_par_site_geo1_geo2_500m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_500m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_500m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_500m_station<-matrix(0,500,26)
Indice_rarete_geo1_geo2_500m_station<-NULL

    for (m in 1:26){
  geo1_geo2_500m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_500m$Code_hab[habitat_par_site_geo1_geo2_500m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_500m$Code_hab[habitat_par_site_geo1_geo2_500m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_500m_station[m]=sum(Indice_rarete_geo1_geo2_500m[names(Indice_rarete_geo1_geo2_500m)%in%geo1_geo2_500m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_500m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_500m=sum(Indice_rarete_geo1_geo2_500m_station)                #somme des indices de rareté  de ttes les sp

     
#750m

  b=unique(habitat_par_site_geo1_geo2_750m$Code_hab)

    Indice_rarete_geo1_geo2_750m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo1_geo2_750m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_750m$Code_Site[habitat_par_site_geo1_geo2_750m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo1_geo2_750m$Code_Site[habitat_par_site_geo1_geo2_750m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_750m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_750m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_750m_station<-matrix(0,500,26)
Indice_rarete_geo1_geo2_750m_station<-NULL

    for (m in 1:26){
  geo1_geo2_750m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_750m$Code_hab[habitat_par_site_geo1_geo2_750m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_750m$Code_hab[habitat_par_site_geo1_geo2_750m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_750m_station[m]=sum(Indice_rarete_geo1_geo2_750m[names(Indice_rarete_geo1_geo2_750m)%in%geo1_geo2_750m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_750m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_750m=sum(Indice_rarete_geo1_geo2_750m_station)                #somme des indices de rareté  de ttes les sp

      
#1000m

  b=unique(habitat_par_site_geo1_geo2_1000m$Code_hab)

    Indice_rarete_geo1_geo2_1000m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo1_geo2_1000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_1000m$Code_Site[habitat_par_site_geo1_geo2_1000m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo1_geo2_1000m$Code_Site[habitat_par_site_geo1_geo2_1000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_1000m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_1000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_1000m_station<-matrix(0,500,26)
Indice_rarete_geo1_geo2_1000m_station<-NULL

    for (m in 1:26){
  geo1_geo2_1000m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_1000m$Code_hab[habitat_par_site_geo1_geo2_1000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_1000m$Code_hab[habitat_par_site_geo1_geo2_1000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_1000m_station[m]=sum(Indice_rarete_geo1_geo2_1000m[names(Indice_rarete_geo1_geo2_1000m)%in%geo1_geo2_1000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_1000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_1000m=sum(Indice_rarete_geo1_geo2_1000m_station)                #somme des indices de rareté  de ttes les sp

                  
#3000m

  b=unique(habitat_par_site_geo1_geo2_3000m$Code_hab)

    Indice_rarete_geo1_geo2_3000m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo1_geo2_3000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_3000m$Code_Site[habitat_par_site_geo1_geo2_3000m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo1_geo2_3000m$Code_Site[habitat_par_site_geo1_geo2_3000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_3000m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_3000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_3000m_station<-matrix(0,500,26)
Indice_rarete_geo1_geo2_3000m_station<-NULL

    for (m in 1:26){
  geo1_geo2_3000m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_3000m$Code_hab[habitat_par_site_geo1_geo2_3000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_3000m$Code_hab[habitat_par_site_geo1_geo2_3000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_3000m_station[m]=sum(Indice_rarete_geo1_geo2_3000m[names(Indice_rarete_geo1_geo2_3000m)%in%geo1_geo2_3000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_3000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_3000m=sum(Indice_rarete_geo1_geo2_3000m_station)                #somme des indices de rareté  de ttes les sp

         
#5000m

  b=unique(habitat_par_site_geo1_geo2_5000m$Code_hab)

    Indice_rarete_geo1_geo2_5000m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo1_geo2_5000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_5000m$Code_Site[habitat_par_site_geo1_geo2_5000m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo1_geo2_5000m$Code_Site[habitat_par_site_geo1_geo2_5000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_5000m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_5000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_5000m_station<-matrix(0,500,26)
Indice_rarete_geo1_geo2_5000m_station<-NULL

    for (m in 1:26){
  geo1_geo2_5000m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_5000m$Code_hab[habitat_par_site_geo1_geo2_5000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_5000m$Code_hab[habitat_par_site_geo1_geo2_5000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_5000m_station[m]=sum(Indice_rarete_geo1_geo2_5000m[names(Indice_rarete_geo1_geo2_5000m)%in%geo1_geo2_5000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_5000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_5000m=sum(Indice_rarete_geo1_geo2_5000m_station)                #somme des indices de rareté  de ttes les sp


source('script.rarete.complementarite.r')

#Selection des dites selon un scenario de rareté complementarité basé sur les habitats

#habitat géo1
selection.rarete.complementarite.geo1_geo2.30m<-run.bin(habitat_par_site_geo1_geo2_30m_pres_abs,Indice_rarete_geo1_geo2_30m_station)
selection.rarete.complementarite.geo1_geo2.60m<-run.bin(habitat_par_site_geo1_geo2_60m_pres_abs,Indice_rarete_geo1_geo2_60m_station)
selection.rarete.complementarite.geo1_geo2.100m<-run.bin(habitat_par_site_geo1_geo2_100m_pres_abs,Indice_rarete_geo1_geo2_100m_station)
selection.rarete.complementarite.geo1_geo2.250m<-run.bin(habitat_par_site_geo1_geo2_250m_pres_abs,Indice_rarete_geo1_geo2_250m_station)
selection.rarete.complementarite.geo1_geo2.500m<-run.bin(habitat_par_site_geo1_geo2_500m_pres_abs,Indice_rarete_geo1_geo2_500m_station)
selection.rarete.complementarite.geo1_geo2.750m<-run.bin(habitat_par_site_geo1_geo2_750m_pres_abs,Indice_rarete_geo1_geo2_750m_station)
selection.rarete.complementarite.geo1_geo2.1000m<-run.bin(habitat_par_site_geo1_geo2_1000m_pres_abs,Indice_rarete_geo1_geo2_1000m_station)
selection.rarete.complementarite.geo1_geo2.3000m<-run.bin(habitat_par_site_geo1_geo2_3000m_pres_abs,Indice_rarete_geo1_geo2_3000m_station)
selection.rarete.complementarite.geo1_geo2.5000m<-run.bin(habitat_par_site_geo1_geo2_5000m_pres_abs,Indice_rarete_geo1_geo2_5000m_station)

#DIVERSITE SELECTIONNE PAR LE SCRIPT RARETE COMPLEMENTARITE#
 source('script.diversite.rarete.selectionnee.r')
 
 diversite.selection.rarete.complementarite.geo1_geo2.30m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2.30m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1_geo2.60m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2.60m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1_geo2.100m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2.100m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1_geo2.250m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2.250m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1_geo2.500m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2.500m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1_geo2.750m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2.750m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1_geo2.1000m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2.1000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1_geo2.3000m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2.3000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1_geo2.5000m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2.5000m,Indice_rarete_station)

#CALCUL DES SAI


source('script.sai.r')
SAI_rarete.complementarite.geo1_geo2.30m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2.30m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1_geo2.60m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2.60m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1_geo2.100m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2.100m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1_geo2.250m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2.250m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1_geo2.500m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2.500m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1_geo2.750m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2.750m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1_geo2.1000m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2.1000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1_geo2.3000m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2.3000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1_geo2.5000m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2.5000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)



 
############################HABITAT geo1_geo2_geo3###########################################


 ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances


relation.habitats.sites.geo1_geo2_geo3<-read.csv2("fc3_intersect_geo1_geo2_geo3.csv")


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



#Load de la requette "REQ_site_geo1_geo2_geo3_30m.csv"

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

#calcul indice de rareté des habitats

#30m

  b=unique(habitat_par_site_geo1_geo2_geo3_30m$Code_hab)

    Indice_rarete_geo1_geo2_geo3_30m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo1_geo2_geo3_30m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_30m$Code_Site[habitat_par_site_geo1_geo2_geo3_30m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo1_geo2_geo3_30m$Code_Site[habitat_par_site_geo1_geo2_geo3_30m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_geo3_30m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_geo3_30m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_geo3_30m_station<-matrix(0,500,26)
Indice_rarete_geo1_geo2_geo3_30m_station<-NULL

    for (m in 1:26){
  geo1_geo2_geo3_30m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_30m$Code_hab[habitat_par_site_geo1_geo2_geo3_30m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_geo3_30m$Code_hab[habitat_par_site_geo1_geo2_geo3_30m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_geo3_30m_station[m]=sum(Indice_rarete_geo1_geo2_geo3_30m[names(Indice_rarete_geo1_geo2_geo3_30m)%in%geo1_geo2_geo3_30m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_geo3_30m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_geo3_30m=sum(Indice_rarete_geo1_geo2_geo3_30m_station)                #somme des indices de rareté  de ttes les sp


   
#60m

  b=unique(habitat_par_site_geo1_geo2_geo3_60m$Code_hab)

    Indice_rarete_geo1_geo2_geo3_60m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo1_geo2_geo3_60m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_60m$Code_Site[habitat_par_site_geo1_geo2_geo3_60m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo1_geo2_geo3_60m$Code_Site[habitat_par_site_geo1_geo2_geo3_60m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_geo3_60m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_geo3_60m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_geo3_60m_station<-matrix(0,500,26)
Indice_rarete_geo1_geo2_geo3_60m_station<-NULL

    for (m in 1:26){
  geo1_geo2_geo3_60m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_60m$Code_hab[habitat_par_site_geo1_geo2_geo3_60m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_geo3_60m$Code_hab[habitat_par_site_geo1_geo2_geo3_60m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_geo3_60m_station[m]=sum(Indice_rarete_geo1_geo2_geo3_60m[names(Indice_rarete_geo1_geo2_geo3_60m)%in%geo1_geo2_geo3_60m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_geo3_60m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_geo3_60m=sum(Indice_rarete_geo1_geo2_geo3_60m_station)                #somme des indices de rareté  de ttes les sp


    
#100m

  b=unique(habitat_par_site_geo1_geo2_geo3_100m$Code_hab)

    Indice_rarete_geo1_geo2_geo3_100m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo1_geo2_geo3_100m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_100m$Code_Site[habitat_par_site_geo1_geo2_geo3_100m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo1_geo2_geo3_100m$Code_Site[habitat_par_site_geo1_geo2_geo3_100m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_geo3_100m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_geo3_100m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_geo3_100m_station<-matrix(0,500,26)
Indice_rarete_geo1_geo2_geo3_100m_station<-NULL

    for (m in 1:26){
  geo1_geo2_geo3_100m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_100m$Code_hab[habitat_par_site_geo1_geo2_geo3_100m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_geo3_100m$Code_hab[habitat_par_site_geo1_geo2_geo3_100m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_geo3_100m_station[m]=sum(Indice_rarete_geo1_geo2_geo3_100m[names(Indice_rarete_geo1_geo2_geo3_100m)%in%geo1_geo2_geo3_100m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_geo3_100m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_geo3_100m=sum(Indice_rarete_geo1_geo2_geo3_100m_station)                #somme des indices de rareté  de ttes les sp

              
#250m

  b=unique(habitat_par_site_geo1_geo2_geo3_250m$Code_hab)

    Indice_rarete_geo1_geo2_geo3_250m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo1_geo2_geo3_250m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_250m$Code_Site[habitat_par_site_geo1_geo2_geo3_250m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo1_geo2_geo3_250m$Code_Site[habitat_par_site_geo1_geo2_geo3_250m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_geo3_250m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_geo3_250m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_geo3_250m_station<-matrix(0,500,26)
Indice_rarete_geo1_geo2_geo3_250m_station<-NULL

    for (m in 1:26){
  geo1_geo2_geo3_250m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_250m$Code_hab[habitat_par_site_geo1_geo2_geo3_250m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_geo3_250m$Code_hab[habitat_par_site_geo1_geo2_geo3_250m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_geo3_250m_station[m]=sum(Indice_rarete_geo1_geo2_geo3_250m[names(Indice_rarete_geo1_geo2_geo3_250m)%in%geo1_geo2_geo3_250m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_geo3_250m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_geo3_250m=sum(Indice_rarete_geo1_geo2_geo3_250m_station)                #somme des indices de rareté  de ttes les sp

     
#500m

  b=unique(habitat_par_site_geo1_geo2_geo3_500m$Code_hab)

    Indice_rarete_geo1_geo2_geo3_500m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo1_geo2_geo3_500m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_500m$Code_Site[habitat_par_site_geo1_geo2_geo3_500m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo1_geo2_geo3_500m$Code_Site[habitat_par_site_geo1_geo2_geo3_500m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_geo3_500m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_geo3_500m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_geo3_500m_station<-matrix(0,500,26)
Indice_rarete_geo1_geo2_geo3_500m_station<-NULL

    for (m in 1:26){
  geo1_geo2_geo3_500m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_500m$Code_hab[habitat_par_site_geo1_geo2_geo3_500m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_geo3_500m$Code_hab[habitat_par_site_geo1_geo2_geo3_500m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_geo3_500m_station[m]=sum(Indice_rarete_geo1_geo2_geo3_500m[names(Indice_rarete_geo1_geo2_geo3_500m)%in%geo1_geo2_geo3_500m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_geo3_500m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_geo3_500m=sum(Indice_rarete_geo1_geo2_geo3_500m_station)                #somme des indices de rareté  de ttes les sp

     
#750m

  b=unique(habitat_par_site_geo1_geo2_geo3_750m$Code_hab)

    Indice_rarete_geo1_geo2_geo3_750m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo1_geo2_geo3_750m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_750m$Code_Site[habitat_par_site_geo1_geo2_geo3_750m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo1_geo2_geo3_750m$Code_Site[habitat_par_site_geo1_geo2_geo3_750m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_geo3_750m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_geo3_750m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_geo3_750m_station<-matrix(0,500,26)
Indice_rarete_geo1_geo2_geo3_750m_station<-NULL

    for (m in 1:26){
  geo1_geo2_geo3_750m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_750m$Code_hab[habitat_par_site_geo1_geo2_geo3_750m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_geo3_750m$Code_hab[habitat_par_site_geo1_geo2_geo3_750m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_geo3_750m_station[m]=sum(Indice_rarete_geo1_geo2_geo3_750m[names(Indice_rarete_geo1_geo2_geo3_750m)%in%geo1_geo2_geo3_750m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_geo3_750m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_geo3_750m=sum(Indice_rarete_geo1_geo2_geo3_750m_station)                #somme des indices de rareté  de ttes les sp

      
#1000m

  b=unique(habitat_par_site_geo1_geo2_geo3_1000m$Code_hab)

    Indice_rarete_geo1_geo2_geo3_1000m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo1_geo2_geo3_1000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_1000m$Code_Site[habitat_par_site_geo1_geo2_geo3_1000m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo1_geo2_geo3_1000m$Code_Site[habitat_par_site_geo1_geo2_geo3_1000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_geo3_1000m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_geo3_1000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_geo3_1000m_station<-matrix(0,500,26)
Indice_rarete_geo1_geo2_geo3_1000m_station<-NULL

    for (m in 1:26){
  geo1_geo2_geo3_1000m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_1000m$Code_hab[habitat_par_site_geo1_geo2_geo3_1000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_geo3_1000m$Code_hab[habitat_par_site_geo1_geo2_geo3_1000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_geo3_1000m_station[m]=sum(Indice_rarete_geo1_geo2_geo3_1000m[names(Indice_rarete_geo1_geo2_geo3_1000m)%in%geo1_geo2_geo3_1000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_geo3_1000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_geo3_1000m=sum(Indice_rarete_geo1_geo2_geo3_1000m_station)                #somme des indices de rareté  de ttes les sp

                  
#3000m

  b=unique(habitat_par_site_geo1_geo2_geo3_3000m$Code_hab)

    Indice_rarete_geo1_geo2_geo3_3000m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo1_geo2_geo3_3000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_3000m$Code_Site[habitat_par_site_geo1_geo2_geo3_3000m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo1_geo2_geo3_3000m$Code_Site[habitat_par_site_geo1_geo2_geo3_3000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_geo3_3000m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_geo3_3000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_geo3_3000m_station<-matrix(0,500,26)
Indice_rarete_geo1_geo2_geo3_3000m_station<-NULL

    for (m in 1:26){
  geo1_geo2_geo3_3000m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_3000m$Code_hab[habitat_par_site_geo1_geo2_geo3_3000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_geo3_3000m$Code_hab[habitat_par_site_geo1_geo2_geo3_3000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_geo3_3000m_station[m]=sum(Indice_rarete_geo1_geo2_geo3_3000m[names(Indice_rarete_geo1_geo2_geo3_3000m)%in%geo1_geo2_geo3_3000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_geo3_3000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_geo3_3000m=sum(Indice_rarete_geo1_geo2_geo3_3000m_station)                #somme des indices de rareté  de ttes les sp

         
#5000m

  b=unique(habitat_par_site_geo1_geo2_geo3_5000m$Code_hab)

    Indice_rarete_geo1_geo2_geo3_5000m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_geo1_geo2_geo3_5000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_5000m$Code_Site[habitat_par_site_geo1_geo2_geo3_5000m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_geo1_geo2_geo3_5000m$Code_Site[habitat_par_site_geo1_geo2_geo3_5000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_geo1_geo2_geo3_5000m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_geo1_geo2_geo3_5000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

geo1_geo2_geo3_5000m_station<-matrix(0,500,26)
Indice_rarete_geo1_geo2_geo3_5000m_station<-NULL
    for (m in 1:26){
  geo1_geo2_geo3_5000m_station[,m]<-c(as.vector(habitat_par_site_geo1_geo2_geo3_5000m$Code_hab[habitat_par_site_geo1_geo2_geo3_5000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_geo1_geo2_geo3_5000m$Code_hab[habitat_par_site_geo1_geo2_geo3_5000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_geo1_geo2_geo3_5000m_station[m]=sum(Indice_rarete_geo1_geo2_geo3_5000m[names(Indice_rarete_geo1_geo2_geo3_5000m)%in%geo1_geo2_geo3_5000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_geo1_geo2_geo3_5000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_geo1_geo2_geo3_5000m=sum(Indice_rarete_geo1_geo2_geo3_5000m_station)                #somme des indices de rareté  de ttes les sp


source('script.rarete.complementarite.r')

#Selection des dites selon un scenario de rareté complementarité basé sur les habitats

#habitat géo1
selection.rarete.complementarite.geo1_geo2_geo3.30m<-run.bin(habitat_par_site_geo1_geo2_geo3_30m_pres_abs,Indice_rarete_geo1_geo2_geo3_30m_station)
selection.rarete.complementarite.geo1_geo2_geo3.60m<-run.bin(habitat_par_site_geo1_geo2_geo3_60m_pres_abs,Indice_rarete_geo1_geo2_geo3_60m_station)
selection.rarete.complementarite.geo1_geo2_geo3.100m<-run.bin(habitat_par_site_geo1_geo2_geo3_100m_pres_abs,Indice_rarete_geo1_geo2_geo3_100m_station)
selection.rarete.complementarite.geo1_geo2_geo3.250m<-run.bin(habitat_par_site_geo1_geo2_geo3_250m_pres_abs,Indice_rarete_geo1_geo2_geo3_250m_station)
selection.rarete.complementarite.geo1_geo2_geo3.500m<-run.bin(habitat_par_site_geo1_geo2_geo3_500m_pres_abs,Indice_rarete_geo1_geo2_geo3_500m_station)
selection.rarete.complementarite.geo1_geo2_geo3.750m<-run.bin(habitat_par_site_geo1_geo2_geo3_750m_pres_abs,Indice_rarete_geo1_geo2_geo3_750m_station)
selection.rarete.complementarite.geo1_geo2_geo3.1000m<-run.bin(habitat_par_site_geo1_geo2_geo3_1000m_pres_abs,Indice_rarete_geo1_geo2_geo3_1000m_station)
selection.rarete.complementarite.geo1_geo2_geo3.3000m<-run.bin(habitat_par_site_geo1_geo2_geo3_3000m_pres_abs,Indice_rarete_geo1_geo2_geo3_3000m_station)
selection.rarete.complementarite.geo1_geo2_geo3.5000m<-run.bin(habitat_par_site_geo1_geo2_geo3_5000m_pres_abs,Indice_rarete_geo1_geo2_geo3_5000m_station)

#DIVERSITE SELECTIONNE PAR LE SCRIPT RARETE COMPLEMENTARITE#
 source('script.diversite.rarete.selectionnee.r')
 
 diversite.selection.rarete.complementarite.geo1_geo2_geo3.30m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2_geo3.30m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1_geo2_geo3.60m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2_geo3.60m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1_geo2_geo3.100m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2_geo3.100m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1_geo2_geo3.250m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2_geo3.250m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1_geo2_geo3.500m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2_geo3.500m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1_geo2_geo3.750m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2_geo3.750m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1_geo2_geo3.1000m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2_geo3.1000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1_geo2_geo3.3000m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2_geo3.3000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.geo1_geo2_geo3.5000m<-species.rarete.evol(selection.rarete.complementarite.geo1_geo2_geo3.5000m,Indice_rarete_station)

#CALCUL DES SAI


source('script.sai.r')
SAI_rarete.complementarite.geo1_geo2_geo3.30m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2_geo3.30m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1_geo2_geo3.60m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2_geo3.60m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1_geo2_geo3.100m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2_geo3.100m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1_geo2_geo3.250m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2_geo3.250m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1_geo2_geo3.500m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2_geo3.500m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1_geo2_geo3.750m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2_geo3.750m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1_geo2_geo3.1000m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2_geo3.1000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1_geo2_geo3.3000m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2_geo3.3000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.geo1_geo2_geo3.5000m<-SAI(diversite.selection.rarete.complementarite.geo1_geo2_geo3.5000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)



 
############################HABITAT g1_g2_g3_b1###########################################


 ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances


relation.habitats.sites.g1_g2_g3_b1<-read.csv2("fc3_intersect_g1_g2_g3_b1.csv")


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



#Load de la requette "REQ_site_g1_g2_g3_b1_30m.csv"

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

#calcul indice de rareté des habitats

#30m

  b=unique(habitat_par_site_g1_g2_g3_b1_30m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_30m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_g1_g2_g3_b1_30m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_30m$Code_Site[habitat_par_site_g1_g2_g3_b1_30m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_g1_g2_g3_b1_30m$Code_Site[habitat_par_site_g1_g2_g3_b1_30m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_30m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_30m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_30m_station<-matrix(0,500,26)
Indice_rarete_g1_g2_g3_b1_30m_station<-NULL

    for (m in 1:26){
  g1_g2_g3_b1_30m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_30m$Code_hab[habitat_par_site_g1_g2_g3_b1_30m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_30m$Code_hab[habitat_par_site_g1_g2_g3_b1_30m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_30m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_30m[names(Indice_rarete_g1_g2_g3_b1_30m)%in%g1_g2_g3_b1_30m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_30m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_30m=sum(Indice_rarete_g1_g2_g3_b1_30m_station)                #somme des indices de rareté  de ttes les sp


   
#60m

  b=unique(habitat_par_site_g1_g2_g3_b1_60m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_60m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_g1_g2_g3_b1_60m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_60m$Code_Site[habitat_par_site_g1_g2_g3_b1_60m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_g1_g2_g3_b1_60m$Code_Site[habitat_par_site_g1_g2_g3_b1_60m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_60m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_60m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_60m_station<-matrix(0,500,26)
Indice_rarete_g1_g2_g3_b1_60m_station<-NULL

    for (m in 1:26){
  g1_g2_g3_b1_60m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_60m$Code_hab[habitat_par_site_g1_g2_g3_b1_60m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_60m$Code_hab[habitat_par_site_g1_g2_g3_b1_60m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_60m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_60m[names(Indice_rarete_g1_g2_g3_b1_60m)%in%g1_g2_g3_b1_60m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_60m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_60m=sum(Indice_rarete_g1_g2_g3_b1_60m_station)                #somme des indices de rareté  de ttes les sp


    
#100m

  b=unique(habitat_par_site_g1_g2_g3_b1_100m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_100m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_g1_g2_g3_b1_100m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_100m$Code_Site[habitat_par_site_g1_g2_g3_b1_100m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_g1_g2_g3_b1_100m$Code_Site[habitat_par_site_g1_g2_g3_b1_100m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_100m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_100m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_100m_station<-matrix(0,500,26)
Indice_rarete_g1_g2_g3_b1_100m_station<-NULL

    for (m in 1:26){
  g1_g2_g3_b1_100m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_100m$Code_hab[habitat_par_site_g1_g2_g3_b1_100m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_100m$Code_hab[habitat_par_site_g1_g2_g3_b1_100m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_100m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_100m[names(Indice_rarete_g1_g2_g3_b1_100m)%in%g1_g2_g3_b1_100m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_100m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_100m=sum(Indice_rarete_g1_g2_g3_b1_100m_station)                #somme des indices de rareté  de ttes les sp

              
#250m

  b=unique(habitat_par_site_g1_g2_g3_b1_250m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_250m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_g1_g2_g3_b1_250m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_250m$Code_Site[habitat_par_site_g1_g2_g3_b1_250m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_g1_g2_g3_b1_250m$Code_Site[habitat_par_site_g1_g2_g3_b1_250m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_250m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_250m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_250m_station<-matrix(0,500,26)
Indice_rarete_g1_g2_g3_b1_250m_station<-NULL

    for (m in 1:26){
  g1_g2_g3_b1_250m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_250m$Code_hab[habitat_par_site_g1_g2_g3_b1_250m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_250m$Code_hab[habitat_par_site_g1_g2_g3_b1_250m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_250m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_250m[names(Indice_rarete_g1_g2_g3_b1_250m)%in%g1_g2_g3_b1_250m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_250m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_250m=sum(Indice_rarete_g1_g2_g3_b1_250m_station)                #somme des indices de rareté  de ttes les sp

     
#500m

  b=unique(habitat_par_site_g1_g2_g3_b1_500m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_500m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_g1_g2_g3_b1_500m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_500m$Code_Site[habitat_par_site_g1_g2_g3_b1_500m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_g1_g2_g3_b1_500m$Code_Site[habitat_par_site_g1_g2_g3_b1_500m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_500m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_500m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_500m_station<-matrix(0,500,26)
Indice_rarete_g1_g2_g3_b1_500m_station<-NULL

    for (m in 1:26){
  g1_g2_g3_b1_500m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_500m$Code_hab[habitat_par_site_g1_g2_g3_b1_500m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_500m$Code_hab[habitat_par_site_g1_g2_g3_b1_500m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_500m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_500m[names(Indice_rarete_g1_g2_g3_b1_500m)%in%g1_g2_g3_b1_500m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_500m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_500m=sum(Indice_rarete_g1_g2_g3_b1_500m_station)                #somme des indices de rareté  de ttes les sp

     
#750m

  b=unique(habitat_par_site_g1_g2_g3_b1_750m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_750m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_g1_g2_g3_b1_750m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_750m$Code_Site[habitat_par_site_g1_g2_g3_b1_750m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_g1_g2_g3_b1_750m$Code_Site[habitat_par_site_g1_g2_g3_b1_750m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_750m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_750m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_750m_station<-matrix(0,500,26)
Indice_rarete_g1_g2_g3_b1_750m_station<-NULL

    for (m in 1:26){
  g1_g2_g3_b1_750m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_750m$Code_hab[habitat_par_site_g1_g2_g3_b1_750m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_750m$Code_hab[habitat_par_site_g1_g2_g3_b1_750m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_750m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_750m[names(Indice_rarete_g1_g2_g3_b1_750m)%in%g1_g2_g3_b1_750m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_750m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_750m=sum(Indice_rarete_g1_g2_g3_b1_750m_station)                #somme des indices de rareté  de ttes les sp

      
#1000m

  b=unique(habitat_par_site_g1_g2_g3_b1_1000m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_1000m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_g1_g2_g3_b1_1000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_1000m$Code_Site[habitat_par_site_g1_g2_g3_b1_1000m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_g1_g2_g3_b1_1000m$Code_Site[habitat_par_site_g1_g2_g3_b1_1000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_1000m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_1000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_1000m_station<-matrix(0,500,26)
Indice_rarete_g1_g2_g3_b1_1000m_station<-NULL

    for (m in 1:26){
  g1_g2_g3_b1_1000m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_1000m$Code_hab[habitat_par_site_g1_g2_g3_b1_1000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_1000m$Code_hab[habitat_par_site_g1_g2_g3_b1_1000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_1000m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_1000m[names(Indice_rarete_g1_g2_g3_b1_1000m)%in%g1_g2_g3_b1_1000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_1000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_1000m=sum(Indice_rarete_g1_g2_g3_b1_1000m_station)                #somme des indices de rareté  de ttes les sp

                  
#3000m

  b=unique(habitat_par_site_g1_g2_g3_b1_3000m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_3000m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_g1_g2_g3_b1_3000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_3000m$Code_Site[habitat_par_site_g1_g2_g3_b1_3000m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_g1_g2_g3_b1_3000m$Code_Site[habitat_par_site_g1_g2_g3_b1_3000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_3000m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_3000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_3000m_station<-matrix(0,500,26)
Indice_rarete_g1_g2_g3_b1_3000m_station<-NULL

    for (m in 1:26){
  g1_g2_g3_b1_3000m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_3000m$Code_hab[habitat_par_site_g1_g2_g3_b1_3000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_3000m$Code_hab[habitat_par_site_g1_g2_g3_b1_3000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_3000m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_3000m[names(Indice_rarete_g1_g2_g3_b1_3000m)%in%g1_g2_g3_b1_3000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_3000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_3000m=sum(Indice_rarete_g1_g2_g3_b1_3000m_station)                #somme des indices de rareté  de ttes les sp

         
#5000m

  b=unique(habitat_par_site_g1_g2_g3_b1_5000m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_5000m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_g1_g2_g3_b1_5000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_5000m$Code_Site[habitat_par_site_g1_g2_g3_b1_5000m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_g1_g2_g3_b1_5000m$Code_Site[habitat_par_site_g1_g2_g3_b1_5000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_5000m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_5000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_5000m_station<-matrix(0,500,26)
Indice_rarete_g1_g2_g3_b1_5000m_station<-NULL

    for (m in 1:26){
  g1_g2_g3_b1_5000m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_5000m$Code_hab[habitat_par_site_g1_g2_g3_b1_5000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_5000m$Code_hab[habitat_par_site_g1_g2_g3_b1_5000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_5000m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_5000m[names(Indice_rarete_g1_g2_g3_b1_5000m)%in%g1_g2_g3_b1_5000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_5000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_5000m=sum(Indice_rarete_g1_g2_g3_b1_5000m_station)                #somme des indices de rareté  de ttes les sp


source('script.rarete.complementarite.r')

#Selection des dites selon un scenario de rareté complementarité basé sur les habitats

#habitat géo1
selection.rarete.complementarite.g1_g2_g3_b1.30m<-run.bin(habitat_par_site_g1_g2_g3_b1_30m_pres_abs,Indice_rarete_g1_g2_g3_b1_30m_station)
selection.rarete.complementarite.g1_g2_g3_b1.60m<-run.bin(habitat_par_site_g1_g2_g3_b1_60m_pres_abs,Indice_rarete_g1_g2_g3_b1_60m_station)
selection.rarete.complementarite.g1_g2_g3_b1.100m<-run.bin(habitat_par_site_g1_g2_g3_b1_100m_pres_abs,Indice_rarete_g1_g2_g3_b1_100m_station)
selection.rarete.complementarite.g1_g2_g3_b1.250m<-run.bin(habitat_par_site_g1_g2_g3_b1_250m_pres_abs,Indice_rarete_g1_g2_g3_b1_250m_station)
selection.rarete.complementarite.g1_g2_g3_b1.500m<-run.bin(habitat_par_site_g1_g2_g3_b1_500m_pres_abs,Indice_rarete_g1_g2_g3_b1_500m_station)
selection.rarete.complementarite.g1_g2_g3_b1.750m<-run.bin(habitat_par_site_g1_g2_g3_b1_750m_pres_abs,Indice_rarete_g1_g2_g3_b1_750m_station)
selection.rarete.complementarite.g1_g2_g3_b1.1000m<-run.bin(habitat_par_site_g1_g2_g3_b1_1000m_pres_abs,Indice_rarete_g1_g2_g3_b1_1000m_station)
selection.rarete.complementarite.g1_g2_g3_b1.3000m<-run.bin(habitat_par_site_g1_g2_g3_b1_3000m_pres_abs,Indice_rarete_g1_g2_g3_b1_3000m_station)
selection.rarete.complementarite.g1_g2_g3_b1.5000m<-run.bin(habitat_par_site_g1_g2_g3_b1_5000m_pres_abs,Indice_rarete_g1_g2_g3_b1_5000m_station)

#DIVERSITE SELECTIONNE PAR LE SCRIPT RARETE COMPLEMENTARITE#
 source('script.diversite.rarete.selectionnee.r')
 
 diversite.selection.rarete.complementarite.g1_g2_g3_b1.30m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1.30m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1.60m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1.60m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1.100m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1.100m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1.250m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1.250m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1.500m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1.500m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1.750m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1.750m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1.1000m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1.1000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1.3000m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1.3000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1.5000m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1.5000m,Indice_rarete_station)

#CALCUL DES SAI


source('script.sai.r')
SAI_rarete.complementarite.g1_g2_g3_b1.30m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1.30m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1.60m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1.60m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1.100m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1.100m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1.250m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1.250m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1.500m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1.500m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1.750m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1.750m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1.1000m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1.1000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1.3000m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1.3000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1.5000m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1.5000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)



 
############################HABITAT g1_g2_g3_b1_b2###########################################


 ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances


relation.habitats.sites.g1_g2_g3_b1_b2<-read.csv2("fc3_intersect_g1_g2_g3_b1_b2.csv")


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



#Load de la requette "REQ_site_g1_g2_g3_b1_b2_30m.csv"

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

#calcul indice de rareté des habitats

#30m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_30m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_30m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_g1_g2_g3_b1_b2_30m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_30m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_30m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_g1_g2_g3_b1_b2_30m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_30m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_30m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_30m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_30m_station<-matrix(0,500,26)
Indice_rarete_g1_g2_g3_b1_b2_30m_station<-NULL

    for (m in 1:26){
  g1_g2_g3_b1_b2_30m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_30m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_30m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_30m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_30m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_30m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_30m[names(Indice_rarete_g1_g2_g3_b1_b2_30m)%in%g1_g2_g3_b1_b2_30m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_30m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_30m=sum(Indice_rarete_g1_g2_g3_b1_b2_30m_station)                #somme des indices de rareté  de ttes les sp


   
#60m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_60m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_60m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_g1_g2_g3_b1_b2_60m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_60m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_60m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_g1_g2_g3_b1_b2_60m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_60m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_60m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_60m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_60m_station<-matrix(0,500,26)
Indice_rarete_g1_g2_g3_b1_b2_60m_station<-NULL

    for (m in 1:26){
  g1_g2_g3_b1_b2_60m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_60m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_60m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_60m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_60m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_60m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_60m[names(Indice_rarete_g1_g2_g3_b1_b2_60m)%in%g1_g2_g3_b1_b2_60m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_60m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_60m=sum(Indice_rarete_g1_g2_g3_b1_b2_60m_station)                #somme des indices de rareté  de ttes les sp


    
#100m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_100m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_100m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_g1_g2_g3_b1_b2_100m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_100m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_100m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_g1_g2_g3_b1_b2_100m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_100m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_100m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_100m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_100m_station<-matrix(0,500,26)
Indice_rarete_g1_g2_g3_b1_b2_100m_station<-NULL

    for (m in 1:26){
  g1_g2_g3_b1_b2_100m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_100m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_100m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_100m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_100m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_100m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_100m[names(Indice_rarete_g1_g2_g3_b1_b2_100m)%in%g1_g2_g3_b1_b2_100m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_100m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_100m=sum(Indice_rarete_g1_g2_g3_b1_b2_100m_station)                #somme des indices de rareté  de ttes les sp

              
#250m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_250m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_250m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_g1_g2_g3_b1_b2_250m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_250m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_250m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_g1_g2_g3_b1_b2_250m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_250m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_250m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_250m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_250m_station<-matrix(0,500,26)
Indice_rarete_g1_g2_g3_b1_b2_250m_station<-NULL

    for (m in 1:26){
  g1_g2_g3_b1_b2_250m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_250m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_250m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_250m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_250m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_250m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_250m[names(Indice_rarete_g1_g2_g3_b1_b2_250m)%in%g1_g2_g3_b1_b2_250m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_250m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_250m=sum(Indice_rarete_g1_g2_g3_b1_b2_250m_station)                #somme des indices de rareté  de ttes les sp

     
#500m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_500m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_500m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_g1_g2_g3_b1_b2_500m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_500m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_500m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_g1_g2_g3_b1_b2_500m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_500m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_500m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_500m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_500m_station<-matrix(0,500,26)
Indice_rarete_g1_g2_g3_b1_b2_500m_station<-NULL

    for (m in 1:26){
  g1_g2_g3_b1_b2_500m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_500m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_500m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_500m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_500m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_500m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_500m[names(Indice_rarete_g1_g2_g3_b1_b2_500m)%in%g1_g2_g3_b1_b2_500m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_500m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_500m=sum(Indice_rarete_g1_g2_g3_b1_b2_500m_station)                #somme des indices de rareté  de ttes les sp

     
#750m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_750m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_750m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_g1_g2_g3_b1_b2_750m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_750m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_750m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_g1_g2_g3_b1_b2_750m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_750m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_750m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_750m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_750m_station<-matrix(0,500,26)
Indice_rarete_g1_g2_g3_b1_b2_750m_station<-NULL

    for (m in 1:26){
  g1_g2_g3_b1_b2_750m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_750m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_750m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_750m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_750m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_750m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_750m[names(Indice_rarete_g1_g2_g3_b1_b2_750m)%in%g1_g2_g3_b1_b2_750m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_750m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_750m=sum(Indice_rarete_g1_g2_g3_b1_b2_750m_station)                #somme des indices de rareté  de ttes les sp

      
#1000m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_1000m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_1000m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_g1_g2_g3_b1_b2_1000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_1000m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_1000m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_g1_g2_g3_b1_b2_1000m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_1000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_1000m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_1000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_1000m_station<-matrix(0,500,26)
Indice_rarete_g1_g2_g3_b1_b2_1000m_station<-NULL

    for (m in 1:26){
  g1_g2_g3_b1_b2_1000m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_1000m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_1000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_1000m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_1000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_1000m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_1000m[names(Indice_rarete_g1_g2_g3_b1_b2_1000m)%in%g1_g2_g3_b1_b2_1000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_1000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_1000m=sum(Indice_rarete_g1_g2_g3_b1_b2_1000m_station)                #somme des indices de rareté  de ttes les sp

                  
#3000m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_3000m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_3000m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_g1_g2_g3_b1_b2_3000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_3000m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_3000m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_g1_g2_g3_b1_b2_3000m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_3000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_3000m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_3000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_3000m_station<-matrix(0,500,26)
Indice_rarete_g1_g2_g3_b1_b2_3000m_station<-NULL

    for (m in 1:26){
  g1_g2_g3_b1_b2_3000m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_3000m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_3000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_3000m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_3000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_3000m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_3000m[names(Indice_rarete_g1_g2_g3_b1_b2_3000m)%in%g1_g2_g3_b1_b2_3000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_3000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_3000m=sum(Indice_rarete_g1_g2_g3_b1_b2_3000m_station)                #somme des indices de rareté  de ttes les sp

         
#5000m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_5000m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_5000m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_g1_g2_g3_b1_b2_5000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_5000m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_5000m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_g1_g2_g3_b1_b2_5000m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_5000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_5000m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_5000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_5000m_station<-matrix(0,500,26)
Indice_rarete_g1_g2_g3_b1_b2_5000m_station<-NULL

    for (m in 1:26){
  g1_g2_g3_b1_b2_5000m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_5000m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_5000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_5000m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_5000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_5000m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_5000m[names(Indice_rarete_g1_g2_g3_b1_b2_5000m)%in%g1_g2_g3_b1_b2_5000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_5000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_5000m=sum(Indice_rarete_g1_g2_g3_b1_b2_5000m_station)                #somme des indices de rareté  de ttes les sp


source('script.rarete.complementarite.r')

#Selection des dites selon un scenario de rareté complementarité basé sur les habitats

#habitat géo1
selection.rarete.complementarite.g1_g2_g3_b1_b2.30m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_30m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_30m_station)
selection.rarete.complementarite.g1_g2_g3_b1_b2.60m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_60m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_60m_station)
selection.rarete.complementarite.g1_g2_g3_b1_b2.100m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_100m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_100m_station)
selection.rarete.complementarite.g1_g2_g3_b1_b2.250m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_250m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_250m_station)
selection.rarete.complementarite.g1_g2_g3_b1_b2.500m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_500m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_500m_station)
selection.rarete.complementarite.g1_g2_g3_b1_b2.750m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_750m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_750m_station)
selection.rarete.complementarite.g1_g2_g3_b1_b2.1000m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_1000m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_1000m_station)
selection.rarete.complementarite.g1_g2_g3_b1_b2.3000m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_3000m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_3000m_station)
selection.rarete.complementarite.g1_g2_g3_b1_b2.5000m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_5000m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_5000m_station)

#DIVERSITE SELECTIONNE PAR LE SCRIPT RARETE COMPLEMENTARITE#
 source('script.diversite.rarete.selectionnee.r')
 
 diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.30m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2.30m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.60m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2.60m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.100m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2.100m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.250m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2.250m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.500m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2.500m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.750m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2.750m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.1000m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2.1000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.3000m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2.3000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.5000m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2.5000m,Indice_rarete_station)

#CALCUL DES SAI


source('script.sai.r')
SAI_rarete.complementarite.g1_g2_g3_b1_b2.30m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.30m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1_b2.60m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.60m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1_b2.100m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.100m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1_b2.250m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.250m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1_b2.500m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.500m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1_b2.750m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.750m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1_b2.1000m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.1000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1_b2.3000m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.3000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1_b2.5000m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2.5000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)



 
############################HABITAT g1_g2_g3_b1_b2_b3###########################################


 ############################CREATION D'UN FICHIER PAR DISTANCE#######################################

#load des fichiers GIS contenant les surfaces de chaque habitats, mais sans distinction des distances


relation.habitats.sites.g1_g2_g3_b1_b2_b3<-read.csv2("fc3_intersect_g1_g2_g3_b1_b2_b3.csv")


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



#Load de la requette "REQ_site_g1_g2_g3_b1_b2_b3_30m.csv"

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

#calcul indice de rareté des habitats

#30m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_b3_30m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_b3_30m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_g1_g2_g3_b1_b2_b3_30m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_30m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_30m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_g1_g2_g3_b1_b2_b3_30m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_30m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_b3_30m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_b3_30m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_b3_30m_station<-matrix(0,500,26)
Indice_rarete_g1_g2_g3_b1_b2_b3_30m_station<-NULL

    for (m in 1:26){
  g1_g2_g3_b1_b2_b3_30m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_30m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_30m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_b3_30m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_30m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_b3_30m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_30m[names(Indice_rarete_g1_g2_g3_b1_b2_b3_30m)%in%g1_g2_g3_b1_b2_b3_30m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_b3_30m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_b3_30m=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_30m_station)                #somme des indices de rareté  de ttes les sp


   
#60m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_b3_60m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_b3_60m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_g1_g2_g3_b1_b2_b3_60m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_60m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_60m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_g1_g2_g3_b1_b2_b3_60m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_60m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_b3_60m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_b3_60m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_b3_60m_station<-matrix(0,500,26)
Indice_rarete_g1_g2_g3_b1_b2_b3_60m_station<-NULL

    for (m in 1:26){
  g1_g2_g3_b1_b2_b3_60m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_60m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_60m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_b3_60m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_60m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_b3_60m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_60m[names(Indice_rarete_g1_g2_g3_b1_b2_b3_60m)%in%g1_g2_g3_b1_b2_b3_60m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_b3_60m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_b3_60m=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_60m_station)                #somme des indices de rareté  de ttes les sp


    
#100m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_b3_100m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_b3_100m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_g1_g2_g3_b1_b2_b3_100m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_100m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_100m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_g1_g2_g3_b1_b2_b3_100m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_100m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_b3_100m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_b3_100m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_b3_100m_station<-matrix(0,500,26)
Indice_rarete_g1_g2_g3_b1_b2_b3_100m_station<-NULL

    for (m in 1:26){
  g1_g2_g3_b1_b2_b3_100m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_100m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_100m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_b3_100m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_100m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_b3_100m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_100m[names(Indice_rarete_g1_g2_g3_b1_b2_b3_100m)%in%g1_g2_g3_b1_b2_b3_100m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_b3_100m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_b3_100m=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_100m_station)                #somme des indices de rareté  de ttes les sp

              
#250m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_b3_250m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_b3_250m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_g1_g2_g3_b1_b2_b3_250m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_250m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_250m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_g1_g2_g3_b1_b2_b3_250m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_250m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_b3_250m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_b3_250m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_b3_250m_station<-matrix(0,500,26)
Indice_rarete_g1_g2_g3_b1_b2_b3_250m_station<-NULL

    for (m in 1:26){
  g1_g2_g3_b1_b2_b3_250m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_250m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_250m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_b3_250m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_250m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_b3_250m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_250m[names(Indice_rarete_g1_g2_g3_b1_b2_b3_250m)%in%g1_g2_g3_b1_b2_b3_250m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_b3_250m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_b3_250m=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_250m_station)                #somme des indices de rareté  de ttes les sp

     
#500m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_b3_500m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_b3_500m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_g1_g2_g3_b1_b2_b3_500m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_500m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_500m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_g1_g2_g3_b1_b2_b3_500m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_500m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_b3_500m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_b3_500m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_b3_500m_station<-matrix(0,500,26)
Indice_rarete_g1_g2_g3_b1_b2_b3_500m_station<-NULL

    for (m in 1:26){
  g1_g2_g3_b1_b2_b3_500m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_500m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_500m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_b3_500m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_500m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_b3_500m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_500m[names(Indice_rarete_g1_g2_g3_b1_b2_b3_500m)%in%g1_g2_g3_b1_b2_b3_500m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_b3_500m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_b3_500m=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_500m_station)                #somme des indices de rareté  de ttes les sp

     
#750m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_b3_750m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_b3_750m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_g1_g2_g3_b1_b2_b3_750m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_750m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_750m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_g1_g2_g3_b1_b2_b3_750m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_750m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_b3_750m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_b3_750m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_b3_750m_station<-matrix(0,500,26)
Indice_rarete_g1_g2_g3_b1_b2_b3_750m_station<-NULL

    for (m in 1:26){
  g1_g2_g3_b1_b2_b3_750m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_750m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_750m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_b3_750m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_750m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_b3_750m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_750m[names(Indice_rarete_g1_g2_g3_b1_b2_b3_750m)%in%g1_g2_g3_b1_b2_b3_750m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_b3_750m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_b3_750m=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_750m_station)                #somme des indices de rareté  de ttes les sp

      
#1000m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_b3_1000m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_b3_1000m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_g1_g2_g3_b1_b2_b3_1000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_1000m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_1000m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_g1_g2_g3_b1_b2_b3_1000m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_1000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_b3_1000m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_b3_1000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_b3_1000m_station<-matrix(0,500,26)
Indice_rarete_g1_g2_g3_b1_b2_b3_1000m_station<-NULL

    for (m in 1:26){
  g1_g2_g3_b1_b2_b3_1000m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_1000m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_1000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_b3_1000m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_1000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_b3_1000m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_1000m[names(Indice_rarete_g1_g2_g3_b1_b2_b3_1000m)%in%g1_g2_g3_b1_b2_b3_1000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_b3_1000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_b3_1000m=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_1000m_station)                #somme des indices de rareté  de ttes les sp

                  
#3000m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_b3_3000m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_b3_3000m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_g1_g2_g3_b1_b2_b3_3000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_3000m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_3000m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_g1_g2_g3_b1_b2_b3_3000m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_3000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_b3_3000m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_b3_3000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_b3_3000m_station<-matrix(0,500,26)
Indice_rarete_g1_g2_g3_b1_b2_b3_3000m_station<-NULL

    for (m in 1:26){
  g1_g2_g3_b1_b2_b3_3000m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_3000m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_3000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_b3_3000m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_3000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_b3_3000m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_3000m[names(Indice_rarete_g1_g2_g3_b1_b2_b3_3000m)%in%g1_g2_g3_b1_b2_b3_3000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_b3_3000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_b3_3000m=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_3000m_station)                #somme des indices de rareté  de ttes les sp

         
#5000m

  b=unique(habitat_par_site_g1_g2_g3_b1_b2_b3_5000m$Code_hab)

    Indice_rarete_g1_g2_g3_b1_b2_b3_5000m=NULL
  site_habitat<-matrix(0,26,length(habitat_par_site_g1_g2_g3_b1_b2_b3_5000m$Code_hab))

    
    for (n in 1:length(b)){ 
    site_habitat[,n]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_5000m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_5000m$Code_hab==b[n]]),rep(0,26-length(habitat_par_site_g1_g2_g3_b1_b2_b3_5000m$Code_Site[habitat_par_site_g1_g2_g3_b1_b2_b3_5000m$Code_hab==b[n]])))
    names(site_habitat[,n])=b[n]
    Indice_rarete_g1_g2_g3_b1_b2_b3_5000m[n]<-2^(26-length(site_habitat[,n][site_habitat[,n]!=0]))
    names(Indice_rarete_g1_g2_g3_b1_b2_b3_5000m)[n]=as.vector(b[n])
    }
    
#calcul indice de rareté des stations    

g1_g2_g3_b1_b2_b3_5000m_station<-matrix(0,500,26)
Indice_rarete_g1_g2_g3_b1_b2_b3_5000m_station<-NULL

    for (m in 1:26){
  g1_g2_g3_b1_b2_b3_5000m_station[,m]<-c(as.vector(habitat_par_site_g1_g2_g3_b1_b2_b3_5000m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_5000m$Code_Site==sites.etudies[m]]),rep(0,500-length(habitat_par_site_g1_g2_g3_b1_b2_b3_5000m$Code_hab[habitat_par_site_g1_g2_g3_b1_b2_b3_5000m$Code_Site==sites.etudies[m]])))
  Indice_rarete_g1_g2_g3_b1_b2_b3_5000m_station[m]=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_5000m[names(Indice_rarete_g1_g2_g3_b1_b2_b3_5000m)%in%g1_g2_g3_b1_b2_b3_5000m_station[,m]])
  }
   
                                     
  names(Indice_rarete_g1_g2_g3_b1_b2_b3_5000m_station)<-sites.etudies
  
   #indice de rarete stations contiens les indices de rareté des stations (ah ben _ça alors !)
    tot_g1_g2_g3_b1_b2_b3_5000m=sum(Indice_rarete_g1_g2_g3_b1_b2_b3_5000m_station)                #somme des indices de rareté  de ttes les sp


source('script.rarete.complementarite.r')

#Selection des dites selon un scenario de rareté complementarité basé sur les habitats

#habitat géo1
selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.30m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_b3_30m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_b3_30m_station)
selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.60m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_b3_60m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_b3_60m_station)
selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.100m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_b3_100m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_b3_100m_station)
selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.250m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_b3_250m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_b3_250m_station)
selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.500m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_b3_500m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_b3_500m_station)
selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.750m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_b3_750m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_b3_750m_station)
selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.1000m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_b3_1000m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_b3_1000m_station)
selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.3000m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_b3_3000m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_b3_3000m_station)
selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.5000m<-run.bin(habitat_par_site_g1_g2_g3_b1_b2_b3_5000m_pres_abs,Indice_rarete_g1_g2_g3_b1_b2_b3_5000m_station)

#DIVERSITE SELECTIONNE PAR LE SCRIPT RARETE COMPLEMENTARITE#
 source('script.diversite.rarete.selectionnee.r')
 
 diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.30m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.30m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.60m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.60m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.100m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.100m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.250m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.250m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.500m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.500m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.750m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.750m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.1000m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.1000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.3000m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.3000m,Indice_rarete_station)
diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.5000m<-species.rarete.evol(selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.5000m,Indice_rarete_station)

#CALCUL DES SAI


source('script.sai.r')
SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.30m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.30m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.60m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.60m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.100m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.100m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.250m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.250m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.500m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.500m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.750m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.750m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.1000m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.1000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.3000m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.3000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)
SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.5000m<-SAI(diversite.selection.rarete.complementarite.g1_g2_g3_b1_b2_b3.5000m,diversite.rarete.moyenne.selection.aleatoire,diversite.selection.rarete.complementarite.poisson)



#STOCKAGE DES RESULTATS

targeted.station.results[,p+1]=c(SAI_rarete.complementarite.geo1.30m,SAI_rarete.complementarite.geo1.60m,SAI_rarete.complementarite.geo1.100m,SAI_rarete.complementarite.geo1.250m,
SAI_rarete.complementarite.geo1.500m,SAI_rarete.complementarite.geo1.750m,SAI_rarete.complementarite.geo1.1000m,SAI_rarete.complementarite.geo1.3000m,SAI_rarete.complementarite.geo1.5000m,
SAI_rarete.complementarite.geo2.30m,SAI_rarete.complementarite.geo2.60m,SAI_rarete.complementarite.geo2.100m,SAI_rarete.complementarite.geo2.250m,
SAI_rarete.complementarite.geo2.500m,SAI_rarete.complementarite.geo2.750m,SAI_rarete.complementarite.geo2.1000m,SAI_rarete.complementarite.geo2.3000m,SAI_rarete.complementarite.geo2.5000m,
SAI_rarete.complementarite.geo3.100m,SAI_rarete.complementarite.geo3.250m,
SAI_rarete.complementarite.geo3.500m,SAI_rarete.complementarite.geo3.750m,SAI_rarete.complementarite.geo3.1000m,SAI_rarete.complementarite.geo3.3000m,SAI_rarete.complementarite.geo3.5000m,
SAI_rarete.complementarite.bent1.30m,SAI_rarete.complementarite.bent1.60m,SAI_rarete.complementarite.bent1.100m,SAI_rarete.complementarite.bent1.250m,
SAI_rarete.complementarite.bent1.500m,SAI_rarete.complementarite.bent1.750m,SAI_rarete.complementarite.bent1.1000m,SAI_rarete.complementarite.bent1.3000m,SAI_rarete.complementarite.bent1.5000m,
SAI_rarete.complementarite.bent2.30m,SAI_rarete.complementarite.bent2.60m,SAI_rarete.complementarite.bent2.100m,SAI_rarete.complementarite.bent2.250m,
SAI_rarete.complementarite.bent2.500m,SAI_rarete.complementarite.bent2.750m,SAI_rarete.complementarite.bent2.1000m,SAI_rarete.complementarite.bent2.3000m,SAI_rarete.complementarite.bent2.5000m,
SAI_rarete.complementarite.bent3.30m,SAI_rarete.complementarite.bent3.60m,SAI_rarete.complementarite.bent3.100m,SAI_rarete.complementarite.bent3.250m,
SAI_rarete.complementarite.bent3.500m,SAI_rarete.complementarite.bent3.750m,SAI_rarete.complementarite.bent3.1000m,SAI_rarete.complementarite.bent3.3000m,SAI_rarete.complementarite.bent3.5000m,
SAI_rarete.complementarite.geo1_geo2.30m,SAI_rarete.complementarite.geo1_geo2.60m,SAI_rarete.complementarite.geo1_geo2.100m,SAI_rarete.complementarite.geo1_geo2.250m,
SAI_rarete.complementarite.geo1_geo2.500m,SAI_rarete.complementarite.geo1_geo2.750m,SAI_rarete.complementarite.geo1_geo2.1000m,SAI_rarete.complementarite.geo1_geo2.3000m,SAI_rarete.complementarite.geo1_geo2.5000m,
SAI_rarete.complementarite.geo1_geo2_geo3.30m,SAI_rarete.complementarite.geo1_geo2_geo3.60m,SAI_rarete.complementarite.geo1_geo2_geo3.100m,SAI_rarete.complementarite.geo1_geo2_geo3.250m,
SAI_rarete.complementarite.geo1_geo2_geo3.500m,SAI_rarete.complementarite.geo1_geo2_geo3.750m,SAI_rarete.complementarite.geo1_geo2_geo3.1000m,SAI_rarete.complementarite.geo1_geo2_geo3.3000m,SAI_rarete.complementarite.geo1_geo2_geo3.5000m,
SAI_rarete.complementarite.g1_g2_g3_b1.30m,SAI_rarete.complementarite.g1_g2_g3_b1.60m,SAI_rarete.complementarite.g1_g2_g3_b1.100m,SAI_rarete.complementarite.g1_g2_g3_b1.250m,
SAI_rarete.complementarite.g1_g2_g3_b1.500m,SAI_rarete.complementarite.g1_g2_g3_b1.750m,SAI_rarete.complementarite.g1_g2_g3_b1.1000m,SAI_rarete.complementarite.g1_g2_g3_b1.3000m,SAI_rarete.complementarite.g1_g2_g3_b1.5000m,
SAI_rarete.complementarite.g1_g2_g3_b1_b2.30m,SAI_rarete.complementarite.g1_g2_g3_b1_b2.60m,SAI_rarete.complementarite.g1_g2_g3_b1_b2.100m,SAI_rarete.complementarite.g1_g2_g3_b1_b2.250m,
SAI_rarete.complementarite.g1_g2_g3_b1_b2.500m,SAI_rarete.complementarite.g1_g2_g3_b1_b2.750m,SAI_rarete.complementarite.g1_g2_g3_b1_b2.1000m,SAI_rarete.complementarite.g1_g2_g3_b1_b2.3000m,SAI_rarete.complementarite.g1_g2_g3_b1_b2.5000m,
SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.30m,SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.60m,SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.100m,SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.250m,
SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.500m,SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.750m,SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.1000m,SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.3000m,SAI_rarete.complementarite.g1_g2_g3_b1_b2_b3.5000m)
}

 targeted.station.results
