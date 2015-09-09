#a partir de la table access : surface_geo1_30m

#repertoire<-"C:/Documents and Settings/simon/Bureau/SIM2011/R/"
#setwd(repertoire)
#input<-"INPUT/"
#output<-"OUTPUT/"
#setwd(input)
#attribut=read.csv2("shannon_surface_geo1_30m.csv")
#attribut.unique<-aggregate(attribut$surface,list(Code_Site=attribut$Code_Site,geomorpho1=attribut$geomorpho1),sum,na.rm=TRUE)
#colnames(attribut.unique)=c("Code_Site","geomorpho1","surface")
 #IDEE : FAIR EUN TABLO QU'ON A PAS AGGREGé PAR SURFACE, LUI APPLIQUER UNIQUE POUR VIRER LES POLYGONES EN PLUSIEURS EXEMPLAIRES, et S'EN SERVIR PR CALCULER LA SURFACE TOTALE.
 # FAIRE UN AUTRE TABLO OU ON AGGREGE PAR SURFACE, PUIS UNIQUE (EVENTUELLMT SI BESOIN) COOL POUR SHANNON

#attribut.surface<-attribut$surface
#attribut.surface.unique<-unique(attribut.surface)
#habitat.surface.totale<-sum(attribut.surface.unique)


#transformation de la requette pour avoir tableau coool pr shannon
#tablo<-xtabs(formula=attribut.unique$surface~attribut.unique$Code_Site+attribut.unique$geomorpho1,attribut.unique)


#definition de la fonction
shannon<-function(tablo,first.station=NULL){

#transforme les 0 en 0.001 pour pas qu'ils fassent chier
for(i in 1:dim(tablo)[1]){
  for (j in 1:dim(tablo)[2]){
    if (tablo[i,j]==0){
    tablo[i,j]<-0.001}
    }}
              
#initialisation

if(is.null(first.station))
{
unselected.sites=rownames(tablo) # Initialisation des sites non-selectionnés : tous les sites du tableau initial
selected.sites=NULL # Initialisation de la liste de sites selectionnés : NULL
habitat.surface.imaginaire=rep(0,dim(tablo)[2]) 
}

 else
{
  unselected.sites=rownames(tablo)[rownames(tablo)!=first.station] # Initialisation des sites non-selectionnés : tous les sites du tableau initial
selected.sites=first.station # Initialisation de la liste de sites selectionnés : NULL
habitat.surface.imaginaire=rep(0,dim(tablo)[2])     

}

       

#tant qu'il reste plus d'un site
while(length(unselected.sites)>1)
{
     H=rep(0,length(unselected.sites))
     
for (i in 1:length(unselected.sites)){
    imaginaire=c(selected.sites,unselected.sites[i])        #on simule un vecteur site qui contient les sites selectionné plus chaque site non selectioné à tour de role
    
#création d'une matrice tablo.selection équivalente à tablo, mais que pour les sites selectionné ds imaginaires. 

tablo.selection<-tablo[rownames(tablo)%in%imaginaire,]

#création d'un vecteur contenant toute les surface differentes uniquement, qui sont ds tablo.selection
if (length(c(selected.sites,1))==1){
  surface.presente=NULL
  for (n in 1:length(tablo.selection)){
    if (tablo.selection[n]%in%surface.presente!=TRUE){
      surface.presente=c(surface.presente,tablo.selection[n])}
      }
}
else{

surface.presente=NULL
for (k in 1:dim(tablo.selection)[2]){    #mettre dim au lieu de length si on a mis une statino au départ
  for (l in 1:dim(tablo.selection)[1]){
    if (tablo.selection[l,k]%in%surface.presente!=TRUE){
       surface.presente=c(surface.presente,tablo.selection[l,k])
       }
       }}
}        
#Calcul de la surface totale d'habitat selectionnés

 habitat.surface.totale<-sum(surface.presente)
  
#on calcul shannon
     
    for (j in 1:dim(tablo)[2]){ 
        habitat.surface.imaginaire[j]=sum(tablo[,j][rownames(tablo)%in%imaginaire])    #vecteur qui contient la somme des surfaces de chaques habitat contenu dans le pool de site imagiaire
        #habitat.surface.imaginaire[j]=aggregate(tablo.imaginaire$hab1,list(site=tablo.imaginaire$Code_Site,habitat=tablo.imaginaire$code_geo2),sum,na.rm=TRUE) #habitat surface imaginaire c koi ? un vecteur? une matrice?
#fait la somm des valeur de surface de chaq habitat pour les sites qui st présent ds imaginaire
#j commence à 2 parcque la colonne 1 c'est les sites
        }
    
        
        habitat.surface.percent=habitat.surface.imaginaire/habitat.surface.totale
        habitat.surface.percent.notnull=habitat.surface.percent[!habitat.surface.percent==0]   #on vire les 0 pour calculer shannon car log(0) existe pas
    #on obtient : habitat.surface.percent=[phab1, phab2...]
    a=log(habitat.surface.percent.notnull,2) #on obtient : a=[log2(phabi1); log2(phab2)...]
    b=habitat.surface.percent.notnull*a  #on obtien b=[phab1*log2(phab1); phab2*log2(phab2)....]
    H[i]=-sum(b)
    if (H[i]==0){
    H[i]<-0.001}
    }
    #on obtien : H =[Havecsite1; Havecsite2....]
    
    imaxi=1
    for (i in 2:length(H)){
        if (H[i]>H[imaxi]){
          imaxi=i}
          }
    #on obtien : imaxi = i pour lequel H est max
    
    bestsite=unselected.sites[imaxi]
    selected.sites=c(selected.sites,bestsite)
    unselected.sites=unselected.sites[unselected.sites!=bestsite]
    } #fin du while :
    
#Ajout du dernier site restant
selected.sites=c(selected.sites,unselected.sites)
 selected.sites   
    } #fin de la fonction
    
#stations.list=shannon(tablo)        