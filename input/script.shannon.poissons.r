#En entrée, avoir talo d'abondance (style presence absence, mais avec abondances) avec site, et code poisson.

#definition de la fonction
shannon.poisson<-function(tablo,first.station=NULL){

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
nombre.poissons.imaginaire=rep(0,dim(tablo)[2])
}

 else
{
unselected.sites=rownames(tablo)[rownames(tablo)!=first.station] # Initialisation des sites non-selectionnés : tous les sites du tableau initial
selected.sites=first.station # Initialisation de la liste de sites selectionnés : NULL
nombre.poissons.imaginaire=rep(0,dim(tablo)[2])

}


#tant qu'il reste plus d'un site
while(length(unselected.sites)>1)
{
     H=rep(0,length(unselected.sites))

for (i in 1:length(unselected.sites)){
    imaginaire=c(selected.sites,unselected.sites[i])        #on simule un vecteur site qui contient les sites selectionné plus chaque site non selectioné à tour de role


    for (j in 1:dim(tablo)[2]){
        nombre.poissons.imaginaire[j]=sum(tablo[,j][rownames(tablo)%in%imaginaire])    #vecteur qui contient la somme des surfaces de chaques habitat contenu dans le pool de site imagiaire
        #habitat.surface.imaginaire[j]=aggregate(tablo.imaginaire$hab1,list(site=tablo.imaginaire$Code_Site,habitat=tablo.imaginaire$code_geo2),sum,na.rm=TRUE) #habitat surface imaginaire c koi ? un vecteur? une matrice?
#fait la somm des valeur de surface de chaq habitat pour les sites qui st présent ds imaginaire
#j commence à 2 parcque la colonne 1 c'est les sites
        }
        nombre.total.poissons.imaginaire=sum(nombre.poissons.imaginaire)

        nombre.de.poissons.percent=nombre.poissons.imaginaire/nombre.total.poissons.imaginaire
        nombre.de.poissons.percent.notnull=nombre.de.poissons.percent[!nombre.de.poissons.percent==0]   #on vire les 0 pour calculer shannon car log(0) existe pas
    #on obtient : habitat.surface.percent=[phab1, phab2...]
    a=log(nombre.de.poissons.percent.notnull,2) #on obtient : a=[log2(phabi1); log2(phab2)...]
    b=nombre.de.poissons.percent.notnull*a  #on obtien b=[phab1*log2(phab1); phab2*log2(phab2)....]
    H[i]=-sum(b)
     if (H[i]==0){
    H[i]<-0.001
    }
    }
    
    #on obtien : H =[Havecsite1; Havecsite2....]

    imaxi=1
    for (i in 2:length(H)){
        if (H[i]>H[imaxi]){
          imaxi=i
          }
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