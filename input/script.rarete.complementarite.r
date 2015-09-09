run.bin<-function(tablo,Indice_rarete_stations,first.station=NULL){


#initialisation

if(is.null(first.station))
{
unselected.sites=rownames(tablo) # Initialisation des sites non-selectionnés : tous les sites du tableau initial
selected.sites=NULL # Initialisation de la liste de sites selectionnés : NULL
habitat.imaginaire=rep(0,dim(tablo)[2])
}

 else
{
  unselected.sites=rownames(tablo)[rownames(tablo)!=first.station] # Initialisation des sites non-selectionnés : tous les sites du tableau initial
selected.sites=first.station # Initialisation de la liste de sites selectionnés : NULL
habitat.imaginaire=rep(0,dim(tablo)[2])

}


#tant qu'il reste plus d'un site
while(length(unselected.sites)>1)
{

    Indice_rarete_pool_selectionne<-rep(0,length(unselected.sites))
for (i in 1:length(unselected.sites)){
    imaginaire=c(selected.sites,unselected.sites[i])        #on simule un vecteur site qui contient les sites selectionné plus chaque site non selectioné à tour de role

#création d'une matrice tablo.selection équivalente à tablo, mais que pour les sites selectionné ds imaginaires.

tablo.selection<-tablo[rownames(tablo)%in%imaginaire,]

#calcul de l'indice de rareté du pool de site selectionné

Indice_rarete_pool_selectionne[i]<-sum(Indice_rarete_stations[names(Indice_rarete_stations)%in%imaginaire])  #A améliorer pr prendre en compte la complementarite

}


    #on obtien : indice.pool.selected.imaginaire =[Idavecsite1; Idavecsite2....]

    imaxi=1
    for (l in 2:length(Indice_rarete_pool_selectionne)){
        if (Indice_rarete_pool_selectionne[l]>Indice_rarete_pool_selectionne[imaxi]){
          imaxi=l}
          }
    #on obtien : imaxi = i pour lequel indice.pool.selected.imaginaire est max

    bestsite=unselected.sites[imaxi]
    selected.sites=c(selected.sites,bestsite)
    unselected.sites=unselected.sites[unselected.sites!=bestsite]
    } #fin du while :

#Ajout du dernier site restant
selected.sites=c(selected.sites,unselected.sites)
 selected.sites
    } #fin de la fonction

#stations.list=shannon(tablo)
