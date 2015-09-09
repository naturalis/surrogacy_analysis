################################################################################



##########################
##   FONCTION RUN.RICH  ##
##########################

#run.rich<-function(table.site,first.station="T_S265",complementarity=T){
   run.rich<-function(table.site,first.station=NULL,complementarity=T){
                   # notre fonction run.rich est fonction de notre tableau "table.site" (tableau de presence absence 
                   #des taxon par site); 
                   #de notre premiere station selectionnée; et de si on utilise la complementarité ou pas dans l'algorithme (on a le choix)

## Conditions initiales

#########################################################################################
if(is.null(first.station))
{
unselected.sites=rownames(table.site) # Initialisation des sites non-selectionnés : tous les sites du tableau initial au départ
selected.sites=NULL # Initialisation de la liste de sites selectionnés : NULL au depart 
habitats.counts=rep(0,dim(table.site)[2]) # Compte le nombre de fois ou chaque habitat est inclu dans le réseau
}

else
{
unselected.sites=rownames(table.site)[rownames(table.site)!=first.station] # Initialisation des sites non-selectionnés : tous les sites du tableau initial
selected.sites=first.station # Initialisation de la liste de sites selectionnés : NULL
habitats.counts=table.site[first.station,] # Compte le nombre de fois ou chaque espèce (habitat) est présente dans le réseau de sites sélectionnés
}

names(habitats.counts)=colnames(table.site) # Noms des colonnes = noms des habitats
#########################################################################################
first.time=T # Booléen pour noter la permière fois ou tous les habitats sont représentés


while(length(unselected.sites)>1)#tant qu'il reste plus d'un site
{

## RICHNESS

habitats.richness.unselected=apply(t(table.site[unselected.sites,]),1,sum) # Comptage du nb d'occurences pr chq habitat (ou gpe fctionnel ou esp) par site
habitats.unselected.notnull=habitats.richness.unselected[habitats.richness.unselected!=0] # Ce n'est pas parce qu'un site n'est pas sélectionné qu'il ne le sera pas par la suite (not null)

## COMPLEMENTARITY

# création d'un tableau ordonné pour le calcul des scores
habitats.richness.selected=sort(habitats.counts[names(habitats.unselected.notnull)]) # Comptages du nb d'occurences pr chq habitats dans les sites déja sélectionnés
# + on supprime les espèces (habitats) non présents dans les sites restants!
# classement par ordre de richesse

table.unselected.sites=table.site[unselected.sites,names(habitats.richness.selected)]

# Creation d'un score de "RICHESSE/COMPLEMENTARITE" pour chaque habitats
# (chaque site donne une ligne de 0 et de 1, dont la somme représente le nombre d'habitats différents par site)
sites.scores=transform(apply(as.matrix(table.unselected.sites),1,function(x) convert2score.rich(x,habitats.richness.selected,complementarity)))

if(dim(sites.scores)[2]==1){
sites.scores.vec=sites.scores[,1]
names(sites.scores.vec)=rownames(sites.scores)
sites.scores=sort(sites.scores.vec)
the.site=names(sites.scores)[length(sites.scores)] # On prend le site qui a le score maximum
}

else{
sites.scores=sites.scores[,do.call(order,transform(t(sites.scores)))] # on ordonne les sites selon score de 0 (plus complémentaires), puis de 1, 2, 3... (moins complémentaires)
# The.site = le score total le plus élevé (0 + 1 + 2 + ...)
the.site=names(sites.scores)[dim(sites.scores)[2]] # Repérage du site ayant le score maximum
}

selected.sites=c(selected.sites,the.site) # Dans la liste des sites sélectionnés, ajout de ce site

# Dans la liste des habitats sélectionnés, chaque représentation d'un habitat inclu dans le réseau est incrémentatée de 1
habitats.counts=habitats.counts+table.site[the.site,] # Dans la liste des habitats sélectionnés, ajout des habitats appartenant au voisinage du site sélectionné pour la richesse en habitats
unselected.sites=unselected.sites[unselected.sites!=the.site] # Dans la liste des sites non sélectionnés, suppression du site ayant le score maximum

} #fin du while : tant qu'il reste des sites à sélectionner


## Ajout du dernier site restant

selected.sites=c(selected.sites,unselected.sites)

# selected.sites
#habitats.representation(table.site,selected.sites)

} # Fin de la fonction run.rich




########################
##   SOUS FONCTION    ##
########################

## convert2score : chaque ligne de présence/absence d'habitat par site est
## convertie en un nombre binaires. Ces nombres binaires sont ensuite convertis
## en nombre décimaux pour la sélection des sites les plus rares. L'habitat qui
## a le score le plus élevé est le plus rare.

convert2score.rich<-function(vec,habitats.richness.selected,complementarity){
if(complementarity==T){
tapply(as.numeric(vec),habitats.richness.selected,sum)}
else{
sum(as.numeric(vec))}
}


################################################################################