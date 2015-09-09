################################################################################



##########################
##   FONCTION RUN.RICH  ##
##########################

#run.rich<-function(table.site,first.station="T_S265",complementarity=T){
   run.rich<-function(table.site,first.station=NULL,complementarity=T){
                   # notre fonction run.rich est fonction de notre tableau "table.site" (tableau de presence absence 
                   #des taxon par site); 

     
## Conditions initiales

#########################################################################################
if(is.null(first.station))
{
unselected.sites=rownames(table.site)
selected.sites=NULL
habitats.counts=rep(0,dim(table.site)[2])
}

else
{
unselected.sites=rownames(table.site)[rownames(table.site)!=first.station]
selected.sites=first.station # Initialisation de la liste de sites selectionn?s : NULL
habitats.counts=table.site[first.station,] # Compte le nombre de fois ou chaque esp?ce (habitat) est pr?sente dans le r?seau de sites s?lectionn?s
}

names(habitats.counts)=colnames(table.site) # Noms des colonnes = noms des habitats
#########################################################################################
first.time=T # Bool?en pour noter la permi?re fois ou tous les habitats sont repr?sent?s


while(length(unselected.sites)>1)#tant qu'il reste plus d'un site
{

## RICHNESS

habitats.richness.unselected=apply(t(table.site[unselected.sites,]),1,sum) # Comptage du nb d'occurences pr chq habitat (ou gpe fctionnel ou esp) par site
habitats.unselected.notnull=habitats.richness.unselected[habitats.richness.unselected!=0] # Ce n'est pas parce qu'un site n'est pas s?lectionn? qu'il ne le sera pas par la suite (not null)

## COMPLEMENTARITY

# cr?ation d'un tableau ordonn? pour le calcul des scores
habitats.richness.selected=sort(habitats.counts[names(habitats.unselected.notnull)]) # Comptages du nb d'occurences pr chq habitats dans les sites d?ja s?lectionn?s
# + on supprime les esp?ces (habitats) non pr?sents dans les sites restants!
# classement par ordre de richesse

table.unselected.sites=table.site[unselected.sites,names(habitats.richness.selected)]

# Creation d'un score de "RICHESSE/COMPLEMENTARITE" pour chaque habitats
# (chaque site donne une ligne de 0 et de 1, dont la somme repr?sente le nombre d'habitats diff?rents par site)
sites.scores=transform(apply(as.matrix(table.unselected.sites),1,function(x) convert2score.rich(x,habitats.richness.selected,complementarity)))

if(dim(sites.scores)[2]==1){
sites.scores.vec=sites.scores[,1]
names(sites.scores.vec)=rownames(sites.scores)
sites.scores=sort(sites.scores.vec)
the.site=names(sites.scores)[length(sites.scores)] # On prend le site qui a le score maximum
}

else{
sites.scores=sites.scores[,do.call(order,transform(t(sites.scores)))] # on ordonne les sites selon score de 0 (plus compl?mentaires), puis de 1, 2, 3... (moins compl?mentaires)
# The.site = le score total le plus ?lev? (0 + 1 + 2 + ...)
the.site=names(sites.scores)[dim(sites.scores)[2]] # Rep?rage du site ayant le score maximum
}

selected.sites=c(selected.sites,the.site) # Dans la liste des sites s?lectionn?s, ajout de ce site

# Dans la liste des habitats s?lectionn?s, chaque repr?sentation d'un habitat inclu dans le r?seau est incr?mentat?e de 1
habitats.counts=habitats.counts+table.site[the.site,] # Dans la liste des habitats s?lectionn?s, ajout des habitats appartenant au voisinage du site s?lectionn? pour la richesse en habitats
unselected.sites=unselected.sites[unselected.sites!=the.site] # Dans la liste des sites non s?lectionn?s, suppression du site ayant le score maximum

} #fin du while : tant qu'il reste des sites ? s?lectionner


## Ajout du dernier site restant

selected.sites=c(selected.sites,unselected.sites)

# selected.sites
#habitats.representation(table.site,selected.sites)

} # Fin de la fonction run.rich




########################
##   SOUS FONCTION    ##
########################

## convert2score : chaque ligne de pr?sence/absence d'habitat par site est
## convertie en un nombre binaires. Ces nombres binaires sont ensuite convertis
## en nombre d?cimaux pour la s?lection des sites les plus rares. L'habitat qui
## a le score le plus ?lev? est le plus rare.

convert2score.rich<-function(vec,habitats.richness.selected,complementarity){
if(complementarity==T){
tapply(as.numeric(vec),habitats.richness.selected,sum)}
else{
sum(as.numeric(vec))}
}


################################################################################