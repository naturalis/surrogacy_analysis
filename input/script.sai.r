##############################
##     FONCTION : SAI       ##
##############################

# Calcul de l'aire sous la courbe (Area Under Curve = AUC)

AUC=function(path){
AUC=sum(sapply(2:length(path),function(x) path[x]-(path[x]-path[x-1])/2))
AUC
}

# Calcul de l'indice d'accumulation pour les espèces (Species Accumulation Index = SAI)

SAI=function(surrogate,random,optimum){
SAI=(AUC(surrogate)-AUC(random))/(AUC(optimum)-AUC(random))
SAI
}

#SAI=function(surrogate.richcomp,random,species.richcomp){
#SAI=(AUC(surrogate.richcomp)-AUC(random))/(AUC(species.richcomp)-AUC(random))
#SAI
#}

# surrogate.rarcomp = H1 (rareté/complémentarité de l'habitat)
# surrogate.richcomp = H2 (richesse/complémentarité de l'habitat)
# surrogate.worst = F# (pire scénario pour l'habitat, SAI=0)
# species.rarcomp = F (optimum, SAI=1)
# species.richcomp
# species.worst
# random = A (aléatoire)

 
SAI.calc=function(surrogate.rarcomp,surrogate.richcomp,surrogate.worst,species.rarcomp,species.richcomp,species.worst,random)
{

SAIs=c(
SAI(surrogate.rarcomp,random,species.rarcomp),
SAI(surrogate.richcomp,random,species.rarcomp),
SAI(surrogate.worst,random,species.rarcomp),
SAI(species.worst,random,species.rarcomp)
)

names(SAIs)=c("surrogate.rarcomp","surrogate.richcomp","surrogate.worst","species.worst")

sort(SAIs)
}
