##############################
##     FONCTION : SAI       ##
##############################

# Calcul de l'aire sous la courbe (Area Under Curve = AUC)

AUC=function(path){
AUC=sum(sapply(2:length(path),function(x) path[x]-(path[x]-path[x-1])/2))
AUC
}



SAI=function(surrogate,random,optimum){
SAI=(AUC(surrogate)-AUC(random))/(AUC(optimum)-AUC(random))
SAI
}

#SAI=function(surrogate.richcomp,random,species.richcomp){
#SAI=(AUC(surrogate.richcomp)-AUC(random))/(AUC(species.richcomp)-AUC(random))
#SAI
#}




# species.rarcomp = F (optimum, SAI=1)
# species.richcomp
# species.worst


 
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
