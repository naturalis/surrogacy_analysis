################################################################################



########################################
##   FONCTION species.richness.evol   ##
########################################

# Calcule un indice de diversité spécifique pour une liste de stations donnéé
species.richness=function(stations.list,observations,method="cumul.sp")
  {
  # On ne garde que les stations et les attributs susceptible de nous intéresser
  observations=observations[observations$Code_Site %in% stations.list,] 
  # Richesse spécifique (cumul.sp = nombre d'espèces)
  if(method=="cumul.sp"){fonc.indice=nlevels(factor(observations$Code_CT_SC_M))}
  } # Fin de la fonction species.richness

#: Calcule l'ensemble des indices de richesse spécifique pas à pas
species.richness.evol=function(stations.list,observations,method="cumul.sp")
  {
  # Création d'un vecteur contenant le nombre d'espèces inclues dans le réseau pour chaque site choisi dans l'ordre
  richness.path=sapply(1:length(stations.list),function(x) species.richness(stations.list[1:x],observations,method))
  richness.path
  } # Fin de la fonction species.richness.evol



################################################################################