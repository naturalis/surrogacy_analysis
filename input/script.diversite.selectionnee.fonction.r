################################################################################



########################################
##   FONCTION species.richness.evol   ##
########################################

# Calcule un indice de diversit� sp�cifique pour une liste de stations donn��
species.richness=function(stations.list,observations,method="cumul.sp")
  {
  # On ne garde que les stations et les attributs susceptible de nous int�resser
  observations=observations[observations$Code_Site %in% stations.list,] 
  # Richesse sp�cifique (cumul.sp = nombre d'esp�ces)
  if(method=="cumul.sp"){fonc.indice=nlevels(factor(observations$Code_CT_SC_M))}
  } # Fin de la fonction species.richness

#: Calcule l'ensemble des indices de richesse sp�cifique pas � pas
species.richness.evol=function(stations.list,observations,method="cumul.sp")
  {
  # Cr�ation d'un vecteur contenant le nombre d'esp�ces inclues dans le r�seau pour chaque site choisi dans l'ordre
  richness.path=sapply(1:length(stations.list),function(x) species.richness(stations.list[1:x],observations,method))
  richness.path
  } # Fin de la fonction species.richness.evol



################################################################################