################################################################################



########################################
##   FONCTION species.richness.evol   ##
########################################


species.richness=function(stations.list,observations,method="cumul.sp")
  {

  observations=observations[observations$Code_Site %in% stations.list,] 

  if(method=="cumul.sp"){fonc.indice=nlevels(factor(observations$Code_poisson))}
  } # Fin de la fonction species.richness


species.richness.evol=function(stations.list,observations,method="cumul.sp")
  {

  richness.path=sapply(1:length(stations.list),function(x) species.richness(stations.list[1:x],observations,method))
  richness.path
  } # Fin de la fonction species.richness.evol



################################################################################