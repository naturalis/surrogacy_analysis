
         # Calcule un indice de diversit� sp�cifique pour une liste de stations donn��

species.rarete=function(stations.list,Indice_rarete_stations)
  {


Indice_rarete_selectionne<-sum(Indice_rarete_stations[names(Indice_rarete_stations)%in%stations.list])


  } # Fin de la fonction species.richness

#: Calcule l'ensemble des indices de richesse sp�cifique pas � pas
species.rarete.evol=function(stations.list,Indice_rarete_stations)
  {
  # Cr�ation d'un vecteur contenant le nombre d'esp�ces inclues dans le r�seau pour chaque site choisi dans l'ordre
  richness.path=sapply(1:length(stations.list),function(x) species.rarete(stations.list[1:x],Indice_rarete_stations))
  richness.path
  } # Fin de la fonction species.richness.evol