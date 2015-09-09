
         # Calcule un indice de diversité spécifique pour une liste de stations donnéé

species.rarete=function(stations.list,Indice_rarete_stations)
  {


Indice_rarete_selectionne<-sum(Indice_rarete_stations[names(Indice_rarete_stations)%in%stations.list])


  } # Fin de la fonction species.richness

#: Calcule l'ensemble des indices de richesse spécifique pas à pas
species.rarete.evol=function(stations.list,Indice_rarete_stations)
  {
  # Création d'un vecteur contenant le nombre d'espèces inclues dans le réseau pour chaque site choisi dans l'ordre
  richness.path=sapply(1:length(stations.list),function(x) species.rarete(stations.list[1:x],Indice_rarete_stations))
  richness.path
  } # Fin de la fonction species.richness.evol