


species.shannon=function(stations.list,observations,method="cumul.sp")
{
  # On ne garde que les stations et les attributs susceptible de nous intéresser
  observations=observations[observations$Code_Site %in% stations.list,]
  # diversite d'sp shannon (cumul.sp = nombre d'espèces)
       Npi=rep(0,length=dim(observations)[1])
  somme=sum(observations$Abondance)
 for (a in 1:dim(observations)[1]){
       Npi[a]=observations$Abondance[a]/somme
       }
       logNpi=log(Npi,2)
       produit=Npi*logNpi
       H=-sum(produit)
  #if(method=="cumul.sp"){fonc.indice=nlevels(factor(observations$Code_poisson))}
  
  } # Fin de la fonction species.richness

#: Calcule l'ensemble des indices de richesse spécifique pas à pas
species.shannon.evol=function(stations.list,observations,method="cumul.sp")
  {
  # Création d'un vecteur contenant le nombre d'espèces inclues dans le réseau pour chaque site choisi dans l'ordre
  shannon.path=sapply(1:length(stations.list),function(x) species.shannon(stations.list[1:x],observations,method))
  shannon.path
  } # Fin de la fonction species.richness.evol


