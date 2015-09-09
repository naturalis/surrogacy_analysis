################################################################################

##########################
##   FONCTION RUN.ALEA  ##
##########################

run.alea=function(table.site){


unselected.sites=rownames(table.site)
selected.sites=NULL #On initialise la liste de sites selectionnés à NULL
selected.sites=sample(unselected.sites,length(unselected.sites))
selected.sites

# habitats.representation(table.site,selected.sites)

} # Fin de la fonction run.alea

################################################################################