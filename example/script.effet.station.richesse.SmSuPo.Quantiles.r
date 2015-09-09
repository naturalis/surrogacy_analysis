spe <- read.csv("SmSuPoSp.csv")
recensement<-spe

# remove duplicates
recensement_unique=unique(recensement)

# create presence/absence table
tablo_pres_abs<-table(recensement_unique$Code_Site,recensement_unique$Code_poisson)

# calculate number of species per site
vec_poisson<-apply(tablo_pres_abs,1,sum)

# import randomization script
source('script.aleatoire.r')

# initialize matrix with 0 values
selection.aleatoire<-matrix(0,dim(tablo_pres_abs)[1],999)
for (i in 1:999){
  
  # populate the matrix, 1000 bootstraps
  selection.aleatoire[,i]<-run.alea(tablo_pres_abs)
}

#Selection des dites selon un scenario de richesse complementarit? bas? sur les poissons
source('script.richesse.complementarite.r')

selection.richesse.complementarite.poisson<-run.rich(tablo_pres_abs)

#Richesse de poisson incluse dans les sites selectionn?s (al?atoirement et selon diversit? de poisson)
source('script.diversite.selectionnee.r')
diversite.selection.aleatoire<-matrix(0,dim(tablo_pres_abs)[1],999)    #999
diversite.selection.aleatoire<-apply(as.matrix(selection.aleatoire),2,function(x) species.richness.evol(x,recensement_unique))
diversite.moyenne.selection.aleatoire<-apply(diversite.selection.aleatoire,1,mean)
diversite.selection.richesse.complementarite.poisson<-species.richness.evol(selection.richesse.complementarite.poisson,recensement_unique)

#Calcul des courbes al?atore moyenne, uper, et low
diversite.moyenne.selection.aleatoire
diversite.borne.sup.selection.aleatoire=apply(diversite.selection.aleatoire,1,function(x) sort(x)[ceiling(0.95*length(x))])
diversite.borne.inf.selection.aleatoire=apply(diversite.selection.aleatoire,1,function(x) sort(x)[ceiling(0.05*length(x))])

# total number of obsered species
nombre.poissons.total <- nlevels(factor(recensement_unique$Code_poisson))

# function to calculate percentages
percentage=function(total,truc){
  truc * 100 / total
}
  
# percentages, RANDOM, upper, low
diversite.moyenne.selection.aleatoire.pourcentage <- c(0,percentage(nombre.poissons.total,diversite.moyenne.selection.aleatoire))
diversite.borne.sup.selection.aleatoire.pourcentage <- c(0,percentage(nombre.poissons.total,diversite.borne.sup.selection.aleatoire))
diversite.borne.inf.selection.aleatoire.pourcentage <- c(0,percentage(nombre.poissons.total,diversite.borne.inf.selection.aleatoire))
diversite.selection.richesse.complementarite.poisson.pourcentage<- c(0,percentage(nombre.poissons.total,diversite.selection.richesse.complementarite.poisson))

gen <- read.csv("SmSuPoGen.csv")
# fam <- read.csv("SmSuPoFam.csv")
# sal<-read.csv("SmSuPoSal.csv")
# tho<-read.csv("SmSuPoTho.csv")
# ara<-read.csv("SmSuPoAra.csv")
# the<-read.csv("SmSuPoThe.csv")
# oxy<-read.csv("SmSuPoOxy.csv")

recensement1<-gen
# recensement.fam<-fam
# recensement.sal<-sal
# recensement.tho<-tho
# recensement.ara<-ara
# recensement.the<-the
# recensement.oxy<-oxy

recensement_unique1=unique(recensement1)
# recensement_unique_fam=unique(recensement.fam)
# recensement_unique_sal=unique(recensement.sal)
# recensement_unique_tho=unique(recensement.tho)
# recensement_unique_ara=unique(recensement.ara)
# recensement_unique_the=unique(recensement.the)
# recensement_unique_oxy=unique(recensement.oxy)

tablo_pres_abs1<-table(recensement_unique1$Code_Site,recensement_unique1$Code_poisson)
# tablo_pres_abs_fam<-table(recensement_unique_fam$Code_Site,recensement_unique_fam$Code_poisson)
# tablo_pres_abs_sal<-table(recensement_unique_sal$Code_Site,recensement_unique_sal$Code_poisson)
# tablo_pres_abs_tho<-table(recensement_unique_tho$Code_Site,recensement_unique_tho$Code_poisson)
# tablo_pres_abs_ara<-table(recensement_unique_ara$Code_Site,recensement_unique_ara$Code_poisson)
# tablo_pres_abs_the<-table(recensement_unique_the$Code_Site,recensement_unique_the$Code_poisson)
# tablo_pres_abs_oxy<-table(recensement_unique_oxy$Code_Site,recensement_unique_oxy$Code_poisson)

vec_poisson1<-apply(tablo_pres_abs1,1,sum)
# vec_poisson_fam<-apply(tablo_pres_abs_fam,1,sum)
# vec_poisson_sal<-apply(tablo_pres_abs_sal,1,sum)
# vec_poisson_tho<-apply(tablo_pres_abs_tho,1,sum)
# vec_poisson_ara<-apply(tablo_pres_abs_ara,1,sum)
# vec_poisson_the<-apply(tablo_pres_abs_the,1,sum)
# vec_poisson_oxy<-apply(tablo_pres_abs_oxy,1,sum)

selection.richesse.complementarite.gen<-run.rich(tablo_pres_abs1)
# selection.richesse.complementarite.fam<-run.rich(tablo_pres_abs_fam)
# selection.richesse.complementarite.sal<-run.rich(tablo_pres_abs_sal)
# selection.richesse.complementarite.tho<-run.rich(tablo_pres_abs_tho)
# selection.richesse.complementarite.ara<-run.rich(tablo_pres_abs_ara)
# selection.richesse.complementarite.the<-run.rich(tablo_pres_abs_the)
# selection.richesse.complementarite.oxy<-run.rich(tablo_pres_abs_oxy)

diversite.selection.richesse.complementarite.gen<-species.richness.evol(selection.richesse.complementarite.gen,recensement_unique)
# diversite.selection.richesse.complementarite.fam<-species.richness.evol(selection.richesse.complementarite.fam,recensement_unique)
# diversite.selection.richesse.complementarite.sal<-species.richness.evol(selection.richesse.complementarite.sal,recensement_unique)
# diversite.selection.richesse.complementarite.tho<-species.richness.evol(selection.richesse.complementarite.tho,recensement_unique)
# diversite.selection.richesse.complementarite.ara<-species.richness.evol(selection.richesse.complementarite.ara,recensement_unique)
# diversite.selection.richesse.complementarite.the<-species.richness.evol(selection.richesse.complementarite.the,recensement_unique)
# diversite.selection.richesse.complementarite.oxy<-species.richness.evol(selection.richesse.complementarite.oxy,recensement_unique)

# calculate SAI
source('script.sai.r')
SAI_richesse.complementarite.gen<-SAI(diversite.selection.richesse.complementarite.gen,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
# SAI_richesse.complementarite.fam<-SAI(diversite.selection.richesse.complementarite.fam,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
# SAI_richesse.complementarite.sal<-SAI(diversite.selection.richesse.complementarite.sal,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
# SAI_richesse.complementarite.tho<-SAI(diversite.selection.richesse.complementarite.tho,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
# SAI_richesse.complementarite.ara<-SAI(diversite.selection.richesse.complementarite.ara,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
# SAI_richesse.complementarite.the<-SAI(diversite.selection.richesse.complementarite.the,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
# SAI_richesse.complementarite.oxy<-SAI(diversite.selection.richesse.complementarite.oxy,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
# SAI_richesse.complementarite.sup<-SAI(diversite.borne.sup.selection.aleatoire,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)
# SAI_richesse.complementarite.inf<-SAI(diversite.borne.inf.selection.aleatoire,diversite.moyenne.selection.aleatoire,diversite.selection.richesse.complementarite.poisson)

diversite.selection.richesse.complementarite.gen.pourcentage<- c(0,percentage(nombre.poissons.total,diversite.selection.richesse.complementarite.gen))
# diversite.selection.richesse.complementarite.fam.pourcentage<- c(0,percentage(nombre.poissons.total,diversite.selection.richesse.complementarite.fam))
# diversite.selection.richesse.complementarite.sal.pourcentage<- c(0,percentage(nombre.poissons.total,diversite.selection.richesse.complementarite.sal))
# diversite.selection.richesse.complementarite.tho.pourcentage<- c(0,percentage(nombre.poissons.total,diversite.selection.richesse.complementarite.tho))
# diversite.selection.richesse.complementarite.oxy.pourcentage<- c(0,percentage(nombre.poissons.total,diversite.selection.richesse.complementarite.oxy))

# plot line graphs
plot(diversite.moyenne.selection.aleatoire.pourcentage,type="l", col="red",xlab="Total Samples",ylab="SAI")
lines(diversite.borne.sup.selection.aleatoire.pourcentage,lty=3, col="red")
lines(diversite.borne.inf.selection.aleatoire.pourcentage,lty=3, col="red")
lines(diversite.selection.richesse.complementarite.poisson.pourcentage,col="blue")
lines(diversite.selection.richesse.complementarite.gen.pourcentage,lty=1,col="black")
lines(diversite.selection.richesse.complementarite.sal.pourcentage,lty=2,col="black")
lines(diversite.selection.richesse.complementarite.tho.pourcentage,lty=3,col="black")
#lines(diversite.selection.richesse.complementarite.fam.pourcentage,lty=3, col="black")
lines(diversite.selection.richesse.complementarite.oxy.pourcentage,lty=4,col="black")


#lines(diversite.selection.richesse.complementarite.gna.pourcentage,lty=5,col="black")

#legend(9,80,c("Random curve","Random sup","Random inf","Optimal","Genus","Woody Vegetation","Family","Salticidae","Thomisidae","Gnaphosidae"),
#       cex=0.8,col=c("red","red","red","blue","green","green","black","black","black","black"),lty=c(1,3,3,1,1,3,1,2,3,4),bty="n")

legend(12,80,c("Genus","Salticidae","Thomisidae","Oxyopidae"),
       cex=0.8,col=c("black","black","black","black"),lty=c(1,2,3,4),bty="n")


