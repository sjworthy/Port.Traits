setwd("~/Desktop/Invasive.Signal/Invasive.Signal/phylos/New.phylo/No.fern.gymno.both.seq.Ginkgo")

# Tree with only species that have rbcL and matK

# Building the phylogenies

# export .fasta files from Geneious

# run MAFFT alignment in Terminal
# place .fasta files into home/username folder
# Terminal code = mafft new.rbcl.fasta > new.rbcl.aln.fasta
# Terminal code = mafft new.matk.fasta > new.matk.aln.fasta

# Alignment Strategy Used = FFT-NS-2

library(ape)
library(picante)
library(geiger)
library(vegan)
library(seqinr)
library(phangorn)
library(phytools)

# Read in alignments and concatenate them
rbcl.aln=read.dna("new.rbcl.aln.fasta", format="fasta")
matk.aln=read.dna("new.matk.aln.fasta", format="fasta")

new.concat=cbind(rbcl.aln, matk.aln, fill.with.gaps=TRUE, check.names=FALSE)
write.dna(new.concat, "new.concat.fasta", format="fasta")

both.seq.concat=read.phyDat("new.concat.fasta", format = "fasta", type="DNA")


# Modeltest

modelTest(both.seq.concat) # GTR + G + I is the best

# Maximum likelihood tree

dist.both.seq=dist.logDet(both.seq.concat)
both.seq.nj.tree=NJ(dist.both.seq)
both.seq.ml.model=pml(both.seq.nj.tree, both.seq.concat)
both.seq.ml.tree=optim.pml(both.seq.ml.model, model="GTR", optGamma = TRUE, optInv = TRUE)
both.seq.optim.ml.tree=optim.pml(both.seq.ml.model, model="GTR", optGamma = TRUE, optInv = TRUE, optNni = TRUE)
both.seq.optim.ml.tree.rooted=root(both.seq.optim.ml.tree$tree, outgroup="Ginkgo biloba|rbcLa|KX283354", resolve.root=T)
plot(both.seq.optim.ml.tree.rooted, cex=0.3)

# Bootstrap trees
# Remove NA node labels in text editor

both.seq.boot=bootstrap.pml(both.seq.optim.ml.tree, bs=1000, optNni=TRUE)
both.seq.bootsrooted <- lapply(both.seq.boot, function(x) root(x, resolve.root=TRUE, outgroup="Ginkgo biloba|rbcLa|KX283354"))
class(both.seq.bootsrooted) <- "multiPhylo"
both.seqtree=plotBS(both.seq.optim.ml.tree$tree, both.seq.bootsrooted, p=75, type="phylo")
both.seqtree.rooted=root(both.seqtree, outgroup="Ginkgo biloba|rbcLa|KX283354", resolve.root=T)
both.seqtree.rooted$tip.label=all.tips
write.tree(both.seqtree.rooted, file="both.seq.tree.bs.txt")

plot(both.seqtree.rooted, cex=0.3)

# Tree with all species

# Building the phylogenies

# export .fasta files from Geneious

# run MAFFT alignment in Terminal
# place .fasta files into home/username folder
# Terminal code = mafft final.rbcl.fasta > new.rbcl.aln.fasta
# Terminal code = mafft final.matk.fasta > new.matk.aln.fasta


# Add sequences that only had (2) rbcL or (2)matk to the concatenated alignment in Terminal
# first add rbcL sequences
# second add matk sequences

# mafft --auto --addfragments only.rbcl.fasta new.concat.fasta>plus.rbcl.aln.fasta
# mafft --auto --addfragments only.matk.fasta plus.rbcl.aln.fasta>new.all.seq.aln.fasta


# Alignment Strategy Used = FFT-NS-2

# Read in alignments and concatenate them

all.seq.aln=read.phyDat("new.all.seq.aln.fasta", format = "fasta", type="DNA")

# Modeltest

modelTest(all.seq.aln) # GTR + G + I is the best

# Maximum likelihood tree

dist.all.seq=dist.logDet(all.seq.aln)
all.seq.nj.tree=NJ(dist.all.seq)
all.seq.ml.model=pml(all.seq.nj.tree, all.seq.aln)
all.seq.ml.tree=optim.pml(all.seq.ml.model, model="GTR", optGamma = TRUE, optInv = TRUE)
all.seq.optim.ml.tree=optim.pml(all.seq.ml.model, model="GTR", optGamma = TRUE, optInv = TRUE, optNni = TRUE)
all.seq.optim.ml.tree.rooted=root(all.seq.optim.ml.tree$tree, outgroup="Ginkgo biloba|rbcLa|KX283354", resolve.root=T)

plot(all.seq.optim.ml.tree.rooted, cex=0.3)

write.tree(all.seq.optim.ml.tree.rooted, file="test.tree.txt")

# Bootstrap trees
# Remove NA node labels in text editor

new.all.seq.boot=bootstrap.pml(all.seq.optim.ml.tree, bs=1000, optNni=TRUE)
new.all.seq.bootsrooted <- lapply(new.all.seq.boot, function(x) root(x, resolve.root=TRUE, outgroup="Ginkgo biloba|rbcLa|KX283354"))
class(new.all.seq.bootsrooted) <- "multiPhylo"

write.tree(new.all.seq.bootsrooted, file = "new.all.seq.bootsrooted.txt")

new.all.seqtree=plotBS(all.seq.optim.ml.tree$tree, new.all.seq.bootsrooted, p=75, type="phylo")

new.all.seqtree.rooted=root(new.all.seqtree, outgroup="Ginkgo biloba|rbcLa|KX283354", resolve.root=T)

new.all.seqtree.rooted$tip.label=all.tips
write.tree(new.all.seqtree.rooted, file="new.all.seq.tree.bs.txt")



# Phylogenetic Signal
setwd("~/Desktop/Invasive.Signal/Invasive.Signal/phylos")

library(phangorn)

# Nate Swenson's github: https://github.com/NGSwenson/SwensonSESYNCWorkshop2017/blob/master/day.3.am/signal.R

# DISPERSAL

# no trait data for species in object names so prune them from tree plus outgroup

names=c("Chaerophyllum tainturieri", "Cyperus compressus", "Cyperus planifolius", "Cyperus surinamensis",
        "Eragrostis secundiflora", "Galium tinctorium", "Glandularia tenera", "Jacquemontia tamnifolia",
        "Ludwigia bonariensis", "Oenothera laciniata", "Panicum scoparium", "Phalaris caroliniana",
        "Portulaca amilis", "Scutellaria racemosa", "Sisyrinchium rosulatum", "Solidago leavenworthii",
        "Stachys floridana", "Verbena brasiliensis", "Vernonia altissima", "Wahlenbergia marginata", "Ginkgo biloba")

# prune phylogeny
dispersal.tree=drop.tip(both.seqtree.rooted,names)
write.tree(dispersal.tree, file="dispersal.tree.txt")

# make matrix of traits
dispersal.matrix=matrix(data=NA, ncol=1, nrow=139)
row.names(dispersal.matrix)=dispersal.tree$tip.label

write.csv(dispersal.matrix, file="dispersal.matrix.csv")

dispersal.matrix=read.csv("dispersal.matrix.csv", row.names = 1, header=T)

dispersal.OBS.null <- c(NA)
for(i in 1:999){
  rand.tree = tipShuffle(dispersal.tree)
  obs=t(data.frame(dispersal.matrix))
  obs2<-phyDat(t(obs),type="USER",levels=attributes(factor(obs))$levels)
  dispersal.OBS.null[i]<-parsimony(rand.tree,obs2,method="sankoff")
}

disperse.obs.real = parsimony(dispersal.tree,obs2,method="sankoff")
# 68

p.value = (rank(c(disperse.obs.real,dispersal.OBS.null))[1])/1000
# p = 0.002


library(picante)

dispersal.trait.df=read.csv("dispersal.trait.df.csv", row.names = 1, header = T)
colnames(dispersal.trait.df)=dispersal.tree$tip.label

dispersal.ses.mpd=ses.mpd(dispersal.trait.df, cophenetic(dispersal.tree), 
                          null.model = "taxa.labels", runs = 999, iterations = 1000)

# Polychory is overdispersed p = 0.999, higher than expected mpd, overdispersion
# Anemochory is clustered p = 0.002, lower than expected mpd, clustered

write.csv(dispersal.ses.mpd, file="dispersal.ses.mpd.csv")

dispersal.ses.mntd=ses.mntd(dispersal.trait.df, cophenetic(dispersal.tree), 
                          null.model = "taxa.labels", runs = 999, iterations = 1000)

# Zoochory is clustered p = 0.013, lower than expected mntd, clustered

# Plot Dispersal

dispersalmode=as.factor(setNames(dispersal.matrix[,1], row.names(dispersal.matrix)))
dotTree(dispersal.tree, dispersalmode, colors = setNames(c("blue","red", "black","gray", "orange"),
                                                         c("Anemochory", "Autochory", "Hydrochory", "Polychory", "Zoochory")),
        ftype="i", fsize=0.7)
library(phytools)



# POLLINATION

# no trait data for species in object names so prune them from tree plus outgroup

names.pollination=c("Crotalaria incana", "Cyclospermum leptophyllum", "Cyperus planifolius", "Cyperus surinamensis", 
        "Cyperus virens", "Eragrostis secundiflora","Fimbristylis miliacea", "Gamochaeta purpurea", 
        "Glandularia tenera", "Jacquemontia tamnifolia","Juncus validus",  "Ludwigia bonariensis", 
        "Ludwigia decurrens", "Panicum scoparium", "Phalaris caroliniana", "Smilax smallii",
        "Strophostyles umbellata", "Ginkgo biloba")

# prune phylogeny
pollination.tree=drop.tip(both.seqtree.rooted,names.pollination)
write.tree(pollination.tree, file="pollination.tree.txt")

# remove Juncus validus because only species with hydrophily

# make matrix of traits
pollination.matrix=matrix(data=NA, ncol=1, nrow=142)
row.names(pollination.matrix)=pollination.tree$tip.label

write.csv(pollination.matrix, file="pollination.matrix.csv")

pollination.matrix=read.csv("pollination.matrix.csv", row.names = 1, header=T)

pollination.OBS.null <- c(NA)
for(i in 1:999){
  rand.tree = tipShuffle(pollination.tree)
  obs=t(data.frame(pollination.matrix))
  obs2<-phyDat(t(obs),type="USER",levels=attributes(factor(obs))$levels)
  pollination.OBS.null[i]<-parsimony(rand.tree,obs2,method="sankoff")
}

pollination.obs.real = parsimony(pollination.tree,obs2,method="sankoff")
# 50

pollination.p.value = (rank(c(pollination.obs.real,pollination.OBS.null))[1])/1000
# p = 0.001

pollination.trait.df=read.csv("pollination.trait.df.csv", row.names = 1, header = T)
colnames(pollination.trait.df)=pollination.tree$tip.label

pollination.ses.mpd=ses.mpd(pollination.trait.df, cophenetic(pollination.tree), 
                          null.model = "taxa.labels", runs = 999, iterations = 1000)


# Zoophily is clumped p = 0.001, lower than expected mpd
# Selfing is overdispersed p = 0.981, higher than expected mpd

write.csv(pollination.ses.mpd, file="pollination.ses.mpd.csv")

pollination.ses.mntd=ses.mntd(pollination.trait.df, cophenetic(pollination.tree), 
                            null.model = "taxa.labels", runs = 999, iterations = 1000)
# Anemophily is clumped p = 0.005, lower than expected mntd

write.csv(pollination.ses.mntd, file="pollination.ses.mntd.csv")


# DURATION

# prune root

names.duration="Ginkgo biloba"

# prune phylogeny

duration.tree=drop.tip(both.seqtree.rooted, names.duration)
write.tree(duration.tree, file="duration.tree.txt")

# make matrix of traits
duration.matrix=matrix(data=NA, ncol=1, nrow=159)
row.names(duration.matrix)=duration.tree$tip.label

write.csv(duration.matrix, file="duration.matrix.csv")

duration.matrix=read.csv("duration.matrix.csv", row.names = 1, header=T)

duration.OBS.null <- c(NA)
for(i in 1:999){
  rand.tree = tipShuffle(duration.tree)
  obs=t(data.frame(duration.matrix))
  obs2<-phyDat(t(obs),type="USER",levels=attributes(factor(obs))$levels)
  duration.OBS.null[i]<-parsimony(rand.tree,obs2,method="sankoff")
}

duration.obs.real = parsimony(duration.tree,obs2,method="sankoff")
# 58

duration.p.value = (rank(c(duration.obs.real,duration.OBS.null))[1])/1000
# p = 0.001

duration.trait.df=read.csv("duration.trait.df.csv", row.names = 1, header = T)
colnames(duration.trait.df)=duration.tree$tip.label

duration.ses.mpd=ses.mpd(duration.trait.df, cophenetic(duration.tree), 
                            null.model = "taxa.labels", runs = 999, iterations = 1000)
# Nothing Significant

duration.ses.mntd=ses.mntd(duration.trait.df, cophenetic(duration.tree), 
                              null.model = "taxa.labels", runs = 999, iterations = 1000)

# Annual is clumped, p = 0.011, mntd lower than expected

# GROWTH HABIT

# prune root and 1 subshrub species

names.ghabit=c("Ginkgo biloba", "Rubus argutus")

# prune phylogeny

ghabit.tree=drop.tip(both.seqtree.rooted, names.ghabit)
write.tree(ghabit.tree, file="ghabit.tree.txt")

# make matrix of traits
ghabit.matrix=matrix(data=NA, ncol=1, nrow=159)
row.names(ghabit.matrix)=ghabit.tree$tip.label

write.csv(ghabit.matrix, file="ghabit.matrix.csv")

ghabit.matrix=read.csv("ghabit.matrix.csv", row.names = 1, header=T)

ghabit.OBS.null <- c(NA)
for(i in 1:999){
  rand.tree = tipShuffle(ghabit.tree)
  obs=t(data.frame(ghabit.matrix))
  obs2<-phyDat(t(obs),type="USER",levels=attributes(factor(obs))$levels)
  ghabit.OBS.null[i]<-parsimony(rand.tree,obs2,method="sankoff")
}

ghabit.obs.real = parsimony(ghabit.tree,obs2,method="sankoff")
# 30

ghabit.p.value = (rank(c(ghabit.obs.real,ghabit.OBS.null))[1])/1000
# p = 0.001

ghabit.trait.df=read.csv("ghabit.trait.df.csv", row.names = 1, header = T)
colnames(ghabit.trait.df)=ghabit.tree$tip.label

ghabit.ses.mpd=ses.mpd(ghabit.trait.df, cophenetic(ghabit.tree), 
                         null.model = "taxa.labels", runs = 999, iterations = 1000)
# Forb, clumped p = 0.001, lower than expected mpd
# Graminoid, clumped p = 0.001, lower than expected mpd
# Multiple, clumped p = 0.003, lower than expected mpd
# Tree, clumped p = 0.014, lower than expected mpd
# Vine, clumped p = 0.019, lower than expected mpd

write.csv(ghabit.ses.mpd, file="ghabit.ses.mpd.csv")

ghabit.ses.mntd=ses.mntd(ghabit.trait.df, cophenetic(ghabit.tree), 
                           null.model = "taxa.labels", runs = 999, iterations = 1000)

# Forb, clumped p = 0.008, lower than expected mntd
# Graminoid, clumped p = 0.001, lower than expected mntd

write.csv(ghabit.ses.mntd, file="ghabit.ses.mntd.csv")

# STATUS
# binary trait so use phylo.d
# Native = 0, Introduced = 1
# prune root

names.status="Ginkgo biloba"

# prune phylogeny

status.tree=drop.tip(all.seqtree.rooted, names.status)

status.data=matrix(data=NA, nrow = 170, ncol=2)

colnames(status.data)=c("binomial", "status")

status.data=as.data.frame(status.data)

status.data$binomial=status.tree$tip.label

write.csv(status.data, file="status.matrix.csv")

status.data=read.csv("status.matrix.csv", row.names = 1, header=T)

# remove node labels from status.tree

status.tree.2=status.tree
status.tree.2$node.label=NULL

status.compar=comparative.data(status.tree.2, status.data, binomial)

statusPhyloD=phylo.d(status.compar, binvar=status, permut = 1000)

print(statusPhyloD)
plot(statusPhyloD)

#Estimated D = 0.7977201
#Probability resulting from no (random) phylogenetic structure = 0.018
#Probability resulting from Brownian phylogenetic strucutre = 0
#Observed = 61.818
#MeanRandom = 69.60332
#MeanBrownian = 31.11549
# simulation tests indicated that the phylogenetic pattern differed significantly from both the 
#Brownian  expectation and from 1 (randomonness). Weber et al. 2013 (weak phylogenetic signal)





