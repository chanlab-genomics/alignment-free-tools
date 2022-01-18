library(ape)
library(phytools)
library(Hmisc)

ref_tree <- read.tree("D:/2020_SS/BioInfo/reference_tree/michael_ref_tree.txt")
ref_tree <- makeNodeLabel(ref_tree,prefix="")
#list_ref <- list()
#for (i in ref_tree$tip){
#  list_ref[[length(list_ref)+1]] <- drop.tip(ref_tree,i)
#  plot(list_ref[[length(list_ref)]])
#}
#list_ref

path_pseudo = "D:/2020_SS/BioInfo/reference_tree/jk_trees"
list_files <- list.files(path_pseudo, full.names = TRUE)
list_files
trees <- list()
for (i in list_files){
  i <- read.tree(i)
  trees[[length(trees)+1]] <- i
}
trees
#count = 0
#list_comp <- list()
#for (i in list_ref){
#  for (k in trees){
#    for (j in i$tip){
#      for (l in k$tip){
#        test <- grepl(j,l)
#        if (test){B
#          inc(count) <- 1
#        }
#      }
#    }
#    if (count==length(i$tip)){
#      B <- prop.clades(i,k)
#      list_comp[[length(list_comp)+1]] <- B
#    }
#    count = 0
#  }
#}
#list_comp

plot(ref_tree)
nodelabels()
tiplabels()

A <- makeNodeLabel(ref_tree,prefix="")
A$node.label <- prop.clades(ref_tree,trees)
list = prop.clades(ref_tree,trees)
list
sum(list)/length(list) 
plot(A,show.node.label=TRUE)
write.tree(A, file="D:/2020_SS/BioInfo/reference_tree/jk_tree_100.txt")
