############################################################################
## This file contains the R commands used to summarize the FP index 
## values for all of the data sets.
## Commands to produce the figures in the manuscript are also 
## included.
############################################################################


##### Plot settings #######

library(ggplot2)
mynamestheme <- theme(
	axis.title = element_text(size=(18),face="bold"),
	axis.text = element_text(size=(14),face="bold")
)



#########################
# 1. Hu-Ch-Go (primate)

primate = read.table("human_chimp_gorilla_fp_combined_table.txt",header=T)

ncol = dim(primate)[2]

primate_ranks = matrix(rank(-1*primate[,2],ties.method="min"),ncol=1)
for (i in 3: ncol) primate_ranks = cbind(primate_ranks,rank(-1*primate[,i],ties.method="min"))

dim(primate_ranks)[2]
#54

primate_iqr=apply(primate_ranks[,1:52],1,IQR)
summary(primate_iqr)

# transpose input matrix, omit the last 2 columns
temp = primate_ranks[,1:52]
primate_ranks_t = data.frame(t(temp))

primate_plot <- ggplot(stack(primate_ranks_t),aes(x=ind,y=values)) + 
    geom_boxplot(aes(fill=ind),notch=T,lwd=1.2) + 
    geom_point(data=data.frame(x=factor(c("X1","X2","X3","X4")),y=primate_ranks[,53]),aes(x=x,y=y,pch="Mean over gene trees"),color='darkred',size=(6),show.legend=T) + 
    geom_point(data=data.frame(x=factor(c("X1","X2","X3","X4")),y=primate_ranks[,54]),aes(x=x,y=y,pch="Rank from species tree"),color='darkred',size=(6),show.legend=T) +
    theme_minimal() + 
    labs(x="Species",y="Rank")  +
    scale_x_discrete(labels=c("X1"="H. sapiens","X2" = "P. troglodytes","X3"="G. gorilla","X4" = "P. pygmaeus")) 
            
print(primate_plot + mynamestheme + scale_shape_manual(values=c('Mean over gene trees' = "circle", 'Rank from species tree' = "triangle"), name = ' ') + guides(fill = "none") + theme(legend.position="top", legend.text = element_text(face = "bold")) + theme(legend.background = element_rect(fill="NA",size = 0.5, linetype="solid",colour="darkblue")))

primate_cor = cor(primate_ranks[,1:52],method="kendall")
pairwise = lower.tri(primate_cor,diag=F)
primate_cormat = primate_cor[pairwise]

primate_gt_st_cor = matrix(0,ncol=1,nrow=52)
for (i in 1:52) primate_gt_st_cor[i,1] = cor(primate_ranks[,i],primate_ranks[,54],method="kendall")

cor(primate_ranks[,53],primate_ranks[,54],method="kendall")



############################
2. Rattlesnake

rattle = read.table("rattlesnake_fp_combined_table.txt",header=T)
ncol = dim(rattle)[2]

rattle_ranks = matrix(rank(-1*rattle[,2],ties.method="min"),ncol=1)
for (i in 3: ncol) rattle_ranks = cbind(rattle_ranks,rank(-1*rattle[,i],ties.method="min"))

dim(rattle_ranks)[2]
#18

rattle_iqr=apply(rattle_ranks[,1:16],1,IQR)
summary(rattle_iqr)

# transpose input matrix, omit the last 2 columns
temp_rattle = rattle_ranks[,1:16]
rattle_ranks_t = data.frame(t(temp_rattle))

rattle_plot <- ggplot(stack(rattle_ranks_t),aes(x=ind,y=values)) + 
    geom_boxplot(aes(fill=ind),notch=T,lwd=1.2) + 
    geom_point(data=data.frame(x=factor(c("X1","X2","X3","X4","X5","X6","X7")),y=rattle_ranks[,17]),aes(x=x,y=y,pch="Rank from mean over gene trees"),color='darkred',size=(6),show.legend=T) + 
    geom_point(data=data.frame(x=factor(c("X1","X2","X3","X4","X5","X6","X7")),y=rattle_ranks[,18]),aes(x=x,y=y,pch="Rank from species tree"),color='darkred',size=(6),show.legend=T) +
    theme_minimal() + 
    labs(x="Species",y="Rank")  +
    scale_x_discrete(labels=c("X1"="S. c. edwardsii","X2" = "S. c. catenatus","X3"="S. c. tergeminus","X4" = "S .m. streckeri ","X5" = "S. m. miliarius","X6" = "S. m. barbouri ","X7" = "A. contortrix")) +
    scale_y_continuous(limits=c(0,8))
            
print(rattle_plot + mynamestheme + scale_shape_manual(values=c('Rank from mean over gene trees' = "circle", 'Rank from species tree' = "triangle"), name = ' ') + guides(fill = "none") + theme(legend.position="top", legend.text = element_text(face = "bold")) + theme(legend.background = element_rect(fill="NA",size = 0.5, linetype="solid",colour="darkblue")))


rattle_cor = cor(rattle_ranks[,1:16],method="kendall")
pairwise = lower.tri(rattle_cor,diag=F)
rattle_cormat = rattle_cor[pairwise]
mean(rattle_cormat)

rattle_gt_st_cor = matrix(0,ncol=1,nrow=16)
for (i in 1:16) rattle_gt_st_cor[i,1] = cor(rattle_ranks[,i],rattle_ranks[,18],method="kendall")
mean(rattle_gt_st_cor)

cor(rattle_ranks[,17],rattle_ranks[,18],method="kendall")


############################
# 3. Yeast

yeast = read.table("yeast_fp_combined_table.txt",header=T)
ncol = dim(yeast)[2]
#109

yeast_ranks = matrix(rank(-1*yeast[,2],ties.method="min"),ncol=1)
for (i in 3: ncol) yeast_ranks = cbind(yeast_ranks,rank(-1*yeast[,i],ties.method="min"))

dim(yeast_ranks)[2]
#108

yeast_iqr=apply(yeast_ranks[,1:106],1,IQR)
summary(yeast_iqr)


temp_yeast = yeast_ranks[,1:106]
yeast_ranks_t = data.frame(t(temp_yeast))

yeast_plot <- ggplot(stack(yeast_ranks_t),aes(x=ind,y=values)) + 
    geom_boxplot(aes(fill=ind),notch=T,lwd=1.2) + 
    geom_point(data=data.frame(x=factor(c("X1","X2","X3","X4","X5","X6","X7","X8")),y=yeast_ranks[,107]),aes(x=x,y=y,pch="Rank from mean over gene trees"),color='darkred',size=(6),show.legend=T) + 
    geom_point(data=data.frame(x=factor(c("X1","X2","X3","X4","X5","X6","X7","X8")),y=yeast_ranks[,108]),aes(x=x,y=y,pch="Rank from species tree"),color='darkred',size=(6),show.legend=T) +
    theme_minimal() + 
    labs(x="Species",y="Rank")  +
    scale_x_discrete(labels=c("X1"="S. kluyveri", "X2"="S. cerevisiae", "X3"="S. paradoxus", "X4"="S. mikatae", "X5"="S. kudriavzevii", "X6"="S. bayanus", "X7"="S. castellii", "X8"="C. albicans")) +
    scale_y_continuous(limits=c(0,9))
                
print(yeast_plot + mynamestheme + scale_shape_manual(values=c('Rank from mean over gene trees' = "circle", 'Rank from species tree' = "triangle"), name = ' ') + guides(fill = "none") + theme(legend.position="top", legend.text = element_text(face = "bold")) + theme(legend.background = element_rect(fill="NA",size = 0.5, linetype="solid",colour="darkblue")))

yeast_cor = cor(yeast_ranks[,1:106],method="kendall")
pairwise = lower.tri(yeast_cor,diag=F)
yeast_cormat = yeast_cor[pairwise]
mean(yeast_cormat)

yeast_gt_st_cor = matrix(0,ncol=1,nrow=106)
for (i in 1:106) yeast_gt_st_cor[i,1] = cor(yeast_ranks[,i],yeast_ranks[,108],method="kendall")
mean(yeast_gt_st_cor)

cor(yeast_ranks[,107],yeast_ranks[,108],method="kendall")




##################################
# 4. Dophin

dolphin = read.table("dolphin_fp_combined_table.txt",header=T)
ncol = dim(dolphin)[2]
#25

dolphin_ranks = matrix(rank(-1*dolphin[,2],ties.method="min"),ncol=1)
for (i in 3: ncol) dolphin_ranks = cbind(dolphin_ranks,rank(-1*dolphin[,i],ties.method="min"))

dim(dolphin_ranks)[2]
# 24

dolphin_iqr=apply(dolphin_ranks[,1:22],1,IQR)
summary(dolphin_iqr)

temp_dolphin = dolphin_ranks[,1:22]
dolphin_ranks_t = data.frame(t(temp_dolphin))

dolphin_plot <- ggplot(stack(dolphin_ranks_t),aes(x=ind,y=values)) + 
    geom_boxplot(aes(fill=ind),notch=T,lwd=1.2) + 
    geom_point(data=data.frame(x=factor(c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19","X20","X21","X22","X23","X24","X25","X26","X27","X28")),y=dolphin_ranks[,23]),aes(x=x,y=y,pch="Rank from mean over gene trees"),color='darkred',size=(6),show.legend=T) +    geom_point(data=data.frame(x=factor(c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19","X20","X21","X22","X23","X24","X25","X26","X27","X28")),y=dolphin_ranks[,24]),aes(x=x,y=y,pch="Rank from species tree"),color='darkred',size=(6),show.legend=T) +
    theme_minimal() + 
    labs(x="Species",y="Rank")  +
    scale_x_discrete(labels=c("X1"="S. bredanensis", "X2"="G. griseus", "X3"="G. macrorhynchus", "X4"="P. electra", "X5"="F. attenuata", "X6"="P. crassidens","X7"="D. capensis", "X8"="D. delphis", "X9"="S. longirostris", "X10"="T. aduncus", "X11"="T. truncatus", "X12"="S. coeruleoalba", "X13"="S. chinensis", "X14"="L. hosei", "X15"="S. attenuata", "X16"="S. fluviatilis", "X17"="C. commersonii", "X18"="L. borealis", "X19"="L. obliquidens", "X20"="L. obscurus", "X21"="O. orca", "X22"="L. albirostris", "X23"="L. acutus", "X24"="N. phocaenoides","X25"="P. phocoena", "X26"="D. leucas", "X27"="M. monoceros", "X28"="P. macrocephalus")) +
    scale_y_continuous(limits=c(0,29))
                
print(dolphin_plot + mynamestheme + scale_shape_manual(values=c('Rank from mean over gene trees' = "circle", 'Rank from species tree' = "triangle"), name = ' ') + guides(fill = "none") + theme(legend.position="top", legend.text = element_text(face = "bold",size=(15)),axis.text.x=element_text(size=(15),angle=45,hjust=1,vjust=1)) + theme(legend.background = element_rect(fill="NA",size = 0.5, linetype="solid",colour="darkblue")))


dolphin_cor = cor(dolphin_ranks[,1:22],method="kendall")
pairwise = lower.tri(dolphin_cor,diag=F)
dolphin_cormat = dolphin_cor[pairwise]
mean(dolphin_cormat)

dolphin_gt_st_cor = matrix(0,ncol=1,nrow=22)
for (i in 1:22) dolphin_gt_st_cor[i,1] = cor(dolphin_ranks[,i],dolphin_ranks[,24],method="kendall")
mean(dolphin_gt_st_cor)

cor(dolphin_ranks[,23],dolphin_ranks[,24],method="kendall")






##############################
# 5. Mammals

mammals = read.table("mammals_fp_combined_table.txt",header=T)
ncol = dim(mammals)[2]
# 450

mammals_ranks = matrix(rank(-1*mammals[,2],ties.method="min"),ncol=1)
for (i in 3: ncol) mammals_ranks = cbind(mammals_ranks,rank(-1*mammals[,i],ties.method="min"))

dim(mammals_ranks)[2]
# 449

mammals_iqr=apply(mammals_ranks[,1:447],1,IQR)
summary(mammals_iqr)

temp_mammals = mammals_ranks[,1:447]
mammals_ranks_t = data.frame(t(temp_mammals))

mammals_plot <- ggplot(stack(mammals_ranks_t),aes(x=ind,y=values)) + 
    geom_boxplot(aes(fill=ind),notch=T,lwd=1.2) +     geom_point(data=data.frame(x=factor(c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19","X20","X21","X22","X23","X24","X25","X26","X27","X28","X29","X30","X31","X32","X33","X34","X35","X36","X37")),y=mammals_ranks[,448]),aes(x=x,y=y,pch="Rank from mean over gene trees"),color='darkred',size=(6),show.legend=T) +    geom_point(data=data.frame(x=factor(c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19","X20","X21","X22","X23","X24","X25","X26","X27","X28","X29","X30","X31","X32","X33","X34","X35","X36","X37")),y=mammals_ranks[,449]),aes(x=x,y=y,pch="Rank from species tree"),color='darkred',size=(6),show.legend=T) +
    theme_minimal() + 
    labs(x="Species",y="Rank")  +
    scale_x_discrete(labels=c("X1"="G. gallus", "X2"="M. eugenii", "X3"="M. domestica", "X4"="M. mulatta", "X5"="G. gorilla", "X6"="H. sapiens", "X7"="P. troglodytes", "X8"="P. pygmaeus", "X9"="C. jacchus", "X10"="T. syrichta", "X11"="O. garnettii", "X12"="M. murinus", "X13"="T. belangeri", "X14"="S. tridecemlineatus", "X15"="D. ordii", "X16"="M. musculus", "X17"="Rattus sp.", "X18"="C. porcellus", "X19"="O. cuniculus", "X20"="O. princeps", "X21"="S. araneus", "X22"="E. europaeus", "X23"="T. truncatus", "X24"="B. taurus", "X25"="S. scrofa", "X26"="V. pacos", "X27"="F. catus", "X28"="C. familiaris", "X29"="E. caballus", "X30"="M. lucifugus", "X31"="P. vampyrus", "X32"="L. africana", "X33"="P. capensis", "X34"="E. telfairi",  "X35"="D. novemcinctus", "X36"="C. hoffmanni", "X37"="O. anatinus")) 
    #scale_y_continuous(limits=c(0,9))
                
print(mammals_plot + mynamestheme + scale_shape_manual(values=c('Rank from mean over gene trees' = "circle", 'Rank from species tree' = "triangle"), name = ' ') + guides(fill = "none") + theme(legend.position="top", legend.text = element_text(face = "bold"),axis.text.x=element_text(size=(12),angle=45,hjust=1,vjust=1)) + theme(legend.background = element_rect(fill="NA",size = 0.5, linetype="solid",colour="darkblue")))


mammals_cor = cor(mammals_ranks[,1:447],method="kendall")
pairwise = lower.tri(mammals_cor,diag=F)
mammals_cormat = mammals_cor[pairwise]
mean(mammals_cormat)

mammals_gt_st_cor = matrix(0,ncol=1,nrow=447)
for (i in 1:447) mammals_gt_st_cor[i,1] = cor(mammals_ranks[,i],mammals_ranks[,449],method="kendall")
mean(mammals_gt_st_cor)

cor(mammals_ranks[,448],mammals_ranks[,449],method="kendall")




############################
# 6. Snake


# Note: I had to read this one in differently -- so the dimensions are different then the other cases because I already removed the first column
mysnakes = matrix(scan("snake_fp_combined_table_try2.csv",skip=1,sep=','),nrow=33,byrow=T)
snakes = data.frame(mysnakes)
ncol = dim(snakes)[2]
# 335

snakes_ranks = matrix(rank(-1*snakes[,1],ties.method="min"),ncol=1)
for (i in 2: ncol) snakes_ranks = cbind(snakes_ranks,rank(-1*snakes[,i],ties.method="min"))

dim(snakes_ranks)[2]
# 335

snakes_iqr=apply(snakes_ranks[,1:333],1,IQR)
summary(snakes_iqr)

temp_snakes = snakes_ranks[,1:333]
snakes_ranks_t = data.frame(t(temp_snakes))


snakes_plot <- ggplot(stack(snakes_ranks_t),aes(x=ind,y=values)) + 
    geom_boxplot(aes(fill=ind),notch=T,lwd=1.2) +     geom_point(data=data.frame(x=factor(c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19","X20","X21","X22","X23","X24","X25","X26","X27","X28","X29","X30","X31","X32","X33")),y=snakes_ranks[,334]),aes(x=x,y=y,pch="Rank from mean over gene trees"),color='darkred',size=(6),show.legend=T) +    geom_point(data=data.frame(x=factor(c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19","X20","X21","X22","X23","X24","X25","X26","X27","X28","X29","X30","X31","X32","X33")),y=snakes_ranks[,335]),aes(x=x,y=y,pch="Rank from species tree"),color='darkred',size=(6),show.legend=T) +
    theme_minimal() + 
    labs(x="Species",y="Rank")  +
    scale_x_discrete(labels=c("X1"="A. carolinensis", "X2"="P. molurus", "X3"="A. granulatus", "X4"="X. javanicus", "X5"="P. carinatus", "X6"="E. carinatus", "X7"="A. feae", "X8"="A. contortrix", "X9"="C. adamanteus", "X10"="C. horridus", "X11"="G. prevostiana", "X12"="L. curtus", "X13"="C. melanurus", "X14"="L. inornatus", "X15"="A. capensis", "X16"="A. bibroni", "X17"="P. cana", "X18"="P. frontali", "X19"="P. notostictus", "X20"="H. modestus", "X21"="G. klingi", "X22"="S. bistrigatus", "X23"="G. smithii", "X24"="M. taeniatum", "X25"="C. pavimentata", "X26"="P. macrops", "X27"="D. punctatus", "X28"="H. nasicus", "X29"="N. rhombifer", "X30"="R. rigida",  "X31"="S. storerioides", "X32"="S. dekayi", "X33"="S. occipitomaculata")) 
    #scale_y_continuous(limits=c(0,9))
                
print(snakes_plot + mynamestheme + scale_shape_manual(values=c('Rank from mean over gene trees' = "circle", 'Rank from species tree' = "triangle"), name = ' ') + guides(fill = "none") + theme(legend.position="top", legend.text = element_text(face = "bold"),axis.text.x=element_text(size=(12),angle=45,hjust=1,vjust=1)) + theme(legend.background = element_rect(fill="NA",size = 0.5, linetype="solid",colour="darkblue")))


snakes_cor = cor(snakes_ranks[,1:333],method="kendall")
pairwise = lower.tri(snakes_cor,diag=F)
snakes_cormat = snakes_cor[pairwise]
mean(snakes_cormat)

snakes_gt_st_cor = matrix(0,ncol=1,nrow=333)
for (i in 1:333) snakes_gt_st_cor[i,1] = cor(snakes_ranks[,i],snakes_ranks[,335],method="kendall")
mean(snakes_gt_st_cor)

cor(snakes_ranks[,334],snakes_ranks[,335],method="kendall")



###########################
# 7. Fungi

# Note: same issue as with snakes above
myfungi = matrix(scan("fungi_fp_combined_table_try2.csv",skip=1,sep=','),nrow=29,byrow=T)
fungi = data.frame(myfungi)
ncol = dim(fungi)[2]
# 685

fungi_ranks = matrix(rank(-1*fungi[,1],ties.method="min",na.last="keep"),ncol=1)
for (i in 2: ncol) fungi_ranks = cbind(fungi_ranks,rank(-1*fungi[,i],ties.method="min",na.last="keep"))

dim(fungi_ranks)[2]
# 685

fungi_iqr=apply(fungi_ranks[,1:683],1,IQR)
summary(fungi_iqr)

temp_fungi = fungi_ranks[,1:683]
fungi_ranks_t = data.frame(t(temp_fungi))

fungi_plot <- ggplot(stack(fungi_ranks_t),aes(x=ind,y=values)) + 
    geom_boxplot(aes(fill=ind),notch=T,lwd=1.2) +     geom_point(data=data.frame(x=factor(c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19","X20","X21","X22","X23","X24","X25","X26","X27","X28","X29")),y=fungi_ranks[,684]),aes(x=x,y=y,pch="Rank from mean over gene trees"),color='darkred',size=(6),show.legend=T) +    geom_point(data=data.frame(x=factor(c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19","X20","X21","X22","X23","X24","X25","X26","X27","X28","X29")),y=fungi_ranks[,685]),aes(x=x,y=y,pch="Rank from species tree"),color='darkred',size=(6),show.legend=T) +
    theme_minimal() + 
    labs(x="Species",y="Rank")  +
    scale_x_discrete(labels=c("X1"="C. jadinii", "X2"="W. anomalus", "X3"="K. marxianus", "X4"="S. cerevisiae", "X5"="Hanseniaspora sp.", "X6"="H. valbyensis", "X7"="H. singularis",
"X8"="H. jakobsenii", "X9"="H. hatyaiensis", "X10"="H. lachancei","X11"="H. pseudoguilliermondii", "X12"="H. opuntiae", "X13"="H. guilliermondii UTAD222", "X14"="H. guilliermondii CBS465", "X15"="H. nectarophila", "X16"="H. thailandica", "X17"="H. meyeri","X18"="H. clermontiae", "X19"="H. uvarum CBS314", "X20"="H. uvarum AWRI3580", "X21"="H. uvarum DSM2768", "X22" = "H. uvarum 34-9", "X23"="H. gamundiae", "X24"="H. occidentalis var. citrica", "X25"="H. occidentalis var. occidentalis", "X26"="H. osmophila CBS313", "X27"="H. osmophila AWRI3579", "X28"="H. vineae T0219AF", "X29"="H. vineae NRRLY-1626")) 
    #scale_y_continuous(limits=c(0,9))
                
print(fungi_plot + mynamestheme + scale_shape_manual(values=c('Rank from mean over gene trees' = "circle", 'Rank from species tree' = "triangle"), name = ' ') + guides(fill = "none") + theme(legend.position="top", legend.text = element_text(face = "bold"),axis.text.x=element_text(size=(12),angle=45,hjust=1,vjust=1)) + theme(legend.background = element_rect(fill="NA",size = 0.5, linetype="solid",colour="darkblue")))

fungi_cor = cor(fungi_ranks[,1:683],method="kendall")
pairwise = lower.tri(fungi_cor,diag=F)
fungi_cormat = fungi_cor[pairwise]
mean(fungi_cormat)

fungi_gt_st_cor = matrix(0,ncol=1,nrow=683)
for (i in 1:683) fungi_gt_st_cor[i,1] = cor(fungi_ranks[,i],fungi_ranks[,685],method="kendall")
mean(fungi_gt_st_cor)

cor(fungi_ranks[,684],fungi_ranks[,685],method="kendall")


###########################
# 8. Animal

animals = read.table("animal_fp_combined_table.txt",header=T)
ncol = dim(animals)[2]
# 764

animals_ranks = matrix(rank(-1*animals[,2],ties.method="min",na.last="keep"),ncol=1)
for (i in 3: ncol) animals_ranks = cbind(animals_ranks,rank(-1*animals[,i],ties.method="min",na.last="keep"))

dim(animals_ranks)[2]
# 763

animals_iqr=apply(animals_ranks[,1:761],1,IQR)
summary(animals_iqr)

temp_animals = animals_ranks[,1:761]
animals_ranks_t = data.frame(t(temp_animals))

animals_plot <- ggplot(stack(animals_ranks_t),aes(x=ind,y=values)) + 
    geom_boxplot(aes(fill=ind),notch=T,lwd=1.2) +     geom_point(data=data.frame(x=factor(c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19","X20","X21","X22","X23","X24","X25","X26","X27","X28","X29","X30","X31","X32","X33","X34","X35","X36","X37")),y=animals_ranks[,762]),aes(x=x,y=y,pch="Rank from mean over gene trees"),color='darkred',size=(6),show.legend=T) +    geom_point(data=data.frame(x=factor(c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19","X20","X21","X22","X23","X24","X25","X26","X27","X28","X29","X30","X31","X32","X33","X34","X35","X36","X37")),y=animals_ranks[,763]),aes(x=x,y=y,pch="Rank from species tree"),color='darkred',size=(6),show.legend=T) +
    theme_minimal() + 
    labs(x="Species",y="Rank")  +
    scale_x_discrete(labels=c("X1"="R. norvegicus", "X2"="R. fuscipes", "X3"="M. musculus", "X4"="C. gliroides", "X5"="H. minahassae", "X6"="M. major", "X7"="L. nouhuysi",
"X8"="P. macrourus", "X9"="H. goliath", "X10"="C. ruemmleri", "X11"="X. barbatus", "X12"="M. lanosus", "X13"="M. rothschildi", "X14"="P. sp", 
"X15"="L. elegans", "X16"="X. myoides", "X17"="P. ellermani", "X18"="P. asper", "X19"="H. chrysogaster", "X20"="C. moncktoni", "X21"="Z. argurus",
"X22"="L. forresti", "X23"="M. fuscus", "X24"="P. oralis", "X25"="P. fumeus", "X26"="P. desertor", "X27"="P. australis", "X28"="M. fuscus",
"X29"="P. albocinereus", "X30"="M. gouldii", "X31"="L. conditor", "X32"="C. penicillatus", "X33"="U. caudimaculatus", "X34"="P. levipes",
"X35"= "S. salebrosus", "X36"="M. rufescens", "X37"="A. imitator")) 
    #scale_y_continuous(limits=c(0,9))
                
print(animals_plot + mynamestheme + scale_shape_manual(values=c('Rank from mean over gene trees' = "circle", 'Rank from species tree' = "triangle"), name = ' ') + guides(fill = "none") + theme(legend.position="top", legend.text = element_text(face = "bold"),axis.text.x=element_text(size=(12),angle=45,hjust=1,vjust=1)) + theme(legend.background = element_rect(fill="NA",size = 0.5, linetype="solid",colour="darkblue")))


animals_cor = cor(animals_ranks[,1:761],method="kendall")
pairwise = lower.tri(animals_cor,diag=F)
animals_cormat = animals_cor[pairwise]
mean(animals_cormat)

animals_gt_st_cor = matrix(0,ncol=1,nrow=761)
for (i in 1:761) animals_gt_st_cor[i,1] = cor(animals_ranks[,i],animals_ranks[,763],method="kendall")
mean(animals_gt_st_cor)

cor(animals_ranks[,762],animals_ranks[,763],method="kendall")




###########################
# 9. plant data

plants = read.table("plant_fp_combined_table_NA.txt",header=T)

ncol = dim(plants)[2]

plants_ranks = matrix(rank(-1*plants[,2],ties.method="min",na.last="keep"),ncol=1)
for (i in 3: ncol) plants_ranks = cbind(plants_ranks,rank(-1*plants[,i],ties.method="min",na.last="keep"))

dim(plants_ranks)[2]
#320

plants_iqr=apply(plants_ranks[,1:318],1,IQR,na.rm=T)
summary(plants_iqr)

temp_plants = plants_ranks[,1:318]
plants_ranks_t = data.frame(t(temp_plants))

plants_plot <- ggplot(stack(plants_ranks_t),aes(x=ind,y=values)) + 
    geom_boxplot(aes(fill=ind),notch=T,lwd=1.2) + 
    geom_point(data=data.frame(x=factor(c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19","X20","X21","X22","X23","X24","X25","X26","X27","X28","X29","X30","X31","X32","X33","X34","X35","X36","X37","X38","X39","X40","X41","X42","X43","X44","X45","X46","X47","X48","X49","X50","X51","X52")),y=plants_ranks[,319]),aes(x=x,y=y,pch="Rank from mean over gene trees"),color='darkred',size=(6),show.legend=T) + 
    geom_point(data=data.frame(x=factor(c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19","X20","X21","X22","X23","X24","X25","X26","X27","X28","X29","X30","X31","X32","X33","X34","X35","X36","X37","X38","X39","X40","X41","X42","X43","X44","X45","X46","X47","X48","X49","X50","X51","X52")),y=plants_ranks[,320]),aes(x=x,y=y,pch="Rank from species tree"),color='darkred',size=(6),show.legend=T) +
    theme_minimal() + 
    labs(x="Species",y="Rank")  +
    scale_x_discrete(labels=c("X1" = "A. pectinata","X2" = "L. tibetica","X3" = "P. volubilis","X4" = "P. tomentosa","X5" = "C. americana","X6" = "P. lasianthos","X7" = "W. fruticosa","X8" = "C. canadensis","X9" = "P. frutescens","X10" = "H. suaveolens", "X11" = "P. barbatus","X12" = "O. basilicum","X13" = "L. angustifolia","X14" = "M. officinalis","X15" = "P. atriplicifolia","X16" = "R. officinalis","X17" = "S. officinalis","X18" = "S. hispanica","X19" = "L. americanus","X20" = "M. piperita", "X21" = "M. spicata","X22" = "M. didyma","X23" = "O. majorana","X24" = "O. vulgare","X25" = "T. vulgaris","X26" = "P. vulgaris","X27" = "N. cataria","X28" = "N. mussinii","X29" = "G. hederacea","X30" = "H. officinalis", "X31" = "A. foeniculum","X32" = "A. reptans","X33" = "T. canadense","X34" = "C. bungei","X35" = "R. myricoides", "X36" = "B. pseudodictamnus","X37" = "M. vulgare","X38" = "L. album","X39" = "L. leonurus","X40" = "L. cardiaca","X41" = "P. fruticosa","X42" = "B. officinalis","X43" = "P. cablin","X44" = "P. bambusetorum","X45" = "H. sanguinea", "X46" = "S. baicalensis","X47" = "C. pyramidata","X48" = "C. tomentosa","X49" = "V. a. castus","X50" = "G. philippensis", "X51" = "P. microphylla", "X52" = "T. grandis")) 
            
print(plants_plot + mynamestheme + scale_shape_manual(values=c('Rank from mean over gene trees' = "circle", 'Rank from species tree' = "triangle"), name = ' ') + guides(fill = "none") + theme(legend.position="top", legend.text = element_text(face = "bold"),axis.text.x=element_text(size=(12),angle=45,hjust=1,vjust=1)) + theme(legend.background = element_rect(fill="NA",size = 0.5, linetype="solid",colour="darkblue")))


plants_cor = cor(plants_ranks[,1:318],method="kendall",use="pairwise.complete.obs")
pairwise = lower.tri(plants_cor,diag=F)
plants_cormat = plants_cor[pairwise]

plants_gt_st_cor = matrix(0,ncol=1,nrow=318)
for (i in 1:318) plants_gt_st_cor[i,1] = cor(plants_ranks[,i],plants_ranks[,320],method="kendall",use="pairwise.complete.obs")

cor(plants_ranks[,319],plants_ranks[,320],method="kendall",use="pairwise.complete.obs")




##################################

# Plot of all data together

all_gt_st_cor = rbind(primate_gt_st_cor,rattle_gt_st_cor,yeast_gt_st_cor,dolphin_gt_st_cor,mammals_gt_st_cor,snakes_gt_st_cor,fungi_gt_st_cor,animals_gt_st_cor,plants_gt_st_cor)
allcors = data.frame(all_gt_st_cor)
allcors$data = c(rep("Primate",length(primate_gt_st_cor)),rep("Rattlesnake",length(rattle_gt_st_cor)),rep("Yeast",length(yeast_gt_st_cor)),rep("Dolphin",length(dolphin_gt_st_cor)),rep("Mammal",length(mammals_gt_st_cor)),rep("Snake",length(snakes_gt_st_cor)),rep("Fungi",length(fungi_gt_st_cor)),
rep("Rodent",length(animals_gt_st_cor)),rep("Plant",length(plants_gt_st_cor)))
	
	
	
primate_cormat_use = matrix(primate_cormat,ncol=1)	
rattle_cormat_use = matrix(rattle_cormat,ncol=1)	
yeast_cormat_use = matrix(yeast_cormat,ncol=1)	
dolphin_cormat_use = matrix(dolphin_cormat,ncol=1)
mammals_cormat_use = matrix(mammals_cormat,ncol=1)		
snakes_cormat_use = matrix(snakes_cormat,ncol=1)
fungi_cormat_use = matrix(fungi_cormat,ncol=1)	
animals_cormat_use = matrix(animals_cormat,ncol=1)	
plants_cormat_use = matrix(plants_cormat,ncol=1)	
all_cormat = rbind(primate_cormat_use,rattle_cormat_use,yeast_cormat_use,dolphin_cormat_use,mammals_cormat_use,snakes_cormat_use,fungi_cormat_use,animals_cormat_use,plants_cormat_use)	
allcormat = data.frame(all_cormat)
allcormat$data = c(rep("Primate",length(primate_cormat_use)),rep("Rattlesnake",length(rattle_cormat_use)),rep("Yeast",length(yeast_cormat_use)),rep("Dolphin",length(dolphin_cormat_use)),rep("Mammal",length(mammals_cormat_use)),rep("Snake",length(snakes_cormat_use)),rep("Fungi",length(fungi_cormat_use)),rep("Rodent",length(animals_cormat_use)),rep("Plant",length(plants_cormat_use)))


alldat = rbind(as.matrix(allcors),as.matrix(allcormat))
alldat = data.frame(alldat)
alldat$which = c(rep("Gene and species",length(all_gt_st_cor)),rep("Pairwise genes",length(all_cormat)))
alldat$data = as.factor(alldat$data)
alldat$which = as.factor(alldat$which)
alldat$all_gt_st_cor = as.numeric(alldat$all_gt_st_cor)



plot1 <- ggplot(alldat,aes(x=data,y=all_gt_st_cor,fill=which,colour=which)) +
  geom_boxplot(color="black",notch=T,lwd=1.2, position = position_dodge(width = 0.78)) +
  scale_x_discrete(name = "", label=c("Dolphin","Fungi","Mammal","Plant","Primate","Rattlesnake","Rodent","Snake","Yeast")) +
  theme_minimal() +
  labs(x="",y="Correlation")
  
   
print(plot1 + mynamestheme + theme(axis.text.x=element_text(size=(18))) + theme(legend.position=c(0.15,0.2), legend.text = element_text(face = "bold",size=(16)), legend.title=element_blank()) + theme(legend.background = element_rect(fill="NA",size = 0.5, linetype="solid",colour="darkblue")))
  

