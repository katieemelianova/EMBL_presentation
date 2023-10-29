library(dplyr)
library(readr)
library(magrittr)

# prior to this I ran the Snakemake pipeline dyno to find and extract LTR sequences from Copia gypsy elements
# the resulting files *all were concatenated and moved locally in files all_hist_impolita and all_hist_sandwicensis

##########################################
#       draw histogram of LTR ages        #
##########################################

# read in the files
all_hist_impo <- read.table("all_hist_impolita", col.names = c("scaffold", "family", "ltr")) %>% mutate(species="D. impolita")
all_hist_sand <- read.table("all_hist_sandwicensis", col.names = c("scaffold", "family", "ltr")) %>% mutate(species="D. sandwicensis")

# read in the files
all_gypsy_impo <- read.table("all_gypsy_impolita", col.names = c("scaffold", "family", "ltr")) %>% mutate(species="D. impolita")
all_gypsy_sand <- read.table("all_gypsy_sandwicensis", col.names = c("scaffold", "family", "ltr")) %>% mutate(species="D. sandwicensis")

# bind them together, here I filter out any LTRs which have a distance of more than 0.2
all_copia<-rbind(all_hist_impo, all_hist_sand) %>% filter(ltr < 0.2) %>% mutate(Superfamily="LTR Copia")
all_gypsy <-rbind(all_gypsy_impo, all_gypsy_sand) %>% filter(ltr < 0.5) %>% mutate(Superfamily="LTR Gypsy")

# bind all that again by superfamily
all<-rbind(all_copia, all_gypsy)

# make a hostogram of LTR distances for both species, 
png(file="EMBL_histogram_insertion_time.png", width = 1000, height = 1000)
all %>%
  ggplot( aes(x=ltr, fill=species)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity', bins=60) +
  #scale_fill_manual(values=c("#69b3a2", "#404080")) +
  scale_fill_manual(values=c("red", "blue")) +
  labs(fill="") +
  xlab("5' and 3' LTR divergence") +
  ylab("Number of TEs") +
  theme(text = element_text(size = 35),
        legend.text = element_text(face = "bold"),
        #axis.title.y=element_blank(),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  facet_wrap(~Superfamily, ncol = 1)
dev.off()


##################################################
#     draw rooted highltighted species tree      #
################################################## 

species_tree<-ape::read.tree("/Users/katieemelianova/Desktop/Diospyros/diospyros_plots/lib1234_speciestree.nwk")
species_tree.rooted <- root(species_tree, which(species_tree$tip.label == "D.sandwicensis"))
tip<-c("D.olen", "D.fasciculosa", "D.macrocarpa", "D.ferrea")
species_tree.rooted<-drop.tip(species_tree.rooted, tip)
species_tree.rooted$species <- species_tree.rooted$tip.label


colours_tips <- case_when(species_tree.rooted$tip.label == "D.sandwicensis" ~ "red",
                          species_tree.rooted$tip.label == "D.impolita" ~"blue",
                          !(species_tree.rooted$tip.label %in% c("D.sandwicensis", "D.impolita")) ~ "black")


dd <- data.frame(taxa=species_tree.rooted$tip.label, tipcols=colours_tips)
p<-ggtree(species_tree.rooted, size=3)
p <- p %<+% dd

png("EMBL_species_highlight_tree.png", width = 1500, height = 1500)
p + geom_tiplab(size=15, aes(color=tipcols)) + 
  scale_colour_manual(values = c("black", "red", "blue")) +
  expand_limits(x = 0.07) + 
  theme_tree2(plot.margin=margin(6, 120, 6, 6)) +
  theme(text=element_text(size=10), 
        legend.position = "none",
        axis.line.x.bottom=element_line(size=3),
        axis.text.x=element_text(size=40, vjust=-1),
        plot.margin = margin(, , 3, , "cm")) + 
  geom_treescale(x=0, y=25, fontsize=15, linesize=1)
dev.off()

# these refer to he plot.margin argument in theme above: which margins refer to which numbers
#up, #right, #down, #left


####################################################################
#     draw scatterplot of TE family size vs LTR insertion time     #
####################################################################

# aggregate by distance (insertion time) and by size of TE family
aggregated_family_dist<-aggregate(all$ltr, list(family = all$family, species=all$species, superfamily=all$Superfamily), mean)
aggregated_family_size<-aggregate(all$ltr, list(family = all$family, species=all$species, superfamily=all$Superfamily), length)

# join them together and do some filtering
aggregated_family<-inner_join(aggregated_family_dist, aggregated_family_size, by = "family") %>% 
  set_colnames(c("family", "species", "superfamily", "distance", "to_delete1", "to_delete2", "size")) %>% 
  dplyr::select(-c(to_delete1, to_delete2)) %>% filter(size < 20 & distance < 0.1)

# plot a facet plot of average insertion time as a function of family size
png("EMBL_TE_vs_FamSize.png", width = 500, height = 1000)
ggplot(aggregated_family, aes(x=distance, y=size, color=species)) +
  geom_point(size=2, shape=19) + geom_smooth(method="lm", se=TRUE, fullrange=TRUE, alpha = .15, aes(fill = species)) +
  xlab("TE subfamily average 5' and 3' LTR divergence") +
  ylab("TE subfamily size") +
  scale_colour_manual(values = c("red", "blue")) +
  facet_wrap(~superfamily, ncol = 1) +
  theme(strip.text = element_text(size = 20),
        legend.text=element_text(size=20),
        legend.title=element_blank(),
        text = element_text(size = 20),
        legend.position = c(.80, .95),
        panel.spacing = unit(3, "lines"),
        legend.background = element_rect(fill='transparent'),
        plot.margin = margin(, 1, , 0.2, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
dev.off()

##################################################
#                   plot TE summary              #
################################################## 

# I here transcribe the numbers from the repeatmasker summary files stored in Diospyros/prsentations/EMBL*summary

san<-data.frame(species="D. sandwicensis",
  total_length=1007252419, 
  total_masked=726814292,
  Ty1Copia=94093519,
  Gypsy=402006824)

# subtract TE counts from categories which contain them
san$total_length <-(san$total_length - (san$total_masked))
san$total_masked <- (san$total_masked - (san$Ty1Copia + san$Gypsy))

imp<-data.frame(species="D. impolita",
  total_length=1834533962,
  total_masked=1388789270,
  Ty1Copia=224869633, 
  Gypsy=833189331)

# subtract TE counts from categories which contain them
imp$total_length <-(imp$total_length - (imp$total_masked))
imp$total_masked <- (imp$total_masked - (imp$Ty1Copia + imp$Gypsy))

# both both species together and reshape the df for plotting
# also convert units to Mb
both<-rbind(san, imp) %>% melt()
both$value<-both$value/1000000

# plot barplot showing proportions of TEs in genome
png("EMBL_TE_summary_barplot.png", width = 1000, height = 1000)
ggplot(both, aes(x=species, y=`value`, fill=variable)) + 
  geom_bar(stat="identity", position="stack") +
  scale_fill_manual(labels= c("Total Genome", 
                              "Total Masked",
                              "LTR Ty1 Copia", 
                              "LTR Gypsy"), 
                    values = c("royalblue1", "hotpink", "orange1", "palegreen3")) +
  theme(legend.title=element_blank(),
        axis.title.x=element_blank(),
        text = element_text(size = 40)) +
  ylab("Total bases (Mbp)") 
dev.off()



######################################################################################
#     Got the Newick tree fron Teerna's tree and trying to plot it with ggtree       #
######################################################################################

# read in original tree
nwk<-ape::read.tree("/Users/katieemelianova/Desktop/Diospyros/diospyros_plots/RadiatingSpeciesDiospyros_ingrp-inds.nwk")

# split the tip name column by BT to get the population names
# I dont understand the format of this function and I dont like it at all
# but it works so whatevs
# also do some data cleaning
species<-sapply(strsplit(nwk$tip.label,"BT"), `[`, 1)
nwk$tip.label[nwk$tip.label =="impolitaBt103"] = "impolitaBT103"
sample_id<-sapply(strsplit(nwk$tip.label,"BT"), `[`, 2)
sample_id<-paste("BT", sample_id, sep = "")
sample_id[sample_id == "BTNA"] = nwk$tip.label[sample_id =="BTNA"]
sample_id<-sapply(strsplit(sample_id, "_"), `[`, 1)
sample_id<-sapply(strsplit(sample_id, "-"), `[`, 1)

# set the species in the newick object
nwk$species <- species
nwk$sample_id<-sample_id



# take the original data frame from the map part and take the rad and BT sample name columns only
mapping_loc<-diospyros_localities %>% 
  dplyr::select(`RAD localities`, `sequenced samples`, Soil) %>% 
  data.frame()

# this is a really dumb way to do it but whatever it works
# apply over each row, split the BTXXX column by ", ", and then add it to a data frame, where the other column is the RAD number
# this gives you a list of data frames, which are bound into one using bind_rows
bt_soil_mapping<-apply(mapping_loc, 1, function(x) strsplit(x[2], ", ") %>% data.frame(rad=x[3])) %>% bind_rows()

# remove the weird rownames the apply function gives it
rownames(bt_soil_mapping) <-NULL
colnames(bt_soil_mapping)<-c("sample_id", "soil")

# remove duplicate samples
bt_soil_mapping_noduplicates <- bt_soil_mapping %>% filter(sample_id != "BT147" & sample_id != "BT296")

# join tip names and soil types
bt_soil_mapping_joined<-left_join(data.frame(sample_id=nwk$sample_id), bt_soil_mapping_noduplicates) %>% dplyr::select(soil)

#set rownames to corresponding tip labels
rownames(bt_soil_mapping_joined) <- nwk$tip.label


# make a list where item name is species name and objects within are the tip labels belonging to that species
groupInfo<-split(nwk$tip.label, nwk$species)

# use groupOTU to group the tips by species and plot
nwk_grouped<-groupOTU(nwk, groupInfo, group_name = "species")

# take the same tree as above but make the tip labels the species name only (no population info)
nwk_grouped_speciesonly<-nwk_grouped
nwk_grouped_speciesonly$tip.label <- nwk_grouped_speciesonly$species


# first make one tree which is grouped and coloured by species
circ<-ggtree(nwk_grouped, aes(color=species), layout='circular') + geom_tiplab(size=5) 
circ_notiplab<-ggtree(nwk_grouped, aes(color=species), layout='circular', size=2.3) + geom_tiplab()
circ_notiplab<-ggtree(nwk_grouped, aes(color=species), layout='circular', size=2.3)


#offset=.012, width=.05
circ_soil<-gheatmap(circ_notiplab, 
                    bt_soil_mapping_joined, 
                    offset=0.0000005, width=.05, colnames=FALSE) + 
  scale_fill_viridis_d(option = "C", name = "Clade", na.value = "gray94") + 
  theme(legend.position="none")

# make the legend for the soiltypes
circ_nolegend<-ggtree(nwk)
soiltype_values<-bt_soil_mapping_joined %>% dplyr::select(soil)

soiltype_values$soil[soiltype_values$soil == "Ironcrust"] <- NA

soiltype_values$soil %>% unique()

# repeat this for soiltype
soiltype_legend<-gheatmap(circ_nolegend, soiltype_values, offset = .0001, width = 2,
                          colnames = FALSE,
                          colnames_offset_y = 1) +
  scale_fill_viridis_d(option = "C", na.value = "gray94", name = "soiltype", 
                       labels = c("Limestone", "Calcareous", "Schist", "Serpentine", "Ultramafic", "Volcanic", "Unknown")
                       #scale_fill_viridis_d(option = "C", na.value = "gray94", name = "Soiltype", labels = c("Ironcrust", "Kalcarious", "Serpentine", "Ultramafic", "Volcano-Sedimentary", "Unknown"))
                       #scale_fill_manual(values=mycolors, name="Soiltype"
  ) + theme(legend.text=element_text(size=50), 
            legend.title= element_blank())

# get the soiltype legend
soiltype_legend <- get_legend(soiltype_legend)

png("EMBL_soiltype_tree.png", width = 1400, height = 1000)
ggarrange(circ_soil, soiltype_legend, widths = c(8,4))
dev.off()



