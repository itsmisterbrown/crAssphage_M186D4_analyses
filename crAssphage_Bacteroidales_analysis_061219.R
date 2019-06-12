library(ggplot2)
library(ggsignif)
library(phyloseq)
library(cowplot)
source("/Users/bpb/Documents/GitHub/microfiltR/microfiltR_source_code.R")

#load data
inf.vs <- phyloseq::import_biom("/Users/bpb/Desktop/InFANT_virome/Bacterial_16S/analysis/Virome_subsample/dataset_files/Complete_virome_subset_013119_ASV_table_w_tax_md.biom")
#update colnames
colnames(tax_table(inf.vs)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#change theme
theme_set(theme_bw())

#update sample data 
sample_data(inf.vs)$Time_point <- factor(sample_data(inf.vs)$Time_point, levels = c("D4.7", "WK.04", "WK.15", "WK.36"), labels = c("Week 1", "Week 4", "Week 15", "Week 36"))
sample_data(inf.vs)$exposure <- c("Infected/exposed", "Infected/exposed", "Infected/exposed", "Infected/exposed", "Infected/exposed", "Infected/exposed", "Infected/exposed", "Uninfected/unexposed", "Uninfected/unexposed", "Uninfected/unexposed", "Uninfected/unexposed")
sample_data(inf.vs)$dyad <- c("186", "197", "519", "521", "521", "526", "526", "701", "720", "720", "720")
sample_data(inf.vs)$FC <- c(3.75, NA, 0.05, 11.26, 19.15, 0.084, 0.088, NA, 0.36, NA, 0.07)
sample_data(inf.vs)$FCbin <- ifelse(test = sample_data(inf.vs)$FC > 1, yes = "> 1", no = "< 1")

#subset to remove samples with no cAssphage data
inf.vs.crass <- microfilter(ps = inf.vs, mdCAT = "crass_cov", mdFACTOR = "NA")
#create clr transformed object
inf.vs.crass.clr <- transform_sample_counts(physeq = inf.vs.crass, fun = function(x) x+1)
#transpose
otu_table(inf.vs.crass.clr) <- clr(t(otu_table(inf.vs.crass.clr)))
#subset to Bacteroidales
inf.vs.crass.clrB <- subset_taxa(inf.vs.crass.clr, Order=="Bacteroidales")

#create alias
phy.obj <- inf.vs.crass.clrB

#add annotation to unannotated bugs
tax_table(phy.obj)[23,6] <- "Bacteroidales"
tax_table(phy.obj)[23,7] <- "sp."
tax_table(phy.obj)[47,7] <- "sp."
tax_table(phy.obj)[41,7] <- "sp."

#melt df
po.bm <- psmelt(phy.obj)

#convert to separate dfs
asv.tab <- otu_table(phy.obj)
tax.tab <- tax_table(phy.obj)
sampledf <- sample_data(phy.obj)

#check correlations
#fold coverage
#convert to numeric
sampledf$FC <- as.numeric(paste0(sampledf$FC))
casv.tab <- cbind.data.frame(asv.tab, log10(sampledf$FC))
cor.tab <- cor(casv.tab, method = "pearson")
View(cor.tab)

#subset to txa with cor coefs > 0.4 and don't include crass coverage column
pc <- rownames(cor.tab[which(cor.tab[,58] > 0.4),])[1:6] #don't include log10(sampledf$FC)
nc <- rownames(cor.tab[which(cor.tab[,58] < -0.4),])

#create combined df
cor.combined <- c(pc, nc)
#subset
po.bm.cor4 <- po.bm[which(po.bm$OTU %in% cor.combined),]

#update factor levels
namevec <- c("Alistipes\ninops", "Parabacteroides\nmerdae", "Paraprevotella\nsp.", "Alistipes\nputredinis", "Parabacteroides\nmerdae", "Bacteroidales\nsp.",  "Parabacteroides\ndistasonis",  "Parabacteroides\nsp.", "Parabacteroides\ndistasonis", "Bacteroides\nuniformis", "Bacteroides\nthetaiotaomicron")
po.bm.cor4$OTU <- factor(po.bm.cor4$OTU, levels = c("ASV224", "ASV186", "ASV258", "ASV40", "ASV55", "ASV140",  "ASV287",  "ASV301", "ASV71", "ASV70", "ASV82"))
po.bm.cor4$annotatedASVs <- paste(po.bm.cor4$Genus, po.bm.cor4$Species, sep = "\n")

#plot the CLR abundances
ggplot(po.bm.cor4, aes(x = OTU, y = Abundance, fill = FCbin)) +
  guides(ncol = 1) +
  geom_boxplot(size=1.5) +
  theme(axis.text.y = element_text(size = 20, colour = "black"),
        axis.title.y = element_text(size = 20, colour = "black"),
        axis.text.x = element_text(size = 10, colour = "black", angle = 315),
        axis.title.x = element_text(size = 0, colour = "black"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 0, colour = "black"),
        legend.position = "right",
        strip.background =element_rect(fill="white"), 
        strip.text = element_text(size = 20, colour = "black"))  +
  labs(y="CLR abundance")


#create ILR balance for simplified model fitting
#add pseudocount
inf.vs.crass.p1 <- transform_sample_counts(physeq = inf.vs.crass, fun = function(x) x+1)

B.bal <- create.balance(df = t(otu_table(inf.vs.crass.p1)), num.tax = "ASV82", den.tax = "ASV55") 

#create df for pairs plots and regression
cdf <- cbind.data.frame(asv.tab, sampledf, B.bal)

#statistical tests on ASV abundance by crassphage FC
#wilcox
wilcox.test(data=cdf, ASV55~FCbin)

#t.test
#subset to hi/lo
hiFC <- subset(cdf, FCbin=="> 1")
loFC <- subset(cdf, FCbin=="< 1")

#check variance
var.test(hiFC$ASV55, loFC$ASV55) #looks good

#perform unpaired two-sample t-test
t.test(ASV55~FCbin, data = cdf, var.equal=T)

#GLM crassphage FC
lmt <- glm(formula = sigvec~log10(FC), data = cdf, family = gaussian(link = identity))

lmt
summary(lmt)

#plot model fit
#gaussian
ggplot(cdf, aes(x=FC, y=sigvec)) +
  geom_point(size=5, aes(color=exposure)) +
  geom_point(color = "grey90", size = 2) +
  theme(axis.text.y = element_text(size = 25, colour = "black"),
        axis.title.y = element_text(size = 20, colour = "black", face = "italic"),
        axis.text.x = element_text(size = 25, colour = "black"),
        axis.title.x = element_text(size = 20, colour = "black"),
        legend.text = element_text(size = 20, colour = "black"),
        legend.title = element_text(size = 0),
        legend.position=c(0.25, 0.75)) +
  scale_x_log10() +
  geom_smooth(method = "glm", data = cdf, formula = y~x, method.args = list(family = gaussian(link = identity)), se = TRUE) + 
  labs(y=expression(paste(log(italic("P.merdae")/italic("B.thetaiotaomicron")))), x="crAssphage FC")

#put it all together in single plot
#with stat comparisons
plot1 <- ggplot(po.bm.cor4, aes(x = OTU, y = Abundance, fill = FCbin)) +
  guides(ncol = 1, fill=guide_legend(title = "crAssphage coverage")) +
  geom_boxplot(size=1.5) +
  theme(axis.text.y = element_text(size = 20, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 10, colour = "black", angle = 315, vjust = 0.9),
        axis.title.x = element_text(size = 0, colour = "black"),
        legend.text = element_text(size = 12, colour = "black"),
        legend.title = element_text(size = 14, colour = "black"),
        legend.position = "right",
        legend.key.size = unit(2, "line"))  +
  labs(y="Bacteroidales abundance (CLR)") + scale_x_discrete(labels=namevec) + 
  geom_signif(y_position=c(6.5, 8.25), xmin=c(4.7, 10.7), xmax=c(5.3, 11.3),
              annotation=c("*", "*"), textsize = 10, size=2, color="black", vjust = 0.5) +
  lims(y=c(-4,8.5))

plot2 <- ggplot(cdf, aes(x=FC, y=sigvec)) +
  geom_point(size=5, aes(color=exposure)) +
  geom_point(color = "grey90", size = 2) +
  theme(axis.text.y = element_text(size = 20, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black", face = "italic"),
        axis.text.x = element_text(size = 20, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black"),
        legend.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 14),
        legend.position="right") +
  scale_x_log10() +
  scale_color_manual(values = c("forestgreen", "orange")) +
  geom_smooth(method = "glm", data = cdf, formula = y~x, method.args = list(family = gaussian(link = identity)), se = TRUE) + 
  labs(y=expression(paste(log(italic(over("B.thetaiotaomicron","P.merdae"))))), x="crAssphage coverage") +
  guides(color=guide_legend(title="HIV status"))

#plot it
cowplot::plot_grid(plot1, plot2, labels = "AUTO", label_size = 25, ncol=1, align = "v")

