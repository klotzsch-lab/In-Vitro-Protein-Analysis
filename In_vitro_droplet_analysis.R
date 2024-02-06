### R script written by Simon Mergenthal
### Humboldt-Universität zu Berlin
### Mechanobiology (Prof. Klotzsch)
### 13.05.2023

### last updated: 06.02.2024

# Script to plot granule data after particle analysis
# Version 2.0 working with all analyzed data with different amounts of parameters

# getting work environment clear
rm(list=ls())
while (dev.cur() !=1) {dev.off()}

### Please enter path to folder, where "Results_particles.csv" after running In_vitro_droplet_analysis.py is located.
path <- "..."


setwd(path)

### loading of packages
pacman::p_load(openxlsx, tidyverse, scales)

options(stringsAsFactors = FALSE) # changes the option, to automatically structure strings as factors

aggregates <- read.csv('Results_particles.csv', sep = "\t", fileEncoding = "utf-8")

# extracting parameter information from image label
aggregates$Protein <- str_extract(aggregates$Label, "FLJ|FLV|FL|DPLDFJ|DPLD|PLDJ|PLDo|PLD|mNG") #for all variants
aggregates$Protein_conc <- str_extract(aggregates$Label, "12.5(?=uM)|\\d+(?=uM)")
aggregates$Protein_conc <- factor(aggregates$Protein_conc , levels=c("1", "5", "12.5", "25", "50", "100"))
aggregates$NaCl_conc <- str_extract(aggregates$Label, "\\d+(?=mM)")
aggregates$NaCl_conc <- factor(aggregates$NaCl_conc, levels = c("50", "100", "150", "250"))
aggregates$Temperature <- str_extract(aggregates$Label, "RT|cold")
aggregates$Temperature[aggregates$Temperature == "cold"] <- "4C"
aggregates$RNA_conc <- str_extract(aggregates$Label, "12.5(?=ng)|\\d+(?=ng)")
aggregates$RNA_conc <- factor(aggregates$RNA_conc, levels=c("0","1","5","12.5","25","50","100","200"))

aggregates <- aggregates[aggregates$Mean>5,]

### Temperature

data_temperature <- aggregates[!is.na(aggregates[,"Temperature"]),]

data_temperature$Temperature <- factor(data_temperature$Temperature, levels = c("RT", "4C"))

conditions_Temperature <- data_temperature %>% group_by(Protein, Temperature) %>% summarise(Amount=n())

means_temp <- aggregate(data_temperature$Feret, FUN = mean, by = list(data_temperature$Protein, data_temperature$Temperature))
names(means_temp) <- c("Protein","Temperature","mean")
means_temp$sd <- aggregate(data_temperature$Feret, FUN = sd, by = list(data_temperature$Protein, data_temperature$Temperature))[,3]
means_temp$median <- aggregate(data_temperature$Feret, FUN = median, by = list(data_temperature$Protein, data_temperature$Temperature))[,3]


ggplot(data = data_temperature, aes(x = Temperature, y = Feret, fill = Temperature, color = Temperature))+
   geom_violin(trim = TRUE)+
   facet_wrap(~Protein)+
   scale_y_continuous(name = "droplet diameter [µm]", breaks = seq(0,20,0.2))+
   coord_cartesian(ylim = c(0.40,8.60))+
   scale_fill_manual(values = c("#70AD45","#5C9BD3"))+
   scale_color_manual(values = c("#5A8A3D","#4A7DAA"))+
   stat_summary(fun = median, fun.min = median, fun.max = median,
                geom = "crossbar", 
                width = 0.4,
                position = position_dodge(width = 0.4),
                show.legend = FALSE)+
   theme(panel.background = element_rect(fill = NA, color = "grey"),
         panel.grid.major.x = element_blank(),
         panel.grid.major.y = element_line(color = c("grey", NA, NA, NA, NA)),
         panel.grid.minor = element_blank(),
         axis.title = element_text(size = 18),
         axis.title.x = element_blank(),
         axis.text = element_text(size = 15),
         axis.text.y = element_text(color = c("#8c8c8c", NA, NA, NA, NA)),
         axis.ticks.y =element_line(color = "grey"),
         plot.title = element_text(size=20),
         legend.text = element_text(size=15),
         legend.title = element_text(size=18),
         strip.text = element_text(size=18),
         legend.position = "bottom",)

### significance test
t_temp_DPLD_4C <- subset(data_temperature, Protein == "DPLD" & Temperature == "4C")
t_temp_DPLD_RT <- subset(data_temperature, Protein == "DPLD" & Temperature == "RT")

t_temp_FL_4C <- subset(data_temperature, Protein == "FL" & Temperature == "4C")
t_temp_FL_RT <- subset(data_temperature, Protein == "FL" & Temperature == "RT")

t_temp_PLD_4C <- subset(data_temperature, Protein == "PLD" & Temperature == "4C")
t_temp_PLD_RT <- subset(data_temperature, Protein == "PLD" & Temperature == "RT")


df_pval_Temperature <- data.frame(matrix(ncol = 3, nrow = 0)) # starting with an empty data frame with 3 columns
colnames(df_pval_Temperature) <-c("condition 1", "condition 2", "significance")

df_pval_Temperature[nrow(df_pval_Temperature) + 1,] <- c("DPLD_4C", "DPLD_RT", wilcox.test(t_temp_DPLD_4C$Feret,t_temp_DPLD_RT$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Temperature[nrow(df_pval_Temperature) + 1,] <- c("FL_4C", "FL_RT", wilcox.test(t_temp_FL_4C$Feret,t_temp_FL_RT$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Temperature[nrow(df_pval_Temperature) + 1,] <- c("PLD_4C", "PLD_RT", wilcox.test(t_temp_PLD_4C$Feret,t_temp_PLD_RT$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)

df_pval_Temperature$significance <- as.numeric(df_pval_Temperature$significance)
df_pval_Temperature <- df_pval_Temperature %>% mutate(sign = case_when(significance == 0 ~ "NA",
                                                                       significance < 0.0001 ~ "****",
                                                                       significance < 0.001 ~ "***",
                                                                       significance < 0.01 ~ "**",
                                                                       significance < 0.05 ~ "*"))

write.csv2(data_temperature, "data_maturation_temperature.csv", row.names=FALSE)
write.csv2(means_temp, "data_maturation_temperature_desciptives.csv", row.names=FALSE)
write.csv2(conditions_Temperature, "data_maturation_temperature_n.csv", row.names=FALSE)
write.csv(df_pval_Temperature, "data_maturation_temperature_significances.csv", row.names=FALSE)

### Protein concentration

data_Protein <- aggregates[!is.na(aggregates[,"Protein_conc"]),]
data_Protein$Protein[data_Protein$Protein == "FLV"] <- "FL"

conditions_Protein <- data_Protein %>% group_by(Protein, Protein_conc) %>% summarise(Amount=n())

means_Protein <- aggregate(data_Protein$Feret, FUN = mean, by = list(data_Protein$Protein_conc, data_Protein$Protein))
names(means_Protein) <- c("Protein_conc","Protein","mean")
means_Protein$sd <- aggregate(data_Protein$Feret, FUN = sd, by = list(data_Protein$Protein_conc, data_Protein$Protein))[,3]
means_Protein$median <- aggregate(data_Protein$Feret, FUN = median, by = list(data_Protein$Protein_conc, data_Protein$Protein))[,3]

ggplot(data = data_Protein, aes(x = Protein_conc, y = Feret))+
   geom_violin(trim=TRUE)+
   facet_wrap(~Protein)+
   scale_y_continuous(name = "droplet diameter [µm]", breaks = seq(0,20,1))+
   stat_summary(fun = median, fun.min = median, fun.max = median,
                geom = "crossbar", 
                width = 0.4,
                position = position_dodge(width = 0.4),
                show.legend = FALSE)+
   coord_cartesian(ylim = c(0.77,16.23))+
   scale_x_discrete(name= "Protein concentration [µM]")+
   theme(panel.background = element_rect(fill = NA, color = "grey"),
         panel.grid.major.x = element_blank(),
         panel.grid.major.y = element_line(color = c("grey", NA, NA, NA, NA)),
         panel.grid.minor = element_blank(),
         axis.title = element_text(size = 18),
         axis.text = element_text(size = 13),
         axis.text.y = element_text(color = c("#8c8c8c", NA, NA, NA, NA)),
         axis.ticks.y =element_line(color = "grey"),
         plot.title = element_text(size=20),
         legend.text = element_text(size=15),
         legend.title = element_text(size=18),
         strip.text = element_text(size=16),
         legend.position = "bottom",)

### significances

t_Protein_1uM_FL <- subset(data_Protein, Protein == "FL" & Protein_conc == "1")
t_Protein_5uM_FL <- subset(data_Protein, Protein == "FL" & Protein_conc == "5")
t_Protein_12.5uM_FL <- subset(data_Protein, Protein == "FL" & Protein_conc == "12.5")
t_Protein_25uM_FL <- subset(data_Protein, Protein == "FL" & Protein_conc == "25")
t_Protein_50uM_FL <- subset(data_Protein, Protein == "FL" & Protein_conc == "50")
t_Protein_100uM_FL <- subset(data_Protein, Protein == "FL" & Protein_conc == "100")

t_Protein_1uM_PLDo <- subset(data_Protein, Protein == "PLDo" & Protein_conc == "1")
t_Protein_5uM_PLDo <- subset(data_Protein, Protein == "PLDo" & Protein_conc == "5")
t_Protein_12.5uM_PLDo <- subset(data_Protein, Protein == "PLDo" & Protein_conc == "12.5")
t_Protein_25uM_PLDo <- subset(data_Protein, Protein == "PLDo" & Protein_conc == "25")
t_Protein_50uM_PLDo <- subset(data_Protein, Protein == "PLDo" & Protein_conc == "50")
t_Protein_100uM_PLDo <- subset(data_Protein, Protein == "PLDo" & Protein_conc == "100")

t_Protein_1uM_DPLD <- subset(data_Protein, Protein == "DPLD" & Protein_conc == "1")
t_Protein_5uM_DPLD <- subset(data_Protein, Protein == "DPLD" & Protein_conc == "5")
t_Protein_12.5uM_DPLD <- subset(data_Protein, Protein == "DPLD" & Protein_conc == "12.5")
t_Protein_25uM_DPLD <- subset(data_Protein, Protein == "DPLD" & Protein_conc == "25")
t_Protein_50uM_DPLD <- subset(data_Protein, Protein == "DPLD" & Protein_conc == "50")
t_Protein_100uM_DPLD <- subset(data_Protein, Protein == "DPLD" & Protein_conc == "100")

t_Protein_1uM_mNG <- subset(data_Protein, Protein == "mNG" & Protein_conc == "1")
t_Protein_5uM_mNG <- subset(data_Protein, Protein == "mNG" & Protein_conc == "5")
t_Protein_12.5uM_mNG <- subset(data_Protein, Protein == "mNG" & Protein_conc == "12.5")
t_Protein_25uM_mNG <- subset(data_Protein, Protein == "mNG" & Protein_conc == "25")
t_Protein_50uM_mNG <- subset(data_Protein, Protein == "mNG" & Protein_conc == "50")
t_Protein_100uM_mNG <- subset(data_Protein, Protein == "mNG" & Protein_conc == "100")

df_pval_Protein <- data.frame(matrix(ncol = 3, nrow = 0)) # starting with an empty data frame with 3 columns
colnames(df_pval_Protein) <-c("condition 1", "condition 2", "significance")

df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("FL_1uM", "FL_5uM", wilcox.test(t_Protein_1uM_FL$Feret,t_Protein_5uM_FL$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("FL_1uM", "FL_12.5uM", wilcox.test(t_Protein_1uM_FL$Feret,t_Protein_12.5uM_FL$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("FL_1uM", "FL_25uM", wilcox.test(t_Protein_1uM_FL$Feret,t_Protein_25uM_FL$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("FL_1uM", "FL_50uM", wilcox.test(t_Protein_1uM_FL$Feret,t_Protein_50uM_FL$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("FL_1uM", "FL_100uM", wilcox.test(t_Protein_1uM_FL$Feret,t_Protein_100uM_FL$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("FL_5uM", "FL_12.5uM", wilcox.test(t_Protein_5uM_FL$Feret,t_Protein_12.5uM_FL$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("FL_5uM", "FL_25uM", wilcox.test(t_Protein_5uM_FL$Feret,t_Protein_25uM_FL$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("FL_5uM", "FL_50uM", wilcox.test(t_Protein_5uM_FL$Feret,t_Protein_50uM_FL$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("FL_5uM", "FL_100uM", wilcox.test(t_Protein_5uM_FL$Feret,t_Protein_100uM_FL$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("FL_12.5uM", "FL_25uM", wilcox.test(t_Protein_12.5uM_FL$Feret,t_Protein_25uM_FL$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("FL_12.5uM", "FL_50uM", wilcox.test(t_Protein_12.5uM_FL$Feret,t_Protein_50uM_FL$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("FL_12.5uM", "FL_100uM", wilcox.test(t_Protein_12.5uM_FL$Feret,t_Protein_100uM_FL$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("FL_25uM", "FL_50uM", wilcox.test(t_Protein_25uM_FL$Feret,t_Protein_50uM_FL$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("FL_25uM", "FL_100uM", wilcox.test(t_Protein_25uM_FL$Feret,t_Protein_100uM_FL$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("FL_50uM", "FL_100uM", wilcox.test(t_Protein_50uM_FL$Feret,t_Protein_100uM_FL$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)

df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("PLDo_1uM", "PLDo_5uM", wilcox.test(t_Protein_1uM_PLDo$Feret,t_Protein_5uM_PLDo$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("PLDo_1uM", "PLDo_12.5uM", wilcox.test(t_Protein_1uM_PLDo$Feret,t_Protein_12.5uM_PLDo$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("PLDo_1uM", "PLDo_25uM", wilcox.test(t_Protein_1uM_PLDo$Feret,t_Protein_25uM_PLDo$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("PLDo_1uM", "PLDo_50uM", wilcox.test(t_Protein_1uM_PLDo$Feret,t_Protein_50uM_PLDo$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("PLDo_1uM", "PLDo_100uM", wilcox.test(t_Protein_1uM_PLDo$Feret,t_Protein_100uM_PLDo$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("PLDo_5uM", "PLDo_12.5uM", wilcox.test(t_Protein_5uM_PLDo$Feret,t_Protein_12.5uM_PLDo$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("PLDo_5uM", "PLDo_25uM", wilcox.test(t_Protein_5uM_PLDo$Feret,t_Protein_25uM_PLDo$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("PLDo_5uM", "PLDo_50uM", wilcox.test(t_Protein_5uM_PLDo$Feret,t_Protein_50uM_PLDo$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("PLDo_5uM", "PLDo_100uM", wilcox.test(t_Protein_5uM_PLDo$Feret,t_Protein_100uM_PLDo$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("PLDo_12.5uM", "PLDo_25uM", wilcox.test(t_Protein_12.5uM_PLDo$Feret,t_Protein_25uM_PLDo$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("PLDo_12.5uM", "PLDo_50uM", wilcox.test(t_Protein_12.5uM_PLDo$Feret,t_Protein_50uM_PLDo$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("PLDo_12.5uM", "PLDo_100uM", wilcox.test(t_Protein_12.5uM_PLDo$Feret,t_Protein_100uM_PLDo$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("PLDo_25uM", "PLDo_50uM", wilcox.test(t_Protein_25uM_PLDo$Feret,t_Protein_50uM_PLDo$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("PLDo_25uM", "PLDo_100uM", wilcox.test(t_Protein_25uM_PLDo$Feret,t_Protein_100uM_PLDo$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("PLDo_50uM", "PLDo_100uM", wilcox.test(t_Protein_50uM_PLDo$Feret,t_Protein_100uM_PLDo$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)

df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("DPLD_1uM", "DPLD_5uM", wilcox.test(t_Protein_1uM_DPLD$Feret,t_Protein_5uM_DPLD$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("DPLD_1uM", "DPLD_12.5uM", wilcox.test(t_Protein_1uM_DPLD$Feret,t_Protein_12.5uM_DPLD$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("DPLD_1uM", "DPLD_25uM", wilcox.test(t_Protein_1uM_DPLD$Feret,t_Protein_25uM_DPLD$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("DPLD_1uM", "DPLD_50uM", wilcox.test(t_Protein_1uM_DPLD$Feret,t_Protein_50uM_DPLD$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("DPLD_1uM", "DPLD_100uM", wilcox.test(t_Protein_1uM_DPLD$Feret,t_Protein_100uM_DPLD$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("DPLD_5uM", "DPLD_12.5uM", wilcox.test(t_Protein_5uM_DPLD$Feret,t_Protein_12.5uM_DPLD$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("DPLD_5uM", "DPLD_25uM", wilcox.test(t_Protein_5uM_DPLD$Feret,t_Protein_25uM_DPLD$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("DPLD_5uM", "DPLD_50uM", wilcox.test(t_Protein_5uM_DPLD$Feret,t_Protein_50uM_DPLD$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("DPLD_5uM", "DPLD_100uM", wilcox.test(t_Protein_5uM_DPLD$Feret,t_Protein_100uM_DPLD$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("DPLD_12.5uM", "DPLD_25uM", wilcox.test(t_Protein_12.5uM_DPLD$Feret,t_Protein_25uM_DPLD$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("DPLD_12.5uM", "DPLD_50uM", wilcox.test(t_Protein_12.5uM_DPLD$Feret,t_Protein_50uM_DPLD$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("DPLD_12.5uM", "DPLD_100uM", wilcox.test(t_Protein_12.5uM_DPLD$Feret,t_Protein_100uM_DPLD$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("DPLD_25uM", "DPLD_50uM", wilcox.test(t_Protein_25uM_DPLD$Feret,t_Protein_50uM_DPLD$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("DPLD_25uM", "DPLD_100uM", wilcox.test(t_Protein_25uM_DPLD$Feret,t_Protein_100uM_DPLD$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("DPLD_50uM", "DPLD_100uM", wilcox.test(t_Protein_50uM_DPLD$Feret,t_Protein_100uM_DPLD$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)

df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("mNG_5uM", "mNG_12.5uM", wilcox.test(t_Protein_5uM_mNG$Feret,t_Protein_12.5uM_mNG$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("mNG_5uM", "mNG_25uM", wilcox.test(t_Protein_5uM_mNG$Feret,t_Protein_25uM_mNG$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("mNG_5uM", "mNG_50uM", wilcox.test(t_Protein_5uM_mNG$Feret,t_Protein_50uM_mNG$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("mNG_5uM", "mNG_100uM", wilcox.test(t_Protein_5uM_mNG$Feret,t_Protein_100uM_mNG$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("mNG_12.5uM", "mNG_25uM", wilcox.test(t_Protein_12.5uM_mNG$Feret,t_Protein_25uM_mNG$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("mNG_12.5uM", "mNG_50uM", wilcox.test(t_Protein_12.5uM_mNG$Feret,t_Protein_50uM_mNG$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("mNG_12.5uM", "mNG_100uM", wilcox.test(t_Protein_12.5uM_mNG$Feret,t_Protein_100uM_mNG$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("mNG_25uM", "mNG_50uM", wilcox.test(t_Protein_25uM_mNG$Feret,t_Protein_50uM_mNG$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("mNG_25uM", "mNG_100uM", wilcox.test(t_Protein_25uM_mNG$Feret,t_Protein_100uM_mNG$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_Protein[nrow(df_pval_Protein) + 1,] <- c("mNG_50uM", "mNG_100uM", wilcox.test(t_Protein_50uM_mNG$Feret,t_Protein_100uM_mNG$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)

df_pval_Protein$significance <- as.numeric(df_pval_Protein$significance)
df_pval_Protein <- df_pval_Protein %>% mutate(stars = case_when(significance == 0 ~ "NA",
                                                               significance < 0.0001 ~ "****",
                                                               significance < 0.001 ~ "***",
                                                               significance < 0.01 ~ "**",
                                                               significance < 0.05 ~ "*"))

write.csv2(data_Protein, "data_maturation_protein.csv", row.names=FALSE)
write.csv2(means_Protein, "data_maturation_protein_desciptives.csv", row.names=FALSE)
write.csv2(conditions_Protein, "data_maturation_protein_n.csv", row.names=FALSE)
write.csv(df_pval_Protein, "data_maturation_protein_significances.csv", row.names=FALSE)

### NaCl concentration

data_NaCl <- subset(aggregates, Protein_conc == "50" & Temperature == "RT")

conditions_NaCl <- data_NaCl %>% group_by(Protein, Protein_conc, NaCl_conc, Temperature) %>% summarise(Amount=n())

means_NaCl <- aggregate(data_NaCl$Feret, FUN = mean, by = list(data_NaCl$NaCl_conc, data_NaCl$Protein))
names(means_NaCl) <- c("NaCl_conc","Protein","mean")
means_NaCl$sd <- aggregate(data_NaCl$Feret, FUN = sd, by = list(data_NaCl$NaCl_conc, data_NaCl$Protein))[,3]
means_NaCl$median <- aggregate(data_NaCl$Feret, FUN = median, by = list(data_NaCl$NaCl_conc, data_NaCl$Protein))[,3]

t_NaCl_50mM_FL <- subset(data_NaCl, Protein == "FL" & NaCl_conc == "50")
t_NaCl_100mM_FL <- subset(data_NaCl, Protein == "FL" & NaCl_conc == "100")
t_NaCl_150mM_FL <- subset(data_NaCl, Protein == "FL" & NaCl_conc == "150")
t_NaCl_250mM_FL <- subset(data_NaCl, Protein == "FL" & NaCl_conc == "250")

t_NaCl_50mM_PLD <- subset(data_NaCl, Protein == "PLD" & NaCl_conc == "50")
t_NaCl_100mM_PLD <- subset(data_NaCl, Protein == "PLD" & NaCl_conc == "100")
t_NaCl_150mM_PLD <- subset(data_NaCl, Protein == "PLD" & NaCl_conc == "150")
t_NaCl_250mM_PLD <- subset(data_NaCl, Protein == "PLD" & NaCl_conc == "250")

t_NaCl_50mM_DPLD <- subset(data_NaCl, Protein == "DPLD" & NaCl_conc == "50")
t_NaCl_100mM_DPLD <- subset(data_NaCl, Protein == "DPLD" & NaCl_conc == "100")
t_NaCl_150mM_DPLD <- subset(data_NaCl, Protein == "DPLD" & NaCl_conc == "150")
t_NaCl_250mM_DPLD <- subset(data_NaCl, Protein == "DPLD" & NaCl_conc == "250")

df_pval_NaCl <- data.frame(matrix(ncol = 3, nrow = 0)) # starting with an empty data frame with 3 columns
colnames(df_pval_NaCl) <-c("condition 1", "condition 2", "significance")

df_pval_NaCl[nrow(df_pval_NaCl) + 1,] <- c("FL_50mM_NaCl", "FL_100mM_NaCl", wilcox.test(t_NaCl_50mM_FL$Feret,t_NaCl_100mM_FL$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_NaCl[nrow(df_pval_NaCl) + 1,] <- c("FL_50mM_NaCl", "FL_150mM_NaCl", wilcox.test(t_NaCl_50mM_FL$Feret,t_NaCl_150mM_FL$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_NaCl[nrow(df_pval_NaCl) + 1,] <- c("FL_50mM_NaCl", "FL_250mM_NaCl", wilcox.test(t_NaCl_50mM_FL$Feret,t_NaCl_250mM_FL$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_NaCl[nrow(df_pval_NaCl) + 1,] <- c("FL_100mM_NaCl", "FL_150mM_NaCl", wilcox.test(t_NaCl_100mM_FL$Feret,t_NaCl_150mM_FL$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_NaCl[nrow(df_pval_NaCl) + 1,] <- c("FL_100mM_NaCl", "FL_250mM_NaCl", wilcox.test(t_NaCl_100mM_FL$Feret,t_NaCl_250mM_FL$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_NaCl[nrow(df_pval_NaCl) + 1,] <- c("FL_150mM_NaCl", "FL_250mM_NaCl", wilcox.test(t_NaCl_150mM_FL$Feret,t_NaCl_250mM_FL$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)

df_pval_NaCl[nrow(df_pval_NaCl) + 1,] <- c("PLD_50mM_NaCl", "PLD_100mM_NaCl", wilcox.test(t_NaCl_50mM_PLD$Feret,t_NaCl_100mM_PLD$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_NaCl[nrow(df_pval_NaCl) + 1,] <- c("PLD_50mM_NaCl", "PLD_150mM_NaCl", wilcox.test(t_NaCl_50mM_PLD$Feret,t_NaCl_150mM_PLD$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_NaCl[nrow(df_pval_NaCl) + 1,] <- c("PLD_50mM_NaCl", "PLD_250mM_NaCl", wilcox.test(t_NaCl_50mM_PLD$Feret,t_NaCl_250mM_PLD$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_NaCl[nrow(df_pval_NaCl) + 1,] <- c("PLD_100mM_NaCl", "PLD_150mM_NaCl", wilcox.test(t_NaCl_100mM_PLD$Feret,t_NaCl_150mM_PLD$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_NaCl[nrow(df_pval_NaCl) + 1,] <- c("PLD_100mM_NaCl", "PLD_250mM_NaCl", wilcox.test(t_NaCl_100mM_PLD$Feret,t_NaCl_250mM_PLD$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_NaCl[nrow(df_pval_NaCl) + 1,] <- c("PLD_150mM_NaCl", "PLD_250mM_NaCl", wilcox.test(t_NaCl_150mM_PLD$Feret,t_NaCl_250mM_PLD$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)

df_pval_NaCl[nrow(df_pval_NaCl) + 1,] <- c("DPLD_50mM_NaCl", "DPLD_100mM_NaCl", wilcox.test(t_NaCl_50mM_DPLD$Feret,t_NaCl_100mM_DPLD$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_NaCl[nrow(df_pval_NaCl) + 1,] <- c("DPLD_50mM_NaCl", "DPLD_150mM_NaCl", wilcox.test(t_NaCl_50mM_DPLD$Feret,t_NaCl_150mM_DPLD$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_NaCl[nrow(df_pval_NaCl) + 1,] <- c("DPLD_50mM_NaCl", "DPLD_250mM_NaCl", wilcox.test(t_NaCl_50mM_DPLD$Feret,t_NaCl_250mM_DPLD$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_NaCl[nrow(df_pval_NaCl) + 1,] <- c("DPLD_100mM_NaCl", "DPLD_150mM_NaCl", wilcox.test(t_NaCl_100mM_DPLD$Feret,t_NaCl_150mM_DPLD$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_NaCl[nrow(df_pval_NaCl) + 1,] <- c("DPLD_100mM_NaCl", "DPLD_250mM_NaCl", wilcox.test(t_NaCl_100mM_DPLD$Feret,t_NaCl_250mM_DPLD$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_NaCl[nrow(df_pval_NaCl) + 1,] <- c("DPLD_150mM_NaCl", "DPLD_250mM_NaCl", wilcox.test(t_NaCl_150mM_DPLD$Feret,t_NaCl_250mM_DPLD$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)

df_pval_NaCl$significance <- as.numeric(df_pval_NaCl$significance)
df_pval_NaCl <- df_pval_NaCl %>% mutate(stars = case_when(significance == 0 ~ "NA",
                                                          significance < 0.0001 ~ "****",
                                                          significance < 0.001 ~ "***",
                                                          significance < 0.01 ~ "**",
                                                          significance < 0.05 ~ "*"))

ggplot(data = data_NaCl, aes(x = NaCl_conc, y = Feret))+
   geom_violin(trim=TRUE)+
   facet_wrap(~Protein)+
   scale_y_continuous(name = "droplet diameter [µm]", breaks = seq(0,12,0.4))+
   coord_cartesian(ylim = c(0.54,11.46))+
   scale_x_discrete(name= "NaCl concentration [mM]")+
   stat_summary(fun = median, fun.min = median, fun.max = median,
                geom = "crossbar", 
                width = 0.4,
                position = position_dodge(width = 0.4),
                show.legend = FALSE)+
   theme(panel.background = element_rect(fill = NA, color = "grey"),
         panel.grid.major.x = element_blank(),
         panel.grid.major.y = element_line(color = c("grey", NA, NA, NA, NA)),
         panel.grid.minor = element_blank(),
         axis.title = element_text(size = 18),
         axis.text = element_text(size = 13),
         axis.text.y = element_text(color = c("#8c8c8c", NA, NA, NA, NA)),
         axis.ticks.y =element_line(color = "grey"),
         plot.title = element_text(size=20),
         legend.text = element_text(size=15),
         legend.title = element_text(size=18),
         strip.text = element_text(size=16),
         legend.position = "bottom",)

ggplot(data = data_NaCl[data_NaCl$Protein=="FL",], aes(x = NaCl_conc, y = Feret))+
   geom_violin(trim=TRUE)+
   scale_y_continuous(name = "droplet diameter [µm]", breaks = seq(0,12,0.4))+
   coord_cartesian(ylim = c(0.54,11.46))+
   scale_x_discrete(name= "NaCl concentration [mM]")+
   stat_summary(fun = median, fun.min = median, fun.max = median,
                geom = "crossbar", 
                width = 0.4,
                position = position_dodge(width = 0.4),
                show.legend = FALSE)+
   theme(panel.background = element_rect(fill = NA, color = "grey"),
         panel.grid.major.x = element_blank(),
         panel.grid.major.y = element_line(color = c("grey", NA, NA, NA, NA)),
         panel.grid.minor = element_blank(),
         axis.title = element_text(size = 18),
         axis.text = element_text(size = 13),
         axis.text.y = element_text(color = c("#8c8c8c", NA, NA, NA, NA)),
         axis.ticks.y =element_line(color = "grey"),
         plot.title = element_text(size=20),
         legend.text = element_text(size=15),
         legend.title = element_text(size=18),
         strip.text = element_text(size=16),
         legend.position = "bottom",)

write.csv2(data_NaCl, "data_maturation_salt.csv", row.names=FALSE)
write.csv2(means_NaCl, "data_maturation_salt_desciptives.csv", row.names=FALSE)
write.csv2(conditions_NaCl, "data_maturation_salt_n.csv", row.names=FALSE)
write.csv(df_pval_NaCl, "data_maturation_salt_significances.csv", row.names=FALSE)


### RNA concentration

data_RNA <- aggregates[!is.na(aggregates[,"RNA_conc"]),]

conditions_RNA <- data_RNA %>% group_by(Protein, RNA_conc) %>% summarise(Amount=n())

means_RNA <- aggregate(data_RNA$Feret, FUN = mean, by = list(data_RNA$Protein, data_RNA$RNA_conc))
names(means_RNA) <- c("Protein","RNA_conc","mean")
means_RNA$sd <- aggregate(data_RNA$Feret, FUN = sd, by = list(data_RNA$Protein, data_RNA$RNA_conc))[,3]
means_RNA$median <- aggregate(data_RNA$Feret, FUN = median, by = list(data_RNA$Protein, data_RNA$RNA_conc))[,3]


ggplot(data = data_RNA, aes(x = RNA_conc, y = Feret))+
   geom_violin(trim = TRUE)+
#   facet_wrap(~Protein)+
   scale_y_continuous(name = "droplet diameter [µm]", breaks = seq(0,20,0.2))+
   xlab("RNA concentration [ng]")+
   coord_cartesian(ylim = c(0.75,15.75))+
   stat_summary(fun = median, fun.min = median, fun.max = median,
                geom = "crossbar", 
                width = 0.4,
                position = position_dodge(width = 0.4),
                show.legend = FALSE)+
   theme(panel.background = element_rect(fill = NA, color = "grey"),
         panel.grid.major.x = element_blank(),
         panel.grid.major.y = element_line(color = c("grey", NA, NA, NA, NA)),
         panel.grid.minor = element_blank(),
         axis.title = element_text(size = 18),
         axis.title.x = element_text(),
         axis.text = element_text(size = 15),
         axis.text.y = element_text(color = c("#8c8c8c", NA, NA, NA, NA)),
         axis.ticks.y =element_line(color = "grey"),
         plot.title = element_text(size=20),
         legend.text = element_text(size=15),
         legend.title = element_text(size=18),
         strip.text = element_text(size=18),
         legend.position = "bottom",)


ggplot(data = data_RNA[!(data_RNA$RNA_conc=="200" | data_RNA$RNA_conc=="100"),], aes(x = RNA_conc, y = Feret))+
   geom_violin(trim = TRUE)+
   #   facet_wrap(~Protein)+
   scale_y_continuous(name = "droplet diameter [µm]", breaks = seq(0,20,0.2))+
   xlab("RNA concentration [ng]")+
   coord_cartesian(ylim = c(0.55,11.65))+
   stat_summary(fun = median, fun.min = median, fun.max = median,
                geom = "crossbar", 
                width = 0.4,
                position = position_dodge(width = 0.4),
                show.legend = FALSE)+
   theme(panel.background = element_rect(fill = NA, color = "grey"),
         panel.grid.major.x = element_blank(),
         panel.grid.major.y = element_line(color = c("grey", NA, NA, NA, NA)),
         panel.grid.minor = element_blank(),
         axis.title = element_text(size = 18),
         axis.title.x = element_text(),
         axis.text = element_text(size = 15),
         axis.text.y = element_text(color = c("#8c8c8c", NA, NA, NA, NA)),
         axis.ticks.y =element_line(color = "grey"),
         plot.title = element_text(size=20),
         legend.text = element_text(size=15),
         legend.title = element_text(size=18),
         strip.text = element_text(size=18),
         legend.position = "bottom",)

### significance test
t_RNA_FL_0 <- subset(data_RNA, Protein == "FL" & RNA_conc == "0")
t_RNA_FL_1 <- subset(data_RNA, Protein == "FL" & RNA_conc == "1")
t_RNA_FL_5 <- subset(data_RNA, Protein == "FL" & RNA_conc == "5")
t_RNA_FL_12.5 <- subset(data_RNA, Protein == "FL" & RNA_conc == "12.5")
t_RNA_FL_25 <- subset(data_RNA, Protein == "FL" & RNA_conc == "25")
t_RNA_FL_50 <- subset(data_RNA, Protein == "FL" & RNA_conc == "50")
t_RNA_FL_100 <- subset(data_RNA, Protein == "FL" & RNA_conc == "100")
t_RNA_FL_200 <- subset(data_RNA, Protein == "FL" & RNA_conc == "200")

df_pval_RNA <- data.frame(matrix(ncol = 3, nrow = 0)) # starting with an empty data frame with 3 columns
colnames(df_pval_RNA) <-c("condition 1", "condition 2", "significance")

df_pval_RNA[nrow(df_pval_RNA) + 1,] <- c("FL_0ng_RNA", "FL_1ng_RNA", wilcox.test(t_RNA_FL_0$Feret,t_RNA_FL_1$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_RNA[nrow(df_pval_RNA) + 1,] <- c("FL_0ng_RNA", "FL_5ng_RNA", wilcox.test(t_RNA_FL_0$Feret,t_RNA_FL_5$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_RNA[nrow(df_pval_RNA) + 1,] <- c("FL_0ng_RNA", "FL_12.5ng_RNA", wilcox.test(t_RNA_FL_0$Feret,t_RNA_FL_12.5$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_RNA[nrow(df_pval_RNA) + 1,] <- c("FL_0ng_RNA", "FL_25ng_RNA", wilcox.test(t_RNA_FL_0$Feret,t_RNA_FL_25$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_RNA[nrow(df_pval_RNA) + 1,] <- c("FL_0ng_RNA", "FL_50ng_RNA", wilcox.test(t_RNA_FL_0$Feret,t_RNA_FL_50$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_RNA[nrow(df_pval_RNA) + 1,] <- c("FL_0ng_RNA", "FL_100ng_RNA", wilcox.test(t_RNA_FL_0$Feret,t_RNA_FL_100$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_RNA[nrow(df_pval_RNA) + 1,] <- c("FL_0ng_RNA", "FL_200ng_RNA", wilcox.test(t_RNA_FL_0$Feret,t_RNA_FL_200$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)

df_pval_RNA[nrow(df_pval_RNA) + 1,] <- c("FL_1ng_RNA", "FL_5ng_RNA", wilcox.test(t_RNA_FL_1$Feret,t_RNA_FL_5$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_RNA[nrow(df_pval_RNA) + 1,] <- c("FL_1ng_RNA", "FL_12.5ng_RNA", wilcox.test(t_RNA_FL_1$Feret,t_RNA_FL_12.5$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_RNA[nrow(df_pval_RNA) + 1,] <- c("FL_1ng_RNA", "FL_25ng_RNA", wilcox.test(t_RNA_FL_1$Feret,t_RNA_FL_25$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_RNA[nrow(df_pval_RNA) + 1,] <- c("FL_1ng_RNA", "FL_50ng_RNA", wilcox.test(t_RNA_FL_1$Feret,t_RNA_FL_50$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_RNA[nrow(df_pval_RNA) + 1,] <- c("FL_1ng_RNA", "FL_100ng_RNA", wilcox.test(t_RNA_FL_1$Feret,t_RNA_FL_100$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_RNA[nrow(df_pval_RNA) + 1,] <- c("FL_1ng_RNA", "FL_200ng_RNA", wilcox.test(t_RNA_FL_1$Feret,t_RNA_FL_200$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)

df_pval_RNA[nrow(df_pval_RNA) + 1,] <- c("FL_5ng_RNA", "FL_12.5ng_RNA", wilcox.test(t_RNA_FL_5$Feret,t_RNA_FL_12.5$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_RNA[nrow(df_pval_RNA) + 1,] <- c("FL_5ng_RNA", "FL_25ng_RNA", wilcox.test(t_RNA_FL_5$Feret,t_RNA_FL_25$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_RNA[nrow(df_pval_RNA) + 1,] <- c("FL_5ng_RNA", "FL_50ng_RNA", wilcox.test(t_RNA_FL_5$Feret,t_RNA_FL_50$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_RNA[nrow(df_pval_RNA) + 1,] <- c("FL_5ng_RNA", "FL_100ng_RNA", wilcox.test(t_RNA_FL_5$Feret,t_RNA_FL_100$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_RNA[nrow(df_pval_RNA) + 1,] <- c("FL_5ng_RNA", "FL_200ng_RNA", wilcox.test(t_RNA_FL_5$Feret,t_RNA_FL_200$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)

df_pval_RNA[nrow(df_pval_RNA) + 1,] <- c("FL_12.5ng_RNA", "FL_25ng_RNA", wilcox.test(t_RNA_FL_12.5$Feret,t_RNA_FL_25$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_RNA[nrow(df_pval_RNA) + 1,] <- c("FL_12.5ng_RNA", "FL_50ng_RNA", wilcox.test(t_RNA_FL_12.5$Feret,t_RNA_FL_50$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_RNA[nrow(df_pval_RNA) + 1,] <- c("FL_12.5ng_RNA", "FL_100ng_RNA", wilcox.test(t_RNA_FL_12.5$Feret,t_RNA_FL_100$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_RNA[nrow(df_pval_RNA) + 1,] <- c("FL_12.5ng_RNA", "FL_200ng_RNA", wilcox.test(t_RNA_FL_12.5$Feret,t_RNA_FL_200$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)

df_pval_RNA[nrow(df_pval_RNA) + 1,] <- c("FL_25ng_RNA", "FL_50ng_RNA", wilcox.test(t_RNA_FL_25$Feret,t_RNA_FL_50$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_RNA[nrow(df_pval_RNA) + 1,] <- c("FL_25ng_RNA", "FL_100ng_RNA", wilcox.test(t_RNA_FL_25$Feret,t_RNA_FL_100$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_RNA[nrow(df_pval_RNA) + 1,] <- c("FL_25ng_RNA", "FL_200ng_RNA", wilcox.test(t_RNA_FL_25$Feret,t_RNA_FL_200$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)

df_pval_RNA[nrow(df_pval_RNA) + 1,] <- c("FL_50ng_RNA", "FL_100ng_RNA", wilcox.test(t_RNA_FL_50$Feret,t_RNA_FL_100$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)
df_pval_RNA[nrow(df_pval_RNA) + 1,] <- c("FL_50ng_RNA", "FL_200ng_RNA", wilcox.test(t_RNA_FL_50$Feret,t_RNA_FL_200$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)

df_pval_RNA[nrow(df_pval_RNA) + 1,] <- c("FL_100ng_RNA", "FL_200ng_RNA", wilcox.test(t_RNA_FL_100$Feret,t_RNA_FL_200$Feret, exact = FALSE, correct = FALSE, conf.int = FALSE)$p.value)


df_pval_RNA$significance <- as.numeric(df_pval_RNA$significance)
df_pval_RNA <- df_pval_RNA %>% mutate(sign = case_when(significance == 0 ~ "NA",
                                                       significance < 0.0001 ~ "****",
                                                       significance < 0.001 ~ "***",
                                                       significance < 0.01 ~ "**",
                                                       significance < 0.05 ~ "*"))

write.csv2(data_RNA, "data_RNA_concentration.csv", row.names=FALSE)
write.csv2(means_RNA, "data_RNA_concentration_desciptives.csv", row.names=FALSE)
write.csv2(conditions_RNA, "data_RNA_concentration_n.csv", row.names=FALSE)
write.csv(df_pval_RNA, "data_RNA_concentration_significances.csv", row.names=FALSE)
