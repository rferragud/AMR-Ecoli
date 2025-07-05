
###################################################################################################################################################
###                                                       LOADING AND PREPARING METADATA                                                       ###                     
###################################################################################################################################################
library(dplyr)
library(svglite)

file = "C:/Users/rferragud/Documents/Projecte/Tables/amrf_gen_fen_ns4.csv"
table = read.delim(file, sep = ",", header = T) %>% filter(QC == "included")
dim(table)
# [1] 169001     13

str(table)

valores_excluir <- c(0.047, 0.094, 0.19, 0.38, 0.75, 1.5, 3)

# Filtrar el data frame quitando esas filas
table <- table %>%
  filter(!resultados_mic %in% valores_excluir)

###################################################################################################################################################
###                                                             TRANSFORMING MICS                                                               ###                     
###################################################################################################################################################

# To change for list of antibiotics of interest

selected_antibiotics = c("meropenem", "gentamicin", "ampicillin", "cefotaxime")

# ciprofloxacin ECOFF = 0.06 
# meropenem ECOFF = 0.06 
# gentamicin ECOFF = 2

# These vectors are used to specify the floor (<=) and ceiling (>=) MIC values to plot for

MIC_range_from = c(0.06, 0.5, 2, 0.03)

MIC_range_to = c(16, 8, 32, 8)

# NOTE: the low range is decided as the highest ≤ value. E.g. for vancomyin, ≤1 is chosen (MIC values: <2  >16 >256  >32 ≤0.5   ≤1 ≥256  ≥32  0.5    1   16    2  256    4   64    8)
# NOTE: the high range is decided as the lowest ≥ value or the lowest > value multiplied by 2. E.g. for vancomyin, ≥32 is chosen (MIC values: <2  >16 >256  >32 ≤0.5   ≤1 ≥256  ≥32  0.5    1   16    2  256    4   64    8)

# To check with Roser

#aaa = which(table$resultados_antibiotic=="ciprofloxacin")
aaa = which(table$resultados_antibiotic=="meropenem")
aaa = which(table$resultados_antibiotic=="gentamicin")
aaa = which(table$resultados_antibiotic=="ampicillin")
aaa = which(table$resultados_antibiotic=="cefotaxime")

length(aaa)

table(table$resultados_mic[aaa], table$resultados_ECOFF[aaa])


# Transform function
transform_MIC = function(x)
{
  z = x
  if(grepl("<",x)){ y = gsub("<", "", x);  z = as.numeric(y)/2; }
  if(grepl(">",x)){ y = gsub(">", "", x);  z = as.numeric(y)*2; }
  if(grepl("<=",x)){ y = gsub("<=", "", x);  z = as.numeric(y); }
  if(grepl(">=",x)){ y = gsub(">=", "", x);  z = as.numeric(y); }
  z = as.numeric(z)
  return(z)
}

extract_MIC_value = function(x)
{
  y = x
  if(grepl("<",x)){ y = gsub("<", "", x); }
  if(grepl(">",x)){ y = gsub(">", "", x); }
  if(grepl("<=",x)){ y = gsub("<=", "", x); }
  if(grepl(">=",x)){ y = gsub(">=", "", x); }
  y = as.numeric(y)
  return(y)
}

# New table with modified MICs
# New variable used to store MIC to be plotted (resultados_mic_modified)
table$resultados_mic_modified = NA

for(a in 1:length(selected_antibiotics))
{
  print(a)
  antibiotic = selected_antibiotics[a]
  print(antibiotic)
  aaa = which(table$resultados_antibiotic==antibiotic); # table row indices corresponding to antibiotic
  mics = table$resultados_mic[aaa]
  # MIC are transformed
  mics_mod = sapply(mics, transform_MIC)
  mics_values = sapply(mics, extract_MIC_value)
  mics_mod[which(mics_mod<=MIC_range_from[a])] = MIC_range_from[a]
  mics_mod[which(mics_mod>=MIC_range_to[a])] = MIC_range_to[a]
  # dealing with <MICs greater (e.g. <0.25) than MIC_range_from (<=0.015)
  tmp = which(grepl("<",mics) & mics_values>MIC_range_from[a])
  if(length(tmp)>0){ mics_mod[tmp] = NA; }
  # dealing with >MICs lower (e.g. >1) than MIC_range_to (>=2)
  tmp = which(grepl(">",mics) & mics_values<MIC_range_to[a])
  if(length(tmp)>0){ mics_mod[tmp] = NA; }
  # revming odd MIC > 
  
  # saving results
  table$resultados_mic_modified[aaa] = mics_mod
  print(table(mics))
  print(table(mics_mod))
  print(table(mics, mics_mod))
}




###################################################################################################################################################
###                                                PLOTTING MICS DISTRIBUTION ALONG WITH GENETIC DETERMINANTS                                  ###                     
###################################################################################################################################################

### Plotting MIC distributions with genotypes

require(ggplot2)
library(RColorBrewer)

getPalette = colorRampPalette(brewer.pal(9, "Paired"))

## plotting variables
size_dot = 1;
size_axis_lines = 0.3;
text_y_offset = 4;
font = "Times";
dot_color = "dimgray";
axis_text_size = 10;
axis_title_size = 20;
ann_text_size = 5;

for(a in 1:length(selected_antibiotics))
{
  antibiotic = selected_antibiotics[a]
  aaa = which(table$resultados_antibiotic==antibiotic); # table row indices corresponding to antibiotic
  table_ab = table[aaa,c("resultados_mic_modified","amrf")]
  colnames(table_ab) = c("MIC","genotype")
  tmp = which(is.na(table_ab[,"MIC"]))
  if(length(tmp)>0){ table_ab = table_ab[-tmp,]; }
  tmp = which(is.na(table_ab[,"genotype"]))
  if(length(tmp)>0){ table_ab[tmp,"genotype"] = "00-No determinant"; }
  table_ab = as.data.frame(table_ab)
  table_ab$MIC_tmp = as.numeric(as.character(table_ab$MIC))
  ## Plotting MIC with genotypes
  gp <- ggplot(table_ab, aes(x = reorder(MIC, MIC_tmp), fill = genotype))
    # gp = gp + scale_fill_brewer(palette = "Paired")
  gp = gp + scale_fill_manual(values = getPalette(length(unique(table_ab$genotype))))
  gp = gp + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
                  axis.line = element_line(colour = "black", size = size_axis_lines), axis.ticks = element_line(size = size_axis_lines))
  gp = gp + geom_bar()
  gp = gp + xlab("MIC") + ylab("Number of Isolates") + ggtitle("MIC distribution") + theme(text = element_text(family = font)) + theme(axis.text = element_text(size=axis_text_size, color="black"), axis.title=element_text(size=axis_title_size), title=element_text(size=axis_title_size))
  plot_file = paste(antibiotic, "_MIC_distribution.png", sep = "")
  ggsave(plot_file, plot = gp, device = "png",width = 10, height = 5, dpi = 300, units = "in")
  plot_file = paste(antibiotic, "_MIC_distribution.no_legend.png", sep = "")
  gp = gp + theme(legend.position = "none")
  # command below only used for a few antibiotics
  # gp = gp + scale_y_continuous(name="Number of Isolates", breaks=seq(0,3000,100), labels=seq(0,3000,100))
  ggsave(plot_file, plot = gp, device = "png",width = 5, height = 5, dpi = 300, units = "in")
  # three commands below only used to zoom in for a few antibiotics
  # gp = gp + coord_cartesian(ylim=c(0, 500)) # NOTE: using only gp + ylim=c(0, 500), will remove data, see url below
  # https://stackoverflow.com/questions/25685185/limit-ggplot2-axes-without-removing-data-outside-limits-zoom
  # plot_file = paste(antibiotic, "_MIC_distribution.no_legend.y0_500.png", sep = "")
  # ggsave(plot_file, plot = gp, device = "png",width = 5, height = 5, dpi = 300, units = "in")
}


getPalette = colorRampPalette(brewer.pal(9, "Paired"))

## plotting variables
size_dot = 1;
size_axis_lines = 0.3;
text_y_offset = 4;
font = "Times";
dot_color = "dimgray";
axis_text_size = 10;
axis_title_size = 20;
ann_text_size = 5;

# amrf

for(a in 1:length(selected_antibiotics))
{
  antibiotic = selected_antibiotics[a]
  aaa = which(table$resultados_antibiotic==antibiotic)  # filas del antibiótico
  table_ab = table[aaa,c("resultados_mic_modified","amrf")]
  colnames(table_ab) = c("MIC","genotype")
  tmp = which(is.na(table_ab[,"MIC"]))
  if(length(tmp)>0){ table_ab = table_ab[-tmp,]; }
  tmp = which(is.na(table_ab[,"genotype"]))
  if(length(tmp)>0){ table_ab[tmp,"genotype"] = "00-No determinant"; }
  table_ab = as.data.frame(table_ab)
  table_ab$MIC_tmp = as.numeric(as.character(table_ab$MIC))
  
  ## —— extraer valores de corte ECOFF y EUCAST
  eco_val <- unique(table$resultados_S_ECOFF[aaa])
  euc_val <- unique(table$resultados_S_EUCAST[aaa])
  
  ## Plotting MIC with genotypes
  gp <- ggplot(table_ab, aes(x = reorder(MIC, MIC_tmp), fill = genotype))
  gp = gp + scale_fill_manual(values = getPalette(length(unique(table_ab$genotype))))
  gp = gp + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
                  axis.line = element_line(colour = "black", size = size_axis_lines), axis.ticks = element_line(size = size_axis_lines))
  gp = gp + geom_bar()
  
  ## —— calcular posición en el eje x para las líneas
  level_order <- levels(reorder(table_ab$MIC, table_ab$MIC_tmp))
  eco_pos <- match(as.character(eco_val), level_order)
  euc_pos <- match(as.character(euc_val), level_order)
  
  ## —— añadir líneas verticales discontinuas
  gp = gp +
    geom_vline(xintercept = eco_pos, linetype = "dashed", size = 0.5, color = "gray") +
    geom_vline(xintercept = euc_pos, linetype = "dashed", size = 0.5, color = "gray")
  
  gp = gp + xlab("MIC") + ylab("Number of Isolates") + ggtitle("MIC distribution") +
    theme(text = element_text(family = font)) +
    theme(axis.text = element_text(size=axis_text_size, color="black"),
          axis.title=element_text(size=axis_title_size),
          title=element_text(size=axis_title_size))
  
  plot_file = paste(antibiotic, "_MIC_distribution_amrf2.svg", sep = "")
  ggsave(plot_file, plot = gp, width = 10, height = 5, dpi = 300, units = "in")
  
  plot_file = paste(antibiotic, "_MIC_distribution.no_legend_amrf2.svg", sep = "")
  gp = gp + theme(legend.position = "none")
  ggsave(plot_file, plot = gp, device = "svg",width = 5, height = 5, dpi = 300, units = "in")
}


# ResFinder

for(a in 1:length(selected_antibiotics))
{
  antibiotic = "ampicillin"
  aaa = which(table$resultados_antibiotic==antibiotic)  # filas del antibiótico
  table_ab = table[aaa,c("resultados_mic_modified","resf")]
  colnames(table_ab) = c("MIC","genotype")
  tmp = which(is.na(table_ab[,"MIC"]))
  if(length(tmp)>0){ table_ab = table_ab[-tmp,]; }
  tmp = which(is.na(table_ab[,"genotype"]))
  if(length(tmp)>0){ table_ab[tmp,"genotype"] = "00-No determinant"; }
  table_ab = as.data.frame(table_ab)
  table_ab$MIC_tmp = as.numeric(as.character(table_ab$MIC))
  
  ## —— extraer valores de corte ECOFF y EUCAST
  eco_val <- unique(table$resultados_S_ECOFF[aaa])
  euc_val <- unique(table$resultados_S_EUCAST[aaa])
  
  ## Plotting MIC with genotypes
  gp <- ggplot(table_ab, aes(x = reorder(MIC, MIC_tmp), fill = genotype))
  gp = gp + scale_fill_manual(values = getPalette(length(unique(table_ab$genotype))))
  gp = gp + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
                  axis.line = element_line(colour = "black", size = size_axis_lines), axis.ticks = element_line(size = size_axis_lines))
  gp = gp + geom_bar()
  
  ## —— calcular posición en el eje x para las líneas
  level_order <- levels(reorder(table_ab$MIC, table_ab$MIC_tmp))
  eco_pos <- match(as.character(eco_val), level_order)
  euc_pos <- match(as.character(euc_val), level_order)
  
  ## —— añadir líneas verticales discontinuas
  gp = gp +
    geom_vline(xintercept = eco_pos, linetype = "dashed", size = 0.5, color = "gray") +
    geom_vline(xintercept = euc_pos, linetype = "dashed", size = 0.5, color = "gray")
  
  gp = gp + xlab("MIC") + ylab("Number of Isolates") + ggtitle("MIC distribution") +
    theme(text = element_text(family = font)) +
    theme(axis.text = element_text(size=axis_text_size, color="black"),
          axis.title=element_text(size=axis_title_size),
          title=element_text(size=axis_title_size))
  
  plot_file = paste(antibiotic, "_MIC_distribution_resf2.svg", sep = "")
  ggsave(plot_file, plot = gp, device = "svg",width = 10, height = 5, dpi = 300, units = "in")
  
  plot_file = paste(antibiotic, "_MIC_distribution.no_legend_resf2.svg", sep = "")
  gp = gp + theme(legend.position = "none")
  ggsave(plot_file, plot = gp, device = "svg",width = 5, height = 5, dpi = 300, units = "in")
}





for(a in 1:length(selected_antibiotics))
{
  antibiotic = "meropenem"
  aaa = which(table$resultados_antibiotic==antibiotic)  # filas del antibiótico
  table_ab = table[aaa,c("resultados_mic_modified","resf")]
  colnames(table_ab) = c("MIC","genotype")
  tmp = which(is.na(table_ab[,"MIC"]))
  if(length(tmp)>0){ table_ab = table_ab[-tmp,]; }
  tmp = which(is.na(table_ab[,"genotype"]))
  if(length(tmp)>0){ table_ab[tmp,"genotype"] = "00-No determinant"; }
  table_ab = as.data.frame(table_ab)
  table_ab$MIC_tmp = as.numeric(as.character(table_ab$MIC))
  
  ## —— extraer valores de corte ECOFF y EUCAST
  eco_val <- unique(table$resultados_S_ECOFF[aaa])
  euc_val <- unique(table$resultados_S_EUCAST[aaa])
  
  ## Plotting MIC with genotypes
  gp <- ggplot(table_ab, aes(x = reorder(MIC, MIC_tmp), fill = genotype))
  gp = gp + scale_fill_manual(values = getPalette(length(unique(table_ab$genotype))))
  gp = gp + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
                  axis.line = element_line(colour = "black", size = size_axis_lines), axis.ticks = element_line(size = size_axis_lines))
  gp = gp + geom_bar()
  
  ## —— calcular posición en el eje x para las líneas
  level_order <- levels(reorder(table_ab$MIC, table_ab$MIC_tmp))
  eco_pos <- match(as.character(eco_val), level_order)
  euc_pos <- match(as.character(euc_val), level_order)
  
  ## —— añadir líneas verticales discontinuas
  gp = gp +
    geom_vline(xintercept = eco_pos, linetype = "dashed", size = 0.5, color = "gray") +
    geom_vline(xintercept = euc_pos, linetype = "dashed", size = 0.5, color = "gray")
  
  gp = gp + xlab("MIC") + ylab("Number of Isolates") + ggtitle("MIC distribution") +
    theme(text = element_text(family = font)) +
    theme(axis.text = element_text(size=axis_text_size, color="black"),
          axis.title=element_text(size=axis_title_size),
          title=element_text(size=axis_title_size))
  
  plot_file = paste(antibiotic, "_MIC_distribution_resf.png", sep = "")
  ggsave(plot_file, plot = gp, device = "png",width = 10, height = 5, dpi = 300, units = "in")
  
  plot_file = paste(antibiotic, "_MIC_distribution.no_legend_resf.png", sep = "")
  gp = gp + theme(legend.position = "none")# command below only used for a few antibiotics
  # gp = gp + scale_y_continuous(name="Number of Isolates", breaks=seq(0,3000,100), labels=seq(0,3000,100))
  ggsave(plot_file, plot = gp, device = "png",width = 5, height = 5, dpi = 300, units = "in")
  # three commands below only used to zoom in for a few antibiotics
  gp = gp + coord_cartesian(ylim=c(0, 500)) # NOTE: using only gp + ylim=c(0, 500), will remove data, see url below
  # https://stackoverflow.com/questions/25685185/limit-ggplot2-axes-without-removing-data-outside-limits-zoom
  plot_file = paste(antibiotic, "_MIC_distribution.no_legend.y0_500.png", sep = "")
  ggsave(plot_file, plot = gp, device = "png",width = 5, height = 5, dpi = 300, units = "in")
}
