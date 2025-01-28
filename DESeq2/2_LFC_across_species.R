# QDR RNAseq Solanum species
# LFC distribution of all DEGs
# contrast: all samples Inf vs Mock
# Severin Einspanier 
# 2024_09_13
library(tidyverse)
library(ggtext)
library(ggdist)
library(glue)
library(patchwork)
rm(list=ls())

setwd()

data_LFC_sus <- read.csv("DeSeq_OUT/combined_SUS_inf_mock.csv") %>% 
  mutate(res="SUSCEPTIBLE") %>% 
  select(-X)

data_LFC <- read.csv("DeSeq_OUT/combined_RES_inf_mock.csv") %>% 
  mutate(res="RESISTANT") %>% 
  select(-X) %>% 
  bind_rows(data_LFC_sus) %>% 
  rename(res="Resistance-Level")
  



bg_color <- "grey97"
font_family <- "Fira Sans"
svg(paste0("C:figures/fig_2/", Sys.Date(),"_LFC_INF_MOCK.svg"), 
    width=12, height=9, bg="transparent")
(p1 <- data_LFC %>% 
    ggplot(aes(x=log2FoldChange,y=species, fill=`Resistance-Level`,col=`Resistance-Level`))+
    #geom_density(alpha=.75, na.rm = T,
    #             adjust=1.5, 
    #             bounds=c(-8,8)
    #             )+
    stat_slab(fill_type = "segments",  
                   normalize="groups", show.legend = F,
                 geom = "slabinterval", 
              alpha=.7) +
    #stat_interval(show.legend = F) +
    stat_summary( geom = "point", 
                 size=3.4, fun = median, show.legend = T)+
    xlim(-5, 5)+
    # adjust legend title here please:
    theme_bw()+
    theme(
      plot.background = element_rect(color = NA, fill = bg_color),
      panel.grid = element_blank(),
      panel.grid.major.x = element_line(linewidth = 0.1, color = "grey75"),
      
      plot.title.position = "plot",
      plot.subtitle = element_textbox_simple(
        margin = margin(t = 4, b = 16), size = 10),
      plot.caption = element_textbox_simple(
        margin = margin(t = 12), size = 7
      ),
      plot.caption.position = "plot",
      #axis.text.y=element_blank(), 
      axis.text.y = element_text(hjust = 0, margin = margin(r = -1), family = "Fira Sans SemiBold", 
                                 colour = "black", size=15, face="italic"),
      axis.text.x = element_text(colour = "black", size=15),
      #plot.margin = margin(4, 4, 4, 4), 
      axis.title.y = element_blank(),
      axis.title.x = element_text(color="black", size=15),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10), # Add space at the top
      strip.text = element_text(size=15, face="italic"),
      legend.position = "top"
    )
)

dev.off()
ggsave("figures/2024_09_13_LFC_INF_MOCK.svg", 
       p1, width=7, height=8, dpi=900)

