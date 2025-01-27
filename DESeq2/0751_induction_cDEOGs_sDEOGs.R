rm(list=ls())
# Check whether cDEOGs are induced when infected

## DEPRECATED!!! look at OGs/cDEOGs_induction.R

library(tidyverse)
# RES

DEGs_res <- read.csv("C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/DeSeq/data/combined_Res_Inf_vs_Mock_results_DEGs.csv")
DEGs_sus <- read.csv("C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/DeSeq/data/combined_Sus_Inf_vs_Mock_results_DEGs.csv")

core_DEOGs <- read.csv("C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/OGs/data/OG_cDEOGs_IDs.csv", check.names = F)

degs_data_wide_res <- DEGs_res %>% 
  pivot_wider(names_from = species, values_from = GeneID) %>% 
  unnest()


induced_cDEOGs <- core_DEOGs %>% 
  dplyr::group_by(species) %>% 
  filter(
    (species == "S. chilense" & GeneID %in% na.omit(degs_data_wide_res$`S. chilense`)) |
      (species == "S. pennellii" & GeneID %in% na.omit(degs_data_wide_res$`S. pennellii`)) |
      (species == "S. habrochaites" & GeneID %in% na.omit(degs_data_wide_res$`S. habrochaites`)) |
      (species == "S. lycopersicoides" & GeneID %in% na.omit(degs_data_wide_res$`S. lycopersicoides`)) |
      (species == "S. pimpinellifolium" & GeneID %in% na.omit(degs_data_wide_res$`S. pimpinellifolium`))
  ) %>%
  ungroup() %>% 
  group_by(OG) %>% 
  summarise(num_species_with_DEG = n_distinct(species))

non_induced <- core_DEOGs %>% 
  anti_join(induced_cDEOGs, by="OG") %>% 
  summarise(OG=OG, 
            num_species_with_DEG=0)

(p1 <- induced_cDEOGs %>% 
    bind_rows(non_induced) %>% 
    ggplot(aes(x=as.factor(num_species_with_DEG)))+
    geom_bar()+
    ggtitle("DEGs in res. INF-Mock (also cDEOG)")+
    labs(x="Number of species with DEG \n Inf. vs. Mock", y="Number of DEOGs")+
    theme_bw()+
    theme(axis.text = element_text(size=11, color="Black"),
          #axis.text.x = element_text(size=11, angle=45,
          #                           hjust=1, vjust=1, color="black"),
          axis.title = element_text(size=11, color="Black"),
          legend.title = element_text(size=13, margin = margin(t = 0, r = 10, b = 0, l = 0)),
          legend.text = element_text(size=11), 
          legend.spacing.x = unit(.5, 'cm'),
          panel.grid.major.x =  element_blank(),
          panel.grid.major.y = element_line(size=.4, color="grey50", linetype = "dashed"),
          panel.grid.minor.y = element_blank(),
          #axis.text.x = element_text(hjust = 0.5, face="italic", color="black"),
          axis.title.x = element_blank()
    ))


# sus


degs_data_wide_sus <- DEGs_sus %>% 
  pivot_wider(names_from = species, values_from = GeneID) %>% 
  unnest()


induced_cDEOGs_sus <- core_DEOGs %>% 
  dplyr::group_by(species) %>% 
  filter(
    (species == "S. chilense" & GeneID %in% na.omit(degs_data_wide_sus$`S. chilense`)) |
      (species == "S. pennellii" & GeneID %in% na.omit(degs_data_wide_sus$`S. pennellii`)) |
      (species == "S. habrochaites" & GeneID %in% na.omit(degs_data_wide_sus$`S. habrochaites`)) |
      (species == "S. lycopersicoides" & GeneID %in% na.omit(degs_data_wide_sus$`S. lycopersicoides`)) |
      (species == "S. pimpinellifolium" & GeneID %in% na.omit(degs_data_wide_sus$`S. pimpinellifolium`))
  ) %>%
  ungroup() %>% 
  group_by(OG) %>% 
  summarise(num_species_with_DEG = n_distinct(species))

non_induced_sus <- core_DEOGs %>% 
  anti_join(induced_cDEOGs_sus, by="OG") %>% 
  summarise(OG=OG, 
            num_species_with_DEG=0)

(p2 <- induced_cDEOGs_sus %>% 
    bind_rows(non_induced_sus) %>% 
    ggplot(aes(x=as.factor(num_species_with_DEG)))+
    geom_bar()+
    ggtitle("DEGs in SUS INF-Mock (also cDEOG)")+
    labs(x="Number of species with DEG \n Inf. vs. Mock", y="Number of DEOGs")+
    theme_bw()+
    theme(axis.text = element_text(size=11, color="Black"),
          #axis.text.x = element_text(size=11, angle=45,
          #                           hjust=1, vjust=1, color="black"),
          axis.title = element_text(size=11, color="Black"),
          legend.title = element_text(size=13, margin = margin(t = 0, r = 10, b = 0, l = 0)),
          legend.text = element_text(size=11), 
          legend.spacing.x = unit(.5, 'cm'),
          panel.grid.major.x =  element_blank(),
          panel.grid.major.y = element_line(size=.4, color="grey50", linetype = "dashed"),
          panel.grid.minor.y = element_blank(),
          #axis.text.x = element_text(hjust = 0.5, face="italic", color="black"),
          axis.title.x = element_blank()
    ))


#
p_tot <- ggpubr::ggarrange(p1,p2, nrow=1)

ggsave("C:/Users/suaph281/Nextcloud/ResiDEvo/sequencing/figures/fig_3/induction_cDEOGs_sDEOGs.svg", p_tot, width=10, height=5)

# Nicer plot 
cDEOGs_res <- induced_cDEOGs %>% 
  bind_rows(non_induced) %>% 
  mutate(resistance="resistant")

complete_data <- induced_cDEOGs_sus %>% 
  bind_rows(non_induced_sus) %>% 
  mutate(resistance="susceptible") %>% 
  bind_rows(cDEOGs_res)

(p_3 <- complete_data %>% 
    group_by(resistance) %>%
    summarise(count = n()) %>%
    mutate(percentage = count / sum(count) * 100) 
  )

%>%
  ggplot(aes(x=as.factor(num_species_with_DEG), 
             fill=resistance, y=percentage))+
  geom_col(position="dodge", stat="count")


### Species specific OG induction 



# Define the function to analyze and plot data for a given species
analyze_species <- function(species_name, DEOGs, orthogroups_genid_edited_long, degs_data_wide_sus, degs_data_wide_res) {
  DEOGs_species <- DEOGs %>% 
    dplyr::filter(!!sym(species_name) == 1 & sum == 1) %>%
    select(Orthogroup) %>% 
    left_join(orthogroups_genid_edited_long, by = "Orthogroup") %>% 
    filter(species == species_name) %>% 
    mutate(induced_sus = ifelse(genelist %in% na.omit(degs_data_wide_sus[[species_name]]), 1, 0),
           induced_res = ifelse(genelist %in% na.omit(degs_data_wide_res[[species_name]]), 1, 0)) %>% 
    mutate(ratio_sus=sum(induced_sus)/n(),
           ratio_res=sum(induced_res)/n()) %>%
    pivot_longer(cols = c(ratio_sus, ratio_res), names_to = "RESISTANCE", values_to = "ratio") %>% 
    select(species, RESISTANCE, ratio) %>% 
    unique()
  
  # Plot the results
  ggplot(DEOGs_species, aes(x = as.factor(RESISTANCE), y=ratio*100, fill=as.factor(RESISTANCE))) +
    geom_col(col="black", size=1, show.legend = F) +
    labs(title = paste("Induction of sDEOGs", species_name),
         x = "DEG in Mock-Infected",
         y = "Induced sDEOGs [%]") +
    theme_bw()+
    ylim(0,100)+
    scale_fill_manual(values = c("#044993", "#FFB632"))+
    theme(axis.text = element_text(size=11, color="Black"),
          #axis.text.x = element_text(size=11, angle=45,
          #                           hjust=1, vjust=1, color="black"),
          axis.title = element_text(size=11, color="Black"),
          legend.title = element_text(size=13, margin = margin(t = 0, r = 10, b = 0, l = 0)),
          legend.text = element_text(size=11), 
          legend.spacing.x = unit(.5, 'cm'),
          panel.grid.major.x =  element_blank(),
          panel.grid.major.y = element_line(size=.4, color="grey50", linetype = "dashed"),
          panel.grid.minor.y = element_blank(),
          #axis.text.x = element_text(hjust = 0.5, face="italic", color="black"),
          axis.title.x = element_blank(),
          strip.text = element_text(size=13, color="black")
    )
}

# Example usage for multiple species
species_list <- c("S. pennellii", "S. lycopersicoides", "S. habrochaites", "S. pimpinellifolium", "S. chilense")  # Add more species as needed

plot_list <- list()
orthogroups_genid_edited_long <- read.csv("C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/DeSeq/data/2024_09_09_ORTHOGORUPS_name_edited.csv", stringsAsFactors = T, header=T, check.names = F)
DEOGs <- read.csv("C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/DeSeq/data/2024_09_09_DEOGs.csv", stringsAsFactors = T, header=T, check.names = F)
colnames(DEOGs)[1]="wayne"

# Loop through each species and apply the function, storing the plots in the list
for (species in species_list) {
  plot_list[[species]] <- analyze_species(species, DEOGs, orthogroups_genid_edited_long, degs_data_wide_sus, degs_data_wide_res)
}

# Combine all plots into one central figure using patchwork
combined_plot <- patchwork::wrap_plots(plot_list)

#
combined_plot
# Save the combined plot as a single figure
ggsave(paste0("C:/Users/suaph281/Nextcloud/ResiDEvo/sequencing/figures/fig_3/", Sys.Date(), "induction_species.svg"), combined_plot, width = 20, height = 10)



