#First Vis for Founderevents

rm(list=ls())
library(tidyverse)

# Verzeichnis, in dem sich die Dateien befinden
directory <- "C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/TAI/data/GenEra/"

# Liste aller Dateien im Verzeichnis (z.B. CSV-Dateien)
file_list <- list.files(path = directory, pattern = "*_founder_summary.tsv", full.names = TRUE)


# Funktion zum Einlesen einer Datei und Entfernen des '#' in der Header-Zeile
read_custom_table <- function(file) {
  # Lese die Datei ein
  df <- read.table(file, sep = "\t", header = TRUE, comment.char = "")
  
  # Entferne den '#' aus den Spaltennamen, falls vorhanden
  names(df) <- sub("^#", "", names(df))
  
  # Extrahiere den Dateinamen ohne Pfad
  filename <- basename(file)
  # Extrahiere die Nummer aus dem Dateinamen (alles vor dem ersten '_')
  file_id <- gsub("_.*", "", filename)
  
  # Füge den Dateinamen als neue Spalte hinzu
  df$taxid <- as.integer(file_id)
  
  
  return(df)
}

combined_data <- do.call(rbind, lapply(file_list, read_custom_table))

# Zeige die ersten Zeilen des kombinierten Dataframes
head(combined_data)

colnames(combined_data)[1]="number_of_genes_families"

ncbi_key <- read.csv("C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/TAI/data/GenEra/ncbi_taxid.csv")

combined_data_edit <- combined_data %>% 
  left_join(ncbi_key, by="taxid") %>% 
  select(!X)



phylorank = data.frame(phylostratum_general=c("Cellular organisms","Eukaryota","Viridiplantae",
                                              "Streptophyta","Streptophytina","Embryophyta",
                                              "Tracheophyta","Euphyllophyta","Spermatophyta",
                                              "Magnoliopsida","Mesangiospermae","Pentapetalae",
                                              "Asterids","Lamiids","Solanales","Solanaceae",
                                              "Solanoideae","Solanum","Lycopersicoides/Lycopersicon",
                                              "Species"),
                       phylorank=seq(1:20)
)

library(ggsci)

cols = pal_simpsons(palette = c("springfield"), alpha = .9)(5)


(p1 <- combined_data_edit %>% 
    left_join(phylorank, by="phylorank") %>%
    arrange(phylorank) %>% 
    mutate(phylostratum_general = fct_reorder(phylostratum_general, phylorank)) %>% 
    ggplot(aes(x=as.factor(phylostratum_general), y=(number_of_genes_families), col=species, group=species))+
    geom_rect(xmin=6 - 0.5, xmax=6 + 0.5, 
              ymin=-Inf, ymax=Inf, fill="grey90", alpha=0.1, inherit.aes = FALSE) +
    geom_rect(xmin=16 - 0.5, xmax=16 + 0.5, 
              ymin=-Inf, ymax=Inf, fill="grey90", alpha=0.1, inherit.aes = FALSE) +
    geom_line(size=2, show.legend = F)+
    geom_point(size=4, alpha=.8)+
    theme_bw()+
    scale_color_manual(values = rev(cols))+
    labs(x="Phylorank",
         y="#Gene of families", 
         color="Species")+
    scale_y_continuous(trans='log10')+
    theme(legend.position = c(.8, .8))+
    theme(axis.text.x = element_text(angle=45, hjust=1),
          axis.text=element_text(color="black", size=13),
          axis.title = element_text(color="black", size=15),
          legend.title = element_text(size=17, face="bold"),
          legend.text = element_text(size=14, face="italic")))

ggsave("C:/Users/suaph281//Nextcloud/ResiDEvo/sequencing/figures/fig_1/2024_09_08_gene_founders.png", 
       p1, dpi=900, width=8, height=7)


svg(paste0("C:/Users/suaph281//Nextcloud/ResiDEvo/sequencing/figures/fig_1/2024_09_08_gene_founders.svg"), 
    width=8, height=7, bg="transparent")

(p1 <- combined_data_edit %>% 
    left_join(phylorank, by="phylorank") %>%
    arrange(phylorank) %>% 
    mutate(phylostratum_general = fct_reorder(phylostratum_general, phylorank)) %>% 
    ggplot(aes(x=as.factor(phylostratum_general), y=(number_of_genes), col=species, group=species))+
    geom_line(size=2, show.legend = F)+
    geom_point(size=4, alpha=.8)+
    theme_bw()+
    scale_color_manual(values = rev(cols))+
    labs(x="Phylorank",
         y="#Genes", 
         color="Species")+
    scale_y_continuous(trans='log10')+
    theme(legend.position = c(.8, .8))+
    theme(axis.text.x = element_text(angle=45, hjust=1),
          axis.text=element_text(color="black", size=13),
          axis.title = element_text(color="black", size=15),
          legend.title = element_text(size=17, face="bold"),
          legend.text = element_text(size=14, face="italic")))

dev.off()


# Try to merge with gene age 

file_list <- list.files(path = directory, pattern = "*gene_age_summary.tsv", full.names = TRUE)


# Funktion zum Einlesen einer Datei und Entfernen des '#' in der Header-Zeile
read_custom_table <- function(file) {
  # Lese die Datei ein
  df <- read.table(file, sep = "\t", header = TRUE, comment.char = "")
  
  # Entferne den '#' aus den Spaltennamen, falls vorhanden
  names(df) <- sub("^#", "", names(df))
  
  # Extrahiere den Dateinamen ohne Pfad
  filename <- basename(file)
  # Extrahiere die Nummer aus dem Dateinamen (alles vor dem ersten '_')
  file_id <- gsub("_.*", "", filename)
  
  # Füge den Dateinamen als neue Spalte hinzu
  df$taxid <- as.integer(file_id)
  
  
  return(df)
}

combined_data_age <- do.call(rbind, lapply(file_list, read_custom_table))

# Zeige die ersten Zeilen des kombinierten Dataframes
head(combined_data_age)

colnames(combined_data_age)[1]="number_of_genes"

#ncbi_key <- read.csv("C:/Users/suaph281/Desktop/GitLab/2024_solanum_ldt_rnaseq/TAI/data/GenEra/ncbi_taxid.csv")

combined_data_edit_age <- combined_data_age %>% 
  left_join(ncbi_key, by="taxid") %>% 
  select(!X)

 

fulldata<- combined_data_edit %>% 
  left_join(phylorank, by="phylorank") %>%
  left_join(combined_data_edit_age, by=c("species", "phylorank", "taxid", "phylostratum")) %>%
  rename("Gene family\nfounder events"=number_of_genes_families, "Gene age"=number_of_genes) %>% 
  pivot_longer(cols=c("Gene family\nfounder events", "Gene age"), names_to="gene_type", 
               values_to="number_of_genes") 


svg(paste0("C:/Users/suaph281//Nextcloud/ResiDEvo/sequencing/figures/fig_1/", Sys.Date(),"_AGE_all.svg"), 
    width=12, height=10, bg="transparent")

(p3 <- fulldata %>% 
    #left_join(phylorank, by="phylorank") %>%
    arrange(phylorank) %>% 
    mutate(phylostratum_general = fct_reorder(phylostratum_general, phylorank)) %>% 
    ggplot(aes(x=as.factor(phylostratum_general), y=(number_of_genes), col=species, group=species))+
    geom_rect(xmin=6 - 0.5, xmax=6 + 0.5, 
              ymin=-Inf, ymax=Inf, fill="grey90", alpha=0.1, inherit.aes = FALSE) +
    geom_rect(xmin=16 - 0.5, xmax=16 + 0.5, 
              ymin=-Inf, ymax=Inf, fill="grey90", alpha=0.1, inherit.aes = FALSE) +geom_line(size=2, show.legend = F)+
    geom_point(size=4, alpha=.95)+
    theme_bw()+
    scale_color_manual(values = rev(cols))+
    labs(x="Phylorank",
         y="Counts", 
         color="Species")+
    scale_y_continuous(trans='log10')+
    theme(legend.position = c(.85, .9))+
    theme(axis.text.x = element_text(angle=45, hjust=1),
          axis.text=element_text(color="black", size=14),
          axis.title = element_text(color="black", size=17),
          legend.title = element_text(size=17, face="bold"),
          legend.text = element_text(size=14, face="italic"),
          strip.text = element_text(size=14, face="bold", color="black"))+
    facet_grid(rows = vars(gene_type), scales = "free_y")
  )
dev.off()

ggsave("C:/Users/suaph281/Nextcloud/ResiDEvo/sequencing/figures/fig_1/2024_09_19_AGE_summary.png", 
       p3, dpi=900, width=12, height=9)
