if(!require('remotes')) install.packages('remotes')
remotes::install_github("nevrome/ggpointgrid")

library(spData)
library(maps)
library(sf)
library(data.table)
library(tidyverse)

# Define levels and colours
os_colors <- c("Calculus - Ancient" = "#EF7C12",
               "Calculus - Industrial" = "#F4B95A",
               "Plaque - Industrial" = "#FC6882",
               "Plaque - Non-industrial" = "#C70E7B",
               "Buccal mucosa - Industrial" = "#8FDA04",
               "Buccal mucosa - Non-industrial" = "#009F3F",
               "Saliva - Industrial" = "#54BCD1",
               "Saliva - Non-industrial" = "#007BC3",
               "Tongue - Industrial" = "#78d9ca",
               "Hard palate - Industrial" = "#cab2d6",
               "Keratinized gingiva - Industrial" = "#F7E690",
               "RS - Industrial" = "gray80")

oral_site_levels <- c("Calculus - Ancient",
                      "Calculus - Industrial",
                      "Plaque - Non-industrial",
                      "Plaque - Industrial",
                      "Buccal mucosa - Non-industrial",
                      "Buccal mucosa - Industrial",
                      "Saliva - Non-industrial",
                      "Saliva - Industrial",
                      "Tongue - Industrial",
                      "Hard palate - Industrial",
                      "Keratinized gingiva - Industrial",
                      "RS - Industrial")

# Read in the cluster info file
cluster_info <- fread("./00-documentation/cluster_info_all_metadata_26_07_2024.tsv") %>%
  mutate(Oral_site = str_replace_all(Oral_site, "anc_calc$", "AC"),
         Oral_site = str_replace_all(Oral_site, "anc_calculus", "AC"),
         Oral_site = str_replace_all(Oral_site, "mod_calculus", "MC"),
         Oral_site = str_remove_all(Oral_site, "sub"),
         Oral_site = str_remove_all(Oral_site, "sup")) |>
  left_join(fread("./countries_coordinates.txt"), by = c("Country" = "country"))

# Read in ancient environmental bins file for filtering
ancient_envt_bins <- fread("./00-documentation/ancient_envt_clusters_26_07_2024.tsv")

# Prepare metadata for countries
country_os_meta <- cluster_info |>
    select(sample_host, country = Country, latitude, longitude,
           Continent, Secondary_sample_accession, primary_cluster, Oral_site, Industrialization) |>
    mutate(Oral_site = recode(Oral_site,
                              `AC` = "Calculus",
                              `MC` = "Calculus",
                              `PQ` = "Plaque",
                              `BM` = "Buccal mucosa",
                              `SL` = "Saliva",
                              `TN` = "Tongue",
                              `HP` = "Hard palate",
                              `KG` = "Keratinized gingiva",
                              `RefSeq` = "RS"),
           Industrialization = recode(Industrialization,
                                      `Non_industrial` = "Non-industrial")) |>
    unite(os, Oral_site:Industrialization, sep = " - ") |>
    mutate(sample_host = if_else(str_detect(sample_host, "Homo"),
                                  "Human", "Other"),
           Continent = recode(Continent,
                              `Australia` = "Oceania",
                              `Carribean` = "South America")) |>
    as_tibble()

# Count the number of samples per country
country_n_samples <- country_os_meta |>
  summarise(No_samples = length(unique(Secondary_sample_accession)),
            .by = c(Continent, country, os, latitude, longitude, sample_host))

# Count the number of MAGs per country
country_n_mags <- country_os_meta |>
  anti_join(ancient_envt_bins) |>
  summarise(No_MAGs = n(), .by = c(country, os))

# Join the number of samples and the number of MAGs per country
country_data <- full_join(country_n_samples, country_n_mags) |>
  mutate(os = fct_relevel(os, oral_site_levels),
         No_MAGs = if_else(country == "Croatia", 0, No_MAGs))

arrange(country_data, Continent, country, os, sample_host) |>
(\(x) fwrite(x, file = "country_nsamples_nmags.tsv", sep = "\t"))()

os_mags_world <- ggplot(world) +
  geom_sf(aes(geometry = geom), fill = "gray90", color = "gray90") +
  ggpointgrid::geom_pointrect(data = country_data,
                              aes(longitude, latitude, color = os, size = No_MAGs, shape = sample_host),
                              scale_x = 10.0, scale_y = 10.0, stroke = 0.8) +
  scale_color_manual(values = os_colors[oral_site_levels]) + # ,
  scale_size_continuous(range  = c(1, 6), 
                        limits = c(0, 4000), 
                        breaks = c(0, 500, 1000, 2000, 3000)) +
  scale_shape_manual(values = c(19, 17)) +
  theme_void() +
  labs(color = "Sample type", size = "No. MAGs", shape = "Sample host") +
  theme(panel.grid.major = element_blank(),
        legend.position = "right")
os_mags_world

# ggsave("os_mags_world_shapes_allsamples.pdf", plot = os_mags_world, device = "pdf",
#        scale = 1, width = 10, height = 5, units = c("in"), dpi = 300, useDingbats = FALSE)
