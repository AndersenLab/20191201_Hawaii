# load packages
library(tidyverse)
library(ggmap)
library(memoise)
library(lubridate)
library(cowplot)
library(pals)
library(grid)
library(gridExtra)
library(scales)
library(kableExtra)
library(leaflet)
library(htmlwidgets)
library(htmltools)

# Set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# load fulcrum data
load('data/fulcrum/fulcrum_dat.Rda')

# define functions
#################################################
### set color palettes                        ###
#################################################
species_palette <- c("C. elegans" = "#BE0032", #7
                     "C. kamaaina" =  "#875692", #4 
                     "C. tropicalis" = "#F38400", #5 
                     "Panagrolaimus sp." = "#C2B280", #8
                     "Oscheius sp." = "#F3C300", #3
                     "C. briggsae" = "#A1CAF1", #6 
                     "Other PCR +" = "#008856", #10
                     "PCR -" = "#848482", #9
                     "Not genotyped" = "#F2F3F4", #1
                     "Tracks only" = "#b3b3b3", #manual
                     "No Nematode" = "#222222")  #2

island_palette <- c("Kauai" = "#E69F00",
                    "Oahu" = "#56B4E9",
                    "Molokai" = "#009E73",
                    "Maui" = "#F0E442",
                    "Big Island" = "#D55E00")

#################################################
### Define Functions                          ###
#################################################
# function used in map_overview to set as numeric
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

# Function for setting islands 
filter_box <- function(longitude, latitude, coords) {
  between(longitude, coords[1], coords[3]) &
    between(latitude, coords[2], coords[4]) &
    !is.na(longitude)
}

# coordinates for islands. Note need to comment out islands not collected from in this list
islands = list(  "Kauai" = c(-159.830818,21.750571,-159.230003,22.350076),
                 "Oahu" = c(-158.323116,21.112767,-157.623081,21.814254),
                 #"Molokai" = c(-157.3515,20.793,-156.6515,21.4956),
                 "Maui" = c(-156.745977,20.405495,-155.942774,21.207099),
                 "Big Island" = c(-156.3651,18.8049,-154.765,20.4064)
)

#gtmap function for map_overview function
gtmap <- function(loc) {
  get_map(location = loc,
          maptype = "terrain-background",
          source = "stamen",
          scale = "auto")
}
mget_map <- memoise(gtmap)

# NOTE: You must configure set_islands to what was sampled in this project for this script to work (see lines150 - 163).
# Overview map plotting function
map_overview <- function(F, label, cso, geoms, face = "plain") {
  island_set = lapply(names(islands), function(i) {
    
    l_position = "none"
    island_size = 2
    imap <- cso %>% dplyr::filter(collection_island == i)
    rects = element_blank()
    map = mget_map(islands[[i]])
    
    # Calculate scalebar
    bb <- attr(map,"bb")
    sbar <- data.frame(lon.start = c(bb$ll.lon + 0.1*(bb$ur.lon - bb$ll.lon)),
                       lon.end = c(bb$ll.lon + 0.25*(bb$ur.lon - bb$ll.lon)),
                       lat.start = c(bb$ll.lat + 0.1*(bb$ur.lat - bb$ll.lat)),
                       lat.end = c(bb$ll.lat + 0.1*(bb$ur.lat - bb$ll.lat)))
    
    sbar$distance <- geosphere::distVincentyEllipsoid(c(sbar$lon.start,sbar$lat.start),
                                                      c(sbar$lon.end,sbar$lat.end))
    
    scalebar.length <- 20
    sbar$lon.end <- sbar$lon.start +
      ((sbar$lon.end-sbar$lon.start)/sbar$distance)*scalebar.length*1000
    ptspermm <- 2.83464567
    
    base_map <- ggplot(imap) +
      ggmap::inset_ggmap(map) +
      #rects +
      geoms +
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
            axis.ticks.y = element_blank(), axis.text.y = element_blank(),
            axis.title.x = element_blank(), axis.title.y = element_blank(),
            panel.background = element_blank(),
            panel.spacing = unit(c(0,0,0,0), "lines"),
            axis.line = element_blank(),
            plot.title = element_text(lineheight=.8, face="bold", vjust=1),
            plot.margin = unit(c(0,0,0,0), "lines"),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = l_position,
            legend.background = element_rect(fill="white"),
            legend.text=element_text(size=12, color = "black", face = face)) +
      coord_equal(ratio=1) +
      scale_x_continuous(limits = islands[[i]][c(1,3)], expand = c(0, 0)) +
      scale_y_continuous(limits = islands[[i]][c(2,4)], expand = c(0, 0)) +
      geom_segment(data = sbar,
                   aes(x = lon.start,
                       xend = lon.end,
                       y = lat.start,
                       yend = lat.end),
                   arrow=arrow(angle = 90, length = unit(0.1, "cm"),
                               ends = "both", type = "open")) +
      geom_text(data = sbar,
                aes(x = (lon.start + lon.end)/2,
                    y = lat.start + 0.025*(bb$ur.lat - bb$ll.lat),
                    label = paste(format(scalebar.length),
                                  'km')),
                hjust = 0.5,
                vjust = 0,
                size = 8/ptspermm)  +
      coord_map(projection = "mercator",
                xlim=c(bb$ll.lon, bb$ur.lon),
                ylim=c(bb$ll.lat, bb$ur.lat)) +
      scale_radius(range = c(island_size, island_size), guide = "none") #+
    #scale_shape_manual(values = shape)
    
    base_map
    
  })
  
  island_set[[6]] <- F
  
  without_label <- plot_grid(plotlist = island_set,
                             labels = c("A - Kauai",
                                        "B - O'ahu",
                                        #"C - Moloka'i",
                                        "C - Maui",
                                        "D - Island of Hawai'i"
                                        #"D", ""
                             ),
                             label_y = 0.98,
                             hjust = 0,
                             label_x = 0.06,
                             align = "vh")
  
  cowplot::plot_grid(without_label, label, nrow = 2, rel_heights = c(1, .05))
}

# Map_collection function for reviewing collection locations
map_collection <- function(df, color_use) {
  
  icos <- iconList(
    red = makeIcon(
      iconUrl = paste0("https://storage.googleapis.com/andersenlab.org/img/red.svg"),
      iconWidth = 20, iconHeight = 20,
      popupAnchorX = 0.001, popupAnchorY = -20,
      iconAnchorX = 20/2, iconAnchorY = 20
    ),
    green = makeIcon(
      iconUrl = paste0("https://storage.googleapis.com/andersenlab.org/img/green.svg"),
      iconWidth = 10, iconHeight = 10,
      popupAnchorX = 0.001, popupAnchorY = -10,
      iconAnchorX = 10/2, iconAnchorY = 10
    )
  )
  df <- dplyr::filter(df, !is.na(df[[color_use]])) %>%
    dplyr::mutate(substrate=ifelse(is.na(substrate), "", substrate)) %>%
    dplyr::arrange(collection_lat_long_method)
  
  #print(df)
  
  
  attach(df)
  leaflet::leaflet(data = df, width = "100%", options = list(zoomControl = T)) %>% 
    leaflet::addTiles( 
      paste0( 
        "https://stamen-tiles-{s}.a.ssl.fastly.net/terrain/{z}/{x}/{y}.png",
        jsonlite::read_json("data/thunderforest.json")$key)  
    ) %>%
    leaflet::addMarkers(~collection_longitude,
                        ~collection_latitude,
                        popup = glue::glue("<h2>{collection_id}</h2><hr />
                                           <strong>collection uplaoded by:</strong> {collection_by}<br />
                                           <strong>latitidue, longitude:</strong> {format(round(collection_latitude, 6), nsmall = 6)}, {format(round(collection_longitude, 6), nsmall = 6)}<br />
                                           <strong>postion method used:</strong> {collection_lat_long_method}<br />
                                           <strong>local time:</strong> {collection_local_time}<br />
                                           <strong>altitude:</strong> {altitude} meters<br />
                                           <strong>landscape:</strong> {landscape}<br /><br />"),
                        popupOptions(maxWidth = 500),
                        icon = icos[ df[[color_use]] ] )
  
  #htmlwidgets::saveWidget(m, tempfile(), selfcontained = FALSE)
  #webshot::webshot("temp.html", file = "map.png",
  #        cliprect = "viewport", vwidth = 1000, vheight = 1000)
}

# plot islands
overview_plot_df <- fulcrum_dat %>%
  dplyr::mutate(collection_type = ifelse(worms_on_sample %in% c("No", "?"), "No Nematode",
                                         ifelse(worms_on_sample == "Tracks", "Tracks only",
                                                ifelse(worms_on_sample == "Yes" & is.na(ITS2_pcr_product), "Not genotyped",
                                                       ifelse(worms_on_sample == "Yes" & ITS2_pcr_product == 0, "PCR -",
                                                              ifelse(species_id %in% c("Chabertia ovina",
                                                                                       "Choriorhabditis cristata",
                                                                                       "Choriorhabditis sp.",
                                                                                       "Heterhabditis zealandica",
                                                                                       "Mesorhabditis sp.",
                                                                                       "no match",
                                                                                       "C. kamaaina",
                                                                                       "Rhabditis terricola",
                                                                                       "Rhanditis tericola",
                                                                                       "Teratorhabditis sp.",
                                                                                       "Unknown",
                                                                                       "unknown",
                                                                                       "Oscheius sp.",
                                                                                       "Panagrolaimus sp.",
                                                                                       "-",
                                                                                       NA),
                                                                     "Other PCR +", species_id)))))) %>%
  dplyr::select(collection_id, collection_type, collection_island, collection_location, species_id, substrate, ITS2_pcr_product, worms_on_sample, collection_longitude, collection_latitude) %>%
  dplyr::distinct(collection_id, collection_type, .keep_all=T) %>% 
  dplyr::group_by(collection_id) %>%
  dplyr::mutate(multiple_type = ifelse(n() > 1, "yes", "no")) %>%
  dplyr::ungroup() %>%
  # optional rename
  dplyr::mutate(collection_type = ifelse(collection_type == "Caenorhabditis kamaaina", "C. kamaaina", collection_type)) %>%
  dplyr::mutate(collection_type = forcats::as_factor(collection_type),
                collection_type = forcats::fct_relevel(collection_type,
                                                       "C. elegans",
                                                       "C. kamaaina",
                                                       "C. tropicalis",
                                                       "C. briggsae",
                                                       "Other PCR +",
                                                       "PCR -",
                                                       "Not genotyped",
                                                       "Tracks only",
                                                       "No Nematode")) %>%
  dplyr::arrange(collection_type) %>% # arrange sets order for collection_id with multiples so highest priority collection type is on top
  dplyr::distinct(collection_id, .keep_all = T) %>% # selects highest priority collection type from a c-label with multiple collection types on it
  dplyr::arrange(desc(collection_type)) %>% # reorders collection_id so highest priority collections are plotted on top
  dplyr::filter(!is.na(collection_id)) %>% # remove any NAs in collection _id
  dplyr::filter(collection_type %in% c("C. elegans", "C. tropicalis", "C. briggsae"))

####################################################
#  Bar chart inset                                 # 
####################################################
# Plot stacked bar chart
bar_chart <- overview_plot_df %>%
  dplyr::group_by(collection_type, collection_island) %>%
  dplyr::mutate(collections_per_island = n()) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(collection_island) %>%
  dplyr::mutate(total_collections = n(), perc_class_island = collections_per_island / total_collections * 100) %>%
  dplyr::arrange(total_collections) %>%
  #dplyr::select(fixed_substrate, plot_class, total_substrates, perc_worm_sub, worm_per_substrate) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(collection_type, collection_island, .keep_all = T) %>%
  dplyr::mutate(collection_island = factor(collection_island, levels = names(island_palette))) %>%
  dplyr::mutate(collection_type = factor(collection_type, levels = c("C. tropicalis", "C. briggsae", "C. elegans")))

# Bar chart for label
plot_bar_chart <- ggplot(data = bar_chart) +
  geom_bar(stat = "identity", aes(x = factor(collection_island), y = perc_class_island, fill = collection_type), colour = "black") +
  scale_fill_manual(values=c(species_palette)) +
  theme_bw() +
  theme(axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        legend.position = "bottom",
        plot.margin = unit(c(0,0,0,0), units = "cm")) + 
  #legend.text = element_text(size = 8, color = "black")) +
  labs(fill = "", x = "", y = "Percentage of all collections") +
  geom_text(aes(x=collection_island, y=102, label=paste0("n=",total_collections)), 
            position = position_dodge(width=1), size = 2.5) +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_y_continuous(breaks = c(25, 50, 75, 100), limits = c(0, 102))

plot_bar_chart_no_legend <- plot_bar_chart +
  theme(legend.position="none",
        plot.margin = unit(c(0,0,0,0), units = "cm"))

plot_bar_chart_legend <- cowplot::plot_grid(cowplot::get_legend(plot_bar_chart))


# Make map overview plot
overview_plot <- map_overview(F = NULL, plot_bar_chart_legend, cso = overview_plot_df, 
                              geoms = c(geom_point(aes(x=collection_longitude,
                                               y=collection_latitude,
                                               fill=collection_type,
                                               size = 1),
                                           color="black",
                                           shape=21,
                                           stroke = 0.5
                              ),
                              scale_fill_manual("species", values = species_palette)
                              ),
                              face="italic"
)
overview_plot
ggsave('plots/DLNR_collection_map.png', width = 6.65, height = 7.5)

# Export DLNR dataframe for elegans, briggsae, tropicalis isolations
DLNR_isolations <- fulcrum_dat %>%
  dplyr::filter(!is.na(collection_id)) %>%
  dplyr::mutate(substrate = ifelse(is.na(substrate_other), substrate,
                                   ifelse(substrate_other == "bark", "Vegetation", substrate))) %>%
  dplyr::filter(species_id %in% c("C. tropicalis", "C. briggsae", "C. elegans")) %>%
  dplyr::select(species_id, 
                collection_id,
                isolation_id,
                substrate,
                collection_island,
                collection_trail,
                collection_longitude,
                collection_latitude) %>%
  dplyr::arrange(species_id, desc(collection_latitude))
  
readr::write_csv(DLNR_isolations, "data/DLNR_isolation_report.csv")

# Export loactions of collections with elegans, briggsae, or tropicalis
readr::write_csv( fulcrum_dat %>%
                    dplyr::filter(!is.na(collection_id)) %>%
                    dplyr::filter(species_id %in% c("C. tropicalis", "C. briggsae", "C. elegans")) %>%
                    dplyr::distinct(collection_id, species_id, .keep_all=T) %>%
                    dplyr::group_by(collection_id) %>%
                    dplyr::summarise(species = paste(species_id, collapse = ", ")) %>%
                    dplyr::left_join(DLNR_isolations %>% dplyr::distinct(collection_id, .keep_all = T)) %>%
                    dplyr::select(species, 
                                  collection_id,
                                  substrate,
                                  collection_island,
                                  collection_trail,
                                  collection_longitude,
                                  collection_latitude) %>%
                    dplyr::arrange(species, desc(collection_latitude)), "data/DLNR_collection_report.csv")

# numbers for DLNR report
collections <- fulcrum_dat %>%
  dplyr::filter(!is.na(collection_id)) %>%
  dplyr::distinct(collection_id) %>%
  dplyr::mutate(num_collections = n()) %>%
  dplyr::distinct(num_collections) %>%
  dplyr::pull(num_collections)

isolations <- fulcrum_dat %>%
  dplyr::filter(!is.na(collection_id)) %>%
  dplyr::filter(!is.na(isolation_id)) %>%
  dplyr::distinct(isolation_id) %>%
  dplyr::mutate(num_isolations = n()) %>%
  dplyr::distinct(num_isolations) %>%
  dplyr::pull(num_isolations)

`SSU + and ITS2 -` <- fulcrum_dat %>%
  dplyr::filter(!is.na(collection_id)) %>%
  dplyr::filter(ITS2_pcr_product == 0 & rhabditid_pcr_product == 1) %>%
  dplyr::mutate(num = n()) %>%
  dplyr::distinct(num) %>%
  dplyr::pull(num)

`SSU + and ITS2 +` <- fulcrum_dat %>%
  dplyr::filter(!is.na(collection_id)) %>%
  dplyr::filter(ITS2_pcr_product == 1 & rhabditid_pcr_product == 1) %>%
  dplyr::mutate(num = n()) %>%
  dplyr::distinct(num) %>%
  dplyr::pull(num)

`C. elegans` <- fulcrum_dat %>%
  dplyr::filter(!is.na(collection_id)) %>%
  dplyr::filter(species_id == "C. elegans") %>%
  dplyr::mutate(num = n()) %>%
  dplyr::distinct(num) %>%
  dplyr::pull(num)

`C. briggsae` <- fulcrum_dat %>%
  dplyr::filter(!is.na(collection_id)) %>%
  dplyr::filter(species_id == "C. briggsae") %>%
  dplyr::mutate(num = n()) %>%
  dplyr::distinct(num) %>%
  dplyr::pull(num)

`C. tropicalis` <- fulcrum_dat %>%
  dplyr::filter(!is.na(collection_id)) %>%
  dplyr::filter(species_id == "C. tropicalis") %>%
  dplyr::mutate(num = n()) %>%
  dplyr::distinct(num) %>%
  dplyr::pull(num)

DLNR_table <- tibble(collections,
                     isolations,
                     `SSU + and ITS2 -`,
                     `SSU + and ITS2 +`,
                     `C. elegans`,
                     `C. briggsae`,
                     `C. tropicalis`) %>%
  gather("class", "count")

# export DLNR table
readr::write_csv(DLNR_table, "data/DLNR_counts_table.csv")
