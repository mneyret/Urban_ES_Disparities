# Urban_ES_Disparities
These scripts allow to reproduce the results from the paper Neyret et al. 2026. 

- The file config.R lists the directories in which data, intermediate results, or figures should be read or saved.
- The file UrbanES_Recreationfunctions.R contains all the functions necessary in UrbanESDisparities_recmodel.R
- UrbanESDisparities_recmodel.R gathers all raw spatial data to generate rasters of park quality, and distance-weighted park quality
- The file UrbanESDisparities_datamerge.R combine these rasters (as well as a raster for heat mitigation) with INSEE socio-demographic data and outputs cleaned up rasters and vector files summarising data at the square (200mx200m) level.
- The file UrbanESDisparities_datamerge.R runs the actual analysis and outputs figures.

!!! Important.
The raw data used to run UrbanESDisparities_recmodel.R are NOT provided here and should be gathered independently; sources are listed below. Instead, only the raster and vector data (needed to run the final script, UrbanESDisparities_datamerge.R) are provided.
If you only wish to reproduce the final results (rather than applying the recreation model to another study site for instance), you should only run this last file.

### Details ###

## config.R.template
Fill this with correct paths in your setup, then save as config.R

## UrbanES_Recreationfunctions.R
Contains all functions to generate the recreation layers.

## UrbanESDisparities_recmodel.R
Relies on UrbanES_Recreationfunctions.R to compute the recreation layers.

Input: 
These should be acquired for the study area (you may extract the extent of the study area from the final recreation maps)

- communes_perimeter.shp

- Cosia/cosia_10m.tif
https://cartes.gouv.fr/rechercher-une-donnee/dataset/IGNF_COSIA. IGN, 2023. CoSIA [WWW Document]. (accessed 12.3.24).

- GlobalHumanSettlement/GHS_2154.tif. 
https://human-settlement.emergency.copernicus.eu/datasets.php. Schiavina, M., Melchiorri, M., Freire, S., 2023. GHS-DUC R2023A - GHS Degree of
Urbanisation Classification, application of the Degree of Urbanisation methodology
(stage II) to GADM 4.1 layer, multitemporal (1975-2030). European Commission,
Joint Research Centre (JRC).
https://10.2905/DC0EB21D-472C-4F5A-8846-823C50836305

INSEE/200m_selection_urb.shp. 
https://www.insee.fr/fr/statistiques/7655515.
Insee, 2019. Revenus, pauvreté et niveau de vie - Données carroyées. Dispositif Fichier
localisé social et fiscal (FiLoSoFi).

- Amenities/osm_parks.shp
- Amenities/osm_fountains.shp
- Amenities/osm_playgrounds.shp
- Amenities/osm_fitness.shp
- Amenities/osm_drinking_water.shp
-All acquired from Open Street Map data.

Amenities/road_noise.tif

Output:
Data/Temp_data/attractiveness.tif
Data/Temp_data/opportunities.tif 
Data/Temp_data/pollution.tif
Data/Temp_data/parks_quality.tif
Data/Temp_data/all_serv_parks.tif
Data/Temp_data/recre_decay_300_whole_region.tif