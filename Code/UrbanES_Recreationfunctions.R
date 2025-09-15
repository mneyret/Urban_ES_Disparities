# slu = Segregated LULC map. Each layer is a landuse. Raster map (spatRaster)
# parks = Patches of the green areas of interest, they are numbered and each associated to its number. Raster map (SparRaster)
# bfp = Same as "parks", but with a 10m buffer around each patch/park. Raster map (SparRaster)
# npatches = number of parks/patches. Numeric
# 
# luH = Land-use numbers for herbacious vegetation. Vector of characters
# luT = Land-use numbers for tree-filled landuses. Vector of characters
# luW = Land-use for waterscapes. Vector of characters
# luB = land-use number corresponding to buildings. Vector of characters
# 
# 
# ftn = Map for the fountains. Point Layer (SpatVector)
# am_int = point layer maps list for the ammenities of interest. List of  Point Layers (SpatVectors)
# 
# luna = land use numbers for the noise-abatting land uses. Vector of characters
# matw = reclassification matrix linking land-uses (those in luna) abd their weight for noise-abattment. 2-column dataframe
# 
# lvl_road = map of the roads. Roads' values are their noise level. Raster map (spatRaster)
# lvls = list of sound levels. Must be the same length as the number of categories in lvl-road. Vector of numerics
# 
# all_serv = Either file name for the model result, or the raster layer straight from the r environment
# outdir = Directory to save the results in

decay = function(d, half){
  x = exp(-log(2)*d/half)
}

scaleraster = function(x){
  nx <- minmax(x)    
  x <- (x - nx[1,]) / (nx[2,] - nx[1,])
}

draw_window_decay <- function(radius_in_pixels, resolution, decay_half_dist) {
  grid <- matrix(NA, 2*radius_in_pixels+1, 2*radius_in_pixels+1)
  grid[1*radius_in_pixels+1, 1*radius_in_pixels+1] = 1
  r = terra::rast(grid)
  dist = terra::distance(r) * resolution
  decay_d = decay(dist, decay_half_dist)
  # Remove further than 1k
  decay_d = decay_d*as.numeric(dist<1*radius_in_pixels*resolution)
  decay_mat = as.matrix(decay_d, wide = T)
  return(decay_mat)
} # Function that draws a circle


visual_quality_distance = function(slu, parks, bfp, luH, luT, luW, ftn, grid, buff) {
  ### Attractiveness of the greenery within the buffer ###
  
  # Raster layer for the grass cover
  herb = c()
  for(i in luH) {
    herb = c(herb,slu[[i]]) # Grouping all the landuse for grasses, in a loop for flexibility
  }
  herb = sprc(herb) # Making it a SpatRaster collection that can be merged
  herb = terra::merge(herb) 
  herb[!is.na(herb)] = 1
  #herb = mask(parks, herb) #Keeping only the herbs in the parks, and associating them with the patch ID
  
  # Raster layer for the tree cover
  tree = c()
  for(i in luT) {
    tree = c(tree,slu[[i]])
  }
  tree = sprc(tree)
  tree = terra::merge(tree)
  tree[!is.na(tree)] = 1
  
  # trick to have 0 where all parks are
  parks0 = parks
  parks0[parks0>0] = 0
  
  print("a")
  #tree = mask(parks, tree)
  cover = sum(parks0, herb, 2*tree, na.rm = T)
  cover = terra::mask(cover, parks)

  compo_fun = function(window){
    counts = table(window)
    counts = counts/sum(counts)
    counts[counts>0.5] = 0.5
    return(sum(counts[names(counts)!=0]))
  }
  

  buffer = focalMat(cover, 300, type=c('circle'), fillNA=FALSE) # A focal window with a 300m radius
  
  score_vege <- focal(cover, buffer, fun=compo_fun, na.policy = "omit") 

  print("b")
  #### Try loop to get greenery within 300m in each park without overlap
  #attract = parks
  #values(attract) = NA
  #
  #for (i in unique(values(parks), na.rm = T)){
  #  print(i)
  #  Mask = mask(minicov, parks == i, maskvalues = c(NA, 0))
  #  f <- focal(Mask, w=49, fun=compo_fun, na.policy = "omit") 
  #  attract = sum(attract,f, na.rm = T)
  #}
  
  
  ### Parc size ###
 # alpha = cellSize(parks,unit="ha") |> zonal(parks, sum) #Gives a dataset linking patches to their area
#  names(alpha) = c("patches","park")
#  alpha$park[is.na(alpha$park)] = 0
  
  
  ### Waterscape-related attractiveness ###
  fountains = rasterizeGeom(ftn, grid, fun = "count") # The points related to fountains are rasterized 
  
  waters = c() # Waters will be the compiled raster of all the LU codes describing waterscapes
  for(i in luW) {
    waters = c(waters,slu[[i]])
  }
  waters = sprc(waters)
  waters = terra::merge(waters)
  
  print("c")
  
  waters = merge(waters, fountains) # Fountains and waterscapes are merged in a single raster layer
  waters[waters == 0] = NA

  DistWater = distance(waters)

  DecayWater = decay(DistWater, half = 300)
  plot(DecayWater)
  score_waters = terra::mask(DecayWater, parks)
  print("d")
  # Total attractiveness score
  score_att = (score_vege + score_waters) /2 # Global attractiveness score
  
  # attractiveness is raster in which each layer represents one of the parameters
  attractiveness = c(score_vege, score_waters, score_att)
  names(attractiveness) = c("score_vege", "score_waters", "score_att")
  
  return(attractiveness)
} 

recre_opp_distance = function(parks, npatches, am_int, grid) {
  ### Parc size ###
  alpha = cellSize(parks,unit="ha") |> zonal(parks, sum) #Gives a dataset linking patches to their area
  names(alpha) = c("patches","park")
  alpha$park[is.na(alpha$park)] = 0
  alpha$park[alpha$park > 2] = 2 #The areas are maxed at 1, giving the area-related score
  park_size = classify(parks, alpha)/2 #Each park's value will be its "area score"
  
  
  ## The amenities
  
  ## Uncomment this section for ammenity-related score limited to 300m ##
  pa = parks 
  pa[!is.na(pa)] = 2 # A single value is attributed to all the parks
  
  resl = rast() # the raster that will store the results
  
  for (i in 1:length(am_int)){ # one loop for each amenity
    rastam = rasterizeGeom(am_int[[i]], grid, fun = "count") # The points of each amenity are rasterized
    rastam[rastam == 0] = NA
    rastam[!is.na(rastam)] = 1 # The pixel values will be either 1 or 0, depending on the presence of the amenity in them
    mg = terra::mask(rastam, parks) # The parks' limits and the amenity's presence are put on a single map
    mg = distance(mg) # Distance to the amenity is calculated in each park
    Decay = decay(mg, 300)
    add(resl) = Decay # Storage of the decay map
  }
  
  am = sum(resl)/3 # Each pixel will show the number of amenities contained within 500m
  am = terra::mask(am, parks)
  am[am == 0] = NA
  
  tot_recre = (park_size + am)/2 # Recreation opportunities score is the mean of the size and the amenities
  
  opportunities = c(park_size, am, tot_recre) #Raster for which each layer is a parameter for the recreation opportunities
  names(opportunities) = c("parc_size", "score_am", "tot_recre")
  
  return(opportunities)
}

pollut = function(parks, slu, luna, luB, lvls, matw, lvl_road) {
  #Visual Pollution
  parbuil = slu[[luB]] # Raster map for the buildings
 # parbuil[is.na(parbuil)] = 0
  
  buildis = distance(parbuil#, #target = as.numeric(luB)
                     ) # distance to buildings
  
  DecayBuild = decay(buildis, half = 300)
  
  # buildis[buildis > 500] = 500
  # buildis = 1 - (buildis / 500) # Distance is maxed at 500m and brought back to a score from 0 to 1

    # High score = high visual pollution (close to building)
  
  # iota = zonal(buildis, parks, fun = "min") # uncomment for whole-park values, comment for for per-pixel values
  # buildis = classify(parks, iota) # uncomment for whole-park values, comment for for per-pixel values
  
  buildis = terra::mask(DecayBuild, parks) # uncomment for per-pixel values, comment for whole-park values
  
  #Noise Pollution
  wght = c() # empty vector that will store the land use that will impact sound propagation
  for(i in luna) {
    wght = c(wght,slu[[i]])
  }
  wght = sprc(wght) # Making it a SpatRaster collection that can be merged
  wght = terra::merge(wght)
  wght = classify(wght, matw) # reclassifying those land use with their weights
  
  
  # Sound levels associated with "road types"
  
  sources = segregate(lvl_road, classes = NULL, keep = FALSE, other = NA) # Separating each "road type"
  sources[!is.na(sources)] = 100 # Roads will be the target, their weight doesn't matter
  sources = merge(sources, wght, first = T) # For each "road type", making a map with this road (valued at 100) and the weighted land uses
  sources[is.na(sources)] = 1 
  
  noise = rast()  # One loop for each road type, results will be stored here
  for (i in 1:length(lvls)){
    print(paste0("boucle ",i))
    nadd = costDist(sources[[i]], target = 100, scale = 1, overwrite = T) # Very heavy function.
    nadd[nadd<1] = 1 # to avoid Inf values after the log()
    nadd = lvls[i] - 10*log10(nadd) # distance-based sound attenuatiion
    add(noise) = nadd
    rm(nadd)
    gc() #small function to clear the unused R memory
  }
  
  park_noise = max(noise) / max(lvls) # Taking the maximum sound level on every point, and making it a ratio on the max sound level
  # High score = high noise pollution
  
  
  # lambda = zonal(park_noise, parks, fun = "min") # uncomment for whole-park values, comment for per-pixel values
  # park_noise = classify(parks, lambda) # uncomment for whole-park values, comment for per-pixel values
  
  park_noise = terra::mask(park_noise, parks) # uncomment for per-pixel values, comment for whole-park values
  
  # The final score
  score_pol = (park_noise + buildis)/2 # Mean of the two pollution scores. High score = high pollution
  pollutions = c(buildis, park_noise, score_pol) # pollution is the raster for which each layer is a parameter of the pollution index
  names(pollutions) = c("dist_build", "noise", "tot_poll")
  
  return(pollutions)
}


serv_tot = function(attractiveness, opportunities, pollutions) {
  # Simply combining all services in a single raster
  # Each layer is one of a parameter
  serv_cult = (scaleraster(attractiveness[[3]]) + 
                 scaleraster(opportunities[[3]]) +
                 (1 - scaleraster(pollutions[[3]])))/3
  
  all_serv = c(attractiveness[[3]], opportunities[[3]], pollutions[[3]], serv_cult)
  names(all_serv) = c("attractiveness", "opportunities", "pollutions", "total_recre")
  
  return(all_serv)
}

qual_decay = function(parks_quality, radius_in_pixels, decay_half_dist){ 
  serv = parks_quality
  buff = buffer(serv, radius_in_pixels*res(serv)[1])
  serv[is.na(serv)] = 0
  serv = terra::mask(serv, buff, maskvalues=FALSE)
  window = draw_window_decay(radius_in_pixels, res(serv)[1], decay_half_dist)
  quali_decaydist = focal(serv, w=window, 
                                  na.policy="omit", 
                                  fun=function(x){quantile(x[x!=0], 0.99, na.rm = T)})
  quali_decaydist = crop(quali_decaydist, parks_quality)
  #terra::writeRaster(quali_decaydist, outdir, datatype = "FLT8S", overwrite = TRUE)
  return(quali_decaydist)
}


