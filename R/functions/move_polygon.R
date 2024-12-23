movePolygon <- function(ogPoly, potSpace, maxAttempts = 1000, seed = 161) {
 
   # Extract the geometry of the input polygon
  set.seed(seed)
  
  newPoly <- sf::st_geometry(ogPoly)
  
  # Calculate the centroid of the bounding polygon and the input polygon
  polyCenter <- sf::st_centroid(newPoly)
  potSpaceCenter <- sf::st_centroid(potSpace)
  
  centeredPoly <- newPoly - polyCenter
  
  ## Try to place the polygon shape in the potential space. If it hasn't worked in a 1000 attempts, stop it
  attempt <- 0
  placedPoly <- NULL
  
  while (attempt < maxAttempts) {
    attempt <- attempt + 1
    print(attempt)
    
    # Generate random translation within the bounding polygon's bounding box
    boundingBox <- sf::st_bbox(potSpace)
    xShift <- runif(1, boundingBox["xmin"], boundingBox["xmax"])
    yShift <- runif(1, boundingBox["ymin"], boundingBox["ymax"])
    
    # Translate the polygon to the new random location
    translatedPolyRaw <- centeredPoly + c(xShift, yShift)
    translatedPoly <- st_as_sf(translatedPolyRaw, crs = st_crs(potSpace))
    
    # 
    if (all(sf::st_within(translatedPoly, potSpace, sparse = FALSE))) {
      placedPoly <- translatedPoly
      break
    }
  }
  
  # Return the placed polygon if successful, otherwise return original
  if (!is.null(placedPoly)) {
    return(placedPoly)
  } else {
    warning("Could not place the polygon within the bounding polygon after maximum attempts.")
    return(NULL)
  }
}
