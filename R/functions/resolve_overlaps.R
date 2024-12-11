# function to remove overlapping polygons 

resolve_overlaps <- function(polygons = NA, keep_largest = TRUE, seed = 161){
  
  sf_use_s2(FALSE)
  overlap_list <- st_intersects(polygons, sparse = TRUE)
  
  #to make sure we only keep the largest
  if(keep_largest == TRUE){
    polygons$area <- as.numeric(st_area(polygons))  #
  }
 
  # Create a vector to track polygons to keep
  keep <- rep(TRUE, nrow(polygons))
  
  for (i in 1:length(overlap_list)) {
    
    if (keep[i]) {
      
      overlapping_indices <- overlap_list[[i]]
      
      # Exclude the current polygon from its overlaps
      overlapping_indices <- overlapping_indices[overlapping_indices != i]
      
      
      if(keep_largest == TRUE) {
        # Compare areas and mark smaller polygons to drop
        for (j in overlapping_indices) {
          if (polygons$area[i] >= polygons$area[j]) {
            keep[j] <- FALSE
          } else {
            keep[i] <- FALSE
            break
          }
        }
        
      }else{
        
        if(length(overlapping_indices) == 0){next}
        set.seed(seed)
        # Mark all as not kept initially
        keep[overlapping_indices] <- FALSE
        # Randomly select one to keep
        random_index <- sample(overlapping_indices, 1)
        keep[random_index] <- TRUE
      }
    }
 #   print(paste0(i, " out of ", nrow(polygons), " done"))
    
  }  
  
  sf_use_s2(TRUE)
  
  # Filter the polygons
  new_polygons <- polygons[keep, ] %>% dplyr::select(-area)
  
  print(paste0(
    "removed ", nrow(polygons) - nrow(new_polygons), " polgyons (",
    round((1 - (nrow(new_polygons)/nrow(polygons)))*100, 1), "% of the data) - ",
    nrow(new_polygons), " are left"
  ))
  return(new_polygons)
}
