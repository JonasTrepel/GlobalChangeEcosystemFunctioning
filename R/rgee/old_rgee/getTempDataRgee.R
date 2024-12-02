#install.packages("geojsonio")
library(rgee)
library(data.table)
library(tidyverse)
ee_Initialize(project = "ee-jonastrepel", drive = TRUE)


years <- c(1950:2023)

annual_temps <- data.table()

for(year in years){
  
  print(year)
  start <- paste0(year, "-01-01")
  end <- paste0(year, "-12-31")
  
  annualComp <- ee$
    #ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR')$
    ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR')$
    select('temperature_2m')$
    filterDate(start, end)$
    mean()$subtract(273.15)

  
  world_ext <- ee$Geometry$Rectangle(
    coords = c(-179.99999, -89.9999, 179.99999, 89.9999),
    proj = "EPSG:4326",
    geodesic = FALSE
  )
  
  acR <- ee_as_rast(image = annualComp,
                    region = world_ext,
                    container = "RGEE",
                    dsn = "annCompMat", 
                   # scale = 500,
                    timePrefix = FALSE,
                    maxPixels = 1e13,
                    via = "drive")
 
  plot(acR)
  
  drive_rm("annCompMat.tif")
  googledrive::drive_empty_trash()
  
  acDT <- as.data.frame(acR, xy = T)
  names(acDT) <- c("x", "y",  paste0("temperature_2m_", year))
  
  if(year == min(years)){
    annual_temps <- acDT
  }else{
    annual_temps <- left_join(acDT, annual_temps)
  }

}
  
fwrite(annual_temps, "data/ProcessedData/annual_mat_era5.csv")

