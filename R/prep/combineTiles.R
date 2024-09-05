## combine elevation 

files.raw <- list.files(path = "../../../../resources/spatial/Elevation/small_tiles/", pattern = ".tif", full.names = T)

r1 <- rast(files.raw[1])
r2 <- rast(files.raw[2])
r3 <- rast(files.raw[3])
r4 <- rast(files.raw[4])
r5 <- rast(files.raw[5])
r6 <- rast(files.raw[6])
r7 <- rast(files.raw[7])
r8 <- rast(files.raw[8])
# r9 <- rast(files.raw[9])
# r10 <- rast(files.raw[10])
# r11 <- rast(files.raw[11])
# r12 <- rast(files.raw[12])



ele.global <- merge(r1, r2, r3, r4, r5, r6, r7, r8, 
                filename = "../../../../resources/spatial/Elevation/global_elevation_450m_SRTM.tif", 
                overwrite = TRUE)

plot(ele.global)


