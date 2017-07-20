#' Extract bands from multiple ASTER L1T hdf files
#'
#' Directly extract all bands in ASTER L1T hdf file/files in the input folder.
#'
#' @param directory Path to ASTER hdf file/ files folder. Assign as string (inside double quotes).
#' Either provide the path inside function or set up satellite image folder as the current working directory
#' before running the function.
#' To define current working directory, either use shortcut key Ctrl+Shift+H  or use \code{\link{setwd}} funtion.
#' @param ext2crop,crop Same as mentioned in \code{\link[ASIP]{arvi}}.
#' @return File named "arvi_'date of satellite image acqisition'.tif" in the input folder.
#' @note 1. Windows users users should be careful while assigning directory. Use "/" to seperate folders not "\\".
#'
#' 2. The base of this function is provided by Land processes distributed active archive center (LP DAAC).
#' Raw code is customized to produce this function with additional functionalities and more ease of use. Authors are thankful to LP DAAC, NASA and USGS.
#' @export
#' @importFrom raster raster writeRaster extent mask crop plot drawExtent
#' @import gdalUtils
aster_bands <- function(directory = getwd(), crop = "n", ext2crop = "none"){
# Create a list of ASTER L1T HDF files in the directory

file_list <- list.files (path = directory, pattern = 'AST_L1T_.*hdf$')
if (length (file_list) == 0)
  stop("Define your satellite image folder path properly")
# Defining the crop extent
ext <- 0
if (crop=="y" || crop=="f")
{
  if (typeof(ext2crop)== "character")
  {
    if (ext2crop!="none")
    {
      shape <- raster::shapefile(ext2crop)
      ext <- raster::extent(shape)
    }
    if (ext2crop=="none" & (crop=="y" || crop=="f")) {print("The extent shapefile is not defined properly")}
  }
  if (typeof(ext2crop)== "S4")
  {
    ext <- raster::extent(ext2crop)
    shape <- ext2crop
  }
}
# Create and set output directory
#-------------------------------------------------------------------------------
for (i in 1 : length (file_list)){
  # Maintains the original filename
  path_file <- paste0 (directory, "/", file_list [i])
  file_name <- file_list [i]
  # Read in the metadata
  meta_data <- gdalinfo (path_file)

  line_date <- meta_data [grep ("CALENDARDATE", meta_data)]
  splits_line <- strsplit (line_date, "=") [[1]]
  date <- splits_line [2]
  date <- paste0 (substr (date, 1, 4), "-", substr (date, 5, 6), "-", substr (date, 7, 8))

  line_dir <- meta_data [grep ("FLYINGDIRECTION", meta_data)]
  splits_line <- strsplit(line_dir, "=") [[1]]
  dirn <- splits_line [2]

  out_dir <- paste0 (directory, "/", "output_bands/")
  suppressWarnings (dir.create (out_dir))

  out_dir <- paste0 (out_dir, date, "_", dirn, "/")
  suppressWarnings (dir.create (out_dir))
  # Define CRS
  # Define Upper left and lower right--need for x, y min/max
  # For offset (pixel size / 2), needs to be defined for 3 ASTER pixel
  # resolutions (15, 30, 90)

  # Grab LR and UL values
  lr_row <- grep ('LOWERRIGHTM', meta_data)
  ul_row <- grep ('UPPERLEFTM', meta_data)
  lr <- substr (meta_data [lr_row [1]], 15, 50)
  ul <- substr (meta_data [ul_row [1]], 14, 50)
  clip4 <- regexpr (', ' , ul)
  clip5 <- regexpr (', ', lr)

  # Define LR and UL x and y values for 15m VNIR Data
  ul_y <- as.numeric ((substr (ul, 1, (clip4 - 1)))) + 7.5
  ul_x <- as.numeric ((substr (ul, (clip4 + 2), 10000))) - 7.5
  lr_y <- as.numeric ((substr (lr, 1, (clip5 - 1)))) - 7.5
  lr_x <- as.numeric ((substr (lr, (clip5 + 2) , 10000))) + 7.5

  # Define LR and UL x and y values for 30m SWIR Data
  ul_y_30m <- as.numeric ((substr (ul, 1, (clip4 - 1)))) + 15
  ul_x_30m <- as.numeric ((substr (ul, (clip4 + 2), 10000))) - 15
  lr_y_30m <- as.numeric ((substr (lr, 1, (clip5 - 1)))) - 15
  lr_x_30m <- as.numeric ((substr (lr, (clip5 + 2) , 10000))) + 15

  # Define LR and UL x and y values for 90m TIR Data
  ul_y_90m <- as.numeric ((substr (ul, 1, (clip4 - 1)))) + 45
  ul_x_90m <- as.numeric ((substr (ul, (clip4 + 2), 10000))) - 45
  lr_y_90m <- as.numeric ((substr (lr, 1, (clip5 - 1)))) - 45
  lr_x_90m <- as.numeric ((substr (lr, (clip5 + 2) , 10000))) + 45

  # Define UTM zone
  utm_row <- grep ('UTMZONECODE', meta_data)
  utm_zone <- substr (meta_data [utm_row [1]], 1, 50)
  clip6 <- regexpr ('=', utm_zone)
  utm_zone <- substr (utm_zone, clip6 + 1, 50)

  # Configure extent properties (15m VNIR)
  y_min <- min(ul_y, lr_y); y_max <- max(ul_y, lr_y)
  x_max <- max(ul_x, lr_x); x_min <- min(ul_x, lr_x)

  # Configure extent properties (30m SWIR)
  y_min_30m <- min(ul_y_30m, lr_y_30m); y_max_30m <- max(ul_y_30m, lr_y_30m)
  x_max_30m <- max(ul_x_30m, lr_x_30m); x_min_30m <- min(ul_x_30m, lr_x_30m)

  # Configure extent properties (90m TIR)
  y_min_90m <- min(ul_y_90m, lr_y_90m); y_max_90m <- max(ul_y_90m, lr_y_90m)
  x_max_90m <- max(ul_x_90m, lr_x_90m); x_min_90m <- min(ul_x_90m, lr_x_90m)

  raster_dims_15m <- raster::extent(x_min, x_max, y_min, y_max)
  raster_dims_30m <- raster::extent(x_min_30m, x_max_30m, y_min_30m, y_max_30m)
  raster_dims_90m <- raster::extent(x_min_90m, x_max_90m, y_min_90m, y_max_90m)

  # Compile Cordinate Reference System string to attach projection information
  crs_string <- paste('+proj=utm +zone=', utm_zone, ' +datum=WGS84 +units=m
                      +no_defs +ellps=WGS84 +towgs84=0,0,0', sep = '')

  # Remove unneccessary variables
  rm(clip4, clip5, clip6, lr, lr_x, lr_y, meta_data, ul, ul_x, ul_y, utm_zone,
     x_min, x_max, y_min, y_max, lr_row, ul_row, utm_row)

  # Get a list of sds names
  sds <- get_subdatasets(path_file)

  # Limit loop to SDS that contain VNIR/SWIR/TIR data (14 max)
  match_vnir <- grep('VNIR_Swath', sds)
  match_swir <- grep('SWIR_Swath', sds)
  match_tir <- grep('TIR_Swath', sds)

  #-----------------------------------------------------------------------------
  if (length(match_vnir)> 0) {
    for (k in min(match_vnir):max(match_vnir)){
      # Isolate the name of the first sds
      sub_dataset<- sds[k]

      # Get the name of the specific SDS
      clip2 <- max(unlist((gregexpr(':', sub_dataset))))

      # Generate output name for tif
      new_file_name <- strsplit(file_name, '.hdf')
      tif_name <- paste(out_dir, substr(sub_dataset,
                                                            (clip2 + 1), 10000),'.tif', sep='')
      sd_name <- paste(new_file_name, substr(sub_dataset, (clip2 + 1), 10000),
                       sep = '_')
      sub_name <- paste(new_file_name, 'ImageData', sep = '_')
      ast_band_name <- gsub(sub_name, '', sd_name)

      # Extract specified SDS and export as Geotiff
      gdal_translate(path_file, tif_name, sd_index=k, output_Raster = FALSE)

      # Open geotiff and add projection (CRS)
      aster_file <- raster(tif_name, crs = crs_string)
      raster::extent(aster_file) <- raster_dims_15m

      # Convert to large raster layer
      aster_file <- raster::calc(aster_file, fun =function(x){x} )

      # Cropping according to the user defined extent
      if (crop == "u" & typeof(ext) == "double")
      {
        plot (aster_file)
        #plot (aster_file, col = grey (0:255/255))
        print("Please define your extent from the map in plot preview for further processing")
        print("You can click on the top left of custom subset region followed by the bottom right")
        ext <- drawExtent()
      }

      if (crop == "y" || crop == "f" || crop == "u")
      {
        aster_file <- raster::crop (aster_file, ext)
        if (crop == "f")
        {
          aster_file <- raster::mask (aster_file, shape)
        }
      }

      # Export the raster layer file (Geotiff format) to the output directory
      writeRaster(aster_file, filename = tif_name,  options = 'INTERLEAVE=BAND',
                  format = 'GTiff', datatype = 'INT1U', overwrite = TRUE,
                  NAflag = 0)
      # Remove unneccessary variables
      rm(aster_file, sub_dataset, sd_name, sub_name, tif_name, new_file_name)
    }
  }
  #-----------------------------------------------------------------------------
  if (length(match_swir) > 0) {
    for (k in min(match_swir):max(match_swir)){
      # Isolate the name of the first sds
      sub_dataset<- sds[k]

      # Get the name of the specific SDS
      clip2 <- max(unlist((gregexpr(':', sub_dataset))))

      # Generate output name for tif
      new_file_name <- strsplit(file_name, '.hdf')
      tif_name <- paste(out_dir, substr(sub_dataset,
                                                            (clip2 + 1), 10000),'.tif', sep='')
      sd_name <- paste(new_file_name, substr(sub_dataset, (clip2 + 1), 10000),
                       sep = '_')
      sub_name <- paste(new_file_name, 'ImageData', sep = '_')
      ast_band_name <- gsub(sub_name, '', sd_name)

      # Extract specified SDS and export as Geotiff
      gdal_translate(path_file, tif_name, sd_index=k, output_Raster = FALSE)

      # Open geotiff and add projection (CRS)
      aster_file <- raster(tif_name, crs = crs_string)
      raster::extent(aster_file) <- raster_dims_30m

      # Convert to large raster layer
      aster_file <- raster::calc(aster_file, fun =function(x){x} )

      # Cropping according to the user defined extent
      if (crop == "u" & ext == 0)
      {
        #plot (aster_file, col = grey (0:255/255))
        plot (aster_file)
        print("Please define your extent from the map in plot preview for further processing")
        print("You can click on the top left of custom subset region followed by the bottom right")
        ext <- drawExtent()
      }

      if (crop == "y" || crop == "f" || crop == "u")
      {
        aster_file <- raster::crop (aster_file, ext)
        if (crop == "f")
        {
          aster_file <- raster::mask (aster_file, shape)
        }
      }

      # Export the raster layer file (Geotiff format) to the output directory
      writeRaster(aster_file, filename = tif_name,  options = 'INTERLEAVE=BAND',
                  format = 'GTiff', datatype = 'INT1U', overwrite = TRUE,
                  NAflag = 0)
      # Remove unneccessary variables
      rm(aster_file, sub_dataset, sd_name, sub_name, tif_name, new_file_name)
    }
  }
  #-----------------------------------------------------------------------------
  if (length(match_tir) > 0) {
    for (k in min(match_tir):max(match_tir)){
      # Isolate the name of the first sds
      sub_dataset<- sds[k]

      # Get the name of the specific SDS
      clip2 <- max(unlist((gregexpr(':', sub_dataset))))

      # Generate output name for tif
      new_file_name <- strsplit(file_name, '.hdf')
      tif_name <- paste(out_dir, substr(sub_dataset,
                                                            (clip2 + 1), 10000),'.tif', sep='')
      sd_name <- paste(new_file_name, substr(sub_dataset, (clip2 + 1), 10000),
                       sep = '_')
      sub_name <- paste(new_file_name, 'ImageData', sep = '_')
      ast_band_name <- gsub(sub_name, '', sd_name)

      # Extract specified SDS and export as Geotiff
      gdal_translate(path_file, tif_name, sd_index=k, output_Raster = FALSE)

      # Open geotiff and add projection (CRS)
      aster_file <- raster(tif_name, crs = crs_string)
      raster::extent(aster_file) <- raster_dims_90m

      # Convert to large raster layer
      aster_file <- raster::calc(aster_file, fun =function(x){x} )

      # Cropping according to the user defined extent
      if (crop == "u" & typeof(ext) == "double")
      {
        plot (aster_file)
        #plot (aster_file, col = grey (0:255/255))
        print("Please define your extent from the map in plot preview for further processing")
        print("You can click on the top left of custom subset region followed by the bottom right")
        ext <- drawExtent()
      }

      if (crop == "y" || crop == "f" || crop == "u")
      {
        aster_file <- raster::crop (aster_file, ext)
        if (crop == "f")
        {
          aster_file <- raster::mask (aster_file, shape)
        }
      }

      # Export the raster layer file (Geotiff format) to the output directory
      writeRaster(aster_file, filename = tif_name,  options = 'INTERLEAVE=BAND',
                  format = 'GTiff', datatype = 'INT2U', overwrite = TRUE,
                  NAflag = 0)
      # Remove unneccessary variables
      rm(aster_file, sub_dataset, sd_name, sub_name, tif_name, new_file_name)
    }
  }
  #-----------------------------------------------------------------------------
  # Remove unneccessary variables
  rm(ast_band_name, clip2, crs_string, file_name,sds, lr_x_30m, lr_x_90m,
     lr_y_30m, lr_y_90m, match_swir, match_tir, match_vnir, raster_dims_15m,
     raster_dims_30m, raster_dims_90m, ul_x_30m, ul_x_90m, ul_y_30m, ul_y_90m,
     x_max_30m, x_max_90m, x_min_30m, x_min_90m, y_max_30m, y_max_90m,y_min_30m,
     y_min_90m)
}}
#-------------------------------------------------------------------------------
