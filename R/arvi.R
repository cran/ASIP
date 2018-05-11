#' Atmospherically Resistant Vegetation Index
#'
#' Atmospherically Resistant Vegetation Index (ARVI) is a vegetation based index which
#' minimizes the effects of atmospheric scattering in comparison to NDVI.
#'
#' @param directory Path to Satellite image folder. Assign as string (inside double quotes).
#' Either assing inside function or set up satellite image folder as the current working directory
#' before running the function.
#' To define current working directory, either use shortcut key Ctrl+Shift+H  or use \code{\link{setwd}} funtion.
#' @param crop Defines the method of cropping outputs to custom extent.
#'
#' "n" <- No cropping required (Default).
#'
#' "u" <- Satellite image will be plotted in the plot window and user can choose the extent by clicking on the top left maximum followed by bottom right maximum.
#'
#' "y" <- Crop to the maximum and minimum extent of the shapefile.
#'
#' "f" <- Crop to exact shapefile boundary.
#' @param ext2crop Path to the shapefile (*.shp) which will be used for cropping. Shapefile should have SAME CORDINATE SYSTEM as the satellite image.
#' Either provide the full path of .shp file or provide the name of the shapefile variable which is already opened.
#' @param gamma It is an aerosol dependant factor. For more details please refer Kaufman and Tanre (1992). By default the value is 1.
#' @return Computed ARVI product
#' @note 1. ARVI = (r_nir - rb)/(r_nir + rb), where
#'
#' rb = r_red - gamma (r_blue - r_red)'  and "r_" denotes Top Of Atmoshpere (TOA) reflection, 'gamma' value is 1 by default as recommended if information about the aerosol type is not available.
#' Please refer Kaufman and Tanre (1992) for more details.
#'
#' Other important notes are mentioned in \code{\link[ASIP]{custom.eqn}}.
#' @export
#' @importFrom raster raster writeRaster extent mask crop plotRGB drawExtent
#' @importFrom utils tail
#' @importFrom rgdal readOGR
#' @references \href{http://ieeexplore.ieee.org/document/134076/?arnumber=134076&tag=1}{Kaufman, Y. J. and D. Tanre (1992) Atmospherically resistant vegetation index (ARVI) for EOS-MODIS, IEEE Transactions on Geoscience and Remote Sensing, 30 (2). doi:10.1109/36.134076.}
#' @examples
#' library (raster)
#' library (rgdal)
#' # Finding the path of the sample satellite image directory.
#' # User may define paths directly like "/home/ur_folder" or "C:/ur_folder"
#' path <- system.file ("TM_sample", package = "ASIP")
#' shapefil <- paste0 (path, "/test.shp")
#' op <- arvi (directory = path, crop = "y", ext2crop = shapefil)
arvi <- function (directory = getwd(), crop = "n", ext2crop = "none", gamma = 1)
{
  # If the directory is not set
  bands <- length (list.files (directory, pattern = "*TIF"))
  if (bands == 0)
  stop("Define your satellite image folder path properly")
  # Finding out which satellite sensor data & name of satellite image data
  files <- list.files (directory)
  for (i in 1: length (files))
  {
    file <- files [i]
    broke_name <- strsplit (file, "_B1.TI")
    broke_name <- broke_name [[1]]
    if (utils::tail (broke_name, 1) == "F")
    {
      sat_fold <- broke_name [1]
      satellite <- substr (sat_fold, 1, 2)
      break ()
    }
  }

  # Defining the crop extent
  if (crop != "n" && crop != "y" && crop !="u" && crop != "f")
    stop ("Define argument 'crop' properly. Use either n, y, f or u in double quotes. Type ?arvi in console to read more about the function")
  if (crop != "n" && ext2crop == "none")
    {
    if (crop != "u")
    stop ("Define argument 'ext2crop' properly if croppping is required, otherwise choose argument 'crop' as n in double quotes")
    }
  if (crop == "y" || crop == "f")
  {
    if (typeof (ext2crop) == "character")
    {
      shape <- rgdal::readOGR (ext2crop)
      ext <- raster::extent (shape)
    }
    if (typeof (ext2crop) == "S4")
    {
      ext <- raster::extent (ext2crop)
      shape <- ext2crop
    }
  }

  meta_data <- readLines (paste0 (directory, "/", sat_fold, "_MTL.txt"))
  count_i <- length (meta_data)
  if (count_i == 0) { print ("ERROR: MTL file not found") }
  ######### Landsat 8 starting###############
  if (satellite == "LC")
  {
    if (crop == "u")
    {
      b5 <- raster (paste0 (directory, "/", sat_fold, "_B5.TIF"))
      b4 <- raster (paste0 (directory, "/", sat_fold, "_B4.TIF"))
      b3 <- raster (paste0 (directory, "/", sat_fold, "_B3.TIF"))
      stak <- raster::stack (c (b5, b4, b3))
      plotRGB (stak, scale = 65536)
      print ("Please define your extent from the map in plot preview for further processing")
      print ("You can click on the top left of custom subset region followed by the bottom right")
      ext <- drawExtent()
    }

    # Extracting values from meta data
    for (i in 1: count_i)
    {
      line <- meta_data [i]
      line_splited <- strsplit (line, " ")
      words <- line_splited [[1]]
      counts <- length (words)
      for (j in 1: counts)
      {
        if (words [j] == "REFLECTANCE_ADD_BAND_5") { nir_refl_add <- as.double (words [j+2])}
        if (words [j] == "REFLECTANCE_MULT_BAND_5") { nir_refl_mult <- as.double (words [j+2])}
        if (words [j] == "REFLECTANCE_ADD_BAND_4") { red_refl_add <- as.double (words [j+2])}
        if (words [j] == "REFLECTANCE_MULT_BAND_4") { red_refl_mult <- as.double (words [j+2])}
        if (words [j] == "REFLECTANCE_ADD_BAND_2"){ blu_refl_add <- as.double (words [j+2])}
        if (words [j] == "REFLECTANCE_MULT_BAND_2"){ blu_refl_mult <- as.double (words [j+2])}

        if (words [j] == "DATE_ACQUIRED"){ data_aq <- as.character (words [j+2])}
        if (words [j] == "SUN_ELEVATION"){ sun_ele <- as.double (words [j+2])}
      }
    }

    nir <- as.integer (raster::raster (paste0 (directory, "/", sat_fold, "_B5.TIF")))
    if (is.null(crop) == FALSE && (crop == "y" || crop == "f" || crop == "u"))
    {
      nir <- raster::crop (nir, ext)
      if (crop == "f")
      {
        nir <- raster::mask (nir, shape)
      }
    }
    toa_nir <- ((nir * nir_refl_mult) + nir_refl_add) /sin (sun_ele* (pi/180))

    red <- as.integer (raster (paste0 (directory, "/", sat_fold, "_B4.TIF")))
    if (crop == "y" || crop == "f" || crop == "u")
    {
      red <- raster::crop (red, ext)
      if (crop == "f")
      {
        red <- raster::mask (red, shape)
      }
    }
    toa_red <- ((red * red_refl_mult) + red_refl_add) / sin (sun_ele* (pi/180))

    blu <- as.integer (raster (paste0 (directory, "/", sat_fold, "_B2.TIF")))
    if (crop == "y" || crop == "f" || crop == "u")
    {
      blu <- raster::crop (blu, ext)
      if (crop == "f")
      {
        blu <- raster::mask (blu, shape)
      }
    }
    toa_blu <- ((blu * blu_refl_mult) + blu_refl_add) / sin (sun_ele* (pi/180))
  }
  ########### Landsat-8 ending ##############
  ########### Landsat-7 starting ##############
  if (satellite == "LE")
  {
    if (crop == "u")
    {
      b4 <- raster (paste0 (directory, "/", sat_fold, "_B4.TIF"))
      b3 <- raster (paste0 (directory, "/", sat_fold, "_B3.TIF"))
      b2 <- raster (paste0 (directory, "/", sat_fold, "_B2.TIF"))
      stak <- stack (c (b4, b3, b2))
      plotRGB (stak)
      print ("Please define your extent from the map in plot preview for further processing")
      print ("You can click on the top left of custom subset region followed by the bottom right")
      ext <- drawExtent()
    }
    qcal_max <- 255
    d <- 0
    for (i in 1: count_i)
    {
      line <- meta_data [i]
      line_splited <- strsplit (line, " ")
      words <- line_splited [[1]]
      counts <- length (words)
      for (j in 1: counts)
      {
        if (words [j] == "QUANTIZE_CAL_MIN_BAND_1") { qcal_min <- as.double (words [j+2])}
        if (words [j] == "EARTH_SUN_DISTANCE") { d <- as.double (words [j+2])}
        if (words [j] == "DATE_ACQUIRED"){ data_aq <- as.character (words [j+2])}
        if (words [j] == "SUN_ELEVATION"){ sun_ele <- as.double (words [j+2])}

        if (words [j] == "RADIANCE_MAXIMUM_BAND_1") { lmax1 <- as.double (words [j+2])}
        if (words [j] == "RADIANCE_MAXIMUM_BAND_3") { lmax3 <- as.double (words [j+2])}
        if (words [j] == "RADIANCE_MAXIMUM_BAND_4"){ lmax4 <- as.double (words [j+2])}

        if (words [j] == "RADIANCE_MINIMUM_BAND_1"){ lmin1 <- as.double (words [j+2])}
        if (words [j] == "RADIANCE_MINIMUM_BAND_3"){ lmin3 <- as.double (words [j+2])}
        if (words [j] == "RADIANCE_MINIMUM_BAND_4"){ lmin4 <- as.double (words [j+2])}
      }
    }

    if (d == 0)
    {
      dat_tex <- as.Date (data_aq)
      jul_ful <- julian (dat_tex)
      yr_rmv_num <- jul_ful%/%365.25
      jul_day <- jul_ful - (yr_rmv_num * 365.25)+ 1.5
      jul_flot <- jul_day%%365.25
      jul_final <- as.integer (jul_flot)
      d <- 1 + (0.0167 * sin ((pi/180) *2 *pi *(jul_final- 93.5) / 365))
    }

    red <- raster (paste0 (directory, "/", sat_fold, "_B3.TIF"))
    if (crop == "y" || crop == "f" || crop == "u")
    {
      red <- raster::crop (red, ext)
      if (crop == "f")
      {
        red <- raster::mask (red, shape)
      }
    }
    rad_b3 <- ((lmax3 -lmin3) / (qcal_max- qcal_min)) * (red- qcal_min) + lmin3
    toa_red <- pi * rad_b3 * d^2  / 1547 * sin (sun_ele* (pi/180))

    nir <- raster (paste0 (directory, "/", sat_fold, "_B4.TIF"))
    if (crop == "y" || crop == "f" || crop == "u")
    {
      nir <- raster::crop (nir, ext)
      if (crop == "f")
      {
        nir <- raster::mask (nir, shape)
      }
    }
    rad_b4 <- ((lmax4- lmin4) / (qcal_max- qcal_min)) * (nir- qcal_min) + lmin4
    toa_nir <- pi * rad_b4 * d^2  / 1044 * sin (sun_ele* (pi/180))

    blu <- raster (paste0 (directory, "/", sat_fold, "_B1.TIF"))
    if (crop == "y" || crop == "f" || crop == "u")
    {
      blu <- raster::crop (blu, ext)
      if (crop == "f")
      {
        blu <- raster::mask (blu, shape)
      }
    }
    rad_b1 <- ((lmax1- lmin1) / (qcal_max- qcal_min)) * (blu- qcal_min) + lmin1
    toa_blu <- pi * rad_b1 * d^2  / 1970 * sin (sun_ele* (pi/180))
  }
  ############## Landsat ETM ending ##################
  ############## Landsat TM starting ##################
  if (satellite == "LT")
  {
    if (crop == "u")
    {
      b4 <- raster (paste0 (directory, "/", sat_fold, "_B4.TIF"))
      b3 <- raster (paste0 (directory, "/", sat_fold, "_B3.TIF"))
      b2 <- raster (paste0 (directory, "/", sat_fold, "_B2.TIF"))
      stak <- stack (c (b4, b3, b2))
      plotRGB (stak)
      print ("Please define your extent from the map in plot preview for further processing")
      print ("You can click on the top left of custom subset region followed by the bottom right")
      ext <- drawExtent()
    }
    qcal_max <- 255
    d <- 0
    for (i in 1: count_i)
    {
      line <- meta_data [i]
      line_splited <- strsplit (line," ")
      words <- line_splited [[1]]
      counts <- length (words)
      for (j in 1: counts)
      {
        if (words [j] == "QUANTIZE_CAL_MIN_BAND_1") { qcal_min <- as.double (words [j+2])}
        if (words [j] == "DATE_ACQUIRED") { data_aq <- as.character (words [j+2])}
        if (words [j] == "SUN_ELEVATION"){ sun_ele <- as.double (words [j+2])}
        if (words [j] == "EARTH_SUN_DISTANCE") { d <- as.double (words [j+2])}
        if (words [j] == "SPACECRAFT_ID") { tm_id <- as.character (words [j+2])}

        if (words [j] == "RADIANCE_MAXIMUM_BAND_1"){ lmax1 <- as.double (words [j+2])}
        if (words [j] == "RADIANCE_MAXIMUM_BAND_3"){ lmax3 <- as.double (words [j+2])}
        if (words [j] == "RADIANCE_MAXIMUM_BAND_4"){ lmax4 <- as.double (words [j+2])}

        if (words [j] == "RADIANCE_MINIMUM_BAND_1"){ lmin1 <- as.double (words [j+2])}
        if (words [j] == "RADIANCE_MINIMUM_BAND_3"){ lmin3 <- as.double (words [j+2])}
        if (words [j] == "RADIANCE_MINIMUM_BAND_4"){ lmin4 <- as.double (words [j+2])}
      }
    }


    if (tm_id == "\"LANDSAT_5\"")
    {
      esun3 <- 1551
      esun4 <- 1036
    }

    if (tm_id != "\"LANDSAT_5\"")
    {
      esun3 <- 1554
      esun4 <- 1033
    }

    if (d == 0)
    {
      dat_tex <- as.Date (data_aq)
      jul_ful <- julian (dat_tex)
      yr_rmv_num <- jul_ful%/%365.25
      jul_day <- jul_ful -(yr_rmv_num *365.25) + 1.5
      jul_flot <- jul_day%%365.25
      jul_final <- as.integer (jul_flot)
      d <- 1 + (0.0167 * sin ((pi/180) *2 *pi *(jul_final- 93.5) / 365))
    }

    red <- raster (paste0 (directory, "/", sat_fold, "_B3.TIF"))
    if (crop == "y" || crop == "f" || crop == "u")
    {
      red <- raster::crop (red, ext)
      if (crop == "f")
      {
        red <- raster::mask (red, shape)
      }
    }
    rad_b3 <- ((lmax3- lmin3) / (qcal_max- qcal_min)) * (red- qcal_min) + lmin3
    toa_red <- pi * rad_b3 * d^2  / esun3 * sin (sun_ele *(pi/180))

    nir <- raster (paste0 (directory, "/", sat_fold, "_B4.TIF"))
    if (crop == "y" || crop == "f" || crop == "u")
    {
      nir <- raster::crop (nir, ext)
      if (crop == "f")
      {
        nir <- raster::mask (nir, shape)
      }
    }
    rad_b4 <- ((lmax4- lmin4) / (qcal_max- qcal_min)) * (nir- qcal_min) + lmin4
    toa_nir <- pi * rad_b4 * d^2  / esun4 * sin (sun_ele* (pi/180))

    blu <- raster (paste0 (directory, "/", sat_fold, "_B1.TIF"))
    if (crop == "y" || crop == "f" || crop == "u")
    {
      blu <- raster::crop (blu, ext)
      if (crop == "f")
      {
        blu <- raster::mask (blu, shape)
      }
    }
    rad_b1 <- ((lmax1- lmin1) / (qcal_max- qcal_min)) * (blu- qcal_min) + lmin1
    toa_blu <- pi * rad_b1 * d^2  / 1958 * sin (sun_ele* (pi/180))
  }
  ######## Landsat TM ending ############
  rb <- toa_red - gamma * (toa_blu - toa_red)
  arvi <- (toa_nir - rb) / (toa_nir + rb)
  return(arvi)
  #raster::writeRaster (arvi, paste0 (directory, "/", "arvi_", data_aq), format = "GTiff", overwrite = TRUE)
  cat ("\nProgram completed, output is produced as a variable named 'arvi'")
}
