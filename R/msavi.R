#' Modified soil adjusted vegetation index
#'
#' Modified Soil Adjusted Vegetation Index (MSAVI) is a vegetation index.
#' Advantage of this index is that, it increases the dynamic range of the vegetation signal while further minimizing the soil background influences,
#' resulting in greater vegetation sensitivity as defined by a 'vegetation signal' to 'soil noise' ratio.
#' @param ext2crop,crop,directory Same as mentioned in \code{\link[ASIP]{arvi}}.
#' @return Computed MSAVI product
#' @note 1. MSAVI=((2r_nir + 1) - ((2r_nir + 1)^2 - 8(r_nir - r_red))^0.5)/2
#'
#' where, "r_" denotes TOA reflectance band.
#'
#' Other important notes are mentioned in \code{\link[ASIP]{custom.eqn}}.
#' @export
#' @references \href{http://www.sciencedirect.com/science/article/pii/0034425794901341}{Qi J, Chehbouni A, Huete A R, Kerr Y, Sorooshian S (1994) A modified soil adjusted vegetation index. Remote Sensing of Environment, 48 (2), pp: 119-126. doi:10.1016/0034-4257(94)90134-1.}
#' @importFrom raster raster writeRaster extent mask crop
#' @importFrom utils tail
#' @examples
#' library (raster)
#' library (rgdal)
#' # Finding the path of the sample satellite image directory.
#' # User may define paths directly like "/home/ur_folder" or "C:/ur_folder"
#' path <- system.file ("TM_sample", package = "ASIP")
#' shapefil <- paste0 (path, "/test.shp")
#' op <- msavi (directory = path, crop = "y", ext2crop = shapefil)
msavi <- function(directory = getwd(), crop = "n", ext2crop = "none")
{
  # If the directory is not set
  bands <- length(list.files(directory,pattern = "*TIF"))
  if (bands == 0)
    stop("Define your satellite image folder path properly")
  # Finding out which satellite sensor data & name of satellite image data
  files <- list.files(directory)
  for (i in 1:length(files))
  {
    file <- files[i]
    broke_name <- strsplit(file, "_B1.TI")
    broke_name <- broke_name[[1]]
    if (utils::tail(broke_name,1) == "F")
    {
      sat_fold <- broke_name[1]
      satellite <- substr(sat_fold,1,2)
      break()
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
      shape <- raster::shapefile (ext2crop)
      ext <- raster::extent (shape)
    }
    if (typeof (ext2crop) == "S4")
    {
      ext <- raster::extent (ext2crop)
      shape <- ext2crop
    }
  }

  meta_data <- readLines(paste0(directory,"/",sat_fold,"_MTL.txt"))
  count_i <- length(meta_data)
  if (count_i==0){print("ERROR: MTL file not found")}
  ######### Landsat 8 starting###############
  if (satellite=="LC")
  {
    if (crop == "u")
    {
      b5 <- raster (paste0 (directory, "/", sat_fold, "_B5.TIF"))
      b4 <- raster (paste0 (directory, "/", sat_fold, "_B4.TIF"))
      b3 <- raster (paste0 (directory, "/", sat_fold, "_B3.TIF"))
      stak <- raster::stack(c(b5,b4,b3))
      plotRGB(stak, scale = 65536)
      print("Please define your extent from the map in plot preview for further processing")
      print("You can click on the top left of custom subset region followed by the bottom right")
      ext <- drawExtent()
    }
    # Extracting values from meta data
    for (i in 1:count_i)
    {
      line <- meta_data[i]
      line_splited <- strsplit(line," ")
      words <- line_splited[[1]]
      counts <- length(words)
      for (j in 1:counts)
      {
        if (words[j]=="REFLECTANCE_ADD_BAND_5"){ nir_refl_add <- as.double(words[j+2])}
        if (words[j]=="REFLECTANCE_MULT_BAND_5"){ nir_refl_mult <- as.double(words[j+2])}
        if (words[j]=="REFLECTANCE_ADD_BAND_4"){ red_refl_add <- as.double(words[j+2])}
        if (words[j]=="REFLECTANCE_MULT_BAND_4"){ red_refl_mult <- as.double(words[j+2])}

        if (words[j]=="DATE_ACQUIRED"){ data_aq <- as.character(words[j+2])}
        if (words[j]=="SUN_ELEVATION"){ sun_ele <- as.double(words[j+2])}
      }
    }

    nir <- as.integer(raster(paste0(directory,"/",sat_fold,"_B5.TIF")))
    if (crop == "y" || crop == "f" || crop == "u")
    {
      nir <- crop(nir, ext)
      if (crop=="f")
      {
        nir <- mask(nir,shape)
      }
    }
    toa_nir <- ((nir * nir_refl_mult) + nir_refl_add)/sin(sun_ele*(pi/180))

    red <- as.integer(raster(paste0(directory,"/",sat_fold,"_B4.TIF")))
    if (crop == "y" || crop == "f" || crop == "u")
    {
      red <- crop(red, ext)
      if (crop=="f")
      {
        red <- mask(red,shape)
      }
    }
    toa_red <- ((red * red_refl_mult) + red_refl_add)/sin(sun_ele*(pi/180))
  }
  ########### Landsat-8 ending ##############
  ########### Landsat-7 starting ##############
  if (satellite=="LE")
  {
    if (crop == "u")
    {
      b4 <- raster (paste0 (directory, "/", sat_fold, "_B4.TIF"))
      b3 <- raster (paste0 (directory, "/", sat_fold, "_B3.TIF"))
      b2 <- raster (paste0 (directory, "/", sat_fold, "_B2.TIF"))
      stak <- stack(c (b4, b3, b2))
      plotRGB(stak)
      print("Please define your extent from the map in plot preview for further processing")
      print("You can click on the top left of custom subset region followed by the bottom right")
      ext <- drawExtent()
    }
    qcal_max <- 255
    d <- 0
    for (i in 1:count_i)
    {
      line <- meta_data[i]
      line_splited <- strsplit(line," ")
      words <- line_splited[[1]]
      counts <- length(words)
      for (j in 1:counts)
      {
        if (words[j]=="QUANTIZE_CAL_MIN_BAND_1"){ qcal_min <- as.double(words[j+2])}
        if (words[j]=="EARTH_SUN_DISTANCE"){ d <- as.double(words[j+2])}
        if (words[j]=="DATE_ACQUIRED"){ data_aq <- as.character(words[j+2])}
        if (words[j]=="SUN_ELEVATION"){ sun_ele <- as.double(words[j+2])}

        if (words[j]=="RADIANCE_MAXIMUM_BAND_3"){ lmax3 <- as.double(words[j+2])}
        if (words[j]=="RADIANCE_MAXIMUM_BAND_4"){ lmax4 <- as.double(words[j+2])}

        if (words[j]=="RADIANCE_MINIMUM_BAND_3"){ lmin3 <- as.double(words[j+2])}
        if (words[j]=="RADIANCE_MINIMUM_BAND_4"){ lmin4 <- as.double(words[j+2])}
      }
    }

    if (d==0)
    {
      dat_tex <- as.Date(data_aq)
      jul_ful <- julian(dat_tex)
      yr_rmv_num <- jul_ful%/%365.25
      jul_day <- jul_ful-(yr_rmv_num*365.25)+1.5
      jul_flot <- jul_day%%365.25
      jul_final <- as.integer(jul_flot)
      d <- 1 + (0.0167 * sin ((pi/180) *2 *pi *(jul_final- 93.5) / 365))
    }

    red <- raster(paste0(directory,"/",sat_fold,"_B3.TIF"))
    if (crop == "y" || crop == "f" || crop == "u")
    {
      red <- crop(red,ext)
      if (crop=="f")
      {
        red <- mask(red,shape)
      }
    }
    rad_b3 <- ((lmax3-lmin3)/(qcal_max-qcal_min)) * (red-qcal_min) + lmin3
    toa_red <- pi * rad_b3 * d^2  / 1547 * sin(sun_ele*(pi/180))

    nir <- raster(paste0(directory,"/",sat_fold,"_B4.TIF"))
    if (crop == "y" || crop == "f" || crop == "u")
    {
      nir <- crop(nir,ext)
      if (crop=="f")
      {
        nir <- mask(nir,shape)
      }
    }
    rad_b4 <- ((lmax4-lmin4)/(qcal_max-qcal_min)) * (nir-qcal_min) + lmin4
    toa_nir <- pi * rad_b4 * d^2  / 1044 * sin(sun_ele*(pi/180))
  }
  ############## Landsat ETM ending ##################
  ############## Landsat TM starting ##################
  if (satellite=="LT")
  {
    if (crop == "u")
    {
      b4 <- raster (paste0 (directory, "/", sat_fold, "_B4.TIF"))
      b3 <- raster (paste0 (directory, "/", sat_fold, "_B3.TIF"))
      b2 <- raster (paste0 (directory, "/", sat_fold, "_B2.TIF"))
      stak <- stack(c (b4, b3, b2))
      plotRGB(stak)
      print("Please define your extent from the map in plot preview for further processing")
      print("You can click on the top left of custom subset region followed by the bottom right")
      ext <- drawExtent()
    }
    qcal_max <- 255
    d <- 0
    for (i in 1:count_i)
    {
      line <- meta_data[i]
      line_splited <- strsplit(line," ")
      words <- line_splited[[1]]
      counts <- length(words)
      for (j in 1:counts)
      {
        if (words[j]=="QUANTIZE_CAL_MIN_BAND_1"){ qcal_min <- as.double(words[j+2])}
        if (words[j]=="DATE_ACQUIRED"){ data_aq <- as.character(words[j+2])}
        if (words[j]=="SUN_ELEVATION"){ sun_ele <- as.double(words[j+2])}
        if (words[j]=="EARTH_SUN_DISTANCE"){ d <- as.double(words[j+2])}
        if (words[j]=="SPACECRAFT_ID"){ tm_id <- as.character(words[j+2])}

        if (words[j]=="RADIANCE_MAXIMUM_BAND_3"){ lmax3 <- as.double(words[j+2])}
        if (words[j]=="RADIANCE_MAXIMUM_BAND_4"){ lmax4 <- as.double(words[j+2])}

        if (words[j]=="RADIANCE_MINIMUM_BAND_3"){ lmin3 <- as.double(words[j+2])}
        if (words[j]=="RADIANCE_MINIMUM_BAND_4"){ lmin4 <- as.double(words[j+2])}
      }
    }

    if (tm_id=="\"LANDSAT_5\"")
    {
      esun3 <- 1551
      esun4 <- 1036
    }

    if (tm_id!="\"LANDSAT_5\"")
    {
      esun3 <- 1554
      esun4 <- 1033
    }

    if (d==0)
    {
      dat_tex <- as.Date(data_aq)
      jul_ful <- julian(dat_tex)
      yr_rmv_num <- jul_ful%/%365.25
      jul_day <- jul_ful-(yr_rmv_num*365.25)+1.5
      jul_flot <- jul_day%%365.25
      jul_final <- as.integer(jul_flot)
      d <- 1 + (0.0167 * sin ((pi/180) *2 *pi *(jul_final- 93.5) / 365))
    }

    red <- raster(paste0(directory,"/",sat_fold,"_B3.TIF"))
    if (crop == "y" || crop == "f" || crop == "u")
    {
      red <- crop(red,ext)
      if (crop=="f")
      {
        red <- mask(red,shape)
      }
    }
    rad_b3 <- ((lmax3-lmin3)/(qcal_max-qcal_min)) * (red-qcal_min) + lmin3
    toa_red <- pi * rad_b3 * d^2  / esun3 * sin(sun_ele*(pi/180))

    nir <- raster(paste0(directory,"/",sat_fold,"_B4.TIF"))
    if (crop == "y" || crop == "f" || crop == "u")
    {
      nir <- crop(nir,ext)
      if (crop=="f")
      {
        nir <- mask(nir,shape)
      }
    }
    rad_b4 <- ((lmax4-lmin4)/(qcal_max-qcal_min)) * (nir-qcal_min) + lmin4
    toa_nir <- pi * rad_b4 * d^2  / esun4 * sin(sun_ele*(pi/180))
  }
  ######## Landsat TM ending ############

  msavi_c1= (2*toa_nir)+1
  msavi_c2= ((2*toa_nir)+1)^2
  msavi_c3= (msavi_c2 - (8 * (toa_nir-toa_red)))^0.5
  msavi <- (msavi_c1 - msavi_c3)/2
  #writeRaster(msavi,paste0(directory,"/","msavi_",data_aq),format="GTiff", overwrite=TRUE)
  return(msavi)
  cat ("\nProgram completed, output is produced as a variable named 'msavi'")
}
