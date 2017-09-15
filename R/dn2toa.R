#' DN to TOA conversion of optical bands
#'
#' Digital number (DN) bands to Top of Atmosphere (TOA) conversion.
#' @param ext2crop,crop,directory Same as mentioned in \code{\link[ASIP]{arvi}}.
#' @param b1 By default Band1 will be processed to TOA reflectance. To cancel production of this band assign value 0.
#' @param b2 By default Band2 will be processed to TOA reflectance. To cancel production of this band assign value 0.
#' @param b3 By default Band3 will be processed to TOA reflectance. To cancel production of this band assign value 0.
#' @param b4 By default Band4 will be processed to TOA reflectance. To cancel production of this band assign value 0.
#' @param b5 By default Band5 will be processed to TOA reflectance. To cancel production of this band assign value 0.
#' @param b6 By default Band6 will be processed to TOA reflectance. To cancel production of this band assign value 0.
#' @param b7 By default Band7 will be processed to TOA reflectance. To cancel production of this band assign value 0.
#' @return Each bands selected will produce corresponding image in *.tif format in the input directory.
#' @note 1. This function followed by \code{\link[ASIP]{multi.indices}} is recommended only if user is intended to produce multiple indices
#' like ndvi & gemi other than running seperate function for each product to save processing time and resources.
#'
#' Other important notes are mentioned in \code{\link[ASIP]{custom.eqn}}.
#' @export
#' @importFrom raster raster writeRaster extent mask crop
#' @importFrom utils tail
#' @references \href{https://landsat.usgs.gov/sites/default/files/documents/Landsat8DataUsersHandbook.pdf}{USGS (2016) Landsat 8 (L8) data users handbook, version 2.}
#'
#' Landsat 7 science data users handbook, NASA. Available at "https://landsat.gsfc.nasa.gov/wp-content/uploads/2016/08/Landsat7_Handbook.pdf".
#' @examples
#' library (raster)
#' library (rgdal)
#' # Finding the path of the sample satellite image directory.
#' # User may define paths directly like "/home/ur_folder" or "C:/ur_folder"
#' path <- system.file ("TM_sample", package = "ASIP")
#' shapefil <- paste0 (path, "/test.shp")
#' # Assign 0 values to band names which are not required
#' dn2toa (path, crop = "f", ext2crop = shapefil, b3=0, b4=0, b5=0, b6 = 0, b7 = 0)

# TOA
dn2toa <- function (directory= getwd(), crop = "n", ext2crop = "none",b1=1,b2=1,b3=1,b4=1,b5=1,b6=1,b7=1)
{
  # If the directory is not set
  bands <- length (list.files(directory,pattern = "*TIF"))
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
        if (words[j]=="REFLECTANCE_ADD_BAND_7"){ swir2_refl_add <- as.double(words[j+2])}
        if (words[j]=="REFLECTANCE_MULT_BAND_7"){ swir2_refl_mult <- as.double(words[j+2])}
        if (words[j]=="REFLECTANCE_ADD_BAND_6"){ swir1_refl_add <- as.double(words[j+2])}
        if (words[j]=="REFLECTANCE_MULT_BAND_6"){ swir1_refl_mult <- as.double(words[j+2])}
        if (words[j]=="REFLECTANCE_ADD_BAND_5"){ nir_refl_add <- as.double(words[j+2])}
        if (words[j]=="REFLECTANCE_MULT_BAND_5"){ nir_refl_mult <- as.double(words[j+2])}
        if (words[j]=="REFLECTANCE_ADD_BAND_4"){ red_refl_add <- as.double(words[j+2])}
        if (words[j]=="REFLECTANCE_MULT_BAND_4"){ red_refl_mult <- as.double(words[j+2])}
        if (words[j]=="REFLECTANCE_ADD_BAND_3"){ green_refl_add <- as.double(words[j+2])}
        if (words[j]=="REFLECTANCE_MULT_BAND_3"){ green_refl_mult <- as.double(words[j+2])}
        if (words[j]=="REFLECTANCE_ADD_BAND_2"){ blu_refl_add <- as.double(words[j+2])}
        if (words[j]=="REFLECTANCE_MULT_BAND_2"){ blu_refl_mult <- as.double(words[j+2])}
        if (words[j]=="REFLECTANCE_ADD_BAND_1"){ aero_refl_add <- as.double(words[j+2])}
        if (words[j]=="REFLECTANCE_MULT_BAND_1"){ aero_refl_mult <- as.double(words[j+2])}
        if (words[j]=="DATE_ACQUIRED"){ data_aq <- as.character(words[j+2])}
        if (words[j]=="SUN_ELEVATION"){ sun_ele <- as.double(words[j+2])}
      }
    }
    # Defining bands & toa calculation

    if (b7==1)
    {
      swir2 <- as.integer(raster(paste0(directory,"/",sat_fold,"_B7.TIF")))
      if (crop == "y" || crop == "f" || crop == "u")
      {
        swir2 <- crop(swir2, ext)
        if (crop=="f")
        {
          swir2 <- mask(swir2,shape)
        }
      }
      toa_swir2 <- ((swir2 * swir2_refl_mult) + swir2_refl_add)/sin(sun_ele*(pi/180))
      writeRaster(toa_swir2,paste0(directory,"/","toa_swir2"),format="GTiff",overwrite=TRUE)
    }
    if (b6==1)
    {
      swir1 <- as.integer(raster(paste0(directory,"/",sat_fold,"_B6.TIF")))
      if (crop == "y" || crop == "f" || crop == "u")
      {
        swir1 <- crop(swir1, ext)
        if (crop=="f")
        {
          swir1 <- mask(swir1,shape)
        }
      }
      toa_swir1 <- ((swir1 * swir1_refl_mult) + swir1_refl_add)/sin(sun_ele*(pi/180))
      writeRaster(toa_swir1,paste0(directory,"/","toa_swir1"),format="GTiff",overwrite=TRUE)
    }
    if (b5==1)
    {
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
      writeRaster(toa_nir,paste0(directory,"/","toa_nir"),format="GTiff",overwrite=TRUE)
    }
    if (b4==1)
    {
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
      writeRaster(toa_red,paste0(directory,"/","toa_red"),format="GTiff",overwrite=TRUE)
    }
    if (b3==1)
    {
      green <- as.integer(raster(paste0(directory,"/",sat_fold,"_B3.TIF")))
      if (crop == "y" || crop == "f" || crop == "u")
      {
        green <- crop(green, ext)
        if (crop=="f")
        {
          green <- mask(green,shape)
        }
      }
      toa_green <- ((green * green_refl_mult) + green_refl_add)/sin(sun_ele*(pi/180))
      writeRaster(toa_green,paste0(directory,"/","toa_green"),format="GTiff",overwrite=TRUE)
    }
    if (b2==1)
    {
      blu <- as.integer(raster(paste0(directory,"/",sat_fold,"_B2.TIF")))
      if (crop == "y" || crop == "f" || crop == "u")
      {
        blu <- crop(blu, ext)
        if (crop=="f")
        {
          blu <- mask(blu,shape)
        }
      }
      toa_blu <- ((blu * blu_refl_mult) + blu_refl_add)/sin(sun_ele*(pi/180))
      writeRaster(toa_blu,paste0(directory,"/","toa_blu"),format="GTiff",overwrite=TRUE)
    }
    if (b1==1)
    {
      aero <- as.integer(raster(paste0(directory,"/",sat_fold,"_B1.TIF")))
      if (crop == "y" || crop == "f" || crop == "u")
      {
        aero <- crop(aero, ext)
        if (crop=="f")
        {
          aero <- mask(aero,shape)
        }
      }
      toa_aero <- ((aero * aero_refl_mult) + aero_refl_add)/sin(sun_ele*(pi/180))
      writeRaster(toa_aero,paste0(directory,"/","toa_aero"),format="GTiff",overwrite=TRUE)
    }
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

        if (words[j]=="RADIANCE_MAXIMUM_BAND_1"){ lmax1 <- as.double(words[j+2])}
        if (words[j]=="RADIANCE_MAXIMUM_BAND_2"){ lmax2 <- as.double(words[j+2])}
        if (words[j]=="RADIANCE_MAXIMUM_BAND_3"){ lmax3 <- as.double(words[j+2])}
        if (words[j]=="RADIANCE_MAXIMUM_BAND_4"){ lmax4 <- as.double(words[j+2])}
        if (words[j]=="RADIANCE_MAXIMUM_BAND_5"){ lmax5 <- as.double(words[j+2])}
        if (words[j]=="RADIANCE_MAXIMUM_BAND_7"){ lmax7 <- as.double(words[j+2])}

        if (words[j]=="RADIANCE_MINIMUM_BAND_1"){ lmin1 <- as.double(words[j+2])}
        if (words[j]=="RADIANCE_MINIMUM_BAND_2"){ lmin2 <- as.double(words[j+2])}
        if (words[j]=="RADIANCE_MINIMUM_BAND_3"){ lmin3 <- as.double(words[j+2])}
        if (words[j]=="RADIANCE_MINIMUM_BAND_4"){ lmin4 <- as.double(words[j+2])}
        if (words[j]=="RADIANCE_MINIMUM_BAND_5"){ lmin5 <- as.double(words[j+2])}
        if (words[j]=="RADIANCE_MINIMUM_BAND_7"){ lmin7 <- as.double(words[j+2])}
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

    if (b1==1)
    {
      blu <- raster(paste0(directory,"/",sat_fold,"_B1.TIF"))
      if (crop == "y" || crop == "f" || crop == "u")
      {
        blu <- crop(blu,ext)
        if (crop=="f")
        {
          blu <- mask(blu,shape)
        }
      }
      rad_b1 <- ((lmax1-lmin1)/(qcal_max-qcal_min)) * (blu-qcal_min) + lmin1
      toa_blu <- pi * rad_b1 * d^2  / 1970 * sin(sun_ele*(pi/180))
      writeRaster(toa_blu,paste0(directory,"/","toa_blu"),format="GTiff",overwrite=TRUE)
    }
    if (b2==1)
    {
      green <- raster(paste0(directory,"/",sat_fold,"_B2.TIF"))
      if (crop == "y" || crop == "f" || crop == "u")
      {
        green <- crop(green,ext)
        if (crop=="f")
        {
          green <- mask(green,shape)
        }
      }
      rad_b2 <- ((lmax2-lmin2)/(qcal_max-qcal_min)) * (green-qcal_min) + lmin2
      toa_green <- pi * rad_b2 * d^2  / 1842 * sin(sun_ele*(pi/180))
      writeRaster(toa_green,paste0(directory,"/","toa_green"),format="GTiff",overwrite=TRUE)
    }

    if (b3==1)
    {
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
      writeRaster(toa_red,paste0(directory,"/","toa_red"),format="GTiff",overwrite=TRUE)
    }
    if (b4==1)
    {
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
      writeRaster(toa_nir,paste0(directory,"/","toa_nir"),format="GTiff",overwrite=TRUE)
    }
    if (b5==1)
    {
      swir1 <- raster(paste0(directory,"/",sat_fold,"_B5.TIF"))
      if (crop == "y" || crop == "f" || crop == "u")
      {
        swir <- crop(swir1,ext)
        if (crop=="f")
        {
          swir1 <- mask(swir1,shape)
        }
      }
      rad_b5 <- ((lmax5-lmin5)/(qcal_max-qcal_min)) * (swir1-qcal_min) + lmin5
      toa_swir1 <- pi * rad_b5 * d^2  / 225.7 * sin(sun_ele*(pi/180))
      writeRaster(toa_swir1,paste0(directory,"/","toa_swir1"),format="GTiff",overwrite=TRUE)
    }
    if (b7==1)
    {
      swir2 <- raster(paste0(directory,"/",sat_fold,"_B7.TIF"))
      if (crop == "y" || crop == "f" || crop == "u")
      {
        swir2 <- crop(swir2,ext)
        if (crop=="f")
        {
          swir2 <- mask(swir2,shape)
        }
      }
      rad_b7 <- ((lmax7-lmin7)/(qcal_max-qcal_min)) * (swir2-qcal_min) + lmin7
      toa_swir2 <- pi * rad_b7 * d^2  / 82.06 * sin(sun_ele*(pi/180))
      writeRaster(toa_swir2,paste0(directory,"/","toa_swir2"),format="GTiff",overwrite=TRUE)
    }
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

        if (words[j]=="RADIANCE_MAXIMUM_BAND_1"){ lmax1 <- as.double(words[j+2])}
        if (words[j]=="RADIANCE_MAXIMUM_BAND_2"){ lmax2 <- as.double(words[j+2])}
        if (words[j]=="RADIANCE_MAXIMUM_BAND_3"){ lmax3 <- as.double(words[j+2])}
        if (words[j]=="RADIANCE_MAXIMUM_BAND_4"){ lmax4 <- as.double(words[j+2])}
        if (words[j]=="RADIANCE_MAXIMUM_BAND_5"){ lmax5 <- as.double(words[j+2])}
        if (words[j]=="RADIANCE_MAXIMUM_BAND_7"){ lmax7 <- as.double(words[j+2])}

        if (words[j]=="RADIANCE_MINIMUM_BAND_1"){ lmin1 <- as.double(words[j+2])}
        if (words[j]=="RADIANCE_MINIMUM_BAND_2"){ lmin2 <- as.double(words[j+2])}
        if (words[j]=="RADIANCE_MINIMUM_BAND_3"){ lmin3 <- as.double(words[j+2])}
        if (words[j]=="RADIANCE_MINIMUM_BAND_4"){ lmin4 <- as.double(words[j+2])}
        if (words[j]=="RADIANCE_MINIMUM_BAND_5"){ lmin5 <- as.double(words[j+2])}
        if (words[j]=="RADIANCE_MINIMUM_BAND_7"){ lmin7 <- as.double(words[j+2])}
      }
    }

    if (tm_id=="\"LANDSAT_5\"")
    {
      esun2 <- 1827
      esun3 <- 1551
      esun4 <- 1036
      esun5 <- 214.9
      esun7 <- 80.65
    }

    if (tm_id!="\"LANDSAT_5\"")
    {
      esun2 <- 1826
      esun3 <- 1554
      esun4 <- 1033
      esun5 <- 214.7
      esun7 <- 80.70
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

    if (b1==1)
    {
      blu <- raster(paste0(directory,"/",sat_fold,"_B1.TIF"))
      if (crop == "y" || crop == "f" || crop == "u")
      {
        blu <- crop(blu,ext)
        if (crop=="f")
        {
          blu <- mask(blu,shape)
        }
      }
      rad_b1 <- ((lmax1-lmin1)/(qcal_max-qcal_min)) * (blu-qcal_min) + lmin1
      toa_blu <- pi * rad_b1 * d^2  / 1958 * sin(sun_ele*(pi/180))
      writeRaster(toa_blu,paste0(directory,"/","toa_blu"),format="GTiff",overwrite=TRUE)
    }
    if (b2==1)
    {
      green <- raster(paste0(directory,"/",sat_fold,"_B2.TIF"))
      if (crop == "y" || crop == "f" || crop == "u")
      {
        green <- crop(green,ext)
        if (crop=="f")
        {
          green <- mask(green,shape)
        }
      }
      rad_b2 <- ((lmax2-lmin2)/(qcal_max-qcal_min)) * (green-qcal_min) + lmin2
      toa_green <- pi * rad_b2 * d^2  / esun2 * sin(sun_ele*(pi/180))
      writeRaster(toa_green,paste0(directory,"/","toa_green"),format="GTiff",overwrite=TRUE)
    }

    if (b3==1)
    {
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
      writeRaster(toa_red,paste0(directory,"/","toa_red"),format="GTiff",overwrite=TRUE)
    }
    if (b4==1)
    {
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
      writeRaster(toa_nir,paste0(directory,"/","toa_nir"),format="GTiff",overwrite=TRUE)
    }
    if (b5==1)
    {
      swir1 <- raster(paste0(directory,"/",sat_fold,"_B5.TIF"))
      if (crop == "y" || crop == "f" || crop == "u")
      {
        swir <- crop(swir1,ext)
        if (crop=="f")
        {
          swir1 <- mask(swir1,shape)
        }
      }
      rad_b5 <- ((lmax5-lmin5)/(qcal_max-qcal_min)) * (swir1-qcal_min) + lmin5
      toa_swir1 <- pi * rad_b5 * d^2  / esun5 * sin(sun_ele*(pi/180))
      writeRaster(toa_swir1,paste0(directory,"/","toa_swir1"),format="GTiff",overwrite=TRUE)
    }
    if (b7==1)
    {
      swir2 <- raster(paste0(directory,"/",sat_fold,"_B7.TIF"))
      if (crop == "y" || crop == "f" || crop == "u")
      {
        swir2 <- crop(swir2,ext)
        if (crop=="f")
        {
          swir2 <- mask(swir2,shape)
        }
      }
      rad_b7 <- ((lmax7-lmin7)/(qcal_max-qcal_min)) * (swir2-qcal_min) + lmin7
      toa_swir2 <- pi * rad_b7 * d^2  / esun7 * sin(sun_ele*(pi/180))
      writeRaster(toa_swir2,paste0(directory,"/","toa_swir2"),format="GTiff",overwrite=TRUE)
    }
  }
  ######## Landsat TM ending ############
print("Program finished, results will be located in the satellite folder")
}
############ TOA FINISHED #######################
# ESUN values are obtained from https://landsat.usgs.gov/esun
