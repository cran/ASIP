#' TIR bands to at satellite brightness temperature conversion
#'
#' Identifies Thermal Infra-Red (TIR) bands and converts them to  at satellite brightness temperature images.
#' @param ext2crop,crop,directory Same as mentioned in \code{\link[ASIP]{arvi}}.
#' @param unit By default the temperature image will be produced in Degree Kelvin. To produce the thermal image in Degree celcius, assign vale "c".
#' To produce the thermal image in Degree celcius, assign vale "c".
#' @return At Satellite Brightness Temperature images in .tif format in input directory.
#' @note 1. FILENAMES OF ANY BAND FILES (*.TIF files) SHOULDN'T CHANGED.
#'
#' 2. Windows users should be careful while assigning directory. Use "/" to seperate folders.
#'
#' 3. Emissivity value used is 1.
#' @export
#' @importFrom raster raster writeRaster extent mask crop
#' @importFrom utils tail
#' @examples
#' library (raster)
#' library (rgdal)
#' # Finding the path of the sample satellite image directory.
#' # User may define paths directly like "/home/ur_folder" or "C:/ur_folder"
#' path <- system.file ("TM_sample", package = "ASIP")
#' shapefil <- paste0 (path, "/test.shp")
#' thermal (directory = path, crop = "y", ext2crop = shapefil, unit = "c")
# DN to thermal
thermal <- function(directory = getwd(), crop = "n", ext2crop = "none",unit = "Deg Kel")
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
        if (words[j]=="RADIANCE_MULT_BAND_10"){ tir1_rad_mult <- as.double(words[j+2])}
        if (words[j]=="RADIANCE_ADD_BAND_10"){ tir1_rad_add <- as.double(words[j+2])}
        if (words[j]=="RADIANCE_MULT_BAND_11"){ tir2_rad_mult <- as.double(words[j+2])}
        if (words[j]=="RADIANCE_ADD_BAND_11"){ tir2_rad_add <- as.double(words[j+2])}
        if (words[j]=="K1_CONSTANT_BAND_10"){ tir1_k1 <- as.double(words[j+2])}
        if (words[j]=="K2_CONSTANT_BAND_10"){ tir1_k2 <- as.double(words[j+2])}
        if (words[j]=="K1_CONSTANT_BAND_11"){ tir2_k1 <- as.double(words[j+2])}
        if (words[j]=="K2_CONSTANT_BAND_11"){ tir2_k2 <- as.double(words[j+2])}

        if (words[j]=="DATE_ACQUIRED"){ data_aq <- as.character(words[j+2])}
        if (words[j]=="SUN_ELEVATION"){ sun_ele <- as.double(words[j+2])}
      }
    }
    # Defining bands & toa calculation
   tir2 <- as.integer(raster(paste0(directory,"/",sat_fold,"_B11.TIF")))
   if (crop == "y" || crop == "f" || crop == "u")
      {
        tir2 <- crop(tir2, ext)
        if (crop=="f")
        {
          tir2 <- mask(tir2,shape)
        }
      }
   rad_tir2 <- (tir2 * tir2_rad_mult) + tir2_rad_add
   temp_tir2 <- tir2_k2/ log((tir2_k1/rad_tir2)+1)
   if (unit=="c")
   {temp_tir2 <- temp_tir2-273.15}
   writeRaster(temp_tir2,paste0(directory,"/","temp_tir2"),format="GTiff",overwrite=TRUE)

   tir1 <- as.integer(raster(paste0(directory,"/",sat_fold,"_B10.TIF")))
   if (crop == "y" || crop == "f" || crop == "u")
      {
        tir1 <- crop(tir1, ext)
        if (crop=="f")
        {
          tir1 <- mask(tir1,shape)
        }
      }
   rad_tir1 <- (tir1 * tir1_rad_mult) + tir1_rad_add
   temp_tir1 <- tir1_k2/ log((tir1_k1/rad_tir1)+1)
   if (unit=="c")
   {temp_tir1 <- temp_tir1-273.15}
   writeRaster(temp_tir1,paste0(directory,"/","temp_tir1"),format="GTiff",overwrite=TRUE)
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
    tir1_k1 <- tir2_k1 <- tir1_k2 <- tir2_k2 <- 999
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

        if (words[j]=="RADIANCE_MAXIMUM_BAND_6_VCID_1"){ lmax61 <- as.double(words[j+2])}
        if (words[j]=="RADIANCE_MAXIMUM_BAND_6_VCID_2"){ lmax62 <- as.double(words[j+2])}

        if (words[j]=="RADIANCE_MINIMUM_BAND_6_VCID_1"){ lmin61 <- as.double(words[j+2])}
        if (words[j]=="RADIANCE_MINIMUM_BAND_6_VCID_2"){ lmin62 <- as.double(words[j+2])}

        if (words[j]=="K1_CONSTANT_BAND_6_VCID_1"){ tir1_k1 <- as.double(words[j+2])}
        if (words[j]=="K2_CONSTANT_BAND_6_VCID_1"){ tir1_k2 <- as.double(words[j+2])}
        if (words[j]=="K1_CONSTANT_BAND_6_VCID_2"){ tir2_k1 <- as.double(words[j+2])}
        if (words[j]=="K2_CONSTANT_BAND_6_VCID_2"){ tir2_k2 <- as.double(words[j+2])}
      }
    }

    if (tir1_k1==999){tir1_k1 <- 666.09}
    if (tir1_k2==999){tir1_k2 <- 1282.71}
    if (tir2_k1==999){tir2_k1 <- 666.09}
    if (tir2_k2==999){tir2_k2 <- 1282.71}

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

   tir1 <- raster(paste0(directory,"/",sat_fold,"_B6_VCID_1.TIF"))
   if (crop == "y" || crop == "f" || crop == "u")
   {
     tir1 <- crop(tir1,ext)
        if (crop=="f")
        {
          tir1 <- mask(tir1,shape)
        }
   }
    rad_b61 <- ((lmax61-lmin61)/(qcal_max-qcal_min)) * (tir1-qcal_min) + lmin61
    tir1 <- tir1_k2/(log((tir1_k1/rad_b61)+1))
    if (unit=="c")
    {tir1 <- tir1-273.15}
    writeRaster(tir1,paste0(directory,"/","tir1"),format="GTiff",overwrite=TRUE)

    tir2 <- raster(paste0(directory,"/",sat_fold,"_B6_VCID_2.TIF"))
    if (crop == "y" || crop == "f" || crop == "u")
    {
      tir2 <- crop(tir2,ext)
      if (crop=="f")
      {
        tir2 <- mask(tir2,shape)
      }
    }
    rad_b62 <- ((lmax62-lmin62)/(qcal_max-qcal_min)) * (tir2-qcal_min) + lmin62
    tir2 <- tir2_k2/(log((tir2_k1/rad_b62)+1))
    if (unit=="c")
    {tir2 <- tir2-273.15}
    writeRaster(tir2,paste0(directory,"/","tir2"),format="GTiff",overwrite=TRUE)
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
    k1 <- k2 <- 999
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

        if (words[j]=="RADIANCE_MAXIMUM_BAND_6"){ lmax6 <- as.double(words[j+2])}
        if (words[j]=="RADIANCE_MINIMUM_BAND_6"){ lmin6 <- as.double(words[j+2])}

        if (words[j]=="K1_CONSTANT_BAND_6"){ k1 <- as.double(words[j+2])}
        if (words[j]=="K2_CONSTANT_BAND_6"){ k2 <- as.double(words[j+2])}
      }
    }

    if (k1==999){k1 <- 607.76}
    if (k2==999){k2 <- 1260.56}

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

    tir <- raster(paste0(directory,"/",sat_fold,"_B6.TIF"))
    if (crop == "y" || crop == "f" || crop == "u")
    {
    tir <- crop(tir,ext)
    if (crop=="f")
      {
        tir <- mask(tir,shape)
      }
    }
    rad_b6 <- ((lmax6-lmin6)/(qcal_max-qcal_min)) * (tir-qcal_min) + lmin6
    tir <- k2/(log((k1/rad_b6)+1))
    if (unit=="c")
    {tir <- tir-273.15}
    writeRaster(tir,paste0(directory,"/","tir"),format="GTiff",overwrite=TRUE)
  }
  print("Program finished, results will be located in the satellite folder")
  }
  ######## Landsat TM ending ############
