#' Normalized Difference Water Index
#'
#' Normalized Difference Water Index (NDWI) is used to monitor changes related to water content in water bodies proposed by McFeeters (1996).
#' @param ext2crop,crop,directory Same as mentioned in \code{\link[ASIP]{arvi}}.
#' @return File named ndwi_'date of satellite image acqisition'.tif in the input folder
#' @note 1. NDWI = (r_green - r_nir) / (r_nir + r_green)
#'
#' where, "r_" denotes TOA reflectance band.
#'
#' 2. There is another NDWI to monitor changes in water content of leaves proposed by Gao (1996).
#' User should understand the requirements and run accordingly.
#'
#' Other important notes are mentioned in \code{\link[ASIP]{custom.eqn}}.
#' @export
#' @importFrom raster raster writeRaster extent mask crop
#' @importFrom utils tail
#' @references \href{http://www.tandfonline.com/doi/abs/10.1080/01431169608948714}{McFeeters, S.K. (1996) The Use of the Normalized Difference Water Index (NDWI) in the Delineation of Open Water Features. International Journal of Remote Sensing, 17, 1425-1432.
#' doi:10.1080/01431169608948714.}
#'
#' \href{http://www.sciencedirect.com/science/article/pii/S0034425796000673}{Gao
#' Bo-cai (1996) NDWI-A normalized difference water index for remote sensing of vegetation liquid water from space. Remote Sensing of Environment, 58 (3), 257-266. doi:10.1016/S0034-4257(96)00067-3.}
#' @examples
#' library (raster)
#' library (rgdal)
#' # Finding the path of the sample satellite image directory.
#' # User may define paths directly like "/home/ur_folder" or "C:/ur_folder"
#' path <- system.file ("TM_sample", package = "ASIP")
#' shapefil <- paste0 (path, "/test.shp")
#' ndwi (directory = path, crop = "y", ext2crop = shapefil)
# ndwi function
ndwi <- function(directory = getwd(), crop = "n", ext2crop = "none")
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
      stak <- stack(c(b5,b4,b3))
      plotRGB(stak)
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
        if (words[j]=="REFLECTANCE_ADD_BAND_3"){ green_refl_add <- as.double(words[j+2])}
        if (words[j]=="REFLECTANCE_MULT_BAND_3"){ green_refl_mult <- as.double(words[j+2])}

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

        if (words[j]=="RADIANCE_MAXIMUM_BAND_2"){ lmax2 <- as.double(words[j+2])}
        if (words[j]=="RADIANCE_MAXIMUM_BAND_4"){ lmax4 <- as.double(words[j+2])}

        if (words[j]=="RADIANCE_MINIMUM_BAND_2"){ lmin2 <- as.double(words[j+2])}
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

        if (words[j]=="RADIANCE_MAXIMUM_BAND_2"){ lmax2 <- as.double(words[j+2])}
        if (words[j]=="RADIANCE_MAXIMUM_BAND_4"){ lmax4 <- as.double(words[j+2])}

        if (words[j]=="RADIANCE_MINIMUM_BAND_2"){ lmin2 <- as.double(words[j+2])}
        if (words[j]=="RADIANCE_MINIMUM_BAND_4"){ lmin4 <- as.double(words[j+2])}
      }
    }

    if (tm_id=="\"LANDSAT_5\"")
    {
      esun2 <- 1827
      esun4 <- 1036
    }

    if (tm_id!="\"LANDSAT_5\"")
    {
      esun2 <- 1826
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
  ndwi <- (toa_green-toa_nir)/(toa_nir+toa_green)
  writeRaster(ndwi,paste0(directory,"/","ndwi_",data_aq),format="GTiff", overwrite=TRUE)
  print("Program completed, output is named as 'ndwi_[date of data acquisition].tif' in satellite image folder")
}
