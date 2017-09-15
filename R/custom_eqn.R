#' Make your own custom satellite image product
#'
#' If any product or index is not available in this package, you don't need to do it manually.
#' This function intakes a custom formula & produced new product according to the formula.
#' This function converts DN bands to corresponding TOA reflectance prior to the computation of user defined formula.
#' @param ext2crop,crop,directory Same as mentioned in \code{\link[ASIP]{arvi}}.
#' @param cus.formula Assign custom formula to be computed AS TEXT input (inside double quotes).
#' To assign bands, ONLY USE BELOW DEFINED WORDS to indicate different bands in the formula.
#'
#' nir for NIR (Near Infra-red) Top Of Atmosphere (TOA) reflectance band.
#'
#' red for Red TOA reflectance band.
#'
#' green for Green TOA reflectance band.
#'
#' blue for Blue TOA reflectance band.
#'
#' swir1 for SWIR-1 (Short Wave Infra-red -1)
#'
#' swir2 for SWIR-2 (Short Wave Infra-red -2)
#'
#' aero for Aerosol/coastal band (Only on Landsat OLI images)
#' @return File named ur raster_'date of satellite image acqisition'.tif in the input folder
#' @note 1. FILENAMES OF ANY BAND FILES (*.TIF files) SHOULDN'T CHANGED.
#'
#' 2. Windows users should be careful while assigning directory. Use "/" to seperate folders not "\\".
#'
#' 3. Earth-sun distance is calculated according to Epema (1992) if the value is not mentioned in the meta data (*MTL.txt) file.
#'
#' 4. Currently recommended ESUN values provided by \href{https://landsat.usgs.gov/esun}{USGS} is used.
#' @export
#' @references \href{http://www.tandfonline.com/doi/ref/10.1080/01431169208904159}{Epema G F (1992) Atmospheric condition and its influence on reflectance
#' of bare soil surfaces in southern Tunisia. International Journal of Remote Sensing, 13(5), pp:853-868. doi:10.1080/01431169208904159.}
#' @importFrom raster raster writeRaster extent mask crop
#' @importFrom utils tail
#' @examples
#' library (raster)
#' library (rgdal)
#' # Finding the path of the sample satellite image directory.
#' # User may define paths directly like "/home/ur_folder" or "C:/ur_folder"
#' path <- system.file ("TM_sample", package = "ASIP")
#' # Input equation should be as text (inside double quotes)
#' eqn <- "(2* (nir^2)+ (2.20 + green))/ (blue / (2 * pi))"
#' shapefil <- paste0 (path, "/test.shp")
#' custom.eqn (directory = path, cus.formula = eqn, crop = "y", ext2crop = shapefil)
custom.eqn <- function (directory=getwd(), cus.formula = "none", crop = "n", ext2crop = "none" )
{
  # If the directory is not set
  bands <- length (list.files(directory,pattern = "*TIF"))
  if (bands == 0)
    stop("Define your satellite image folder path properly")

  if (typeof(cus.formula)!= "character")
    stop("Please assign your formula as text in double quotes")

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
      stak <- raster::stack (c (b5, b4, b3))
      plotRGB (stak, scale = 65536)
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

    nir <- as.integer(raster::raster(paste0(directory,"/",sat_fold,"_B5.TIF")))
    if (crop == "y" || crop == "f" || crop == "u")
    {
      nir <- raster::crop(nir, ext)
      if (crop=="f")
      {
        nir <- raster::mask(nir,shape)
      }
    }
    nir <- ((nir * nir_refl_mult) + nir_refl_add)/sin(sun_ele*(pi/180))

    red <- as.integer(raster(paste0(directory,"/",sat_fold,"_B4.TIF")))
    if (crop == "y" || crop == "f" || crop == "u")
    {
      red <- raster::crop(red, ext)
      if (crop=="f")
      {
        red <- raster::mask(red,shape)
      }
    }
    red <- ((red * red_refl_mult) + red_refl_add)/sin(sun_ele*(pi/180))

    blu <- as.integer(raster(paste0(directory,"/",sat_fold,"_B2.TIF")))
    if (crop == "y" || crop == "f" || crop == "u")
    {
      blu <- raster::crop(blu, ext)
      if (crop=="f")
      {
        blu <- raster::mask(blu,shape)
      }
    }
    blue <- ((blu * blu_refl_mult) + blu_refl_add)/sin(sun_ele*(pi/180))

    green <- as.integer(raster(paste0(directory,"/",sat_fold,"_B3.TIF")))
    if (crop == "y" || crop == "f" || crop == "u")
    {
      green <- raster::crop(green, ext)
      if (crop=="f")
      {
        green <- raster::mask(green,shape)
      }
    }
    green <- ((green * green_refl_mult) + green_refl_add)/sin(sun_ele*(pi/180))

    swir1 <- as.integer(raster(paste0(directory,"/",sat_fold,"_B6.TIF")))
    if (crop == "y" || crop == "f" || crop == "u")
    {
      swir1 <- raster::crop(swir1, ext)
      if (crop=="f")
      {
        swir1 <- raster::mask(swir1,shape)
      }
    }
    swir1 <- ((swir1 * swir1_refl_mult) + swir1_refl_add)/sin(sun_ele*(pi/180))

    swir2 <- as.integer(raster(paste0(directory,"/",sat_fold,"_B7.TIF")))
    if (crop == "y" || crop == "f" || crop == "u")
    {
      swir2 <- raster::crop(swir2, ext)
      if (crop=="f")
      {
        swir2 <- raster::mask(swir2,shape)
      }
    }
    swir2 <- ((swir2 * swir2_refl_mult) + swir2_refl_add)/sin(sun_ele*(pi/180))
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

    red <- raster(paste0(directory,"/",sat_fold,"_B3.TIF"))
    if (crop == "y" || crop == "f" || crop == "u")
    {
      red <- raster::crop(red,ext)
      if (crop=="f")
      {
        red <- raster::mask(red,shape)
      }
    }
    rad_b3 <- ((lmax3-lmin3)/(qcal_max-qcal_min)) * (red-qcal_min) + lmin3
    red <- pi * rad_b3 * d^2  / 1547 * sin(sun_ele*(pi/180))

    nir <- raster(paste0(directory,"/",sat_fold,"_B4.TIF"))
    if (crop == "y" || crop == "f" || crop == "u")
    {
      nir <- raster::crop(nir,ext)
      if (crop=="f")
      {
        nir <- raster::mask(nir,shape)
      }
    }
    rad_b4 <- ((lmax4-lmin4)/(qcal_max-qcal_min)) * (nir-qcal_min) + lmin4
    nir <- pi * rad_b4 * d^2  / 1044 * sin(sun_ele*(pi/180))

    blu <- raster(paste0(directory,"/",sat_fold,"_B1.TIF"))
    if (crop == "y" || crop == "f" || crop == "u")
    {
      blu <- raster::crop(blu,ext)
      if (crop=="f")
      {
        blu <- raster::mask(blu,shape)
      }
    }
    rad_b1 <- ((lmax1-lmin1)/(qcal_max-qcal_min)) * (blu-qcal_min) + lmin1
    blue <- pi * rad_b1 * d^2  / 1970 * sin(sun_ele*(pi/180))


    green <- raster(paste0(directory,"/",sat_fold,"_B2.TIF"))
    if (crop == "y" || crop == "f" || crop == "u")
    {
      green <- raster::crop(green,ext)
      if (crop=="f")
      {
        green <- raster::mask(green,shape)
      }
    }
    rad_b2 <- ((lmax2-lmin2)/(qcal_max-qcal_min)) * (green-qcal_min) + lmin2
    green <- pi * rad_b2 * d^2  / 1842 * sin(sun_ele*(pi/180))

    swir1 <- raster(paste0(directory,"/",sat_fold,"_B5.TIF"))
    if (crop == "y" || crop == "f" || crop == "u")
    {
      swir1 <- raster::crop(swir1,ext)
      if (crop=="f")
      {
        swir1 <- raster::mask(swir1,shape)
      }
    }
    rad_b5 <- ((lmax5-lmin5)/(qcal_max-qcal_min)) * (swir1-qcal_min) + lmin5
    swir1 <- pi * rad_b5 * d^2  / 225.7 * sin(sun_ele*(pi/180))

    swir2 <- raster(paste0(directory,"/",sat_fold,"_B7.TIF"))
    if (crop == "y" || crop == "f" || crop == "u")
    {
      swir2 <- raster::crop(swir2,ext)
      if (crop=="f")
      {
        swir2 <- raster::mask(swir2,shape)
      }
    }
    rad_b7 <- ((lmax7-lmin7)/(qcal_max-qcal_min)) * (swir2-qcal_min) + lmin7
    swir2 <- pi * rad_b7 * d^2  / 82.06 * sin(sun_ele*(pi/180))
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

    red <- raster(paste0(directory,"/",sat_fold,"_B3.TIF"))
    if (crop == "y" || crop == "f" || crop == "u")
    {
      red <- raster::crop(red,ext)
      if (crop=="f")
      {
        red <- raster::mask(red,shape)
      }
    }
    rad_b3 <- ((lmax3-lmin3)/(qcal_max-qcal_min)) * (red-qcal_min) + lmin3
    red <- pi * rad_b3 * d^2  / esun3 * sin(sun_ele*(pi/180))

    nir <- raster(paste0(directory,"/",sat_fold,"_B4.TIF"))
    if (crop == "y" || crop == "f" || crop == "u")
    {
      nir <- raster::crop(nir,ext)
      if (crop=="f")
      {
        nir <- raster::mask(nir,shape)
      }
    }
    rad_b4 <- ((lmax4-lmin4)/(qcal_max-qcal_min)) * (nir-qcal_min) + lmin4
    nir <- pi * rad_b4 * d^2  / esun4 * sin(sun_ele*(pi/180))

    blu <- raster(paste0(directory,"/",sat_fold,"_B1.TIF"))
    if (crop == "y" || crop == "f" || crop == "u")
    {
      blu <- raster::crop(blu,ext)
      if (crop=="f")
      {
        blu <- raster::mask(blu,shape)
      }
    }
    rad_b1 <- ((lmax1-lmin1)/(qcal_max-qcal_min)) * (blu-qcal_min) + lmin1
    blue <- pi * rad_b1 * d^2  / 1958 * sin(sun_ele*(pi/180))

    green <- raster(paste0(directory,"/",sat_fold,"_B2.TIF"))
    if (crop == "y" || crop == "f" || crop == "u")
    {
      green <- raster::crop(green,ext)
      if (crop=="f")
      {
        green <- raster::mask(green,shape)
      }
    }
    rad_b2 <- ((lmax2-lmin2)/(qcal_max-qcal_min)) * (green-qcal_min) + lmin2
    green <- pi * rad_b2 * d^2  / esun2 * sin(sun_ele*(pi/180))

    swir1 <- raster(paste0(directory,"/",sat_fold,"_B5.TIF"))
    if (crop == "y" || crop == "f" || crop == "u")
    {
      swir1 <- raster::crop(swir1,ext)
      if (crop=="f")
      {
        swir1 <- raster::mask(swir1,shape)
      }
    }
    rad_b5 <- ((lmax5-lmin5)/(qcal_max-qcal_min)) * (swir1-qcal_min) + lmin5
    swir1 <- pi * rad_b5 * d^2  / esun5 * sin(sun_ele*(pi/180))

    swir2 <- raster(paste0(directory,"/",sat_fold,"_B7.TIF"))
    if (crop == "y" || crop == "f" || crop == "u")
    {
      swir2 <- raster::crop(swir2,ext)
      if (crop=="f")
      {
        swir2 <- raster::mask(swir2,shape)
      }
    }
    rad_b7 <- ((lmax7-lmin7)/(qcal_max-qcal_min)) * (swir2-qcal_min) + lmin7
    swir2 <- pi * rad_b7 * d^2  / esun7 * sin(sun_ele*(pi/180))
  }
  ######## Landsat TM ending ############
  op <- eval(parse(text = cus.formula))
  raster::writeRaster(op,paste0(directory,"/","ur raster_",data_aq),format="GTiff", overwrite=TRUE)
  print("Program completed, output is named as 'ur raster_[date of data acquisition].tif' in satellite image folder")
}
