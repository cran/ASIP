#' Cropping of satellite image bands
#'
#' Crop desired satellite image bands either using a shapefile or draw custom extent from a plot image while running the function.
#' @param ext2crop,crop,directory Same as mentioned in \code{\link[ASIP]{arvi}}.
#' @param b1 By default Band1 will be cropped. To cancel cropping of this band assign value 0.
#' @param b2 By default Band2 will be cropped. To cancel cropping of this band assign value 0.
#' @param b3 By default Band3 will be cropped. To cancel cropping of this band assign value 0.
#' @param b4 By default Band4 will be cropped. To cancel cropping of this band assign value 0.
#' @param b5 By default Band5 will be cropped. To cancel cropping of this band assign value 0.
#' @param b6 By default Band6 will be cropped. To cancel cropping of this band assign value 0.
#' @param b7 By default Band7 will be cropped. To cancel cropping of this band assign value 0.
#' @return Each bands selected will cropped and produce corresponding <bandname>_crop.tif format in the input directory.
#' @note 1. FILENAMES OF ANY BAND FILES (*.TIF files) SHOULDN'T CHANGED.
#'
#' 2. Windows users should be careful while assigning directory. Use "/" to seperate folders not "\\".
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
#' # Assign 0 values to band names which are not required
#' crop.bands (path, crop = "f", ext2crop = shapefil, b3=0, b4=0, b5=0, b6 = 0, b7 = 0)
crop.bands <- function (directory= getwd(), crop = "n", ext2crop = "none",b1=1,b2=1,b3=1,b4=1,b5=1,b6=1,b7=1)
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
  ######### Landsat 8 starting###############
  if (satellite=="LC")
  {
    nir <- as.integer(raster(paste0(directory,"/",sat_fold,"_B5.TIF")))
    red <- as.integer(raster(paste0(directory,"/",sat_fold,"_B4.TIF")))
    green <- as.integer(raster(paste0(directory,"/",sat_fold,"_B3.TIF")))
    if (crop == "u")
    {
      stak <- raster::stack(c(nir,red,green))
      plotRGB(stak, scale = 65536)
      print("Please define your extent from the map in plot preview for further processing")
      print("You can click on the top left of custom subset region followed by the bottom right")
      ext <- drawExtent()
    }
    # Defining bands

    if (b7==1)
    {
      swir2 <- as.integer(raster(paste0(directory,"/",sat_fold,"_B7.TIF")))
      swir2_crop <- crop(swir2, ext)
      if (crop=="f")
      {
        swir2_crop <- mask(swir2_crop,shape)
      }
      writeRaster(swir2_crop,paste0(directory,"/","swir2_crop"),format="GTiff",overwrite=TRUE)
    }

    if (b6==1)
    {
      swir1 <- as.integer(raster(paste0(directory,"/",sat_fold,"_B6.TIF")))
      swir1_crop <- crop(swir1, ext)
      if (crop=="f")
      {
        swir1_crop <- mask(swir1_crop,shape)
      }
      writeRaster(swir1_crop,paste0(directory,"/","swir1_crop"),format="GTiff",overwrite=TRUE)
    }

    if (b5==1)
    {
      nir_crop <- crop(nir, ext)
      if (crop=="f")
      {
        nir_crop <- mask(nir_crop,shape)
      }
      writeRaster(nir_crop,paste0(directory,"/","nir_crop"),format="GTiff",overwrite=TRUE)
    }

    if (b4==1)
    {
      red_crop <- crop(red, ext)
      if (crop=="f")
      {
        red_crop <- mask(red_crop,shape)
      }
      writeRaster(red_crop,paste0(directory,"/","red_crop"),format="GTiff",overwrite=TRUE)
    }

    if (b3==1)
    {
      green_crop <- crop(green, ext)
      if (crop=="f")
      {
        green_crop <- mask(green_crop,shape)
      }
      writeRaster(green_crop,paste0(directory,"/","green_crop"),format="GTiff",overwrite=TRUE)
    }

    if (b2==1)
    {
      blu <- as.integer(raster(paste0(directory,"/",sat_fold,"_B2.TIF")))
      blu_crop <- crop(blu, ext)
      if (crop=="f")
      {
        blu_crop <- mask(blu_crop,shape)
      }
      writeRaster(blu_crop,paste0(directory,"/","blu_crop"),format="GTiff",overwrite=TRUE)
    }

    if (b1==1)
    {
      aero <- as.integer(raster(paste0(directory,"/",sat_fold,"_B1.TIF")))
      aero_crop <- crop(aero, ext)
      if (crop=="f")
      {
        aero_crop <- mask(aero_crop,shape)
      }
      writeRaster(aero_crop,paste0(directory,"/","aero_crop"),format="GTiff",overwrite=TRUE)
    }
  }
  ########### Landsat-8 ending ##############
  ########### Landsat-ETM+ and TM starting ##############
  if (satellite=="LE" || satellite=="LT")
  {
    nir <- raster(paste0(directory,"/",sat_fold,"_B4.TIF"))
    red <- raster(paste0(directory,"/",sat_fold,"_B3.TIF"))
    green <- raster(paste0(directory,"/",sat_fold,"_B2.TIF"))
    if (crop == "u")
    {
      stak <- stack(c (nir, red, green))
      plotRGB(stak)
      print("Please define your extent from the map in plot preview for further processing")
      print("You can click on the top left of custom subset region followed by the bottom right")
      ext <- drawExtent()
    }

    if (b1==1)
    {
      blu <- raster(paste0(directory,"/",sat_fold,"_B1.TIF"))
      blu_crop <- crop(blu,ext)
      if (crop=="f")
      {
        blu_crop <- mask(blu_crop,shape)
      }
      writeRaster(blu_crop,paste0(directory,"/","blu_crop"), format= "GTiff", overwrite=TRUE)
    }

    if (b2==1)
    {
      green_crop <- crop(green,ext)
      if (crop == "f")
      {
        green_crop <- mask(green_crop,shape)
      }
      writeRaster(green_crop,paste0(directory,"/","green_crop"),format="GTiff",overwrite=TRUE)
    }

    if (b3==1)
    {
      red_crop <- crop(red,ext)
      if (crop == "f")
      {
        red_crop <- mask(red_crop,shape)
      }
      writeRaster (red_crop,paste0(directory,"/","red_crop"),format="GTiff",overwrite=TRUE)
    }

    if (b4==1)
    {
      nir_crop <- crop(nir,ext)
      if (crop=="f")
      {
        nir_crop <- mask(nir_crop,shape)
      }
      writeRaster(nir_crop,paste0(directory,"/","nir_crop"),format="GTiff",overwrite=TRUE)
    }

    if (b5==1)
    {
      swir1 <- raster(paste0(directory,"/",sat_fold,"_B5.TIF"))
      swir1_crop <- crop(swir1,ext)
      if (crop=="f")
      {
        swir1_crop <- mask(swir1_crop,shape)
      }
      writeRaster(swir1_crop,paste0(directory,"/","swir1_crop"),format="GTiff",overwrite=TRUE)
    }

    if (b7==1)
    {
      swir2 <- raster(paste0(directory,"/",sat_fold,"_B7.TIF"))
      swir2_crop <- crop(swir2,ext)
      if (crop=="f")
      {
        swir2_crop <- mask(swir2_crop,shape)
      }
      writeRaster(swir2_crop,paste0(directory,"/","swir2_crop"),format="GTiff",overwrite=TRUE)
    }
  }
}
  ############## Landsat ETM & TM ending ##################
