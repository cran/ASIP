#' Make your own custom Sentinel-2 MSI satellite image products
#'
#' This function is dedicated to Sentinel-2 MSI satellite image processing. Provide your custom equation to
#' produce the desired results (Tested only Sentinel 2 L1C products).
#' @param ext2crop,crop,directory Same as mentioned in \code{\link[ASIP]{arvi}}.
#' @param cus.formula Assign custom formula to be computed AS TEXT input (inside double quotes).
#' To assign bands, ONLY USE BAND NUMBERS (b1, b2,....,b12) to indicate different bands in the formula.
#'
#' @return Computed Sentinel 2 custom product
#' @note 1. FILENAMES OF ANY BAND FILES (*.jp2 files) SHOULDN'T CHANGED.
#'
#' 2. Bands with same resolution can only be computed.
#'
#' 2. Windows users should be careful while assigning directory. Use "/" to seperate folders not "\\".
#' @export
#' @importFrom raster raster writeRaster extent mask crop
#' @importFrom rgdal readGDAL
#' @importFrom stringr str_sub
#' @examples
#' library (raster)
#' library (rgdal)
#' # Finding the path of the sample satellite image directory.
#' # User may define paths directly like "/home/ur_folder" or "C:/ur_folder"
#' ##path <- system.file ("TM_sample", package = "ASIP")
#' # Input equation should be as text (inside double quotes)
#' eqn <- "((2 * b4)+ (b3+pi+b8))/(b3+b4+b8)"
#' ##shapefil <- paste0 (path, "/test.shp")
#' ##op <- custom.eqn (directory = path, cus.formula = eqn, crop = "y", ext2crop = shapefil)
sen2_msi <- function (directory=getwd(), cus.formula = "none", crop = "n", ext2crop = "none" )
{
  if (typeof(cus.formula)!= "character")
    stop("Please assign your formula as text in double quotes")
 # Directory check
  band_names <- list.files(directory, pattern = "*jp2")
  band_count <- length(band_names)
  if (band_count == 0)
    stop("Define your satellite image folder path properly")

  band_name <- band_names[1]
  if (stringr::str_sub(band_name,-3) == "jp2")
  {
    broke_name <- strsplit(band_name, "_B01.jp2")
    broke_name <- broke_name[[1]]
  }
  rm (band_count, band_name, band_names)
  b1 = b2 = b3 = b4 = b5= b6= b7= b8= b9 =b10= b11= b12 =0

  if (crop == "u")
  {
    b8= b4 = b3 =1
  }
  ## Getting the bands required #####
  for (i in 1:12)
  {
    i=gettext(i)
    if (grepl(paste0("b",i),cus.formula)==TRUE)
    {
      assign(paste0("b",i),1)
    }
  }
  ## Getting the bands required #####

  ##### Crop Module ######
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
  ###### Crop module ######
  i <- 0
  for (i in 1:9)
  {
    i <- gettext (i)
    if (eval(parse(text= paste0("b",i))) == 1)
    {
      s2a <- raster(paste0(directory, "/",broke_name,"_B0",i,".jp2"))
      assign(paste0("b",i), s2a)
    }
  }

  for (i in 0:2)
  {
    i <- gettext(i)
    if (eval(parse(text= paste0("b1",i))) == 1)
    {
      s2a <- raster(paste0(directory, "/",broke_name,"_B1",i,".jp2"))
      assign(paste0("b1",i), s2a)
    }
  }

  rm (i,s2a)

  # Making lists of bands and band names
  band_na = list()
  bands = list()
  j=0
  for (i in 1 :12)
  {
    if (typeof (eval(parse(text= paste0("b",i)))) == "S4")
    {
      j=j+1 #Number of bands
      band_na[[j]]=paste0("b",i)
      bands [[j]]= eval(parse(text= paste0("b",i)))
    }
  }
  # Custom extent defining
  if (crop == "u")
  {
    stak <- raster::stack(c(b8,b4,b3))
    plotRGB(stak, scale = 7000)
    cat("\nPlease define your extent from the map in plot preview for further processing\n")
    cat("\nYou can click on the top left of custom subset region followed by the bottom right\n")
    ext <- raster::drawExtent()
    rm (stak)
  }

  # Croping

  if (crop == "y" || crop == "f" || crop == "u")
  {
    for (i in 1:j)
    {
      bands[[i]] = raster::crop (bands[[i]],ext)
      if (crop == "f")
      {
       bands[[i]] <- raster::mask(bands[[i]],shape)
      }
      assign(band_na[[i]], bands[[i]])
    }
  }
  if (crop == "u") {rm (ext)}
  if (crop == "f") {rm (shape)}
## Band calculator
  sen2_cus_eqn <- eval(parse(text = cus.formula))

  rm (b1,b2,b3, b4, b5,b6,b7,b8, b9,b10,b11,b12, cus.formula, band_na, bands, broke_name, crop, ext2crop, i ,j)
  #raster::writeRaster(op, paste0 (directory,"/output_ASIP"),format="GTiff", overwrite=TRUE)
  return (sen2_cus_eqn)
  rm (directory)
  cat("\nProgram completed, output is produced as a variable named 'sen2_cus_eqn'")
}
