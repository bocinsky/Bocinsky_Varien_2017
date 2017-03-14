distill <- function(file, gray = T){
  system(paste0("gs -sDEVICE=pdfwrite -dNOPAUSE -dQUIET -dBATCH -dCompatibilityLevel=1.7 -sColorConversionStrategy=/CMYK -dPDFSETTINGS=/prepress -dEmbedAllFonts=true -dSubsetFonts=false -dAutoRotatePages=/None -sOutputFile=./temp.pdf ",file))
  system(paste0("rm ",file))
  system(paste0("mv ./temp.pdf ",file))
  
  if(gray){
    system(paste0("gs -sDEVICE=pdfwrite -dNOPAUSE -dQUIET -dBATCH -dCompatibilityLevel=1.7 -sProcessColorModel=DeviceGray -sColorConversionStrategy=Gray -dPDFSETTINGS=/prepress -dEmbedAllFonts=true -dSubsetFonts=false -dAutoRotatePages=/None -sOutputFile=",gsub(".pdf","_gray.pdf",file)," ",file))
  }
  
  system(paste0("convert -density 600 ",file," ",gsub(".pdf",".png",file)))
  if(gray){
    system(paste0("convert -colorspace gray -density 600 ",gsub(".pdf","_gray.pdf",file)," ",gsub(".pdf","_gray.png",file)))
  }
  
}