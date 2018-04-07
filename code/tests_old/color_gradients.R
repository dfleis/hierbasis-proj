parula_cols <-
  c(
    rgb(53,  42,  135, maxColorValue = 255),
    rgb(15,  92,  221, maxColorValue = 255),
    rgb(18,  125, 216, maxColorValue = 255),
    rgb(7,   156, 207, maxColorValue = 255),
    rgb(21,  177, 180, maxColorValue = 255),
    rgb(89,  189, 140, maxColorValue = 255),
    rgb(165, 190, 107, maxColorValue = 255),
    rgb(225, 185, 82,  maxColorValue = 255),
    rgb(252, 206, 46,  maxColorValue = 255),
    rgb(249, 251, 14,  maxColorValue = 255)
  )
parula <- colorRampPalette(parula_cols)

jet_cols <- c("#00008F", "#00009F", "#0000AF", "#0000BF", "#0000CF", 
              "#0000DF", "#0000EF", "#0000FF", "#0010FF", "#0020FF", 
              "#0030FF", "#0040FF", "#0050FF", "#0060FF", "#0070FF", 
              "#0080FF", "#008FFF", "#009FFF", "#00AFFF", "#00BFFF", 
              "#00CFFF", "#00DFFF", "#00EFFF", "#00FFFF", "#10FFEF", 
              "#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF", "#60FF9F", 
              "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60", "#AFFF50", 
              "#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10", "#FFFF00", 
              "#FFEF00", "#FFDF00", "#FFCF00", "#FFBF00", "#FFAF00", 
              "#FF9F00", "#FF8F00", "#FF8000", "#FF7000", "#FF6000", 
              "#FF5000", "#FF4000", "#FF3000", "#FF2000", "#FF1000", 
              "#FF0000", "#EF0000", "#DF0000", "#CF0000", "#BF0000", 
              "#AF0000", "#9F0000", "#8F0000", "#800000")
jet <- colorRampPalette(jet_cols)
