# Read and do some basic cleaning on the Colitis data set (Burczynski et al., 2006)
# (Note: See 
# https://mdozmorov.github.io/BIOS567/assets/presentation_Bioconductor/GEO.pdf
# for a quick primer on this format)

#===== libraries =====#
#source("http://www.bioconductor.org/biocLite.R")
#biocLite("GEOquery")
library(Biobase)
library(GEOquery)

#===== load file ====#
GDS1615 <- getGEO(filename = "data/GDS1615_full.soft.gz")#, destdir = "data/")

#==== extract and restructure ====#
dat <- Table(GDS1615)
dat2 <- dat[,3:129]

x <- matrix(unlist(dat2), ncol = ncol(dat2))
y <- factor(GDS1615@dataTable@columns$disease.state, 
            levels = c("normal", "ulcerative colitis", "Crohn's disease"))
y <- as.numeric(y)

#==== write to disk =====#
colitis <- list(x = x, y = y)
save(colitis, file = "data/colitis.rda")






