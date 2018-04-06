setwd("~/Shared/Data-Science/Data-Source-Model-Repository/MedGen/scripts/")

library(XML)
library(parallel)

##
mc.cores <- 55
sdir <- "../sources"
ddir <- "../data"

###############################################################################@
## Source information ----
###############################################################################@

sfi <- read.table(
   file.path(sdir, "ARCHIVES/ARCHIVES.txt"),
   sep="\t",
   header=T,
   stringsAsFactors=FALSE
)
MedGen_sourceFiles <- sfi[which(sfi$inUse), c("url", "current")]

###############################################################################@
## Data from medgen ----
###############################################################################@
## decompress gz
for(f in sfi$file){
  gzf <- file.path(sdir,f)
  system(paste0("gzip -d ", file.path(sdir,f)))
  }

###########################
## Basic information

## crossId
MedGen_conso <- read.table(file.path(sdir,"MGCONSO.RFF"), sep = "|", header = TRUE, comment.char = "", quote = "", 
                      fill = TRUE, colClasses = c("character"))
MedGen_conso$SAB[MedGen_conso$SAB == "MSH"] <- "MeSH"
MedGen_conso$SAB[MedGen_conso$SAB == "NCI"] <- "NCIt"
MedGen_conso$SAB[MedGen_conso$SAB == "SNOMEDCT_US"] <- "SNOMEDCT"
MedGen_conso$SAB[MedGen_conso$SAB == "ORDO"] <- "ORPHA"
MedGen_conso$DB <- "MedGen"

MedGen_crossId <- MedGen_conso[,c("DB","X.CUI","SAB","SDUI")]
names(MedGen_crossId) <- c("DB1","id1","DB2","id2")

## entryId
MedGen_entryId <- MedGen_conso[!duplicated(MedGen_conso$X.CUI), c("DB","X.CUI")]
names(MedGen_entryId) <- c("DB","id")

MedGen_def <- read.table(file.path(sdir,"MGDEF.RFF"), sep = "|", header = TRUE, comment.char = "", quote = "", 
                           fill = TRUE, colClasses = c("character"))
MedGen_entryId$definition <- MedGen_def$DEF[match(MedGen_entryId$id, MedGen_def$X.CUI)]

## idNames
MedGen_idNames <- read.table(file.path(sdir,"NAMES.RGG"), sep = "|", header = TRUE, comment.char = "", quote = "", 
                         fill = TRUE, colClasses = c("character"))
MedGen_idNames$DB <- "MedGen"
MedGen_idNames <- MedGen_idNames[,c("DB","X.CUI","name")]
names(MedGen_idNames) <- c("DB","id","name")
MedGen_idNames$canonical <- TRUE

rm(MedGen_conso, MedGen_def)

###############################################################################@
## Writing tables ----
###############################################################################@
message("Writing tables...")
message(Sys.time())
toSave <- grep("^MedGen[_]", ls(), value=T)
for(f in toSave){
  message(paste("Saving", f))
  ## Ensure unicity
  assign(f, get(f))
  if(length(names(f))==0){
    f <- unique(f)
  }
  ##
  write.table(
    get(f),
    file=file.path(ddir, paste(f, ".txt", sep="")),
    sep="\t",
    row.names=FALSE, col.names=TRUE,
    quote=FALSE
  )
}
message(Sys.time())
message("... Done\n")
