rm(list = ls())
gc()

# setwd("~/Shared/Data-Science/Data-Source-Model-Repository/MedGen/scripts/")

library(here)
library(XML)
library(parallel)
library(here)
source(here("..", "00-Utils/writeLastUpdate.R"))
library(dplyr)
library(tidyr)
library(tibble)
library(ReDaMoR)
##
mc.cores <- 55
sdir <- "~/Shared/Data-Science/Data-Source-Model-Repository/MedGen/sources"
ddir <- here("data")

###############################################################################@
## Data model ----
###############################################################################@
load(here("model", "MedGen.rda"))
# dm <- model_relational_data()
save(dm, file = here("model", "MedGen.rda"))

###############################################################################@
## Source information ----
###############################################################################@

sfi <- read.table(
   file.path(sdir, "ARCHIVES/ARCHIVES.txt"),
   sep="\t",
   header=T,
   stringsAsFactors=FALSE
)
MedGen_sourceFiles <- sfi[which(sfi$inUse), c("file", "url", "current")]

###############################################################################@
## Data from medgen ----
###############################################################################@
## decompress gz
for(f in MedGen_sourceFiles$file){
  gzf <- file.path(sdir,f)
  print(gzf)
  system(paste0("gzip -df ", gzf))
  # system(paste0("rm ",file.path(sdir,gsub(".gz","",f))))
}

###########################
## Basic information

## crossId
MedGen_conso <- read.table(file.path(sdir,"MGCONSO.RFF"), sep = "|", header = TRUE, comment.char = "", quote = "", 
                      fill = TRUE, colClasses = c("character"))
table(MedGen_conso$SAB)
MedGen_conso$SAB[MedGen_conso$SAB == "MSH"] <- "MeSH"
MedGen_conso$SAB[MedGen_conso$SAB == "NCI"] <- "NCIt"
MedGen_conso$SAB[MedGen_conso$SAB == "SNOMEDCT_US"] <- "SNOMEDCT"
MedGen_conso$SAB[MedGen_conso$SAB == "ORDO"] <- "ORPHA"
MedGen_conso$SAB[MedGen_conso$SAB == "HPO"] <- "HP"
MedGen_conso <- MedGen_conso[!MedGen_conso$SAB == "",]
table(MedGen_conso$SAB)
##
MedGen_conso$DB <- "MedGen"
MedGen_conso$SDUI <- gsub(".*:","",MedGen_conso$SDUI)
MedGen_conso$SDUI <- gsub(".*_","",MedGen_conso$SDUI)
MedGen_conso$id <- ifelse(MedGen_conso$SAB %in% c("NCIt","SNOMEDCT"),MedGen_conso$SCUI,MedGen_conso$SDUI)
MedGen_conso$canonical <- FALSE

MedGen_crossId <- unique(MedGen_conso[,c("DB","X.CUI","SAB","id")])
names(MedGen_crossId) <- c("DB1","id1","DB2","id2")
table(MedGen_crossId$DB1)
table(MedGen_crossId$DB2)
dim(MedGen_crossId)

## When removing prefix, an integer is a correct disease ID
table(!is.na(as.numeric(sub("^[^[:digit:]]*", "", MedGen_crossId$id2))))
table(!is.na(as.numeric(sub("^[^[:digit:]]*", "", MedGen_crossId$id1))))
toKeep <- MedGen_crossId[which(!is.na(as.numeric(sub("^[^[:digit:]]*", "", MedGen_crossId$id2))) &
                                !is.na(as.numeric(sub("^[^[:digit:]]*", "", MedGen_crossId$id1)))),]
dim(toKeep)
toCheck <- MedGen_crossId[-which(!is.na(as.numeric(sub("^[^[:digit:]]*", "", MedGen_crossId$id2))) &
                            !is.na(as.numeric(sub("^[^[:digit:]]*", "", MedGen_crossId$id1)))),]
dim(toCheck)

## Remove any DBs that are not disease DBs and DB1 can only be "EFO" or "Orphanet"
## check wrong IDs, remove weird ones still
table(toKeep$DB2)
table(toKeep$DB1)

## Remove self references
MedGen_crossId[which(MedGen_crossId$id1 == MedGen_crossId$id2),]
## MedGen_crossId <- MedGen_crossId[-which(MedGen_crossId$id1 == MedGen_crossId$id2),]

## Keep only disease
MedGen_dis <- read.table(file.path(sdir,"MGSTY.RFF"), sep = "|", header = TRUE, comment.char = "", quote = "", 
                         fill = TRUE, colClasses = c("character"))
MedGen_type <- MedGen_dis %>% 
  mutate(db = "MedGen") %>% 
  select(db, 
         id = X.CUI,
         type = STY)
# MedGen_dis <- MedGen_dis[grep(paste("Disease or Syndrome","Acquired Abnormality",
#                                     "Anatomical Abnormality","Congenital Abnormality",sep = "|"),
#                               MedGen_dis$STY),]
MedGen_crossId <- MedGen_crossId[which(MedGen_crossId$id1 %in% MedGen_dis$X.CUI),]
dim(MedGen_crossId)

## idNames
MedGen_idNames <- read.table(file.path(sdir,"NAMES.RGG"), sep = "|", header = TRUE, 
                             comment.char = "", quote = "", 
                             fill = TRUE, colClasses = c("character")) %>%
  mutate(DB = "MedGen",
         canonical = TRUE) %>%
  select(DB, 
         id = X.CUI, 
         syn = name, 
         canonical)
# MedGen_idNames$DB <- "UMLS"
# MedGen_idNames <- MedGen_idNames[,c("DB","X.CUI","name")]
# MedGen_idNames$canonical <- FALSE
##
MedGen_idNames <- bind_rows(MedGen_idNames,
                            MedGen_conso %>%
                              select(DB, 
                                     id  = X.CUI, 
                                     syn = STR, 
                                     canonical)) %>%
  distinct()
dim(MedGen_idNames)
## * Finding non-numeric ID 
id <-  (MedGen_idNames$id)
numidsuf <- as.numeric(sub("^[^[:digit:]]*", "", id))
MedGen_idNames[which(is.na(numidsuf)),]
# MedGen_idNames <- MedGen_idNames[-which(is.na(numidsuf)),] 

## Check characters for \t, \n, \r and put to ASCII
MedGen_idNames$syn <- iconv(x = MedGen_idNames$syn,to="ASCII//TRANSLIT")
MedGen_idNames$syn <- gsub(paste("\n","\t","\r", sep = "|")," ",MedGen_idNames$syn)
MedGen_idNames$syn <- gsub("\"","'",MedGen_idNames$syn)
table(unlist(sapply(MedGen_idNames$syn, strsplit, split = "")))

## Remove NA
table(is.na(MedGen_idNames$syn))
MedGen_idNames <- MedGen_idNames[!is.na(MedGen_idNames$syn),]
## Remove duplicated (keep canonical)
dim(MedGen_idNames)
dim(unique(MedGen_idNames))
# MedGen_idNames <- MedGen_idNames[order(MedGen_idNames$canonical,decreasing = T),]
# MedGen_idNames <- unique(MedGen_idNames)
dim(MedGen_idNames)

## all idNames in entryId
table(MedGen_idNames$id %in% MedGen_dis$X.CUI)
MedGen_idNames <- MedGen_idNames[which(MedGen_idNames$id %in% MedGen_dis$X.CUI),]

## Remove empty names, ifany
nc <- nchar(MedGen_idNames$syn)
head(table(nc))
MedGen_idNames[which(nc < 4),]
## Remove names of 0 or 1 character long
MedGen_idNames[which(nc == 0),]
MedGen_idNames[which(nc == 1),]
# MedGen_idNames <- MedGen_idNames[-which(nc == 1),]

## entryId
MedGen_entryId <- unique(data.frame(DB = c(MedGen_crossId$DB1,MedGen_idNames$DB),
                             id = c(MedGen_crossId$id1,MedGen_idNames$id),
                             stringsAsFactors = F))
MedGen_def <- read.table(file.path(sdir,"MGDEF.RFF"), sep = "|", header = TRUE, comment.char = "", quote = "", 
                           fill = TRUE, colClasses = c("character"))
MedGen_entryId$def <- MedGen_def$DEF[match(MedGen_entryId$id, MedGen_def$X.CUI)]
## Check characters for \t, \n, \r and put to ASCII
MedGen_entryId$def <- iconv(x = MedGen_entryId$def,to="ASCII//TRANSLIT")
MedGen_entryId$def <- gsub(paste("\n","\t","\r", sep = "|")," ",MedGen_entryId$def)
MedGen_entryId$def <- gsub("\"","'",MedGen_entryId$def)
MedGen_entryId$def <- gsub("\\\\","",MedGen_entryId$def)
MedGen_entryId$level <- NA
table(unlist(sapply(MedGen_entryId$def, strsplit, split = "")))

## check NA
table(is.na(MedGen_entryId$def))
## Remove duplicated (keep canonical)
dim(MedGen_entryId)
dim(unique(MedGen_entryId))

## all idNames in entryId
table(MedGen_entryId$id %in% MedGen_dis$X.CUI)
# MedGen_entryId <- MedGen_entryId[which(MedGen_entryId$id %in% MedGen_dis$X.CUI),]
## Remove empty names, ifany
nc <- nchar(MedGen_entryId$def)
head(table(nc))
MedGen_entryId[which(nc < 3),]
## Remove names of 0 or 1 character long
MedGen_entryId[which(nc == 0),]
MedGen_entryId[which(nc == 1),]
MedGen_entryId$def[which(nc < 3)] <- NA

## If definition is NA, label is used (canonical)
tmp <- MedGen_idNames %>% filter(canonical)
MedGen_entryId <- MedGen_entryId %>% 
  mutate(def = case_when(is.na(def) ~ tmp$syn[match(id,tmp$id)],
                         TRUE ~ def))

## Not all ids have a canonical label, use the non-canonical one
tmp <- MedGen_idNames %>% filter(id %in% MedGen_entryId[is.na(MedGen_entryId$def),"id"])
MedGen_entryId <- MedGen_entryId %>% 
  mutate(def = case_when(is.na(def) ~ tmp$syn[match(id,tmp$id)],
                         TRUE ~ def))
table(MedGen_type$id %in% MedGen_entryId$id)

rm(MedGen_conso, MedGen_def, MedGen_dis)

## 
MedGen_HPO <- read.table(file.path(sdir,"MedGen_HPO_Mapping.txt"), sep = "|", header = TRUE, comment.char = "", quote = "", 
                         fill = TRUE, colClasses = c("character")) %>%
  filter(X.CUI %in% MedGen_entryId$id) %>%
  select(id = X.CUI,
         hp = SDUI) %>%
  mutate(DB = "MedGen", 
         hp = gsub(".*:", "", hp))
table(MedGen_HPO$id %in% MedGen_entryId$id)

head(MedGen_HPO)

MedGen_OMIM_HPO <- read.table(file.path(sdir,"MedGen_HPO_OMIM_Mapping.txt"), sep = "|", header = TRUE, comment.char = "", quote = "", 
                         fill = TRUE, colClasses = c("character")) 
MedGen_OMIM_entryId <- MedGen_OMIM_HPO %>%
  mutate(DB = "OMIM",
         level = NA) %>%
  select(DB,
         id = MIM_number,
         def = OMIM_name,
         level) %>%
  distinct()
MedGen_OMIM_entryId$def <- iconv(x = MedGen_OMIM_entryId$def,to="ASCII//TRANSLIT")
MedGen_OMIM_entryId$def <- gsub(paste("\n","\t","\r", sep = "|")," ",MedGen_OMIM_entryId$def)
MedGen_OMIM_entryId$def <- gsub("\"","'",MedGen_OMIM_entryId$def)
# table(unlist(sapply(MedGen_OMIM_entryId$def, strsplit, split = "")))

MedGen_OMIM_idNames <- MedGen_OMIM_entryId %>%
  select(DB, 
         id, 
         syn = def) %>%
  mutate(canonical = TRUE)
MedGen_OMIM_HPO <- MedGen_OMIM_HPO %>%
  mutate(DB = "OMIM", 
         phenoDB = "HP") %>%
  select(DB, 
         id = MIM_number,
         hp = HPO_ID)
head(MedGen_OMIM_HPO)

table(MedGen_OMIM_HPO$id %in% MedGen_OMIM_entryId$id)


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
    quote=TRUE,
    qmethod = "double"
  )
}
message(Sys.time())
message("... Done\n")

##############################################################
## Check model
# source("../00-Utils/autoCheckModel.R")
