############## Relative Read Abundance from molecular scatology ###############

## Sarah Brown | WDFW | December 2021
### Modified by: 
#### Adrianne Akmajian | Makah Fisheries Management | January 2022
#### Zoë K. Lewis | Western Washington University | February 23, 2022

#Load packages
library("plyr") 
library(tidyr)
library(tidyverse)
library(ggplot2)
library(car)

setwd()

############################### Read in data #################################

codes <- read.csv("species_codes_common_99perc.csv", 
                  header = T) #Read in species codes 

metadata <- read.csv("metadata_sealion_01122022.csv", 
                     header = TRUE) # Read in sample metadata

raw16s <- read.delim("99%_02042022_Makah Tribe Sea Lion diet_16S_Read_Number_per_Sample.txt", sep="\t", header=TRUE, stringsAsFactors = TRUE) # Read in 16S data and taxonomy assignments

rawCOI<- read.delim("99%_02042022_Makah Tribe Sea Lion diet_COIsal_Read_Number_per_Sample.txt", sep="\t", header=TRUE, stringsAsFactors = TRUE) # Read in COI data and taxonomy assignments

############################### 16s data #######################################

  ############################# Format 16s data  #############################
  ## creates a table with sample, species, and read counts 

raw16sord<- raw16s[order(raw16s$sample),] #order by sample id  

clean16s<- subset(raw16sord, speciesID !="None"& speciesID !="PHOLAR" 
                  & speciesID !="HOMSAP" & speciesID !="EUMJUB" 
                  & speciesID !="ZALCAL" & speciesID != "CALURS") 
  # remove reads with no Taxonomic ID and IDs w/ no sequences 
  # remove non-prey reads (Human, sea lion and seal DNA) 


  ####################### Loop 1:  Initial 16s percentages #####################
  ## calculates prey percentage i.e. prey species reads per total sample reads 
  ## create table with sample, species, read count, and % each species/sample 

### create empty variable
temp1 <- NULL 

### create loop for all values of Sample
for (k in levels(clean16s$sample))  
{  
  sub1 <- subset(clean16s,sample == k, 
                 select = c(sample, speciesID, freq))# subset each sample
  Totseqs1 <- sum(sub1$freq)                         # sum all seqs/sample
  sums1 <- rep(Totseqs1,length(sub1$freq))           # create vector of sums 
  mat1 <- data.frame(sub1$sample,sub1$speciesID,
                     sub1$freq,sums1)                # create dataframe with sums
  temp1 <- rbind(temp1,mat1)                         # append dataframe to temp
  print(k)}                                          # print each sample name

### calculate %species/sample 
perc1 <- (temp1$sub1.freq/temp1$sums1)*100      

### create dataframe with % species
data16s1 <- data.frame(sample =temp1$sub1.sample, 
                       species =temp1$sub1.speciesID, 
                       count = temp1$sub1.freq, perc1)  
                      

### remove all species where prey is less that 1% within a sample
filter16s1<- subset(data16s1, perc1>1, select = c(sample, species, count)) 

  ################# Loop 2: 16s percentages w/o prey species <1% ###############
  ## Recalculate percentages without low frequency species ##
  ## create table with sample, species, read count, and % of each species/sample 

temp2 <- NULL # empty variable

for (k in levels(filter16s1$sample))  #creates a loop for all values of Sample
{  
  sub2 <- subset(filter16s1,sample == k, 
                 select = c(sample, species, count)) # subset each sample
  Totseqs2 <- sum(sub2$count)                        # sum all seqs/sample
  sums2 <- rep(Totseqs2, length(sub2$count))         # create vector of sums 
  mat2 <- data.frame(sub2$sample,
                     sub2$species,sub2$count,sums2)  # create dataframe with sums
  temp2 <- rbind(temp2,mat2)                         # append dataframe to temp
  print(k)}                                          # print sample name 

### calculate new percentages of each species per sample
perc16s2 <- (temp2$sub2.count/temp2$sums2)*100                                     

### create dataframe including percentages calculated without species <%1
data16s2 <- data.frame(sample =temp2$sub2.sample, species =temp2$sub2.species, 
                       count = temp2$sub2.count, perc16s2)

  ###################### 16s % data with metadata ##############################

### add species labels
final16s <- merge(data16s2, codes, by = "species", sort = F) 

### merge metadata and sample averages
merged16sdata <- merge(final16s,metadata,
                       by.x ="sample", 
                       by.y = "WDFW_code") 

write.csv(merged16sdata,"SL_16S_DNA_sample_averages_ZKL.csv")

############################## COI prey percentage#############################
  
  ############################### Format COI data ###############################
  ## creates a table with sample, species, and read counts 
  
  rawCOItab <- as.data.frame(rawCOI) # Tabulates sequence Ids 
  rawCOIord<- rawCOI[order(rawCOI$sample),] #order by sample id  
  cleanCOI<- subset(rawCOIord, speciesID !="None"& speciesID !="PHOLAR" 
                    & speciesID !="HOMSAP" & speciesID !="EUMJUB" 
                    & speciesID !="ZALCAL" & speciesID != "CALURS") 
  # remove reads with no Taxonomic ID and IDs w/ no sequences 
  # remove non-prey reads (Human, sea lion and seal DNA) 

  ####################### Loop 3:  Initial COI percentages #####################
  ## calculates prey percentage i.e. prey species reads per total sample reads 
  ## create table with sample, species, read count, and % each species/sample 

### create empty variable
temp3 <- NULL 

### create loop for all values of Sample
for (k in levels(cleanCOI$sample))  
{  
  sub3 <- subset(cleanCOI,sample == k, 
                 select = c(sample, speciesID, freq))# subset each sample
  Totseqs3 <- sum(sub3$freq)                         # sum all seqs/sample
  sums3 <- rep(Totseqs3,length(sub3$freq))           # create vector of sums 
  mat3 <- data.frame(sub3$sample,sub3$speciesID,
                     sub3$freq,sums3)                # create dataframe with sums
  temp3 <- rbind(temp3,mat3)                         # append dataframe to temp
  print(k)}                                          # print each sample name

### calculate %species/sample 
percCOI1 <- (temp3$sub3.freq/temp3$sums3)*100           

### create dataframe with % species
dataCOI1 <- data.frame(sample =temp3$sub3.sample, 
                       species =temp3$sub3.speciesID, 
                       count = temp3$sub3.freq, percCOI1)

#remove all species where prey is less that 1% within a sample
filterCOI1<- subset(dataCOI1, percCOI1>1, select = c(sample, species, count)) 


  ################# Loop 4: COI percentages w/o prey species <1% ###############
  ## Recalculate percentages without low frequency species ##
  ## create table with sample, species, read count, and % of each species/sample 

temp4 <- NULL # empty variable

for (k in levels(filterCOI1$sample))  #creates a loop for all values of Sample
{  
  sub4 <- subset(filterCOI1,sample == k, 
                 select = c(sample, species, count)) # subset each sample
  Totseqs4 <- sum(sub4$count)                        # sum all seqs/sample
  sums4 <- rep(Totseqs4, length(sub4$count))         # create vector of sums 
  mat4 <- data.frame(sub4$sample,
                     sub4$species,sub4$count,sums4)  # create dataframe with sums
  temp4 <- rbind(temp4,mat4)                         # append dataframe to temp
  print(k)}                                          # print sample name 

### calculate new percentages of each species per sample
percCOI2 <- (temp4$sub4.count/temp4$sums4)*100          

### create dataframe including percentages calculated without species <%1
dataCOI2 <- data.frame(sample =temp4$sub4.sample, 
                       species =temp4$sub4.species, 
                       count = temp4$sub4.count, 
                       percCOI2)

###################### COI % data with metadata ##############################

### add species labels
finalCOI <- merge(dataCOI2, codes, by = "species", sort = F)

### merge metadata and sample averages
mergedCOIdata <- merge(finalCOI,metadata,
                       by.x ="sample", 
                       by.y = "WDFW_code") #merge metadata and sample averages
  
  write.csv(mergedCOIdata,"SL_COI_DNA_sample_averages_ZKL.csv")
  
###################### 16s scaled to COI data ##############################
 
  rawdna16 <- read.csv('SL_16S_DNA_sample_averages_ZKL.csv')
  rawdna16$sample <- rawdna16$sample.y
dna16 <- rawdna16 %>% select(-sample.y)
  
  rawdnaCoi <- read.csv('SL_COI_DNA_sample_averages_ZKL.csv')
  rawdnaCoi$sample <- rawdnaCoi$sample.y
dnaCoi <- rawdnaCoi %>% select(-sample.y)
  
#### May choose to subset only portion of the data #########3
  
  dna16.10plus <- subset(dna16,
                         tapply(count , sample , sum , na.rm=TRUE)[sample] >= 10)
  dna16.10plus$species <- factor(dna16.10plus$species)
  dna16.10plus$sample <- factor(dna16.10plus$sample)
  nlevels(dna16$sample) - nlevels(dna16.10plus$sample)
  # Prints number of samples removed if greater than 10 species
  
  # Salmon spp
  sspp <- c('ONCGOR', 'ONCOR', 'ONCCLA', 'ONCKET', 'ONCNER', 'ONCTSH', 'ONCMYK', 'ONCKIS', 'SALSAL', 'SALALP', 'SALMAL','SALTRU')
  dna16.10plus$Sspp <- dna16.10plus$species %in% sspp
  with(dna16.10plus , table(species,Sspp,useNA='ifany'))
  
  # is there Salmon in sample
  SlmSamp.temp <- with(dna16.10plus,
                       tapply(Sspp , sample , function(x) sum(x , na.rm = TRUE) > 0 )
  )
  table(SlmSamp.temp) # X samples with no Salmon
  dna16.10plus$SlmSamp <- SlmSamp.temp[as.character(dna16.10plus$sample)]
  
  # proportion of Sspp in sample
  samsampprop.temp <- by(dna16.10plus , dna16.10plus$sample, function(x) sum(x$perc16s2[x$Sspp])/100)
  dna16.10plus$samsampprop <- samsampprop.temp[as.character(dna16.10plus$sample)]
  
  ############## done with 16
  
  # Subset the Coi data set >= 5 AND only for those samples with Salmon in 16S
  dnaCoi.10plus <- subset(dnaCoi,
                          (tapply(count , sample , sum , na.rm=TRUE)[as.character(sample)] >= 5) & SlmSamp.temp[as.character(sample)])
  dnaCoi.10plus$species <- factor(dnaCoi.10plus$species)
  dnaCoi.10plus$sample <- factor(dnaCoi.10plus$sample)
  nlevels(dnaCoi$sample) - nlevels(dnaCoi.10plus$sample) # Prints number samples filtered
  str(dnaCoi.10plus) 
  
  ## rescale % in Coi to 16
  dnaCoi.10plus$precT <- samsampprop.temp[as.character(dnaCoi.10plus$sample)] * dnaCoi.10plus$percCOI2/100
  
  ## put 16 % on 0-1 scale
  dna16.10plus$precT <- dna16.10plus$perc16s2/100
  
  # merge
  dna.new.10plus <- merge(dna16.10plus , dnaCoi.10plus,
                          by =  c("sample"  ,  "species"  , "day" , "month"   ,  "year"   , "site" ), all = TRUE)
  
  # condition "is the sample in Coi ?"
  cond <- as.character(dna.new.10plus$sample) %in% levels(dnaCoi.10plus$sample)
  table(cond,useNA='ifany')
  
  # apply pcnt as a heirarchy of conditions
  dna.new.10plus$PRCNT <- with(dna.new.10plus,
                               ifelse(!cond , precT.x , # if no Coi for sample - prec16
                                      ifelse( !(is.na(Sspp) | Sspp)  , precT.x , # if Coi but not Sspp or extra Sspp - prec16
                                              ifelse(is.na(precT.y),0,precT.y) # if Sspp not in Coi, 0, otherwise Coi
                                      )
                               )
  )
  # check all PRCNT sum to 1 within a sample
  summary(aggregate(PRCNT ~ sample , sum , data = dna.new.10plus))
  
  # order dna new and remove the extra columns
  dna.new <- dna.new.10plus[,!(names(dna.new.10plus) %in% c('X.x','X.y','X.1.x','X.1.y'))]
  dna.new <- dna.new[order(dna.new$sample,dna.new$Sspp,dna.new$species),]
  
  
  dna.new <- subset(dna.new,PRCNT != 0)   # remove zero diet values
  
  dna.new$PRCNT100 <- dna.new$PRCNT*100
  
  ##write.csv(dna.new , 'dna_new.csv' , row.names= FALSE)
  
  
  ###################################################################################################################
  #Lastly, merge the data to a final spreadsheet
  #replace .x with "16S" and .y with "COI"
  
  merged_data <- merge(dna.new, codes, by = "species", sort = F)
  ##write.csv(merged_data,"SSL_Final_DNA_sample_averages_long.csv")
  
##### COnvert to wide to add zeros for RRA data
  
  RRA <- merged_data  %>% select(species, sample, PRCNT100)
  
  RRA.wide <- pivot_wider(RRA, names_from=sample, values_from=PRCNT100) %>%  
    replace(is.na(.), 0)

  RRA.long <- gather(RRA.wide, sample, PRCNT100, EJ_001:PV_012, factor_key=TRUE)
  
  mergedRRA <- merge(RRA.long,metadata,
                        by.x ="sample")

  write.csv(mergedRRA, "SSL_RRA.csv")  
  
  