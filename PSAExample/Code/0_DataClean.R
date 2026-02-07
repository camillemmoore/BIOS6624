#############################################################
## Program Name: DataClean.R							                 ##
## Purpose: To read in the dataset and label the variables ##
## 			   Create any new variables.					             ##
## Output: A Clean .csv file is exported for reading in to ##
##         all other analysis files.					             ##
## Created by:  Kevin Josey and updated by Camille Moore   ##
#############################################################

# READ IN THE RAW DATA AND CREATE VARIABLE TRANSFORMS
# Note: This is the place where I would come back to
#       for any data cleaning, outlier removal more
#       transforms.  This keeps all data cleaning   
#		    and variable decisions in one place.


# Dependencies
library(dplyr)

###You will have to set your working directory####
setwd("/Users/mooreca/Documents/BIOS6624_Repositories/PSAExample/")

#Read in the raw data
psa <- read.csv("DataRaw/prostate.csv")

colnames(psa) <- tolower(colnames(psa))

#Add indicator variable format for our categorical gleason score variable
psaclean <- psa %>%
  mutate(lpsa = log(psa),
         grade6 = ifelse(gleason == 6, 1, 0),
         grade7 = ifelse(gleason == 7, 1, 0),
         grade8 = ifelse(gleason == 8, 1, 0)) %>%
  mutate(gd6ageint = age*grade6,
         gd7ageint = age*grade7,
         gd8ageint = age*grade8) %>%
  mutate(svi_factor = ifelse(svi == 1, 'Yes', 'No'))
         
# export the clean dataset for use in other analysis programs
write.csv(psaclean, "DataProcessed/psaclean.csv", row.names = F)

# further data cleaning based on reviewing descriptive statistics
psaclean2 <- psaclean %>%
  filter(wt <= 400) %>%
  mutate(cavolsviint = cavol*svi,
         wtsviint = wt*svi,
         agesviint = age*svi,
         bphsviint = bph*svi,
         cappensviint = cappen*svi,
         gd6sviint = grade6*svi,
         gd7sviint = grade7*svi,
         gd8sviint = grade8*svi)

# export the clean dataset for use in other analysis programs
write.csv(psaclean2, "DataProcessed/psaclean2.csv", row.names = F)
