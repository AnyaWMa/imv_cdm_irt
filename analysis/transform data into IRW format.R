library(tidyr)
library(dplyr)

setwd("C:/Users/jeffr/OneDrive/Desktop/stanford/CDM")

load("CDM_AlcoholApplication_Data.Rdata")

df <- data0[,41:80]

df <- tibble::rownames_to_column(df, var = "id")

# Reshape the data into long format:
df_long <- df %>%
  pivot_longer(
    cols = -id,         
    names_to = "item",  
    values_to = "response"
  )

save(df_long, file = "CDM_AlcoholApplication_Data_IRW.RData")

