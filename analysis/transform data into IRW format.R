library(tidyr)
library(dplyr)
library(tibble)

setwd("C:/Users/jeffr/OneDrive/Desktop/stanford/CDM")

load("CDM_AlcoholApplication_Data.Rdata")

# same items used in wenchao's study including all dichotomous response
df <- data0[,41:80]

df <- tibble::rownames_to_column(df, var = "id")

# Reshape the data into long format:
df_long <- df %>%
  pivot_longer(
    cols = -id,         
    names_to = "item",  
    values_to = "response"
  )

# add Q-matrix
item_names <- colnames(df)[colnames(df) != "id"]

# Q-matrix revised by the experts
rownames(Q0) <- item_names

q_mat_df <- as.data.frame(Q0)
q_mat_df$item <- rownames(q_mat_df)

df_long_q <- df_long %>%
  left_join(q_mat_df, by = "item")

# rename columns
q_names <- c("Qmatrix_skill1", "Qmatrix_skill2", "Qmatrix_skill3", "Qmatrix_skill4")

colnames(df_long_q)[4:7] <- q_names

# output data
df <- df_long_q
save(df, file = "CDM_AlcoholApplication_Data_IRW.RData")


