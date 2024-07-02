df_svnt <- read.csv(here::here("dataraw", "nih", "paired B and sVNT_for_DHodgson.csv"))
df_svnt_long <- read.csv(here::here("dataraw", "nih", "COVAX_longitudinal_analysis.csv"))
df_bio_nih <- read.csv(here::here("dataraw","nih",  "CoVax_probe_analysis.csv"))

# Deanonymise the data for df_svnt
require(tidyverse)
pids <- df_svnt$PID 
PID_new <- pids %>% gsub("ALF", "S1", .) %>% gsub("CHW", "S2", .) %>% gsub("JHH", "S3", .) %>%  gsub("QCH", "S4", .) %>% gsub("PCH", "S5", .) %>%  gsub("WCH", "S6", .)
df_svnt$PID <- PID_new
sampleid <- df_svnt$SampleID 
SampleID_new <- sampleid %>% gsub("ALF", "S1", .) %>% gsub("CHW", "S2", .) %>% gsub("JHH", "S3", .) %>%  gsub("QCH", "S4", .) %>% gsub("PCH", "S5", .) %>%  gsub("WCH", "S6", .)
df_svnt$SampleID <- SampleID_new



# Deanonymise the data for df_svnt_long
pids <- df_svnt_long$PID 
PID_new <- pids %>% gsub("ALF", "S1", .) %>% gsub("CHW", "S2", .) %>% gsub("JHH", "S3", .) %>%  gsub("QCH", "S4", .) %>% gsub("PCH", "S5", .) %>%  gsub("WCH", "S6", .)
df_svnt_long$PID <- PID_new
sampleid <- df_svnt_long$Sample_ID 
Sample_ID_new <- sampleid%>% gsub("ALF", "S1", .) %>% gsub("CHW", "S2", .) %>% gsub("JHH", "S3", .) %>%  gsub("QCH", "S4", .) %>% gsub("PCH", "S5", .) %>%  gsub("WCH", "S6", .)
df_svnt_long$Sample_ID <- Sample_ID_new
df_svnt_long$Sample_ID <- Sample_ID_new

# Deanonymise the data for df_bio_nih
pids <- df_bio_nih$PID 
PID_new <- pids %>% gsub("ALF", "S1", .) %>% gsub("CHW", "S2", .) %>% gsub("JHH", "S3", .) %>%  gsub("QCH", "S4", .) %>% gsub("PCH", "S5", .) %>%  gsub("WCH", "S6", .)
df_bio_nih$PID <- PID_new
sampleid <- df_bio_nih$SampleID 
SampleID_new <- sampleid %>% gsub("ALF", "S1", .) %>% gsub("CHW", "S2", .) %>% gsub("JHH", "S3", .) %>%  gsub("QCH", "S4", .) %>% gsub("PCH", "S5", .) %>%  gsub("WCH", "S6", .)
df_bio_nih$SampleID <- SampleID_new

df_bio_nih_new <- df_bio_nih %>% select(!StudySite)

write.csv(df_svnt, here::here("data", "nih", "paired B and sVNT_for_DHodgson.csv"))
write.csv(df_svnt_long, here::here("data", "nih", "COVAX_longitudinal_analysis.csv"))
write.csv(df_bio_nih_new, here::here("data","nih",  "CoVax_probe_analysis.csv"))