library(RSQLite)
library(reshape2)
library(FactoMineR)
library(factoextra)
library(corrplot)
library(psych)
library(dplyr)
library(DescTools)
library(tidyr)
library(waffle)
library(ggplot2)
library(ggpubr)
library(Rlof)

#********************************************************************************************************************************************************************************************************
# Code to investigate business model change (BMC) activity in firms. It imports historical financial data from the project's SQLite database, uses dimensionality reduction with a collection
# of variables, uses outlier detection to create a measure of BMC, and test association of BMC activity patterns and firm characteristics
# STEP 1:
#   1A: Get data from db, re-sample to obtain all firms with full 12 years of data (with no NA), and prepare data for MFA in next step
#   1B: Inspect the data for systematic missing values, systematic zeros and extreme outliers, to avoid issues with factor analysis).  Remove variables and treat extreme outlier data points accordingly
# STEP 2: Run MFA (multiple factor analysis, a weigthed PCA across multiple timepoints as dimensions) for all ratios. This is to choose ratios per value creation, delivery and capture dimensions
# STEP 3: Run LOF (local outlier factor analysis) per company
# STEP 4: Label change events with "1" (per value dimension that is) using a cut-off that I need to calculate based on percentiles
# STEP 5: Aggregate events by BM across 1, 2 or 3 years, and laben the BMC with intensity 1 if done in 3yrs, intensity 2 if done in 2yrs, or intensity 3 if done in 3 yrs
# STEP 6: Based on BMCs per firm over time from step 5, group firms into 4 BMC pattern groups (based on intensity=mean intensity over 12yrs and frequency=no.BMCs over 12 yrs): 
#         LoFre-LoMag; LoFre-HiMag; HiFre-LoMag; HiFre-HiMag
# STEP 7: For each firm, get the:
#     7A: Mean of Log(Total assets) over the 12 years, and group firms into 3 size groups: Lo; Mid; Hi
#     7B: Industry (6-digit GICS code), and group firms by them
#     7C: Age (year of incorporation), and group firms into 3 age groups
# STEP 8: Run Chi-Square test of independence using the groups in STEP 6 and STEP 7
# STEP 9: Visualization of no. value dim change events, and no. BMC events (by intensity)
#********************************************************************************************************************************************************************************************************

#===STEP 1 (import, inspect(1A), filter(1B) and prepare data for MFA) ==============================================================================================================================================
#get imputed 2006+ data from db (1A):
dbconn <- dbConnect(dbDriver("SQLite"), dbname = "~/Documents/SQLite/BMC_activity_data.sqlite")
query_1 <- "SELECT ID_company, parameter, timepoint, value FROM ratios_imputed_table WHERE timepoint  >= 2006 AND value IS NOT NULL"
ratios_raw <- data.frame(dbGetQuery(dbconn, query_1), stringsAsFactors = F)

query_2 <- "SELECT ID_company, parameter, timepoint, value FROM financials_imputed_table WHERE parameter = 'Tota' AND timepoint  >= 2006 AND value IS NOT NULL"
financials_raw <- data.frame(dbGetQuery(dbconn, query_2), stringsAsFactors = F)

query_3 <- "SELECT ID_company, name, gics_code, incorp_date FROM companies_table"
companies_raw <- data.frame(dbGetQuery(dbconn, query_3), stringsAsFactors = F)

dbDisconnect(dbconn)

#merge all:
raw <- rbind(ratios_raw, financials_raw)
raw <- merge(raw, companies_raw, by = "ID_company")
raw <- raw[,c(1,5,6,7,2,3,4)]

#Inspect data / data wrangling (1B):
temp <- raw
temp$ID_company_timep <- paste(temp$ID_company, temp$timepoint, sep = "_") # new column with ID company and timepoint
temp2 <- dcast(temp, ID_company_timep ~ parameter, value.var = "value") #transform data from long to wide format
temp3 <- data.frame(psych::describe(temp2[,-1])) # descriptive analysis
temp3$vars <- row.names(temp3)
#additional descriptives for NAs, zeroes, etc:
no_obs <- length(unique(temp$ID_company))*12
temp4 <- raw %>% 
  group_by(parameter) %>% 
  summarise(n = n(), 
            NAs = no_obs-n(), 
            NAs_perc = (no_obs-n())*100/no_obs, 
            Zeros = sum(value==0), 
            Zeros_perc = sum(value==0)*100/no_obs)
colnames(temp4)[1] <- "vars"
raw_summary <- merge(temp4, temp3[,-2], by = "vars")
rm(temp,temp2,temp3,temp4) #clear temporary dataframes from memory

#Result:
#rd.sa and rd.ot have too many zeros (%25), while others have >3%. This means contribution to these variable is significantly lower --> Drop them for consistency.
data <- raw[raw$parameter!="rd.sa" & raw$parameter!="rd.ot",] #remove rd.sa and rd.ot
#Also, there are extreme outliers in most variables --> Winsorise all variables
data <- data %>% 
  dplyr::group_by(parameter) %>% 
  mutate(value = DescTools::Winsorize(value, probs = c(0.01, 0.99), na.rm = T)) #winzorize extreme outliers per variable (1% level)

#Inspect data again:
temp <- data
temp$ID_company_timep <- paste(temp$ID_company, temp$timepoint, sep = "_") # new column with ID company and timepoint
temp2 <- dcast(temp, ID_company_timep ~ parameter, value.var = "value") #transform data from long to wide format
temp3 <- data.frame(psych::describe(temp2[,-1])) # descriptive analysis
temp3$vars <- row.names(temp3)
#additional descriptives for NAs, zeroes, etc:
no_obs <- length(unique(temp$ID_company))*12
temp4 <- raw %>% 
  group_by(parameter) %>% 
  summarise(n = n(), 
            NAs = no_obs-n(), 
            NAs_perc = (no_obs-n())*100/no_obs, 
            Zeros = sum(value==0), 
            Zeros_perc = sum(value==0)*100/no_obs)
colnames(temp4)[1] <- "vars"
raw_summary <- merge(temp4, temp3[,-2], by = "vars")
rm(temp,temp2,temp3,temp4) #clear temporary dataframes from memory
#...Result: skewness significantly reduced (from ~100 in some cases to ~5)
#===End STEP 1=================================================================================================================================================================================


#===STEP 2 (run MFA and identify ratios) ======================================================================================================================================================
# Useful MFA reference: https://stats.stackexchange.com/questions/18617/can-i-do-a-pca-on-repeated-measures-for-data-reduction
df_mfa <- data[data$parameter != "Tota",] #filter out total assets data, leave only bmc financial ratios
df_mfa$parameter_timepoint <- paste0(df_mfa$parameter, df_mfa$timepoint)
df_mfa <- dcast(df_mfa, ID_company ~ parameter_timepoint, value.var = "value") #transform from long to wide
df_mfa <- df_mfa[complete.cases(df_mfa),] #Remove firms with missing values, because some variables with more NAs might be unfavoured over others with less NAs, affecting their contribution to the dimensions)

#calculate MFA with 10 underlying dimensions, so then I can choose between them
mfa <- FactoMineR::MFA(df_mfa[,-1], group=rep((2017-min_year+1),21), ncp = 10, graph = F,
                       name.group=c("cg.in","cg.op","co.pe","de.op","de.sa","em.ca","fi.ot","fr.op","ia.ta","oi.fa","ot.em","ro.ic","sa.ap","sa.ar","sa.nu","sa.ta","se.co","se.sa","st.tu","ta.ia","ta.ta"))
#visual analytics:
fviz_screeplot(mfa)
fviz_mfa_var(mfa, palette="igv","group", axes = c(1,2))
fviz_mfa_var(mfa, palette="igv","group", axes = c(2,3))

fviz_contrib(mfa, choice = "group", axes = 1, top = 10) #dimension /princ comp 1
fviz_contrib(mfa, choice = "group", axes = 2, top = 10) #dimension /princ comp 2
fviz_contrib(mfa, choice = "group", axes = 3, top = 10) #dimension /princ comp 3
write.csv(data.frame(mfa$group$contrib),"~/Desktop/variable_contrib_to_dims.csv", row.names = T) #export contributions (in %) of each ratio to the underlying dimensions

# Result (after assessing variable_contrib_to_dims.csv and visual analytics):
  #The choice of ratios is:
    # Vcre: de.op, de.sa and ta.ta (combined contribution of 55%. Individual contributors contributing to more than 10%, and none of them contributes more than 4% to other dimensions)
    # Vdel: cg.op, ot.em and se.co (combined contribution of 40.5%. Individual contributors contributing to more than 10%, and none of them contributes more than 3% to other dimensions)
#    Vcap: co.pe and sa.ar (combined contribution of 58%. Individual contributors contributing to more than 10%, and none of them contributes more than 4% to other dimensions)
#===End STEP 2=================================================================================================================================================================================




#===STEP 3 (run LOF and store LOF values per company, timepoint and variable) =================================================================================================================
list_of_companies <- unique(df_mfa$ID_company) # this is the final sample -  firms with no NAs (the same used for MFA)

#Prepare datasets for LOF - transform original data (with variables selected in STEP 2) from long to semi-wide (ID_company and timepoints in separated columns):
data_vcre <- data[data$parameter %in% c("de.op", "de.sa", "ta.ta"), ]
data_vcre <- dcast(data_vcre, ID_company + timepoint ~ parameter, value.var = "value")
data_vdel <- data[data$parameter %in% c("cg.op", "ot.em", "se.co"), ]
data_vdel <- dcast(data_vdel, ID_company + timepoint ~ parameter, value.var = "value")
data_vcap <- data[data$parameter %in% c("co.pe", "sa.ar"), ]
data_vcap <- dcast(data_vcap, ID_company + timepoint ~ parameter, value.var = "value")
result <- data.frame(c()) # this will store resulting LOF values
pb <- txtProgressBar(min = 0, max = length(list_of_companies), initial = 0) #progress bar

#loop through individual companies, get data from previously built datasets and estimate LOF:
# Note on choice of LOF"s K: After inspection of ~6 random cases, K=3 seems to be the most accurate, but sometimes is 2, and also 4. To aggregate them, I tested with mean of the three, median, and max. Median is the most accurate.
for (i in 1:length(list_of_companies)) {
  temp_vcre <- data_vcre[data_vcre$ID_company==list_of_companies[i], -1]
  res_lof_vcre <- data.frame(Rlof::lof(temp_vcre[,-1], 2:4)) #LOF calculation with a range of Ks (from 2-4)
  res_lof_vcre <- rapply(res_lof_vcre, f=function(x) ifelse(is.nan(x) | is.infinite(x),0,x), how="replace") # replace "Inf" and "NaN" by 0, so that further calculation of change events is not affected
  temp_vcre$vcre_LOF <- apply(res_lof_vcre[1:ncol(res_lof_vcre)],1,function(x) median(as.numeric(x),na.rm = T)) #Median of all LOFs (calculated with different Ks). This is per timepoint
  
  temp_vdel <- data_vdel[data_vdel$ID_company==list_of_companies[i], -1]
  res_lof_vdel <- data.frame(Rlof::lof(temp_vdel[,-1], 2:4)) #LOF calculation with a range of Ks (from 2-4)
  res_lof_vdel <- rapply(res_lof_vdel, f=function(x) ifelse(is.nan(x) | is.infinite(x),0,x), how="replace") # replace "Inf" and "NaN" by 0, so that further calculation of change events is not affected
  temp_vdel$vdel_LOF <- apply(res_lof_vdel[1:ncol(res_lof_vdel)],1,function(x) median(as.numeric(x),na.rm = T)) #Median of all LOFs (calculated with different Ks). This is per timepoint
  
  temp_vcap <- data_vcap[data_vcap$ID_company==list_of_companies[i], -1]
  res_lof_vcap <- data.frame(Rlof::lof(temp_vcap[,-1], 2:4)) #LOF calculation with a range of Ks (from 2-4)
  res_lof_vcap <- rapply(res_lof_vcap, f=function(x) ifelse(is.nan(x) | is.infinite(x),0,x), how="replace") # replace "Inf" and "NaN" by 0, so that further calculation of change events is not affected
  temp_vcap$vcap_LOF <- apply(res_lof_vcap[1:ncol(res_lof_vcap)],1,function(x) median(as.numeric(x),na.rm = T)) #Median of all LOFs (calculated with different Ks). This is per timepoint
  
  temp1 <- merge(temp_vcre[,c(1,ncol(temp_vcre))], temp_vdel[,c(1,ncol(temp_vdel))], by = "timepoint")
  temp2 <- merge(temp1, temp_vcap[,c(1,ncol(temp_vcap))], by = "timepoint")
  temp2$ID_company <- list_of_companies[i]
  result <- rbind(result, temp2[,c(5,1,2:4)])
  setTxtProgressBar(pb,i)
}
#===End STEP 3=================================================================================================================================================================================




#===STEP 4 (label change events with 1 based on cut-off value) ===============================================================================================================================
# Notes: Cut-offs calculated at the full-sample level, not sub-industry. How percentage values (i.e. 66%, 95% or 99%) are chosen?
#   - upper bound at 99% allows to remove extremely outlying data points potentially resulting from incongruencies in data reported
#   - lower bound at 66% means that significant changes in value dimensions associated with BMC are those in the upper third of the distribution - consistent with org. change literature
#   - smaller range for value delivery allows to maintain symmetry in the resulting number of change events across dimensions
lb_vcre <- as.numeric(quantile(result$vcre_LOF, .66, na.rm = T))
ub_vcre <- as.numeric(quantile(result$vcre_LOF, .99, na.rm = T))

lb_vdel <- as.numeric(quantile(result$vdel_LOF, .66, na.rm = T))
ub_vdel <- as.numeric(quantile(result$vdel_LOF, .95, na.rm = T))

lb_vcap <- as.numeric(quantile(result$vcap_LOF, .66, na.rm = T))
ub_vcap <- as.numeric(quantile(result$vcap_LOF, .99, na.rm = T))

#Label change events per value dimension. 1 if falls between the lower and upper bounds, 0 otherwise :
result$vcre_change <- ifelse(result$vcre_LOF >= lb_vcre & result$vcre_LOF <= ub_vcre, 1, 0)
result$vdel_change <- ifelse(result$vdel_LOF >= lb_vdel & result$vdel_LOF <= ub_vdel, 1, 0)
result$vcap_change <- ifelse(result$vcap_LOF >= lb_vcap & result$vcap_LOF <= ub_vcap, 1, 0)
#===End STEP 4 ================================================================================================================================================================================




#===STEP 5 (Aggregate events at the BM level) =================================================================================================================================================
# Notes: BMC events aggregated in time periods of:
#     - 1 year (fast changes where all dims are changed simultaneously)
#     - 2 years (semi-fast changes where two dims are changed in a single year, and the remaining one the next year)
#     - 3 years (slow changes where each individual dim is changed per year)

dtfm <- result[,c(1,2,6,7,8)] #get only change event columns

#-1yr changes (speed=3):
dtfm_result <- data.frame(ID_company = as.character(), timepoint = as.character(), intensity = as.numeric(), stringsAsFactors = F)
for (c in 1:length(list_of_companies)) { #treat one company at a time
  dtfm_temporal <- dtfm[dtfm$ID_company==list_of_companies[c],]
  for (t in 2:(nrow(dtfm_temporal)-1)) { # this range avoids first and last rows, so that it doesn't throw error (reduces N but it guarantees the events within the time window are not part of a larger change outside the time window)
    vcre_pre <- dtfm_temporal[t-1, 3]
    vdel_pre <- dtfm_temporal[t-1, 4]
    vcap_pre <- dtfm_temporal[t-1, 5]
    vcre <- dtfm_temporal[t, 3]
    vdel <- dtfm_temporal[t, 4]
    vcap <- dtfm_temporal[t, 5]
    vcre_post <- dtfm_temporal[t+1, 3]
    vdel_post <- dtfm_temporal[t+1, 4]
    vcap_post <- dtfm_temporal[t+1, 5]
    #if ((vcre_pre+vdel_pre+vcap_pre==0) & (vcre_post+vdel_post+vcap_post==0)) { #if it's bounded by 0's, then...
    #case 1:
    if (vcre>=1 & vdel>=1 & vcap>=1) { dtfm_result[nrow(dtfm_result)+1, ] <- list(dtfm_temporal[t,1], dtfm_temporal[t,2], 3) }
    #}
  }
}
#-end

#-2yr changes (speed=2):
dtfm_result2 <- data.frame(ID_company = as.character(), timepoint = as.character(), intensity = as.numeric(), stringsAsFactors = F)
for (c in 1:length(list_of_companies)) { #treat one company at a time
  dtfm_temporal <- dtfm[dtfm$ID_company==list_of_companies[c],]
  if (nrow(dtfm_temporal)>3) { #companies need at least 4 data points to go through this process (2yr changes plus 0,0,0 before and 0,0,0 after)
    for (t in 2:(nrow(dtfm_temporal)-2)) { # this range avoids first, last and second last rows, so that it doesn't throw error
      vcre_pre <- dtfm_temporal[t-1, 3]
      vdel_pre <- dtfm_temporal[t-1, 4]
      vcap_pre <- dtfm_temporal[t-1, 5]
      vcre <- dtfm_temporal[t, 3]
      vdel <- dtfm_temporal[t, 4]
      vcap <- dtfm_temporal[t, 5]
      vcre_post <- dtfm_temporal[t+1, 3]
      vdel_post <- dtfm_temporal[t+1, 4]
      vcap_post <- dtfm_temporal[t+1, 5]
      vcre_post2 <- dtfm_temporal[t+2, 3]
      vdel_post2 <- dtfm_temporal[t+2, 4]
      vcap_post2 <- dtfm_temporal[t+2, 5]
      #if ((vcre_pre+vdel_pre+vcap_pre==0) & (vcre_post2+vdel_post2+vcap_post2==0)) { #if it's bounded by 0's, then...
      #case 2:
      if (vcre>=1 & vdel>=1 & vcap>=0 & vcre_post>=0 & vdel_post>=0 & vcap_post>=1) { dtfm_result2[nrow(dtfm_result2)+1, ] <- list(dtfm_temporal[t,1], dtfm_temporal[t+1,2], 2) }
      #case 3:
      if (vcre>=1 & vdel>=0 & vcap>=1 & vcre_post>=0 & vdel_post>=1 & vcap_post>=0) { dtfm_result2[nrow(dtfm_result2)+1, ] <- list(dtfm_temporal[t,1], dtfm_temporal[t+1,2], 2) }
      #case 4:
      if (vcre>=0 & vdel>=1 & vcap>=1 & vcre_post>=1 & vdel_post>=0 & vcap_post>=0) { dtfm_result2[nrow(dtfm_result2)+1, ] <- list(dtfm_temporal[t,1], dtfm_temporal[t+1,2], 2) }
      #case 5:
      if (vcre>=1 & vdel>=0 & vcap>=0 & vcre_post>=0 & vdel_post>=1 & vcap_post>=1) { dtfm_result2[nrow(dtfm_result2)+1, ] <- list(dtfm_temporal[t,1], dtfm_temporal[t+1,2], 2) }
      #case 6:
      if (vcre>=0 & vdel>=1 & vcap>=0 & vcre_post>=1 & vdel_post>=0 & vcap_post>=1) { dtfm_result2[nrow(dtfm_result2)+1, ] <- list(dtfm_temporal[t,1], dtfm_temporal[t+1,2], 2) }
      #case 7:
      if (vcre>=0 & vdel>=0 & vcap>=1 & vcre_post>=1 & vdel_post>=1 & vcap_post>=0) { dtfm_result2[nrow(dtfm_result2)+1, ] <- list(dtfm_temporal[t,1], dtfm_temporal[t+1,2], 2) }
      #}
    }
  }
}
#-end

#-3yr changes (speed=1):
dtfm_result3 <- data.frame(ID_company = as.character(), timepoint = as.character(), intensity = as.numeric(), stringsAsFactors = F)
for (c in 1:length(list_of_companies)) { #treat one company at a time
  dtfm_temporal <- dtfm[dtfm$ID_company==list_of_companies[c],]
  if (nrow(dtfm_temporal)>4) { #companies need at least 5 data points to go through this process (3yr changes plus 0,0,0 before and 0,0,0 after)
    for (t in 2:(nrow(dtfm_temporal)-3)) { # this range avoids first, last and second last rows, so that it doesn't throw error
      vcre_pre <- dtfm_temporal[t-1, 3]
      vdel_pre <- dtfm_temporal[t-1, 4]
      vcap_pre <- dtfm_temporal[t-1, 5]
      vcre <- dtfm_temporal[t, 3]
      vdel <- dtfm_temporal[t, 4]
      vcap <- dtfm_temporal[t, 5]
      vcre_post <- dtfm_temporal[t+1, 3]
      vdel_post <- dtfm_temporal[t+1, 4]
      vcap_post <- dtfm_temporal[t+1, 5]
      vcre_post2 <- dtfm_temporal[t+2, 3]
      vdel_post2 <- dtfm_temporal[t+2, 4]
      vcap_post2 <- dtfm_temporal[t+2, 5]
      vcre_post3 <- dtfm_temporal[t+3, 3]
      vdel_post3 <- dtfm_temporal[t+3, 4]
      vcap_post3 <- dtfm_temporal[t+3, 5]
      #if ((vcre_pre+vdel_pre+vcap_pre==0) & (vcre_post3+vdel_post3+vcap_post3==0)) { #if it's bounded by 0's, then...
      #case 8:
      if (vcre>=1 & vdel>=0 & vcap>=0 & vcre_post>=0 & vdel_post>=1 & vcap_post>=0 & vcre_post2>=0 & vdel_post2>=0 & vcap_post2>=1) { dtfm_result3[nrow(dtfm_result3)+1, ] <- list(dtfm_temporal[t,1], dtfm_temporal[t+2,2], 1) }
      #case 9:
      if (vcre>=0 & vdel>=1 & vcap>=0 & vcre_post>=1 & vdel_post>=0 & vcap_post>=0 & vcre_post2>=0 & vdel_post2>=0 & vcap_post2>=1) { dtfm_result3[nrow(dtfm_result3)+1, ] <- list(dtfm_temporal[t,1], dtfm_temporal[t+2,2], 1) }
      #case 10:
      if (vcre>=0 & vdel>=0 & vcap>=1 & vcre_post>=1 & vdel_post>=0 & vcap_post>=0 & vcre_post2>=0 & vdel_post2>=1 & vcap_post2>=0) { dtfm_result3[nrow(dtfm_result3)+1, ] <- list(dtfm_temporal[t,1], dtfm_temporal[t+2,2], 1) }
      #case 11:
      if (vcre>=0 & vdel>=1 & vcap>=0 & vcre_post>=0 & vdel_post>=0 & vcap_post>=1 & vcre_post2>=1 & vdel_post2>=0 & vcap_post2>=0) { dtfm_result3[nrow(dtfm_result3)+1, ] <- list(dtfm_temporal[t,1], dtfm_temporal[t+2,2], 1) }
      #case 12:
      if (vcre>=1 & vdel>=0 & vcap>=0 & vcre_post>=0 & vdel_post>=0 & vcap_post>=1 & vcre_post2>=0 & vdel_post2>=1 & vcap_post2>=0) { dtfm_result3[nrow(dtfm_result3)+1, ] <- list(dtfm_temporal[t,1], dtfm_temporal[t+2,2], 1) }
      #case 13:
      if (vcre>=0 & vdel>=0 & vcap>=1 & vcre_post>=0 & vdel_post>=1 & vcap_post>=0 & vcre_post2>=1 & vdel_post2>=0 & vcap_post2>=0) { dtfm_result3[nrow(dtfm_result3)+1, ] <- list(dtfm_temporal[t,1], dtfm_temporal[t+2,2], 1) }
      #}
    }
  }
}
#-end

result_bmc <- rbind(unique(dtfm_result), unique(dtfm_result2), unique(dtfm_result3)) #Important: unique is used because the >=0 clause makes the same event to fall multiple times in the conditionals, leading to duplicates
result_bmc <- result_bmc[order(result_bmc$ID_company, result_bmc$timepoint), ]
result_bmc_2 <- aggregate(intensity ~ ., max, data = result_bmc) # collapse timepoints with multiple change events, by leaving those with the highest intensity

result_bmc_summary <- result_bmc_2 %>% 
  group_by(ID_company) %>% 
  summarise(n=n())

results_bmc_summary_per_change <- result_bmc_2 %>% 
  group_by(intensity) %>% 
  summarise(n=n())
#===End STEP 5 ==================================================================================================================================================================================




#===STEP 6 (group firms by BMC pattern - in terms of intensity and frequency of changes ==========================================================================================================
result_bmc_3 <- result_bmc_2 %>% 
  group_by(ID_company) %>% 
  summarise(intensity = as.numeric(mean(intensity)), frequency = as.numeric(n()))

result_bmc_3 <- data.frame(result_bmc_3)
result_bmc_3$intensity_group <- cut(result_bmc_3$intensity, breaks=c(0,2,3), labels=c("LoInt", "HiInt")) #categorise cutting into specific values, so that firms with same intensity are not split into different grous
result_bmc_3$frequency_group <- cut(result_bmc_3$frequency, breaks=c(0,1,10), labels=c("LoFre", "HiFre")) #categorise cutting into specific values, so that firms with same frequency are not split into different grous
result_bmc_3$BMC_pattern_group <- factor(paste(result_bmc_3$intensity_group, result_bmc_3$frequency_group, sep = "-"))
#===End STEP 6 ====================================================================================================================================================================================




#===STEP 7 (group firms by firm characteristics)===================================================================================================================================================
#7A. Size (log(total assets)) - results in a dataframe with companies and sizes
result_size <- dcast(data, ID_company + timepoint ~ parameter, value.var = "value")
result_size <- result_size[,c("ID_company", "timepoint", "Tota")]
result_size <- result_size[result_size$ID_company %in% unique(result_bmc_3$ID_company),] # filter companies to align with bmc data
result_size$Tota <- log(result_size$Tota)
result_size_2 <- result_size %>% group_by(ID_company) %>% summarise(Tota_mean = mean(Tota, na.rm = T))
result_size_2$Tota_group <- ntile(result_size_2$Tota_mean, 5) #categorise by cutting into 4 groups of equal size (doesn't have the problem of intensity and frequency as it's truly continuous)
result_size_2$Tota_group <- factor(result_size_2$Tota_group, labels = c("Sm", "SmMi", "Mi", "MiLa", "La"))

#7B. Industry - results in a dataframe with companies and industries
result_industry <- aggregate(gics_code ~ ., max, data = data[,c(1,3)]) #collapse all multiple records of companies into one (and corresponding GICS code)
result_industry <- result_industry[result_industry$ID_company %in% unique(result_bmc_3$ID_company),] # filter companies to align with bmc data
result_industry$gics_code <- substr(result_industry$gics_code, 1, 6) #level of industry analysis: Industry (6-digit GICS code)
result_industry$gics_code <- factor(result_industry$gics_code)

#7C. Age (year of incorporation) - results in a dataframe with companies and ages
result_age <- aggregate(incorp_date ~ ., max, data = data[,c(1,4)]) #collapse all multiple records of companies into one (and corresponding incorp_date)
result_age <- result_age[result_age$ID_company %in% unique(result_bmc_3$ID_company),] # filter companies to align with bmc data
result_age$incorp_date <- as.numeric(substr(result_age$incorp_date, 1, 4)) #get only the year
result_age$age_group <- cut(result_age$incorp_date, breaks=c(0,1969,1984,1991,1997,2007), labels = c("Mature", "Mid-Mature", "Mid", "Young-Mid", "Young")) #categorise cutting into specific values, so that firms with same incorp year are not split into different grous
#=== End STEP 7 ====================================================================================================================================================================================




#===STEP 8 (Run Chi_Square test, print results and contigency tables) ==============================================================================================================================
# merge all firm BMC activity patterns and firm characteristics dataframes into a single one:
result_total <- merge(result_bmc_3[,c(1,6)], result_size_2[,c(1,3)], by = "ID_company")
result_total <- merge(result_total, result_industry, by = "ID_company")
result_total <- merge(result_total, result_age[,c(1,3)], by = "ID_company")

# calculate Chi Square tests for BMC activity patterns ~ size!
chisq_size <- chisq.test(result_total[,2], result_total[,3])
chisq_size

#build desc. and frequency table for BMC activity patterns ~ size!
size_table <- result_total %>% 
  group_by(Tota_group) %>% 
  summarise(n=n())

size_table2 <- result_total %>%
  group_by(Tota_group, BMC_pattern_group) %>% 
  summarise(n=n()) %>% 
  tidyr::spread(BMC_pattern_group, n)

# calculate Chi Square tests for BMC activity patterns ~ age!
chisq_age <- chisq.test(result_total[,2], result_total[,5])
chisq_age

#build desc. and frequency table for BMC activity patterns ~ age!
age_table <- result_total %>% 
  group_by(age_group) %>% 
  summarise(n=n())

age_table2 <- result_total %>% 
  group_by(age_group, BMC_pattern_group) %>% 
  summarise(n=n()) %>% 
  tidyr::spread(BMC_pattern_group, n)

# calculate Chi Square tests for BMC activity patterns ~ industry!
chisq_industry <- chisq.test(result_total[,2], result_total[,4])
chisq_industry

#build desc. and frequency table for BMC activity patterns ~ industry!
industry_table <- result_total %>% 
  group_by(gics_code) %>% 
  summarise(n=n())

industry_table2 <- result_total %>%
  group_by(gics_code, BMC_pattern_group) %>% 
  summarise(n=n()) %>% 
  tidyr::spread(BMC_pattern_group, n)
#===End STEP 8 ====================================================================================================================================================================================




#===Step 9 (visualization) ========================================================================================================================================================================
#set dataframe with value dimension data for plotting:
df_waffle_yrs <- dtfm %>% #dtfm is original dataframe with companies, timepoints and value dim changes
  group_by(timepoint) %>% 
  summarise(vcre=sum(vcre_change, na.rm = T), vdel=sum(vdel_change, na.rm = T), vcap=sum(vcap_change, na.rm = T))

df_waffle_yrs <- melt(df_waffle_yrs, id.vars="timepoint") # transform to long format
df_waffle_yrs <- data.frame(df_waffle_yrs)
df_waffle_yrs <- df_waffle_yrs[df_waffle_yrs$timepoint!="2006" & df_waffle_yrs$timepoint!="2017",] #remove first and last year to keep consistency with the BMC plot (as BMCs are  calculated from 2007 to 2016)

#set df with bmc data for plotting
df_waffle_yrs2 <- result_bmc_2 %>% 
  group_by(timepoint, intensity) %>% 
  summarise(n=n())
df_waffle_yrs2 <- data.frame(df_waffle_yrs2)

#Good reference (https://awesomeopensource.com/project/hrbrmstr/waffle)

#Plot 1 - Calculated value dim. change events across timepoints:
p1 <- ggplot(df_waffle_yrs, aes(fill = variable, values = value)) +
  geom_waffle(color = "white", size = .3, n_rows = 25, flip = T) +
  facet_wrap(~timepoint, nrow = 1, strip.position = "bottom") +
  scale_x_discrete() + 
  scale_y_continuous(labels = function(x) x * 25, # make this multiplyer the same as n_rows
                     expand = c(0,0)) +
  scale_fill_manual(name = "Value dimension", labels = c("Value creation", "Value delivery", "Value capture"), values = c("#d95f02", "#7570b3", "#1b9e77")) +
  coord_equal() +
  labs(title = "A. Change events by individual value dimension", x = "Year", y = "Count") +
  theme_minimal(base_family = "sans", base_size = 22) +
  theme(panel.grid = element_blank(), axis.ticks.y = element_line(), plot.title = element_text(size=27, face="bold", hjust = 0.5), plot.margin=unit(c(0,0,-1,0), "cm")) +
  guides(fill = guide_legend(reverse = T))

#Plot 2 - Calculated BMC events by intesity, across timepoints (BMC requency variable not included here because firms have ~1-2 changes over the 10 years...requires a different plot to visualise)
p2 <- ggplot(df_waffle_yrs2, aes(fill = intensity, values = n)) +
  geom_waffle(color = "white", size = .3, n_rows = 15, flip = T) +
  facet_wrap(~timepoint, nrow = 1, strip.position = "bottom") +
  scale_x_discrete() + 
  scale_y_continuous(labels = function(x) x * 15, # make this multiplyer the same as n_rows
                     expand = c(0,0)) +
  scale_fill_manual(name = "BMC intensity", labels = c("Low: 3-yr change", "Mid: 2-yr change", "High: 1-yr change"), values = c("#deebf7", "#9ecae1", "#3182bd")) +
  coord_equal() +
  labs(title = "B. Business model change events by intensity", x = "Year", y = "Count") +
  theme_minimal(base_family = "sans", base_size = 22) +
  theme(panel.grid = element_blank(), axis.ticks.y = element_line(), plot.title = element_text(size=27, face="bold", hjust = 0.5),  plot.margin=unit(c(-1,0,0,0), "cm")) +
  guides(fill = guide_legend(reverse = T))

ggarrange(p1, p2, nrow = 2, ncol = 1)
#===End STEP 9 ====================================================================================================================================================================================


