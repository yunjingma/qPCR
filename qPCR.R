library(dplyr)
library(tidyr)
# Read the data into a data frame and specify "#VALUE!" as NA
# Replace the file path "path_to_your_file.xlsx" with the actual file path
df <- readxl::read_excel(path = "~/Desktop/qPCR/qPCR/raw_data.xlsx", na = "#VALUE!")

# Clean the data and convert Ct values to numeric
df$ct <- gsub(",", ".", df$ct)
df$ct <- as.numeric(df$ct)

# Calculate the average Ct values for each gene in each sample
df_avg <- df %>%
  group_by(Sample, Gene) %>%
  summarise(avg_ct = mean(ct, na.rm = TRUE))

# Pivot the dataframe to have genes as columns
df_pivot <- pivot_wider(df_avg, names_from = Gene, values_from = avg_ct)

# Calculate ΔCt for ERF114 and JAZ10
df_pivot$delta_Ct_ERF114 <- df_pivot$ERF114 - df_pivot$UBQ
df_pivot$delta_Ct_JAZ10 <- df_pivot$JAZ10 - df_pivot$UBQ

write.csv(df_pivot, "df_pivot.csv", row.names = FALSE)

# Replace sample names with "Col-0 Mock" where applicable
for (i in 22:24) {
  df_pivot$Sample <- gsub(paste0("cDNA", i), "the1_EGCG", df_pivot$Sample)
}

# Now 'df_pivot' should have "Col-0 Mock" as sample names instead of "cDNA1", "cDNA2", etc.

# Assuming the control group is Col-0 mock
control_group <- subset(df_pivot, grepl("Col-0_Mock", Sample))

# Calculate ΔΔCt for ERF114 and JAZ10
df_pivot <- merge(df_pivot, control_group, by = "Sample", suffixes = c("", "_control"))
df_pivot$delta_delta_Ct_ERF114 <- df_pivot$delta_Ct_ERF114 - df_pivot$delta_Ct_ERF114_control
df_pivot$delta_delta_Ct_JAZ10 <- df_pivot$delta_Ct_JAZ10 - df_pivot$delta_Ct_JAZ10_control

# Calculate fold change (FC) for ERF114 and JAZ10
df_pivot$FC_ERF114 <- 2^(-df_pivot$delta_delta_Ct_ERF114)
df_pivot$FC_JAZ10 <- 2^(-df_pivot$delta_delta_Ct_JAZ10)

# Perform statistical analysis if needed
# For example, t-test comparing treatment groups to the control group

# Assuming treatment group is EGCG treatment for both ERF114 and JAZ10
treatment_group <- subset(df_pivot, grepl("EGCG", Sample))

# Perform t-test for ERF114
ttest_ERF114 <- t.test(treatment_group$FC_ERF114, control_group$FC_ERF114)

# Perform t-test for JAZ10
ttest_JAZ10 <- t.test(treatment_group$FC_JAZ10, control_group$FC_JAZ10)

# Print the results
print("ERF114 t-test:")
print(ttest_ERF114)
print("JAZ10 t-test:")
print(ttest_JAZ10)

# Replace commas with decimal points in the 'ct' column
data$ct <- gsub(",", ".", data$ct)

# Convert 'ct' column to numeric
data$ct <- as.numeric(data$ct)

# Custom function to calculate mean while ignoring NA values
mean_without_na <- function(x) {
  if (all(is.na(x))) {
    return(NA)  # Return NA if all values are NA
  } else {
    return(mean(x, na.rm = TRUE))  # Calculate mean without NA values
  }
}
# Create a new column to indicate the index within each group
avg_ct$Index <- rep(1:(nrow(data)/3), each = 3)

# Define the groups manually, considering the missing sample
data$Group <- ifelse(data$Index %in% c(5, 11, 14), data$Index - 1, data$Index)

# Remove unnecessary columns
data <- data[, !(names(data) %in% c("Index"))]

# Calculate average ct values while considering only non-NA values
avg_ct <- aggregate(ct ~ Sample + Gene, data = data, FUN = mean_without_na)

cDNA24_avg <- data.frame(Sample = "cDNA24", Gene = "UBQ", ct = 22.71200)

# Bind the new data frame to avg_ct
avg_ct <- rbind(avg_ct, cDNA24_avg)


# Step 3: Normalize Ct values (optional if you have a reference gene)
# Assuming UBQ is the reference gene
ref_gene <- "UBQ"
avg_ct$Norm_ct <- ifelse(avg_ct$Gene == ref_gene, avg_ct$ct, avg_ct$ct - avg_ct$ct[avg_ct$Gene == ref_gene])

# Step 4: Calculate relative expression levels
# Assuming ΔΔCt method with the first 3 samples as reference
reference_samples <- avg_ct$Norm_ct[1:3]

avg_ct$Rel_Expression <- ifelse(avg_ct$Sample %in% avg_ct$Sample[1:3], 1, 2^-(avg_ct$Norm_ct - reference_samples))





# Now you can proceed with your analysis using the 'Group' column






# Calculate the standard deviation of relative expression for each group
stdev <- aggregate(Relative_Expression ~ Group, data = data, FUN = sd)

# Print the result
print(stdev)

# Step 5: Perform statistical analysis if needed
# For example, you can perform t-tests between groups
t_test_result <- t.test(Rel_Expression ~ Sample, data = avg_ct)

# View the results
print(avg_ct)
print(t_test_result)
