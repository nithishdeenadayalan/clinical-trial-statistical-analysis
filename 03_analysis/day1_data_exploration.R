# -----------------------------------------------
# Clinical Trial Statistical Analysis
# Day 1 — Data Exploration and Documentation
# Nithish Deenadayalan
# -----------------------------------------------

library(survival)
library(dplyr)

# Load data
data(lung, package = "survival")

# -----------------------------------------------
# STEP 1 — Basic exploration
# -----------------------------------------------

# Dimensions
cat("Rows:", nrow(lung), "\n")
cat("Columns:", ncol(lung), "\n")

# Variable types
str(lung)

# Summary statistics
summary(lung)

# -----------------------------------------------
# STEP 2 — Missing value analysis
# -----------------------------------------------

cat("\n--- Missing values per variable ---\n")
missing_summary <- data.frame(
  variable    = names(lung),
  missing_n   = colSums(is.na(lung)),
  missing_pct = round(colSums(is.na(lung)) / nrow(lung) * 100, 1)
)
print(missing_summary)

# -----------------------------------------------
# STEP 3 — Key distributions
# -----------------------------------------------

cat("\n--- Age distribution ---\n")
cat("Mean age:", round(mean(lung$age, na.rm = TRUE), 1), "\n")
cat("Range:", min(lung$age, na.rm = TRUE), "to", max(lung$age, na.rm = TRUE), "\n")

cat("\n--- Sex breakdown ---\n")
print(table(lung$sex))
cat("1 = Male, 2 = Female\n")

cat("\n--- ECOG performance score breakdown ---\n")
print(table(lung$ph.ecog, useNA = "always"))
cat("0 = Fully active, 1 = Restricted, 2 = Ambulatory, 3 = Limited self-care\n")

cat("\n--- Survival status ---\n")
print(table(lung$status))
cat("1 = Censored (alive at study end), 2 = Dead\n")

cat("\n--- Median survival time ---\n")
cat("Median survival:", median(lung$time, na.rm = TRUE), "days\n")
cat("Range:", min(lung$time), "to", max(lung$time), "days\n")




# -----------------------------------------------
# STEP 4 — Build the Data Dictionary
# -----------------------------------------------

data_dictionary <- data.frame(
  variable = c("inst","time","status","age","sex",
               "ph.ecog","ph.karno","pat.karno","meal.cal","wt.loss"),
  type = c("Numeric","Numeric","Numeric","Numeric","Numeric",
           "Numeric","Numeric","Numeric","Numeric","Numeric"),
  description = c(
    "Institution code — hospital where patient was treated",
    "Survival time in days from study entry to death or last follow-up",
    "Patient status at end of study — 1 = censored (alive), 2 = dead",
    "Patient age in years at study entry",
    "Patient sex — 1 = Male, 2 = Female",
    "ECOG performance score rated by physician — 0 = fully active, 1 = restricted but ambulatory, 2 = ambulatory but unable to work, 3 = limited self-care",
    "Karnofsky performance score rated by physician — 0 to 100 scale",
    "Karnofsky performance score rated by patient — 0 to 100 scale",
    "Calories consumed at meals",
    "Weight loss in pounds in last 6 months — negative values indicate weight gain"
  ),
  missing_n = colSums(is.na(lung)),
  missing_pct = paste0(round(colSums(is.na(lung)) / nrow(lung) * 100, 1), "%"),
  missing_strategy = c(
    "Exclude from analysis — not a clinical variable",
    "No missing values — primary outcome variable",
    "No missing values — primary outcome variable",
    "No missing values — complete case analysis",
    "No missing values — complete case analysis",
    "1 missing value — exclude from ECOG subgroup analysis",
    "1 missing value — exclude from Karnofsky analysis",
    "3 missing values — exclude from patient-rated analysis",
    "47 missing values — exclude variable from primary analysis",
    "14 missing values — use available cases for weight loss analysis"
  )
)

print(data_dictionary)

# Save to CSV for your documentation folder
write.csv(data_dictionary, "data_dictionary.csv", row.names = FALSE)
cat("\nData dictionary saved to data_dictionary.csv\n")