## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----eval=FALSE---------------------------------------------------------------
# library(starling)
# 
# # Load example datasets
# data(dx_data)      # Diagnosis/case data
# data(vax_data)     # Vaccination records
# 
# # Clean the case data
# cases <- clean_the_nest(
#   data = dx_data,
#   data_type = "cases",
#   id_var = "identity",
#   diagnosis = "disease_name",
#   lettername1 = "first_name",
#   lettername2 = "surname",
#   dob = "date_of_birth",
#   medicare = "medicare_no",
#   gender = "gender",
#   postcode = "postcode",
#   onset_date = "diagnosis_date"
# )
# 
# # Clean the vaccination data (convert to wide format)
# vaccines <- clean_the_nest(
#   data = vax_data,
#   data_type = "vaccination",
#   lie_nest_flat = TRUE,  # Convert from long to wide format
#   id_var = "patient_id",
#   lettername1 = "firstname",
#   lettername2 = "last_name",
#   dob = "birth_date",
#   medicare = "medicare_number",
#   gender = "gender",
#   postcode = "postcode",
#   vax_type = "vaccine_delivered",
#   vax_date = "service_date"
# )

## ----eval=FALSE---------------------------------------------------------------
# # Link cases to vaccination history
# linked_data <- murmuration(
#   df1 = cases,
#   df2 = vaccines,
#   linkage_type = "v2c",
#   blocking_var = "gender",
#   compare_vars = c("lettername1", "lettername2", "dob"),
#   threshold_value = 12,
#   days_allowed_before_event = 7,
#   clean_eggs = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# # Create analysis-ready dataset with labels and categories
# final_data <- preening(linked_data)
# 
# # Now ready for analysis!
# library(gtsummary)
# final_data %>%
#   select(age5cat, gender, vaccination_status, admission_outcome) %>%
#   tbl_summary(by = vaccination_status)

## ----eval=FALSE---------------------------------------------------------------
# data(dx_data)
# data(hosp_data)
# data(vax_data)
# 
# # Clean all datasets
# cases <- clean_the_nest(
#   data = dx_data,
#   data_type = "cases",
#   id_var = "identity",
#   diagnosis = "disease_name",
#   lettername1 = "first_name",
#   lettername2 = "surname",
#   dob = "date_of_birth",
#   medicare = "medicare_no",
#   gender = "gender",
#   postcode = "postcode",
#   onset_date = "diagnosis_date"
# )
# 
# hospitals <- clean_the_nest(
#   data = hosp_data,
#   data_type = "hospital",
#   id_var = "patient_id",
#   lettername1 = "firstname",
#   lettername2 = "last_name",
#   dob = "birth_date",
#   medicare = "medicare_number",
#   gender = "sex",
#   postcode = "zip_codes",
#   icd_code = "icd_codes",
#   admission_date = "date_of_admission",
#   discharge_date = "date_of_discharge"
# )
# 
# vaccines <- clean_the_nest(
#   data = vax_data,
#   data_type = "vaccination",
#   lie_nest_flat = TRUE,
#   id_var = "patient_id",
#   lettername1 = "firstname",
#   lettername2 = "last_name",
#   dob = "birth_date",
#   medicare = "medicare_number",
#   gender = "gender",
#   postcode = "postcode",
#   vax_type = "vaccine_delivered",
#   vax_date = "service_date"
# )
# 
# # Link cases to vaccinations
# cases_vax <- murmuration(
#   df1 = cases,
#   df2 = vaccines,
#   linkage_type = "v2c",
#   blocking_var = "gender",
#   compare_vars = c("lettername1", "lettername2", "dob", "medicare10"),
#   threshold_value = 12,
#   clean_eggs = TRUE
# )
# 
# # Link the combined dataset to hospitalizations
# complete_data <- murmuration(
#   df1 = cases_vax,
#   df2 = hospitals,
#   linkage_type = "v2h",
#   blocking_var = "gender",
#   compare_vars = c("lettername1", "lettername2", "dob", "medicare10"),
#   days_allowed_before_event = 7,
#   days_allowed_after_event = 30,
#   one_row_per_person = TRUE,
#   clean_eggs = TRUE
# )
# 
# # Prepare for analysis
# analysis_data <- preening(complete_data)

## ----eval=FALSE---------------------------------------------------------------
# data(manifest_data)
# data(vax_data)
# 
# # Clean manifest data
# manifest <- clean_the_nest(
#   data = manifest_data,
#   data_type = "cases",
#   id_var = "passenger_id",
#   lettername1 = "first_name",
#   lettername2 = "surname",
#   dob = "date_of_birth",
#   gender = "gender"
# )
# 
# # Clean vaccination data
# vaccines <- clean_the_nest(
#   data = vax_data,
#   data_type = "vaccination",
#   lie_nest_flat = TRUE,
#   id_var = "patient_id",
#   lettername1 = "firstname",
#   lettername2 = "last_name",
#   dob = "birth_date",
#   gender = "gender",
#   vax_type = "vaccine_delivered",
#   vax_date = "service_date"
# )
# 
# # Link to determine vaccination status at time of flight
# flight_vax <- murmuration(
#   df1 = manifest,
#   df2 = vaccines,
#   linkage_type = "v2e",
#   event_date = as.Date("2024-03-15"),  # Flight date
#   blocking_var = "gender",
#   compare_vars = c("lettername1", "lettername2", "dob"),
#   days_allowed_before_event = 14,  # Valid if vaccinated ≥14 days before flight
#   clean_eggs = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# data(linelist_data)
# data(vax_data)
# 
# # Clean outbreak linelist
# outbreak <- clean_the_nest(
#   data = linelist_data,
#   data_type = "cases",
#   id_var = "case_id",
#   lettername1 = "first_name",
#   lettername2 = "surname",
#   dob = "date_of_birth",
#   medicare = "medicare_no",
#   gender = "gender",
#   postcode = "postcode",
#   onset_date = "onset_date"
# )
# 
# # Clean vaccination data
# vaccines <- clean_the_nest(
#   data = vax_data,
#   data_type = "vaccination",
#   lie_nest_flat = TRUE,
#   id_var = "patient_id",
#   lettername1 = "firstname",
#   lettername2 = "last_name",
#   dob = "birth_date",
#   medicare = "medicare_number",
#   gender = "gender",
#   postcode = "postcode",
#   vax_type = "vaccine_delivered",
#   vax_date = "service_date"
# )
# 
# # Link to determine vaccination status at time of exposure
# outbreak_vax <- murmuration(
#   df1 = outbreak,
#   df2 = vaccines,
#   linkage_type = "v2e",
#   event_date = as.Date("2024-06-01"),  # Festival date
#   blocking_var = "postcode",
#   compare_vars = c("lettername1", "lettername2", "dob", "medicare10"),
#   days_allowed_before_event = 7,
#   clean_eggs = TRUE
# )

## ----eval=FALSE---------------------------------------------------------------
# # Check date format
# class(data$dob)  # Should return "Date"
# 
# # Convert if needed
# data$dob <- as.Date(data$dob, format = "%Y-%m-%d")
# # or using lubridate
# data$dob <- lubridate::ymd(data$dob)

## ----eval=FALSE---------------------------------------------------------------
# citation("starling")

