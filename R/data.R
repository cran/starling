#' Example Diagnosis Dataset
#'
#' A sample dataset containing diagnosis records for 10 individuals with
#' respiratory illnesses (Influenza A/B, COVID-19, RSV).
#'
#' @format A data frame with 10 rows and 10 variables:
#' \describe{
#'   \item{identity}{Patient identifier}
#'   \item{disease_name}{Name of diagnosed disease}
#'   \item{first_name}{Patient's first name}
#'   \item{surname}{Patient's surname}
#'   \item{date_of_birth}{Patient's date of birth}
#'   \item{medicare_no}{Medicare number}
#'   \item{gender}{Gender (M/F)}
#'   \item{postcode}{Postcode of residence}
#'   \item{indigenous_status}{Indigenous status}
#'   \item{diagnosis_date}{Date of diagnosis}
#' }
#' @examples
#' data(dx_data)
#' head(dx_data)
"dx_data"

#' Example Hospitalization Dataset
#'
#' A sample dataset containing hospitalization records for 10 individuals,
#' with 3 overlapping with dx_data.
#'
#' @format A data frame with 10 rows and 11 variables:
#' \describe{
#'   \item{patient_id}{Patient identifier}
#'   \item{firstname}{Patient's first name}
#'   \item{last_name}{Patient's last name}
#'   \item{birth_date}{Patient's date of birth}
#'   \item{medicare_number}{Medicare number with check digit}
#'   \item{sex}{Sex (1=Male, 2=Female)}
#'   \item{zip_codes}{Postcode}
#'   \item{cultural_heritage}{Indigenous status}
#'   \item{icd_codes}{ICD-10 diagnosis code}
#'   \item{date_of_admission}{Hospital admission date}
#'   \item{date_of_discharge}{Hospital discharge date}
#' }
#' @examples
#' data(hosp_data)
#' head(hosp_data)
"hosp_data"

#' Example Vaccination Dataset
#'
#' A sample dataset containing vaccination records in long format for multiple
#' individuals, with some having multiple vaccination entries.
#'
#' @format A data frame with 15 rows and 9 variables:
#' \describe{
#'   \item{patient_id}{Patient identifier}
#'   \item{firstname}{Patient's first name}
#'   \item{last_name}{Patient's last name}
#'   \item{birth_date}{Patient's date of birth}
#'   \item{medicare_number}{Medicare number}
#'   \item{gender}{Gender (Male/Female)}
#'   \item{postcode}{Postcode of residence}
#'   \item{vaccine_delivered}{Type of vaccine administered}
#'   \item{service_date}{Date of vaccination}
#' }
#' @examples
#' data(vax_data)
#' head(vax_data)
"vax_data"

#' Example Flight Manifest Dataset
#'
#' A sample dataset containing passenger manifest for flight QF123 on
#' 2024-03-15, used to demonstrate event-based vaccination linkage.
#'
#' @format A data frame with 8 rows and 8 variables:
#' \describe{
#'   \item{passenger_id}{Passenger identifier}
#'   \item{first_name}{Passenger's first name}
#'   \item{surname}{Passenger's surname}
#'   \item{date_of_birth}{Passenger's date of birth}
#'   \item{gender}{Gender (M/F)}
#'   \item{seat_number}{Seat assignment}
#'   \item{flight_number}{Flight number}
#'   \item{flight_date}{Date of flight}
#' }
#' @examples
#' data(manifest_data)
#' head(manifest_data)
"manifest_data"

#' Example Outbreak Linelist Dataset
#'
#' A sample dataset containing measles outbreak cases associated with a music
#' festival on 2024-06-01.
#'
#' @format A data frame with 10 rows and 10 variables:
#' \describe{
#'   \item{case_id}{Case identifier}
#'   \item{first_name}{Patient's first name}
#'   \item{surname}{Patient's surname}
#'   \item{date_of_birth}{Patient's date of birth}
#'   \item{medicare_no}{Medicare number}
#'   \item{gender}{Gender (M/F)}
#'   \item{postcode}{Postcode of residence}
#'   \item{onset_date}{Date of symptom onset}
#'   \item{notification_date}{Date case was notified}
#'   \item{case_classification}{Case classification (Confirmed/Probable)}
#' }
#' @examples
#' data(linelist_data)
#' head(linelist_data)
"linelist_data"


