#' Marketing Data
#'
#' A dataset containing information on the annual household income along 
#' with 13 other demographic factors of shopping mall customers in the 
#' San Francisco Bay Area.
#' 
#' @format A data frame with 8993 rows and 14 variables:
#' \describe{
#'   \item{Income}{ANNUAL INCOME OF HOUSEHOLD (PERSONAL INCOME IF SINGLE)}
#'   \item{Sex}{Sex, 1= male, 2 = female}
#'   \item{Marital}{Marital status}
#'   \item{Age}{Age (categorized)}
#'   \item{Edu}{Education:  1. Grade 8 or less
#'                          2. Grades 9 to 11
#'                          3. Graduated high school
#'                          4. 1 to 3 years of college
#'                          5. College graduate
#'                          6. Grad Study}
#'    \item{Occupation}{Occupation}
#'    \item{HowLong}{HOW LONG HAVE YOU LIVED IN THE SAN FRAN./OAKLAND/SAN JOSE AREA?}
#'    \item{Dual}{DUAL INCOMES (IF MARRIED)}
#'    \item{Persons}{PERSONS IN YOUR HOUSEHOLD}
#'    \item{Persons<18}{PERSONS IN HOUSEHOLD UNDER 18}
#'    \item{Status}{HOUSEHOLDER STATUS}
#'    \item{Type}{TYPE OF HOME}
#'    \item{Ethnic}{ETHNIC CLASSIFICATION}
#'    \item{Language}{WHAT LANGUAGE IS SPOKEN MOST OFTEN IN YOUR HOME?
#'                    1. English
#'                    2. Spanish
#'                    3. Other}
#' }
#' 
#' @usage data(marketing)
#' 
#' @source Impact Resources, Inc., Columbus, OH (1987). Dataset available from
#' \url{https://hastie.su.domains/ElemStatLearn/data.html}, formerly part of the
#'  ElemStatLearn package.
#' 
"marketing"