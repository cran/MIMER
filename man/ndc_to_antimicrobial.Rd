% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_antibiotics.R
\name{ndc_to_antimicrobial}
\alias{ndc_to_antimicrobial}
\title{Convert 'ndc' code to corresponding Antibiotic code.}
\usage{
ndc_to_antimicrobial(ndc,
 class_names,
 include_missing_NDCs = TRUE)
}
\arguments{
\item{ndc}{A vector containing ndc codes. Will be coerced to character.}

\item{class_names}{A vector containing antibacterial class names - eg: c("antimicrobial", "antibacterial").}

\item{include_missing_NDCs}{includes a hardcoded database of NDCs that are present in MIMIC-IV but not in NDC database.}
}
\value{
Vector of antimicrobials in antibiotic class from AMR package.
}
\description{
Function to convert 'ndc' code to corresponding Antibiotic code.
}
\details{
ndc_to_antimicrobial
}
