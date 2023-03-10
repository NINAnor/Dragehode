% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funksjoner.R
\name{beregn.popstruktur}
\alias{beregn.popstruktur}
\title{Calculate population structure}
\usage{
beregn.popstruktur(
  utvalgt.lokalitet,
  lokalitetsdata,
  transektdata,
  rutedata,
  rutebredde = 1,
  fert = "Fert.planter",
  veg = "Veg.planter",
  sma = "Småplanter",
  tot = "Ant.DR",
  forekomst_transekt = "Dragehode",
  forekomst_avstand = "Forekomst dragehode (m)",
  quantiles = c(0.025, 0.975)
)
}
\arguments{
\item{utvalgt.lokalitet}{chosen locality}

\item{lokalitetsdata}{data for that locality}

\item{transektdata}{transect data}

\item{rutedata}{quadrate data}

\item{rutebredde}{quadrate width (defaults to 1)}

\item{fert}{Fertile plants column}

\item{veg}{Vegetative plants column}

\item{sma}{seedlins column}

\item{tot}{total plants column}

\item{forekomst_transekt}{occurrence transect data}

\item{forekomst_avstand}{occurrence distance}

\item{quantiles}{quantiles to be calsulated (dafaults to 0.025, 0.975)}
}
\value{
population structure list per location
}
\description{
Calculate population structure
}
