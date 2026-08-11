\name{pohang}
\alias{pohang}
\docType{data}
\title{Pohang Annual r-Largest Daily Precipitation Data}
\description{
Annual r-largest daily precipitation data for Pohang (station 138), Korea.
}
\usage{data(pohang)}
\format{
A numeric matrix with 74 rows and 10 columns. Rows correspond to years from
1949 through 2022. Columns contain the first through tenth largest daily
precipitation observations within each year, ordered from largest to smaller values.
Measurements are in millimeters.
}
\source{
Korea Meteorological Administration, Open MET Data Portal (\url{https://data.kma.go.kr}).
}
\details{
The dataset is included as an example of the block-wise input format expected by
\code{ms.rlos()}. The package uses the annual block structure and assumes columns
are arranged in decreasing order within each row.
}
\keyword{datasets}
