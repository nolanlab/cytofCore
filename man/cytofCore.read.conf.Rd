\name{cytofCore.read.conf}
\alias{cytofCore.read.conf}
\title{
Read CyToF produced conf file.
}
\description{
Read CyToF produced conf file.
}
\usage{
cytofCore.read.conf(file)
}

\arguments{
  \item{file}{
Conf file.
}
}
\details{
Expects CyToF-produced conf file with header line and entries with Symbol,
Mass, R2, Intercept and Slope.
}
\value{
Matrix with the above described columns
}
\author{
Michael Linderman (mlinderm@stanford.edu)
}
\examples{
}
