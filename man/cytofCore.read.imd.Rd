\name{cytofCore.read.imd}
\alias{cytofCore.read.imd}
\title{
Read an CyToF IMD file
}
\description{
Read CyToF IMD binary file and return intensity/pulse data or dual count data
}
\usage{
cytofCore.read.imd(file, analytes = NULL, conf = NULL, pulse_thresh = 3, start_push = 0, num_pushes = NULL)
}
\arguments{
  \item{file}{
IMD file
}
  \item{analytes}{
Vector of analyte names. Length of vector must match number of analytes in IMD file. If NULL, \code{conf} must be provided.
}
  \item{conf}{
CyToF conf file or matrix. If conf data is specified, only dual counts are returned.
}
  \item{pulse_thresh}{
Pulse threshold used in dual count computation
}
  \item{start_push}{
Starting push when reading from middle of file
}
  \item{num_pushes}{
Number of pushes to read. NULL reads to end of file.
}
}
\details{
Reads CyToF IMD binary file.

IMD file contains intensity (I) and pulse data (P) for all analytes and pushes. Format of each push is

I_analyte0 P_analyte0 I_analyte1 ...

where all values are 16-bit unsigned integers. File starts with push 0.

Dual counts are calculated as:

dual = round(intensity * slope + intercept)
if (dual < pulse_thresh)
	dual = max(dual, pulse)
}
\value{
\item{intensity }{Intensity data}
\item{pulse }{Pulse data}
\item{dual }{Dual counts}
}
\references{
Bandura, et al. Analytical Chemistry 81:6813-6822, 2009.
}
\author{
Michael Linderman (mlinderm@stanford.edu)
}
\note{
%%  ~~further notes~~
}
\examples{
}
