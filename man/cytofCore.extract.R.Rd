\name{cytofCore.extract.R}
\alias{cytofCore.extract.R}
\alias{cytofCore.extract.native}
\title{
Extraction of cells from CyToF IMD files
}
\description{
Identifies and extracts cells from CyToF IMD files using simple thresholding algorithm. Two implementations are availabe: in pure R, and with a C++ core.
}
\usage{
cytofCore.extract.R(imd, conf, pulse_thresh = 3, num_pushes = NULL, thresh = 10, sigma = 3, num_sigma = 3, min_length = 10, max_length = 75, noise.subtraction = TRUE, noise.min.length = 30, slope.filter = TRUE, freq = 77000)
cytofCore.extract.native(imd, conf, pulse_thresh = 3, num_pushes = .Machine$integer.max, thresh = 10, sigma = 3, num_sigma = 3, min_length = 10, max_length = 75, noise.subtraction = TRUE, noise.min.length = 30, slope.filter = TRUE, freq = 77000)
}
\arguments{
  \item{imd}{
IMD file
}
  \item{conf}{
CyToF conf file or matrix.
}
  \item{pulse_thresh}{
Pulse threshold used in dual count computation.
}
  \item{num_pushes}{
Number of pushes to read. Default values reads to end of file.
}
  \item{thresh}{
Threshold for start of a cell.
}
  \item{sigma}{
Standard deviation of guassian filter applied to time series prior to thresholding.
}
  \item{num_sigma}{
Number of sigma to extend gaussian filter.
}
  \item{min_length}{
Minimum length of cell in pushes.
}
  \item{max_length}{
Maximum length of cell in pushes.
}
  \item{noise.subtraction}{
Subtract noise from cell counts.
}
  \item{noise.min.length}{
Minimum length for valid noise region in pushes.
}
  \item{slope.filter}{
Set cell quality to zero if slope has more than one zero crossing.
}
  \item{freq}{
Number of pushes per second.
}
}
\details{
Extracts cells from CyToF IMD file with conservative algorithm similar but not identical to that
used by manufacturer's software. ".R" is a native R implementation, while ".native" is a C++-based
implementation that is substantially more performant.
}
\value{
	\item{cells.found }{Matrix of cells}
	\item{quality }{Quality metric, in [0,1] where 0 is lowest quality and 1 is highest quality}
}
\references{
Bandura, et al. Analytical Chemistry 81:6813-6822, 2009.
}
\author{
Michael Linderman (mlinderm@stanford.edu)
}
\examples{
}
