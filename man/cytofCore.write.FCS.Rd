\name{cytofCore.write.FCS}
\alias{cytofCore.write.FCS}
\title{
Write an FCS file
}
\description{
Write and FCS file from a data matrix. Specifically intended for use with cell data extracted from IMD files.
}
\usage{
cytofCore.write.FCS(x, filename, what = "numeric")
}

\arguments{
  \item{x}{
A data matrix
}
  \item{filename}{
A character scalar giving the output file name
}
  \item{what}{
A character scale defining the output data type. One in 'integer', 'numeric', 'double'. Note that forcing the data type to 'integer' may result in considerable loss of precision. We recommend using the default data type.
}
}
\details{
Create FCS 3.0 standard file from the input data matrix.
}
\value{
A character scale of the file name.
}
\references{
\url{http://www.isac-net.org}
}
\author{
Michael Linderman (mlinderm@stanford.edu)
}
\seealso{
\code{\link{write.FCS}}
}
\examples{
}
