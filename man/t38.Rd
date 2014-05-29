\name{t38}
\alias{t38}
\title{
CDO Pricing with the C_t38 copula
}
\description{
't38' compute 5 tranches spreads under the C_t38 copula
}
\usage{
t38(theta1Input, theta2Input, dF, MInput, dateInput)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta1Input}{
a numeric giving the Pearson's correlation for dependence specification
}
  \item{theta2Input}{
a numeric giving the Pearson's correlation for dependence specification
}
  \item{dF}{
a numeric giving the degree of freedom
}
  \item{MInput}{
a numeric giving the Monte Carlo simulation runs
}
  \item{dateInput}{
a character giving the pricing date, e.g. c("2007-10-23")
}
}
\details{
Please make sure that the data sets of "defIntensity.csv" and "payday.csv" have been correctly installed in such paths: "C:/defIntensity.csv", "C:/payday.csv". The both data sets can be downloaded from "https://github.com/YafeiXu/xyfQuantlet".
}
\value{
A vector with 5 numerics will be returned, from left to right: equity, junior mezzanine, senior mezzanine, junior senior, senior.
}
\references{
The master thesis, CDO, HAME Copulas and an R Package 'CDO', can be downloaded from
https://sites.google.com/site/cdowithr/.
}
\author{
Yafei Xu <yafei.xu@hu-berlin.de>
}

\examples{
t38(0.3, 0.3, 13, 100, c("2007-10-23"))
}
