\name{joe}
\alias{joe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CDO Pricing with the Joe copula
}
\description{
'joe' compute 5 tranches spreads under the Joe copula
}
\usage{
joe(rhoInput, MInput, dateInput)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rhoInput}{
a numeric giving the Kendall's tau for dependence specification
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
joe(0.3, 100, c("2007-10-23"))
}
