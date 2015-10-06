\name{SimesBH}
\alias{SimesBH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Benjamini-Bogomolov Selection-Adjusted Procedure (BH-BH)
}
\description{
Applies the Selection-Adjusted Procedure developed in Benjamini and Bogomolov (2014), targeting the 
control of the expected average false discovery proportion over the selected families.  Given a list of 
families of p-values, calculates the Simes p-value for each family and selects families by applying the 
Benjamini-Hochberg procedure on these p-values. Following this, inside each selected family, the 
Benjamini-Hochberg procedure is performed at the level adjusted for selection of families at the 
previous step.
}
\usage{
SimesBH(ListOfFamilies,
        familySelectionAlpha=0.05, 
        insideFamilyAlpha=0.05,
        echo="OFF")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ListOfFamilies}{
  a list containing vectors of families of p-values.}
  \item{familySelectionAlpha}{
  The level at which the Benjamini-Hochberg procedure is applied for selecting the families (the nominal 
  level for control of the FDR at the level of families.)
  }
   \item{insideFamilyAlpha}{
  The nominal level for control of the expected average false discovery proportion over the selected 
families.
  }
   \item{echo}{
  "OFF"/"ON"
  whether to print out the calculated Simes p-value and the members of each family before the selection.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
doesn't return a value.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Uriya Bartal, directed by Yoav Benjamini and Marina Bogomolov
}
\note{
Version 1.0 - 29.1.2015
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
y=list(c(0, 0.06, 0.07), 
       c(0,0.055, 0.065), 
       c(0,0.075,0.08))
SimesBH (y,echo="ON")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
