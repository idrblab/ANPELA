% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Ranking.R
\name{Ranking}
\alias{Ranking}
\title{Ranking}
\usage{
Ranking(data, name = "result", savepath = "./ANPELA_res")
}
\arguments{
\item{data}{Character, the R object resulting from the function "Assess" "CSIassess" or "PTIassess", or obtained by loading from the resulting RData file of these funcitons.}

\item{name}{Character, the filename of the overall ranking data and figure file.}

\item{savepath}{Character, the absolute path of the folder which will store the overall ranking data and figure file.}
}
\value{
A CSV file named "_Ranking_Table.csv", recording the overall ranking and the values of criteria. \if{html}{\out{<br>}}A PDF file named "_Ranking_Figure.pdf", helping users better understand the differences among various data processing workflows, where the different colors represent different performance assessment levels: green indicates "good," and red indicates "poor".
}
\description{
Ranking ranks all workflows assessed by the function "CSIassess" and "PTIassess" based on collective consideration of values and grades (classified by well-defined cutoffs) under each criterion.
}
\examples{
\donttest{
}
}
