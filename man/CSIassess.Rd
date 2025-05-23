% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CSIassess.R
\name{CSIassess}
\alias{CSIassess}
\title{CSI (Cell Subpopulation Identification) Assess}
\usage{
CSIassess(
  name = "result",
  data,
  respath,
  clusteringM = c("FlowSOM", "PhenoGraph", "Mclust"),
  Phenograph_k = 30,
  ncluster = 8,
  Ca_metric = "AUC",
  Cb_metric = "Silhouette coefficient (SC)",
  Cc_metric = "relative weighted consistency (CWrel)",
  ntop = NULL,
  DEP = NULL,
  marker_path = NULL,
  known_celltype_path = NULL,
  save_processed_res = "one_folder",
  savepath = "./ANPELA_res",
  cores = floor(parallel::detectCores()/2),
  ...
)
}
\arguments{
\item{name}{Character, the filename of the RData file in the "assess_res" folder which will store the assessment results.}

\item{data}{Character, the R object resulting from the function "Process", "FCprocess" or "MCprocess", or obtained by loading from the resulting RData file of these funcitons when the \code{save_processed_res} parameter in these functions is set to "one_RData".}

\item{respath}{Character, the absolute path of the folder storing the resulting "info_saved.RData" file and the "process_res" folder of the function "Process", "FCprocess" or "MCprocess" when the \code{save_processed_res} parameter in these functions is set to "one_folder".}

\item{clusteringM}{Character, the method of clustering the processed data prior to performance assessment, including "FlowSOM" and "PhenoGraph".
\if{html}{\out{<br>}}\strong{FlowSOM}: a widely used cluster clustering algorithm designed for high-dimensional cytometry data (the number of clusters needs to be specified). ANPELA uses the function “SOM” in R package "FlowSOM" to implement the algorithm.
\if{html}{\out{<br>}}\strong{PhenoGraph}: a well-designed clustering algorithm developed to define phenotypes in high-dimensional single-cell data (the number of clusters need not be specified). ANPELA uses the function “Rphenograph” in R package "cytofkit".}

\item{Phenograph_k}{Character, the number of nearest neighbours used in PhenoGraph clustering method.
\if{html}{\out{<br>}}Only needed when the argument of "clusteringM" is selected as "PhenoGraph".}

\item{ncluster}{Integer, the number of clusters for meta clustering in FlowSOM.
\if{html}{\out{<br>}}Only needed when the argument of "clusteringM" is selected as "FlowSOM".}

\item{Ca_metric}{Character, the assessing metric under Criterion Ca for the "CSI" study type, including "AUC and "F1 score".
\if{html}{\out{<br>}}\strong{AUC}: a metric that measures the ability of a model to distinguish between classes by calculating the area under the Receiver Operating Characteristic (ROC) curve, where a higher AUC indicates better model performance in terms of classification.
\if{html}{\out{<br>}}\strong{F1 score}: a metric that evaluates a model’s accuracy by calculating the harmonic mean of precision and recall, where a higher F1 Score indicates a better balance between the model’s ability to identify true positives and avoid false positives.}

\item{Cb_metric}{Character, the assessing metric under Criterion Cb for the "CSI" study type, including "Silhouette coefficient (SC)", "Xie-Beni index (XB)", "Calinski-Harabasz index (CH)", "Davies-Bouldin index (DB)", "purity" and "Rand index (RI)".
\if{html}{\out{<br>}}\strong{Silhouette coefficient (SC)}: a metric that measures how similar a sample is to its own cluster compared to other clusters by calculating the difference between the mean intra-cluster distance and the mean nearest-cluster distance.
\if{html}{\out{<br>}}\strong{Xie-Beni index (XB)}: a metric that assesses the compactness and separation of clusters by calculating the ratio of the sum of squared distances of samples to their cluster centroids over the minimum distance squared between any two cluster centroids.
\if{html}{\out{<br>}}\strong{Calinski-Harabasz index (CH)}: a metric that evaluates clustering quality by assessing the ratio of between-cluster dispersion to within-cluster dispersion, aiming to maximize inter-cluster distances while minimizing intra-cluster distances.
\if{html}{\out{<br>}}\strong{Davies-Bouldin index (DB)}: a metric that evaluates clustering results by calculating the average similarity ratio between each cluster and its most similar other cluster, where lower values indicate better clustering performance.
\if{html}{\out{<br>}}\strong{purity}: a metric that measures the extent to which each cluster contains data points primarily from a single class by identifying the dominant class in each cluster and calculating the proportion of these dominant class counts to the total number of samples.
\if{html}{\out{<br>}}\strong{Rand index (RI)}: a metric that evaluates the similarity between clustering results and ground truth labels by calculating the ratio of correctly classified sample pairs (both true positives and true negatives) to the total number of pairs.}

\item{Cc_metric}{Character, the assessing metric under Criterion Cc for the "CSI" study type, including "relative weighted consistency (CWrel)" and "consistency score (CS)".
\if{html}{\out{<br>}}\strong{relative weighted consistency (CWrel)}: a metric that evaluates the consistency of clustering with respect to a reference partition by weighting the consistency of data point assignments based on their relative importance, with a higher CWrel indicating a stronger alignment with the reference clustering structure.
\if{html}{\out{<br>}}\strong{consistency score (CS)}: a metric that assesses the stability of clustering results by measuring how consistently data points are assigned to the same clusters across different runs or variations of the clustering process, with a higher CS indicating more stable and reliable clustering.}

\item{ntop}{Integer, the number of the most differentially expressed markers that are truncated for calculating the CWrel value.
\if{html}{\out{<br>}}Only needed when the argument of "Cc_metric" is selected as "relative weighted consistency (CWrel)". This value must be less than the number of your selected markers.}

\item{DEP}{Character, the absolute filepath of the CSV file including the differentially expressed proteins used as the prior knowledge for the fourth criterion.
\if{html}{\out{<br>}}It is a table of one column without the column name, each table cell includes one protein typically in the format of "channel description (channel name)", for example: "CD20(FITC.A)".}

\item{marker_path}{Character, the absolute file path of the CSV or XLSX file containing the markers for cell type annotation, and detailed format requirements can be found on the website https://github.com/idrblab/ANPELA.}

\item{known_celltype_path}{Character, the absolute file path of the CSV file containing the gold-standard cell type annotation results, with the first column being cell IDs and the second column being cell types, and these cell types should correspond pre-cisely to those in the ‘marker_path’ file.}

\item{save_processed_res}{Character, the format of the data processing output files. "no" denotes that the results would not be saved. "one_folder" denotes that successfully processed results will be saved as separate RData files in the "process_res" folder. "one_RData" denotes that all processed results will be saved as one RData file in the "process_res" folder.}

\item{savepath}{Character, the absolute path of the folder which will store the assessment results.}

\item{cores}{Integer, the number of CPU cores to be employed for performing parallel computing.
\if{html}{\out{<br>}}To avoid memory explosion due to parallel computing, the default is the largest integers not greater than half of the number of CPU cores on the current host.}
}
\value{
The \strong{assess_res} folder stores the assessment output file named \code{name}\strong{_assess.RData}, which contains 2 lists, "table" and "table2", providing the raw scores for different assessment criteria and performance assessment levels categorized by thresholds, respectively.
\if{html}{\out{<br>}}In addition, the file \strong{log.txt} is also generated simultaneously, recording the processing details.
}
\description{
CSIassess() assesses processing performance of all workflows which are used while running the function ""FCprocess" or "MCprocess" based on comprehensive criteria (each with distinct underlying theories) from the perspective of CSI studies.
}
\examples{
\donttest{
}
}
