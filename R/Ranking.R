
#' @title Ranking
#' @description Ranking ranks all workflows assessed by the function "CSIassess" and "PTIassess" based on collective consideration of values and grades (classified by well-defined cutoffs) under each criterion.
#'
#' @param name Character, the filename of the overall ranking data and figure file.
#' @param data Character, the R object resulting from the function "Assess" "CSIassess" or "PTIassess", or obtained by loading from the resulting RData file of these funcitons.
#' @param savepath Character, the absolute path of the folder which will store the overall ranking data and figure file.
#'
#' @return A CSV file named "_Ranking_Table.csv", recording the overall ranking and the values of criteria. <br>A PDF file named "_Ranking_Figure.pdf", helping users better understand the differences among various data processing workflows, where the different colors represent different performance assessment levels: green indicates "good," and red indicates "poor".
#'
#' @export
#'
#' @examples
#' \donttest{
#' }

Ranking <- function(data, name = "result", savepath = "./ANPELA_res") {
  if (any(!is.na(data$table[, 4]))) {
    table <- data$table
    table2 <- data$table2
  } else {
    table <- data$table[, -4]
    table2 <- data$table2[, -4]
  }

  Color_score <- apply(table2, 1, sum)
  Value_score <- apply(table, 1, sum)
  OverallRank <- order(Color_score, Value_score, decreasing = TRUE)
  table_res <- table[OverallRank, ]

  csvresult <- cbind(Rank.OverallRank = 1:nrow(table), Value = table_res)
  write.csv(csvresult, file = paste0(savepath,"/", name, "_Ranking_Table.csv"))
  pheatmapresult <- table2[OverallRank, ]
  rownames(pheatmapresult) <- paste0("Rank", 1:nrow(pheatmapresult), " ", rownames(pheatmapresult))

  pdf(paste0(savepath,"/", name, "_Ranking_Figure.pdf"), height = ceiling(nrow(pheatmapresult)/100*12.5)+2)
  pheatmap::pheatmap(pheatmapresult, cluster_rows = FALSE, cluster_cols = FALSE,
                     angle_col = "90", annotation_legend = FALSE, show_rownames = TRUE, show_colnames = TRUE,
                     border_color = "white", cellheight = 9, cellwidth = 30, fontsize_row = 6, fontsize_col = 8,
                     color = c("#CA2125", "#008237"), legend = FALSE,
                     breaks = c(1,4,10),
                     legend_breaks = c(seq(1, 10, by = 1.5)))
  dev.off()
}
