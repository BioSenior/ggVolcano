#' Add regulate column to the DEG result
#'
#' @param data a data of DEG result.
#' @param log2FC_name the column name corresponding to log2FC
#' @param fdr_name the column name corresponding to fdr(p-adjust)
#' @param log2FC the threshold value of log2FC.
#' @param fdr the threshold value of fdr.
#'
#' @return
#' @export
#'
#' @examples
#' # load the data
#' data(deg_data)
#' head(deg_data)
#'
#' data_new <- add_regulate(deg_data, log2FC_name = "log2FoldChange",
#'                          fdr_name = "padj",log2FC = 1, fdr = 0.05)
#' head(data_new)
add_regulate <- function(data, log2FC_name = "log2FoldChange",
                         fdr_name = "padj",
                         log2FC = 1, fdr = 0.05){
  colnames(data)[colnames(data) == log2FC_name] <- "log2FoldChange"
  colnames(data)[colnames(data) == fdr_name] <- "padj"

  data$regulate <- "Normal"

  loc_up <- intersect(which(data$log2FoldChange > log2FC),
                      which(data$padj<fdr))
  loc_down <- intersect(which(data$log2FoldChange < (-log2FC)),
                        which(data$padj<fdr))

  data$regulate[loc_up] <- "Up"
  data$regulate[loc_down] <- "Down"
  return(data)
}

#' A volcano plot
#'
#' @param data a data of DEG result.
#' @param x the column name corresponding to the x-axis, defualt is "log2FoldChange"
#' @param y the column name corresponding to the y-axis, defualt is "padj".
#' @param pointSize size of the point.
#' @param pointShape shape of the point.
#' @param fills a vector containing the fill color of the point.
#' @param colors a vector containing the stroke color of the point.
#' @param x_lab label of x-axis.
#' @param y_lab label of y-axis.
#' @param legend_title title of the legend.
#' @param add_line a logical value, whether to add a dashed line, defult is TRUE.
#' @param label the column name corresponding to the label.
#' @param label_number how many gene labels you want to show in the plot.
#' @param output a logical value, whether to save the image, defult is TRUE.
#' @param filename if the output = TRUE, please set a filename.
#'
#' @return
#' @export
#'
#' @examples
#' # load the data
#' data(deg_data)
#' data <- add_regulate(deg_data, log2FC_name = "log2FoldChange",
#'                      fdr_name = "padj",log2FC = 1, fdr = 0.05)
#' # plot
#' ggvolcano(data, x = "log2FoldChange", y = "padj",
#'           label = "row", label_number = 10, output = FALSE)
ggvolcano <- function(data,
                      x = "log2FoldChange", y = "padj",
                      pointSize = 1, pointShape = 21,
                      fills = c("#00AFBB","#999999","#FC4E07"),
                      colors = c("#00AFBB","#999999","#FC4E07"),
                      x_lab = NULL, y_lab = NULL, legend_title = NULL,
                      add_line = TRUE, label = "row", label_number = 10,
                      output = TRUE, filename = "volcano_plot"){

  colnames(data)[colnames(data) == x] <- "x"
  colnames(data)[colnames(data) == y] <- "y"
  colnames(data)[colnames(data) == label] <- "geneName"

  data$label <- rep("",nrow(data))
  data$label[order(data$y)[1:label_number]] <- data$geneName[order(data$x)[1:label_number]]

  p <- ggplot(data, aes(x, -log10(y), color = regulate))
  p <- p +
    geom_point(aes(fill = regulate), size = pointSize, shape = pointShape) +
    scale_fill_manual(values=fills) +
    scale_color_manual(values=colors) +
    theme_bw() +
    theme(title = element_text(size = 15),
          text = element_text(size = 15),
          legend.position = c(0.01, 0.99),
          legend.justification = c(0, 1),
          legend.background = element_blank()) +
    # 添加标签：
    geom_text_repel(aes(label = label), size=3, max.overlaps = 100, key_glyph = draw_key_point)+
    # 设置部分图例不显示：
    guides(fill = guide_legend(title = legend_title %||% "Regulate"),
           color = guide_legend(title = legend_title %||% "Regulate"))+
    # 修改坐标轴：
    labs(x = x_lab %||% TeX("$Log_2 \\textit{FC}$"),
         y = y_lab %||% TeX("$-Log_{10} \\textit{FDR} $"))

  if (add_line == TRUE) {
    p <- p +
      geom_vline(xintercept = c(-1,1), linetype ="dashed") +
      geom_hline(yintercept = -log10(0.05), linetype ="dashed")
  }

  # 保存图片：
  if (output == TRUE) {
    ggsave(paste0(filename, ".pdf"), plot = p, height = 5, width = 6)
  }

  return(p)
}


#' A gradual volcano plot
#'
#'
#' @param data a data of DEG result.
#' @param x the column name corresponding to the x-axis, defualt is "log2FoldChange"
#' @param y the column name corresponding to the y-axis, defualt is "padj".
#' @param pointSizeRange a two elements vector specifing the range of the point size.
#' @param fills a vector containing the fill color of the point.
#' @param colors a vector containing the stroke color of the point.
#' @param x_lab label of x-axis.
#' @param y_lab label of y-axis.
#' @param legend_title title of the legend.
#' @param add_line a logical value, whether to add a dashed line, defult is TRUE.
#' @param label the column name corresponding to the label.
#' @param label_number how many gene labels you want to show in the plot.
#' @param output a logical value, whether to save the image, defult is TRUE.
#' @param filename if the output = TRUE, please set a filename.
#'
#' @return
#' @export
#'
#' @examples
#' # load the data
#' data(deg_data)
#' data <- add_regulate(deg_data, log2FC_name = "log2FoldChange",
#'                      fdr_name = "padj",log2FC = 1, fdr = 0.05)
#'
#' # plot
#' gradual_volcano(deg_data, x = "log2FoldChange", y = "padj",
#'              label = "row", label_number = 10, output = FALSE)
gradual_volcano <- function(data,
                            x = "log2FoldChange", y = "padj",
                            pointSizeRange = c(0.5, 4),
                            fills = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"),
                            colors = c("#17194e","#68bfe7","#f9ed36","#a22f27","#211f1f"),
                            x_lab = NULL, y_lab = NULL, legend_title = NULL,
                            add_line = TRUE, label = "row", label_number = 10,
                            output = TRUE, filename = "volcano_plot"
                            ){
  colnames(data)[colnames(data) == x] <- "x"
  colnames(data)[colnames(data) == y] <- "y"
  colnames(data)[colnames(data) == label] <- "geneName"

  data$label <- rep("",nrow(data))
  data$label[order(data$y)[1:label_number]] <- data$geneName[order(data$x)[1:label_number]]

  p <- ggplot(data,aes(x, -log10(y)))

  if (add_line == TRUE) {
    p <- p +
      # 横向水平参考线：
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
      # 纵向垂直参考线：
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#999999")
  }

  p <- p +
    # 散点图:
    geom_point(aes(size=-log10(y), fill = -log10(y), color =  -log10(y)), shape = 21)+
    # 指定颜色渐变模式：
    scale_fill_gradientn(values = seq(0,1,0.2),
                         colors = fills)+
    scale_color_gradientn(values = seq(0,1,0.2),
                          colors = colors)+
    # 指定散点大小渐变模式：
    scale_size_continuous(range = pointSizeRange)+
    # 主题调整：
    theme_bw()+
    # 调整主题和图例位置：
    theme(panel.grid = element_blank(),
          legend.position = c(0.01,0.7),
          legend.justification = c(0,1),
          legend.background = element_blank()
    )+
    # 设置部分图例不显示：
    guides(fill = guide_colourbar(title = legend_title %||% "-Log10_q-value"),
           color = "none",
           size = "none")+
    # 添加标签：
    geom_text_repel(aes(label = label, color = -log10(y)), size=3, max.overlaps = 10)+
    # 修改坐标轴：
    labs(x = x_lab %||% TeX("$Log_2 \\textit{FC}$"),
         y = y_lab %||% TeX("$-Log_{10} \\textit{FDR} $"))

  # 保存图片：
  if (output == TRUE) {
    ggsave(paste0(filename, ".pdf"), plot = p, height = 5, width = 6)
  }
  return(p)
}


#' A term volcano plot
#'
#' @param data a data of DEG result.
#' @param term_data a two columns dataframe containing some genes' GO terms information.
#' @param x the column name corresponding to the x-axis, defualt is "log2FoldChange"
#' @param y the column name corresponding to the y-axis, defualt is "padj".
#' @param normal_point_color the color of normal points. Defult is "#999999".
#' @param normal_point_size the size of normal points. Defult is 1.
#' @param deg_point_color the stroke color of these deg points, defualt is "black".
#' @param deg_point_fill a named vector containing the fill color of these deg points. Make sure that the names of the elements in the vector correspond to the pathways and colors.
#' @param deg_point_size the size of these deg points. Defult is 2.
#' @param legend_background_fill the fill color of legend background.
#' @param legend_title title of the legend.
#' @param x_lab label of x-axis.
#' @param y_lab label of y-axis.
#' @param add_line a logical value, whether to add a dashed line, defult is TRUE.
#' @param label the column name corresponding to the label.
#' @param label_number how many gene labels you want to show in the plot.
#' @param output a logical value, whether to save the image, defult is TRUE.
#' @param filename if the output = TRUE, please set a filename.
#'
#'
#' @return
#' @export
#'
#' @examples
#' # load the data
#' data(deg_data)
#' data(term_data)
#' data <- add_regulate(deg_data, log2FC_name = "log2FoldChange",
#'                      fdr_name = "padj",log2FC = 1, fdr = 0.05)
#'
#' # plot
#' term_volcano(deg_data, term_data, x = "log2FoldChange", y = "padj",
#'              label = "row", label_number = 10, output = FALSE)
term_volcano <- function(data, term_data,
                         x = "log2FoldChange", y = "padj",
                         normal_point_color = "#999999",
                         normal_point_size = 1,
                         deg_point_color = "black",
                         deg_point_fill = c(dendritic="#49c2c6", "ion transport."="#fbcbcc",
                                            metabolic="#eef0ac",myelin="#b1daa7",
                                            synaptic="#d0d0a0"),
                         deg_point_size = 2,
                         legend_background_fill = "#fefde2",legend_title = NULL,
                         x_lab = NULL, y_lab = NULL,
                         add_line = TRUE,
                         label = "row", label_number = 10,
                         output = TRUE, filename = "volcano_plot"){

  colnames(term_data) <- c("geneName", "term")
  colnames(data)[colnames(data) == x] <- "x"
  colnames(data)[colnames(data) == y] <- "y"
  colnames(data)[colnames(data) == label] <- "geneName"

  # 添加GO一列：
  data$GO_term <- "others"
  term_data <- term_data[term_data$geneName %in% data$geneName,]
  data[term_data$geneName,]$GO_term <- term_data$term

  # 计算上调下调数目：
  Down_num <- length(which(data$y < 0.05 & data$x < 0))
  Up_num <- length(which(data$y < 0.05 & data$x > 0))

  # 设定原始散点颜色：
  color <- rep(normal_point_color, nrow(data))

  # 选取p值最显著的label_number个加上标签：
  data$label <- rep("",nrow(data))
  data$label[order(data$y)[1:label_number]] <- data$geneName[order(data$x)[1:label_number]]

  p <- ggplot(data[which(data$GO_term!="others"),],
         aes(x,-log10(y),fill = GO_term))+
    geom_point(data=data[which(data$GO_term=="others"),],
               aes(x,-log10(y)),
               size = normal_point_size, color=normal_point_color) +
    # 彩色散点：
    geom_point(size = deg_point_size, shape=21, color=deg_point_color) +
    scale_fill_manual(values = deg_point_fill) +
    labs(x = x_lab %||% TeX("$Log_2 \\textit{FC}$"),
         y = y_lab %||% TeX("$-Log_{10} \\textit{FDR} $"))+
    theme(title = element_text(size = 15), text = element_text(size = 15)) +
    theme_bw() +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          legend.position = c(0.01, 0.99),
          legend.justification = c(0, 1),
          # 图例大框颜色：
          legend.background = element_rect(
            fill = legend_background_fill, # 填充色
            colour = "black", # 框线色
            # 线条宽度
            size = 0.2),
          # 图例符号颜色：
          legend.key = element_rect(fill = legend_background_fill),
          # 调整图例大小：
          legend.key.size = unit(12, "pt"))+
    # 添加gene标签：部分标签没有显示是因为重叠，可以修改max.overlaps值；
    geom_text_repel(data = data, aes(x,-log10(y), label = label), size=3, max.overlaps = 100)+
    guides(fill = guide_legend(title = legend_title %||% "GO terms"))

  if (add_line == TRUE) {
    p <- p +
      geom_vline(xintercept = 0, linetype ="longdash") +
      geom_hline(yintercept = -log10(0.05), linetype ="longdash")
  }

  # 保存图片：
  if (output == TRUE) {
    ggsave(paste0(filename, ".pdf"), plot = p, height = 5, width = 6)
  }
  return(p)
}
