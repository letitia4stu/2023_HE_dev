rm(list=ls()) 
setwd('/Users/lichaoran/dataprocess/2016-nature-HE/HE_bowtie_new_version')
library('dplyr')
library('myTAI')
library('edgeR')
library(ggplot2)
#Expression genes in every time points with PS value
PhyloExpressionSet=read.csv('/Users/lichaoran/Downloads/he-genage-new/global-dev-exp-genage.csv',check.names = F)
is.ExpressionSet(PhyloExpressionSet)
#modify function
PlotSignaturechaa <- function (ExpressionSet, measure = "TAI", TestStatistic = "FlatLineTest", 
                               modules = NULL, permutations = 1000, lillie.test = FALSE, 
                               p.value = TRUE, shaded.area = FALSE, custom.perm.matrix = NULL, 
                               xlab = "Ontogeny", ylab = "Transcriptome Index", main = "", 
                               lwd = 4, alpha = 0.1, y.ticks = 10) 
{
  if (!is.element(measure, c("TAI", "TDI", "TPI"))) 
    stop("Measure '", measure, "' is not available for this function. Please specify a measure supporting by this function.", 
         call. = FALSE)
  if (!is.element(TestStatistic, c("FlatLineTest", "ReductiveHourglassTest", 
                                   "EarlyConservationTest", "ReverseHourglassTest"))) 
    stop("Please choose a 'TestStatistic' that is supported by this function. E.g. TestStatistic = 'FlatLineTest', TestStatistic = 'ReductiveHourglassTest', TestStatistic = 'EarlyConservationTest', TestStatistic = 'ReverseHourglassTest'.", 
         call. = FALSE)
  cat("Plot signature: '", measure, "' and test statistic: '", 
      TestStatistic, "' running ", permutations, " permutations.")
  cat("\n")
  Stage <- NULL
  if (measure == "TAI") {
    TI <- tibble::tibble(Stage = names(TAI(ExpressionSet)), 
                         TI = TAI(ExpressionSet))
  }
  if (measure == "TDI") {
    TI <- tibble::tibble(Stage = names(TDI(ExpressionSet)), 
                         TI = TDI(ExpressionSet))
  }
  if (measure == "TPI") {
    TI <- tibble::tibble(Stage = names(TPI(ExpressionSet)), 
                         TI = TPI(ExpressionSet))
  }
  bm <- tibble::as_tibble(bootMatrix(ExpressionSet = ExpressionSet, 
                                     permutations = permutations))
  if (TestStatistic == "FlatLineTest") {
    if (is.null(custom.perm.matrix)) {
      resList <- FlatLineTest(ExpressionSet = ExpressionSet, 
                              permutations = permutations)
    }
    else if (!is.null(custom.perm.matrix)) {
      resList <- FlatLineTest(ExpressionSet = ExpressionSet, 
                              custom.perm.matrix = custom.perm.matrix)
    }
  }
  if (TestStatistic == "ReductiveHourglassTest") {
    if (lillie.test) {
      if (is.null(custom.perm.matrix)) {
        resList <- ReductiveHourglassTest(ExpressionSet = ExpressionSet, 
                                          modules = modules, permutations = permutations, 
                                          lillie.test = TRUE)
      }
      else if (!is.null(custom.perm.matrix)) {
        resList <- ReductiveHourglassTest(ExpressionSet = ExpressionSet, 
                                          modules = modules, lillie.test = TRUE, custom.perm.matrix = custom.perm.matrix)
      }
    }
    if (!lillie.test) {
      if (is.null(custom.perm.matrix)) {
        resList <- ReductiveHourglassTest(ExpressionSet = ExpressionSet, 
                                          modules = modules, permutations = permutations, 
                                          lillie.test = FALSE)
      }
      else if (!is.null(custom.perm.matrix)) {
        resList <- ReductiveHourglassTest(ExpressionSet = ExpressionSet, 
                                          modules = modules, lillie.test = FALSE, custom.perm.matrix = custom.perm.matrix)
      }
    }
  }
  if (TestStatistic == "ReverseHourglassTest") {
    if (lillie.test) {
      if (is.null(custom.perm.matrix)) {
        resList <- ReverseHourglassTest(ExpressionSet = ExpressionSet, 
                                        modules = modules, permutations = permutations, 
                                        lillie.test = TRUE)
      }
      else if (!is.null(custom.perm.matrix)) {
        resList <- ReverseHourglassTest(ExpressionSet = ExpressionSet, 
                                        modules = modules, lillie.test = TRUE, custom.perm.matrix = custom.perm.matrix)
      }
    }
    if (!lillie.test) {
      if (is.null(custom.perm.matrix)) {
        resList <- ReverseHourglassTest(ExpressionSet = ExpressionSet, 
                                        modules = modules, permutations = permutations, 
                                        lillie.test = FALSE)
      }
      else if (!is.null(custom.perm.matrix)) {
        resList <- ReverseHourglassTest(ExpressionSet = ExpressionSet, 
                                        modules = modules, lillie.test = FALSE, custom.perm.matrix = custom.perm.matrix)
      }
    }
  }
  if (TestStatistic == "EarlyConservationTest") {
    if (lillie.test) {
      if (is.null(custom.perm.matrix)) {
        resList <- EarlyConservationTest(ExpressionSet = ExpressionSet, 
                                         modules = modules, permutations = permutations, 
                                         lillie.test = TRUE)
      }
      else if (!is.null(custom.perm.matrix)) {
        resList <- EarlyConservationTest(ExpressionSet = ExpressionSet, 
                                         modules = modules, lillie.test = TRUE, custom.perm.matrix = custom.perm.matrix)
      }
    }
    if (!lillie.test) {
      if (is.null(custom.perm.matrix)) {
        resList <- EarlyConservationTest(ExpressionSet = ExpressionSet, 
                                         modules = modules, permutations = permutations, 
                                         lillie.test = FALSE)
      }
      else if (!is.null(custom.perm.matrix)) {
        resList <- EarlyConservationTest(ExpressionSet = ExpressionSet, 
                                         modules = modules, lillie.test = FALSE, custom.perm.matrix = custom.perm.matrix)
      }
    }
  }
  pval <- resList$p.value
  pval <- format(pval, digits = 3)
  TI.ggplot <- ggplot2::ggplot(TI, ggplot2::aes(x = as.numeric(Stage), 
                                                y = TI, group = 1)) + ggplot2::geom_ribbon(ggplot2::aes(ymin = TI - 
                                                                                                          apply(bm, 2, stats::sd), ymax = TI + apply(bm, 2, stats::sd)), 
                                                                                           alpha = alpha) + ggplot2::geom_line(lwd = lwd) + ggplot2::theme_minimal() + 
    ggplot2::labs(x = xlab, y = ylab, title = main) + theme_bw()+theme(legend.background = element_blank()) + theme(legend.title = element_blank())+ggplot2::theme(title = ggplot2::element_text(size = 18, 
                                                                                                   face = "bold"), legend.title = ggplot2::element_text(size = 18, 
                                                                                                                                                        face = "bold"), legend.text = ggplot2::element_text(size = 18, 
                                                                                                                                                                                                            face = "bold"), axis.title = ggplot2::element_text(size = 18, 
                                                                                                                                                                                                                                                               face = "bold"), axis.text.y = ggplot2::element_text(size = 18, 
                                                                                                                                                                                                                                                                                                                   face = "bold"), axis.text.x = ggplot2::element_text(size = 18, 
                                                                                                                                                                                                                                                                                                                                                                       face = "bold"), panel.background = ggplot2::element_blank(), 
                                                                     strip.text.x = ggplot2::element_text(size = 18, colour = "black", 
                                                                                                          face = "bold")) + ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = y.ticks))+ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = y.ticks))
  if ((TestStatistic == "FlatLineTest") && p.value) {
    TI.ggplot <- TI.ggplot + ggplot2::annotate("text", x = 2, 
                                               y = max(TI$TI) + (max(TI$TI)/30), label = paste0("p_flt = ", 
                                                                                                pval), size = 6)
    cat("Significance status of signature: ", ifelse(as.numeric(pval) <= 
                                                       0.05, "significant.", "not significant (= no evolutionary signature in the transcriptome)."))
    cat("\n")
  }
  if (TestStatistic == "ReductiveHourglassTest") {
    if (p.value) {
      TI.ggplot <- TI.ggplot + ggplot2::annotate("text", 
                                                 x = 2, y = max(TI$TI) + (max(TI$TI)/30), label = paste0("p_rht = ", 
                                                                                                         pval), size = 6)
    }
    if (shaded.area) {
      TI.ggplot <- TI.ggplot + ggplot2::geom_rect(data = TI, 
                                                  ggplot2::aes(xmin = modules[[2]][1], xmax = modules[[2]][length(modules[[2]])], 
                                                               ymin = min(TI) - (min(TI)/50), ymax = Inf), 
                                                  fill = "#4d004b", alpha = alpha)
    }
    stage.names <- names(ExpressionSet)[3:ncol(ExpressionSet)]
    cat("Modules: \n early = {", paste0(stage.names[modules[[1]]], 
                                        " "), "}", "\n", "mid = {", paste0(stage.names[modules[[2]]], 
                                                                           " "), "}", "\n", "late = {", paste0(stage.names[modules[[3]]], 
                                                                                                               " "), "}")
    cat("\n")
    cat("Significance status of signature: ", ifelse(as.numeric(pval) <= 
                                                       0.05, "significant.", "not significant (= no evolutionary signature in the transcriptome)."))
  }
  if (TestStatistic == "ReverseHourglassTest") {
    if (p.value) {
      TI.ggplot <- TI.ggplot + ggplot2::annotate("text", 
                                                 x = 2, y = max(TI$TI) + (max(TI$TI)/30), label = paste0("p_reverse_hourglass = ", 
                                                                                                         pval), size = 6)
    }
    if (shaded.area) {
      TI.ggplot <- TI.ggplot + ggplot2::geom_rect(data = TI, 
                                                  ggplot2::aes(xmin = modules[[2]][1], xmax = modules[[2]][length(modules[[2]])], 
                                                               ymin = min(TI) - (min(TI)/50), ymax = Inf), 
                                                  fill = "#4d004b", alpha = alpha)
    }
    stage.names <- names(ExpressionSet)[3:ncol(ExpressionSet)]
    cat("Modules: \n early = {", paste0(stage.names[modules[[1]]], 
                                        " "), "}", "\n", "mid = {", paste0(stage.names[modules[[2]]], 
                                                                           " "), "}", "\n", "late = {", paste0(stage.names[modules[[3]]], 
                                                                                                               " "), "}")
    cat("\n")
    cat("Significance status of signature: ", ifelse(as.numeric(pval) <= 
                                                       0.05, "significant.", "not significant (= no evolutionary signature in the transcriptome)."))
  }
  if (TestStatistic == "EarlyConservationTest") {
    if (p.value) {
      TI.ggplot <- TI.ggplot + ggplot2::annotate("text", 
                                                 x = 2, y = max(TI$TI) + (max(TI$TI)/30), label = paste0("p_ect = ", 
                                                                                                         pval), size = 6)
    }
    if (shaded.area) {
      TI.ggplot <- TI.ggplot + ggplot2::geom_rect(data = TI, 
                                                  ggplot2::aes(xmin = modules[[2]][1], xmax = modules[[2]][length(modules[[2]])], 
                                                               ymin = min(TI) - (min(TI)/50), ymax = Inf), 
                                                  fill = "#4d004b", alpha = alpha * 0.5)
    }
    stage.names <- names(ExpressionSet)[3:ncol(ExpressionSet)]
    cat("Modules: \n early = {", paste0(stage.names[modules[[1]]], 
                                        " "), "}", "\n", "mid = {", paste0(stage.names[modules[[2]]], 
                                                                           " "), "}", "\n", "late = {", paste0(stage.names[modules[[3]]], 
                                                                                                               " "), "}")
    cat("\n")
    cat("Significance status of signature: ", ifelse(as.numeric(pval) <= 
                                                       0.05, "significant.", "not significant (= no evolutionary signature in the transcriptome)."))
  }
  TI.ggplot <- TI.ggplot + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, 
                                                                              vjust = 1, hjust = 1))
  return(TI.ggplot)
}


PlotSignaturechaa( ExpressionSet = PhyloExpressionSet,
                   measure       = "TAI", 
                   permutations  = 1000,
                   modules = list(early=1:9,
                                  mid=10:20,
                                  late=21:29),
                   TestStatistic = "ReductiveHourglassTest",
                   xlab          = "Developmental time point (h)")
ReductiveHourglassTest(PhyloExpressionSet,
                       modules = list(early = 1:12, mid = 13:20, late = 21:29), 
                       permutations = 1000)
pcha=PlotSignaturechaa( ExpressionSet = PhyloExpressionSet,
                        measure       = "TAI", 
                        permutations  = 1000,
                        modules = list(early=1:12,
                                       mid=13:20,
                                       late=21:29),
                        TestStatistic = "ReductiveHourglassTest",
                        xlab          = "Developmental time point (h)")


at=pcha$data
write.csv(at,"29-timepoints-tai.csv")

cowplot::save_plot("55-29_exp_TAI.pdf",
                   pcha,
                   base_height = 10,
                   base_width = 14 )

