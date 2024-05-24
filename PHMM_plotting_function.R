plot_Results <- function(modelResults, numPlots, klen){
  
  # Split by rank for % correct ids
  percentCorrect_group <- modelResults$ResultsData$`CorrectID% - Group` %>%
    group_by(rank) %>%
    group_split()
  
  # Split by rank for TP,FP,TN,FN
  confidenceVal_group <- modelResults$ResultsData$`TP,FP,TN,FN - Group` %>%
    group_by(rank) %>%
    group_split()
  
  # another vector for klen (because dplyr is getting confused)
  klen_vec <- klen
  
  # Creating directory for model results plots
  dirP <- paste0(dir, "/ModelResults")
  dir.create(dirP, showWarnings = FALSE)
  
  # Code segment for plotting % CorrectIDs per group
  
  #*** Changed all if statements - was getting error before
  # If running order level or all levels
  if("Order" %in% subclassRank){
    
    # Finding average percent correct ids across all folds and then sorting by average percent correct ids for plotting
    percentCorrect_groupAvgOrd <- suppressWarnings(percentCorrect_group[[3]] %>% 
                                                     group_by(Group, klen) %>% 
                                                     dplyr::summarise(rank=rank, `CorrectId%` = mean(`CorrectId%`), n=sum(n), fold="Avg") %>% 
                                                     mutate(`Fold&klength` = paste0(fold, ",", klen, sep="")) %>% 
                                                     arrange(desc(`CorrectId%`)) %>% 
                                                     dplyr::distinct() %>%
                                                     dplyr::relocate(any_of(c("rank", "Group", "fold", "klen", "CorrectId%", "n", "Fold&klength"))))
    
    sortOrd <- percentCorrect_groupAvgOrd[,"Group"] %>% rownames_to_column('sortVar')
    
    # percentCorrect_groupAvgOrd <- rbind(percentCorrect_group[[3]], percentCorrect_groupAvgOrd)
    percentCorrect_groupAvgOrd <- percentCorrect_groupAvgOrd %>%
      mutate(Group = factor(Group, levels = rev(as.character(unique(sortOrd$Group)))))
    
    # Divide into two plots
    # Count up which groups fall into one of two groups for family based on the number of orders and 
    # number of individual families for each order
    alphaGroup_ord <- percentCorrect_groupAvgOrd %>%
      distinct() %>%
      group_by(Group) %>%
      tally() 
    alphaGroup_ord$plotNum <- cumsum(alphaGroup_ord$n) %/% floor(sum(alphaGroup_ord$n) / numPlots[1]) + 1
    alphaGroup_ord <- alphaGroup_ord %>% select(!n)
    alphaGroup_ord$plotNum <- ifelse(alphaGroup_ord$plotNum == max(alphaGroup_ord$plotNum), alphaGroup_ord$plotNum - 1, alphaGroup_ord$plotNum)
    
    # Merge with percentCorrect_groupAvgFam
    percentCorrect_groupAvgOrd <- merge(percentCorrect_groupAvgOrd, alphaGroup_ord, by="Group")
    
    # Order Heatmap plot
    ordPlots <- list()
    foreach(i=1:length(unique(percentCorrect_groupAvgOrd$plotNum))) %do% {
      ordPlots[[i]] <- ggplot(as.data.frame(percentCorrect_groupAvgOrd[which(percentCorrect_groupAvgOrd$plotNum == i),]), aes(y = Group, x = `klen`)) + 
        geom_tile(aes(fill=`CorrectId%`)) +
        geom_text(colour="darkorchid1", size=6,aes(label=n), fontface = "bold") + 
        scale_fill_viridis(discrete=FALSE) +
        theme_dark(base_size = 20) +
        guides(fill = guide_colourbar(barwidth = 0.5,
                                      barheight = 20)) +
        # ggtitle(classRank) +
        theme(panel.grid.major = element_blank()) + 
        ylab("Order") +
        xlab("klength")
    }
    foreach(i=1:length(unique(percentCorrect_groupAvgOrd$plotNum))) %do% {
      ggsave(paste0(dirP, "/ord_percent_correct_", i, ".png", sep=""), ordPlots[[i]], width = 10, height = 15, dpi = 500, bg="white")
    }
    # Add filepath to modelResults
    modelResults$PlotFilepaths$`Order - HeatMap` <- unlist(foreach(i=1:length(unique(percentCorrect_groupAvgOrd$plotNum))) %do% paste0(dirP, "/ord_percent_correct_", i, ".png", sep=""))
    
    # ROC plots using cutpointr package
    confPlotOrd <- list()
    foreach(i=1:length(klen_vec)) %do% { confPlotOrd[[i]] <- roc(data = confidenceVal_group[[3]] %>%
                                                                   dplyr::filter(klen == klen_vec[i]) %>%
                                                                   mutate(truth = as.factor(ifelse(assignedTaxa == Group, 1, 0))) %>%
                                                                   select(!confThreshold) %>%
                                                                   distinct(), x = probScore, class = truth,
                                                                 pos_class = 1, neg_class = 0, direction = ">=") 
    confPlotOrd[[i]] <- confPlotOrd[[i]] %>% mutate(klenauc = paste0(klen_vec[i], ", AUC=", round(cutpointr::auc(confPlotOrd[[i]]), 2), sep="")) }
    confPlotOrd <- rbindlist(confPlotOrd)
    nGroupsOrd <- length(unique(confidenceVal_group[[3]]$Group))
    nSeqsOrd <- length(unique(confidenceVal_group[[3]]$recordID))
    
    ggOrd <- ggplot(data=confPlotOrd, aes(x=fpr, y = tpr, col = klenauc)) + 
      geom_line(size=1.2, alpha=0.8) + 
      scale_colour_viridis_d(name = "klength & AUC") + 
      xlab("Average FPR") + ylab("Average TPR")  +
      ggtitle(paste0(classRank, " - Order, numFolds = ", numFolds, ", nGroups = ", nGroupsOrd, ", nSeqs = ", nSeqsOrd, sep="")) +
      theme(text = element_text(size=16)) +
      annotate("segment", x = 0, xend = 1, y = 0, yend = 1, color="black", linetype="dashed") + 
      theme(panel.spacing = unit(1, "lines"))
    # ggsave(paste0(dirP, "/ord_roc.png"), plot=ggOrd, width = 12, height = 8, dpi = 500, bg="white")
    
    # MCC scores
    confPlotOrd$mcc <- ((confPlotOrd$tn*confPlotOrd$tp) - (confPlotOrd$fp*confPlotOrd$fn)) / sqrt((confPlotOrd$tn+confPlotOrd$fn) * (confPlotOrd$fp+confPlotOrd$tp) * (confPlotOrd$tn+confPlotOrd$fp) * (confPlotOrd$fn+confPlotOrd$tp))
    confPlotOrd$klen <- substr(confPlotOrd$klenauc, 1, 1)
    confPlotOrd$mcc <- ifelse(is.nan(confPlotOrd$mcc) & is.infinite(confPlotOrd$x.sorted), 0, confPlotOrd$mcc)
    confPlotOrd$x.sorted <- ifelse(is.infinite(confPlotOrd$x.sorted), 1, confPlotOrd$x.sorted)
    
    ggOrd2 <- ggplot(data=confPlotOrd, aes(x=x.sorted, y = mcc, col = klen)) + 
      geom_line(size=1.2, alpha=0.8) + 
      scale_colour_viridis_d(name = "klength") + 
      xlab("Confidence Threshold") + ylab("Average MCC")  +
      #ggtitle(paste0(classRank, " - Order, numFolds = ", numFolds, ", nGroups = ", nGroupsOrd, ", nSeqs = ", nSeqsOrd, sep="")) +
      theme(text = element_text(size=16)) +
      theme(panel.spacing = unit(1, "lines"))
    # ggsave(paste0(dirP, "/ord_mcc.png"), plot=ggOrd2, width = 12, height = 8, dpi = 500, bg="white")
    
    # F1 scores
    confPlotOrd$f1 <- (2*confPlotOrd$tp) / ((2 * confPlotOrd$tp) + confPlotOrd$fn + confPlotOrd$fp)
    
    ggOrd3 <- ggplot(data=confPlotOrd, aes(x=x.sorted, y = f1, col = klen)) + 
      geom_line(size=1.2, alpha=0.8) + 
      scale_colour_viridis_d(name = "klength") + 
      xlab("Confidence Threshold") + ylab("Average F1 score")  +
      #ggtitle(paste0(classRank, " - Order, numFolds = ", numFolds, ", nGroups = ", nGroupsOrd, ", nSeqs = ", nSeqsOrd, sep="")) +
      theme(text = element_text(size=16)) +
      theme(panel.spacing = unit(1, "lines"))
    # ggsave(paste0(dirP, "/ord_f1.png"), plot=ggOrd3, width = 12, height = 8, dpi = 500, bg="white")
    
    layoutOrd <- rbind(c(1),
                       c(2),
                       c(3))
    ggCombOrd <- grid.arrange(ggOrd, ggOrd2, 
                              ggOrd3, layout_matrix = layoutOrd)
    ggsave(paste0(dirP, "/ord_roc_mcc_f1.png"), plot=ggCombOrd, width = 10, height = 15, dpi = 500, bg="white")
    
    # Add filepath to modelResults
    modelResults$PlotFilepaths$`Order - ROC,MCC,F1` <- paste0(dirP, "/ord_roc_mcc_f1.png")
  }
  
  # If running family level or all levels
  if("Family" %in% subclassRank){
    
    # Find orders for each family listed in the results for sorting
    dfTaxaOrdFam <- dfTaxa[,c("ordID", "famID")] %>% distinct()
    colnames(dfTaxaOrdFam)[2] <- "Group"
    
    ordGroup_fam <- merge(percentCorrect_group[[1]], dfTaxaOrdFam, by="Group")
    
    # Finding average percent correct ids across all folds and then sorting by average percent correct ids per facet for plotting
    percentCorrect_groupAvgFam <- suppressWarnings(ordGroup_fam %>% 
                                                     group_by(Group, klen) %>% 
                                                     dplyr::summarise(rank=rank, `CorrectId%` = mean(`CorrectId%`), n=sum(n), fold="Avg", ordID=ordID) %>% 
                                                     mutate(`Fold&klength` = paste0(fold, ",", klen, sep="")) %>% 
                                                     arrange(desc(`CorrectId%`)) %>% 
                                                     dplyr::distinct() %>%
                                                     dplyr::relocate(any_of(c("rank", "Group", "fold", "klen", "CorrectId%", "n", "Fold&klength", "ordID"))))
    
    sortFam <- percentCorrect_groupAvgFam[,"Group"] %>% rownames_to_column('sortVar')
    
    # percentCorrect_groupAvgFam <- rbind(ordGroup_fam, percentCorrect_groupAvgFam)
    percentCorrect_groupAvgFam <- percentCorrect_groupAvgFam %>%
      mutate(Group = factor(Group, levels = rev(as.character(unique(sortFam$Group)))))
    
    # Divide into two plots
    # Count up which groups fall into one of two groups for family based on the number of orders and 
    # number of individual families for each order
    alphaGroup_fam <- percentCorrect_groupAvgFam %>%
      distinct() %>%
      group_by(ordID) %>%
      tally() 
    alphaGroup_fam$plotNum <- cumsum(alphaGroup_fam$n) %/% floor(sum(alphaGroup_fam$n) / numPlots[2]) + 1
    alphaGroup_fam <- alphaGroup_fam %>% select(!n)
    alphaGroup_fam$plotNum <- ifelse(alphaGroup_fam$plotNum == max(alphaGroup_fam$plotNum), alphaGroup_fam$plotNum - 1, alphaGroup_fam$plotNum)
    
    # Merge with percentCorrect_groupAvgFam
    percentCorrect_groupAvgFam <- merge(percentCorrect_groupAvgFam, alphaGroup_fam, by="ordID")
    
    # Family Heatmap plots
    famPlots <- list()
    foreach(i=1:length(unique(percentCorrect_groupAvgFam$plotNum))) %do% {
      famPlots[[i]] <- ggplot(as.data.frame(percentCorrect_groupAvgFam[which(percentCorrect_groupAvgFam$plotNum == i),]), aes(`klen`, Group)) + 
        geom_tile(aes(fill=`CorrectId%`)) +
        geom_text(colour="darkorchid2", size=6, aes(label=n), fontface = "bold") + 
        scale_fill_viridis(discrete=FALSE) +
        theme_dark(base_size = 19) +
        # theme(panel.grid.major = element_blank(), 
        #       legend.position = "none") +
        facet_col(vars(ordID), scales = "free_y", space = "free") + 
        # ggtitle(paste0("Order (", plotFamH[[i]]$ordID[1], " - ", plotFamT[[i]]$ordID[1], ")", sep="")) + 
        ylab("Family") +
        xlab("klength")
    }
    # yleft <- textGrob("Family", rot = 90, gp = gpar(fontsize = 20))
    # bottom <- textGrob("kLength", gp = gpar(fontsize = 20))
    # famPlots <- grid.arrange(grobs = famPlots, ncol = 2, left = yleft, bottom = bottom) 
    foreach(i=1:length(unique(percentCorrect_groupAvgFam$plotNum))) %do% {
      ggsave(paste0(dirP, "/fam_percent_correct_", i, ".png", sep=""), famPlots[[i]], width = 10, height = 15, dpi = 500, bg="white")
    }
    # Add filepath to modelResults
    modelResults$PlotFilepaths$`Family - HeatMap` <- unlist(foreach(i=1:length(unique(percentCorrect_groupAvgFam$plotNum))) %do% paste0(dirP, "/fam_percent_correct_", i, ".png", sep=""))
    
    # Family ROC plot
    
    # ROC plots using cutpointr package
    
    # TPR vs FPR
    confPlotFam <- list()
    foreach(i=1:length(klen_vec)) %do% { confPlotFam[[i]] <- roc(data = confidenceVal_group[[1]] %>%
                                                                   dplyr::filter(klen == klen_vec[i]) %>%
                                                                   mutate(truth = as.factor(ifelse(assignedTaxa == Group, 1, 0))) %>%
                                                                   select(!confThreshold) %>%
                                                                   distinct(), x = probScore, class = truth,
                                                                 pos_class = 1, neg_class = 0, direction = ">=") 
    
    confPlotFam[[i]] <- confPlotFam[[i]] %>% mutate(klenauc = paste0(klen_vec[i], ", AUC=", round(cutpointr::auc(confPlotFam[[i]]), 2), sep="")) }
    confPlotFam <- rbindlist(confPlotFam)
    
    nGroupsFam <- length(unique(confidenceVal_group[[1]]$Group))
    nSeqsFam <- length(unique(confidenceVal_group[[1]]$recordID))
    
    ggFam <- ggplot(data=confPlotFam, aes(x=fpr, y = tpr, col = klenauc)) + 
      geom_line(size=1.2, alpha=0.8) + 
      scale_colour_viridis_d(name = "klength & AUC") + 
      xlab("Average FPR") + ylab("Average TPR")  +
      ggtitle(paste0(classRank, " - Family, numFolds = ", numFolds, ", nGroups = ", nGroupsFam, ", nSeqs = ", nSeqsFam, sep="")) +
      theme(text = element_text(size=16)) +
      annotate("segment", x = 0, xend = 1, y = 0, yend = 1, color="black", linetype="dashed") + 
      theme(panel.spacing = unit(1, "lines"))
    # ggsave(paste0(dirP, "/fam_roc.png"), plot=ggFam, width = 12, height = 8, dpi = 500, bg="white")
    
    # MCC scores
    confPlotFam$mcc <- ((confPlotFam$tn*confPlotFam$tp) - (confPlotFam$fp*confPlotFam$fn)) / sqrt((confPlotFam$tn+confPlotFam$fn) * (confPlotFam$fp+confPlotFam$tp) * (confPlotFam$tn+confPlotFam$fp) * (confPlotFam$fn+confPlotFam$tp))
    confPlotFam$klen <- substr(confPlotFam$klenauc, 1, 1)
    confPlotFam$mcc <- ifelse(is.nan(confPlotFam$mcc) & is.infinite(confPlotFam$x.sorted), 0, confPlotFam$mcc)
    confPlotFam$x.sorted <- ifelse(is.infinite(confPlotFam$x.sorted), 1, confPlotFam$x.sorted)
    
    ggFam2 <- ggplot(data=confPlotFam, aes(x=x.sorted, y = mcc, col = klen)) + 
      geom_line(size=1.2, alpha=0.8) + 
      scale_colour_viridis_d(name = "klength") + 
      xlab("Confidence Threshold") + ylab("Average MCC")  +
      # ggtitle(paste0(classRank, " - Family, numFolds = ", numFolds, ", nGroups = ", nGroupsFam, ", nSeqs = ", nSeqsFam, sep="")) +
      theme(text = element_text(size=16)) +
      theme(panel.spacing = unit(1, "lines"))
    # ggsave(paste0(dirP, "/fam_mcc.png"), plot=ggFam2, width = 12, height = 8, dpi = 500, bg="white")
    
    # F1 scores
    confPlotFam$f1 <- (2*confPlotFam$tp) / ((2 * confPlotFam$tp) + confPlotFam$fn + confPlotFam$fp)
    
    ggFam3 <- ggplot(data=confPlotFam, aes(x=x.sorted, y = f1, col = klen)) + 
      geom_line(size=1.2, alpha=0.8) + 
      scale_colour_viridis_d(name = "klength") + 
      xlab("Confidence Threshold") + ylab("Average F1 score")  +
      # ggtitle(paste0(classRank, " - Family, numFolds = ", numFolds, ", nGroups = ", nGroupsFam, ", nSeqs = ", nSeqsFam, sep="")) +
      theme(text = element_text(size=16)) +
      theme(panel.spacing = unit(1, "lines"))
    # ggsave(paste0(dirP, "/fam_f1.png"), plot=ggFam3, width = 12, height = 8, dpi = 500, bg="white")
    
    layoutFam <- rbind(c(1),
                       c(2),
                       c(3))
    ggCombFam <- grid.arrange(ggFam, ggFam2, 
                              ggFam3, layout_matrix = layoutFam)
    ggsave(paste0(dirP, "/fam_roc_mcc_f1.png"), plot=ggCombFam, width = 10, height = 15, dpi = 500, bg="white")
    
    # Add filepath to modelResults
    modelResults$PlotFilepaths$`Family - ROC,MCC,F1` <- paste0(dirP, "/fam_roc_mcc_f1.png")
  }
  
  # If running genus level or all levels
  if("Genus" %in% subclassRank){  
    
    # Find orders for each family listed in the results for sorting
    dfTaxaOrdFamGen <- dfTaxa[,c("ordID", "famID", "genID")] %>% distinct()
    colnames(dfTaxaOrdFamGen)[3] <- "Group"
    
    ordFamGroup_gen <- merge(percentCorrect_group[[2]], dfTaxaOrdFamGen, by="Group") %>%
      mutate(ordFamID = paste0(ordID, ", ", famID, sep=""))
    
    # Finding average percent correct ids across all folds and then sorting by average percent correct ids per facet for plotting
    percentCorrect_groupAvgGen <- suppressWarnings(ordFamGroup_gen %>% 
                                                     group_by(Group, klen) %>% 
                                                     dplyr::summarise(rank=rank, `CorrectId%` = mean(`CorrectId%`), n=sum(n), fold="Avg", ordFamID=ordFamID) %>% 
                                                     mutate(`Fold&klength` = paste0(fold, ",", klen, sep="")) %>% 
                                                     arrange(desc(`CorrectId%`)) %>% 
                                                     dplyr::distinct() %>%
                                                     dplyr::relocate(any_of(c("rank", "Group", "fold", "klen", "CorrectId%", "n", "Fold&klength", "ordFamID"))))
    ordFamGroup_gen <- ordFamGroup_gen %>% select(!c("ordID", "famID"))
    
    sortGen <- percentCorrect_groupAvgGen[,"Group"] %>% rownames_to_column('sortVar')
    
    # percentCorrect_groupAvgGen <- rbind(ordFamGroup_gen, percentCorrect_groupAvgGen)
    percentCorrect_groupAvgGen <- percentCorrect_groupAvgGen %>%
      mutate(Group = factor(Group, levels = rev(as.character(unique(sortGen$Group)))))
    
    # Divide into 6 plots
    # Count up which groups fall into one of two groups for family based on the number of orders and 
    # number of individual families for each order
    alphaGroup_gen <- percentCorrect_groupAvgGen %>% 
      distinct() %>%
      group_by(ordFamID) %>%
      tally() 
    alphaGroup_gen$plotNum <- cumsum(alphaGroup_gen$n) %/% (ceiling(sum(alphaGroup_gen$n) / numPlots[3])) + 1
    alphaGroup_gen <- alphaGroup_gen %>% select(!n)
    alphaGroup_gen$plotNum <- ifelse(alphaGroup_gen$plotNum == max(alphaGroup_gen$plotNum), alphaGroup_gen$plotNum - 1, alphaGroup_gen$plotNum)
    
    # Merge with percentCorrect_groupAvgFam
    percentCorrect_groupAvgGen <- merge(percentCorrect_groupAvgGen, alphaGroup_gen, by="ordFamID")
    
    # Genus Heatmaps
    genPlots <- list()
    foreach(i=1:length(unique(percentCorrect_groupAvgGen$plotNum))) %do% {
      genPlots[[i]] <- ggplot(as.data.frame(percentCorrect_groupAvgGen[which(percentCorrect_groupAvgGen$plotNum == i),]), aes(`klen`, Group)) + 
        geom_tile(aes(fill=`CorrectId%`)) +
        geom_text(colour="darkorchid2", size=6, aes(label=n), fontface = "bold") + 
        scale_fill_viridis(discrete=FALSE) +
        theme_dark(base_size = 19) +
        # theme(panel.grid.major = element_blank(), 
        #      legend.position = "none") +
        facet_col(vars(ordFamID), scales = "free_y", space = "free") + 
        ylab("Genus") +
        xlab("klength")
    }
    
    # Save results plots
    foreach(i=1:length(unique(percentCorrect_groupAvgGen$plotNum))) %do% {
      ggsave(paste0(dirP, "/gen_percent_correct_", i, ".png", sep=""), genPlots[[i]], width = 10, height = 15, dpi = 500, bg="white")
    }
    # Add filepath to modelResults
    modelResults$PlotFilepaths$`Genus - HeatMap` <- unlist(foreach(i=1:length(unique(percentCorrect_groupAvgGen$plotNum))) %do% paste0(dirP, "/gen_percent_correct_", i, ".png", sep=""))
    
    # Genus ROC plot
    # TPR vs FPR
    confPlotGen <- list()
    foreach(i=1:length(klen_vec)) %do% { confPlotGen[[i]] <- roc(data = confidenceVal_group[[2]] %>%
                                                                   dplyr::filter(klen == klen_vec[i]) %>%
                                                                   mutate(truth = as.factor(ifelse(assignedTaxa == Group, 1, 0))) %>%
                                                                   select(!confThreshold) %>%
                                                                   distinct(), x = probScore, class = truth,
                                                                 pos_class = 1, neg_class = 0, direction = ">=") 
    
    confPlotGen[[i]] <- confPlotGen[[i]] %>% mutate(klenauc = paste0(klen_vec[i], ", AUC=", round(cutpointr::auc(confPlotGen[[i]]), 2), sep="")) }
    confPlotGen <- rbindlist(confPlotGen)
    
    nGroupsGen <- length(unique(confidenceVal_group[[2]]$Group))
    nSeqsGen <- length(unique(confidenceVal_group[[2]]$recordID))
    
    ggGen <- ggplot(data=confPlotGen, aes(x=fpr, y = tpr, col = klenauc)) + 
      geom_line(size=1.2, alpha=0.8) + 
      scale_colour_viridis_d(name = "klength & AUC") + 
      xlab("Average FPR") + ylab("Average TPR")  +
      ggtitle(paste0(classRank, " - Genus, numFolds = ", numFolds, ", nGroups = ", nGroupsGen, ", nSeqs = ", nSeqsGen, sep="")) +
      theme(text = element_text(size=16)) +
      annotate("segment", x = 0, xend = 1, y = 0, yend = 1, color="black", linetype="dashed") + 
      theme(panel.spacing = unit(1, "lines"))
    # ggsave(paste0(dirP, "/gen_roc.png"), plot=ggGen, width = 12, height = 8, dpi = 500, bg="white")
    
    # MCC scores
    confPlotGen$mcc <- ((confPlotGen$tn*confPlotGen$tp) - (confPlotGen$fp*confPlotGen$fn)) / sqrt((confPlotGen$tn+confPlotGen$fn) * (confPlotGen$fp+confPlotGen$tp) * (confPlotGen$tn+confPlotGen$fp) * (confPlotGen$fn+confPlotGen$tp))
    confPlotGen$klen <- substr(confPlotGen$klenauc, 1, 1)
    confPlotGen$mcc <- ifelse(is.nan(confPlotGen$mcc) & is.infinite(confPlotGen$x.sorted), 0, confPlotGen$mcc)
    confPlotGen$x.sorted <- ifelse(is.infinite(confPlotGen$x.sorted), 1, confPlotGen$x.sorted)
    
    ggGen2 <- ggplot(data=confPlotGen, aes(x=x.sorted, y = mcc, col = klen)) + 
      geom_line(size=1.2, alpha=0.8) + 
      scale_colour_viridis_d(name = "klength") + 
      xlab("Confidence Threshold") + ylab("Average MCC")  +
      # ggtitle(paste0(classRank, " - Genus, numFolds = ", numFolds, ", nGroups = ", nGroupsGen, ", nSeqs = ", nSeqsGen, sep="")) +
      theme(text = element_text(size=16)) +
      theme(panel.spacing = unit(1, "lines"))
    # ggsave(paste0(dirP, "/gen_mcc.png"), plot=ggGen2, width = 12, height = 8, dpi = 500, bg="white")
    
    # F1 scores
    confPlotGen$f1 <- (2*confPlotGen$tp) / ((2 * confPlotGen$tp) + confPlotGen$fn + confPlotGen$fp)
    
    ggGen3 <- ggplot(data=confPlotGen, aes(x=x.sorted, y = f1, col = klen)) + 
      geom_line(size=1.2, alpha=0.8) + 
      scale_colour_viridis_d(name = "klength") + 
      xlab("Confidence Threshold") + ylab("Average F1 score")  +
      # ggtitle(paste0(classRank, " - Genus, numFolds = ", numFolds, ", nGroups = ", nGroupsGen, ", nSeqs = ", nSeqsGen, sep="")) +
      theme(text = element_text(size=16)) +
      theme(panel.spacing = unit(1, "lines"))
    # ggsave(paste0(dirP, "/gen_f1.png"), plot=ggGen3, width = 12, height = 8, dpi = 500, bg="white")
    
    layoutGen <- rbind(c(1),
                       c(2),
                       c(3))
    ggCombGen <- grid.arrange(ggGen, ggGen2, 
                              ggGen3, layout_matrix = layoutGen)
    ggsave(paste0(dirP, "/gen_roc_mcc_f1.png"), plot=ggCombGen, width = 10, height = 15, dpi = 500, bg="white")
    
    # Add filepath to modelResults
    modelResults$PlotFilepaths$`Genus - ROC,MCC,F1` <- paste0(dirP, "/gen_roc_mcc_f1.png")
  }
  
  # If running species level or all levels
  if("Species" %in% subclassRank){
    
    # Find orders for each family listed in the results for sorting
    dfTaxaFamGen <- dfTaxa[,c("famID", "genID", "spID")] %>% distinct()
    colnames(dfTaxaFamGen)[3] <- "Group"
    
    famGenGroup_sp <- merge(percentCorrect_group[[4]], dfTaxaFamGen, by="Group") %>%
      mutate(famGenID = paste0(famID, ", ", genID, sep=""))
    
    # Finding average percent correct ids across all folds and then sorting by average percent correct ids per facet for plotting
    percentCorrect_groupAvgSp <- suppressWarnings(famGenGroup_sp %>% 
                                                    group_by(Group, klen) %>% 
                                                    dplyr::summarise(rank=rank, `CorrectId%` = mean(`CorrectId%`), n=sum(n), fold="Avg", famGenID=famGenID) %>% 
                                                    mutate(`Fold&klength` = paste0(fold, ",", klen, sep="")) %>% 
                                                    arrange(desc(`CorrectId%`)) %>% 
                                                    dplyr::distinct() %>%
                                                    dplyr::relocate(any_of(c("rank", "Group", "fold", "klen", "CorrectId%", "n", "Fold&klength", "famGenID"))))
    famGenGroup_sp <- famGenGroup_sp %>% select(!c("famID", "genID"))
    
    sortSp <- percentCorrect_groupAvgSp[,"Group"] %>% rownames_to_column('sortVar')
    
    # percentCorrect_groupAvgSp <- rbind(famGenGroup_sp, percentCorrect_groupAvgSp)
    percentCorrect_groupAvgSp <- percentCorrect_groupAvgSp %>%
      mutate(Group = factor(Group, levels = rev(as.character(unique(sortSp$Group)))))
    
    # Divide into 12 plots
    # Count up which groups fall into one of two groups for species based on the number of families and 
    # number of individual genuses for each family
    alphaGroup_sp <- percentCorrect_groupAvgSp %>% 
      distinct() %>%
      group_by(famGenID) %>%
      tally() 
    alphaGroup_sp$plotNum <- cumsum(alphaGroup_sp$n) %/% (ceiling(sum(alphaGroup_sp$n) / numPlots[4])) + 1
    alphaGroup_sp <- alphaGroup_sp %>% select(!n)
    alphaGroup_sp$plotNum <- ifelse(alphaGroup_sp$plotNum == max(alphaGroup_sp$plotNum), alphaGroup_sp$plotNum - 1, alphaGroup_sp$plotNum)
    
    # Merge with percentCorrect_groupAvgFam
    percentCorrect_groupAvgSp <- merge(percentCorrect_groupAvgSp, alphaGroup_sp, by="famGenID")
    
    # Species Heatmap plot
    spPlots <- list()
    foreach(i=1:length(unique(percentCorrect_groupAvgSp$plotNum))) %do% {
      spPlots[[i]] <- ggplot(as.data.frame(percentCorrect_groupAvgSp[which(percentCorrect_groupAvgSp$plotNum == i),]), aes(`klen`, Group)) + 
        geom_tile(aes(fill=`CorrectId%`)) +
        geom_text(colour="darkorchid2", size=6, aes(label=n), fontface = "bold") + 
        scale_fill_viridis(discrete=FALSE) +
        theme_dark(base_size = 19) +
        # theme(panel.grid.major = element_blank(), 
        #      legend.position = "none") +
        facet_col(vars(famGenID), scales = "free_y", space = "free") + 
        ylab("Species") +
        xlab("klength")
    }
    
    # Save results plots
    foreach(i=1:length(unique(percentCorrect_groupAvgSp$plotNum))) %do% {
      ggsave(paste0(dirP, "/sp_percent_correct_", i, ".png", sep=""), spPlots[[i]], width = 10, height = 15, dpi = 500, bg="white")
    }
    # Add filepath to modelResults
    modelResults$PlotFilepaths$`Species - HeatMap` <- unlist(foreach(i=1:length(unique(percentCorrect_groupAvgSp$plotNum))) %do% paste0(dirP, "/sp_percent_correct_", i, ".png", sep=""))
    
    # Species ROC plot
    # TPR vs FPR
    confPlotSp <- list()
    foreach(i=1:length(klen_vec)) %do% { confPlotSp[[i]] <- roc(data = confidenceVal_group[[4]] %>%
                                                                  dplyr::filter(klen == klen_vec[i]) %>%
                                                                  mutate(truth = as.factor(ifelse(assignedTaxa == Group, 1, 0))) %>%
                                                                  select(!confThreshold) %>%
                                                                  distinct(), x = probScore, class = truth,
                                                                pos_class = 1, neg_class = 0, direction = ">=") 
    
    confPlotSp[[i]] <- confPlotSp[[i]] %>% mutate(klenauc = paste0(klen_vec[i], ", AUC=", round(cutpointr::auc(confPlotSp[[i]]), 2), sep="")) }
    confPlotSp <- rbindlist(confPlotSp)
    
    nGroupsSp <- length(unique(confidenceVal_group[[4]]$Group))
    nSeqsSp <- length(unique(confidenceVal_group[[4]]$recordID))
    
    ggSp <- ggplot(data=confPlotSp, aes(x=fpr, y = tpr, col = klenauc)) + 
      geom_line(size=1.2, alpha=0.8) + 
      scale_colour_viridis_d(name = "klength & AUC") + 
      xlab("Average FPR") + ylab("Average TPR")  +
      ggtitle(paste0(classRank, " - Species, numFolds = ", numFolds, ", nGroups = ", nGroupsSp, ", nSeqs = ", nSeqsSp, sep="")) +
      theme(text = element_text(size=16)) +
      annotate("segment", x = 0, xend = 1, y = 0, yend = 1, color="black", linetype="dashed") + 
      theme(panel.spacing = unit(1, "lines"))
    # ggsave(paste0(dirP, "/sp_roc.png"), plot=ggSp, width = 12, height = 8, dpi = 500, bg="white")
    
    # MCC scores
    confPlotSp$mcc <- ((confPlotSp$tn*confPlotSp$tp) - (confPlotSp$fp*confPlotSp$fn)) / sqrt((confPlotSp$tn+confPlotSp$fn) * (confPlotSp$fp+confPlotSp$tp) * (confPlotSp$tn+confPlotSp$fp) * (confPlotSp$fn+confPlotSp$tp))
    confPlotSp$klen <- substr(confPlotSp$klenauc, 1, 1)
    confPlotSp$mcc <- ifelse(is.nan(confPlotSp$mcc) & is.infinite(confPlotSp$x.sorted), 0, confPlotSp$mcc)
    confPlotSp$x.sorted <- ifelse(is.infinite(confPlotSp$x.sorted), 1, confPlotSp$x.sorted)
    
    ggSp2 <- ggplot(data=confPlotSp, aes(x=x.sorted, y = mcc, col = klen)) + 
      geom_line(size=1.2, alpha=0.8) + 
      scale_colour_viridis_d(name = "klength") + 
      xlab("Confidence Threshold") + ylab("Average MCC")  +
      # ggtitle(paste0(classRank, " - Species, numFolds = ", numFolds, ", nGroups = ", nGroupsSp, ", nSeqs = ", nSeqsSp, sep="")) +
      theme(text = element_text(size=16)) +
      theme(panel.spacing = unit(1, "lines"))
    # ggsave(paste0(dirP, "/sp_mcc.png"), plot=ggSp2, width = 12, height = 8, dpi = 500, bg="white")
    
    # F1 scores
    confPlotSp$f1 <- (2*confPlotSp$tp) / ((2 * confPlotSp$tp) + confPlotSp$fn + confPlotSp$fp)
    
    ggSp3 <- ggplot(data=confPlotSp, aes(x=x.sorted, y = f1, col = klen)) + 
      geom_line(size=1.2, alpha=0.8) + 
      scale_colour_viridis_d(name = "klength") + 
      xlab("Confidence Threshold") + ylab("Average F1 score")  +
      # ggtitle(paste0(classRank, " - Species, numFolds = ", numFolds, ", nGroups = ", nGroupsSp, ", nSeqs = ", nSeqsSp, sep="")) +
      theme(text = element_text(size=16)) +
      theme(panel.spacing = unit(1, "lines"))
    # ggsave(paste0(dirP, "/sp_f1.png"), plot=ggSp3, width = 12, height = 8, dpi = 500, bg="white")
    
    layoutSp <- rbind(c(1),
                      c(2),
                      c(3))
    ggCombSp <- grid.arrange(ggSp, ggSp2, 
                             ggSp3, layout_matrix = layoutSp)
    ggsave(paste0(dirP, "/sp_roc_mcc_f1.png"), plot=ggCombSp, width = 10, height = 15, dpi = 500, bg="white")
    
    # Add filepath to modelResults
    modelResults$PlotFilepaths$`Species - ROC,MCC,F1` <- paste0(dirP, "/sp_roc_mcc_f1.png")
  }
  
  # Plot paneling by Order if all ranks are being run
  if(length(subclassRank)==4){
    # Find the unique Orders
    uniqOrd <- as.character(unique(percentCorrect_group[[3]]$Group))
    
    # Average the percent correct values across folds
    percentCorrect_OrdPan <- suppressWarnings(percentCorrect_group[[3]] %>%
                                                mutate(ordID = Group) %>%
                                                dplyr::relocate(any_of(c("rank", "Group", "fold", "klen", "CorrectId%", "n", "Fold&klength", "ordID"))) %>% 
                                                group_by(Group, klen) %>% 
                                                dplyr::summarise(rank=rank, `CorrectId%` = mean(`CorrectId%`), n=sum(n), fold="Avg", ordID=ordID) %>% 
                                                dplyr::distinct() %>%
                                                arrange(desc(`CorrectId%`)) %>%
                                                mutate(class = classRank))
    
    # Average the percent correct values across folds
    percentCorrect_FamPan <- suppressWarnings(percentCorrect_group[[1]] %>%
                                                mutate(famID = Group) %>%
                                                dplyr::relocate(any_of(c("rank", "Group", "fold", "klen", "CorrectId%", "n", "Fold&klength", "famID"))) %>% 
                                                group_by(Group, klen) %>% 
                                                dplyr::summarise(rank=rank, `CorrectId%` = mean(`CorrectId%`), n=sum(n), fold="Avg", famID=famID) %>% 
                                                dplyr::distinct()) 
    # Attach order & genus & species to family level results dataframe
    dfTaxaOrdFamPanel <- dfTaxa[,c("ordID", "famID")] %>% 
      distinct()
    percentCorrect_FamPan <- suppressWarnings(merge(percentCorrect_FamPan, dfTaxaOrdFamPanel, by="famID") %>%
                                                arrange(desc(`CorrectId%`)))
    # Sort according to correct ID and group
    sortFam <- as.data.frame(percentCorrect_FamPan[,"Group"]) %>% rownames_to_column('sortVar')
    colnames(sortFam)[2] <- "Group"
    percentCorrect_FamPan <- percentCorrect_FamPan %>%
      mutate(Group = factor(Group, levels = rev(as.character(unique(sortFam$Group)))))
    
    # Average percent correct values across folds
    percentCorrect_GenPan <- suppressWarnings(percentCorrect_group[[2]] %>%
                                                mutate(genID = Group) %>%
                                                dplyr::relocate(any_of(c("rank", "Group", "fold", "klen", "CorrectId%", "n", "Fold&klength", "genID"))) %>% 
                                                group_by(Group, klen) %>% 
                                                dplyr::summarise(rank=rank, `CorrectId%` = mean(`CorrectId%`), n=sum(n), fold="Avg", genID=genID) %>% 
                                                dplyr::distinct())
    # Attach order & family to genus level results dataframe 
    dfTaxaOrdFamPanel <- dfTaxa[,c("ordID", "famID","genID")] %>% 
      distinct()
    percentCorrect_GenPan <- suppressWarnings(merge(percentCorrect_GenPan, dfTaxaOrdFamPanel, by="genID") %>%
                                                arrange(desc(`CorrectId%`)))
    # Sort according to correct ID and group
    sortGen <- as.data.frame(percentCorrect_GenPan[,"Group"]) %>% rownames_to_column('sortVar')
    colnames(sortGen)[2] <- "Group"
    percentCorrect_GenPan <- percentCorrect_GenPan %>%
      mutate(Group = factor(Group, levels = rev(as.character(unique(sortGen$Group)))))
    
    # Attach order & family & genus to species level results dataframe and average percent correct values across folds
    percentCorrect_SpPan <- suppressWarnings(percentCorrect_group[[4]] %>%
                                               mutate(spID = Group) %>%
                                               dplyr::relocate(any_of(c("rank", "Group", "fold", "klen", "CorrectId%", "n", "Fold&klength", "spID"))) %>% 
                                               group_by(Group, klen) %>% 
                                               dplyr::summarise(rank=rank, `CorrectId%` = mean(`CorrectId%`), n=sum(n), fold="Avg", spID=spID) %>% 
                                               dplyr::distinct())
    # Attach order & family & species to genus level results dataframe 
    dfTaxaOrdFamPanel <- dfTaxa[,c("ordID", "famID","genID","spID")] %>% 
      distinct()
    percentCorrect_SpPan <- suppressWarnings(merge(percentCorrect_SpPan, dfTaxaOrdFamPanel, by="spID") %>%
                                               arrange(desc(`CorrectId%`)))
    # Sort according to correct ID and group
    sortSp <- as.data.frame(percentCorrect_SpPan[,"Group"]) %>% rownames_to_column('sortVar')
    colnames(sortSp)[2] <- "Group"
    percentCorrect_SpPan <- percentCorrect_SpPan %>%
      mutate(Group = factor(Group, levels = rev(as.character(unique(sortSp$Group)))))
    
    # Empty lists
    ordPanel <- list()
    famPanel <- list()
    genPanel <- list()
    spPanel <- list()
    
    # Subset each percent correct plot by order and extract out the groups present in that order
    foreach(i=1:length(uniqOrd)) %do% {
      ordPanel[[i]] <- ggplot(as.data.frame(percentCorrect_OrdPan %>% filter(ordID == uniqOrd[i])), aes(`klen`, Group)) + 
        geom_tile(aes(fill=`CorrectId%`)) +
        geom_text(colour="darkorchid2", size=6, aes(label=n), fontface = "bold") + 
        scale_fill_viridis(discrete=FALSE) +
        theme_dark(base_size = 19) +
        ggtitle("Order") +
        theme(panel.grid.major = element_blank(), 
              legend.position = "none") +
        facet_col(vars(class), scales = "free_y", space = "free") + 
        ylab("") +
        xlab("")
      # Depending on how large the species dataset is, create multiple plots for species if greater than 100 rows
      if(nrow(percentCorrect_SpPan %>% filter(ordID == uniqOrd[i]))<60){
        famPanel[[i]] <- ggplot(as.data.frame(percentCorrect_FamPan %>% filter(ordID == uniqOrd[i])), aes(`klen`, Group)) + 
          geom_tile(aes(fill=`CorrectId%`)) +
          geom_text(colour="darkorchid2", size=6, aes(label=n), fontface = "bold") + 
          scale_fill_viridis(discrete=FALSE) +
          theme_dark(base_size = 19) +
          ggtitle("Family") +
          theme(panel.grid.major = element_blank(), 
                legend.position = "none") +
          facet_col(vars(ordID), scales = "free_y", space = "free") + 
          ylab("") +
          xlab("")
        genPanel[[i]] <- ggplot(as.data.frame(percentCorrect_GenPan %>% filter(ordID == uniqOrd[i])), aes(`klen`, Group)) + 
          geom_tile(aes(fill=`CorrectId%`)) +
          geom_text(colour="darkorchid2", size=6, aes(label=n), fontface = "bold") + 
          scale_fill_viridis(discrete=FALSE) +
          theme_dark(base_size = 19) +
          ggtitle("Genus") +
          theme(panel.grid.major = element_blank(), 
                legend.position = "none") +
          facet_col(vars(famID), scales = "free_y", space = "free") + 
          ylab("") +
          xlab("")
        spPanel[[i]] <- ggplot(as.data.frame(percentCorrect_SpPan %>% filter(ordID == uniqOrd[i])), aes(`klen`, Group)) + 
          geom_tile(aes(fill=`CorrectId%`)) +
          geom_text(colour="darkorchid2", size=6, aes(label=n), fontface = "bold") + 
          scale_fill_viridis(discrete=FALSE) +
          theme_dark(base_size = 19) +
          ggtitle("Species") +
          theme(panel.grid.major = element_blank(), 
                legend.position = "none") +
          facet_col(vars(genID), scales = "free_y", space = "free") + 
          ylab("") +
          xlab("klength")
      } else if(nrow(percentCorrect_SpPan %>% filter(ordID == uniqOrd[i]))>= 60 & nrow(percentCorrect_SpPan %>% filter(ordID == uniqOrd[i]))<100) {
        famPanel[[i]] <- ggplot(as.data.frame(percentCorrect_FamPan %>% filter(ordID == uniqOrd[i])), aes(`klen`, Group)) + 
          geom_tile(aes(fill=`CorrectId%`)) +
          geom_text(colour="darkorchid2", size=6, aes(label=n), fontface = "bold") + 
          scale_fill_viridis(discrete=FALSE) +
          theme_dark(base_size = 19) +
          ggtitle("Family") +
          theme(panel.grid.major = element_blank(), 
                legend.position = "none") +
          facet_col(vars(ordID), scales = "free_y", space = "free") + 
          ylab("") +
          xlab("")
        genPanel[[i]] <- ggplot(as.data.frame(percentCorrect_GenPan %>% filter(ordID == uniqOrd[i])), aes(`klen`, Group)) + 
          geom_tile(aes(fill=`CorrectId%`)) +
          geom_text(colour="darkorchid2", size=6, aes(label=n), fontface = "bold") + 
          scale_fill_viridis(discrete=FALSE) +
          theme_dark(base_size = 19) +
          ggtitle("Genus") +
          theme(panel.grid.major = element_blank(), 
                legend.position = "none") +
          facet_col(vars(famID), scales = "free_y", space = "free") + 
          ylab("") +
          xlab("klength")
        spPanel[[i]] <- ggplot(as.data.frame(percentCorrect_SpPan %>% filter(ordID == uniqOrd[i])), aes(`klen`, Group)) + 
          geom_tile(aes(fill=`CorrectId%`)) +
          geom_text(colour="darkorchid2", size=6, aes(label=n), fontface = "bold") + 
          scale_fill_viridis(discrete=FALSE) +
          theme_dark(base_size = 19) +
          ggtitle("Species") +
          theme(panel.grid.major = element_blank(), 
                legend.position = "none") +
          facet_col(vars(genID), scales = "free_y", space = "free") + 
          ylab("") +
          xlab("klength")
      } else if(nrow(percentCorrect_SpPan %>% filter(ordID == uniqOrd[i]))>=100 & nrow(percentCorrect_SpPan %>% filter(ordID == uniqOrd[i]))<=400) {
        famPanel[[i]] <- ggplot(as.data.frame(percentCorrect_FamPan %>% filter(ordID == uniqOrd[i])), aes(`klen`, Group)) + 
          geom_tile(aes(fill=`CorrectId%`)) +
          geom_text(colour="darkorchid2", size=6, aes(label=n), fontface = "bold") + 
          scale_fill_viridis(discrete=FALSE) +
          theme_dark(base_size = 17) +
          ggtitle("Family") +
          theme(panel.grid.major = element_blank(), 
                legend.position = "none") +
          facet_col(vars(ordID), scales = "free_y", space = "free") + 
          ylab("") +
          xlab("klength")
        genPanel[[i]] <- ggplot(as.data.frame(percentCorrect_GenPan %>% filter(ordID == uniqOrd[i])), aes(`klen`, Group)) + 
          geom_tile(aes(fill=`CorrectId%`)) +
          geom_text(colour="darkorchid2", size=6, aes(label=n), fontface = "bold") + 
          scale_fill_viridis(discrete=FALSE) +
          theme_dark(base_size = 17) +
          ggtitle("Genus") +
          theme(panel.grid.major = element_blank(), 
                legend.position = "none") +
          facet_col(vars(famID), scales = "free_y", space = "free") + 
          ylab("") +
          xlab("klength")
        percentCorrect_SpPan1 <- percentCorrect_SpPan %>% filter(ordID == uniqOrd[i])
        # Divide into four plots
        # Count up which groups fall into one of four groups for species based on the number of genuses
        # number of individual families for each order
        alphaGroup <- percentCorrect_SpPan1 %>%
          distinct() %>%
          group_by(genID) %>%
          tally() 
        alphaGroup$plotNum <- cumsum(alphaGroup$n) %/% floor(sum(alphaGroup$n) / 3) + 1
        alphaGroup <- alphaGroup %>% select(!n)
        alphaGroup$plotNum <- ifelse(alphaGroup$plotNum == max(alphaGroup$plotNum), alphaGroup$plotNum - 1, alphaGroup$plotNum)
        percentCorrect_SpPan1 <- suppressWarnings(merge(percentCorrect_SpPan1, alphaGroup, by="genID") %>%
                                                    arrange(desc(`CorrectId%`)))
        
        # spPanel gets broken down into sublists 1:4
        spPanel[[i]] <- list()
        foreach(j=1:length(unique(percentCorrect_SpPan1$plotNum))) %do% {
          if(j==1 | j ==3){
            spPanel[[i]][[j]] <- ggplot(as.data.frame(percentCorrect_SpPan1 %>% filter(plotNum == j)), aes(`klen`, Group)) + 
              geom_tile(aes(fill=`CorrectId%`)) +
              geom_text(colour="darkorchid2", size=6, aes(label=n), fontface = "bold") + 
              scale_fill_viridis(discrete=FALSE) +
              theme_dark(base_size = 19) +
              ggtitle("Species") +
              theme(panel.grid.major = element_blank(), 
                    legend.position = "none") +
              facet_col(vars(genID), scales = "free_y", space = "free") + 
              ylab("") +
              xlab("klength")
          } else {
            spPanel[[i]][[j]] <- ggplot(as.data.frame(percentCorrect_SpPan1 %>% filter(plotNum == j)), aes(`klen`, Group)) + 
              geom_tile(aes(fill=`CorrectId%`)) +
              geom_text(colour="darkorchid2", size=6, aes(label=n), fontface = "bold") + 
              scale_fill_viridis(discrete=FALSE) +
              theme_dark(base_size = 19) +
              ggtitle("") +
              theme(panel.grid.major = element_blank(), 
                    legend.position = "none") +
              facet_col(vars(genID), scales = "free_y", space = "free") + 
              ylab("") +
              xlab("klength")
          }
        }
      } else {
        famPanel[[i]] <- ggplot(as.data.frame(percentCorrect_FamPan %>% filter(ordID == uniqOrd[i])), aes(`klen`, Group)) + 
          geom_tile(aes(fill=`CorrectId%`)) +
          geom_text(colour="darkorchid2", size=6, aes(label=n), fontface = "bold") + 
          scale_fill_viridis(discrete=FALSE) +
          theme_dark(base_size = 17) +
          ggtitle("Family") +
          theme(panel.grid.major = element_blank(), 
                legend.position = "none") +
          facet_col(vars(ordID), scales = "free_y", space = "free") + 
          ylab("") +
          xlab("klength")
        genPanel[[i]] <- ggplot(as.data.frame(percentCorrect_GenPan %>% filter(ordID == uniqOrd[i])), aes(`klen`, Group)) + 
          geom_tile(aes(fill=`CorrectId%`)) +
          geom_text(colour="darkorchid2", size=6, aes(label=n), fontface = "bold") + 
          scale_fill_viridis(discrete=FALSE) +
          theme_dark(base_size = 17) +
          ggtitle("Genus") +
          theme(panel.grid.major = element_blank(), 
                legend.position = "none") +
          facet_col(vars(famID), scales = "free_y", space = "free") + 
          ylab("") +
          xlab("klength")
        percentCorrect_SpPan1 <- percentCorrect_SpPan %>% filter(ordID == uniqOrd[i])
        # Divide into nine plots
        # Count up which groups fall into one of nine groups for species based on the number of genuses
        alphaGroup <- percentCorrect_SpPan1 %>%
          distinct() %>%
          group_by(genID) %>%
          tally() 
        alphaGroup$plotNum <- cumsum(alphaGroup$n) %/% floor(sum(alphaGroup$n) / 6) + 1
        alphaGroup <- alphaGroup %>% select(!n)
        alphaGroup$plotNum <- ifelse(alphaGroup$plotNum == max(alphaGroup$plotNum), alphaGroup$plotNum - 1, alphaGroup$plotNum)
        percentCorrect_SpPan1 <- suppressWarnings(merge(percentCorrect_SpPan1, alphaGroup, by="genID") %>%
                                                    arrange(desc(`CorrectId%`)))
        
        # spPanel gets broken down into sublists 1:9
        spPanel[[i]] <- list()
        foreach(j=1:length(unique(percentCorrect_SpPan1$plotNum))) %do% {
          if(j==1 | j ==4){
            spPanel[[i]][[j]] <- ggplot(as.data.frame(percentCorrect_SpPan1 %>% filter(plotNum == j)), aes(`klen`, Group)) + 
              geom_tile(aes(fill=`CorrectId%`)) +
              geom_text(colour="darkorchid2", size=6, aes(label=n), fontface = "bold") + 
              scale_fill_viridis(discrete=FALSE) +
              theme_dark(base_size = 19) +
              ggtitle("Species") +
              theme(panel.grid.major = element_blank(), 
                    legend.position = "none") +
              facet_col(vars(genID), scales = "free_y", space = "free") + 
              ylab("") +
              xlab("klength")
          } else {
            spPanel[[i]][[j]] <- ggplot(as.data.frame(percentCorrect_SpPan1 %>% filter(plotNum == j)), aes(`klen`, Group)) + 
              geom_tile(aes(fill=`CorrectId%`)) +
              geom_text(colour="darkorchid2", size=6, aes(label=n), fontface = "bold") + 
              scale_fill_viridis(discrete=FALSE) +
              theme_dark(base_size = 19) +
              ggtitle("") +
              theme(panel.grid.major = element_blank(), 
                    legend.position = "none") +
              facet_col(vars(genID), scales = "free_y", space = "free") + 
              ylab("") +
              xlab("klength")
          }
        }
      }
    }
    
    # Make paneled plots based on the size of the dataset (in terms of number of rows)
    foreach(i=1:length(uniqOrd)) %do% {
      if(length(spPanel[[i]])==11){
        if(nrow(spPanel[[i]]$data)<15){
          g1 <- ggplotGrob(ordPanel[[i]])
          g2 <- ggplotGrob(famPanel[[i]])
          g3 <- ggplotGrob(genPanel[[i]])
          g4 <- ggplotGrob(spPanel[[i]])
          g <- rbind(g1, g2, g3, g4, size = "first")
          g$widths <- unit.pmax(g1$widths, g2$widths, g3$widths, g4$widths)
          ggsave(paste0(dirP, "/percent_correct_panel_", uniqOrd[i], ".png"), plot=g, width = 7, height = 15, dpi = 500, bg="white")
        } else if(nrow(spPanel[[i]]$data)>=15 & nrow(spPanel[[i]]$data)<=50){
          g1 <- ggplotGrob(ordPanel[[i]])
          g2 <- ggplotGrob(famPanel[[i]])
          g3 <- ggplotGrob(genPanel[[i]])
          g4 <- ggplotGrob(spPanel[[i]])
          g <- rbind(g1, g2, g3, g4, size = "first")
          g$widths <- unit.pmax(g1$widths, g2$widths, g3$widths, g4$widths)
          ggsave(paste0(dirP, "/percent_correct_panel_", uniqOrd[i], ".png"), plot=g, width = 7, height = 15, dpi = 500, bg="white")
        } else if(nrow(spPanel[[i]]$data)>50){
          g1 <- ggplotGrob(ordPanel[[i]])
          g2 <- ggplotGrob(famPanel[[i]])
          g3 <- ggplotGrob(genPanel[[i]])
          g4 <- ggplotGrob(spPanel[[i]])
          g <- rbind(g1, g2, g3, size = "first")
          fg1 <- gtable_frame(g, debug = TRUE)
          fg2 <- gtable_frame(g4, debug = TRUE)
          fg12 <- gtable_frame(gtable_cbind(fg1, fg2),
                               width = unit(2, "null"),
                               height = unit(1, "null"))
          ggsave(paste0(dirP, "/percent_correct_panel_", uniqOrd[i], ".png"), plot=fg12, width = 15, height = 20, dpi = 500, bg="white")
        }
      } else if(length(spPanel[[i]])==4){
        # Order, family, genus go on one page
        g1 <- ggplotGrob(ordPanel[[i]])
        g2 <- ggplotGrob(famPanel[[i]])
        g3 <- ggplotGrob(genPanel[[i]])
        g <- rbind(g1, g2, size = "first")
        fg1 <- gtable_frame(g, debug = TRUE)
        fg2 <- gtable_frame(g3, debug = TRUE)
        fg12 <- gtable_frame(gtable_cbind(fg1, fg2),
                             width = unit(3, "null"),
                             height = unit(2, "null"))
        ggsave(paste0(dirP, "/percent_correct_panel_", uniqOrd[i], "_1.png"), plot=fg12, width = 20, height = 22, dpi = 500, bg="white")
        
        # Placing species on two additional pages
        
        # First Page
        gs1 <- gtable_frame(ggplotGrob(spPanel[[i]][[1]]))
        gs2 <- gtable_frame(ggplotGrob(spPanel[[i]][[2]]))
        gsp1 <- gtable_frame(gtable_cbind(gs1, gs2, gs3, gs4),
                             width = unit(3, "null"),
                             height = unit(3, "null"))
        ggsave(paste0(dirP, "/percent_correct_panel_", uniqOrd[i], "_2.png"), plot=gsp1, width = 20, height = 22, dpi = 500, bg="white")
        
        # Second Page
        gs3 <- gtable_frame(ggplotGrob(spPanel[[i]][[3]]))
        gs4 <- gtable_frame(ggplotGrob(spPanel[[i]][[4]]))
        gsp2 <- gtable_frame(gtable_cbind(gs3, gs4),
                             width = unit(3, "null"),
                             height = unit(3, "null"))
        ggsave(paste0(dirP, "/percent_correct_panel_", uniqOrd[i], "_3.png"), plot=gsp2, width = 20, height = 22, dpi = 500, bg="white")
        
      } else if(length(spPanel[[i]])==6){
        # Order, family, genus go on one page
        g1 <- ggplotGrob(ordPanel[[i]])
        g2 <- ggplotGrob(famPanel[[i]])
        g3 <- ggplotGrob(genPanel[[i]])
        g <- rbind(g1, g2, size = "first")
        fg1 <- gtable_frame(g, debug = TRUE)
        fg2 <- gtable_frame(g3, debug = TRUE)
        fg12 <- gtable_frame(gtable_cbind(fg1, fg2),
                             width = unit(5, "null"),
                             height = unit(3, "null"))
        ggsave(paste0(dirP, "/percent_correct_panel_", uniqOrd[i], "_1.png"), plot=fg12, width = 23, height = 24, dpi = 500, bg="white")
        
        # Placing species on three additional pages
        
        # First species page
        gs1 <- gtable_frame(ggplotGrob(spPanel[[i]][[1]]))
        gs2 <- gtable_frame(ggplotGrob(spPanel[[i]][[2]]))
        gs3 <- gtable_frame(ggplotGrob(spPanel[[i]][[3]]))
        gsp1 <- gtable_frame(gtable_cbind(gs1, gs2, gs3),
                             width = unit(5, "null"),
                             height = unit(3, "null"))
        ggsave(paste0(dirP, "/percent_correct_panel_", uniqOrd[i], "_2.png"), plot=gsp1, width = 23, height = 24, dpi = 500, bg="white")
        
        # Second species page
        gs4 <- gtable_frame(ggplotGrob(spPanel[[i]][[4]]))
        gs5 <- gtable_frame(ggplotGrob(spPanel[[i]][[5]]))
        gs6 <- gtable_frame(ggplotGrob(spPanel[[i]][[6]]))
        gsp2 <- gtable_frame(gtable_cbind(gs4, gs5, gs6),
                             width = unit(5, "null"),
                             height = unit(3, "null"))
        ggsave(paste0(dirP, "/percent_correct_panel_", uniqOrd[i], "_3.png"), plot=gsp2, width = 23, height = 24, dpi = 500, bg="white")
        
        # Third species page
        #gs7 <- gtable_frame(ggplotGrob(spPanel[[i]][[7]]))
        #gs8 <- gtable_frame(ggplotGrob(spPanel[[i]][[8]]))
        # gs9 <- gtable_frame(ggplotGrob(spPanel[[i]][[9]]))
        #gsp3 <- gtable_frame(gtable_cbind(gs7, gs8),
        #width = unit(4, "null"),
        #height = unit(2, "null"))
        #ggsave(paste0(dirP, "/percent_correct_panel_", uniqOrd[i], "_4.png"), plot=gsp3, width = 20, height = 23, dpi = 500, bg="white")
      }
    }
  }
  return(modelResults)
}
