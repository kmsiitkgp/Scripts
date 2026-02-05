#source("/hpc/home/kailasamms/projects/RNASeq/RNASeq_DESeq2_Functions.R")
source("C:/Users/kailasamms/Documents/GitHub/R-Scripts/RNASeq_DESeq2_Functions.R")

parent_path <- "C:/Users/kailasamms/Box/Saravana@cedars/05. Bioinformatics/CRISPR/CRISPR_Lena/count_results/"
output_path <- "C:/Users/kailasamms/Box/Saravana@cedars/05. Bioinformatics/CRISPR/CRISPR_Lena/"

#for (suffix in c("bowtie", "bowtie2", "mageck", "kms")){
for (library in c("CDH12KO", "CDH12ACT", "ImmuneKO", "ImmuneACT")){
  
  # Read the count summary file and count files
  countsummary <- read.table(paste0(parent_path, "count_results_", library, "/", "mageck.total", ".countsummary.txt"), 
                             header=TRUE)
  
  count_raw <- read.table(paste0(parent_path, "count_results_", library, "/", "mageck.total", ".count.txt"), 
                             header=TRUE)
  
  # count_total1 <- read.table(paste0(parent_path, "count_results_", library, "/", "mageck.total", ".count_normalized.txt"), 
  #                            header=TRUE)
  count_total <- count_raw %>%
    dplyr::mutate(across(.cols=c(everything(), -c(sgRNA, Gene)), ~ (1000000*.)/sum(.)))
  
  # count_median <- read.table(paste0(parent_path, "count_results_", library, "/", "mageck.median", ".count_normalized.txt"), 
  #                            header=TRUE)
  
  count_median <- count_raw %>%
    dplyr::mutate(across(.cols=c(everything(), -c(sgRNA, Gene)), ~ ()))
  
  # Gini index
  p1 <- ggplot(data = countsummary, aes(x=Label, y=GiniIndex,fill=Label)) + 
    geom_col() +              
    theme_classic() +         #display with x and y axis lines and no gridlines
    my_theme +
    labs(x = "Sample", y = "Gini Index", title = stringr::str_wrap("Evenness of sgRNA reads", 30)) +
    coord_cartesian(clip = "off") +
    geom_text(aes(label=GiniIndex), y = 0, hjust = 0, angle = 90)
  
  # Missed sgRNAs
  p2 <- ggplot(data = countsummary, aes(x=Label, y=Zerocounts,fill=Label)) + 
    geom_col() +              
    theme_classic() +         #display with x and y axis lines and no gridlines
    my_theme +
    labs(x = "Sample", y = "Number of Missed sgRNAs", title = stringr::str_wrap("Zero Count sgRNAs", 30)) +
    coord_cartesian(clip = "off") +
    geom_text(aes(label=Zerocounts), y = 0, hjust = 0, angle = 90)
  
  # Read mapping
  data <- countsummary %>%
    dplyr::mutate(Percent_Mapped = round(100*Mapped/Reads, digits=2),
                  Percent_Unmapped = round(100*(Reads-Mapped)/Reads, digits=2)) %>%
    dplyr::select(Label, Percent_Mapped, Percent_Unmapped) %>%
    tidyr::pivot_longer(cols=!Label, names_to = "Mapping", values_to = "Percent")
  
  p3 <- ggplot(data = data, aes(x=Label, y=Percent,fill=Mapping)) + 
    geom_col(position="stack") +              
    theme_classic() +         #display with x and y axis lines and no gridlines
    my_theme +
    labs(x = "Sample", y = "Reads", title = stringr::str_wrap("Mapping Ratio", 30)) +
    coord_cartesian(clip = "off") +
    geom_text(aes(label=Percent), position = position_stack(vjust = .5), hjust = 0, angle = 90)
  
  count_raw_violin <- count_raw %>% 
    dplyr::select(everything(), -Gene) %>%
    tidyr::pivot_longer(cols = !sgRNA, names_to="Sample", values_to="Reads")
  
  p4 <- ggplot(data = count_raw_violin, aes(x=Sample, y=log2(Reads),fill=Sample)) +
    geom_violin() +        
    theme_classic() +
    my_theme +
    labs(x = "Sample", y = "log10(Reads)", title = stringr::str_wrap("Raw Reads Distribution", 30)) +
    coord_cartesian(ylim = c(0, 8), clip = "off")
  
  count_total_violin <- count_total %>% 
    dplyr::select(everything(), -Gene) %>%
    tidyr::pivot_longer(cols = !sgRNA, names_to="Sample", values_to="Reads")
  
  p5 <- ggplot(data = count_total_violin, aes(x=Sample, y=log10(Reads),fill=Sample)) +
    geom_violin() +        
    theme_classic() +
    my_theme +
    labs(x = "Sample", y = "log10(Reads)", title = stringr::str_wrap("Total Normalized Reads Distribution", 30)) +
    coord_cartesian(ylim = c(0, 8), clip = "off")
  
  count_median_violin <- count_median %>% 
    dplyr::select(everything(), -Gene) %>%
    tidyr::pivot_longer(cols = !sgRNA, names_to="Sample", values_to="Reads")
  
  p6 <- ggplot(data = count_median_violin, aes(x=Sample, y=log10(Reads),fill=Sample)) +
    geom_violin() +        
    theme_classic() +
    my_theme +
    labs(x = "Sample", y = "log10(Reads)", title = stringr::str_wrap("Median Normalized Reads Distribution", 30)) +
    coord_cartesian(ylim = c(0, 8), clip = "off")
  
  # Save all plots
  plot_grid(p1, p2, p3, p4, p5, p6,
            align = c("hv"),
            axis = c("tblr"),
            nrow = 6,
            ncol = 1,
            rel_widths = 1,
            rel_heights = 1,
            labels = NULL,
            label_size = 14,
            label_fontfamily = NULL,
            label_fontface = "bold",
            label_colour = NULL,
            label_x = 0,
            label_y = 1,
            hjust = -0.5,
            vjust = 1.5,
            scale = 1,
            greedy = TRUE,
            byrow = TRUE)
  
  ggsave(filename = paste0(library, "_QC.tiff"),
         plot = last_plot(),
         device = "tiff",
         path = output_path,
         width = 11,
         height = 22,
         units = c("in"),	 
         dpi = 300,
         limitsize = TRUE,
         bg = "white")
}


for (suffix in c("bowtie", "bowtie2", "mageck", "kms")){
  for (comparison in c("Day7vs0_2D", "Day14vs0_2D", "Day7vs0_3D", "Day14vs0_3D")){
    
    # Read the common essential genes
    common_essential <- read.xlsx("C:/Users/kailasamms/Box/Saravana@cedars/05. Bioinformatics/Depmap Gene Dependency Profile Summary.xlsx",
                                  sheet = "Common Essential Human Genes")
    common_essential <- common_essential$All.combined
    
    # Read the gene level data and remove common essential genes
    gdata <- read.table(paste0(parent_path, "DEG_results/", suffix, "/", comparison, ".gene_summary.txt"), 
                        header=TRUE) %>%
      dplyr::filter(!(id %in% common_essential)) %>%
      dplyr::select(id, neg.score, neg.lfc, neg.fdr, neg.goodsgrna) %>%
      dplyr::rename(RRAScore = neg.score,
                    LFC = neg.lfc, 
                    FDR = neg.fdr,
                    GoodsgRNA = neg.goodsgrna)
    
    # Read the sgRNA level data and remove common essential genes
    sdata <- read.table(paste0(parent_path, "DEG_results/", suffix, "/", comparison, ".sgrna_summary.txt"), 
                        header=TRUE) %>%
      dplyr::filter(!(Gene %in% common_essential)) %>%
      dplyr::select(sgrna, Gene, LFC, FDR)
    
    
    # Visualization of negative selections and positive selections 
    # Volcano plot
    output_path <- "C:/Users/kailasamms/Desktop/"
    Target <- str_split(string=comparison, pattern="vs")[[1]][1]
    Reference <- "Day 0"
    log2_cutoff <- 0.58
    padj_cutoff <- 0.05
    file_suffix <- paste0(comparison, "_", suffix)
    volcano_df <- gdata %>% dplyr::rename(log2FC = LFC, padj = FDR, SYMBOL = id)
    disp_genes <- "MAPK1"
    
    p1 <- plot_volcano(volcano_df, disp_genes, Target, Reference, log2_cutoff, padj_cutoff, file_suffix, output_path)
    
    # Rank all the genes based on their RRA scores and label genes in the rank plot
    # Rank plot
    gdata$Rank = base::rank(gdata$RRAScore)  # the smallest RRA score gets rank 1
    
    p2 <- ggplot(data = gdata, aes(x=Rank, y=-log10(RRAScore), fill=RRAScore)) +
      geom_point() +
      theme_classic() +         #display with x and y axis lines and no gridlines
      my_theme +
      labs(x = "Rank", y = expression("-log"[10]*"(RRA Score)"), title = stringr::str_wrap("RRA Ranks", 30)) +
      coord_cartesian(clip = "off") +
      geom_text_repel(data = gdata %>% dplyr::filter(id %in% disp_genes),
                      mapping = aes(label = id),
                      force = 0.5,
                      point.size = 1,
                      angle = 0,
                      #vjust = 0,
                      #hjust = 0,
                      #direction = "y",
                      box.padding = 1,  # increases line length somehow
                      point.padding = 0.1,
                      max.overlaps = Inf,
                      xlim = c(NA, NA),
                      ylim = c(-Inf,NA),
                      min.segment.length = 0.2,
                      #min.segment.length = dplyr::if_else(draw_line == "TRUE", 0.2, Inf),
                      #arrow = arrow(length = unit(0.015, "npc")),
                      position = position_quasirandom())
    
    # Visualize negative and positive selected genes separately
    # Dot plot
    # NOTE: We consider genes with LFC more than 2 Standard deviations as 
    # differentially selected
    cutoff_limit <- stats::sd(gdata$LFC)*2
    palette <- c("#808080", "#00BFC4", "#F8766D")
    names(palette) <- c("Not Selected", "Negatively Selected", "Positively Selected")
    
    gdata$RandomIndex = sample(1:nrow(gdata), nrow(gdata))
    gdata = gdata[order(-gdata$LFC), ]
    gdata <- gdata %>% 
      dplyr::mutate(cutoff = dplyr::case_when(abs(LFC) >= cutoff_limit & LFC > 0 ~ "Positively Selected",
                                              abs(LFC) >= cutoff_limit & LFC < 0 ~ "Negatively Selected",
                                              TRUE ~ "Not Selected"))
    
    gg <- gdata[gdata$LFC<0, ]
    
    p3 <- ggplot(data = gg, aes(x=RandomIndex, y=LFC, fill=cutoff)) +
      
      # If you want all points to be fixed size, use size=0.2 within geom_point()
      # If you want to adjust point size based on column in volcano_df, declare in global aes() within ggplot()
      ggplot2::geom_point(col="black", 
                          shape=21,
                          stroke=0.5,
                          position=position_jitter(h=0.01,w=0.01)) +
      
      # Define the theme of plot
      ggplot2::theme_classic() +
      
      # Define the color of the dots
      scale_fill_manual(values = palette) +
      
      # Define the axis, plot headings
      ggplot2::labs(x = "RandomIndex",
                    y = "LFC",
                    #fill = "Direction",
                    fill = "log2FC",
                    size = "-log10(padj)*log2FC",
                    title = stringr::str_wrap("Negatively selected Genes", 30)) +
      
      # Draw line to mark the cutoffs
      geom_hline(yintercept = -cutoff_limit, color = "black", linetype = "dotted", linewidth = 0.5) +
      
      # Define axis
      coord_cartesian(clip = "off") +
      
      # Adjust font size, style
      my_theme +
      
      # Add gene labels
      geom_text_repel(data = gg %>% dplyr::filter(id %in% disp_genes),
                      mapping = aes(label = id),
                      #size = 2,
                      force = 0.5,
                      point.size = 1,
                      angle = 0,
                      #vjust = 0,
                      #hjust = 0,
                      #direction = "y",
                      box.padding = 1,  # increases line length somehow
                      point.padding = 0.1,
                      max.overlaps = Inf,
                      xlim = c(NA, NA),
                      ylim = c(-Inf,NA),
                      min.segment.length = 0.2,
                      #min.segment.length = dplyr::if_else(draw_line == "TRUE", 0.2, Inf),
                      #arrow = arrow(length = unit(0.015, "npc")),
                      position = position_quasirandom())
 
    
    gg <- gdata[gdata$LFC>0, ]
    p4 <- ggplot(data = gg, aes(x=RandomIndex, y=LFC, fill=cutoff)) +
      
      # If you want all points to be fixed size, use size=0.2 within geom_point()
      # If you want to adjust point size based on column in volcano_df, declare in global aes() within ggplot()
      ggplot2::geom_point(col="black", 
                          shape=21,
                          stroke=0.5,
                          position=position_jitter(h=0.01,w=0.01)) +
      
      # Define the theme of plot
      ggplot2::theme_classic() +
      
      # Define the color of the dots
      scale_fill_manual(values = palette) +
      
      # Define the axis, plot headings
      ggplot2::labs(x = "RandomIndex",
                    y = "LFC",
                    #fill = "Direction",
                    fill = "log2FC",
                    size = "-log10(padj)*log2FC",
                    title = stringr::str_wrap("Positively selected Genes", 30)) +
      
      # Draw line to mark the cutoffs
      geom_hline(yintercept = cutoff_limit, color = "black", linetype = "dotted", linewidth = 0.5) +
      
      # Define axis
      coord_cartesian(clip = "off") +
      
      # Adjust font size, style
      my_theme +
      
      # Add gene labels
      geom_text_repel(data = gg %>% dplyr::filter(id %in% disp_genes),
                      mapping = aes(label = id),
                      #size = 2,
                      force = 0.5,
                      point.size = 1,
                      angle = 0,
                      #vjust = 0,
                      #hjust = 0,
                      #direction = "y",
                      box.padding = 1,  # increases line length somehow
                      point.padding = 0.1,
                      max.overlaps = Inf,
                      xlim = c(NA, NA),
                      ylim = c(-Inf,NA),
                      min.segment.length = 0.2,
                      #min.segment.length = dplyr::if_else(draw_line == "TRUE", 0.2, Inf),
                      #arrow = arrow(length = unit(0.015, "npc")),
                      position = position_quasirandom())
    
    
    # Save all plots
    plot_grid(p1,p2,p3,p4,
              align = c("hv"),
              axis = c("tblr"),
              nrow = 2,
              ncol = 2,
              rel_widths = 1,
              rel_heights = 1,
              labels = NULL,
              label_size = 14,
              label_fontfamily = NULL,
              label_fontface = "bold",
              label_colour = NULL,
              label_x = 0,
              label_y = 1,
              hjust = -0.5,
              vjust = 1.5,
              scale = 1,
              greedy = TRUE,
              byrow = TRUE)
    
    ggsave(filename = paste0(comparison, "_", suffix,".tiff"),
           plot = last_plot(),
           device = "tiff",
           path = output_path,
           width = 11,
           height = 11,
           units = c("in"),	 
           dpi = 300,
           limitsize = TRUE,
           bg = "white")
    
  }  
}


# Combine all reads together
