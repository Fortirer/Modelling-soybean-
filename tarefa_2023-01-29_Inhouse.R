library(consensusDE)
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene) 

#file_list <- list.files("/home/rstudio/inhouse",recursive = TRUE, pattern = "*bam$",full=TRUE)
file_list <- list.files("/srv/database/data/NGS_raw_data/transcriptome/rnaSeq_lafieco/soja/inhouse", recursive = TRUE, pattern = "*bam$", full = TRUE)

# download from pc to server
#scp -r  /srv/database/data/NGS_raw_data/transcriptome/rnaSeq_lafieco/soja/Gmax_275_Wm82.a2.v1.gene.gff3 janaina@143.107.247.188:/home/janaina/Downloads/

gff3file <- file.path("/home/janaina/Downloads/Gmax_275_Wm82.a2.v1.gene.gff3") # gene
txdb <- makeTxDbFromGFF(gff3file, format="gff3", circ_seqs=character())

#names <- substr(basename(file_list), start = 1, stop = 3)
#names_1 <- sub(pattern = "15[0-4]", replacement = "Elev", x = names)
#names_2 <- sub(pattern = "25[0-4]", replacement = "Elev.Temp", x = names_1)
#names_3 <- sub(pattern = "35[0-4]", replacement = "Amb", x = names_2)

#names_new <- sub(pattern = "45[0-4]", replacement = "Temp", x = names_3)
#sample_table <- data.frame("file" = basename(file_list),"group" = names_new)

sample_table <- data.frame("file" = basename(file_list),"group" =  c("Amb","Amb","Amb","Amb","Amb","Elev","Elev","Elev","Elev","Elev","Elev.Temp","Elev.Temp","Elev.Temp","Elev.Temp","Elev.Temp","Temp","Temp","Temp","Temp","Temp"))

bam_dir <- as.character(gsub(basename(file_list)[1], "", file_list[1]))

#make.names(summarized_dm3$group) ##### edita as strings para um formato aceitável pelo R 
#summarized_dm3$group <- make.names(summarized_dm3$group)

#nao usei setwd()
summarized_dm3 <- buildSummarized(sample_table = sample_table,
                                  bam_dir = bam_dir,
                                  tx_db = txdb,
                                  read_format = "paired",  # deu 20 warnings paired - colocar single
                                  output_log="/home/janaina/Downloads",                                  
                                  force_build = TRUE)
                              

#all_pairs_airway <- multi_de_pairs(summarized = summarized_dm3, paired = "unpaired" , ruv_correct = T) #### com todos os dados

set.seed(1234) # Para garantir que a sua subamostragem seja a mesma. Uma vez que vc deseje alterar algum parâmetro na all_pairs_airway(), a set.seed() garante que o cunjunto de subamostragem de genes seja o mesmo que na análise anterior.
summarized_dm3_filter <- sample(summarized_dm3, 1000) # faz uma subamostragem n=1000 do seu data.frame
all_pairs_airway <- multi_de_pairs(summarized = summarized_dm3_filter, paired = "unpaired" , ruv_correct = T)


# call multi_de_pairs()
all_pairs_airway_ruv <- multi_de_pairs(summarized = summarized_dm3_filter,
                                       paired = "unpaired",
                                       ruv_correct = TRUE)
                                       
# To view all the comparisons conducted:
names(all_pairs_airway$merged)

# to access data of a particular comparison           
head(all_pairs_airway$merged[["Elev-Amb"]])

# access the summarized experiment (now including the residuals under the "W_1" column)
all_pairs_airway_ruv$summarized@phenoData@data

#Grafico Boxplot
pdf(file = "/home/janaina/Downloads/BoxPlot_FilterPaired.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 7) # The height of the plot in inches

diag_plots(se_in = summarized_dm3_filter, boxplot = TRUE)
# Step 3: Run dev.off() to create the file!
dev.off()

#Grafico mapped reads
pdf(file = "/home/janaina/Downloads/MappedPlotPaired.pdf")   # The directory you want to save the file in
diag_plots(se_in = summarized_dm3_filter,
           name = "airway example data",
           mapped_reads = TRUE)
dev.off()
           
# Relative Log Expression
pdf(file = "/home/janaina/Downloads/RLEPlotPaired.pdf")   # The directory you want to save the file in
diag_plots(se_in = summarized_dm3_filter,
           name = "airway example data",
           rle = TRUE)
dev.off()

#PCA
pdf(file = "/home/janaina/Downloads/PCAPlotPaired.pdf")   # The directory you want to save the file in
diag_plots(se_in = summarized_dm3_filter,pca = TRUE)
dev.off()

############################# COMPARAÇÃO ELEVADO E AMBIENTE ####################################
#MAPlot
# 1. View all the comparisons conducted
names(all_pairs_airway$merged)

comparison <- all_pairs_airway$merged[["Elev-Amb"]]

# this will not work unless in a list and will stop, producing an error. E.g.
#diag_plots(merged_in = comparison,
#           name = "Elev-Amb",
#           ma = TRUE)

# Error message:
#merged_in is not a list. If you want to plot with one comparison only,
#put the single dataframe into a list as follows. my_list <- list("name"=
#merged_in)

# 3. Put into a new list as instructed by the error
comparison_list <- list("Elev-Amb" = comparison)

# this will not work unless the appropriate columns are labelled
# "ID", "AveExpr", and "Adj_PVal"

# 4. Relabel the columns for plotting
# inspecting the column names reveals that the "Adj_PVal" column needs to be specified.


# Here, we will relabel "edger_adj_p" with "Adj_PVal" to use this p-value, using
# the "gsub" command as follows (however, we could also use one of the others or
# the p_max column)

colnames(comparison_list[["Elev-Amb"]]) <- gsub("edger_adj_p", "Adj_PVal",
                                                 colnames(comparison_list[["Elev-Amb"]]))

# after label
colnames(comparison_list[["Elev-Amb"]])

# 5. Plot MA
pdf(file = "/home/janaina/Downloads/ElevAmbMA2Plot_Filter.pdf",   # The directory you want to save the file in
    width = 3, # The width of the plot in inches
    height = 3)
diag_plots(merged_in = comparison_list,
           name = "Elev-Amb",
           ma = TRUE)
dev.off()

#Volcano
pdf(file = "/home/janaina/Downloads/ElevAmbVolcanoPlot.pdf")
diag_plots(merged_in = comparison_list,
           name = "Elev-Amb",
           volcano = TRUE)
dev.off()

#p-value distribution
pdf(file = "/home/janaina/Downloads/ElevAmbPvaluePlot.pdf")
diag_plots(merged_in = comparison_list,
           name = "Elev-Amb",
           p_dist = TRUE)
dev.off()

############################# COMPARAÇÃO ELEVADO E TEMPERATURA ####################################
#MAPlot
# 1. View all the comparisons conducted
names(all_pairs_airway$merged)

comparison <- all_pairs_airway$merged[["Elev-Temp"]]

# this will not work unless in a list and will stop, producing an error. E.g.
diag_plots(merged_in = comparison,
           name = "Elev-Temp",
           ma = TRUE)

# Error message:
#merged_in is not a list. If you want to plot with one comparison only,
#put the single dataframe into a list as follows. my_list <- list("name"=
#merged_in)

# 3. Put into a new list as instructed by the error
comparison_list <- list("Elev-Temp" = comparison)

# this will not work unless the appropriate columns are labelled
# "ID", "AveExpr", and "Adj_PVal"

# 4. Relabel the columns for plotting
# inspecting the column names reveals that the "Adj_PVal" column needs to be specified.


# Here, we will relabel "edger_adj_p" with "Adj_PVal" to use this p-value, using
# the "gsub" command as follows (however, we could also use one of the others or
# the p_max column)

colnames(comparison_list[["Elev-Temp"]]) <- gsub("edger_adj_p", "Adj_PVal",
                                                 colnames(comparison_list[["Elev-Temp"]]))

# after label
colnames(comparison_list[["Elev-Temp"]])

# 5. Plot MA
pdf(file = "/home/janaina/Downloads/ElevTempMA2Plot_Filter.pdf",   # The directory you want to save the file in
    width = 1, # The width of the plot in inches
    height = 1)
diag_plots(merged_in = comparison_list,
           name = "Elev-Temp",
           ma = TRUE)
dev.off()

#Volcano
pdf(file = "/home/janaina/Downloads/EleTempVolcanoPlot.pdf")
diag_plots(merged_in = comparison_list,
           name = "Elev-Temp",
           volcano = TRUE)
dev.off()

#p-value distribution
pdf(file = "/home/janaina/Downloads/ElevTempPvaluePlot.pdf")
diag_plots(merged_in = comparison_list,
           name = "Elev-Temp",
           p_dist = TRUE)
dev.off()

############################# COMPARAÇÃO TEMPERATURA E AMBIENTE ####################################
#MAPlot
# 1. View all the comparisons conducted
names(all_pairs_airway$merged)

comparison <- all_pairs_airway$merged[["Temp-Amb"]]

# this will not work unless in a list and will stop, producing an error. E.g.
diag_plots(merged_in = comparison,
           name = "Temp-Amb",
           ma = TRUE)

# Error message:
#merged_in is not a list. If you want to plot with one comparison only,
#put the single dataframe into a list as follows. my_list <- list("name"=
#merged_in)

# 3. Put into a new list as instructed by the error
comparison_list <- list("Temp-Amb" = comparison)

# this will not work unless the appropriate columns are labelled
# "ID", "AveExpr", and "Adj_PVal"

# 4. Relabel the columns for plotting
# inspecting the column names reveals that the "Adj_PVal" column needs to be specified.


# Here, we will relabel "edger_adj_p" with "Adj_PVal" to use this p-value, using
# the "gsub" command as follows (however, we could also use one of the others or
# the p_max column)

colnames(comparison_list[["Temp-Amb"]]) <- gsub("edger_adj_p", "Adj_PVal",
                                                 colnames(comparison_list[["Temp-Amb"]]))

# after label
colnames(comparison_list[["Temp-Amb"]])

# 5. Plot MA
pdf(file = "/home/janaina/Downloads/TempAmbMA2Plot_Filter.pdf",   # The directory you want to save the file in
    width = 1, # The width of the plot in inches
    height = 1)
diag_plots(merged_in = comparison_list,
           name = "Temp-Amb",
           ma = TRUE)
dev.off()

#Volcano
pdf(file = "/home/janaina/Downloads/TempAmbVolcanoPlot.pdf")
diag_plots(merged_in = comparison_list,
           name = "Temp-Amb",
           volcano = TRUE)
dev.off()

#p-value distribution
pdf(file = "/home/janaina/Downloads/TempAmbPvaluePlot.pdf")
diag_plots(merged_in = comparison_list,
           name = "Temp-Amb",
           p_dist = TRUE)
dev.off()

############################# COMPARAÇÃO ELEVADO.TEMPERATURA - AMBIENTE ##################################
#MAPlot
# 1. View all the comparisons conducted
names(all_pairs_airway$merged)

comparison <- all_pairs_airway$merged[["Elev.Temp-Amb"]]

# this will not work unless in a list and will stop, producing an error. E.g.
diag_plots(merged_in = comparison,
           name = "Elev.Temp-Amb",
           ma = TRUE)

# Error message:
#merged_in is not a list. If you want to plot with one comparison only,
#put the single dataframe into a list as follows. my_list <- list("name"=
#merged_in)

# 3. Put into a new list as instructed by the error
comparison_list <- list("Elev.Temp-Amb" = comparison)

# this will not work unless the appropriate columns are labelled
# "ID", "AveExpr", and "Adj_PVal"

# 4. Relabel the columns for plotting
# inspecting the column names reveals that the "Adj_PVal" column needs to be specified.


# Here, we will relabel "edger_adj_p" with "Adj_PVal" to use this p-value, using
# the "gsub" command as follows (however, we could also use one of the others or
# the p_max column)

colnames(comparison_list[["Elev.Temp-Amb"]]) <- gsub("edger_adj_p", "Adj_PVal",
                                                 colnames(comparison_list[["Elev.Temp-Amb"]]))

# after label
colnames(comparison_list[["Elev.Temp-Amb"]])

# 5. Plot MA
pdf(file = "/home/janaina/Downloads/ElevTemp-AmbMA2Plot_Filter.pdf")   # The directory you want to save the file in
    #width = 10, # The width of the plot in inches
    #height = 10)
diag_plots(merged_in = comparison_list,
           name = "Elev.Temp-Amb",
           ma = TRUE)
dev.off()

#Volcano
pdf(file = "/home/janaina/Downloads/ElevTemp-AmbVolcanoPlot.pdf")
diag_plots(merged_in = comparison_list,
           name = "Elev.Temp-Amb",
           volcano = TRUE)
dev.off()

#p-value distribution
pdf(file = "/home/janaina/Downloads/ElevTemp-AmbPvaluePlot.pdf")
diag_plots(merged_in = comparison_list,
           name = "Elev.Temp-Amb",
           p_dist = TRUE)
dev.off()

############################# COMPARAÇÃO ELEVADO.TEMPERATURA-TEMPERATURA ##################################
#MAPlot
# 1. View all the comparisons conducted
names(all_pairs_airway$merged)

comparison <- all_pairs_airway$merged[["Elev.Temp-Temp"]]

# this will not work unless in a list and will stop, producing an error. E.g.
#diag_plots(merged_in = comparison,
 #          name = "Elev.Temp-Temp",
  #         ma = TRUE)

# Error message:
#merged_in is not a list. If you want to plot with one comparison only,
#put the single dataframe into a list as follows. my_list <- list("name"=
#merged_in)

# 3. Put into a new list as instructed by the error
comparison_list <- list("Elev.Temp-Temp" = comparison)

# this will not work unless the appropriate columns are labelled
# "ID", "AveExpr", and "Adj_PVal"

# 4. Relabel the columns for plotting
# inspecting the column names reveals that the "Adj_PVal" column needs to be specified.


# Here, we will relabel "edger_adj_p" with "Adj_PVal" to use this p-value, using
# the "gsub" command as follows (however, we could also use one of the others or
# the p_max column)

colnames(comparison_list[["Elev.Temp-Temp"]]) <- gsub("edger_adj_p", "Adj_PVal",
                                                 colnames(comparison_list[["Elev.Temp-Temp"]]))

# after label
colnames(comparison_list[["Elev.Temp-Temp"]])

# 5. Plot MA
pdf(file = "/home/janaina/Downloads/ElevTemp-TempMA2Plot_Filter.pdf",   # The directory you want to save the file in
   width = 1, # The width of the plot in inches
    height = 1)
diag_plots(merged_in = comparison_list,
           name = "Elev.Temp-Temp",
           ma = TRUE)
dev.off()

#Volcano
pdf(file = "/home/janaina/Downloads/ElevTemp-TempVolcanoPlot.pdf")
diag_plots(merged_in = comparison_list,
           name = "Elev.Temp-Temp",
           volcano = TRUE)
dev.off()

#p-value distribution
pdf(file = "/home/janaina/Downloads/ElevTemp-TempPvaluePlot.pdf")
diag_plots(merged_in = comparison_list,
           name = "Elev.Temp-Temp",
           p_dist = TRUE)
dev.off()


############################ COMPARAÇÃO ELEVADO.TEMPERATURA-ELEVADO ##################################
#MAPlot
# 1. View all the comparisons conducted
names(all_pairs_airway$merged)

comparison <- all_pairs_airway$merged[["Elev-Elev.Temp"]]

# this will not work unless in a list and will stop, producing an error. E.g.
#diag_plots(merged_in = comparison,
 #          name = "Elev.Temp-Temp",
  #         ma = TRUE)

# Error message:
#merged_in is not a list. If you want to plot with one comparison only,
#put the single dataframe into a list as follows. my_list <- list("name"=
#merged_in)

# 3. Put into a new list as instructed by the error
comparison_list <- list("Elev-Elev.Temp" = comparison)

# this will not work unless the appropriate columns are labelled
# "ID", "AveExpr", and "Adj_PVal"

# 4. Relabel the columns for plotting
# inspecting the column names reveals that the "Adj_PVal" column needs to be specified.


# Here, we will relabel "edger_adj_p" with "Adj_PVal" to use this p-value, using
# the "gsub" command as follows (however, we could also use one of the others or
# the p_max column)

colnames(comparison_list[["Elev-Elev.Temp"]]) <- gsub("edger_adj_p", "Adj_PVal",
                                                 colnames(comparison_list[["Elev-Elev.Temp"]]))

# after label
colnames(comparison_list[["Elev-Elev.Temp"]])

# 5. Plot MA
pdf(file = "/home/janaina/Downloads/ElevTemp-TempMA2Plot_Filter.pdf",   # The directory you want to save the file in
   width = 1, # The width of the plot in inches
    height = 1)
diag_plots(merged_in = comparison_list,
           name = "Elev-Elev.Temp",
           ma = TRUE)
dev.off()

#Volcano
pdf(file = "/home/janaina/Downloads/Elev-Elev.TempVolcanoPlot.pdf")
diag_plots(merged_in = comparison_list,
           name = "Elev-Elev.Temp",
           volcano = TRUE)
dev.off()

#p-value distribution
pdf(file = "/home/janaina/Downloads/ElevTemp-TempPvaluePlot.pdf")
diag_plots(merged_in = comparison_list,
           name = "Elev-Elev.Temp",
           p_dist = TRUE)
dev.off()











