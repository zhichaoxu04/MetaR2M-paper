library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)

# -------------------------------
# ---------------- Correlation Matrix
# -------------------------------
# load("/Users/xu/Library/CloudStorage/OneDrive-InsideMDAnderson/MDACC/PW/Cohort Cleaning/RDA_033122.rData")
load("S:/Rotation/PW/Cohort_Cleaning/RDA_033122.rData")
load("S:/Rotation/PW/Cohort_Cleaning/Top10PCs.rData")

load("/Users/xu/Library/CloudStorage/OneDrive-InsideMDAnderson/MDACC/PW/Cohort Cleaning/RDA_script/RDA_033122.rData")
load("/Users/xu/Library/CloudStorage/OneDrive-InsideMDAnderson/MDACC/PW/Cohort Cleaning/RDA_script/Top10PCs.rData")

outcome <- "SBP_adj" 
exposure <- "Age"
CovarName <- c("Age", "Sex", "Currsmk", "BMI", "Alcohol2")

data <- Save_033122 %>%  
  # dplyr::filter(cohort == "offspring") %>% 
  mutate(cohort = as.numeric(as.factor(cohort))) %>%
  dplyr::select(shareid, 
                all_of(CovarName),
                all_of(outcome), 
                starts_with("ID_")) %>% 
  filter(complete.cases(.))  

X <- as.numeric(data %>% dplyr::select(all_of(exposure)) %>% pull())
Y <- as.numeric(data %>% dplyr::select(all_of(outcome)) %>% pull())
M <- as.matrix(data %>% dplyr::select(starts_with("ID_")) )
PCs <- Top10PCs[[1]]
colnames(PCs) <- paste("PC", 1:10, sep = "")
Covar <- data %>% dplyr::select(all_of(CovarName), 
                                -all_of(exposure))  # Other 5 variables
Covar <- Covar %>% bind_cols(PCs)


# load("/Users/xu/Library/CloudStorage/OneDrive-InsideMDAnderson/MDACC/PW/Cohort Cleaning/101023/Result/CFOLS_1_3_iSIS_SBP_adj_Age_FDR.RData")
# load("S:/Rotation/PW/Cohort_Cleaning/RDA_script/022624/Result/iSIS_PC/CFOLS2_1_3_iSIS_SBP_adj_Age_FDR_PC_10.RData")
load("/Users/xu/Library/CloudStorage/OneDrive-InsideMDAnderson/MDACC/PW/Cohort Cleaning/RDA_script/022624/Result/iSIS_PC/CFOLS2_1_3_iSIS_SBP_adj_Age_FDR_PC_10.RData")
SBP_FDR1 <- result$select1
SBP_FDR2 <- result$select2

# load("/Users/xu/Library/CloudStorage/OneDrive-InsideMDAnderson/MDACC/PW/Cohort Cleaning/101023/Result/CFOLS_1_3_iSIS_SBP_adj_Age_NOFDR.RData")
load("/Users/xu/Library/CloudStorage/OneDrive-InsideMDAnderson/MDACC/PW/Cohort Cleaning/RDA_script/022624/Result/iSIS_PC/CFOLS2_1_3_iSIS_SBP_adj_Age_NOFDR_PC_10.RData")
# load("S:/Rotation/PW/Cohort_Cleaning/RDA_script/022624/Result/iSIS_PC/CFOLS2_1_3_iSIS_SBP_adj_Age_NOFDR_PC_10.RData")
SBP_ALL1 <- result$select1
SBP_ALL2 <- result$select2
NonFDR1 <- SBP_ALL1[!SBP_ALL1 %in% SBP_FDR1]
NonFDR2 <- SBP_ALL2[!SBP_ALL2 %in% SBP_FDR2]

# Regress out X from M
MX_All <- function(mmm){
  data1 <- data.frame(Med = mmm,
                      envir = X,
                      Cov = Covar)
  l <- summary(stats::lm('Med ~.', data = data1))
  res <- stats::residuals(l)
  return(res)
}

MX_Alpha <- function(mmm){
  data1 <- data.frame(Med = mmm,
                      envir = X,
                      Cov = Covar)
  l <- summary(stats::lm('Med ~.', data = data1))
  res <- coef(l)["envir", 1]
  return(res)
}

# Function to perform multivariate regression and extract beta coefficients
YM_Beta_Multivariate <- function(selectedM, subsample = 1) {
  # Subset the mediator matrix using the selected indices
  if(subsample == 1){
    idx1 <- 1:(nrow(data)/ 2) 
  }else{
    idx1 <- ((nrow(data)/ 2)+1):nrow(data)
  }
  
  MedMatrix <- M[, selectedM, drop = FALSE]  # Subset based on indices, ensure result is a data frame
  
  # Create a data frame for the regression model
  data2 <- data.frame(Response = Y, 
                      MedMatrix, 
                      Cov = Covar)
  
  # Fit the multivariate regression model
  l <- stats::lm(Response ~ ., data = data2 %>% slice(idx1))  # Include all selected mediators and covariates
  
  # Extract the beta coefficients for the mediators
  betas <- coef(l)[2:(length(selectedM) + 1)]  # Extract coefficients for mediators (ignoring intercept)
  
  return(betas)
}

MedMatrix1 <- M %>% 
  as.data.frame() %>% 
  dplyr::select(all_of(SBP_FDR1), all_of(NonFDR1))
MedResMatrix1 <- apply(MedMatrix1, 2, MX_All)
CorM1 <- cor(MedResMatrix1)
eigenM1 <- eigen(CorM1)$values
c(min(eigenM1), max(eigenM1))
round(quantile(abs(CorM1[1:length(SBP_FDR1), (length(SBP_FDR1)+1):nrow(CorM1)]), probs = c(0, 0.25, 0.5, 0.75, 1)), 3)
Alpha1 <- apply(MedMatrix1, 2, MX_Alpha)

quantile(abs(result$AllAlpha1)[SBP_FDR1], probs = c(0, 0.25, 0.5, 0.75, 1))
quantile(abs(result$AllAlpha2)[SBP_FDR2], probs = c(0, 0.25, 0.5, 0.75, 1))
quantile(abs(result$AllAlpha1)[NonFDR2], probs = c(0, 0.25, 0.5, 0.75, 1))
quantile(abs(result$AllAlpha2)[NonFDR1], probs = c(0, 0.25, 0.5, 0.75, 1))


quantile(abs(result$AllBeta1)[SBP_FDR1], probs = c(0, 0.25, 0.5, 0.75, 1))
quantile(abs(result$AllBeta2)[SBP_FDR2], probs = c(0, 0.25, 0.5, 0.75, 1))
quantile(abs(result$AllBeta1)[-c(SBP_FDR2, NonFDR2)], probs = c(0, 0.25, 0.5, 0.75, 1))
quantile(abs(result$AllBeta2)[-c(SBP_FDR1, NonFDR1)], probs = c(0, 0.25, 0.5, 0.75, 1))


# Run the multivariate regression and get the beta coefficients
quantile(abs(YM_Beta_Multivariate(selectedM = SBP_FDR1, subsample = 1)), probs = c(0, 0.25, 0.5, 0.75, 1))
quantile(abs(YM_Beta_Multivariate(selectedM = SBP_FDR2, subsample = 2)), probs = c(0, 0.25, 0.5, 0.75, 1))
quantile(abs(YM_Beta_Multivariate(selectedM = NonFDR2, subsample = 1)), probs = c(0, 0.25, 0.5, 0.75, 1))
quantile(abs(YM_Beta_Multivariate(selectedM = NonFDR1, subsample = 2)), probs = c(0, 0.25, 0.5, 0.75, 1))


MedMatrix2 <- M %>% 
  as.data.frame() %>% 
  dplyr::select(all_of(SBP_FDR2), all_of(NonFDR2))
MedResMatrix2 <- apply(MedMatrix2, 2, MX_All)
CorM2 <- cor(MedResMatrix2)
eigenM2 <- eigen(CorM2)$values
c(min(eigenM2), max(eigenM2))
round(quantile(abs(CorM2[1:length(SBP_FDR2), (length(SBP_FDR2)+1):nrow(CorM2)]), probs = c(0, 0.25, 0.5, 0.75, 1)), 3)
Alpha2 <- apply(MedMatrix2, 2, MX_Alpha)
quantile(abs(Alpha2)[1:length(SBP_FDR2)], probs = c(0, 0.25, 0.5, 0.75, 1))
quantile(abs(Alpha2)[(1+length(SBP_FDR2)):length(Alpha2)], probs = c(0, 0.25, 0.5, 0.75, 1))

library(circlize)
# col_fun = colorRamp2(c(-0.4, 0, 0.4), c("blue", "white", "red"))

HeatSBP1 <- Heatmap(CorM1, name = "Correlation", 
        row_split = factor(c(rep(1, length(SBP_FDR1)), rep(2, length(NonFDR1))), levels = c(1,2)), 
        column_split = factor(c(rep(1, length(SBP_FDR1)),rep(2, length(NonFDR1))), levels = c(1,2)),
        show_row_names = FALSE, show_column_names = FALSE, 
        row_title = NULL, column_title = NULL,
        cluster_row_slices = FALSE, cluster_column_slices = FALSE,
        # column_title = paste0(": Selected M in 1st subsample"),
        show_row_dend = FALSE, show_column_dend = F,
        row_gap = unit(c(2), "mm"), column_gap =  unit(c(2), "mm"),
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 3:4),
                                                            labels = c("Selected Mediators", "Non"), 
                                                            labels_gp = gpar(col = "white", fontsize = 15))),
        left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 3:4),
                                                             labels = c(paste("p =", length(SBP_FDR1)), paste("p =", length(NonFDR1))), 
                                                             labels_gp = gpar(col = "white", fontsize = 15))))

HeatSBP2 <- Heatmap(CorM2, name = "Correlation", 
                    row_split = factor(c(rep(1, length(SBP_FDR2)), rep(2, length(NonFDR2))), levels = c(1,2)), 
                    column_split = factor(c(rep(1, length(SBP_FDR2)),rep(2, length(NonFDR2))), levels = c(1,2)),
                    show_row_names = FALSE, show_column_names = FALSE, 
                    row_title = NULL, column_title = NULL,
                    cluster_row_slices = FALSE, cluster_column_slices = FALSE,
                    # column_title = paste0(": Selected M in 1st subsample"),
                    show_row_dend = FALSE, show_column_dend = F,
                    row_gap = unit(c(2), "mm"), column_gap =  unit(c(2), "mm"),
                    top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 3:4),
                                                                        labels = c("Selected Mediators", "Non"), 
                                                                        labels_gp = gpar(col = "white", fontsize = 15))),
                    left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 3:4),
                                                                     labels = c(paste("p =", length(SBP_FDR2)), paste("p =", length(NonFDR2))), 
                                                                     labels_gp = gpar(col = "white", fontsize = 15))))





# HDL VS. SEX
outcome <- "HDL" 
exposure <- "Sex"
data <- Save_033122 %>%  
  # dplyr::filter(cohort == "offspring") %>% 
  mutate(cohort = as.numeric(as.factor(cohort))) %>%
  dplyr::select(shareid, 
                all_of(CovarName),
                all_of(outcome), 
                starts_with("ID_")) %>% 
  filter(complete.cases(.))  

X <- as.numeric(data %>% dplyr::select(all_of(exposure)) %>% pull())
Y <- as.numeric(data %>% dplyr::select(all_of(outcome)) %>% pull())
M <- as.matrix(data %>% dplyr::select(starts_with("ID_")) )
PCs <- Top10PCs[[2]]
colnames(PCs) <- paste("PC", 1:10, sep = "")
Covar <- data %>% dplyr::select(all_of(CovarName), 
                                -all_of(exposure))  # Other 5 variables
Covar <- Covar %>% bind_cols(PCs)

# load("/Users/xu/Library/CloudStorage/OneDrive-InsideMDAnderson/MDACC/PW/Cohort Cleaning/101023/Result/CFOLS_1_3_iSIS_HDL_Sex_FDR.RData")
# load("S:/Rotation/PW/Cohort_Cleaning/RDA_script/022624/Result/iSIS_PC/CFOLS2_1_3_iSIS_HDL_Sex_FDR_PC_10.RData")
load("/Users/xu/Library/CloudStorage/OneDrive-InsideMDAnderson/MDACC/PW/Cohort Cleaning/RDA_script/022624/Result/iSIS_PC/CFOLS2_1_3_iSIS_HDL_Sex_FDR_PC_10.RData")
FDR1 <- result$select1
FDR2 <- result$select2

# load("/Users/xu/Library/CloudStorage/OneDrive-InsideMDAnderson/MDACC/PW/Cohort Cleaning/101023/Result/CFOLS_1_3_iSIS_HDL_Sex_NOFDR.RData")
# load("S:/Rotation/PW/Cohort_Cleaning/RDA_script/022624/Result/iSIS_PC/CFOLS2_1_3_iSIS_HDL_Sex_NOFDR_PC_10.RData")
load("/Users/xu/Library/CloudStorage/OneDrive-InsideMDAnderson/MDACC/PW/Cohort Cleaning/RDA_script/022624/Result/iSIS_PC/CFOLS2_1_3_iSIS_HDL_Sex_NOFDR_PC_10.RData")
ALL1 <- result$select1
ALL2 <- result$select2
NonFDR1 <- ALL1[!ALL1 %in% FDR1]
NonFDR2 <- ALL2[!ALL2 %in% FDR2]

MedMatrix1 <- M %>% 
  as.data.frame() %>% 
  dplyr::select(all_of(FDR1), all_of(NonFDR1))
MedResMatrix1 <- apply(MedMatrix1, 2, MX_All)
CorM1 <- cor(MedResMatrix1)
eigenM1 <- eigen(CorM1)$values
c(min(eigenM1), max(eigenM1))
quantile(abs(CorM1[1:length(FDR1), (length(FDR1)+1):nrow(CorM1)]), probs = c(0, 0.25, 0.5, 0.75, 1))


MedMatrix2 <- M %>% 
  as.data.frame() %>% 
  dplyr::select(all_of(FDR2), all_of(NonFDR2))
MedResMatrix2 <- apply(MedMatrix2, 2, MX_All)
CorM2 <- cor(MedResMatrix2)
eigenM2 <- eigen(CorM2)$values
c(min(eigenM2), max(eigenM2))
quantile(abs(CorM2[1:length(FDR2), (length(FDR2)+1):nrow(CorM2)]), probs = c(0, 0.25, 0.5, 0.75, 1))
quantile(abs(rnorm(300, 0, 0.1)), probs = c(0, 0.25, 0.5, 0.75, 1))


quantile(abs(result$AllAlpha1)[FDR1], probs = c(0, 0.25, 0.5, 0.75, 1))
quantile(abs(result$AllAlpha2)[FDR2], probs = c(0, 0.25, 0.5, 0.75, 1))
quantile(abs(result$AllAlpha1)[NonFDR2], probs = c(0, 0.25, 0.5, 0.75, 1))
quantile(abs(result$AllAlpha2)[NonFDR1], probs = c(0, 0.25, 0.5, 0.75, 1))
quantile(abs(rnorm(150, 0, 1.5)), probs = c(0, 0.25, 0.5, 0.75, 1))

quantile(abs(result$AllBeta1)[FDR1], probs = c(0, 0.25, 0.5, 0.75, 1))
quantile(abs(result$AllBeta2)[FDR2], probs = c(0, 0.25, 0.5, 0.75, 1))
quantile(abs(result$AllBeta1)[-c(FDR2, NonFDR2)], probs = c(0, 0.25, 0.5, 0.75, 1))
quantile(abs(result$AllBeta2)[-c(FDR1, NonFDR1)], probs = c(0, 0.25, 0.5, 0.75, 1))


# --- Simulations setting
p <- 1500
TrueM <- 10
p1 <- 150
p2 <- 150
p3 <- 1050
corMat <- diag(p)
corCoef1 <- 0.2
corCoef2 <- 0.2
corMatDat <- matrix(runif((TrueM+p1)^2, corCoef1, corCoef2), nrow = (TrueM+p1), ncol = (TrueM+p1))
corMatDat <- matrix(rnorm((TrueM+p1)^2, 0, sd = 0.1), nrow = (TrueM+p1), ncol = (TrueM+p1))
corMatDat[upper.tri(corMatDat)] <- t(corMatDat)[upper.tri(corMatDat)]
diag(corMatDat) <- 1
corMatDat <- Matrix::nearPD(corMatDat, corr = TRUE, base.matrix = TRUE, maxit = 300)$mat
summary(corMatDat[upper.tri(corMatDat, diag = F)])
corMat[1:(TrueM+p1), 1:(TrueM+p1)] <- corMatDat
diag(corMat) <- 1
eigenM <- eigen(corMat[1:TrueM, 1:TrueM])$values
c(min(eigenM), max(eigenM))

HeatHDL1 <- Heatmap(CorM1, name = "Correlation", 
                    row_split = factor(c(rep(1, length(FDR1)), rep(2, length(NonFDR1))), levels = c(1,2)), 
                    column_split = factor(c(rep(1, length(FDR1)),rep(2, length(NonFDR1))), levels = c(1,2)),
                    show_row_names = FALSE, show_column_names = FALSE, 
                    row_title = NULL, column_title = NULL,
                    cluster_row_slices = FALSE, cluster_column_slices = FALSE,
                    # column_title = paste0(": Selected M in 1st subsample"),
                    show_row_dend = FALSE, show_column_dend = F,
                    row_gap = unit(c(2), "mm"), column_gap =  unit(c(2), "mm"),
                    top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 3:4),
                                                                        labels = c("Selected Mediators", "Non"), 
                                                                        labels_gp = gpar(col = "white", fontsize = 15))),
                    left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 3:4),
                                                                     labels = c(paste("p =", length(FDR1)), paste("p =", length(NonFDR1))), 
                                                                     labels_gp = gpar(col = "white", fontsize = 15))))


HeatHDL2 <- Heatmap(CorM2, name = "Correlation", 
                    row_split = factor(c(rep(1, length(FDR2)), rep(2, length(NonFDR2))), levels = c(1,2)), 
                    column_split = factor(c(rep(1, length(FDR2)),rep(2, length(NonFDR2))), levels = c(1,2)),
                    show_row_names = FALSE, show_column_names = FALSE, 
                    row_title = NULL, column_title = NULL,
                    cluster_row_slices = FALSE, cluster_column_slices = FALSE,
                    # column_title = paste0(": Selected M in 1st subsample"),
                    show_row_dend = FALSE, show_column_dend = F,
                    row_gap = unit(c(2), "mm"), column_gap =  unit(c(2), "mm"),
                    top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 3:4),
                                                                        labels = c("Selected Mediators", "Non"), 
                                                                        labels_gp = gpar(col = "white", fontsize = 15))),
                    left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 3:4),
                                                                     labels = c(paste("p =", length(FDR2)), paste("p =", length(NonFDR2))), 
                                                                     labels_gp = gpar(col = "white", fontsize = 15))))


p1 = grid.grabExpr(draw(HeatSBP1))
p2 = grid.grabExpr(draw(HeatSBP2))
p3 = grid.grabExpr(draw(HeatHDL1))
p4 = grid.grabExpr(draw(HeatHDL2))
cowplot::plot_grid(p1, p2, p3, p4, nrow = 2, labels = c("A","B","C","D"), label_size = 20)


# Assuming HeatSBP1, HeatSBP2, HeatHDL1, and HeatHDL2 are Heatmap objects
ht_list = HeatmapList(HeatSBP1, HeatSBP2, HeatHDL1, HeatHDL2)
# Draw the heatmaps in a 2 by 2 layout
draw(ht_list, heatmap_legend_side = 'bot', annotation_legend_side = 'bot', padding = unit(2, "cm"))

library(gridExtra)
library(grid)

# Assuming you have plots saved as images or can plot them to devices
# grid.arrange() can be used to arrange grobs, not directly images
# For demonstration, let's assume you have grobs: grob1, grob2, grob3, grob4
grid.arrange(grob1, grob2, grob3, grob4, nrow = 2, ncol = 2)


# Load Point estimate for BMIXED
load("S:/Rotation/PW/Cohort_Cleaning/RDA_script/013023/Result/BMIXPoint_1_3.RData")
load("S:/Rotation/PW/Cohort_Cleaning/RDA_script/013023/Result/BMIXPoint_2_3.RData") # Different seed
result1$output
result2$output


load("S:/Rotation/PW/CoverP_Simulation/R2_Script/020923/Result/Sim_80_3_2023_200_0_1.5_2.RData")
library(openxlsx)
openxlsx::write.xlsx(summary, file = "C:/Users/Zxu7/Desktop/SimTem_091823.xlsx")



library(tidyverse)
library(ggplot2)
load("S:/Rotation/PW/Cohort_Cleaning/RDA_script/013023/Result/BMIX_CI_2022_3_1_SBP_adj_Age.RData")
load("S:/Rotation/PW/Cohort_Cleaning/RDA_script/013023/Result/BMIX_CI_2022_3_2_HDL_Sex.RData")

# ------ Calculate the 95% CI
# tem <- do.call("rbind", result) %>%
#   as.data.frame() %>%
#   slice_head(prop = 0.5) %>%  # Select all numeric values
#   mutate_all(~ as.numeric(.))
# 
# apply(tem, 2, quantile, probs = c(0.025, 0.975),na.rm = T)


tem <- NULL
for(i in 1:dim(result)[1]){
  tem_bind <- result[[i]]
  tem <- rbind(tem_bind, tem)
}

apply(tem, 2, quantile, probs = c(0.025, 0.975),na.rm = T)




library(gtools)
Output <- data.frame(Study=1:12, N=NA, Rsq.med=NA, se=NA, CI_low=NA, CI_upper=NA, pab1=NA, pab1_fdr=NA, pab2=NA, pab2_fdr=NA,
                     RYX=NA,RYXM=NA, RYM=NA, Total1=NA, Total2=NA, Direct1=NA, Direct2=NA)
# load("S:/Rotation/PW/Cohort_Cleaning/RDA_script/092023-DGAI/Result/CFOLS_Infer1_1_3_1_HDL_Iscore_2010.RData")

filelist = list.files(path = "S:/Rotation/PW/Cohort_Cleaning/RDA_script/092023-DGAI/Result", 
                      recursive = FALSE,
                      pattern = "CFOLS_NoCov.*2010.RData",
                      full.names = TRUE)
filelist <- gtools::mixedsort(filelist)

for (i in 1:12) {
  load(filelist[i])
  Output[i, 2:17] <- as.numeric(c(result$output[c(34,1,30,3,4,5,6,7,8,21,12,18)], 
                                   result$M1c, result$M2c, 
                                   result$M1Gamma, result$M2Gamma))
}



library(openxlsx)
openxlsx::write.xlsx(Output, file = "C:/Users/Zxu7/Desktop/RDA+DGAI.xlsx")























