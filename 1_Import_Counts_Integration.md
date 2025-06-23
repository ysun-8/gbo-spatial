Import and data processing for Visium Batch 1 on gGBOs
================
Yusha Sun
2025-06-23

``` r
library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(cowplot)
options(future.globals.maxSize = 1e9)
setwd("~/Library/CloudStorage/Box-Box/GBO_2/Hypoxia_Structure/GBM_Spatial")
source('~/Library/CloudStorage/Box-Box/GBO/SEQ_tools/useful_functions.R')
```

In these sections, we will import the gGBO 10X Visium spatial
transcriptomic sections. To re-run this code, set the directories
accordingly based on the processed data on GEO (link: ).

First, we import FFPE (Batch 2) for UP-10072 and UP-9059 in order.

``` r
GBO <- Load10X_Spatial('~/Library/CloudStorage/Box-Box/GBO_2/Hypoxia_Structure/GBM_Spatial/Visium_FFPE/10072_9059/outs/',
                       filename = 'filtered_feature_bc_matrix.h5')
coords <- GetTissueCoordinates(GBO)
GBO@meta.data$x.coord <- coords[,1]
GBO@meta.data$y.coord  <- coords[,2]
GBO@meta.data$cells <- rownames(GBO@meta.data)

cells <- GBO@meta.data %>%
  dplyr::filter(x.coord < 1350, x.coord > 750,
                y.coord > 500, y.coord < 1300)
combined_10072 <- subset(GBO, cells = cells$cells)
```

    ## Warning: Not validating Centroids objects
    ## Not validating Centroids objects

    ## Warning: Not validating FOV objects
    ## Not validating FOV objects
    ## Not validating FOV objects
    ## Not validating FOV objects
    ## Not validating FOV objects
    ## Not validating FOV objects

    ## Warning: Not validating Seurat objects

``` r
SpatialDimPlot(combined_10072, label = TRUE, label.size = 3, pt.size.factor = 5, image.alpha = 1) + theme(aspect.ratio = 0.8)
```

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](1_Import_Counts_Integration_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
cells <- GBO@meta.data %>%
  dplyr::filter(x.coord < 2200, x.coord > 1590,
                y.coord > 1950, y.coord < 2250)
combined_9059 <- subset(GBO, cells = cells$cells)
```

    ## Warning: Not validating Centroids objects

    ## Warning: Not validating Centroids objects

    ## Warning: Not validating FOV objects
    ## Not validating FOV objects
    ## Not validating FOV objects
    ## Not validating FOV objects
    ## Not validating FOV objects
    ## Not validating FOV objects

    ## Warning: Not validating Seurat objects

``` r
SpatialDimPlot(combined_9059, label = TRUE, label.size = 3, pt.size.factor = 5, image.alpha = 1) + theme(aspect.ratio = 2)
```

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](1_Import_Counts_Integration_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

We next import FFPE (Batch 2) for UP-9121 and UP-7790 in order.

``` r
GBO <- Load10X_Spatial('~/Library/CloudStorage/Box-Box/GBO_2/Hypoxia_Structure/GBM_Spatial/Visium_FFPE/9121_7790/outs/',
                       filename = 'filtered_feature_bc_matrix.h5')
#SpatialDimPlot(GBO, label = TRUE, label.size = 3, pt.size.factor = 2.25) + theme(aspect.ratio = 0.6)
#SpatialDimPlot(GBO, crop = TRUE, label = TRUE, interactive = TRUE)
coords <- GetTissueCoordinates(GBO)
GBO@meta.data$x.coord <- coords[,1]
GBO@meta.data$y.coord  <- coords[,2]
GBO@meta.data$cells <- rownames(GBO@meta.data)

cells <- GBO@meta.data %>%
  dplyr::filter(x.coord < 1950, x.coord > 1300,
                y.coord > 780, y.coord < 1500)
combined_7790 <- subset(GBO, cells = cells$cells)
```

    ## Warning: Not validating Centroids objects
    ## Not validating Centroids objects

    ## Warning: Not validating FOV objects
    ## Not validating FOV objects
    ## Not validating FOV objects
    ## Not validating FOV objects
    ## Not validating FOV objects
    ## Not validating FOV objects

    ## Warning: Not validating Seurat objects

``` r
SpatialDimPlot(combined_7790, label = TRUE, label.size = 3, pt.size.factor = 5, image.alpha = 1) + theme(aspect.ratio = 0.8)
```

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](1_Import_Counts_Integration_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
cells <- GBO@meta.data %>%
  dplyr::filter(x.coord < 1680, x.coord > 1350,
                y.coord > 1700, y.coord < 2050)
combined_9121 <- subset(GBO, cells = cells$cells)
```

    ## Warning: Not validating Centroids objects

    ## Warning: Not validating Centroids objects

    ## Warning: Not validating FOV objects
    ## Not validating FOV objects
    ## Not validating FOV objects
    ## Not validating FOV objects
    ## Not validating FOV objects
    ## Not validating FOV objects

    ## Warning: Not validating Seurat objects

``` r
SpatialDimPlot(combined_9121, label = TRUE, label.size = 3, pt.size.factor = 5, image.alpha = 1) + theme(aspect.ratio = 0.8)
```

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](1_Import_Counts_Integration_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

Next up is FF (Batch 1) for several organoids in slice A.

``` r
GBO <- Load10X_Spatial('~/Library/CloudStorage/Box-Box/GBO_2/Hypoxia_Structure/GBM_Spatial/FF_Visium_10_2024/Visium_A/outs/',
                       filename = 'filtered_feature_bc_matrix.h5')
SpatialDimPlot(GBO, label = TRUE, label.size = 3, pt.size.factor = 2.25) + theme(aspect.ratio = 1.6)
```

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](1_Import_Counts_Integration_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
#SpatialDimPlot(GBO, crop = TRUE, label = TRUE, interactive = TRUE)
coords <- GetTissueCoordinates(GBO)
GBO@meta.data$x.coord <- coords[,1]
GBO@meta.data$y.coord  <- coords[,2]
GBO@meta.data$cells <- rownames(GBO@meta.data)

cells <- GBO@meta.data %>%
  dplyr::filter(x.coord < 3760, x.coord > 950,
                y.coord > 1150, y.coord < 2950)
combined_9121_2 <- subset(GBO, cells = cells$cells)
```

    ## Warning: Not validating Centroids objects
    ## Not validating Centroids objects

    ## Warning: Not validating FOV objects
    ## Not validating FOV objects
    ## Not validating FOV objects
    ## Not validating FOV objects
    ## Not validating FOV objects
    ## Not validating FOV objects

    ## Warning: Not validating Seurat objects

``` r
SpatialDimPlot(combined_9121_2, label = TRUE, label.size = 3, pt.size.factor = 5, image.alpha = 1) + theme(aspect.ratio = 0.8)
```

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](1_Import_Counts_Integration_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
cells <- GBO@meta.data %>%
  dplyr::filter(x.coord < 3340, x.coord > 960,
                y.coord > 3100, y.coord < 4900)
combined_10072_2 <- subset(GBO, cells = cells$cells)
```

    ## Warning: Not validating Centroids objects

    ## Warning: Not validating Centroids objects

    ## Warning: Not validating FOV objects
    ## Not validating FOV objects
    ## Not validating FOV objects
    ## Not validating FOV objects
    ## Not validating FOV objects
    ## Not validating FOV objects

    ## Warning: Not validating Seurat objects

``` r
SpatialDimPlot(combined_10072_2, label = TRUE, label.size = 3, pt.size.factor = 5, image.alpha = 1) + theme(aspect.ratio = 1.2)
```

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](1_Import_Counts_Integration_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->

Next up is FF (Batch 1) for several organoids in slice B. Note that some
organoids on slice B are normal cortical organoids and we not import
those here.

``` r
GBO <- Load10X_Spatial('~/Library/CloudStorage/Box-Box/GBO_2/Hypoxia_Structure/GBM_Spatial/FF_Visium_10_2024/Visium_B/outs/',
                       filename = 'filtered_feature_bc_matrix.h5')
SpatialDimPlot(GBO, label = TRUE, label.size = 3, pt.size.factor = 2.25) + theme(aspect.ratio = 1.6)
```

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](1_Import_Counts_Integration_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
#SpatialDimPlot(GBO, crop = TRUE, label = TRUE, interactive = TRUE)
coords <- GetTissueCoordinates(GBO)
GBO@meta.data$x.coord <- coords[,1]
GBO@meta.data$y.coord  <- coords[,2]
GBO@meta.data$cells <- rownames(GBO@meta.data)

cells <- GBO@meta.data %>%
  dplyr::filter(y.coord < 4200, y.coord > 2500,
                x.coord > 2000, x.coord < 3200)
combined_7790_2 <- subset(GBO, cells = cells$cells)
```

    ## Warning: Not validating Centroids objects
    ## Not validating Centroids objects

    ## Warning: Not validating FOV objects
    ## Not validating FOV objects
    ## Not validating FOV objects
    ## Not validating FOV objects
    ## Not validating FOV objects
    ## Not validating FOV objects

    ## Warning: Not validating Seurat objects

``` r
SpatialDimPlot(combined_7790_2, label = TRUE, label.size = 3, pt.size.factor = 5, image.alpha = 1) + theme(aspect.ratio = 0.8)
```

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](1_Import_Counts_Integration_files/figure-gfm/unnamed-chunk-5-2.png)<!-- -->

``` r
cells <- GBO@meta.data %>%
  dplyr::filter(y.coord < 6300, y.coord > 4400,
                x.coord > 2000, x.coord < 3400)
combined_10072_3 <- subset(GBO, cells = cells$cells)
```

    ## Warning: Not validating Centroids objects

    ## Warning: Not validating Centroids objects

    ## Warning: Not validating FOV objects
    ## Not validating FOV objects
    ## Not validating FOV objects
    ## Not validating FOV objects
    ## Not validating FOV objects
    ## Not validating FOV objects

    ## Warning: Not validating Seurat objects

``` r
SpatialDimPlot(combined_10072_3, label = TRUE, label.size = 3, pt.size.factor = 5, image.alpha = 1) + theme(aspect.ratio = 1.2)
```

    ## Scale for fill is already present.
    ## Adding another scale for fill, which will replace the existing scale.

![](1_Import_Counts_Integration_files/figure-gfm/unnamed-chunk-5-3.png)<!-- -->

We combine all of these Batch 1 and 2 gGBO data and integrate them with
Seurat v5.

``` r
name_list <- c(combined_10072, combined_10072_2, combined_10072_3, combined_7790, combined_7790_2, combined_9059, 
               combined_9121, combined_9121_2)
sliceIDs <- c("slice1", "slice2", "slice3", "slice4", "slice5",
              "slice6", "slice7", "slice8")
line_list <- c('10072', '10072', '10072', '7790', '7790', '9059', '9121', '9121') 
visium.list <- c()
for (i in 1:8){
  print(i)
  seurat <- name_list[[i]]
  seurat$slice <- i
  seurat$sliceid <- sliceIDs[i]
  seurat$line <- line_list[i]
  names(seurat@images) <- sliceIDs[i]
  seurat@images[[sliceIDs[i]]]@key <- sliceIDs[i]
  visium.list[[i]] <- seurat
}
```

    ## [1] 1
    ## [1] 2
    ## [1] 3
    ## [1] 4
    ## [1] 5
    ## [1] 6
    ## [1] 7
    ## [1] 8

``` r
visium.merge <- merge(visium.list[[1]], c(visium.list[[2]], visium.list[[3]], visium.list[[4]], visium.list[[5]], 
                                          visium.list[[6]], visium.list[[7]], visium.list[[8]]))
```

    ## Warning: Some cell names are duplicated across objects provided. Renaming to
    ## enforce unique cell names.

``` r
visium.merge <- visium.merge %>% SCTransform(assay = "Spatial", vars.to.regress = c('nCount_Spatial')) %>% RunPCA(npcs = 30)
```

    ## Running SCTransform on assay: Spatial

    ## Warning: The `slot` argument of `GetAssayData()` is deprecated as of SeuratObject 5.0.0.
    ## ℹ Please use the `layer` argument instead.
    ## ℹ The deprecated feature was likely used in the Seurat package.
    ##   Please report the issue at <https://github.com/satijalab/seurat/issues>.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## vst.flavor='v2' set. Using model with fixed slope and excluding poisson genes.

    ## Calculating cell attributes from input UMI matrix: log_umi

    ## Variance stabilizing transformation of count matrix of size 15969 by 465

    ## Model formula is y ~ log_umi

    ## Get Negative Binomial regression parameters per gene

    ## Using 2000 genes, 465 cells

    ## Found 135 outliers - those will be ignored in fitting/regularization step

    ## Second step: Get residuals using fitted parameters for 15969 genes

    ## Computing corrected count matrix for 15969 genes

    ## Calculating gene attributes

    ## Wall clock passed: Time difference of 1.935733 secs

    ## Determine variable features

    ## Regressing out nCount_Spatial

    ## Centering data matrix

    ## Place corrected count matrix in counts slot

    ## Warning: The `slot` argument of `SetAssayData()` is deprecated as of SeuratObject 5.0.0.
    ## ℹ Please use the `layer` argument instead.
    ## ℹ The deprecated feature was likely used in the Seurat package.
    ##   Please report the issue at <https://github.com/satijalab/seurat/issues>.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## vst.flavor='v2' set. Using model with fixed slope and excluding poisson genes.

    ## Calculating cell attributes from input UMI matrix: log_umi

    ## Variance stabilizing transformation of count matrix of size 13717 by 267

    ## Model formula is y ~ log_umi

    ## Get Negative Binomial regression parameters per gene

    ## Using 2000 genes, 267 cells

    ## Found 151 outliers - those will be ignored in fitting/regularization step

    ## Second step: Get residuals using fitted parameters for 13717 genes

    ## Computing corrected count matrix for 13717 genes

    ## Calculating gene attributes

    ## Wall clock passed: Time difference of 1.252915 secs

    ## Determine variable features

    ## Regressing out nCount_Spatial

    ## Centering data matrix

    ## Place corrected count matrix in counts slot

    ## vst.flavor='v2' set. Using model with fixed slope and excluding poisson genes.

    ## Calculating cell attributes from input UMI matrix: log_umi

    ## Variance stabilizing transformation of count matrix of size 13387 by 255

    ## Model formula is y ~ log_umi

    ## Get Negative Binomial regression parameters per gene

    ## Using 2000 genes, 255 cells

    ## Found 220 outliers - those will be ignored in fitting/regularization step

    ## Second step: Get residuals using fitted parameters for 13387 genes

    ## Computing corrected count matrix for 13387 genes

    ## Calculating gene attributes

    ## Wall clock passed: Time difference of 1.108927 secs

    ## Determine variable features

    ## Regressing out nCount_Spatial

    ## Centering data matrix

    ## Place corrected count matrix in counts slot

    ## vst.flavor='v2' set. Using model with fixed slope and excluding poisson genes.

    ## Calculating cell attributes from input UMI matrix: log_umi

    ## Variance stabilizing transformation of count matrix of size 15677 by 666

    ## Model formula is y ~ log_umi

    ## Get Negative Binomial regression parameters per gene

    ## Using 2000 genes, 666 cells

    ## Found 150 outliers - those will be ignored in fitting/regularization step

    ## Second step: Get residuals using fitted parameters for 15677 genes

    ## Computing corrected count matrix for 15677 genes

    ## Calculating gene attributes

    ## Wall clock passed: Time difference of 1.814265 secs

    ## Determine variable features

    ## Regressing out nCount_Spatial

    ## Centering data matrix

    ## Place corrected count matrix in counts slot

    ## vst.flavor='v2' set. Using model with fixed slope and excluding poisson genes.

    ## Calculating cell attributes from input UMI matrix: log_umi

    ## Variance stabilizing transformation of count matrix of size 13611 by 215

    ## Model formula is y ~ log_umi

    ## Get Negative Binomial regression parameters per gene

    ## Using 2000 genes, 215 cells

    ## Found 225 outliers - those will be ignored in fitting/regularization step

    ## Second step: Get residuals using fitted parameters for 13611 genes

    ## Computing corrected count matrix for 13611 genes

    ## Calculating gene attributes

    ## Wall clock passed: Time difference of 1.080423 secs

    ## Determine variable features

    ## Regressing out nCount_Spatial

    ## Centering data matrix

    ## Place corrected count matrix in counts slot

    ## vst.flavor='v2' set. Using model with fixed slope and excluding poisson genes.

    ## Calculating cell attributes from input UMI matrix: log_umi

    ## Variance stabilizing transformation of count matrix of size 13575 by 168

    ## Model formula is y ~ log_umi

    ## Get Negative Binomial regression parameters per gene

    ## Using 2000 genes, 168 cells

    ## Found 287 outliers - those will be ignored in fitting/regularization step

    ## Second step: Get residuals using fitted parameters for 13575 genes

    ## Computing corrected count matrix for 13575 genes

    ## Calculating gene attributes

    ## Wall clock passed: Time difference of 1.031641 secs

    ## Determine variable features

    ## Regressing out nCount_Spatial

    ## Centering data matrix

    ## Place corrected count matrix in counts slot

    ## vst.flavor='v2' set. Using model with fixed slope and excluding poisson genes.

    ## Calculating cell attributes from input UMI matrix: log_umi

    ## Variance stabilizing transformation of count matrix of size 12739 by 74

    ## Model formula is y ~ log_umi

    ## Get Negative Binomial regression parameters per gene

    ## Using 2000 genes, 74 cells

    ## Found 324 outliers - those will be ignored in fitting/regularization step

    ## Second step: Get residuals using fitted parameters for 12739 genes

    ## Computing corrected count matrix for 12739 genes

    ## Calculating gene attributes

    ## Wall clock passed: Time difference of 0.8203611 secs

    ## Determine variable features

    ## Regressing out nCount_Spatial

    ## Centering data matrix

    ## Place corrected count matrix in counts slot

    ## vst.flavor='v2' set. Using model with fixed slope and excluding poisson genes.

    ## Calculating cell attributes from input UMI matrix: log_umi

    ## Variance stabilizing transformation of count matrix of size 14646 by 340

    ## Model formula is y ~ log_umi

    ## Get Negative Binomial regression parameters per gene

    ## Using 2000 genes, 340 cells

    ## Found 234 outliers - those will be ignored in fitting/regularization step

    ## Second step: Get residuals using fitted parameters for 14646 genes

    ## Computing corrected count matrix for 14646 genes

    ## Calculating gene attributes

    ## Wall clock passed: Time difference of 1.33175 secs

    ## Determine variable features

    ## Regressing out nCount_Spatial

    ## Centering data matrix

    ## Place corrected count matrix in counts slot

    ## Set default assay to SCT

    ## Warning in PrepDR(object = object, features = features, verbose = verbose): The
    ## following 255 features requested have not been scaled (running reduction
    ## without them): MT-CO1, GAPDH, MT-ND3, MALAT1, HIST1H1B, SLC2A3, HMGCS1,
    ## MTRNR2L12, HSPA5, TUBA1B, AL627171.2, PLCG2, BNIP3, PGK1, SNHG6, HSPA1B, RPS13,
    ## TPT1, RPS27, HMGN2, LDHA, EPB41L4A-AS1, RPL32, RPL17, ZFAS1, STMN1, HIST1H1E,
    ## RPL26, RPL10, SFPQ, RPL30, MT-ND5, RPL39, MTRNR2L8, EIF1, HIST1H1D, GAS5,
    ## RPL21, FAU, IER3, IDH1, RPS18, RPL24, CDK2AP1, MYOF, TPI1, RPS19, RPS25, H2AFV,
    ## AC097534.2, HSP90B1, MSMO1, SEZ6, RPL36, SCRG1, HAPLN1, SNHG25, UBE2S, RRAD,
    ## RPL9, NDUFB1, BTBD17, HSPA1A, RPS21, SNHG8, NTRK2, TSPAN12, LMO3, RPL31,
    ## HNRNPAB, HIST2H2BE, STC1, FUT5, MASP1, YWHAB, SBDS, SRSF3, TMEM255A, CRISPLD1,
    ## PCDHGC3, TMSB15A, PABPN1, PANTR1, SNHG1, AQP1, CACNA2D1, IFIT3, GPX3, NCAM2,
    ## ARID3A, NETO2, NONO, SELENOP, RPL18A, CRNDE, HNRNPH1, AEBP1, GABBR2, TNR,
    ## HLA-B, PLPPR1, CHRNA9, RPL22, COL4A5, TATDN1, SUB1, AC092958.1, RPL12, CHORDC1,
    ## SEPTIN7, RSL1D1, GTF3C6, LINC00461, LINC-PINT, TUBB3, RPL3, SREK1, SRGAP2,
    ## OLFML2A, GIHCG, SLCO1C1, MRPS15, NSG1, NT5E, QPCTL, CXXC4, EIF3M, KLHL4, LAMA1,
    ## RBM4, RPS16, MPPED2, SLITRK2, SOX10, SNHG5, SLN, FCRLA, HIST1H4C, HMGB1,
    ## CXCL10, CALM1, HNRNPA1, SNHG7, RPL34, LHFPL3, ALDOA, RPS23, RPS5, AL355075.4,
    ## TSC22D1, ZNF395, CHI3L2, RPS24, TMEM158, AC009133.1, SOCS1, RPL38, ADGRL3,
    ## AC107398.3, MT1E, COL2A1, UBB, RIMS4, LINC00511, AC027031.2, ORC6, HIST1H4J,
    ## AC010319.3, MEST, WTAP, HLA-C, SLC6A11, GABRB1, SPPL2A, EIF5A, G3BP1, KLHL11,
    ## SLC26A7, NRG1, IGFBP7, CTXND1, IFRD1, MIR210HG, CHAC1, POSTN, RPS15A, SNHG9,
    ## CA8, DLGAP1, LHX9, CKS1B, GOLGA4, ELAVL4, FAM181A, CD34, RPL4, SOX11, HPF1,
    ## C5orf56, HIST2H2AA4, MIR99AHG, HIST1H2AK, KCNMA1, MZT1, LRRN3, FKBP1A,
    ## FGD5-AS1, STX3, NT5C3A, TAFA5, ZNF91, ASS1, HIST1H2BD, NSG2, TIMM8B, NSRP1,
    ## NKAIN3, CHRNA1, RRM2B, SCN3A, RORB, RAB34, PMM2, GNAI1, NACA, PALM2-AKAP2,
    ## CDH20, MRPL23, EPHB2, OLFM1, CD59, MN1, RGS11, GRIA4, SNHG19, CPVL, XYLT1,
    ## RBBP4, ANOS1, CHST6, TJP1, ZBTB12, TSPAN15, NUDT4, TNPO1, ACSS1, EIF3E, LRRC3B,
    ## DBF4B, ARL6IP1, MT-ATP8, FIGN, HAP1, STMP1, CALCRL

    ## PC_ 1 
    ## Positive:  S100B, OLIG1, HES6, ID3, MT-CO3, DLL3, MARCKSL1, COL20A1, MT-CO2, ACAT2 
    ##     MAP2, LMNB1, NES, PCNA, FXYD6, SPRY4, FABP7, SPC24, HNRNPA2B1, METRN 
    ##     CTXN1, MT-ND2, UCHL1, PLAT, DHCR24, PMP2, DBI, CCND1, NCAN, FAM181B 
    ## Negative:  IGFBP5, VEGFA, SLC2A1, ENO2, NDRG1, NRN1, DDIT4, PPP1R3C, ERO1A, BHLHE40 
    ##     CRYAB, ANGPTL4, SCD, HILPDA, P4HA1, BNIP3L, SLC3A2, VIM, ARRDC3, ADM 
    ##     DDIT3, CHPF, ANKRD37, WSB1, CEBPD, PER1, FAM162A, ABCA1, ID2, IRS2 
    ## PC_ 2 
    ## Positive:  ACAT2, DHCR24, MYC, OAZ1, ARC, CRYAB, OLIG1, SCD, PPP1R15B, HIST1H1C 
    ##     PIM3, PCNA, NFIL3, KHDC4, CTXN1, SRSF2, FASN, MVD, FDFT1, HES6 
    ##     DNAJB9, SNRPB, HMGCR, UBA52, RACK1, MYBL2, PRPF3, YBX1, EEF2, AHSA1 
    ## Negative:  MT3, FAM107A, VEGFA, AQP4, APOE, OBSL1, RGMA, GFAP, SPARCL1, F3 
    ##     CEBPD, CHI3L1, TNC, EDNRB, GRIA1, SLC1A3, IGFBP5, AGT, FAM162A, LMO2 
    ##     LPL, CST3, NOL3, EHD2, CHPF, WSB1, SLC4A4, HILPDA, SEMA6D, LRIG1 
    ## PC_ 3 
    ## Positive:  MT-ATP6, MT-CO2, MT-CO3, MT-CYB, MT-ND4, TUBA1A, MT-ND1, MARCKSL1, OLIG1, SCD 
    ##     FABP7, HES6, PTPRZ1, H3F3A, CST3, DBI, MYO9B, CAV1, MT-ND2, NCAN 
    ##     ID3, BCAN, PTN, TTYH1, METRN, GPRC5A, TXNIP, OLIG2, ANGPTL4, DUSP15 
    ## Negative:  CRYAB, FTL, IGFBP5, HILPDA, MT2A, HIST1H1C, MT1X, IGFBP3, SLC3A2, SQSTM1 
    ##     DDIT3, CEBPD, VEGFA, SPC24, TRMT112, ST13, CLU, SOD2, CD44, WEE1 
    ##     DNAJB9, SPRY4, HIST1H2BN, ENO1, PPP1R15A, UBC, ERO1A, ARC, DNAJB1, NAMPT 
    ## PC_ 4 
    ## Positive:  GFAP, AQP4, CHI3L1, HSPB1, CST3, AGT, HOPX, LGALS1, FAM107A, FTH1 
    ##     MT2A, LMO2, SLC1A3, TMSB4X, MT-ND2, SOD2, MPC1, MLC1, SPARC, TUBA4A 
    ##     METRN, TXNIP, FAM181B, GAS7, ATF5, LGALS3, ISG15, CRY1, CEBPB, ARF4 
    ## Negative:  VEGFA, ENO2, NRN1, SLC2A1, WSB1, DDIT4, HILPDA, BHLHE40, COL20A1, NDRG1 
    ##     FAM162A, IGFBP2, IGFBP5, ERO1A, P4HA1, IGFBP3, DLL3, EGLN1, ID2, ENO1 
    ##     NKAIN4, CHD7, ARRDC3, TMEM45A, PLOD2, CD9, PDGFRA, COL9A3, DCX, TRIO 
    ## PC_ 5 
    ## Positive:  HES5, APOE, AQP4, NCAN, ALDOC, SPARCL1, PPDPF, ID3, SOX4, OAZ1 
    ##     DLL3, CKB, PIM3, ADCYAP1R1, ID1, P4HA1, H3F3B, CST3, NKAIN4, SLC3A2 
    ##     ID2, PPP1R15B, PCSK1N, PEA15, GADD45G, CPE, CLU, ADM, DCX, MARCKSL1 
    ## Negative:  AKAP12, CAV1, SCD, ABCA1, MYO9B, BHLHE41, ACAT2, PYGB, FASN, MT2A 
    ##     FABP3, DHCR24, MVD, SAT1, ANXA2, LSS, MT1X, ANGPTL4, PPP1R3C, IGFBP3 
    ##     UPP1, FDPS, EPAS1, NDRG1, TNFRSF12A, EHD2, FGFBP3, CCT6A, DUSP15, DNER

``` r
visium.merge <- IntegrateLayers(
  object = visium.merge, method = CCAIntegration, k.weight = 40, normalization.method = 'SCT',
  orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)
visium.merge <- FindNeighbors(visium.merge, reduction = "integrated.cca", dims = 1:30)
```

    ## Computing nearest neighbor graph

    ## Computing SNN

``` r
visium.merge <- FindClusters(visium.merge, resolution = 0.3, cluster.name = "cca_clusters")
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 2450
    ## Number of edges: 95983
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8765
    ## Number of communities: 5
    ## Elapsed time: 0 seconds

``` r
visium.merge <- RunUMAP(visium.merge, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca", min.dist = 0.25, n.neighbors = 27)
```

    ## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    ## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    ## This message will be shown once per session

    ## 19:23:30 UMAP embedding parameters a = 1.121 b = 1.057

    ## 19:23:30 Read 2450 rows and found 30 numeric columns

    ## 19:23:30 Using Annoy for neighbor search, n_neighbors = 27

    ## 19:23:30 Building Annoy index with metric = cosine, n_trees = 50

    ## 0%   10   20   30   40   50   60   70   80   90   100%

    ## [----|----|----|----|----|----|----|----|----|----|

    ## **************************************************|
    ## 19:23:30 Writing NN index file to temp file /var/folders/12/vmh3qpbx26729gr39pv19d9h0000gn/T//Rtmpk2Yszg/file3dc32dd419d6
    ## 19:23:30 Searching Annoy index using 1 thread, search_k = 2700
    ## 19:23:30 Annoy recall = 100%
    ## 19:23:31 Commencing smooth kNN distance calibration using 1 thread with target n_neighbors = 27
    ## 19:23:31 Initializing from normalized Laplacian + noise (using RSpectra)
    ## 19:23:31 Commencing optimization for 500 epochs, with 99936 positive edges
    ## 19:23:31 Using rng type: pcg
    ## 19:23:33 Optimization finished

``` r
DimPlot(visium.merge, reduction = 'umap.cca', group.by = c('slice', 'cca_clusters', 'line'), pt.size = 0.4) & NoAxes()
```

![](1_Import_Counts_Integration_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
ggsave('Figures/Visium_spots_UMAP.pdf', units = 'in', width = 8, height = 3)

SpatialDimPlot(visium.merge, group.by = c('slice')) & NoAxes()
```

![](1_Import_Counts_Integration_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

``` r
#saveRDS(visium.merge, 'GBO_Visium_CCA_Integrated_110524.RDS')
```

``` r
installed.packages()[names(sessionInfo()$otherPkgs), "Version"]
```

    ##       future      viridis  viridisLite      cowplot    patchwork    lubridate 
    ##     "1.40.0"      "0.6.5"      "0.4.2"      "1.1.3"      "1.3.0"      "1.9.4" 
    ##      forcats      stringr        dplyr        purrr        readr        tidyr 
    ##      "1.0.0"      "1.5.1"      "1.1.4"      "1.0.4"      "2.1.5"      "1.3.1" 
    ##       tibble      ggplot2    tidyverse       Seurat SeuratObject           sp 
    ##      "3.2.1"      "3.5.2"      "2.0.0"      "5.3.0"      "5.1.0"      "2.2-0"
