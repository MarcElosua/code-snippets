library(Seurat)
# BiocManager::iBiocManager::install("SpatialExperiment")nstall("SpatialExperiment")
library(dplyr)
library(here)
library(glue)
library(tidyverse)
# library(ggspavis)

source(here::here("misc/paths.R"))

# Load Seurat Spatial Object
se <- "misc/{robj_dir}/20220527_tonsil_atlas_spatial_seurat_obj.rds.rds" %>%
    glue() %>%
    here() %>%
    readRDS(file = .)
rownames(se@meta.data) <- colnames(se)

# Load deconvolution
decon_mtrx_ls <- "{fig1}/{robj_dir}/decon_mtrx_annotation_figure_1.rds" %>%
    glue::glue() %>%
    here::here() %>%
    readRDS(file = .)
se@meta.data <- cbind(se@meta.data, decon_mtrx_ls[[2]])

se <- se[, se$gem_id == "esvq52_nluss5"]
se@images <- se@images[Seurat::Images(se) == "esvq52_nluss5"]

# Exrtact tissue image
img <- SeuratObject::GetImage(se, image = "esvq52_nluss5", mode = "raster")
img <- grid::rasterGrob(
    img,
    interpolate = FALSE,
    width = grid::unit(1, "npc"),
    height = grid::unit(1, "npc"))

# Extract tissue coordinates
coord <- as.matrix(SeuratObject::GetTissueCoordinates(se, image = "esvq52_nluss5"))

# Extract and process metadata
dd <- cbind(se@meta.data,coord) %>%
    pivot_longer(c(Cycling, epithelial, Activated.NBC), names_to = "cell_type", values_to = "proportion") %>%
    arrange(barcode, desc(proportion)) %>% 
    group_by(barcode) %>%
    slice_head(n = 1) %>%
    filter(proportion >= 0.1)

# Make Image plot
p <- ggplot() +
    # Plot underlying image
    annotation_custom(
        grob = img,
        xmin = 0, xmax = ncol(img$raster),
        ymin = 0, ymax = nrow(img$raster)
        ) +
    coord_fixed(
        xlim = c(0, ncol(img$raster)),
        ylim = c(0, nrow(img$raster)))

# Make Y negative
ymax <- max(p$coordinates$limits$y)
dd$coord_y_i <- abs(dd$imagerow - ymax)
    
p + geom_point(
    data = dd,
    aes(
        x = imagecol,
        y = coord_y_i,
        color = cell_type,
        fill = proportion,
        alpha = proportion),
    size = 2) +
    scale_alpha_continuous(range = c(0, 1)) +
    theme_void()
