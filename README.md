# Global analysis suggests nitrogen deposition as an underestimated driver of vegetation greening

This repository contains all scripts used for the analyses in the manuscript:

> *"Global analysis suggests nitrogen deposition as an underestimated driver of vegetation greening"*
> Authors: Jonas Trepel, Joe Atkinson, Oliver Baines, Matthew Kerr, Elizabeth le Roux, Hannah Rubin, Jens-Christian Svenning, Robert Buitenwerf  

---

## Repository Structure:

---

### Data Acquisition and Preparation
All datasets used in this study are publicly available. Atmospheric nitrogen deposition estimates are taken from Rubin et al. 2023 (https://doi.org/10.5194/acp-23-7091-2023). Climate variables are derived from the ERA5 monthly reanalysis dataset (https://developers.google.com/earth-engine/datasets/catalog/ECMWF_ERA5_MONTHLY). Human land-modification data are obtained from Theobald et al. (2020) (https://zenodo.org/records/3963013). Vegetation data (EVI) are sourced from the MODIS MOD13A1 product (https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD13A1). The alternative nitrogen-deposition product from Zhu et al. is available at https://doi.org/10.6084/m9.figshare.26778574. Fire-occurrence data are taken from the MCD64A1 Collection 6.1 product (https://lpdaac.usgs.gov/products/mcd64a1v061/). Temporal N-deposition datasets can be downloaded for the United States from https://nadp.slh.wisc.edu/committees/tdep/#tdep-maps and for Europe from https://thredds.met.no/thredds/catalog/data/EMEP/2025_Reporting/catalog.html. Protected-area boundaries are obtained from the World Database on Protected Areas (https://www.protectedplanet.net/en/search-areas?geo_type=site).

#### **R/rgee/: Remote Sensing Data from Google Earth Engine (via `rgee`)**
- `R/rgee/initialize_rgee.R` – Set up `rgee` for accessing GEE.
- `R/rgee/rgee_annual_<variable>_download.R` – Download annual composites of RS variables (e.g., EVI, temperature).

#### **R/prep/: Data preparation**
- `R/prep/combine_hmi_and_get_change.R` - calculate difference in human land modification between 2015-2000
- `R/prep/calculate_europe_total_n_depo.R` - Aggregate annual N deposition totals for Europe.
- `R/prep/get_pa_shapes.R` – Combine and prepare WDPA shapefiles.
- `R/prep/rasterize_biomes_and_pas.R` – Rasterize PA boundaries and biomes for easy extraction later on.
- `R/prep/create_grid.R` - create global 5x5 km grid
- `R/prep/extract_covs.R` - extract covariates
- `R/prep/find_controls_for_pas.R` - find suitable control sites for protected areas
- `R/prep/extract_time_series_values.R` - extract time series values to calculate changes in EVI and climate 
- `R/prep/calculate_trends.R` - calculate trends in time series 
- `R/prep/combine_data_for_analysis.R` - Combine final datasets to be used in subsequent analysis

### Analysis

#### **Main analysis**
- `R/analysis/sdmtmb_world_grid.R`- run global greening models for global grid dataset
- `R/analysis/sdmtmb_pa_dataset.R`- run global greening models for protected area dataset 

#### **Sensitivity**
- `R/analysis/sensitivitiy_zhu_n_depo.R` - test relationship between greening and N depo using an alternative product 
- `R/analysis/sensitivitiy_n_depo_trend.R` - test relationship between greening and N depo trend in europe and the USA 


### Visualization

- `R/viz/extract_model_preds_world_grid.R` - extract predictions from the global grid dataset models
- `R/viz/extract_model_preds_pa_dataset.R` - extract model predictions from the protected areas dataset model
- `R/viz/grid_figures.R` – Figures based on grid dataset (Figs 1&2, plus supplementary figures).
- `R/viz/pa_figures.R` – Figures based on PA dataset (Fig. 3 + supplementary figures).
- `R/viz/sensitivity_zhu_et_al_estimates.R` – Plot results of model with alternative N depositon product.

### Functions

- `R/functions/move_polygon.R` – Displace polygons (e.g., to find control sites).
- `R/functions/resolve_overlaps.R` – Resolve overlapping PA boundaries.

---

## Requirements

- R ≥ 4.5.0
- Google Earth Engine account and `rgee` setup required for remote sensing extraction.


## Contact

For questions, please don't hesitate to contact: *jonas.trepel[at]bio.au.dk*
