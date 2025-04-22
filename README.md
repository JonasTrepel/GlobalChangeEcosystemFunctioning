# Nitrogen Deposition, Climate Change and Land Modification Are Reshaping Global Ecosystem Functioning

This repository contains all scripts used for the analyses in the manuscript:

> **"Nitrogen Deposition, Climate Change and Land Modification Are Reshaping Global Ecosystem Functioning"**  
> Authors: *Jonas Trepel, Joe Atkinson, Oliver Baines, Matthew Kerr, Elizabeth le Roux, Hannah Rubin, Jens-Christian Svenning, Robert Buitenwerf*  

---

## Repository Structure:

---

### Data Acquisition and Preparation

#### **Remote Sensing Data from Google Earth Engine (via `rgee`)**
- `R/rgee/initialize_rgee.R` – Set up `rgee` for accessing GEE.
- `R/rgee/rgee_annual_<variable>_download.R` – Download annual composites of RS variables (e.g., EVI, temperature).

#### **Protected Areas (PA) Data**
- `R/prep/get_pa_shapes.R` – Combine and prepare WDPA and other PA shapefiles.
- `R/prep/rasterize_pas.R` – Rasterize PA boundaries for analysis.

#### **Nitrogen Deposition**
- `R/prep/calculate_europe_total_n_depo.R` – Aggregate annual N deposition totals for Europe.

#### **Covariate Extraction & Dataset Assembly**
- `R/prep/create_grid_and_extract_covariates.R` – Create global grid and extract environmental covariates.
- `R/prep/find_controls_for_pas.R` – Identify suitable control sites for PAs.
- `R/prep/combine_pas_and_controls.R` – Merge PA and control datasets.
- `R/prep/create_pa_grid_1km.R` – Create 1 km resolution grid within PAs.
- `R/prep/extract_covs_pas_controls.R` – Extract covariates for PA and control datasets.
- `R/prep/get_grid_subset.R` – Subset grid to 400k global, 100k for Europe and USA.
- `R/prep/extract_time_series_data.R` – Extract time series data for key variables.

---

### Analysis

#### **Climate and N Deposition Trends**
- `R/analysis/remotePARTS/get_climate_trends_PARTS_loop.R`

#### **EVI Trend Hypotheses**
- **Global**: `test_mean_evo_hypotheses_grid_PARTS.R`
- **Functional Biomes**: `test_mean_evi_hypotheses_grid_functional_biomes_PARTS.R`
- **Olson Biomes**: `test_mean_evi_hypotheses_grid_olsonbiomes_PARTS.R`

#### **Burned Area Trend Hypotheses**
- `test_burned_area_hypotheses_grid_PARTS.R`
- `test_burned_area_hypotheses_grid_functional_biomes_PARTS.R`
- `test_burned_area_hypotheses_grid_olsonbiomes_PARTS.R`

#### **Green-Up Trend Hypotheses**
- `test_greenup_hypotheses_grid_PARTS.R`
- `test_greenup_hypotheses_grid_functional_biomes_PARTS.R`
- `test_greenup_hypotheses_grid_olsonbiomes_PARTS.R`

#### **PA Effects**
- `test_pa_characteristics_controls_pas_remotePARTS.R` – Analyze protection effect and PA characteristics.

---

### Sensitivity Analyses: Nitrogen Deposition

- `test_zhu_n_depo.R` – Global trends using Zhu et al.'s N deposition dataset.
- `zhu_n_depo_evi_functional_biomes.R` – Fucntional biome specific trends using Zhu et al.'s N deposition dataset.
- `zhu_n_depo_evi_olson_biomes.R` – Olson Biome specific trends using Zhu et al.'s N deposition dataset.
- `n_depo_trend_sensitivity.R` – Test effect of N deposition trend on ecosystem functioning in US and Europe.
- `test_europe_mean_evi_n_depo_syndrome.R` - Test effect of N deposition syndrome (combinantion of low, medium, high trend/mean) on ecosystem functioning in Europe
- `test_usa_mean_evi_n_depo_syndrome.R` - Test effect of N deposition syndrome (combinantion of low, medium, high trend/mean) on ecosystem functioning in the US

---

### Visualization

- `R/viz/grid_figures.R` – Figures based on grid dataset (Figs 1–4).
- `R/viz/pa_figures.R` – Figures based on PA dataset.
- `R/viz/global_change_drivers_comparison_protection.R` – Compare global change drivers in and outside PAs.
- Note that many figure components have to be later on combined in inkscape or a similar software

---

### Functions

- `R/functions/extract_gls_estimates.R` – Extract GLS model estimates.
- `R/functions/move_polygon.R` – Displace polygons (e.g., to find control sites).
- `R/functions/resolve_overlaps.R` – Resolve overlapping PA boundaries.

---

## Requirements

- R ≥ 4.2.0
- Google Earth Engine account and `rgee` setup required for remote sensing extraction.


## Contact

For questions, please contact: *jonas.trepel[at]gmail.com*
