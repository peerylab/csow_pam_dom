# csow_pam_dom

This repository contains anonymized California spotted owl encounter histories (i.e., randomly generated cell IDs and scaled latitudinal information) from the Sierra Nevada Bioacoustic Monitoring Program for 2021--2024. Complete data with locational information are owned by the Department of Forest and Wildlife Ecology at the University of Wisconsin-Madison, and are available to qualified researchers by contacting the monitoring program principal investigator (M. Zachariah Peery; email: mpeery@wisc.edu) and requesting access. Additionally, complete R code and monitoring program shapefiles are available by contacting the monitoring program data scientist (Jay Winiarski; email: jwiniarski@wisc.edu).

# Passive acoustic monitoring can provide insights into occupancy dynamics and impacts of disturbance for at-risk species

## Abstract

Climate and land use change are dramatically altering the frequency, intensity, and extent of ecological disturbances, which threatens the persistence of at-risk species. To curb the pace and scale of disturbances, balance management and conservation priorities, and alleviate associated population declines, managers require high-quality information on species’ responses to disturbance and their population trends across broad spatial scales that challenge the capacity of traditional, local-scale monitoring programs. Passive acoustic monitoring is a scalable approach to obtain occurrence data, but the extent to which it can be used to model occupancy dynamics and their environmental drivers remains uncertain. Here, we demonstrate how passive acoustic surveys can be analyzed within a Bayesian dynamic occupancy modeling framework to robustly estimate occupancy dynamics and responses to disturbance in the California spotted owl (*Strix occidentalis occidentalis*), which is threatened by increasingly large, severe ‘megafires.’ From 2021–2024, we collected ~2 million hours of audio from autonomous recording units deployed across seven national forests in the Sierra Nevada, California, USA. Spotted owls were less likely to initially occupy and colonize sites that were severely burned, and more likely to go locally extinct following high-severity fire. Further, we observed declining post-fire occupancy trajectories, particularly when sites burned ≥50% high-severity. Occupancy trends varied by national forest, but declined by 2% across the entire region. Our findings—which closely align with those from intensive, ‘gold standard’ demographic studies—demonstrate that large-scale passive acoustic monitoring paired with dynamic occupancy models can effectively detect species’ responses to disturbance and estimate population trends, offering valuable insights for management across multiple spatial scales. Finally, we provide specific recommendations to help other passive acoustic monitoring programs successfully detect ecological responses to disturbance and track population changes.

## Repository files

### data

 * [data](./data): Folder containing unmarked multiseason data frame with anonymized encounter histories (`data/anonymized_umf.rds`)

### code

 * [01-dynamic_occupancy_model.R](01-dynamic_occupancy_model.R): Fit the dynamic occupancy model with the `ubms` package
 * [02-model_data_summaries.R](02-model_data_summaries.R): Create figures and model summaries for the main manuscript and appendix
