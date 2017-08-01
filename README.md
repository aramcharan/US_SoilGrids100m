# Soil Property and Class Maps of the Conterminous US at 100 meter Spatial Resolution

Summary: These are results of spatial predictions of soil property and soil classes for the conterminous U.S. at 100 m spatial resolution. Three national U.S. soil point datasets — NCSS Characterization Database, the National Soil Information System (NASIS), and the Rapid Carbon Assessment (RaCA) datasets — were combined with remote sensing images and detailed conventional soil polygon maps, and used to generate complete-coverage gridded predictions of soil properties (percent organic carbon, total nitrogen, bulk density, pH, and percent sand and clay) and classes (taxonomic great group and particle size in the control section). An ensemble model from the machine learning algorithms random forest and gradient boosting, as implemented in R packages ranger and xgboost, were used to generate spatial predictions. The model validation results indicate an average classification accuracy of 60% for great groups, and 66% for texture classes; for soil properties R-square at validation points ranged from 62% for total nitrogen to 87% for pH. This hybrid "SoilGrids+" modelling system that incorporates remote sensing data, global and local predictions of soil properties, traditional soil polygon maps, and machine learning, opens a possibility for combining traditional (soil survey data) with state-of-the-art Machine Learning technology, with an objective to make soil data more accurate, easier to update, more accessible, and easier to use.

![alt text](https://github.com/aramcharan/US_SoilGrids100m/blob/master/Results/US48_SoilGrids100mWorkflow.png)
![alt text](https://github.com/aramcharan/US_SoilGrids100m/blob/master/Results/Fig_USA48_distribution_training_points.png "Training points used to generate spatial predictions for soil-classes and soil properties")

*Please cite as:*

* Ramcharan A.,Hengl T., Nauman T., Brungard C., Waltman S., Wills S., Thompson J. (2017) **[Soil Property and Class Maps of the Conterminous US at 100 meter Spatial Resolution based on a Compilation of National Soil Point Observations and Machine Learning](https://arxiv.org/abs/1705.08323).** Submitted to Soil Science Society of America Journal.

# Download Maps

All maps are available for download under the [Open Database License (ODbl) v1.0](https://opendatacommons.org/licenses/odbl/) and can be downloaded from https://doi.org/10.18113/S1KW2H without restrictions. File names beginning with "TAXgg" are the probability maps (0-100%) of occurrence of the named soil great group for each 100 m pixel; file names with "PSCS" are the probability maps (0-100%) of occurrence of the named modified particle size (inn the control section) for each 100 m pixel. For soil property maps, "sl1", "sl2", "sl3", "sl4", "sl5", "sl6", and "sl7" refers to soil property maps at soil depths 0, 5, 15, 30, 60, 100, and 200 cm respectively. In addition to reported soil properties, preliminary, unvalidated maps of Mg, K, and EC are provided in the dataset at the reported soil depths. 

# Disclaimer

These are results of spatial predictions based on using Machine Learning algorithms attached to the above-listed paper and hence some errors and artifacts are still possible. We aim at updating these maps regularly i.e. as the new training / point data arrives. 
