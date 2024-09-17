The dataset contains two Excel spreadsheets, one csv file and one R file.

Data Contents: meta_data_allo.csv is the complete dataset used in the meta-regression analyses 

meta_analyses_code.R includes all R code used for the meta-regression analyses and producing figures 

TableS1.xlsx reports the full list of included and excluded studies

TableS2.xlsx reports the complete dataset, including moderator variables. (This file is the xlsx version of the csv file)

The meta_data_allo.csv and Table S2 files have the following columns:

Phylum: Taxonomic rank of the study organism 

Class: Taxonomic rank of the study organism 

Order: Taxonomic rank of the study organism 

Family: Taxonomic rank of the study organism 

Species: Species name of the study organism, as written in the paper 

rotl_spp: Species name of the study organism, as written in the Open Tree of Life (OTL) database. If the original species was not included in the OTL database, a substitution was found from the same genus or family. 

Weapon: Type of weapon studied 

Citation: Author and year of the study included in the meta-analysis 

Study: Identifier for each unique study 

Injured: Number of individuals with injured/damaged weapons 

Total: Total number of individuals observed 

Proportion: Proportion of individuals with damaged weapons (calculated as injured/total) 

Weapon_size: Weapon length of a large male of the study species. Units are variable, but same as Body_size. 

Body_size: Body length of a large male of the study species. Units are variable, but same as Weapon_size. 

Rel_size: Relative size of weapon (calculated as Weapon_size/Body_size) 

Size: Categorical variable (small or large) indicating relative weapon size. (Small weapons are less than 1/3 length of body; large weapons are longer than 1/3 length of body) 

Allometry: Continuous variable indicating the allometric slope of the weapon 

Regenerate: Categorical variable (Y or N) indicating whether the study organism can regenerate its weapon 

Ram: Categorical variable (Y or N) indicating whether the study organism uses its weapon to ram opponents during fights 

Bias_large: Categorical variable (Y or N) indicating whether there was evidence from the study that large individuals suffered higher rates of weapon damage than small individuals

Observation: Identifier for each unique observation

