# CompBioFP - Exploring Microbial Community Dynamics

## Louis Ruby Elsalim (2602185431) & Micky Malvino Kusandiwinata (2602174522)

youtube video link for presentation: https://youtu.be/6CjrBZg44CQ

## Short Description


This repository contains a Streamlit app designed to visualize gut microbiome data, and then visualize it. The functionality of the app is to give a more pretty version (intuitive GUI) of comparison between microbial abundance between individuals with Autism Spectrum Disorder (ASD) and healthy individuals. The app provides interactive tools for filtering data, visualizing microbial differences, and exploring insights about specific taxa.

## Features
1. Interactive Filters 🕹️
* Adjust the Minimum Fold Change to display taxa with significant differences in abundance between ASD and Healthy groups. 
* Select the Number of Top Taxa (in Pairs) to control the number of taxa visualized. 

2. Dynamic Data Table 🔢
* Displays filtered data including:
  * Taxa names (cleaned of prefixes like g__ and s__).
  * Mean abundance in ASD and Healthy groups.
  * Fold change values.

3. Bacterial Insights 🦠
* Precompiled insights about specific taxa, including:
  * Effect: Biological role or impact.
  * Affects: Targeted system in the body (example: gut).
  * More in: Prevalence in ASD or Healthy groups.

4. Visualizations 👁️
* Heatmap: Color-coded visualization of microbial abundance in ASD and Healthy groups.
* Bar Plot: Highlights fold changes for selected taxa.
* Shannon Diversity Index:
  * Box plot comparing microbial diversity between groups.
  * Statistical analysis using a Mann-Whitney U test with p-value displayed.

***

To run run the Streamlit app: ```streamlit run app.py```

***

## Dataset
The dataset (ASD meta abundance.csv) contains microbial abundance data for individuals with ASD and healthy individuals. The key columns are:
* Taxonomy which is the names of microbial taxa.
* Columns labled with A* which represents Abundance data for ASD individuals.
* Columns labled with B* which represents data for Healthy individuals.


## How It Works
Preprocessing:

* Cleans taxonomy names by removing prefixes (g__ and s__).
* Calculates mean abundance for ASD and Healthy groups.
* Computes fold change and Shannon diversity index.

User Interaction:
* Sidebar sliders for filtering data dynamically.
* Outputs include a filtered data table, bacterial insights, and visualizations.


Statistical Analysis:
* Uses Mann-Whitney U test to compare Shannon diversity indices between ASD and Healthy groups.
