import streamlit as st
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import mannwhitneyu

# Load the dataset
file_path = 'ASD meta abundance.csv'
data = pd.read_csv(file_path)

# Remove 'g__' and 's__' from bacterial names
if 'Taxonomy' in data.columns:
    data['Taxonomy'] = data['Taxonomy'].str.replace('g__', '', regex=True).str.replace('s__', '', regex=True)

# Precompute columns for ASD and Healthy
asd_columns = [col for col in data.columns if col.startswith('A')]
healthy_columns = [col for col in data.columns if col.startswith('B')]

# Calculate mean abundances
data['ASD_mean'] = data[asd_columns].mean(axis=1)
data['Healthy_mean'] = data[healthy_columns].mean(axis=1)

# Calculate fold change
data['Fold_Change'] = data['ASD_mean'] / data['Healthy_mean']
data['Fold_Change'] = data['Fold_Change'].replace([float('inf'), -float('inf')], float('nan'))

# Shannon Index calculation function
def calculate_shannon_index(abundance_matrix):
    proportions = abundance_matrix.div(abundance_matrix.sum(axis=0), axis=1)  # Normalize abundances to proportions
    shannon_index = -np.sum(proportions * np.log(proportions + 1e-10), axis=0)  # Add small value to avoid log(0)
    return shannon_index

# Calculate Shannon Index for ASD and Healthy groups
asd_shannon_index = calculate_shannon_index(data[asd_columns])
healthy_shannon_index = calculate_shannon_index(data[healthy_columns])

# Add Shannon Index to the dataset (as averages)
data['ASD_Shannon'] = asd_shannon_index.mean()
data['Healthy_Shannon'] = healthy_shannon_index.mean()

# Perform Mann-Whitney U test for Shannon Index comparison
stat, p_value = mannwhitneyu(asd_shannon_index, healthy_shannon_index)

# Function to provide bacterial insights
def bacterial_insights(taxon):
    insights = {
        "Bacteroides;Bacteroides stercoris": {
            "Effect": "Part of normal gut flora, but overgrowth can lead to inflammation.",
            "Affects": "Gut inflammation",
            "More in": "ASD",
        },
        "Unclassified;Firmicutes bacterium CAG:95": {
            "Effect": "Associated with energy harvest and metabolic changes.",
            "Affects": "Gut metabolism",
            "More in": "ASD",
        },
        "Bacteroides;Bacteroides plebeius": {
            "Effect": "Involved in carbohydrate metabolism, overgrowth disrupts homeostasis.",
            "Affects": "Gut homeostasis",
            "More in": "ASD",
        },
        "Bacteroides;Bacteroides plebeius CAG:211": {
            "Effect": "Involved in carbohydrate metabolism, overgrowth disrupts homeostasis.",
            "Affects": "Gut homeostasis",
            "More in": "ASD",
        },
        "Desulfovibrio;Desulfovibrio piger": {
            "Effect": "Produces hydrogen sulfide, leading to gut irritation.",
            "Affects": "Gut lining and inflammation",
            "More in": "ASD",
        },
        "Prevotella;Prevotella sp. CAG:279": {
            "Effect": "Linked to inflammation and altered gut-brain signaling.",
            "Affects": "Gut-brain axis",
            "More in": "ASD",
        },
        "Bacteroides;Bacteroides finegoldii": {
            "Effect": "Identified as a chondroitin sulfate (CS)-degrading bacterium in the human gut microbiota.",
            "Affects": "Metabolism of glycosaminoglycans",
            "More in": "ASD",
        },
        "Bacteroides;Bacteroides salyersiae": {
            "Effect": "Demonstrated potent chondroitin sulfate-degrading activity in the human gut microbiota.",
            "Affects": "Plays a role in the degradation of chondroitin sulfate",
            "More in": "ASD",
        },
        "Sutterella;Sutterella sp. CAG:397": {
            "Effect": "Associated with gastrointestinal distress.",
            "Affects": "Gut health",
            "More in": "ASD",
        },
        "Clostridium;Clostridium sp. CAG:417": {
            "Effect": "Some species produce toxins affecting gut and nervous system.",
            "Affects": "Gut and nervous system",
            "More in": "ASD",
        },
        "Fusobacterium;Fusobacterium mortiferum": {
            "Effect": "Linked to inflammatory diseases.",
            "Affects": "Inflammation",
            "More in": "ASD",
        },
    }
    return insights.get(taxon, {"Effect": "No specific information available.", "Affects": "Unknown", "More in": "Unknown"})

# Streamlit App Title
st.title("Gut Microbiome Analysis: ASD vs Healthy")

# Add a brief description of the project
st.markdown("""
### About this App
This interactive tool compares microbial abundance and diversity between individuals with Autism Spectrum Disorder (ASD) and healthy individuals. It includes Shannon Index calculations to assess the diversity of microbial communities in each group.

### How to Use
1. Use the sliders on the left to adjust:
   - **Minimum Fold Change**: Only show taxa with a fold change above this threshold.
   - **Number of Top Taxa (in Pairs)**: Select the number of top taxa to display as pairs (ASD and Healthy).
2. View the following:
   - **Filtered Data Table**: Shows the selected taxa with their calculated mean abundance, fold change, and Shannon index.
   - **Heatmap**: Visualizes mean abundance of taxa in ASD and Healthy groups.
   - **Bar Plot**: Highlights fold changes for the selected taxa.
   - **Shannon Index Comparison**: Compares diversity (Shannon Index) for ASD and Healthy groups.

### Why This Matters
Studying the gut microbiome can help uncover potential biomarkers or therapeutic targets for improving the quality of life for individuals with ASD.
""")

# Sidebar controls
st.sidebar.header("Filter Parameters")
fold_change_threshold = st.sidebar.slider("Minimum Fold Change", min_value=1.0, max_value=10.0, value=2.0, step=0.1)
top_n_taxa_pairs = st.sidebar.slider("Number of Top Taxa (in Pairs)", min_value=1, max_value=25, value=1, step=1)

# Adjust the dataset for filtering and pairing
filtered_data = data[data['Fold_Change'] >= fold_change_threshold]

# Select top taxa based on ASD mean abundance and ensure pairing with Healthy
top_asd_taxa = filtered_data.nlargest(top_n_taxa_pairs, 'ASD_mean')
top_healthy_taxa = filtered_data.nlargest(top_n_taxa_pairs, 'Healthy_mean')

# Combine the two groups for a paired dataset
paired_taxa = pd.concat([top_asd_taxa, top_healthy_taxa]).drop_duplicates(subset='Taxonomy', keep='first')

# Display filtered data in Streamlit
st.subheader("Filtered Data")
st.dataframe(paired_taxa[['Taxonomy', 'ASD_mean', 'Healthy_mean', 'Fold_Change']])

# Display bacterial insights
st.subheader("Bacterial Insights")
for index, row in paired_taxa.iterrows():
    insights = bacterial_insights(row['Taxonomy'])
    st.markdown(f"**{row['Taxonomy']}**: {insights['Effect']} (Affects: {insights['Affects']}, More in: {insights['More in']})")

# Heatmap visualization
st.subheader("Heatmap of Selected Taxa")
if not paired_taxa.empty:
    # Prepare heatmap data
    heatmap_data = paired_taxa[['Taxonomy', 'ASD_mean', 'Healthy_mean']].set_index('Taxonomy')
    heatmap_data = heatmap_data.fillna(0)  # Ensure no blanks in the heatmap
    fig, ax = plt.subplots(figsize=(12, 8))
    sns.heatmap(
        heatmap_data,
        annot=False,  # Remove numbers from the cells
        cmap='viridis',
        cbar_kws={'label': 'Mean Abundance'}
    )
    st.pyplot(fig)
else:
    st.warning("No data to display. Adjust the filter parameters.")

# Bar plot visualization
st.subheader("Fold Change of Selected Taxa")
if not paired_taxa.empty:
    fig, ax = plt.subplots(figsize=(12, 6))
    sns.barplot(x="Fold_Change", y="Taxonomy", data=paired_taxa, palette="coolwarm")
    ax.set_title("Top Taxa by Fold Change (ASD / Healthy)")
    ax.set_xlabel("Fold Change")
    ax.set_ylabel("Taxa")
    st.pyplot(fig)
else:
    st.warning("No data to display. Adjust the filter parameters.")

# Shannon Index Distribution (Boxplot)
st.subheader("Shannon Diversity Index Distribution")
fig, ax = plt.subplots(figsize=(6, 4))
sns.boxplot(data=[asd_shannon_index, healthy_shannon_index], palette='coolwarm')
ax.set_xticklabels(['ASD', 'Healthy'])
ax.set_title("Shannon Diversity Index Distribution")
ax.set_ylabel("Shannon Index")
st.pyplot(fig)

# Display P-value
st.markdown(f"""
### Statistical Test Result
**P-value for Shannon Index Comparison (ASD vs Healthy):** {p_value:.5f}
""")
