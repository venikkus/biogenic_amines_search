import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import matplotlib.colors as mcolors
import sys

def set_style() -> None:
    """
    Sets the custom styling for all matplotlib plots with transparency
    """
    custom_style = {
        "axes.facecolor": "#333334",      # Dark background for plot area
        "figure.facecolor": "none",       # Transparent figure background
        "savefig.facecolor": "none",      # Transparent save background
        "figure.edgecolor": "none",       # No figure border
        "axes.edgecolor": "lightgray",    # Axis edge color
        # "axes.labelcolor": "white",       # Axis label color
        # "xtick.color": "white",           # X-tick color
        # "ytick.color": "white",           # Y-tick color
        # "text.color": "white",            # Global text color
    }
    plt.rcParams.update(custom_style)

palette = [
    "#213100", "#3A5700", "#4A6E00",
    "#70A700", "#8ACE00", "none"
][::-1]  # Reverse the color order

cmap = mcolors.LinearSegmentedColormap.from_list("custom_cmap", palette)

set_style()

input_tsv = sys.argv[1]
output_png = sys.argv[2]

df = pd.read_csv(input_tsv, sep='\t', index_col=0)

plt.figure(figsize=(max(8, len(df.columns)*0.6), max(6, len(df.index)*0.3)),
           facecolor='none')

ax = sns.heatmap(
    df,
    cmap=cmap,
    linewidths=0.5,
    linecolor='lightgray',
    square=True,
    annot=True,
    fmt="d",
    annot_kws={"fontsize": 6},  # White annotations
    cbar_kws={"label": "Number of copies"}
)

plt.title("Gene copy number", fontsize=14, pad=20, color='white')
plt.xlabel("Genes", fontsize=12)
plt.ylabel("Samples", fontsize=12)
plt.xticks(rotation=90, ha='right', fontsize=9)
plt.yticks(fontsize=9)

cbar = ax.collections[0].colorbar
cbar.set_label('Number of copies', rotation=270, labelpad=20, color='white')
cbar.ax.yaxis.set_tick_params(color='white')  # White ticks for colorbar
plt.setp(cbar.ax.yaxis.get_ticklabels(), color='white')  # White tick labels

plt.tight_layout()
plt.savefig(output_png, dpi=600, bbox_inches='tight', transparent=True)
plt.close()