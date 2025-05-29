import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import matplotlib.colors as mcolors
import sys


def set_style() -> None:
    """
    Sets the custom styling for all matplotlib plots.
    """
    custom_style = {
        "axes.facecolor": "#333334",
        # "font.family": "Montserrat",
        "axes.labelcolor": "white",
        "figure.edgecolor": "white",
    }
    plt.rcParams.update(custom_style)

palette = [
        "white", "#8ACE00", "#70A700", "#4A6E00",
        "#3A5700", "#213100",
    ]
cmap = mcolors.LinearSegmentedColormap.from_list("custom_cmap", palette)

set_style()

# Аргументы: [input_tsv, output_png]
input_tsv = sys.argv[1]
output_png = sys.argv[2]

# Читаем данные
df = pd.read_csv(input_tsv, sep='\t', index_col=0)

# Определяем максимальное значение для нормализации
max_value = df.max().max()

# Создаем фигуру
plt.figure(figsize=(max(8, len(df.columns)*0.6), max(6, len(df.index)*0.3)))

# Строим хитмап с аннотациями
ax = sns.heatmap(
    df,
    cmap=cmap,
    linewidths=0.5,
    linecolor='lightgray',
    fontcolor='white',
    square=True,
    annot=True,  # Показывать значения в ячейках
    fmt="d",     # Формат целых чисел
    annot_kws={"fontsize": 6},
    cbar_kws={"label": "Number of copies"}
)

# Настройки внешнего вида
plt.title("Gene copy number", fontsize=14, pad=20)
plt.xlabel("Genes", fontsize=12)
plt.ylabel("Samples", fontsize=12)
plt.xticks(rotation=90, ha='right', fontsize=9)
plt.yticks(fontsize=9)

# Настройка цветовой шкалы
cbar = ax.collections[0].colorbar
cbar.set_label('Number of copies', rotation=270, labelpad=20)

plt.tight_layout()
plt.savefig(output_png, dpi=600, bbox_inches='tight')
plt.close()