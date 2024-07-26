import pandas as pd
import os
import numpy as np
import seaborn as sns
import sys
import matplotlib.pyplot as plt
import shap


def string_to_float(string):
    return float(string.strip('[]'))


file_path = "/lustrehome/mbossa/Nuses/Analysis/Moliere_electron20000MeV.parquet"

df_e = pd.read_csv(file_path)

df_e['R4'] = df_e['R4'].apply(string_to_float)

plt.figure(figsize=(10, 6))
plt.hist(df_e['R4'], bins=20, color='skyblue')
plt.xlabel('R4 Value (mm)')
#plt.ylabel('Frequency')
plt.title('Histogram of R4 Values')
plt.show()


