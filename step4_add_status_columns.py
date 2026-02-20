import pandas as pd
import numpy as np

def add_methylation_status(file_path):
    """
    Hozzáadja a hyper/hypo státusz oszlopokat az annotált metilációs fájlhoz
    a különbség (Difference) oszlopok alapján.
    """
    print(f"Fájl betöltése: {file_path}")
    df = pd.read_csv(file_path)

    # Oszlopnevek meghatározása a fejléc alapján
    mv9_diff_col = 'Difference (Mv9_7Hadd_Logistic regression p<0.05 after correction. Min obs 10)'
    igri_diff_col = 'Difference (Igri_7Hadd_Logistic regression p<0.05 after correction. Min obs 10)'

    def get_status(diff):
        if pd.isna(diff) or diff == '':
            return np.nan
        try:
            val = float(diff)
            return 'hyper' if val > 0 else 'hypo'
        except:
            return np.nan

    # Státusz oszlopok kiszámítása
    if mv9_diff_col in df.columns:
        df['Mv9_7Hadd_Status'] = df[mv9_diff_col].apply(get_status)
    
    if igri_diff_col in df.columns:
        df['Igri_7Hadd_Status'] = df[igri_diff_col].apply(get_status)

    # Eredmény mentése
    df.to_csv(file_path, index=False)
    print(f"Sikeresen frissítve: {file_path}")

if __name__ == "__main__":
    add_methylation_status('annotated_significant_methylation.csv')
