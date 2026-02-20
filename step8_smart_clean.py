import pandas as pd
import os

CONTEXTS = ['CpG', 'CHG', 'CHH']

def clean_file(context):
    input_file = f"{context}_all.txt"
    output_file = f"clean_{context}_all.txt"
    
    if not os.path.exists(input_file):
        print(f"Skipping {input_file}: Not found.")
        return

    print(f"Processing {input_file}...")
    
    # Beolvasás (feltételezve a tabulátor tagolást a head parancs alapján)
    try:
        df = pd.read_csv(input_file, sep='	', low_memory=False)
    except Exception as e:
        print(f"Hiba a fájl olvasásakor: {e}")
        return

    # Numerikussá alakítás, a hibás értékek (pl. 'null') NaN-ra váltanak
    cols_to_check = ['Mv9', 'Igri', '7Hadd']
    for col in cols_to_check:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    initial_count = len(df)

    # 1. LÉPÉS: Csak azokat a sorokat dobjuk el, ahol MINDHÁROM oszlop NaN.
    # Ha legalább egy mintában van adat, azt megtartjuk.
    df_filtered = df.dropna(subset=cols_to_check, how='all')
    
    # Adatok mentése
    df_filtered.to_csv(output_file, sep='	', index=False)
    
    final_count = len(df_filtered)
    print(f"  Eredeti sorok: {initial_count}")
    print(f"  Megmaradt sorok: {final_count}")
    print(f"  Eltávolítva (üres sorok): {initial_count - final_count}")
    print(f"  Mentve: {output_file}")
    print("-" * 30)

if __name__ == "__main__":
    for context in CONTEXTS:
        clean_file(context)
