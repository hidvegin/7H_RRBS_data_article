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
    
    # Read file (assuming tab separation based on head command output)
    try:
        df = pd.read_csv(input_file, sep='\t', low_memory=False)
    except Exception as e:
        print(f"Error reading file: {e}")
        return
    
    # Convert to numeric, invalid values (e.g. 'null') convert to NaN
    cols_to_check = ['Mv9', 'Igri', '7Hadd']
    for col in cols_to_check:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    
    initial_count = len(df)
    
    # STEP 1: Remove only rows where ALL THREE columns are NaN.
    # If there is data in at least one sample, we keep it.
    df_filtered = df.dropna(subset=cols_to_check, how='all')
    
    # Save data
    df_filtered.to_csv(output_file, sep='\t', index=False)
    
    final_count = len(df_filtered)
    print(f"  Original rows: {initial_count}")
    print(f"  Remaining rows: {final_count}")
    print(f"  Removed (empty rows): {initial_count - final_count}")
    print(f"  Saved: {output_file}")
    print("-" * 30)

if __name__ == "__main__":
    for context in CONTEXTS:
        clean_file(context)
