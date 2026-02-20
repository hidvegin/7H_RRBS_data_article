import pandas as pd
import os

def filter_context_file(input_file, output_file):
    print(f"Processing {input_file}...")
    
    try:
        # Read the file (assuming tab-delimited based on inspection)
        df = pd.read_csv(input_file, sep='	', low_memory=False)
        
        cols_to_check = ['Mv9', 'Igri', '7Hadd']
        
        # Verify columns exist
        missing_cols = [c for c in cols_to_check if c not in df.columns]
        if missing_cols:
            print(f"  Error: Missing columns {missing_cols} in {input_file}. Skipping.")
            return

        initial_count = len(df)
        
        # Count NaNs in the specific columns row-wise
        # axis=1 checks across columns for each row
        nan_counts = df[cols_to_check].isna().sum(axis=1)
        
        # Keep rows where NaN count is less than 2 (i.e., 0 or 1 NaN allowed)
        # This means we drop rows with 2 or 3 NaNs
        df_filtered = df[nan_counts < 2].copy()
        
        filtered_count = len(df_filtered)
        
        # Replace remaining NaNs in these columns (and potentially others if desired, 
        # but prompt specified "if somewhere NaN remains... empty value")
        # We will replace ALL NaNs in the dataframe with empty string for the output 
        # or just let to_csv handle it with na_rep parameter.
        # The user specifically said "if somewhere still remains NaN... stay empty value".
        # Using na_rep='' in to_csv is the cleanest way to output empty strings for missing data.
        
        # Save to new file
        df_filtered.to_csv(output_file, sep='	', index=False, na_rep='')
        
        print(f"  Done.")
        print(f"  Original rows: {initial_count}")
        print(f"  Filtered rows: {filtered_count}")
        print(f"  Removed rows: {initial_count - filtered_count}")
        print(f"  Saved to: {output_file}")
        print("-" * 30)
        
    except Exception as e:
        print(f"  An error occurred processing {input_file}: {e}")

if __name__ == "__main__":
    files_to_process = [
        ("CpG_all.txt", "filtered_CpG_all.txt"),
        ("CHG_all.txt", "filtered_CHG_all.txt"),
        ("CHH_all.txt", "filtered_CHH_all.txt")
    ]
    
    for input_f, output_f in files_to_process:
        if os.path.exists(input_f):
            filter_context_file(input_f, output_f)
        else:
            print(f"Input file not found: {input_f}")
