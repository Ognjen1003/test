import pandas as pd

class LookupTableSingleton:
    _lookup_data = {}

    @staticmethod
    def load_lookup_table(case: str) -> pd.DataFrame:
        # Check if the data for the requested case is already loaded
        if case in LookupTableSingleton._lookup_data:
            print(f"Returning cached data for case: {case}")
            return LookupTableSingleton._lookup_data[case]

        # URL for the lookup table CSV based on the case
        url = f"https://raw.githubusercontent.com/dvulin/lookup/main/lookup_table_{case}.csv"
        
        # Load the data from the CSV and process it
        try:
            df_lookup = pd.read_csv(url)
            df_lookup['mu_L'] = df_lookup['mu_L'] * 0.001  # Convert to Pa·s
            df_lookup['mu_g'] = df_lookup['mu_g'] * 0.001  # Convert to Pa·s
            print(f"Loaded data for case: {case}")

            # Store the data in the cache (persist across requests)
            LookupTableSingleton._lookup_data[case] = df_lookup
            return df_lookup
        except Exception as e:
            print(f"Error loading lookup table for case {case}: {e}")
            raise e
