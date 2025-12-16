import pandas as pd

class Mappings:
    
    FINAL_DWC_SCHEMA = [
        'collectionID',
        'datasetID',
        'basisOfRecord',
        'associatedSequences',
        'materialSampleID',
        'eventDate',
        'year',
        'month',
        'day',
        'decimalLatitude',
        'decimalLongitude',
        'verbatimDepth',
        'measurementID',
        'measurementMethod',
        'measurementType', 
        'measurementUnit', 
        'measurementValue',
        'scientificNameID',
        'scientificName',
        'taxonRank',
        'higherClassification',
        'domain',
        'kingdom',
        'phylum',
        'order',
        'class',
        'family',
        'genus',
        'specificEpithet',
        'occurrenceStatus',
        'occurrenceID'
    ]
    
    @staticmethod
    def parse_tax_string(tax_string):
        """
        Parses a semicolon-separated taxonomic string into a Darwin Core (DwC) dictionary.

        The string uses a prefix__value format (e.g., 'k__Bacteria;p__Firmicutes')
        and maps prefixes to DwC fields, inferring the lowest available 'taxonRank'.

        Parameters:
            tax_string (str): The taxonomic string to parse.

        Returns:
            dict: A dictionary of DwC taxonomic terms and the inferred 'taxonRank'.
        """
        prefix_map = {
            "sk": "domain",   # superkingdom → DwC domain & superkingdom not available yet
            "k": "kingdom",
            "p": "phylum",
            "c": "class",
            "o": "order",
            "f": "family",
            "g": "genus",
            "s": "specificEpithet"
        }
        
        dwc = {
            "taxonRank": None,
            "domain": None,
            "kingdom": None,
            "phylum": None,
            "order": None,
            "class": None,
            "family": None,
            "genus": None,
            "specificEpithet": None
        }

        if not isinstance(tax_string, str):
            return dwc

        for part in tax_string.split(";"):
            if "__" not in part:
                continue

            prefix, value = part.split("__", 1)
            value = value.strip()

            if not value or prefix not in prefix_map:
                continue

            dwc[prefix_map[prefix]] = value

        # Infer rank and scientificName (lowest available rank)
        if dwc["genus"]:
                dwc["taxonRank"] = "genus"

        elif dwc["family"]:
            dwc["taxonRank"] = "family"

        elif dwc["order"]:
            dwc["taxonRank"] = "order"

        elif dwc["class"]:
            dwc["taxonRank"] = "class"

        elif dwc["phylum"]:
            dwc["taxonRank"] = "phylum"

        elif dwc["kingdom"]:
            dwc["taxonRank"] = "kingdom"

        elif dwc["domain"]:
            dwc["taxonRank"] = "domain"

        return dwc
    
    @staticmethod
    def strip_taxonomy_to_dwc(df: pd.DataFrame):
        """
        Applies `parse_tax_string` to the 'annotation_lineage' column of a DataFrame and expands the results into new DwC columns.

        Parameters:
            df (pd.DataFrame): The input DataFrame containing the 'annotation_lineage' column.

        Returns:
            pd.DataFrame: The original DataFrame concatenated with new columns for DwC taxonomic terms.
        """
        # Apply to DataFrame and expand dict into columns
        dwc_df = df["annotation_lineage"].apply(Mappings.parse_tax_string).apply(pd.Series)

        # Concatenate back to original DataFrame
        return pd.concat([df, dwc_df], axis=1)

    @staticmethod
    def add_ymd(df) -> pd.DataFrame:
        """
        Converts the 'collection_date_start' column to datetime objects and extracts 'year', 'month', and 'day' into new columns.

        Parameters:
            df (pd.DataFrame): The input DataFrame containing the 'collection_date_start' column.

        Returns:
            pd.DataFrame: The DataFrame with the date column converted to datetime and new 'year', 'month', and 'day' columns added.
        """
        
        dates = pd.to_datetime(df["collection_date_start"], utc=True, errors="coerce")

        df["collection_date_start"] = dates

        df["year"] = dates.dt.year
        df["month"] = dates.dt.month
        df["day"] = dates.dt.day

        return df
    
    def add_occurrence(df: pd.DataFrame) -> pd.DataFrame:
        """
        Adds the 'occurrenceStatus' column (set to 'present') and generates a unique 'occurrenceID'.

        The 'occurrenceID' is created by joining values from a predefined list of DwC identifier columns using a pipe ('|') separator.

        Parameters:
            df (pd.DataFrame): The input DataFrame.

        Returns:
            pd.DataFrame: The DataFrame with 'occurrenceStatus' and 'occurrenceID' columns added.
        """
        sep = '|'
        columns = ['collectionID','datasetID','associatedSequences','materialSampleID','measurementID','scientificNameID']
        df['occurrenceID'] = (
            df[columns].apply(lambda row: sep.join(row.dropna().astype(str)), axis=1)
        )
        df['occurrenceStatus'] = 'present'
        return df
    
    def add_basis_of_record(df: pd.DataFrame) -> pd.DataFrame:
        """
        Adds the 'basisOfRecord' column to the DataFrame and sets its value to 'digital sequence information' for all rows.

        Parameters:
            df (pd.DataFrame): The input DataFrame.

        Returns:
            pd.DataFrame: The DataFrame with the 'basisOfRecord' column added.
        """
        df['basisOfRecord'] = 'digital sequence information'
        return df
    
    def parse_ena(df: pd.DataFrame) -> pd.DataFrame:
        """
        Maps and renames ENA (European Nucleotide Archive) columns to corresponding Darwin Core (DwC) terms.

        Specifically, it creates new DwC columns like 'collectionID', 'datasetID', 'materialSampleID', 'eventDate',
        'decimalLatitude', 'decimalLongitude', and 'verbatimDepth' from existing ENA columns.

        Parameters:
            df (pd.DataFrame): The input DataFrame containing ENA data.

        Returns:
            pd.DataFrame: The DataFrame with new columns added and populated with DwC mapping values.
        """
        # Executes renames to provide to darwin core
        df['collectionID'] = df['study_accession|secondary_study_accession']
        df['datasetID'] = df['run_accession']
        df['materialSampleID'] = df['sample_accession|secondary_sample_accession']
        df['eventDate'] = df['collection_date_start']
        df['decimalLatitude'] = df['lat']
        df['decimalLongitude'] = df['long']
        df['verbatimDepth'] = df['depth']
        return df
    
    def parse_mgnify(df: pd.DataFrame) -> dict[str, pd.DataFrame]:
        """
        Maps Mgnify-specific columns to Darwin Core (DwC) terms for sequence data and taxonomic annotations.

        Specifically, it creates mappings for 'associatedSequences', 'measurementID', 'measurementMethod', 'scientificNameID',
        'scientificName', and 'higherClassification' from the input DataFrame columns.

        Parameters:
            df (pd.DataFrame): The input DataFrame containing Mgnify metadata columns (e.g., dsi_id, analysis_id, annotation_name).

        Returns:
            pd.DataFrame: The DataFrame with new DwC-compliant columns added and populated.
        """
        df['associatedSequences'] = df['dsi_id']
        df['measurementID'] = df['analysis_id']
        df['measurementMethod'] = df['analysis_name']
        df['scientificNameID'] = df['annotation_id']
        df['scientificName'] = df['annotation_name']
        df['higherClassification'] = df['annotation_lineage']
        
        return df
        
    @staticmethod
    def split_relative_absolute_counts(df: pd.DataFrame) -> dict[str, pd.DataFrame]:
        """
        Splits a DataFrame containing both absolute and relative counts into two separate DataFrames suitable for DwC measurement extension.

        The function creates two DataFrames: one with 'measurementValue' set to 'count_absolute' and 'measurementUnit' set to '#', and one
        with 'measurementValue' set to 'count_relative' and 'measurementUnit' set to '%'.

        Parameters:
            df (pd.DataFrame): The input DataFrame containing 'count_absolute' and 'count_relative' columns.

        Returns:
            dict[str, pd.DataFrame]: A dictionary containing the two DataFrames:
                - "absolute": DataFrame mapped for absolute counts.
                - "relative": DataFrame mapped for relative counts.
        """
        absolute_df = df
        relative_df = df.copy()
        
        absolute_df["measurementValue"] = absolute_df["count_absolute"]
        absolute_df["measurementUnit"] = "#"
        absolute_df["measurementType"] = "absolute count"
        
        relative_df["measurementValue"] = relative_df["count_relative"]
        relative_df["measurementUnit"] = "%"
        relative_df["measurementType"] = "relative count"
        
        return {
            "absolute": absolute_df,
            "relative": relative_df
        }
        
    def template_dataframe() -> pd.DataFrame:
        """
        Creates an empty pandas DataFrame with predefined columns representing a combination of ENA and Mgnify metadata required for DwC mapping.

        Returns:
            pd.DataFrame: An empty DataFrame with the specified column template.
        """
        columns = [
            "study_accession|secondary_study_accession",
            "run_accession",
            "dsi_id",
            "sample_accession|secondary_sample_accession",
            "collection_date_start",
            "lat",
            "long",
            "depth",
            "analysis_id",
            "analysis_name",
            "count_absolute",
            "count_relative",
            "annotation_id",
            "annotation_name",
            "annotation_lineage"
        ]
        
        return pd.DataFrame(columns=columns)
    
    def enforce_schema_order(df: pd.DataFrame) -> pd.DataFrame:
        """
        Selects and reorders columns in the DataFrame according to the FINAL_DWC_SCHEMA.
        """
        
        missing_cols = list(set(Mappings.FINAL_DWC_SCHEMA) - set(df.columns))
        
        if missing_cols:
            print(f"Warning: Adding missing DwC columns: {missing_cols}")
            for col in missing_cols:
                df[col] = None 
        
        return df[Mappings.FINAL_DWC_SCHEMA]
    
    def run_mappings_from_filled_template(df: pd.DataFrame) -> dict[str, pd.DataFrame]:
        """
        Executes the full chain of Darwin Core (DwC) mapping operations on the input DataFrame.

        This involves mapping ENA and Mgnify data, parsing taxonomy strings, extracting date components,
        adding occurrence and basis of record fields, and returning the final DwC-compliant DataFrame.

        Parameters:
            df (pd.DataFrame): The input DataFrame, typically derived from the template and filled with sequence metadata.

        Returns:
            pd.DataFrame: The resulting DataFrame with all DwC fields added and populated.
        """
        df = Mappings.parse_ena(df)
        df = Mappings.parse_mgnify(df)
        df = Mappings.strip_taxonomy_to_dwc(df)
        df = Mappings.add_ymd(df)
        df = Mappings.add_occurrence(df)
        df = Mappings.add_basis_of_record(df)
        
        split_dfs = Mappings.split_relative_absolute_counts(df)
        
        for key,value in split_dfs.items():
            split_dfs[key] = Mappings.enforce_schema_order(value)
        
        return split_dfs