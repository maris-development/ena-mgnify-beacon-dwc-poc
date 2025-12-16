
import io
from mappings import Mappings
from utils import *
from time import sleep
import pandas as pd

ONLY_USE_CACHE = True # If True, only use cached data and do not make new requests
# Note: to build partial cache (for testing), run first with this set to False, and 
# comment out 'blocking' steps once you've gathered enought test data for that step (mgnify api calls)

OUTPUT_DIR = "results"

ENA_QUERY= """
(
    (
        tag = "env_tax:marine"
        AND tag = "env_geo:marine"
    ) AND (
        (
            tag != "env_tax:brackish"
            AND tag != "env_tax:freshwater"
            AND tag != "env_tax:terrestrial"
        ) OR (
            tag != "env_geo:coastal"
            AND tag != "env_geo:freshwater"
            AND tag != "env_geo:terrestrial"
        )
    ) OR (
        tag = "env_geo:marine"
        AND (
            tag = "env_geo:coastal"
            OR tag = "env_geo:freshwater"
            OR tag = "env_geo:terrestrial"
        )
        AND tag = "env_tax:marine"
        AND (
            tag = "env_tax:brackish"
            OR tag = "env_tax:freshwater"
            OR tag = "env_tax:terrestrial"
        )
    )
)
AND collection_date = "*"
AND depth >= 0
AND library_strategy = "amplicon"


"""

ENA_QUERY_FIELDS = [
    "run_accession",
    "study_accession",
    "secondary_study_accession",
    "sample_accession",
    "secondary_sample_accession",
    "collection_date_start",
    "lat",
    "lon",
    "depth"
]

if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Step 1: Query ENA for marine samples
    ena_runs = Utils.ena_api_search_request(ENA_QUERY_FIELDS, ENA_QUERY, save_file="ena_marine_samples.json", use_cache=True)
    
    run_accessions = set([entry['run_accession'] for entry in ena_runs])
    
    print(f"Retrieved {len(run_accessions)} run accessions.")
    
    study_accessions = set([entry['study_accession'] for entry in ena_runs])
    
    print (f"Retrieved {len(study_accessions)} unique study accessions.")
    
    # Step 2: Download study metadata from MGnify for each study accession
    all_study_metadata = {}
    
    STUDIES_OUTPUT_FILE = "all_study_metadata.json"
    
    # Load previously downloaded study metadata if available
    if os.path.exists(STUDIES_OUTPUT_FILE):
        with open(STUDIES_OUTPUT_FILE, "r", encoding="utf-8") as f:
            all_study_metadata = json.load(f)

    # Filter out study accessions that have already been downloaded
    study_accessions = [sa for sa in study_accessions if sa not in all_study_metadata.keys()]
    
    if not ONLY_USE_CACHE:
        total = len(study_accessions)
        print (f"{total} study accessions to download metadata for.")
        
        # Download data for each study accession
        count = 0
        for study_accession in study_accessions:
            count += 1
            print (f"({count}/{total}) Downloading metadata for study {study_accession}...")
            
            all_study_metadata[study_accession] = Utils.download_study_metadata(study_accession)
            
            # Store intermediate results after each download
            with open(STUDIES_OUTPUT_FILE, "w", encoding="utf-8") as f:
                json.dump(all_study_metadata, f, indent=4)
    
    
    # Step 3: Find applicable TSV files in the study data
    tsv_files = []
    
    for study_accession, study_datas in all_study_metadata.items():
        if study_datas is None:
            continue
        for study_data in study_datas:
            if study_data["included"] is not None: 
                _tsv_files = [ item for item in study_data["included"] if item.get("links", {}).get("self", "").endswith(".tsv")]
                
                tsv_files.extend(_tsv_files)
            
    # Store intermediate information about TSV files 
    with open("mgnify_tsv_files.json", "w", encoding="utf-8") as out_f:
        json.dump(tsv_files, out_f, indent=4)
        
    # Filter TSV files to only those matching run accessions
    filtered_tsv_files = {}
    
    for tsv in tsv_files:
        file = tsv.get("links", {}).get("self", "")
        
        basename = os.path.basename(file)
        
        # Only keep merged SSU OTU tables
        if basename.lower().endswith("ssu_otu.tsv"):
            nameparts = basename.split("_")
            run_id = nameparts[0]
            is_merged = nameparts[1].lower() == "merged"
            if run_id in run_accessions and is_merged:
                filtered_tsv_files[run_id] = tsv
            
    print (f"{len(filtered_tsv_files)} TSV files after filtering to matching run accessions.")
        
    # Step 4: Process each run and map data
    count = 0
    done = 0
    total = len(ena_runs)
    
    for run in ena_runs:
        count += 1
        
        study_metadata = all_study_metadata.get(run['study_accession'])
        
        if study_metadata is None:
            # print (f"No study metadata found for study {run['study_accession']}")
            continue
        
        run_tsv = filtered_tsv_files.get(run['run_accession'])
        
        if run_tsv is None:
            # Should not happen when running for the whole desired set
            print (f"No TSV file found for run {run['run_accession']}")
            continue
        
        print(f"({count}/{done}/{total}) Processing run {run['run_accession']}...")
        
        # Find the TSV file link, extract analysis ID and name (from attributes & url), and download the TSV content
        run_tsv_url = run_tsv.get("links", {}).get("self", "")
        
        match = re.search(r"/analyses/(.+?)/file/", run_tsv_url)

        analysis_id = match.group(1) if match else None
        analysis_name = run_tsv.get("attributes", {}).get("group-type", "")
        
        # Download and process the TSV content, skipping the first two header lines
        tsv_content = ""
        
        try:
            tsv_content = Utils.download_file(run_tsv_url, cache_dir="tsv_file_cache", only_use_cache=ONLY_USE_CACHE)
        except FileNotFoundError:
            continue
        except Exception as e:
            print (f"Error downloading TSV for run {run['run_accession']}: {str(e)}")
            continue
        
        tsv_content = "\n".join(tsv_content.splitlines()[2:])
        
        # Load TSV content into a DataFrame
        df = pd.read_csv(io.StringIO(tsv_content), sep="\t", names=["otu_id", "abs_count", "taxonomy", "taxid"])
    
        # Populate additional required fields (MGnify + ENA metadata)
        df["study_accession|secondary_study_accession"] = f"{run['study_accession']}|{run['secondary_study_accession']}"
        df["run_accession"] = run['run_accession']
        df['dsi_id'] = df['otu_id']
        df["sample_accession|secondary_sample_accession"] = f"{run['sample_accession']}|{run['secondary_sample_accession']}"
        df["collection_date_start"] = run['collection_date_start']
        df["lat"] = run['lat']
        df["long"] = run['lon']
        df["depth"] = run['depth']
        df["analysis_id"] = analysis_id
        df["analysis_name"] = analysis_name
        df["count_absolute"] = df["abs_count"]
        df["count_relative"] = df["abs_count"] / df["abs_count"].sum() * 100
        df["annotation_id"] = df["taxid"]
        df["annotation_name"] = (
            df["taxonomy"]
            .str.split(";")
            .apply(lambda x: x[-1][3:] if len(x) > 0 else "")
            .str.replace("_", " ")
        )        
        df["annotation_lineage"] = df["taxonomy"]
        
        # Check against template, keeping only fields that are in the template
        template: pd.DataFrame = Mappings.template_dataframe()
        result_df = df.reindex(columns=template.columns)
        
        # Save the intermediate filled template to CSV
        output_filename = f"mapped_{run['run_accession']}.ebi.csv"    
        result_df.to_csv(os.path.join(OUTPUT_DIR, output_filename), index=False)    
        
        # Apply mappings to the DataFrame, this will add DwC fields
        result_dfs = Mappings.run_mappings_from_filled_template(result_df)

        for key, _df in result_dfs.items():
            # Save each of the split DataFrames to separate CSVs
            output_filename_dwc = f"mapped_{run['run_accession']}.{key}.dwc.csv"    
            _df.to_csv(os.path.join(OUTPUT_DIR, output_filename_dwc), index=False)
        
        done += 1
    
    
    # Step 5: Combine all individual result DwC CSVs into a single Parquet file
    for csv_type in ["absolute", "relative"]:
        result_dfs = []
        
        for csv in os.listdir(OUTPUT_DIR):
            if not csv.endswith(f".{csv_type}.dwc.csv"):
                continue
        
            df = pd.read_csv(os.path.join(OUTPUT_DIR, csv))
            
            result_dfs.append(df)
        
        print(f"Combining {len(result_dfs)} individual result CSVs ({csv_type}) into a single Parquet file...")
        
        # Combine all DataFrames
        combined_df = pd.concat(result_dfs, ignore_index=True)
        
        # Bugfix: Ensure numeric fields are properly typed (input data contains some invalid values, e.g. 2000m, N/A, etc.)
        combined_df["verbatimDepth"] = pd.to_numeric(combined_df["verbatimDepth"], errors="coerce")
        
        # Save combined DataFrame to Parquet
        combined_df.to_parquet(os.path.join(OUTPUT_DIR, f"dwc_combined.{csv_type}.parquet"), index=False)
        
    print("Processing complete.")
        
    