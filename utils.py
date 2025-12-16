import re
import requests
import os
import json
import pandas as pd
from urllib.parse import urlencode
from time import sleep
import hashlib

class Utils:
    @staticmethod
    def download_file(url: str, cache_dir: str = "", only_use_cache: bool = False) -> str:
        """
        Downloads the content of a file from the given URL.

        :param url: The URL of the file to download.
        :return: The content of the file as a string.
        """
        cache_hash = hashlib.md5(url.encode('utf-8')).hexdigest()
        cache_filename = os.path.join(cache_dir, cache_hash + ".cache") 
        
        if cache_dir:
            os.makedirs(cache_dir, exist_ok=True)
            if os.path.exists(cache_filename):
                with open(cache_filename, 'r', encoding='utf-8') as f:
                    return f.read()
                
        if only_use_cache:
            raise FileNotFoundError(f"Cache file not found for URL: {url}")
                
        response = requests.get(url)
        response.raise_for_status()  # Raise an error for bad status codes
        
        if cache_dir:
            with open(cache_filename, 'w', encoding='utf-8') as f:
                f.write(response.text)
                
        return response.text
    
    @staticmethod
    def download_study_metadata(study_accession, page = 1):
        url = f"https://www.ebi.ac.uk/metagenomics/api/v1/analyses?study_accession={study_accession}&include=downloads&pagesize=1000&page={page}"
        
        print(f"Fetching page {page} for study {study_accession}... (URL: {url})")
        
        response = requests.get(url)
        
        response.raise_for_status()
        
        analyses = response.json()
        
        if analyses['meta']['pagination']['count'] == 0:
            print(f"No analyses found for study {study_accession}")
            return None
        
        data = [{
            "data": analyses['data'],
            "included": analyses['included']
        }]
        
        sleep (0.33)  # Be polite and avoid overwhelming the server
        
        if analyses['links']['next'] != None:
            more_data = Utils.download_study_metadata(study_accession, page+1)
            data.extend(more_data)
            
        return data
        
    @staticmethod
    def ena_api_search_request(fields: list, query: str, result="read_run", save_file = "ena_results.json" , save_path = "", use_cache = True):
        filepath = os.path.join(save_path, save_file)

        if use_cache and os.path.exists(filepath):
            print(f"Loading cached ENA API results from {filepath}...")
            with open (filepath, "r", encoding="utf-8") as f:
                return json.load(f)
        
        host = "https://www.ebi.ac.uk/ena/portal/api/search"
        
        query = Utils.clean_text(query)
        
        request_body_data = {
            'fields': ",".join(fields),
            'result': result,
            'format': 'json',
            'query': query
        }
        
        
        # The actual HTTP POST request using the requests library
        try:
            print(f"Sending request to ENA API at {host}...")

            response = requests.post(host, data=request_body_data)
            response.raise_for_status() # Raise an exception for bad status codes (4xx or 5xx)

            print(f"Status Code: {response.status_code}")
            
            json_data = response.json()
            
            # Always store intermediate results to file
            with open(filepath, 'w', encoding='utf-8') as f:
                json.dump(json_data, f, indent=4)
                
            return json_data

        except requests.exceptions.RequestException as e:
            print(f"An error occurred during the API request: {e}")
            
        return None

    @staticmethod
    def clean_text(input_string: str) -> str:
        """
        Cleans a string by replacing multiple spaces with a single space
        and removing all newline characters.
        
        Equivalent to the provided VB.NET CleanText function.
        
        :param input_string: The string to be cleaned.
        :return: The cleaned string.
        """
        # Replace multiple spaces (2 or more) with a single space
        # The \s character class matches any whitespace character (space, tab, newline, etc.)
        result = re.sub(r'\s{2,}', ' ', input_string)
        
        # Remove all newline characters (\n and \r)
        # The original VB was more specific about newlines, but since \s{2,} 
        # already collapses most whitespace, focusing on removing remaining 
        # newlines/carriage returns is usually sufficient, though re.sub(r'[\r\n]', '', result)
        # is the direct equivalent of the VB logic.
        result = re.sub(r'[\r\n]', '', result)
        
        return result