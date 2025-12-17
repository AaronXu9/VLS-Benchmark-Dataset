#!/usr/bin/env python3
"""
DrugBench API Client
====================

A template script for accessing DrugBench web services (https://www.drugbench.org/).

This client provides functionality to:
1. Search by UniProt ID (protein ID) to get affinity data
2. Search by SMILES to get affinity data
3. Download similarity matrix data (CSV and SDF files) for given thresholds

Usage:
    python drugbench_client.py

Requirements:
    pip install requests pandas
    
Author: PLATE-VS Team
Date: December 2025
"""

import os
import requests
import pandas as pd
from typing import Optional, List, Dict, Any, Union
from pathlib import Path
import json
import zipfile
import io
import time
import urllib.parse


class DrugBenchClient:
    """
    Client for interacting with DrugBench web services.
    
    The DrugBench database provides protein-ligand affinity data that can be
    accessed through various search methods and bulk downloads.
    """
    
    BASE_URL = "https://www.drugbench.org"
    API_URL = f"{BASE_URL}/api"
    
    # Available qcov levels for similarity matrix
    QCOV_LEVELS = [50, 70, 95, 100]
    
    def __init__(self, timeout: int = 30, output_dir: str = "./drugbench_data"):
        """
        Initialize the DrugBench client.
        
        Args:
            timeout: Request timeout in seconds
            output_dir: Directory to save downloaded files
        """
        self.timeout = timeout
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'DrugBench-Python-Client/1.0',
            'Accept': 'application/json, text/csv, */*'
        })
    
    # =========================================================================
    # Search by UniProt ID / Protein ID
    # =========================================================================
    
    def search_by_uniprot(self, uniprot_id: str, page: int = 1, limit: int = 100) -> Dict[str, Any]:
        """
        Search for affinity data by UniProt protein ID.
        
        Args:
            uniprot_id: UniProt accession ID (e.g., 'P00533')
            page: Page number for pagination
            limit: Number of results per page
            
        Returns:
            Dictionary containing affinity data for the protein
        """
        endpoint = f"{self.API_URL}/molecules"
        
        params = {
            'protein_id': uniprot_id,
            'page': page,
            'limit': limit
        }
        
        try:
            response = self.session.get(endpoint, params=params, timeout=self.timeout)
            response.raise_for_status()
            return response.json()
                
        except requests.exceptions.RequestException as e:
            print(f"Error searching by UniProt ID '{uniprot_id}': {e}")
            return {'error': str(e), 'uniprot_id': uniprot_id}
    
    def get_protein_ligands(self, uniprot_id: str) -> pd.DataFrame:
        """
        Get all ligands and their affinity values for a given protein.
        
        Args:
            uniprot_id: UniProt accession ID
            
        Returns:
            DataFrame with ligand SMILES and affinity data
        """
        # Use the download endpoint to get all data as CSV
        endpoint = f"{self.API_URL}/molecules/download"
        
        payload = {
            'filters': {
                'protein_id': uniprot_id
            }
        }
        
        try:
            response = self.session.post(endpoint, json=payload, timeout=self.timeout)
            response.raise_for_status()
            
            # The response is a CSV string
            return pd.read_csv(io.StringIO(response.text))
            
        except requests.exceptions.RequestException as e:
            print(f"Error getting ligands for '{uniprot_id}': {e}")
            return pd.DataFrame()
    
    # =========================================================================
    # Search by SMILES
    # =========================================================================
    
    def search_by_smiles(self, smiles: str, exact_match: bool = False) -> pd.DataFrame:
        """
        Search for affinity data by SMILES string.
        
        Args:
            smiles: SMILES string of the compound
            exact_match: If True, performs exact match. If False, performs substructure/similarity search.
            
        Returns:
            DataFrame containing affinity data for matching compounds
        """
        if exact_match:
            # Use the molecules endpoint for exact/ILIKE search
            endpoint = f"{self.API_URL}/molecules"
            params = {'smiles': smiles, 'limit': 100}
            
            try:
                response = self.session.get(endpoint, params=params, timeout=self.timeout)
                response.raise_for_status()
                data = response.json()
                return pd.DataFrame(data)
            except requests.exceptions.RequestException as e:
                print(f"Error searching by SMILES: {e}")
                return pd.DataFrame()
        else:
            # Use the similarity search endpoint
            endpoint = f"{self.API_URL}/search/ligand"
            payload = {'smiles': smiles, 'threshold': 0.7, 'limit': 100}
            
            try:
                response = self.session.post(endpoint, json=payload, timeout=self.timeout)
                response.raise_for_status()
                data = response.json()
                return pd.DataFrame(data)
            except requests.exceptions.RequestException as e:
                print(f"Error searching by SMILES: {e}")
                return pd.DataFrame()

    def download_affinity_data(self, query: str, query_type: str = 'uniprot') -> Optional[Path]:
        """
        Download affinity data for a given query to a CSV file.
        
        Args:
            query: The search query (UniProt ID or SMILES)
            query_type: Type of query - 'uniprot' or 'smiles'
            
        Returns:
            Path to the downloaded file, or None if failed
        """
        endpoint = f"{self.API_URL}/molecules/download"
        
        filters = {}
        if query_type == 'uniprot':
            filters['protein_id'] = query
        elif query_type == 'smiles':
            filters['smiles'] = query
        else:
            print(f"Unsupported query type: {query_type}")
            return None
            
        payload = {'filters': filters}
        
        try:
            response = self.session.post(endpoint, json=payload, timeout=self.timeout)
            response.raise_for_status()
            
            # Generate filename
            safe_query = query.replace('/', '_').replace('\\', '_')[:50]
            filename = f"affinity_{query_type}_{safe_query}.csv"
            filepath = self.output_dir / filename
            
            with open(filepath, 'wb') as f:
                f.write(response.content)
            
            print(f"Downloaded: {filepath}")
            return filepath
            
        except requests.exceptions.RequestException as e:
            print(f"Error downloading affinity data: {e}")
            return None
    
    # =========================================================================
    # Similarity Matrix Downloads
    # =========================================================================
    
    def download_similarity_matrix_csv(
        self, 
        similarity_threshold: float = 0.9,
        qcov_level: int = 100
    ) -> Optional[Path]:
        """
        Download similarity matrix data as CSV for a given threshold.
        
        Args:
            similarity_threshold: Similarity threshold (e.g., 0.9, 0.8, 0.7)
            qcov_level: Query coverage level (50, 70, 95, 100)
            
        Returns:
            Path to the downloaded CSV file, or None if failed
        """
        if qcov_level not in self.QCOV_LEVELS:
            print(f"Warning: qcov_level {qcov_level} not in standard levels {self.QCOV_LEVELS}")
            
        endpoint = f"{self.API_URL}/similarity-matrix/download-uniprot"
        
        params = {
            'threshold': similarity_threshold,
            'qcov_level': qcov_level
        }
        
        try:
            response = self.session.get(endpoint, params=params, timeout=self.timeout)
            response.raise_for_status()
            
            filename = f"similarity_matrix_qcov{qcov_level}_threshold_{similarity_threshold}.csv"
            filepath = self.output_dir / filename
            
            with open(filepath, 'wb') as f:
                f.write(response.content)
            
            print(f"Downloaded similarity matrix CSV: {filepath}")
            return filepath
            
        except requests.exceptions.RequestException as e:
            print(f"Error downloading similarity matrix CSV: {e}")
            return None
    
    def download_similarity_sdf(
        self, 
        similarity_threshold: float = 0.9
    ) -> Optional[Path]:
        """
        Download SDF files (zipped) for a given similarity threshold.
        
        Args:
            similarity_threshold: Similarity threshold (e.g., 0.9, 0.8, 0.7)
            
        Returns:
            Path to the downloaded/extracted SDF files directory, or None if failed
        """
        endpoint = f"{self.API_URL}/similarity-matrix/download-sdf"
        
        params = {
            'threshold': similarity_threshold
        }
        
        try:
            # This endpoint redirects to S3, so we follow redirects
            response = self.session.get(endpoint, params=params, timeout=self.timeout, allow_redirects=True)
            response.raise_for_status()
            
            # Save the tar.gz file
            zip_filename = f"similarity_sdf_threshold_{similarity_threshold}.tar.gz"
            zip_filepath = self.output_dir / zip_filename
            
            with open(zip_filepath, 'wb') as f:
                f.write(response.content)
            
            print(f"Downloaded SDF archive: {zip_filepath}")
            return zip_filepath
            
        except requests.exceptions.RequestException as e:
            print(f"Error downloading SDF files: {e}")
            return None
    
    def download_all_similarity_data(
        self, 
        thresholds: List[float] = [0.9, 0.8, 0.7],
        qcov_level: int = 100
    ) -> Dict[float, Dict[str, Optional[Path]]]:
        """
        Download all similarity data (CSV and SDF) for multiple thresholds.
        
        Args:
            thresholds: List of similarity thresholds to download
            qcov_level: Query coverage level for CSV data
            
        Returns:
            Dictionary mapping thresholds to their downloaded file paths
        """
        results = {}
        
        for threshold in thresholds:
            print(f"\nDownloading data for threshold {threshold}...")
            results[threshold] = {
                'csv': self.download_similarity_matrix_csv(threshold, qcov_level),
                'sdf': self.download_similarity_sdf(threshold)
            }
            # Be respectful of the server
            time.sleep(1)
        
        return results
    
    # =========================================================================
    # Utility Methods
    # =========================================================================
    
    def check_service_status(self) -> Dict[str, bool]:
        """
        Check if DrugBench web services are accessible.
        
        Returns:
            Dictionary with service status
        """
        status = {}
        
        endpoints = {
            'main': self.BASE_URL,
            'api': self.API_URL
        }
        
        for name, url in endpoints.items():
            try:
                response = self.session.head(url, timeout=10)
                # 405 is fine for API root if it doesn't support HEAD/GET
                status[name] = response.status_code in [200, 405]
            except requests.exceptions.RequestException:
                status[name] = False
        
        return status


# =============================================================================
# Example Usage
# =============================================================================

def main():
    """
    Main function demonstrating DrugBench API client usage.
    """
    print("=" * 60)
    print("DrugBench API Client - Template Script")
    print("https://www.drugbench.org/")
    print("=" * 60)
    
    client = DrugBenchClient()
    
    # Check status
    print("\nChecking service status...")
    print(client.check_service_status())
    
    # Example 1: Search by UniProt ID
    uniprot_id = "P00533"
    print(f"\nSearching for UniProt ID: {uniprot_id}")
    df = client.get_protein_ligands(uniprot_id)
    if not df.empty:
        print(f"Found {len(df)} ligands. First 5:")
        print(df.head())
    else:
        print("No ligands found or error occurred.")
        
    # Example 2: Download Similarity Data
    print("\nDownloading similarity data (threshold 0.9)...")
    client.download_similarity_matrix_csv(0.9)
    # client.download_similarity_sdf(0.9) # Uncomment to download large files
    
    print("\nDone.")

if __name__ == "__main__":
    main()
