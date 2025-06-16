import pandas as pd
import requests
from pathlib import Path
import time

def get_gene_name(uniprot_id):
    """Get gene name from UniProt API for a given UniProt ID."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            # Get the first gene name if available
            if 'genes' in data and len(data['genes']) > 0 and 'geneName' in data['genes'][0]:
                return data['genes'][0]['geneName'].get('primary', '')
        return ''
    except Exception as e:
        print(f"Error fetching gene name for {uniprot_id}: {str(e)}")
        return ''

def main():
    # Read the input file
    input_file = Path("data/metadata/markerdb_biomarker.csv")
    df = pd.read_csv(input_file, header=0)
    
    # Create a new column for gene names
    print("Fetching gene names from UniProt...")
    gene_names = []
    
    # Process each UniProt ID
    for uniprot_id in df.iloc[:, 0]:  # Assuming UniProt IDs are in the first column
        gene_name = get_gene_name(uniprot_id)
        gene_names.append(gene_name)
        # Add a small delay to be nice to the API
        time.sleep(0.5)
        
    # Add the gene names as a new column
    df['gene_name'] = gene_names
    
    # Save the updated DataFrame
    output_file = input_file.parent / "markerdb_biomarker_with_genes.csv"
    df.to_csv(output_file, index=False)
    print(f"Updated file saved to: {output_file}")

if __name__ == "__main__":
    main() 