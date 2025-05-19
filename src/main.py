import sys
from Bio import Entrez
import ssl
ssl._create_default_https_context = ssl._create_unverified_context
from dotenv import load_dotenv
import os

load_dotenv()
Entrez.email = os.getenv('NCBI_EMAIL')

def search_pubmed(query, max_results=100):
    """
    Searches PubMed for a given query and returns a list of article IDs.
    """
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()
        id_list = record["IdList"]
        print(f"Found {len(id_list)} articles for query: '{query}'")
        return id_list
    except Exception as e:
        print("An error occurred during PubMed search:", e)
        return []

def fetch_pubmed_data(id_list):
    """
    Retrieves detailed records for a list of PubMed article IDs.
    """
    if not id_list:
        return None
    try:
        ids = ",".join(id_list)
        handle = Entrez.efetch(db="pubmed", id=ids, rettype="xml", retmode="text")
        data = handle.read()
        handle.close()
        return data
    except Exception as e:
        print("An error occurred while fetching PubMed data:", e)
        return None

def main():
    query = "mental health therapy guidelines"
    
    if len(sys.argv) > 1:
        query = " ".join(sys.argv[1:])
    
    id_list = search_pubmed(query, max_results=100)
    
    data = fetch_pubmed_data(id_list)

    if data and data.strip():
        output_file = "data/pubmed/results.xml"
        try:
            with open(output_file, "w", encoding="utf-8") as f:
                f.write(data.decode('utf-8'))
            print(f"Results written to {output_file}")
        except Exception as e:
            print(f"Error writing to file: {e}")
    else:
        print("No data retrieved or data is empty.")

if __name__ == '__main__':
    main()
