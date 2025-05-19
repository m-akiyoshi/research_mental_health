from Bio import Entrez
from pathlib import Path
from tqdm import tqdm
import time
import ssl
ssl._create_default_https_context = ssl._create_unverified_context
from dotenv import load_dotenv
import os

load_dotenv()
Entrez.email = os.getenv('NCBI_EMAIL')

def search_pmc(query, max_results=50):
    handle = Entrez.esearch(
        db="pmc",
        term=f"{query} AND open access[filter]",
        retmax=max_results
    )
    record = Entrez.read(handle)
    pmc_ids = record["IdList"]
    return pmc_ids

def fetch_pmc_article(pmcid):
    handle = Entrez.efetch(
        db="pmc",
        id=pmcid,
        rettype="full",
        retmode="xml"
    )
    data = handle.read()
    return data

def fetch_full_pmc_articles(query, max_results=50):
    print(f"🔍 Searching PMC for: {query}")
    pmc_ids = search_pmc(query, max_results)
    print(f"✅ Found {len(pmc_ids)} open-access articles")

    output_dir = Path("data/pmc")
    output_dir.mkdir(parents=True, exist_ok=True)

    for pmcid in tqdm(pmc_ids):
        xml_data = fetch_pmc_article(pmcid)
        file_path = output_dir / f"{pmcid}.xml"
        with open(file_path, "w", encoding="utf-8") as f:
            f.write(xml_data.decode('utf-8'))
        time.sleep(0.5)

    print(f"✅ Downloaded {len(pmc_ids)} full-text XML articles to {output_dir}")

if __name__ == '__main__':
    query = "mental health therapy guidelines"
    fetch_full_pmc_articles(query)