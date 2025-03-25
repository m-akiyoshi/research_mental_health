# parse_pmc_xml.py

import os
import json
from lxml import etree
from pathlib import Path
from tqdm import tqdm
from bs4 import BeautifulSoup

def parse_pmc_xml(file_path):
    with open(file_path, 'rb') as f:
        xml_content = f.read()

    soup = BeautifulSoup(xml_content, "xml")

    def extract_text(tag):
        return tag.get_text(strip=True) if tag else None

    title = extract_text(soup.find("article-title"))
    doi = extract_text(soup.find("article-id", {"pub-id-type": "doi"}))
    pmid = extract_text(soup.find("article-id", {"pub-id-type": "pmid"}))
    pmc = extract_text(soup.find("article-id", {"pub-id-type": "pmc"}))
    journal = extract_text(soup.find("journal-title"))

    # Metadata
    pub_date_tag = soup.find("pub-date", {"pub-type": "epub"})
    pub_date = {
        "day": extract_text(pub_date_tag.find("day")),
        "month": extract_text(pub_date_tag.find("month")),
        "year": extract_text(pub_date_tag.find("year")),
    } if pub_date_tag else None

    # Extract abstract
    abstract = soup.find("abstract")
    abstract_text = " ".join([extract_text(p) for p in abstract.find_all("p")]) if abstract else None

    # Extract keywords
    keywords = [extract_text(k) for k in soup.find_all("kwd")]

    # Extract authors
    authors = []
    for contrib in soup.find_all("contrib", {"contrib-type": "author"}):
        name = contrib.find("name")
        if name:
            given = extract_text(name.find("given-names"))
            surname = extract_text(name.find("surname"))
            full_name = f"{given} {surname}".strip()
            authors.append(full_name)

    # Extract affiliations
    affiliations = {}
    for aff in soup.find_all("aff"):
        aff_id = aff.get("id")
        institutions = [extract_text(i) for i in aff.find_all("institution")]
        country = extract_text(aff.find("country"))
        email = extract_text(aff.find("email"))
        affiliations[aff_id] = {
            "institutions": institutions,
            "country": country,
            "email": email
        }

    # Combine into structured JSON
    article_data = {
        "title": title,
        "doi": doi,
        "pmid": pmid,
        "pmc": pmc,
        "journal": journal,
        "publication_date": pub_date,
        "abstract": abstract_text,
        "keywords": keywords,
        "authors": authors,
        "affiliations": affiliations
    }


    return json.dumps(article_data, indent=2, ensure_ascii=False)

def parse_directory(input_dir, output_file):
    parsed = []
    input_dir = Path(input_dir)
    
    for xml_file in tqdm(list(input_dir.glob("*.xml"))):
        try:
            data = parse_pmc_xml(xml_file)
            parsed.append(data)
        except Exception as e:
            print(f"❌ Failed to parse {xml_file.name}: {e}")

    # Create output directory if it doesn't exist
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, "w", encoding="utf-8") as f:
        json.dump(parsed, f, indent=2, ensure_ascii=False)

    print(f"\n✅ Parsed {len(parsed)} articles and saved to {output_file}")

# Example usage
parse_directory("data/pmc", "data/parsed/parsed_articles.json")