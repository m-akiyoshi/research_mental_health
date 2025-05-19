import os
import os
import json
import xml.etree.ElementTree as ET
from tqdm import tqdm

def extract_namespace(tag):
    if tag.startswith("{"):
        return tag[1:].split("}")[0]
    return ''

def get_text_from_element(element):
    if element is None:
        return ''
    return ' '.join(element.itertext()).strip()

def parse_article(xml_file_path):
    try:
        tree = ET.parse(xml_file_path)
        root = tree.getroot()

        ns = {
            'mml': 'http://www.w3.org/1998/Math/MathML',
            'xlink': 'http://www.w3.org/1999/xlink'
        }

        article_meta = root.find('.//article-meta')
        if article_meta is None:
            return None

        title_el = article_meta.find('./title-group/article-title')
        abstract_el = article_meta.find('./abstract')
        body_el = root.find('.//body')

        doi_el = article_meta.find('./article-id[@pub-id-type="doi"]')
        pmid_el = article_meta.find('./article-id[@pub-id-type="pmid"]')
        pmc_el = article_meta.find('./article-id[@pub-id-type="pmc"]')

        journal_el = root.find('.//journal-title')

        pub_date_el = article_meta.find('./pub-date')
        pub_date = ''
        if pub_date_el is not None:
            year = pub_date_el.find('year')
            month = pub_date_el.find('month')
            day = pub_date_el.find('day')
            parts = [get_text_from_element(p) for p in [year, month, day] if p is not None]
            pub_date = '-'.join(parts)

        keywords = [get_text_from_element(kw) for kw in article_meta.findall('.//kwd')]

        authors = []
        for contrib in article_meta.findall('.//contrib[@contrib-type="author"]'):
            name_el = contrib.find('./name')
            if name_el is not None:
                surname = get_text_from_element(name_el.find('./surname'))
                given_names = get_text_from_element(name_el.find('./given-names'))
                full_name = f"{given_names} {surname}".strip()
                authors.append(full_name)

        affiliations = [get_text_from_element(aff) for aff in article_meta.findall('.//aff')]

        return {
            'title': get_text_from_element(title_el),
            'abstract': get_text_from_element(abstract_el),
            'body': get_text_from_element(body_el),
            'doi': get_text_from_element(doi_el),
            'pmid': get_text_from_element(pmid_el),
            'pmc': get_text_from_element(pmc_el),
            'journal': get_text_from_element(journal_el),
            'publication_date': pub_date,
            'keywords': keywords,
            'authors': authors,
            'affiliations': affiliations
        }

    except ET.ParseError:
        print(f"Failed to parse {xml_file_path}")
        return None

def parse_directory(input_dir):
    parsed_articles = []

    for filename in tqdm(os.listdir(input_dir)):
        if filename.endswith('.xml'):
            filepath = os.path.join(input_dir, filename)
            article_data = parse_article(filepath)
            if article_data:
                parsed_articles.append(article_data)

    return parsed_articles

if __name__ == '__main__':
    input_directory = 'data/pmc'
    output_file = 'data/parsed/parsed_articles.json'

    results = parse_directory(input_directory)

    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2, ensure_ascii=False)

    print(f"Parsed {len(results)} articles and saved to {output_file}")