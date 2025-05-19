import json
import spacy

print("🔁 Loading models...")
nlp_scibert = spacy.load("en_core_sci_scibert")          
nlp_bionlp = spacy.load("en_ner_bionlp13cg_md")          

def relabel_entity(entity_text, bionlp_doc):
    """Look for entity_text in bionlp_doc and return the most specific label."""
    for ent in bionlp_doc.ents:
        if ent.text.lower() == entity_text.lower():
            return ent.label_
    return "ENTITY"

def extract_entities(text):
    """Extract from sci_scibert and relabel from bionlp."""
    doc_scibert = nlp_scibert(text)
    doc_bionlp = nlp_bionlp(text)

    unique = set()
    entities = []
    for ent in doc_scibert.ents:
        if ent.text.lower() in unique:
            continue
        unique.add(ent.text.lower())
        label = relabel_entity(ent.text, doc_bionlp)
        entities.append({
            "text": ent.text,
            "label": label,
            "start_char": ent.start_char,
            "end_char": ent.end_char
        })
    return entities

with open("data/parsed/parsed_articles.json", "r") as f:
    articles = json.load(f)

output = []
print(f"🔍 Processing {len(articles)} articles...")
for article in articles:
    abstract = article.get("abstract", "")
    entities = extract_entities(abstract)
    output.append({
        "title": article.get("title"),
        "abstract": abstract,
        "entities": entities
    })

with open("data/ner/ner_output_hybrid.json", "w") as f:
    json.dump(output, f, indent=2)

print("✅ Done! Output written to ner_output_hybrid.json")