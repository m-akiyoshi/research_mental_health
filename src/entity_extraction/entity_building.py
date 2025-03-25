import json
import spacy
import os

# Available models:
# en_ner_bionlp13cg_md
# en_core_sci_sm
# en_core_sci_md
# en_core_sci_scibert


# Load the SciSpacy model
model_name = "en_core_sci_scibert"
nlp = spacy.load(model_name)  # or another SciSpacy NER model

# Load the parsed JSO
# Create output file if it doesn't exist
with open(f"data/parsed/parsed_articles.json", "r", encoding="utf-8") as f:
        articles = json.load(f)

# Output structure
ner_results = []

# Process each article
for article in articles:
    title = article.get("title", "")
    abstract = article.get("abstract", "")
    doc = nlp(abstract)

    entities = []
    for ent in doc.ents:
        entities.append({
            "text": ent.text,
            "label": ent.label_,
            "start_char": ent.start_char,
            "end_char": ent.end_char
        })

    ner_results.append({
        "title": title,
        "abstract": abstract,
        "entities": entities
    })

# Save output to JSON
with open(f"data/ner/ner_output_{model_name}.json", "w", encoding="utf-8") as f:
    json.dump(ner_results, f, indent=2, ensure_ascii=False)

print("✅ NER extraction complete! Results saved to ner_output.json.")