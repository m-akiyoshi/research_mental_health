import json
import spacy

# Available models:
# en_ner_bionlp13cg_md
# en_core_sci_sm
# en_core_sci_md
# en_core_sci_scibert


model_name = "en_core_sci_scibert"
nlp = spacy.load(model_name)  

with open(f"data/parsed/parsed_articles.json", "r", encoding="utf-8") as f:
        articles = json.load(f)

ner_results = []

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

with open(f"data/ner/ner_output_{model_name}.json", "w", encoding="utf-8") as f:
    json.dump(ner_results, f, indent=2, ensure_ascii=False)

print("✅ NER extraction complete! Results saved to ner_output.json.")