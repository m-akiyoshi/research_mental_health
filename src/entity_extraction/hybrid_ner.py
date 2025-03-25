import spacy
from spacy.tokens import Span

# Load both models
nlp_wide = spacy.load("en_core_sci_scibert")          # High recall, low specificity
nlp_specific = spacy.load("en_ner_bionlp13cg_md")     # Low recall, high specificity

def hybrid_ner(text):
    # Step 1: Extract entities using wide model
    doc_wide = nlp_wide(text)
    spans = list(doc_wide.ents)

    # Step 2: Use specific model to label each span
    labeled_ents = []
    for span in spans:
        # Run the specific model just on the span text
        subdoc = nlp_specific(span.text)
        if subdoc.ents:
            # If it finds a specific label, use it
            labeled_ent = {
                "text": span.text,
                "label": subdoc.ents[0].label_,
                "start_char": span.start_char,
                "end_char": span.end_char
            }
        else:
            # Otherwise fallback to generic "ENTITY"
            labeled_ent = {
                "text": span.text,
                "label": "ENTITY",
                "start_char": span.start_char,
                "end_char": span.end_char
            }
        labeled_ents.append(labeled_ent)

    return labeled_ents

# Test abstract
abstract = """
Major depressive disorder (MDD) is a common mental health condition characterized by persistent sadness, 
loss of interest in activities, and cognitive impairments. Advances in neuroimaging have enabled researchers 
to identify brain regions such as the prefrontal cortex and amygdala as key areas implicated in MDD. 
This study investigates the correlation between serotonin levels and structural changes in the hippocampus 
of individuals diagnosed with MDD.
"""

# Run hybrid NER
entities = hybrid_ner(abstract)

# Print output
print("🧠 Hybrid Named Entities:\n")
for ent in entities:
    print(f"→ {ent['text']}  |  Label: {ent['label']}")