import scispacy
import spacy

# Load a SciSpaCy model (you can switch to en_core_sci_md or en_ner_bionlp13cg_md, etc.)
nlp = spacy.load("en_ner_bionlp13cg_md")

# Example abstract
abstract = """
Major depressive disorder (MDD) is a common mental health condition characterized by persistent sadness, 
loss of interest in activities, and cognitive impairments. Advances in neuroimaging have enabled researchers 
to identify brain regions such as the prefrontal cortex and amygdala as key areas implicated in MDD. 
This study investigates the correlation between serotonin levels and structural changes in the hippocampus 
of individuals diagnosed with MDD.
"""

# Run NER
doc = nlp(abstract)

# Print named entities
print("\n🧠 Named Entities:\n")
for ent in doc.ents:
    print(f"→ {ent.text}  |  Label: {ent.label_}")
    