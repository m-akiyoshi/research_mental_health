import spacy

nlp = spacy.load("en_ner_bionlp13cg_md")

abstract = """
Major depressive disorder (MDD) is a common mental health condition characterized by persistent sadness, 
loss of interest in activities, and cognitive impairments. Advances in neuroimaging have enabled researchers 
to identify brain regions such as the prefrontal cortex and amygdala as key areas implicated in MDD. 
This study investigates the correlation between serotonin levels and structural changes in the hippocampus 
of individuals diagnosed with MDD.
"""

doc = nlp(abstract)

print("\n🧠 Named Entities:\n")
for ent in doc.ents:
    print(f"→ {ent.text}  |  Label: {ent.label_}")
    