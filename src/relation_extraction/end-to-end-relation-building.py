from transformers import AutoTokenizer, AutoModelForSeq2SeqLM
import torch
import re

model_name = "Babelscape/rebel-large"
# model_name = "dmis-lab/biobert-v1.1"
tokenizer = AutoTokenizer.from_pretrained(model_name)
model = AutoModelForSeq2SeqLM.from_pretrained(model_name)
model.eval()

def extract_triplets(text):
    triplets = []
    relation, subject, relation, object_ = '', '', '', ''
    text = text.strip()
    current = 'x'
    for token in text.replace("<s>", "").replace("<pad>", "").replace("</s>", "").split():
        if token == "<triplet>":
            current = 't'
            if relation != '':
                triplets.append({'head': subject.strip(), 'type': relation.strip(),'tail': object_.strip()})
                relation = ''
            subject = ''
        elif token == "<subj>":
            current = 's'
            if relation != '':
                triplets.append({'head': subject.strip(), 'type': relation.strip(),'tail': object_.strip()})
            object_ = ''
        elif token == "<obj>":
            current = 'o'
            relation = ''
        else:
            if current == 't':
                subject += ' ' + token
            elif current == 's':
                object_ += ' ' + token
            elif current == 'o':
                relation += ' ' + token
    if subject != '' and relation != '' and object_ != '':
        triplets.append({'head': subject.strip(), 'type': relation.strip(),'tail': object_.strip()})
    return triplets

text = "[RE] Serotonin levels are correlated with hippocampus structure in major depressive disorder. The prefrontal cortex and amygdala are also affected in this condition."
triples = []

def generate_triples(texts):
    inputs = tokenizer([text], return_tensors="pt", max_length=512, truncation=True)
    print("inputs", inputs)
    with torch.no_grad():
        outputs = model.generate(
            **inputs,
            max_length=512,
            num_beams=8,
            temperature=0.7,
            num_return_sequences=1,
            early_stopping=True
        )
    print("outputs", outputs)

    decoded = tokenizer.batch_decode(outputs, skip_special_tokens=False)

    print("decoded", decoded)

    print("\n🧠 Raw generated output:\n", decoded)

    for idx, sentence in enumerate(decoded):
        et = extract_triplets(sentence)
        for t in et:
            triples.append((t['head'], t['type'], t['tail']))

    print("triples", triples)
    structured = []
    print("triples", triples)
    for triple in triples:
        parts = re.findall(r'"(.*?)"', triple)
        print("parts", parts)
        if len(parts) == 3:
            structured.append({
                "subject": parts[0],
                "relation": parts[1],
                "object": parts[2]
            })

    print("\n🔗 Extracted Triples:")
    for rel in structured:
        print(f"→ {rel['subject']}  —[{rel['relation']}]->  {rel['object']}")

generate_triples(text)