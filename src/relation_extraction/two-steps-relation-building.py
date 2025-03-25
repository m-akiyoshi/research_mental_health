import json
from transformers import AutoTokenizer, AutoModelForSequenceClassification
import torch

# Load model
model_name = "michiyasunaga/BioLinkBERT-base"
tokenizer = AutoTokenizer.from_pretrained(model_name)
model = AutoModelForSequenceClassification.from_pretrained(model_name)
model.eval()

# Load your NER JSON
with open("data/ner/ner_output_abstract.json") as f:
    data = json.load(f)

def predict_relation(context, head_entity, tail_entity):
    # Format input (BioLinkBERT expects sentence with marked entities)
    input_text = context.replace(head_entity, f"[E1]{head_entity}[/E1]").replace(
        tail_entity, f"[E2]{tail_entity}[/E2]"
    )
    
    inputs = tokenizer(
        input_text,
        return_tensors="pt",
        max_length=512,
        truncation=True,
        padding="max_length"
    )
    with torch.no_grad():
        logits = model(**inputs).logits
        pred = torch.argmax(logits, dim=1).item()
    
    # You may map pred to label if model provides it
    return pred

# Run on all entity pairs
for item in data:
    context = item["abstract"]
    entities = item["entities"]
    
    for i, ent1 in enumerate(entities):
        for j, ent2 in enumerate(entities):
            if i == j:
                continue
            relation = predict_relation(context, ent1["text"], ent2["text"])
            print(f"{ent1['text']} → {ent2['text']} = Predicted relation ID: {relation}")