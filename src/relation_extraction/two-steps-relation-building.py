import json
from transformers import AutoTokenizer, AutoModelForSequenceClassification
import torch
from tqdm import tqdm

model_name = "michiyasunaga/BioLinkBERT-large"
tokenizer = AutoTokenizer.from_pretrained(model_name)
model = AutoModelForSequenceClassification.from_pretrained(model_name)
model.eval()

with open("data/ner/ner_output_hybrid.json") as f:
    data = json.load(f)

def predict_relation(context, head_entity, tail_entity):
    input_text = context.replace(head_entity, f"[E1]{head_entity}[/E1]", 1).replace(
        tail_entity, f"[E2]{tail_entity}[/E2]", 1
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

    return pred

output = []

for item in tqdm(data, desc="🔗 Predicting Relations"):
    context = item["abstract"]
    entities = item["entities"]
    relations = []

    for i, ent1 in enumerate(entities):
        for j, ent2 in enumerate(entities):
            if i == j:
                continue
            try:
                relation_id = predict_relation(context, ent1["text"], ent2["text"])
                print(f"✅ Predicted relation: {ent1['text']} → {ent2['text']} ({relation_id})")
                if relation_id != 0:
                    relations.append({
                        "head": ent1["text"],
                        "tail": ent2["text"],
                        "relation_id": relation_id
                    })
                    print(f"✅ Relation added: {ent1['text']} → {ent2['text']} ({relation_id})")
            except Exception as e:
                print(f"❌ Failed for {ent1['text']} → {ent2['text']}: {e}")

    output.append({
        "title": item.get("title"),
        "abstract": context,
        "entities": entities,
        "relations": relations
    })

with open("data/relations/relation_output.json", "w") as f:
    json.dump(output, f, indent=2)

print("✅ Done! Output saved to data/relations/relation_output.json")