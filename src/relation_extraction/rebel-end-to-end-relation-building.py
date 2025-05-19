import json
from transformers import AutoTokenizer, AutoModelForSeq2SeqLM
from tqdm import tqdm
from datetime import datetime
import os


with open("data/parsed/parsed_articles.json") as f:
    data = json.load(f)


def write_file(output_dir, file_name, output_data):
    print(f"Writing to {output_dir}/{file_name}")
    os.makedirs(output_dir, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = f"{output_dir}/{file_name}_{timestamp}.json"

    try:
        with open(output_file, "x") as f:
            json.dump(output_data, f, indent=2)
        print(f"✅ Created new file: {output_file}")
    except FileExistsError:
        print(f"⚠️ File already exists: {output_file}")


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

# Load model and tokenizer
# model_name = "Babelscape/rebel-large"
# model_name = "alexyamato/finetuned-rebel-bio"
model_name = "results/checkpoint-5"

tokenizer = AutoTokenizer.from_pretrained(model_name)
model = AutoModelForSeq2SeqLM.from_pretrained(model_name)
gen_kwargs = {
    "max_length": 256,
    "length_penalty": 0,
    "num_beams": 3,
    "num_return_sequences": 3,
}

output_data = []

for item in tqdm(data, desc="🔗 Predicting Relations"):
    text = item["abstract"]
    model_inputs = tokenizer(text, max_length=256, padding=True, truncation=True, return_tensors = 'pt')

    generated_tokens = model.generate(
        model_inputs["input_ids"].to(model.device),
        attention_mask=model_inputs["attention_mask"].to(model.device),
        **gen_kwargs,
    )

    decoded_preds = tokenizer.batch_decode(generated_tokens, skip_special_tokens=False)

    all_triples = []
    for sentence in decoded_preds:
        triples = extract_triplets(sentence)
        for t in triples:
            if t not in all_triples:
                all_triples.append(t)

    output_data.append({
        "title": item.get("title", ""),
        "abstract": item["abstract"],
        "triples": all_triples
    })

write_file("data/relations", f"relations_output", output_data)