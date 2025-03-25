from transformers import AutoTokenizer, AutoModelForSeq2SeqLM
import torch
import re

# Load model & tokenizer
# model_name = "Babelscape/rebel-large"
model_name = "dmis-lab/biobert-v1.1"
tokenizer = AutoTokenizer.from_pretrained(model_name)
model = AutoModelForSeq2SeqLM.from_pretrained(model_name)
model.eval()

# Input text
text = "[RE] Serotonin levels are correlated with hippocampus structure in major depressive disorder. The prefrontal cortex and amygdala are also affected in this condition."

# Tokenize
inputs = tokenizer([text], return_tensors="pt", max_length=512, truncation=True)
print("inputs", inputs)
# Generate with beam search
with torch.no_grad():
    outputs = model.generate(
        **inputs,
        max_length=512,
        # num_beams=5,
        num_beams=8,           # More beams for diverse output
        temperature=0.7,
        # num_return_sequences=1,
        early_stopping=True
    )
print("outputs", outputs)
# Decode output
decoded = tokenizer.decode(outputs[0], skip_special_tokens=True)

print("\n🧠 Raw generated output:\n", decoded)

# Extract structured triples using regex
triples = re.findall(r"\(.*?\)", decoded)
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