from transformers import AutoTokenizer, AutoModelForSeq2SeqLM

model = AutoModelForSeq2SeqLM.from_pretrained("Babelscape/rebel-large")
tokenizer = AutoTokenizer.from_pretrained("Babelscape/rebel-large")

input_text = "Langerin, a novel C-type lectin specific to Langerhans cells, is an endocytic receptor that induces the formation of Birbeck granules. We have identified a type II Ca2+-dependent lectin displaying mannose-binding specificity, exclusively expressed by Langerhans cells (LC), and named Langerin. LC are uniquely characterized by Birbeck granules (BG), which are organelles consisting of superimposed and zippered membranes. Here, we have shown that Langerin is constitutively associated with BG and that antibody to Langerin is internalized into these structures. Remarkably, transfection of Langerin cDNA into fibroblasts created a compact network of membrane structures with typical features of BG. Langerin is thus a potent inducer of membrane superimposition and zippering leading to BG formation. Our data suggest that induction of BG is a consequence of the antigen-capture function of Langerin, allowing routing into these organelles and providing access to a nonclassical antigen-processing pathway."

inputs = tokenizer(input_text, return_tensors="pt")
outputs = model.generate(**inputs)
print("base model:", tokenizer.decode(outputs[0]))

model = AutoModelForSeq2SeqLM.from_pretrained("results/checkpoint-5")
tokenizer = AutoTokenizer.from_pretrained("results/checkpoint-5")


inputs = tokenizer(input_text, return_tensors="pt")
outputs = model.generate(**inputs)
print("finetuned model:", tokenizer.decode(outputs[0]))