from transformers import AutoTokenizer, AutoModelForSeq2SeqLM, Seq2SeqTrainer, Seq2SeqTrainingArguments, DataCollatorForSeq2Seq, TrainerCallback,PrinterCallback
from datasets import Dataset
import json

model_name = "Babelscape/rebel-large"
tokenizer = AutoTokenizer.from_pretrained(model_name)
model = AutoModelForSeq2SeqLM.from_pretrained(model_name)

with open("data/BioRED/biored_seq2seq_train.json") as f:
    raw_data = json.load(f)

with open("data/BioRED/biored_seq2seq_dev.json") as f:
    dev_data = json.load(f)


dataset = Dataset.from_list(raw_data)
eval_dataset = Dataset.from_list(dev_data)

def preprocess(example):
    model_inputs = tokenizer(
        example["input_text"],
        max_length=512,
        truncation=True,
        padding="max_length"
    )
    labels = tokenizer(
        example["target_text"],
        max_length=256,
        truncation=True,
        padding="max_length"
    )
    model_inputs["labels"] = [
        (label if label != tokenizer.pad_token_id else -100)
        for label in labels["input_ids"]
    ]
    return model_inputs

tokenized = dataset.map(preprocess, remove_columns=["input_text", "target_text"])
tokenized_eval = eval_dataset.map(preprocess, remove_columns=["input_text", "target_text"])
train_dataset = tokenized.shuffle(seed=42)
print("Train dataset size:", len(tokenized))
print("Sample input_ids length:", len(tokenized[0]['input_ids']))
print("Sample labels length:", len(tokenized[0]['labels']))

args = Seq2SeqTrainingArguments(
    overwrite_output_dir=True,
    output_dir="data/relations",
    evaluation_strategy="steps",
    logging_dir='logs',
    logging_strategy="steps",
    logging_steps=50,
    save_steps=200,
    learning_rate=3e-5,
    per_device_train_batch_size=4,
    num_train_epochs=5,
    save_total_limit=1,
    fp16=False,
    do_train=True,
    do_eval=True,
    disable_tqdm=False,
    remove_unused_columns=False,
)

trainer = Seq2SeqTrainer(
    model=model,
    args=args,
    train_dataset=train_dataset,
    tokenizer=tokenizer,
    data_collator=DataCollatorForSeq2Seq(tokenizer, model=model),
    eval_dataset=tokenized_eval,
    callbacks=[PrinterCallback()]
)

print("trainer", trainer)

trainer.train()

print("training done!")
model.save_pretrained("./finetuned-rebel")
tokenizer.save_pretrained("./finetuned-rebel")