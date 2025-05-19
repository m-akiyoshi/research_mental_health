from torch.utils.data import Dataset

class REBELDataset(Dataset):
    def __init__(self, data, tokenizer, max_input_length=512, max_output_length=256):
        self.data = data
        self.tokenizer = tokenizer
        self.max_input_length = max_input_length
        self.max_output_length = max_output_length

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        item = self.data[idx]
        input_enc = self.tokenizer(
            item["input_text"],
            truncation=True,
            padding="max_length",
            max_length=self.max_input_length,
            return_tensors="pt"
        )
        target_enc = self.tokenizer(
            item["target_text"],
            truncation=True,
            padding="max_length",
            max_length=self.max_output_length,
            return_tensors="pt"
        )
        input_enc = {k: v.squeeze() for k, v in input_enc.items()}
        target_enc = {k: v.squeeze() for k, v in target_enc.items()}
        return {
            "input_ids": input_enc["input_ids"],
            "attention_mask": input_enc["attention_mask"],
            "labels": target_enc["input_ids"]
        }