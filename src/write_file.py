import os
from datetime import datetime
import json

def write_file(output_dir, file_name, output_data):
    print(f"Writing to {output_dir}/{file_name}")
    os.makedirs(output_dir, exist_ok=True)

    # Construct the file path
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_file = f"{output_dir}/{file_name}_{timestamp}.json"

    # Write the output
    try:
        with open(output_file, "x") as f:
            json.dump(output_data, f, indent=2)
        print(f"✅ Created new file: {output_file}")
    except FileExistsError:
        print(f"⚠️ File already exists: {output_file}")
