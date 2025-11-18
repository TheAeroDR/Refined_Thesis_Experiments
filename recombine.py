import os

directory = "."
mapping_file = "file_parts_list.txt"
skip_files = {"split.py", "recombine.py"}  # filenames to skip

directory = os.path.abspath(directory)

with open(mapping_file, "r", encoding="utf-8") as map_file:
    for line in map_file:
        line = line.strip()
        if not line or "\t" not in line:
            continue
        rel_original_file, rel_parts_str = line.split("\t")
        original_file = os.path.join(directory, rel_original_file)
        parts = [os.path.join(directory, p.strip()) for p in rel_parts_str.split(",")]

        if os.path.basename(original_file) in skip_files:
            print(f"Skipping {rel_original_file}")
            continue

        os.makedirs(os.path.dirname(original_file), exist_ok=True)
        with open(original_file, "wb") as outfile:
            for part in parts:
                with open(part, "rb") as infile:
                    outfile.write(infile.read())
        print(f"Recombined {rel_original_file} from {len(parts)} parts")

        # Delete split parts
        for part in parts:
            os.remove(part)
        print(f"Deleted {len(parts)} split parts for {rel_original_file}")

# Delete mapping file
os.remove(mapping_file)
print(f"Deleted mapping file: {mapping_file}")
