import os

# Configuration
directory = "."
chunk_size = 40 * 1024 * 1024  # 40 MB
min_size = 100 * 1024 * 1024   # 100 MB
mapping_file = "file_parts_list.txt"
skip_files = {"split.py", "recombine.py"}  # filenames to skip

directory = os.path.abspath(directory)

with open(mapping_file, "w", encoding="utf-8") as map_file:
    for root, _, files in os.walk(directory):
        if ".git" in os.path.relpath(root, directory).split(os.sep):
            continue
        for f in files:
            if f in skip_files:
                continue
            full_path = os.path.join(root, f)
            size = os.path.getsize(full_path)
            if size > min_size:
                rel_path = os.path.relpath(full_path, directory)
                parts = []
                with open(full_path, "rb") as infile:
                    part_num = 1
                    while True:
                        chunk = infile.read(chunk_size)
                        if not chunk:
                            break
                        part_name = f"{full_path}.part{part_num}"
                        with open(part_name, "wb") as part_file:
                            part_file.write(chunk)
                        parts.append(os.path.relpath(part_name, directory))
                        part_num += 1
                # Record mapping
                map_file.write(f"{rel_path}\t{','.join(parts)}\n")
                print(f"Split {rel_path} into {len(parts)} parts")
                # Delete original file
                os.remove(full_path)
                print(f"Deleted original file: {rel_path}")
