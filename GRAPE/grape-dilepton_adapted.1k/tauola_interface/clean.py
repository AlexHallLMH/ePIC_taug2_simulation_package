#Cleans the HEPMC file and saves

import re
import subprocess

input_filename = "Grape_dup.hepmc"
output_filename = "Grape_done.hepmc"
last_line_fixed = -1  # Tracks how far down we've fixed

def try_read_event_external(filename):
    result = subprocess.run(
        ["python3", "reader_tester.py", filename],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True
    )
    output = result.stdout
    print(output)
    has_error = "ERROR::ReaderAscii:" in output
    return not has_error, output

def parse_error_output(output, last_line_fixed, filename):
    match = re.search(r"(\d+)\s+vs\s+(\d+)\s+expected", output)
    if not match:
        return None, None, last_line_fixed
    correct_vertex = int(match.group(1))

    failed_line_match = re.search(r"Parsing failed at line:\n(.*)", output)
    if not failed_line_match:
        return None, None, last_line_fixed
    failed_line = failed_line_match.group(1).strip()

    with open(filename, 'r') as f:
        lines = f.readlines()

    for i in range(last_line_fixed + 1, len(lines)):
        if lines[i].strip() == failed_line:
            return correct_vertex, failed_line, i

    print("⚠ Failed line not found beyond last fixed line.")
    return None, None, last_line_fixed

def fix_hepmc_file(filename, failed_line, correct_vertex, last_line_fixed):
    with open(filename, 'r') as f:
        lines = f.readlines()

    for i in range(last_line_fixed + 1, len(lines)):
        if lines[i].strip() == failed_line:
            p_index = i
            break
    else:
        print("⚠ Failed line not found in file.")
        return False, last_line_fixed

    for j in range(p_index, -1, -1):
        if lines[j].startswith("E "):
            parts = lines[j].strip().split()
            if len(parts) >= 3:
                print(f"🔧 Fixing event at line {j}: vertex {parts[2]} → {correct_vertex}")
                parts[2] = str(correct_vertex)
                lines[j] = ' '.join(parts) + '\n'
                break
    else:
        print("⚠ Could not find event header line.")
        return False, last_line_fixed

    with open(filename, 'w') as f:
        f.writelines(lines)

    return True, p_index

# Initial cleaning
while True:
    success, output = try_read_event_external(input_filename)
    if success:
        print("✅ Initial file cleaned.")
        break
    correct_vertex, failed_line, new_last_line_fixed = parse_error_output(output, last_line_fixed, input_filename)
    if correct_vertex is None:
        print("❌ Couldn't parse error output.")
        break
    success, last_line_fixed = fix_hepmc_file(input_filename, failed_line, correct_vertex, last_line_fixed)
    if not success:
        print("❌ Failed to fix file.")
        break

# Copy cleaned file to output
with open(input_filename, 'r') as fin, open(output_filename, 'w') as fout:
    fout.writelines(fin.readlines())
print(f"📄 Cleaned file written to {output_filename}")
