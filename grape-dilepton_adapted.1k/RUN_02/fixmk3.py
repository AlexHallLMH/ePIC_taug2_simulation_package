import re
import subprocess
import pyhepmc
from pyhepmc.io import ReaderAscii, WriterAscii
import math

input_filename = "Grape_dup.hepmc"
temp_filename = "Grape_temp.hepmc"
output_filename = "Grape_done.hepmc"
last_line_fixed = -1  # Tracks how far down we've fixed

# Step 1: Try to read the file using pyhepmc via subprocess
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

# Step 2: Parse error message and get correct vertex + failed line, searching from last_line_fixed+1
def parse_error_output(output, last_line_fixed, filename):
    match = re.search(r"(\d+)\s+vs\s+(\d+)\s+expected", output)
    if not match:
        return None, None, last_line_fixed
    correct_vertex = int(match.group(1))

    failed_line_match = re.search(r"Parsing failed at line:\n(.*)", output)
    if not failed_line_match:
        return None, None, last_line_fixed
    failed_line = failed_line_match.group(1).strip()

    # Now find line number of failed line in file starting from last_line_fixed + 1
    with open(filename, 'r') as f:
        lines = f.readlines()

    for i in range(last_line_fixed + 1, len(lines)):
        if lines[i].strip() == failed_line:
            return correct_vertex, failed_line, i

    print("⚠ Failed line not found beyond last fixed line.")
    return None, None, last_line_fixed

# Step 3: Fix the file in-place
def fix_hepmc_file(filename, failed_line, correct_vertex, last_line_fixed):
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Find the index of the failed line AFTER last_line_fixed
    for i in range(last_line_fixed + 1, len(lines)):
        if lines[i].strip() == failed_line:
            p_index = i
            break
    else:
        print("⚠ Failed line not found in file.")
        return False, last_line_fixed

    # Find nearest event declaration line above
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

# Step 4: Iteratively fix input
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

# Step 5: Modify the now-clean file (add photons)
reader = ReaderAscii(input_filename)
writer = WriterAscii(temp_filename)
event = pyhepmc.GenEvent()

def delta_p(p1, p2):
    return math.sqrt((p1.px - p2.px)**2 + (p1.py - p2.py)**2 + (p1.pz - p2.pz)**2)

while not reader.failed():
    reader.read_event(event)
    if reader.failed():
        break

    beams = [p for p in event.particles if p.status == 3 and abs(p.pid) in [2212, 11]]
    finals = [p for p in event.particles if p.status == 1]
    taus = [p for p in event.particles if abs(p.pid) == 15 and p.status == 1]

    for tau in taus:
        closest_beam = min(beams, key=lambda b: delta_p(b.momentum, tau.momentum))
        final_candidates = [p for p in finals if p.pid == closest_beam.pid]
        if not final_candidates:
            continue
        final_beam = min(final_candidates, key=lambda p: delta_p(p.momentum, closest_beam.momentum))

        px = closest_beam.momentum.px - final_beam.momentum.px
        py = closest_beam.momentum.py - final_beam.momentum.py
        pz = closest_beam.momentum.pz - final_beam.momentum.pz
        e  = closest_beam.momentum.e  - final_beam.momentum.e

        # Skip zero-momentum photons
        if px == 0 and py == 0 and pz == 0 and e == 0:
            continue

        photon = pyhepmc.GenParticle(pyhepmc.FourVector(px, py, pz, e), 22, 2)

        vtx1 = pyhepmc.GenVertex()
        vtx1.add_particle_in(closest_beam)
        vtx1.add_particle_out(final_beam)
        vtx1.add_particle_out(photon)
        event.add_vertex(vtx1)

        vtx2 = pyhepmc.GenVertex()
        vtx2.add_particle_in(photon)
        vtx2.add_particle_out(tau)
        event.add_vertex(vtx2)

    writer.write_event(event)
    event.clear()

reader.close()
writer.close()

# Step 6: FINAL CLEANUP of new file → Grape_done.hepmc
print("🔁 Final re-cleaning pass on output file...")
last_line_fixed = -1
while True:
    success, output = try_read_event_external(temp_filename)
    if success:
        print(f"✅ Final output file is clean. Writing to {output_filename}")
        with open(temp_filename, 'r') as fin, open(output_filename, 'w') as fout:
            fout.writelines(fin.readlines())
        break
    correct_vertex, failed_line, new_last_line_fixed = parse_error_output(output, last_line_fixed, temp_filename)
    if correct_vertex is None:
        print("❌ Couldn't parse final error output.")
        break
    success, last_line_fixed = fix_hepmc_file(temp_filename, failed_line, correct_vertex, last_line_fixed)
    if not success:
        print("❌ Failed to fix final file.")
        break
