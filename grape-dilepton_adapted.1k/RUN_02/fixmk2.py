import re
import subprocess

# The input file to fix (using a duplicate of Grape so that we don't overwrite the original)
input_filename = "Grape_dup.hepmc"

# Step 1: Try to read the file using pyhepmc from a subprocess
def try_read_event_external():
    """Runs reader_tester.py and captures error text even if returncode is 0."""
    result = subprocess.run(
        ["python3", "reader_tester.py", input_filename],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,  # Merge stderr into stdout
        text=True
    )
    output = result.stdout
    print(output)
    has_error = "ERROR::ReaderAscii:" in output
    return not has_error, output

# Step 2: Parse the error message to get the correct vertex count and failed line
def parse_error_output(output):
    """Parses output from reader_tester.py to get correct vertex and failed line"""
    # Look for the vertex mismatch line: e.g. "9  vs  8 expected"
    match = re.search(r"(\d+)\s+vs\s+(\d+)\s+expected", output)
    if not match:
        return None, None
    correct_vertex = int(match.group(1))
    failed_vertex = int(match.group(2))

    # Find the line that failed (P line)
    failed_line_match = re.search(r"Parsing failed at line:\n(.*)", output)
    if not failed_line_match:
        return None, None
    failed_line = failed_line_match.group(1).strip()

    return correct_vertex, failed_line

# Step 3: Fix the .hepmc file in-place
def fix_hepmc_file(filename, failed_line, correct_vertex):
    """Fixes the vertex count in the event containing the failed line"""
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Find the index of the failed P line
    try:
        p_index = lines.index(failed_line + "\n")
    except ValueError:
        print(f"⚠ Failed line not found in file: {failed_line}")
        return False

    # Search backwards for the nearest 'E' line (event declaration)
    for i in range(p_index, -1, -1):
        if lines[i].startswith("E "):
            parts = lines[i].strip().split()
            if len(parts) >= 3:
                print(f"🔧 Fixing event at line {i}: vertex {parts[2]} → {correct_vertex}")
                parts[2] = str(correct_vertex)  # Replace vertex count
                lines[i] = ' '.join(parts) + '\n'
                break
    else:
        print("⚠ Could not find event declaration line for failed P line.")
        return False

    # Overwrite the file with the fixed content
    with open(filename, 'w') as f:
        f.writelines(lines)

    return True

# Main loop: iterate until the file reads cleanly
while True:
    success, output = try_read_event_external()
    if success:
        print("✅ HepMC file is now clean.")
        break  # File is fully fixed

    # Parse the error output
    correct_vertex, failed_line = parse_error_output(output)
    if correct_vertex is None or failed_line is None:
        print("❌ Failed to parse error output. Exiting.")
        print(output)
        break

    # Fix the file
    fix_success = fix_hepmc_file(input_filename, failed_line, correct_vertex)
    if not fix_success:
        print("❌ Failed to fix the HepMC file. Exiting.")
        break

# Step 4: Run your original logic over the now-cleaned file
from pyhepmc.io import ReaderAscii, WriterAscii
import pyhepmc
import math

output_filename = "Grape_done.hepmc"
reader = ReaderAscii(input_filename)
writer = WriterAscii(output_filename)
event = pyhepmc.GenEvent()

def delta_p(p1, p2):
    """Euclidean delta-p between two 3-momenta."""
    return math.sqrt((p1.px - p2.px)**2 + (p1.py - p2.py)**2 + (p1.pz - p2.pz)**2)

while not reader.failed():
    reader.read_event(event)
    if reader.failed():
        print('FFS')
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

# Final cleaning pass to fix any vertex count issues introduced by photon addition
print("🔁 Running final cleanup on output file...")

# Now target the final output file
input_filename = output_filename  # i.e. "Grape_done.hepmc"

while True:
    success, output = try_read_event_external()
    if success:
        print("✅ Final HepMC file is clean after photon addition.")
        break  # Done!

    correct_vertex, failed_line = parse_error_output(output)
    if correct_vertex is None or failed_line is None:
        print("❌ Failed to parse error output from final pass.")
        print(output)
        break

    fix_success = fix_hepmc_file(input_filename, failed_line, correct_vertex)
    if not fix_success:
        print("❌ Final cleanup failed.")
        break
