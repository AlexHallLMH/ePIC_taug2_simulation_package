import re
import subprocess
import pyhepmc
from pyhepmc.io import ReaderAscii, WriterAscii

input_filename = "Grape_tauola.hepmc"
temp_filename = "Grape_temp.hepmc"
output_filename = "Grape_final.hepmc"

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

# Step 2: Process the clean file
reader = ReaderAscii(input_filename)
writer = WriterAscii(temp_filename)
event = pyhepmc.GenEvent()

def get_tau_descendants(particle):
    descendants = set()
    stack = list(particle.end_vertex.particles_out if particle.end_vertex else [])
    while stack:
        p = stack.pop()
        if id(p) in descendants:
            continue
        descendants.add(id(p))
        if p.end_vertex:
            stack.extend(p.end_vertex.particles_out)
    return descendants

while not reader.failed():
    reader.read_event(event)
    if reader.failed():
        break

    particles = list(event.particles)

    # Find last tau+ and tau-
    last_tau_plus = None
    last_tau_minus = None
    for p in particles:
        if p.pid == 15:
            last_tau_plus = p
        elif p.pid == -15:
            last_tau_minus = p

    keep_ids = set()
    selected_taus = []

    if last_tau_plus:
        keep_ids.add(id(last_tau_plus))
        selected_taus.append(last_tau_plus)
    if last_tau_minus:
        keep_ids.add(id(last_tau_minus))
        selected_taus.append(last_tau_minus)

    # Add tau descendants
    for tau in selected_taus:
        keep_ids |= get_tau_descendants(tau)

    # Add all final-state particles
    for p in particles:
        if p.status == 1:
            keep_ids.add(id(p))

    # Build cleaned event
    new_event = pyhepmc.GenEvent()
    id_to_new_particle = {}
    vertex_map = {}

    for p in particles:
        if id(p) not in keep_ids:
            continue

        new_p = pyhepmc.GenParticle(p.momentum, p.pid, p.status)
        id_to_new_particle[id(p)] = new_p

        if p.production_vertex:
            old_vtx = p.production_vertex
            if old_vtx not in vertex_map:
                new_vtx = pyhepmc.GenVertex(old_vtx.position)
                vertex_map[old_vtx] = new_vtx
                new_event.add_vertex(new_vtx)
            vertex_map[old_vtx].add_particle_out(new_p)
        else:
            new_event.add_particle(new_p)  # incoming particle

    for old_vtx, new_vtx in vertex_map.items():
        for p_in in old_vtx.particles_in:
            if id(p_in) in id_to_new_particle:
                new_vtx.add_particle_in(id_to_new_particle[id(p_in)])

    writer.write_event(new_event)
    event.clear()

reader.close()
writer.close()

# Step 3: Final clean pass
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
