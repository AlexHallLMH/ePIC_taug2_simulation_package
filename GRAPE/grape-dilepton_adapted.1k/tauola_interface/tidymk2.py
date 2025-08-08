import re
import subprocess
import pyhepmc
from pyhepmc.io import ReaderAscii, WriterAscii
import math

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

# First: run cleaning routine on the input file "Grape_done.hepmc" to ensure no errors
last_line_fixed = -1
while True:
    success, output = try_read_event_external(input_filename)
    if success:
        print("✅ Input file cleaned.")
        break
    correct_vertex, failed_line, new_last_line_fixed = parse_error_output(output, last_line_fixed, input_filename)
    if correct_vertex is None:
        print("❌ Couldn't parse error output during cleaning.")
        break
    success, last_line_fixed = fix_hepmc_file(input_filename, failed_line, correct_vertex, last_line_fixed)
    if not success:
        print("❌ Failed to fix file during cleaning.")
        break

reader = ReaderAscii(input_filename)
writer = WriterAscii(temp_filename)
event = pyhepmc.GenEvent()

def delta_p(p1, p2):
    return math.sqrt((p1.px - p2.px)**2 + (p1.py - p2.py)**2 + (p1.pz - p2.pz)**2)

while not reader.failed():
    reader.read_event(event)
    if reader.failed():
        break

    first_proton = None
    first_electron = None
    last_tau_plus = None
    last_tau_minus = None
    final_state_particles = []

    # Find the first proton and first electron/positron with status 3,
    # collect final state particles with status 1,
    # and track last tau+ and tau- (any status)
    for p in event.particles:
        if p.status == 3:
            if p.pid == 2212 and first_proton is None:
                first_proton = p
            elif abs(p.pid) == 11 and first_electron is None:
                first_electron = p

        if p.status == 1:
            final_state_particles.append(p)

        if p.pid == 15:
            last_tau_plus = p
        elif p.pid == -15:
            last_tau_minus = p

    # Skip event if missing any required particles
    if not (first_proton and first_electron and last_tau_plus and last_tau_minus):
        print("⚠ Skipping event due to missing required particles.")
        event.clear()
        continue

    # Create new event
    new_event = pyhepmc.GenEvent()

    # Map original particles to new particles for parentage handling
    orig_to_new = {}

    # Add first proton and electron to new event
    new_first_proton = pyhepmc.GenParticle(first_proton.momentum, first_proton.pid, first_proton.status)
    new_event.add_particle(new_first_proton)
    orig_to_new[id(first_proton)] = new_first_proton

    new_first_electron = pyhepmc.GenParticle(first_electron.momentum, first_electron.pid, first_electron.status)
    new_event.add_particle(new_first_electron)
    orig_to_new[id(first_electron)] = new_first_electron

    # Add final state particles to new event
    for p in final_state_particles:
        new_p = pyhepmc.GenParticle(p.momentum, p.pid, p.status)
        new_event.add_particle(new_p)
        orig_to_new[id(p)] = new_p

    # Add last tau+ and tau- (could be any status)
    new_last_tau_plus = pyhepmc.GenParticle(last_tau_plus.momentum, last_tau_plus.pid, last_tau_plus.status)
    new_event.add_particle(new_last_tau_plus)
    orig_to_new[id(last_tau_plus)] = new_last_tau_plus

    new_last_tau_minus = pyhepmc.GenParticle(last_tau_minus.momentum, last_tau_minus.pid, last_tau_minus.status)
    new_event.add_particle(new_last_tau_minus)
    orig_to_new[id(last_tau_minus)] = new_last_tau_minus

    # Find which tau is closest to proton and which closest to electron for parentage
    dist_plus_to_proton = delta_p(first_proton.momentum, last_tau_plus.momentum)
    dist_plus_to_electron = delta_p(first_electron.momentum, last_tau_plus.momentum)
    dist_minus_to_proton = delta_p(first_proton.momentum, last_tau_minus.momentum)
    dist_minus_to_electron = delta_p(first_electron.momentum, last_tau_minus.momentum)

    # Assign parents accordingly:
    # Tau closest to proton → proton parent
    # Tau closest to electron → electron parent

    # Determine parent mapping for taus
    if dist_plus_to_proton < dist_plus_to_electron:
        tau_plus_parent = new_first_proton
    else:
        tau_plus_parent = new_first_electron

    if dist_minus_to_proton < dist_minus_to_electron:
        tau_minus_parent = new_first_proton
    else:
        tau_minus_parent = new_first_electron

    # Set up parent-child relations: create vertices for the tau+ and tau- with parents
    # Create vertex for tau+
    vtx_plus = pyhepmc.GenVertex()
    vtx_plus.add_particle_in(tau_plus_parent)
    vtx_plus.add_particle_out(new_last_tau_plus)
    new_event.add_vertex(vtx_plus)

    # Create vertex for tau-
    vtx_minus = pyhepmc.GenVertex()
    vtx_minus.add_particle_in(tau_minus_parent)
    vtx_minus.add_particle_out(new_last_tau_minus)
    new_event.add_vertex(vtx_minus)

    # For final state particles and first proton/electron: preserve original parentage
    # except we already created vertices for the taus

    # Loop over all original particles for parentage reconstruction,
    # skip the ones we already handled (first proton, first electron, last taus)
    for p in event.particles:
        if id(p) not in orig_to_new:
            continue
        new_p = orig_to_new[id(p)]

        # Find parents of p that are also kept in new event
        original_parents = [parent for parent in p.parents if id(parent) in orig_to_new]

        if not original_parents:
            continue

        # Create a new vertex for this parentage
        vtx = pyhepmc.GenVertex()
        for parent in original_parents:
            new_parent = orig_to_new[id(parent)]
            vtx.add_particle_in(new_parent)
        vtx.add_particle_out(new_p)
        new_event.add_vertex(vtx)

    writer.write_event(new_event)
    event.clear()

reader.close()
writer.close()

print("🔁 Final cleaning pass on temp output file...")
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
