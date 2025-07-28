import re
import subprocess
import pyhepmc
from pyhepmc.io import ReaderAscii, WriterAscii
import math

input_filename = "Grape_dup.hepmc"
temp_filename = "Grape_temp.hepmc"
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

reader = ReaderAscii(input_filename)
writer = WriterAscii(temp_filename)
event = pyhepmc.GenEvent()

def delta_p(p1, p2):
    return math.sqrt((p1.px - p2.px)**2 + (p1.py - p2.py)**2 + (p1.pz - p2.pz)**2)

while not reader.failed():
    reader.read_event(event)
    if reader.failed():
        break

    incoming_proton = None
    incoming_electron = None
    final_proton = None
    final_electron = None
    taus = []

    for p in event.particles:
        if p.status == 3:
            if p.pid == 2212 and incoming_proton is None:
                incoming_proton = p
            elif p.pid == -11 and incoming_electron is None:
                incoming_electron = p
        elif p.status == 1:
            if p.pid == 2212 and final_proton is None:
                final_proton = p
            elif p.pid == -11 and final_electron is None:
                final_electron = p
            elif abs(p.pid) == 15:
                taus.append(p)

    if not (incoming_proton and incoming_electron and final_proton and final_electron) or len(taus) < 2:
        print("⚠ Skipping event due to missing particles.")
        event.clear()
        continue

    used_taus = set()
    new_event = pyhepmc.GenEvent()
    new_taus = []

    new_incoming_proton = pyhepmc.GenParticle(incoming_proton.momentum, incoming_proton.pid, incoming_proton.status)
    new_incoming_electron = pyhepmc.GenParticle(incoming_electron.momentum, incoming_electron.pid, incoming_electron.status)
    new_event.add_particle(new_incoming_proton)
    new_event.add_particle(new_incoming_electron)

    new_final_proton = pyhepmc.GenParticle(final_proton.momentum, final_proton.pid, final_proton.status)
    new_final_electron = pyhepmc.GenParticle(final_electron.momentum, final_electron.pid, final_electron.status)
    new_event.add_particle(new_final_proton)
    new_event.add_particle(new_final_electron)

    for p in taus:
        new_p = pyhepmc.GenParticle(p.momentum, p.pid, p.status)
        new_event.add_particle(new_p)
        new_taus.append((p, new_p))

    # --- Proton side photon + tau vertex ---
    delta = incoming_proton.momentum - final_proton.momentum
    photon = pyhepmc.GenParticle(delta, 22, 2)
    new_event.add_particle(photon)

    tau_candidates = [(orig, new) for (orig, new) in new_taus if id(orig) not in used_taus]
    matched_orig, matched_new = min(tau_candidates, key=lambda pair: delta_p(delta, pair[0].momentum))
    used_taus.add(id(matched_orig))

    vtx1 = pyhepmc.GenVertex()
    vtx1.add_particle_in(new_incoming_proton)
    vtx1.add_particle_out(new_final_proton)
    vtx1.add_particle_out(photon)
    new_event.add_vertex(vtx1)

    vtx2 = pyhepmc.GenVertex()
    vtx2.add_particle_in(photon)
    vtx2.add_particle_out(matched_new)

    # Corrected: p_virtual = p_photon - p_tau
    p_virtual_momentum = photon.momentum - matched_orig.momentum
    virtual_tau = pyhepmc.GenParticle(p_virtual_momentum, -15, 2)  # virtual tau, status=2
    new_event.add_particle(virtual_tau)
    vtx2.add_particle_out(virtual_tau)
    new_event.add_vertex(vtx2)

    # --- Electron side photon + tau vertex ---
    delta = incoming_electron.momentum - final_electron.momentum
    photon2 = pyhepmc.GenParticle(delta, 22, 2)
    new_event.add_particle(photon2)

    tau_candidates = [(orig, new) for (orig, new) in new_taus if id(orig) not in used_taus]
    if tau_candidates:
        matched_orig2, matched_new2 = min(tau_candidates, key=lambda pair: delta_p(delta, pair[0].momentum))
        used_taus.add(id(matched_orig2))

        vtx3 = pyhepmc.GenVertex()
        vtx3.add_particle_in(new_incoming_electron)
        vtx3.add_particle_out(new_final_electron)
        vtx3.add_particle_out(photon2)
        new_event.add_vertex(vtx3)

        vtx4 = pyhepmc.GenVertex()
        vtx4.add_particle_in(photon2)
        vtx4.add_particle_out(matched_new2)

        # Corrected: p_virtual = p_photon - p_tau
        p_virtual_momentum2 = photon2.momentum - matched_orig2.momentum
        virtual_tau2 = pyhepmc.GenParticle(p_virtual_momentum2, -15, 2)
        new_event.add_particle(virtual_tau2)
        vtx4.add_particle_out(virtual_tau2)
        new_event.add_vertex(vtx4)

    writer.write_event(new_event)
    event.clear()

reader.close()
writer.close()

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