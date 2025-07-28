from pyhepmc.io import ReaderAscii, WriterAscii
import pyhepmc
import math
import re
import io
import sys
from contextlib import redirect_stderr

input_filename = "Grape_fixed_vertices.hepmc"
output_filename = "Grape_done.hepmc"

def delta_p(p1, p2):
    """Euclidean delta-p between two 3-momenta."""
    return math.sqrt((p1.px - p2.px)**2 + (p1.py - p2.py)**2 + (p1.pz - p2.pz)**2)

def parse_error_log(error_log):
    print('here')
    print(error_log)
    """Extract correct vertex, incorrect vertex, and problematic P line."""
    vertex_match = re.search(r'(\d+)\s+vs\s+(\d+)\s+expected', error_log)
    p_line_match = re.search(r'^P\s.*$', error_log, re.MULTILINE)

    if not vertex_match or not p_line_match:
        return None

    correct = int(vertex_match.group(1))
    wrong = int(vertex_match.group(2))
    p_line = p_line_match.group(0).strip()
    return correct, wrong, p_line

def fix_hepmc_vertex(filename, correct_vertex, problematic_p_line):
    """Locate event containing the problematic P line and fix its vertex count."""
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Locate the P line
    p_idx = next((i for i, line in enumerate(lines) if line.strip() == problematic_p_line), None)
    if p_idx is None:
        raise ValueError("Couldn't find P line: " + problematic_p_line)

    # Find the E line just above
    e_idx = next((j for j in range(p_idx, -1, -1) if lines[j].startswith("E ")), None)
    if e_idx is None:
        raise ValueError("No E line found above P line.")

    parts = lines[e_idx].strip().split()
    parts[2] = str(correct_vertex)
    lines[e_idx] = " ".join(parts) + "\n"

    with open(filename, 'w') as f:
        f.writelines(lines)
    print(f"🔧 Fixed vertex number in event starting on line {e_idx + 1}")

def process_events_with_fix(input_filename, output_filename):
    writer = WriterAscii(output_filename)
    event = pyhepmc.GenEvent()
    keep_trying = True

    while keep_trying:
        # Attempt to open and read through the file
        reader = ReaderAscii(input_filename)
        event.clear()
        error_buffer = io.StringIO()

        with redirect_stderr(error_buffer):
            reader.read_event(event)

        if reader.failed():
            error_output = error_buffer.getvalue()
            result = parse_error_log(error_output)
            reader.close()

            if result is None:
                print("❌ Unexpected failure; no parsable error.")
                break

            correct_vertex, wrong_vertex, p_line = result
            fix_hepmc_vertex(input_filename, correct_vertex, p_line)
            continue  # Try again from the top

        # No read failure – process and write event
        beams = [p for p in event.particles if p.status == 3 and abs(p.pid) in [2212, 11]]
        finals = [p for p in event.particles if p.status == 1]
        taus = [p for p in event.particles if abs(p.pid) == 15 and p.status == 1]

        for tau in taus:
            # Find closest beam particle
            closest_beam = min(beams, key=lambda b: delta_p(b.momentum, tau.momentum), default=None)
            if not closest_beam:
                continue

            # Match final state of same PID
            final_candidates = [p for p in finals if p.pid == closest_beam.pid]
            if not final_candidates:
                continue
            final_beam = min(final_candidates, key=lambda p: delta_p(p.momentum, closest_beam.momentum))

            # Infer photon 4-momentum
            px = closest_beam.momentum.px - final_beam.momentum.px
            py = closest_beam.momentum.py - final_beam.momentum.py
            pz = closest_beam.momentum.pz - final_beam.momentum.pz
            e = closest_beam.momentum.e - final_beam.momentum.e

            photon = pyhepmc.GenParticle(pyhepmc.FourVector(px, py, pz, e), 22, 2)

            # Vertex: beam → beam + photon
            vtx1 = pyhepmc.GenVertex()
            vtx1.add_particle_in(closest_beam)
            vtx1.add_particle_out(final_beam)
            vtx1.add_particle_out(photon)
            event.add_vertex(vtx1)

            # Vertex: photon → tau
            vtx2 = pyhepmc.GenVertex()
            vtx2.add_particle_in(photon)
            vtx2.add_particle_out(tau)
            event.add_vertex(vtx2)

        writer.write_event(event)
        event.clear()

        # Continue reading more events
        while not reader.failed():
            error_buffer = io.StringIO()
            with redirect_stderr(error_buffer):
                reader.read_event(event)

            if reader.failed():
                error_output = error_buffer.getvalue()
                result = parse_error_log(error_output)
                reader.close()
                if result:
                    correct_vertex, wrong_vertex, p_line = result
                    fix_hepmc_vertex(input_filename, correct_vertex, p_line)
                    break  # Restart outer loop
                else:
                    print("❌ Unparsable error during loop.")
                    return

            # Process and write current event
            beams = [p for p in event.particles if p.status == 3 and abs(p.pid) in [2212, 11]]
            finals = [p for p in event.particles if p.status == 1]
            taus = [p for p in event.particles if abs(p.pid) == 15 and p.status == 1]

            for tau in taus:
                closest_beam = min(beams, key=lambda b: delta_p(b.momentum, tau.momentum), default=None)
                if not closest_beam:
                    continue

                final_candidates = [p for p in finals if p.pid == closest_beam.pid]
                if not final_candidates:
                    continue
                final_beam = min(final_candidates, key=lambda p: delta_p(p.momentum, closest_beam.momentum))

                px = closest_beam.momentum.px - final_beam.momentum.px
                py = closest_beam.momentum.py - final_beam.momentum.py
                pz = closest_beam.momentum.pz - final_beam.momentum.pz
                e = closest_beam.momentum.e - final_beam.momentum.e

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
        keep_trying = False  # All done

    writer.close()
    print(f"\n✅ Final output written to {output_filename}")

# Run it
if __name__ == "__main__":
    process_events_with_fix(input_filename, output_filename)
