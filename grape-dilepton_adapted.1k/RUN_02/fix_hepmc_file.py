from pyhepmc.io import ReaderAscii, WriterAscii  # ← Use updated imports
import pyhepmc
import math

input_filename = "Grape_fixed_vertices.hepmc"
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
        # Find closest beam particle to tau
        closest_beam = min(beams, key=lambda b: delta_p(b.momentum, tau.momentum))

        # Find corresponding final state particle of same type
        final_candidates = [p for p in finals if p.pid == closest_beam.pid]
        if not final_candidates:
            continue
        final_beam = min(final_candidates, key=lambda p: delta_p(p.momentum, closest_beam.momentum))

        # Construct inferred photon 4-momentum
        px = closest_beam.momentum.px - final_beam.momentum.px
        py = closest_beam.momentum.py - final_beam.momentum.py
        pz = closest_beam.momentum.pz - final_beam.momentum.pz
        e  = closest_beam.momentum.e  - final_beam.momentum.e

        photon = pyhepmc.GenParticle(pyhepmc.FourVector(px, py, pz, e), 22, 2)
        # photon.id = ... ← No need to set `.id`, it's read-only in pyhepmc

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

reader.close()
writer.close()
print(f"✅ Fixed HepMC file written to {output_filename}")
