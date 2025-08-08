# reader_tester.py — test-read a HepMC file and trigger known parsing errors
import pyhepmc
from pyhepmc.io import ReaderAscii
import sys

filename = sys.argv[1]
reader = ReaderAscii(filename)
event = pyhepmc.GenEvent()

event_count = 0

while not reader.failed():
    reader.read_event(event)
    if reader.failed():
        print(">>> ERROR detected after reading event", event_count)
        break
    event_count += 1
    event.clear()

reader.close()
