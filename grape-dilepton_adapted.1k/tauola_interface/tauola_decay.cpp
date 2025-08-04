#include "../../TAUOLA.1.1.8-LHC/TAUOLA/include/Tauola/Tauola.h"
#include "../../TAUOLA.1.1.8-LHC/TAUOLA/include/Tauola/TauolaHepMC3Event.h"

#include "../../HepMC3-3.2.6/include/HepMC3/ReaderAscii.h"
#include "../../HepMC3-3.2.6/include/HepMC3/WriterAscii.h"
#include "../../HepMC3-3.2.6/include/HepMC3/Print.h"

#include <iostream>

using namespace std;
using namespace HepMC3;
using namespace Tauolapp;

int main() {
    // Initialize Tauola (must be done before processing any events)
    Tauola::initialize();

    // Open the input HepMC file
    ReaderAscii reader("Grape_done.hepmc");

    // Open the output HepMC file
    WriterAscii writer("Grape_tauola.hepmc");

    GenEvent evt(Units::GEV, Units::CM);

    // Read and process events until EOF or failure
    while (!reader.failed()) {
        reader.read_event(evt);
        if (reader.failed()) break;
        if (evt.particles().empty()) break;

        cout << "BEFORE TAU DECAYS:" << endl;
        Print::listing(evt);

        // Decay taus
        TauolaHepMC3Event t_event(&evt);
        t_event.decayTaus();

        cout << "AFTER TAU DECAYS:" << endl;
        Print::listing(evt);

        // Write modified event to output file
        writer.write_event(evt);

        // Clear event for the next read
        evt.clear();
    }

    // Close files
    reader.close();
    writer.close();

    return 0;
}
