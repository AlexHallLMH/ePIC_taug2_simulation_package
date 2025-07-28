#include "../../TAUOLA.1.1.8-LHC/TAUOLA/include/Tauola/Tauola.h"
#include "../../TAUOLA.1.1.8-LHC/TAUOLA/include/Tauola/TauolaHepMC3Event.h"

// HepMC3 headers for reading and printing events
#include "../../HepMC3-3.2.6/include/HepMC3/ReaderAscii.h"
#include "../../HepMC3-3.2.6/include/HepMC3/Print.h"

#include <iostream>

using namespace std;
using namespace HepMC3;
using namespace Tauolapp;

int main(void) {
    int NumberOfEvents = 10; // process up to 10 events from the file

    // Initialize Tauola (must be done before processing any events)
    Tauola::initialize();

    // Open the input HepMC file "Grape.hepmc"
    ReaderAscii reader("Grape_done.hepmc");

    for (int iEvent = 0; iEvent < NumberOfEvents; ++iEvent) {
        // Create a fresh GenEvent for each iteration.
        GenEvent evt(Units::GEV, Units::CM);
        
        // Read the next event from the file.
        // If no event is read, break the loop.
        if (!reader.read_event(evt)) {
            cerr << "No more events or error reading event." << endl;
            break;
        }

        cout << "BEFORE TAU DECAYS:" << endl;
        Print::listing(evt);

        // Wrap the event with Tauola's HepMC3 interface and decay the taus.
        TauolaHepMC3Event t_event(&evt);
        t_event.decayTaus();

        cout << "AFTER TAU DECAYS:" << endl;
        Print::listing(evt);
    }

    return 0;
}
