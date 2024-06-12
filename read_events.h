// FastEEC 0.1. This header file is a part of FastEEC which
// reads the input data file and checks that the file contains 
// the right number of input values for each jet. The inputs for each
// jet are provided in the order : jet_no, p_t, eta and phi of the 
// particles in a jet.

#ifndef READ_EVENT_H
#define READ_EVENT_H

#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include "fastjet/PseudoJet.hh" // Include fastjet header

template<typename Container>
void read_event(Container& events, std::ifstream& inputfile, int nevent) {
    std::string line;
    int event_count = 0;
    while (std::getline(inputfile, line) && event_count < nevent) {
        std::istringstream linestream(line);
        // take substrings to avoid problems when there are extra "pollution"
        // characters (e.g. line-feed).
        if (line.substr(0,4) == "#END") {return;}
        if (line.substr(0,1) == "#") {continue;}
        int event_number;
        double pt, eta, phi;
        linestream >> event_number >> pt >> eta >> phi;
        
        // Check if all variables were successfully read
        if (linestream.fail()) {
            std::cout << "Error: Incorrect input format in line: " << line << std::endl;
            std::exit(0);
        }
        
        double px = pt * cos(phi);
        double py = pt * sin(phi);
        double pz = pt * sinh(eta);
        double E = sqrt(px * px + py * py + pz * pz);
        fastjet::PseudoJet particle(px, py, pz, E);
        if (event_number >= nevent) {
            break;
        } else {
            events[event_number].push_back(particle);
        }
    }
    
    //Check if enough events were provided
    for (int i=0; i<nevent; ++i)
        if (events[i].size()==0){
            std::cout << "Error: Input file does not have enough events." << std::endl;
            std::exit(0);
        }
    
}

#endif // READ_EVENT_H
