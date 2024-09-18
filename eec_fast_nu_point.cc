// FastEEC 0.2 -- This file is part of FastEEC. This code implements
// the fast evaluation method for ν-point energy correlators using C/A clustering
// and using a dynamical subjet radius with fixed number of subjets  1 < nsub < 17.
//
// Copyright (C) 2024, A. Budhraja and W. Waalewijn
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib> // for atoi
#include "fastjet/ClusterSequence.hh"
#include "banner.h"
#include "read_events.h"
#include "eec_nu_point.h"

using namespace std;

//Specify jet definition. Radius R is as defined in the eec_nu_point.h header file
fastjet::JetDefinition jetdef(fastjet::genkt_algorithm, R, 0.0,fastjet::pt_scheme);

int main(int argc, char* argv[]){
    
    //print banner
    printHeader();
    
    //Check syntax
    if (argc != 8){
        cout << "Usage: ./eec_fast_nu_point input_file events N nsub minbin nbins output_file\n"
        << "\n"
        << "The input_file from which events should be read\n"
        << "events > 0 is the number of events \n"
        << "N > 0 specifies which point correlator to compute \n"
        << "17 > number of subjets (nsub) > 1\n"
        << "minbin < 0 specifies the minimum bin value as log10(R_L) needed for sampling output \n"
        << "number of bins > 1\n"
        << "The output_file in which data should be stored\n";
        return 0;
    }
    
    //Convert input parameters to variables
    int nevent = atoi(argv[2]);
    nu = atof(argv[3]);
    const double minbin = atof(argv[5]);
    const int nbins = atoi(argv[6]);
    string fname = argv[7];
    nsub = atoi(argv[4]);
    
    ifstream inputfile(argv[1]);
    
    //Perform checks on input parameters
    if (!inputfile){
        cout << "Error: Cannot open the file.\n";
        return 0;
    }
    
    if (nevent <= 0){
        cout << "Error: Number of events must be a positive integer.\n";
        return 0;
    }
    
    if (nu <= 0 ){
        cout << "Error: N must be positive.\n";
        return 0;
    }
    
    if (nsub < 2 || nsub > 16){
        cout << "Error: Number of subjets must be between 2 and 16.\n";
        return 0;
    }
    
    if (minbin >= 0){
        cout << "Error: The minimum bin value cannot be greater than 0.\n";
        return 0;
    }
    
    if (nbins <= 1 || nbins > 1000){
        cout << "Error: Total number of bins must be larger than 1 and smaller than 1000.\n";
        return 0;
    }
    
    
    vector< vector< fastjet::PseudoJet > > events(nevent);
    fastjet::ClusterSequence clustseq;
    vector<fastjet::PseudoJet> jet;
    
    //Read the input file
    read_event(events,inputfile,nevent);
    
    //Close input file
    inputfile.close();
    
    //Initialize the array that will contain results (defined in eec_nu_point.h) with zeroes
    for (int l = 0; l < nbins; ++l)
        res[l] = 0;
    
    //Sample over events
    for (int ievent = 0; ievent < nevent; ++ievent){
        double total_pt = 0.0; // Reset total_pt for each event
        //Calculate the jet pt
        for (long unsigned int m = 0; m < events[ievent].size(); ++m) {
            total_pt += events[ievent][m].pt();
        }
        
        //Using C/A to cluster the whole event into one jet
        cout.setstate(ios_base::failbit);
        
        clustseq = fastjet::ClusterSequence(events[ievent], jetdef);
        jet = clustseq.inclusive_jets(0);
        
        cout.clear();
        
        //Check that everything is clustered into one jet
        if (jet.size()>1){
            cout << "Error: the provided event contains more than one jet. Please provide one jet or increase the value of R in eec_nu_point.h" << endl;
        }
        
        //Calculate the EEC recursively
        EEC(jet[0],total_pt,minbin,nbins);
    }
    
    //Open output file
    ofstream f;
    f.open(fname);
    
    //Print number of events, bins and bin range
    f << nevent << " " << nbins << " " << minbin << " " << maxbin << "\n";
    
    //Output histogram for all paticles
    for (int l = 0; l < nbins; ++l)
        f << res[l] << " ";
    
    //Close output file
    f.close();
    
    return 0;
}
