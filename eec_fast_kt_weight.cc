// FastEEC 0.1 -- This file is part of FastEEC. This code implements
// the fast evaluation method with higher weights for N point correlators 
// using kt clustering and a jet resolution f = f' kt_min^2, with f' >0.
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
#include "eec_higher_weight.h"

using namespace std;

//Specify jet definition. Radius R is as defined in the eec_higher_weight.h header file
fastjet::JetDefinition jetdef(fastjet::genkt_algorithm, R, 1.0,fastjet::pt_scheme);

int main(int argc, char* argv[]){
    
    //print banner
    printHeader();
    
    //Check syntax
    if (argc != 9){
        cout << "Usage: ./eec_fast_kt_weight input_file events N resolutionfactor kappa minbin nbins output_file\n"
        << "\n"
        << "The input_file from which events should be read\n"
        << "events > 0 is the number of events \n"
        << "1 < N < 9 specifies which point correlator to compute \n"
        << "resolution factor (called f' in paper for k_T clustering) > 0\n"
        << "weight kappa > 0\n"
        << "minbin < 0 specifies the minimum bin value as log10(R_L) needed for sampling output \n"
        << "number of bins > 1\n"
        << "The output_file in which data should be stored\n";
        return 0;
    }
    
    //Convert input parameters to variables
    int nevent = atoi(argv[2]);
    fac = atof(argv[4]);
    double weight = atof(argv[5]);
    const double minbin = atof(argv[6]);
    const int nbins = atoi(argv[7]);
    string fname = argv[8];
    int n = atoi(argv[3]);
    
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
    
    if (n<=1 || n > 8){
        cout << "Error: Specify N value greater than 1.\n";
        return 0;
    }

     if (fac <= 0){
        cout << "Error: Resolution factor must be positive.\n";
        return 0;
    }
    
    if (weight < 0){
        cout << "Error: weight kappa must be positive.\n";
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
    
    read_event(events,inputfile,nevent);
    
    //Close input file
    inputfile.close();
    
    //Initialize array that will contain results (defined in eec_higher_weight.h) with zeroes
    for (int i = 0; i < nbins; ++i)
        res[i] = 0;
    
    //Sample over events
    for (int ievent = 0; ievent < nevent; ++ievent){
        double total_pt = 0.0; // Reset total_pt for each event
        //Calculate the jet pt
        for (long unsigned int i = 0; i < events[ievent].size(); ++i) {
            total_pt += events[ievent][i].pt();
        }
        
        //Using kt to cluster the whole event into one jet
        cout.setstate(ios_base::failbit);
        clustseq = fastjet::ClusterSequence(events[ievent], jetdef);
        jet = clustseq.inclusive_jets(0);
        
        cout.clear();
        
        //Check that everything is clustered into one jet
        if (jet.size()>1){
            cout << "Error: the provided event contains more than one jet. Please provide one jet or increase the value of R in eec_higher_weight.h" << endl;
        }
        
        //Calculate the EEC recursively
        EEC(jet[0],total_pt,n,weight,minbin,nbins);
    }
    
    //Open output file
    ofstream f;
    f.open(fname);
    
    //Print number of events, bins and bin range
    f << nevent << " " << nbins << " " << minbin << " " << maxbin << "\n";
    
    //Output histogram for all paticles
    for (int i = 0; i < nbins; ++i)
        f << res[i] << " ";
    
    //Close output file
    f.close();
    
    return 0;
}
