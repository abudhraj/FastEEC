// FastEEC 0.3 -- This file is part of FastEEC. This code implements
// the computation of N-point projected energy correlators for any (integer
// or non-integer) N, following our new parametrization. Here we compute
// the largest distances R1 with respect to a special particle, which is then
// summed over. This provides a reasonable description of the correlator in the
// perturbative region, without any approximations; though there are some differences
// in the transition to the non-perturbative region. This code achieves a computational 
// scaling of M^2 log(M), for any N value.
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

using namespace std;

//Global variables
const double PI = 3.14159265358979323846;
const double TWOPI = 2 * PI;
const double maxbin = 0;
std::vector<double> res(1000);

// Function to sort the distances
bool sortbydist(const pair<double,double> &a,
               const pair<double,double> &b)
{
    return (a.first < b.first);
}

//Function to convert the distances into corresponding bin positions
int binposcalc(double ang, const double minbin, const int nbins){
    int binpos;
    //Since we take a log, I include this if-statement
    if (ang <= pow(10,minbin))
        binpos=0;
    else
        binpos = std::trunc((log10(ang)-minbin)/(maxbin-minbin)*nbins);
    
    //Handle overflow and underflow
    if (binpos < 0)
        binpos = 0;
    if (binpos >= nbins)
        binpos = nbins-1;
    
    return binpos;
}

int main(int argc, char* argv[]){
    
    //print banner
    printHeader();
    
    //Check syntax
    if (argc != 7){
        cout << "Usage: ./eec_angles input_file events N minbin nbins output_file\n"
        << "\n"
        << "The input_file from which events should be read\n"
        << "events > 0 is the number of events \n"
        << "nu > 0 specifies which point correlator to compute \n"
        << "minbin < 0 specifies the minimum bin value as log10(R_L) needed for sampling output \n"
        << "number of bins > 1\n"
        << "The output_file in which data should be stored\n";
        return 0;
    }
    
    //Convert input parameters to variables
    int nevent = atoi(argv[2]);
    double nu = atof(argv[3]);
    const double minbin = atof(argv[4]);
    const int nbins = atoi(argv[5]);
    string fname = argv[6];
    
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
        cout << "Error: Nu must be positive.\n";
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
    
    read_event(events,inputfile,nevent);
    
    //Close input file
    inputfile.close();
    
    //Initialize array that will contain results (defined in eec_nu_point.h) with zeroes
    for (int l = 0; l < nbins; ++l)
        res[l] = 0;
    
    cout.setstate(ios_base::failbit);
    cout.clear();
    
    // Sample over events
    for (int ievent = 0; ievent < nevent; ++ievent) {
        double total_pt = 0.0; // Reset total_pt for each event
        //Calculate the jet pt
        for (long unsigned int i = 0; i < events[ievent].size(); ++i) {
            total_pt += events[ievent][i].pt();
        }
        // Loop over each particle as the reference particle `s`
        for (long unsigned int s = 0; s < events[ievent].size(); ++s) {
            double E_s = events[ievent][s].pt()/total_pt;
            
            // Create a list of distances and energies
            std::vector< pair<double,double> > ang;
            for (long unsigned int j = 0; j < events[ievent].size(); ++j){
                
                double dphiabs = fabs(events[ievent][s].phi() - events[ievent][j].phi());
                double dphi = dphiabs > PI ? TWOPI - dphiabs : dphiabs;
                double chi = sqrt(pow(fabs(events[ievent][s].rap() - events[ievent][j].rap()), 2) + pow(dphi, 2));
                
                ang.push_back({chi, events[ievent][j].pt()/total_pt});
            }
            // Sort the list by distances
            sort(ang.begin(), ang.end(), sortbydist );
            
            // Handle the contact term contribution (R1 = 0)
            res[0] += pow(E_s, nu);
            
            // Calculate contributions for the sorted list
            double esum = E_s;
            for(long unsigned int k = 1 ; k < ang.size(); ++k){
                double esum_k = esum + ang[k].second;
                double eec = E_s * (pow(esum_k, nu-1) - pow(esum, nu-1));
                
                // Add to the corresponding bin
                res[binposcalc(ang[k].first,minbin,nbins)] += eec;
                esum = esum_k;
            }
        }
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
