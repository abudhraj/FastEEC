// FastEEC 0.2 -- This header file is a part of FastEEC. This code implements
// the fast evaluation method for ν-point energy correlators clustered using a C/A
// clustering algorithm with a fixed number of subjets nsub. This code is an
// extension of the standard N-point projected correlators for integer N values
// provided within FastEEC 0.1. It can be used to compute the ν-point energy correlator
// for any value of ν. The number of subjets that can be utilized for the computation
// is currently limited to nsub = 16.
//
// Copyright (C) 2024, A. Budhraja and W. Waalewijn
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, version 2 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#ifndef EEC_H
#define EEC_H

#include <vector>
#include <cmath>
#include "fastjet/PseudoJet.hh"

const double PI = 3.14159265358979323846;
const double TWOPI = 2 * PI;

//Global variables that can be accessed directly in the main code eec_fast_nu_point.cc
const double maxbin = 0;
std::vector<double> res(1000);
double nu;
int nsub;
//We have chosen a jet radius R that assumes that the provided input is clustered into one jet
const double R = 1.5;
const int ppow2[17] = {1,2,4,8,16,32,64,128,256,512,1024,2048, 4096, 8192, 16384, 32768, 65536};

//Function to compute distance in (eta,phi)
double ang(const fastjet::PseudoJet& v1, const fastjet::PseudoJet& v2) {
    double dphiabs = std::fabs(v1.phi() - v2.phi());
    double dphi = dphiabs > PI ? TWOPI - dphiabs : dphiabs;
    double dy = v1.rap() - v2.rap();
    return std::sqrt(dy*dy + dphi*dphi);
}

//Function to compute the bin position
int binposcalc(double chi, const double minbin, const int nbins){
    int binpos;
    //Since we take a log, we include this if-statement
    if (chi <= pow(10,minbin))
        binpos=0;
    else
        binpos = std::trunc((log10(chi)-minbin)/(maxbin-minbin)*nbins);
    
    //Handle overflow and underflow
    if (binpos < 0)
        binpos = 0;
    if (binpos >= nbins)
        binpos = nbins-1;
    
    return binpos;
}

class EEC {
public:
    EEC(fastjet::PseudoJet tree, double jpt, const double minbin, const int nbins) {
        compute(tree, jpt, minbin, nbins);
    }
    
    void compute(const fastjet::PseudoJet& tree, double jpt, const double minbin, const int nbins)
    {
        // Pass by reference
        fastjet::PseudoJet p1, p2;
        
        if (tree.has_parents(p1, p2))
        {
            {
                // Variables used in recursion that require space. Here
                //   partdist stores the distance between all possible pairs of pesudojets, and
                //   dist is for storing the maximum separation for each subset of particles in the recursion, 
                //      labelling the subset by the corresponding binary number.
                //   w stores the weights for each subset, separating the contribution by the number of particles (to use in recursion)
                //   ptfrac stores the momentum fractions of the pseudojets, 
                
                // As dist and w require large memory, we dynamically allocate the variables based on actual number of subjets created
                int ppow2npt = std::pow(2,nsub);
                std::vector<unsigned long> i(nsub);
                std::vector<int> partdist(nsub * nsub);
                std::vector<int> dist(ppow2npt);
                std::vector<double> w(ppow2npt * nsub);
                std::vector<double> ptfrac(nsub);
                
                std::vector<fastjet::PseudoJet> b1, b2;
                
                //Resolve as much of each of the branches as specified by nsub
                b1 = p1.exclusive_subjets_up_to(nsub/2);
                b2 = p2.exclusive_subjets_up_to(nsub/2);
                
                //Determine how many particles are in each branch
                int nb1 = b1.size(), nb2 = b2.size();
                
                //Total number of subjets should be at most 16.
                if(nb1+nb2>16){                
                    std::cout << "The number subjets must be between 2 and 16" << std::endl;
                    return;
                }
                
                //Append b2 to b1 for convenience
                b1.insert(b1.end(), b2.begin(), b2.end());
                
                //Initialize arrays with zeroes
                for (int j = 0; j<ppow2npt; ++j)
                    for (int k=0; k<nsub; ++k)
                        w[j * nsub + k]=0;                
                for (int j = 0; j<ppow2npt; ++j)
                    dist[j]=0;
                
                //Calculate the momentum fractions
                for (unsigned long j =0; j<b1.size(); ++j)
                    ptfrac[j] = b1[j].pt()/jpt ;
                
                //Precalculate bin position for all angles between pair of pseudojets
                for (unsigned long j=0; j<b1.size(); ++j)
                    for (unsigned long k=0; k<b1.size(); ++k)
                        partdist[j * nsub + k] = binposcalc(ang(b1[j],b1[k]),minbin,nbins);
                
                //Calculate w and dist
                
                //Calculate all 1 particle weights, i.e.
                //       z_1^nu, z_2^nu, ... , z_n^nu.
                // for which the distance is 0.
                for (i[0]=0; i[0]<b1.size(); ++i[0]){                    
                    w[ppow2[i[0]] * nsub + 0] = std::pow(ptfrac[i[0]],nu);
                    //The binary number corresponding to the single particle i[0] is obtained using powers of 2 (ppow2).
                    dist[ppow2[i[0]]] = 0;
                }
                
                //Calculate all possible 2 particle weights and subtract the corresponding 1 particle weights from this, i.e. 
                //       (z_1 + z_2)^nu - z_1^nu - z_2^nu, (z_1 + z_3)^nu - z_1^nu - z_3^nu, ...
                //For n>=2, we also account for the appropriate combinatorial factors in the weight, to avoid double counting.
                //This is kept track of by separating the contribution to the weight by the number of particles involved.
                //Finally, the separation between the particles is also determined. 
                if (b1.size()>=2){
                    int step = 2;
                    for (i[0]=0; i[0]<b1.size()-1; ++i[0])
                        for (i[1]=i[0]+1; i[1]<b1.size(); ++i[1]){
                            //Calculate the binary number corresponding to the subset with particles i[0] and i[1] and their total energy. 
                            int ind = 0;
                            double en = 0;
                            for (int j=0; j<step; ++j){
                                ind += ppow2[i[j]];
                                en += ptfrac[i[j]];
                            }
                            
                            //Calculate the direct (2-particle) contribution to the weight
                            w[ind * nsub + step-1]=std::pow(en,nu);
                            
                            //Include lower-point contributions recursively.
                            for (int j=0; j<step; ++j)
                                for (int k=0; k<step-1; ++k)
                                    w[ind * nsub +k] -= w[(ind - ppow2[i[j]]) * nsub + k]/(step-1-k);
                            
                            //Determine the distance using the precalculate partdist
                            dist[ind] = partdist[i[0] * nsub + i[1]];
                        }
                }
                
                //Calculate all possible 3 particle weights and subtract the corresponding lower-points weights from this, i.e.
                //     (z_1 + z_2 + z_3)^nu - (z_1 + z_2)^nu - (z_1 + z_3)^nu - (z_2 + z_3)^nu + z_1^nu + z_2^nu + z_3^nu , ....
                //Finally, the maximum separation between the particles is also determined.
                if (b1.size()>=3){
                    int step = 3;
                    for (i[0]=0; i[0]<b1.size()-2; ++i[0])
                        for (i[1]=i[0]+1; i[1]<b1.size()-1; ++i[1])
                            for (i[2]=i[1]+1; i[2]<b1.size(); ++i[2]){
                                //Calculate the binary number corresponding to the subset with particles i[0], i[1] and i[2] and their total energy.
                                int ind = 0;
                                double en = 0;
                                for (int j=0; j<step; ++j){
                                    ind+=ppow2[i[j]];
                                    en+=ptfrac[i[j]];
                                }
                                
                                //Calculate the direct (3-particle) contribution.
                                w[ind * nsub + step-1]=std::pow(en,nu);
                                
                                //Include lower-point contributions recusrively.
                                for (int j=0; j<step; ++j)
                                    for (int k=0; k<step-1; ++k)
                                        w[ind * nsub + k] -= w[(ind - ppow2[i[j]]) * nsub + k]/(step-1-k);
                                
                                //Calculate the maximum distance recursively from the subset without i[0].
                                dist[ind] = dist[ind-ppow2[i[0]]];
                                for (int j=1; j<step; ++j)
                                    dist[ind] = std::max(dist[ind], partdist[i[0] * nsub + i[j]]);
                            }
                }
                
                //4 particles
                if (b1.size()>=4){
                    int step = 4;
                    for (i[0]=0; i[0]<b1.size()-3; ++i[0])
                        for (i[1]=i[0]+1; i[1]<b1.size()-2; ++i[1])
                            for (i[2]=i[1]+1; i[2]<b1.size()-1; ++i[2])
                                for (i[3]=i[2]+1; i[3]<b1.size(); ++i[3]){
                                    //Calculate the index
                                    int ind = 0;
                                    double en = 0;
                                    for (int j=0; j<step; ++j){
                                        ind+=ppow2[i[j]];
                                        en+=ptfrac[i[j]];
                                    }
                                    
                                    //Calculate the direct contribution
                                    w[ind * nsub + step-1]=std::pow(en,nu);
                                    
                                    //Include lower-point contributions
                                    for (int j=0; j<step; ++j)
                                        for (int k=0; k<step-1; ++k)
                                            w[ind * nsub + k] -= w[(ind - ppow2[i[j]]) * nsub + k]/(step-1-k);
                                    
                                    //Distance
                                    dist[ind] = dist[ind-ppow2[i[0]]];
                                    for (int j=1; j<step; ++j)
                                        dist[ind] = std::max(dist[ind], partdist[i[0] * nsub + i[j]]);
                                }
                }
                
                //5 particles
                if (b1.size()>=5){
                    int step = 5;
                    for (i[0]=0; i[0]<b1.size()-4; ++i[0])
                        for (i[1]=i[0]+1; i[1]<b1.size()-3; ++i[1])
                            for (i[2]=i[1]+1; i[2]<b1.size()-2; ++i[2])
                                for (i[3]=i[2]+1; i[3]<b1.size()-1; ++i[3])
                                    for (i[4]=i[3]+1; i[4]<b1.size(); ++i[4]){
                                        //Calculate the index
                                        int ind = 0;
                                        double en = 0;
                                        for (int j=0; j<step; ++j){
                                            ind+=ppow2[i[j]];
                                            en+=ptfrac[i[j]];
                                        }
                                        
                                        //Calculate the direct contribution
                                        w[ind * nsub + step-1]=std::pow(en,nu);
                                        
                                        //Include lower-point contributions
                                        for (int j=0; j<step; ++j)
                                            for (int k=0; k<step-1; ++k)
                                                w[ind * nsub + k] -= w[(ind - ppow2[i[j]]) * nsub + k]/(step-1-k);
                                        
                                        //Distance
                                        dist[ind] = dist[ind-ppow2[i[0]]];
                                        for (int j=1; j<step; ++j)
                                            dist[ind] = std::max(dist[ind], partdist[i[0] * nsub + i[j]]);
                                    }
                }
                
                
                //6 particles
                if (b1.size()>=6){
                    int step = 6;
                    for (i[0]=0; i[0]<b1.size()-5; ++i[0])
                        for (i[1]=i[0]+1; i[1]<b1.size()-4; ++i[1])
                            for (i[2]=i[1]+1; i[2]<b1.size()-3; ++i[2])
                                for (i[3]=i[2]+1; i[3]<b1.size()-2; ++i[3])
                                    for (i[4]=i[3]+1; i[4]<b1.size()-1; ++i[4])
                                        for (i[5]=i[4]+1; i[5]<b1.size(); ++i[5]){
                                            //Calculate the index
                                            int ind = 0;
                                            double en = 0;
                                            for (int j=0; j<step; ++j){
                                                ind+=ppow2[i[j]];
                                                en+=ptfrac[i[j]];
                                            }
                                            
                                            //Calculate the direct contribution
                                            w[ind * nsub + step-1]=std::pow(en,nu);
                                            
                                            //Include lower-point contributions
                                            for (int j=0; j<step; ++j)
                                                for (int k=0; k<step-1; ++k)
                                                    w[ind * nsub + k] -= w[(ind - ppow2[i[j]]) * nsub + k]/(step-1-k);
                                            
                                            //Distance
                                            dist[ind] = dist[ind-ppow2[i[0]]];
                                            for (int j=1; j<step; ++j)
                                                dist[ind] = std::max(dist[ind], partdist[i[0] * nsub + i[j]]);
                                        }
                }
                
                //7 particles
                if (b1.size()>=7){
                    int step = 7;
                    for (i[0]=0; i[0]<b1.size()-6; ++i[0])
                        for (i[1]=i[0]+1; i[1]<b1.size()-5; ++i[1])
                            for (i[2]=i[1]+1; i[2]<b1.size()-4; ++i[2])
                                for (i[3]=i[2]+1; i[3]<b1.size()-3; ++i[3])
                                    for (i[4]=i[3]+1; i[4]<b1.size()-2; ++i[4])
                                        for (i[5]=i[4]+1; i[5]<b1.size()-1; ++i[5])
                                            for (i[6]=i[5]+1; i[6]<b1.size(); ++i[6]){
                                                //Calculate the index
                                                int ind = 0;
                                                double en = 0;
                                                for (int j=0; j<step; ++j){
                                                    ind+=ppow2[i[j]];
                                                    en+=ptfrac[i[j]];
                                                }
                                                
                                                //Calculate the direct contribution
                                                w[ind * nsub + step-1]=std::pow(en,nu);
                                                
                                                //Include lower-point contributions
                                                for (int j=0; j<step; ++j)
                                                    for (int k=0; k<step-1; ++k)
                                                        w[ind * nsub + k] -= w[(ind - ppow2[i[j]]) * nsub + k]/(step-1-k);
                                                
                                                //Distance
                                                dist[ind] = dist[ind-ppow2[i[0]]];
                                                for (int j=1; j<step; ++j)
                                                    dist[ind] = std::max(dist[ind], partdist[i[0] * nsub + i[j]]);
                                            }
                }
                
                //8 particles
                if (b1.size()>=8){
                    int step = 8;
                    for (i[0]=0; i[0]<b1.size()-7; ++i[0])
                        for (i[1]=i[0]+1; i[1]<b1.size()-6; ++i[1])
                            for (i[2]=i[1]+1; i[2]<b1.size()-5; ++i[2])
                                for (i[3]=i[2]+1; i[3]<b1.size()-4; ++i[3])
                                    for (i[4]=i[3]+1; i[4]<b1.size()-3; ++i[4])
                                        for (i[5]=i[4]+1; i[5]<b1.size()-2; ++i[5])
                                            for (i[6]=i[5]+1; i[6]<b1.size()-1; ++i[6])
                                                for (i[7]=i[6]+1; i[7]<b1.size(); ++i[7]){
                                                    //Calculate the index
                                                    int ind = 0;
                                                    double en = 0;
                                                    for (int j=0; j<step; ++j){
                                                        ind+=ppow2[i[j]];
                                                        en+=ptfrac[i[j]];
                                                    }
                                                    
                                                    //Calculate the direct contribution
                                                    w[ind * nsub + step-1]=std::pow(en,nu);
                                                    
                                                    //Include lower-point contributions
                                                    for (int j=0; j<step; ++j)
                                                        for (int k=0; k<step-1; ++k)
                                                            w[ind * nsub + k] -= w[(ind - ppow2[i[j]]) * nsub + k]/(step-1-k);
                                                    
                                                    //Distance
                                                    dist[ind] = dist[ind-ppow2[i[0]]];
                                                    for (int j=1; j<step; ++j)
                                                        dist[ind] = std::max(dist[ind], partdist[i[0] * nsub + i[j]]);
                                                }
                }
                
                //9 particles
                if (b1.size()>=9){
                    int step = 9;
                    for (i[0]=0; i[0]<b1.size()-8; ++i[0])
                        for (i[1]=i[0]+1; i[1]<b1.size()-7; ++i[1])
                            for (i[2]=i[1]+1; i[2]<b1.size()-6; ++i[2])
                                for (i[3]=i[2]+1; i[3]<b1.size()-5; ++i[3])
                                    for (i[4]=i[3]+1; i[4]<b1.size()-4; ++i[4])
                                        for (i[5]=i[4]+1; i[5]<b1.size()-3; ++i[5])
                                            for (i[6]=i[5]+1; i[6]<b1.size()-2; ++i[6])
                                                for (i[7]=i[6]+1; i[7]<b1.size()-1; ++i[7])
                                                    for (i[8]=i[7]+1; i[8]<b1.size(); ++i[8]){
                                                        //Calculate the index
                                                        int ind = 0;
                                                        double en = 0;
                                                        for (int j=0; j<step; ++j){
                                                            ind+=ppow2[i[j]];
                                                            en+=ptfrac[i[j]];
                                                        }
                                                        
                                                        //Calculate the direct contribution
                                                        w[ind * nsub + step-1]=std::pow(en,nu);
                                                        
                                                        //Include lower-point contributions
                                                        for (int j=0; j<step; ++j)
                                                            for (int k=0; k<step-1; ++k)
                                                                w[ind * nsub + k] -= w[(ind - ppow2[i[j]]) * nsub + k]/(step-1-k);
                                                        
                                                        //Distance
                                                        dist[ind] = dist[ind-ppow2[i[0]]];
                                                        for (int j=1; j<step; ++j)
                                                            dist[ind] = std::max(dist[ind], partdist[i[0] * nsub + i[j]]);
                                                    }
                }
                
                //10 particles
                if (b1.size()>=10){
                    int step = 10;
                    for (i[0]=0; i[0]<b1.size()-9; ++i[0])
                        for (i[1]=i[0]+1; i[1]<b1.size()-8; ++i[1])
                            for (i[2]=i[1]+1; i[2]<b1.size()-7; ++i[2])
                                for (i[3]=i[2]+1; i[3]<b1.size()-6; ++i[3])
                                    for (i[4]=i[3]+1; i[4]<b1.size()-5; ++i[4])
                                        for (i[5]=i[4]+1; i[5]<b1.size()-4; ++i[5])
                                            for (i[6]=i[5]+1; i[6]<b1.size()-3; ++i[6])
                                                for (i[7]=i[6]+1; i[7]<b1.size()-2; ++i[7])
                                                    for (i[8]=i[7]+1; i[8]<b1.size()-1; ++i[8])
                                                        for (i[9]=i[8]+1; i[9]<b1.size(); ++i[9]){
                                                            
                                                            //Calculate the index
                                                            int ind = 0;
                                                            double en = 0;
                                                            for (int j=0; j<step; ++j){
                                                                ind+=ppow2[i[j]];
                                                                en+=ptfrac[i[j]];
                                                            }
                                                            
                                                            //Calculate the direct contribution
                                                            w[ind * nsub + step-1]=std::pow(en,nu);
                                                            
                                                            //Include lower-point contributions
                                                            for (int j=0; j<step; ++j)
                                                                for (int k=0; k<step-1; ++k)
                                                                    w[ind * nsub + k] -= w[(ind - ppow2[i[j]]) * nsub + k]/(step-1-k);
                                                            
                                                            //Distance
                                                            dist[ind] = dist[ind-ppow2[i[0]]];
                                                            for (int j=1; j<step; ++j)
                                                                dist[ind] = std::max(dist[ind], partdist[i[0] * nsub + i[j]]);
                                                        }
                }
                
                
                //11 particles
                if (b1.size()>=11){
                    int step = 11;
                    for (i[0]=0; i[0]<b1.size()-10; ++i[0])
                        for (i[1]=i[0]+1; i[1]<b1.size()-9; ++i[1])
                            for (i[2]=i[1]+1; i[2]<b1.size()-8; ++i[2])
                                for (i[3]=i[2]+1; i[3]<b1.size()-7; ++i[3])
                                    for (i[4]=i[3]+1; i[4]<b1.size()-6; ++i[4])
                                        for (i[5]=i[4]+1; i[5]<b1.size()-5; ++i[5])
                                            for (i[6]=i[5]+1; i[6]<b1.size()-4; ++i[6])
                                                for (i[7]=i[6]+1; i[7]<b1.size()-3; ++i[7])
                                                    for (i[8]=i[7]+1; i[8]<b1.size()-2; ++i[8])
                                                        for (i[9]=i[8]+1; i[9]<b1.size()-1; ++i[9])
                                                            for (i[10]=i[9]+1; i[10]<b1.size(); ++i[10]){
                                                                
                                                                //Calculate the index
                                                                int ind = 0;
                                                                double en = 0;
                                                                for (int j=0; j<step; ++j){
                                                                    ind+=ppow2[i[j]];
                                                                    en+=ptfrac[i[j]];
                                                                }
                                                                
                                                                //Calculate the direct contribution
                                                                w[ind * nsub + step-1]=std::pow(en,nu);
                                                                
                                                                //Include lower-point contributions
                                                                for (int j=0; j<step; ++j)
                                                                    for (int k=0; k<step-1; ++k)
                                                                        w[ind * nsub + k] -= w[(ind - ppow2[i[j]]) * nsub + k]/(step-1-k);
                                                                
                                                                //Distance
                                                                dist[ind] = dist[ind-ppow2[i[0]]];
                                                                for (int j=1; j<step; ++j)
                                                                    dist[ind] = std::max(dist[ind], partdist[i[0] * nsub + i[j]]);
                                                            }
                }
                
                //12 particles
                if (b1.size()>=12){
                    int step = 12;
                    for (i[0]=0; i[0]<b1.size()-11; ++i[0])
                        for (i[1]=i[0]+1; i[1]<b1.size()-10; ++i[1])
                            for (i[2]=i[1]+1; i[2]<b1.size()-9; ++i[2])
                                for (i[3]=i[2]+1; i[3]<b1.size()-8; ++i[3])
                                    for (i[4]=i[3]+1; i[4]<b1.size()-7; ++i[4])
                                        for (i[5]=i[4]+1; i[5]<b1.size()-6; ++i[5])
                                            for (i[6]=i[5]+1; i[6]<b1.size()-5; ++i[6])
                                                for (i[7]=i[6]+1; i[7]<b1.size()-4; ++i[7])
                                                    for (i[8]=i[7]+1; i[8]<b1.size()-3; ++i[8])
                                                        for (i[9]=i[8]+1; i[9]<b1.size()-2; ++i[9])
                                                            for (i[10]=i[9]+1; i[10]<b1.size()-1; ++i[10])
                                                                for (i[11]=i[10]+1; i[11]<b1.size(); ++i[11]){
                                                                    
                                                                    //Calculate the index
                                                                    int ind = 0;
                                                                    double en = 0;
                                                                    for (int j=0; j<step; ++j){
                                                                        ind+=ppow2[i[j]];
                                                                        en+=ptfrac[i[j]];
                                                                    }
                                                                    
                                                                    //Calculate the direct contribution
                                                                    w[ind * nsub + step-1]=std::pow(en,nu);
                                                                    
                                                                    //Include lower-point contributions
                                                                    for (int j=0; j<step; ++j)
                                                                        for (int k=0; k<step-1; ++k)
                                                                            w[ind * nsub + k] -= w[(ind - ppow2[i[j]]) * nsub + k]/(step-1-k);
                                                                    
                                                                    //Distance
                                                                    dist[ind] = dist[ind-ppow2[i[0]]];
                                                                    for (int j=1; j<step; ++j)
                                                                        dist[ind] = std::max(dist[ind], partdist[i[0] * nsub + i[j]]);
                                                                }
                }
                
                //13 particles
                if (b1.size()>=13){
                    int step = 13;
                    for (i[0]=0; i[0]<b1.size()-12; ++i[0])
                        for (i[1]=i[0]+1; i[1]<b1.size()-11; ++i[1])
                            for (i[2]=i[1]+1; i[2]<b1.size()-10; ++i[2])
                                for (i[3]=i[2]+1; i[3]<b1.size()-9; ++i[3])
                                    for (i[4]=i[3]+1; i[4]<b1.size()-8; ++i[4])
                                        for (i[5]=i[4]+1; i[5]<b1.size()-7; ++i[5])
                                            for (i[6]=i[5]+1; i[6]<b1.size()-6; ++i[6])
                                                for (i[7]=i[6]+1; i[7]<b1.size()-5; ++i[7])
                                                    for (i[8]=i[7]+1; i[8]<b1.size()-4; ++i[8])
                                                        for (i[9]=i[8]+1; i[9]<b1.size()-3; ++i[9])
                                                            for (i[10]=i[9]+1; i[10]<b1.size()-2; ++i[10])
                                                                for (i[11]=i[10]+1; i[11]<b1.size()-1; ++i[11])
                                                                    for (i[12]=i[11]+1; i[12]<b1.size(); ++i[12]){
                                                                        
                                                                        //Calculate the index
                                                                        int ind = 0;
                                                                        double en = 0;
                                                                        for (int j=0; j<step; ++j){
                                                                            ind+=ppow2[i[j]];
                                                                            en+=ptfrac[i[j]];
                                                                        }
                                                                        
                                                                        //Calculate the direct contribution
                                                                        w[ind * nsub + step-1]=std::pow(en,nu);
                                                                        
                                                                        //Include lower-point contributions
                                                                        for (int j=0; j<step; ++j)
                                                                            for (int k=0; k<step-1; ++k)
                                                                                w[ind * nsub + k] -= w[(ind - ppow2[i[j]]) * nsub + k]/(step-1-k);
                                                                        
                                                                        //Distance
                                                                        dist[ind] = dist[ind-ppow2[i[0]]];
                                                                        for (int j=1; j<step; ++j)
                                                                            dist[ind] = std::max(dist[ind], partdist[i[0] * nsub + i[j]]);
                                                                    }
                }
                
                //14 particles
                if (b1.size()>=14){
                    int step = 14;
                    for (i[0]=0; i[0]<b1.size()-13; ++i[0])
                        for (i[1]=i[0]+1; i[1]<b1.size()-12; ++i[1])
                            for (i[2]=i[1]+1; i[2]<b1.size()-11; ++i[2])
                                for (i[3]=i[2]+1; i[3]<b1.size()-10; ++i[3])
                                    for (i[4]=i[3]+1; i[4]<b1.size()-9; ++i[4])
                                        for (i[5]=i[4]+1; i[5]<b1.size()-8; ++i[5])
                                            for (i[6]=i[5]+1; i[6]<b1.size()-7; ++i[6])
                                                for (i[7]=i[6]+1; i[7]<b1.size()-6; ++i[7])
                                                    for (i[8]=i[7]+1; i[8]<b1.size()-5; ++i[8])
                                                        for (i[9]=i[8]+1; i[9]<b1.size()-4; ++i[9])
                                                            for (i[10]=i[9]+1; i[10]<b1.size()-3; ++i[10])
                                                                for (i[11]=i[10]+1; i[11]<b1.size()-2; ++i[11])
                                                                    for (i[12]=i[11]+1; i[12]<b1.size()-1; ++i[12])
                                                                        for (i[13]=i[12]+1; i[13]<b1.size(); ++i[13]){
                                                                            
                                                                            //Calculate the index
                                                                            int ind = 0;
                                                                            double en = 0;
                                                                            for (int j=0; j<step; ++j){
                                                                                ind+=ppow2[i[j]];
                                                                                en+=ptfrac[i[j]];
                                                                            }
                                                                            
                                                                            //Calculate the direct contribution
                                                                            w[ind * nsub + step-1]=std::pow(en,nu);
                                                                            
                                                                            //Include lower-point contributions
                                                                            for (int j=0; j<step; ++j)
                                                                                for (int k=0; k<step-1; ++k)
                                                                                    w[ind * nsub + k] -= w[(ind - ppow2[i[j]]) * nsub + k]/(step-1-k);
                                                                            
                                                                            //Distance
                                                                            dist[ind] = dist[ind-ppow2[i[0]]];
                                                                            for (int j=1; j<step; ++j)
                                                                                dist[ind] = std::max(dist[ind], partdist[i[0] * nsub + i[j]]);
                                                                        }
                }
                
                //15 particles
                if (b1.size()>=15){
                    int step = 15;
                    for (i[0]=0; i[0]<b1.size()-14; ++i[0])
                        for (i[1]=i[0]+1; i[1]<b1.size()-13; ++i[1])
                            for (i[2]=i[1]+1; i[2]<b1.size()-12; ++i[2])
                                for (i[3]=i[2]+1; i[3]<b1.size()-11; ++i[3])
                                    for (i[4]=i[3]+1; i[4]<b1.size()-10; ++i[4])
                                        for (i[5]=i[4]+1; i[5]<b1.size()-9; ++i[5])
                                            for (i[6]=i[5]+1; i[6]<b1.size()-8; ++i[6])
                                                for (i[7]=i[6]+1; i[7]<b1.size()-7; ++i[7])
                                                    for (i[8]=i[7]+1; i[8]<b1.size()-6; ++i[8])
                                                        for (i[9]=i[8]+1; i[9]<b1.size()-5; ++i[9])
                                                            for (i[10]=i[9]+1; i[10]<b1.size()-4; ++i[10])
                                                                for (i[11]=i[10]+1; i[11]<b1.size()-3; ++i[11])
                                                                    for (i[12]=i[11]+1; i[12]<b1.size()-2; ++i[12])
                                                                        for (i[13]=i[12]+1; i[13]<b1.size()-1; ++i[13])
                                                                            for (i[14]=i[13]+1; i[14]<b1.size(); ++i[14]){
                                                                                
                                                                                //Calculate the index
                                                                                int ind = 0;
                                                                                double en = 0;
                                                                                for (int j=0; j<step; ++j){
                                                                                    ind+=ppow2[i[j]];
                                                                                    en+=ptfrac[i[j]];
                                                                                }
                                                                                
                                                                                //Calculate the direct contribution
                                                                                w[ind * nsub + step-1]=std::pow(en,nu);
                                                                                
                                                                                //Include lower-point contributions
                                                                                for (int j=0; j<step; ++j)
                                                                                    for (int k=0; k<step-1; ++k)
                                                                                        w[ind * nsub + k] -= w[(ind - ppow2[i[j]]) * nsub + k]/(step-1-k);
                                                                                
                                                                                //Distance
                                                                                dist[ind] = dist[ind-ppow2[i[0]]];
                                                                                for (int j=1; j<step; ++j)
                                                                                    dist[ind] = std::max(dist[ind], partdist[i[0] * nsub + i[j]]);
                                                                            }
                }
                
                //16 particles
                if (b1.size()>=16){
                    int step = 16;
                    for (i[0]=0; i[0]<b1.size()-15; ++i[0])
                        for (i[1]=i[0]+1; i[1]<b1.size()-14; ++i[1])
                            for (i[2]=i[1]+1; i[2]<b1.size()-13; ++i[2])
                                for (i[3]=i[2]+1; i[3]<b1.size()-12; ++i[3])
                                    for (i[4]=i[3]+1; i[4]<b1.size()-11; ++i[4])
                                        for (i[5]=i[4]+1; i[5]<b1.size()-10; ++i[5])
                                            for (i[6]=i[5]+1; i[6]<b1.size()-9; ++i[6])
                                                for (i[7]=i[6]+1; i[7]<b1.size()-8; ++i[7])
                                                    for (i[8]=i[7]+1; i[8]<b1.size()-7; ++i[8])
                                                        for (i[9]=i[8]+1; i[9]<b1.size()-6; ++i[9])
                                                            for (i[10]=i[9]+1; i[10]<b1.size()-5; ++i[10])
                                                                for (i[11]=i[10]+1; i[11]<b1.size()-4; ++i[11])
                                                                    for (i[12]=i[11]+1; i[12]<b1.size()-3; ++i[12])
                                                                        for (i[13]=i[12]+1; i[13]<b1.size()-2; ++i[13])
                                                                            for (i[14]=i[13]+1; i[14]<b1.size()-1; ++i[14])
                                                                                for (i[15]=i[14]+1; i[15]<b1.size(); ++i[15]){
                                                                                    
                                                                                    //Calculate the index
                                                                                    int ind = 0;
                                                                                    double en = 0;
                                                                                    for (int j=0; j<step; ++j){
                                                                                        ind+=ppow2[i[j]];
                                                                                        en+=ptfrac[i[j]];
                                                                                    }
                                                                                    
                                                                                    //Calculate the direct contribution
                                                                                    w[ind * nsub + step-1]=std::pow(en,nu);
                                                                                    
                                                                                    //Include lower-point contributions
                                                                                    for (int j=0; j<step; ++j)
                                                                                        for (int k=0; k<step-1; ++k)
                                                                                            w[ind * nsub + k] -= w[(ind - ppow2[i[j]]) * nsub + k]/(step-1-k);
                                                                                    
                                                                                    //Distance
                                                                                    dist[ind] = dist[ind-ppow2[i[0]]];
                                                                                    for (int j=1; j<step; ++j)
                                                                                        dist[ind] = std::max(dist[ind], partdist[i[0] * nsub + i[j]]);
                                                                                }
                }
                
                //Add to histogram
                for (int j = 1; j<ppow2npt; ++j)
                    if (j>=ppow2[nb1] && (j % ppow2[nb1])>0)
                        for (int k=0; k<nsub; ++k)
                            res[dist[j]]+=w[j * nsub + k];
                
            } //Ensure that the memory is deallocated before every recursion
            
            //Recurse on the parents
            EEC(p1, jpt, minbin, nbins);
            EEC(p2, jpt, minbin, nbins);
        }
        else
        {   // Compute the correlator when all transverse momenta weights are kept on the same particle
            res[0] += std::pow(tree.pt()/jpt,nu);
        }
    }
};

#endif // EEC_H
