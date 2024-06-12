// FastEEC 0.1 -- This header file is a part of FastEEC. This code implements
// the fast evaluation method for N point correlators clustered with C/A or kt
// clustering with a dynamical jet resolution f. The code contains N point
// correlators up to N = 8.
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
#include <fstream>
#include <sstream>
#include <cmath>
#include "fastjet/PseudoJet.hh"

const double PI = 3.14159265358979323846;
const double TWOPI = 2 * PI;

//Global variables that can be accessed directly in the main code eec_fast.cc/eec_fast_kt.cc
const double maxbin = 0;
std::vector<double> res(1000);
double fac;
//We have chosen a jet radius R that assumes that the provided input is clustered into one jet
const double R = 1.5;

//Function to compute distance in (eta,phi)
double ang(const std::vector<fastjet::PseudoJet>& v, int m, int n) {
    double dphiabs = std::fabs(v[m].phi() - v[n].phi());
    double dphi = dphiabs > PI ? TWOPI - dphiabs : dphiabs;
    double dy = v[m].rap() - v[n].rap();
    return std::sqrt(dy*dy + dphi*dphi);
}

double ang(const fastjet::PseudoJet& v, const fastjet::PseudoJet& w) {
    double dphiabs = std::fabs(v.phi() - w.phi());
    double dphi = dphiabs > PI ? TWOPI - dphiabs : dphiabs;
    double dy = v.rap() - w.rap();
    return std::sqrt(dy*dy + dphi*dphi);
}

class EEC {
public:
    EEC(fastjet::PseudoJet tree, double jpt, int n, const double minbin, const int nbins) {
        compute(tree, jpt, n, minbin, nbins);
    }
    
    void compute(const fastjet::PseudoJet& tree, double jpt, int n, const double minbin, const int nbins)
    { // Pass by reference
        fastjet::PseudoJet p1, p2;
        
        if (tree.has_parents(p1, p2))
        {
            //Calculate the distance in (eta,phi) between the parents
            const double th = ang(p1,p2);
            
            std::vector<fastjet::PseudoJet> b1, b2;
            
            const double dist = ((th * th)/(R*R));
            
            //Resolve as much of each of the branches as needed
            b1 = p1.exclusive_subjets(dist / fac);
            b2 = p2.exclusive_subjets(dist / fac);
            
            //Determine how many particles in each branch
            const int nb1 = b1.size(), nb2 = b2.size();
            const int factorial[10] = {1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800};
            
            //Append b2 to b1 for convenience
            b1.insert(b1.end(), b2.begin(), b2.end());
            //Precalculate all bin positions
            double bindist[200][200];
            
            for (int j1=0; j1<nb1+nb2; ++j1)
                for (int j2=j1+1; j2<nb1+nb2; ++j2)
                {
                    double chi = ang(b1[j1],b1[j2]);
                    
                    //Convert this to a bin position
                    int binpos = chi <= pow(10,minbin) ? 0 : std::trunc((log10(chi)-minbin)/(maxbin-minbin)*nbins);
                    
                    //Handle overflow and underflow
                    if (binpos < 0)
                        binpos = 0;
                    if (binpos >= nbins)
                        binpos = nbins-1;
                    
                    bindist[j1][j2]=binpos;
                }
            switch(n)
            {
                case 2:
                { 
                    const int npt = 2;
                    int i[npt];
                    //Loop over all pseudojets, with i in branch 1 and j in branch 2
                    for (i[0] = 0; i[0] < nb1; ++i[0])
                        for (i[1] = std::max(i[0], nb1); i[1] < nb1 + nb2; ++i[1])
                        {
                            int binpos = 0;
                            //Calculate the max of all pairs of distances
                            for (int j1 = 0; j1 < npt; ++j1)
                                for (int j2 = j1 + 1; j2 < npt; ++j2)
                                    binpos = (binpos > bindist[i[j1]][i[j2]] ? binpos : bindist[i[j1]][i[j2]]);
                            
                            //Number of permutations
                            int perms = 2;
                            
                            //Add to histogram, including permutation factor and weight of particles transverse momenta divided by the total jet pt
                            double efac = 1;
                            for (int j = 0; j < npt; ++j)
                                efac *= b1[i[j]].pt() / jpt;
                            
                            res[binpos] += perms * efac;
                        }
                }
                    break;
                case 3:
                {
                    const int npt = 3;
                    int i[npt];
                    //Loop over all pseudojets, with i in branch 1 and m in branch 2
                    for (i[0] = 0; i[0] < nb1; ++i[0])
                        for (i[1] = i[0]; i[1] < nb1 + nb2; ++i[1])
                            for (i[2] = std::max(i[1], nb1); i[2] < nb1 + nb2; ++i[2])
                            {
                                int binpos = 0;
                                //Calculate the max of all pairs of distances
                                for (int j1 = 0; j1 < npt; ++j1)
                                    for (int j2 = j1 + 1; j2 < npt; ++j2)
                                        binpos = (binpos > bindist[i[j1]][i[j2]] ? binpos : bindist[i[j1]][i[j2]]);
                                
                                //Calculate number of permutations
                                int perm = 6;
                                if (i[0] == i[1] && i[1] == i[2]) perm = 1;
                                if (i[0] == i[1] || i[1] == i[2] || i[0] == i[2]) perm = 3;
                                
                                //Add to histogram, including permutation factor and weight of particles transverse momenta divided by the total jet pt
                                double efac = 1;
                                for (int j = 0; j < npt; ++j)
                                    efac *= b1[i[j]].pt() / jpt;
                                
                                res[binpos] += perm * efac;
                            }
                }
                    break;
                case 4:
                {
                    const int npt = 4;
                    int i[npt];
                    //Loop over all pseudojets, with i in branch 1 and l in branch 2
                    for (i[0] = 0; i[0] < nb1; ++i[0])
                        for (i[1] = i[0]; i[1] < nb1 + nb2; ++i[1])
                            for (i[2] = i[1]; i[2] < nb1 + nb2; ++i[2])
                                for (i[3] = std::max(i[2], nb1); i[3] < nb1 + nb2; ++i[3])
                                {
                                    int binpos = 0;
                                    //Calculate the max of all pairs of distances
                                    for (int j1 = 0; j1 < npt; ++j1)
                                        for (int j2 = j1 + 1; j2 < npt; ++j2)
                                            binpos = (binpos > bindist[i[j1]][i[j2]] ? binpos : bindist[i[j1]][i[j2]]);
                                    
                                    //Calculate number of permutations
                                    int perm = factorial[npt-1];
                                    
                                    int oldval=i[0], cnt=1;
                                    for (int j=1; j<npt; ++j)
                                        if (i[j]!=oldval)
                                        {
                                            perm/=factorial[cnt-1];
                                            
                                            oldval=i[j];
                                            cnt=1;
                                        } else {
                                            cnt++;
                                        }
                                    perm/=factorial[cnt-1];
                                    
                                    //Add to histogram, including permutation factor and weight of particles transverse momenta divided by the total jet pt
                                    double efac = 1;
                                    for (int j = 0; j < npt; ++j)
                                        efac *= b1[i[j]].pt() / jpt;
                                    
                                    res[binpos] += perm * efac;
                                }
                }
                    break;
                case 5:
                {
                    const int npt = 5;
                    int i[npt];
                    //Loop over all pseudojets, with i in branch 1 and m in branch 2
                    for (i[0] = 0; i[0] < nb1; ++i[0])
                        for (i[1] = i[0]; i[1] < nb1+nb2; ++i[1])
                            for (i[2] = i[1]; i[2] < nb1+nb2; ++i[2])
                                for (i[3] = i[2]; i[3] < nb1+nb2; ++i[3])
                                    for (i[4] = std::max(i[3],nb1); i[4] < nb1+nb2; ++i[4])
                                    {
                                        int binpos = 0;
                                        //Calculate the max of all pairs of distances
                                        for (int j1 = 0; j1 < npt; ++j1)
                                            for (int j2 = j1 + 1; j2 < npt; ++j2)
                                                binpos = (binpos > bindist[i[j1]][i[j2]] ? binpos : bindist[i[j1]][i[j2]]);
                                        
                                        //Calculate number of permutations
                                        int perm = factorial[npt-1];
                                        
                                        int oldval=i[0], cnt=1;
                                        for (int j=1; j<npt; ++j)
                                            if (i[j]!=oldval)
                                            {
                                                perm/=factorial[cnt-1];
                                                
                                                oldval=i[j];
                                                cnt=1;
                                            } else {
                                                cnt++;
                                            }
                                        perm/=factorial[cnt-1];
                                        
                                        //Add to histogram, including permutation factor and weight of particles transverse momenta divided by the total jet pt
                                        double efac = 1;
                                        for (int j = 0; j < npt; ++j)
                                            efac *= b1[i[j]].pt() / jpt;
                                        
                                        res[binpos] += perm * efac;
                                    }
                }
                    break;
                case 6:
                {
                    const int npt = 6;
                    int i[npt];
                    //Loop over all pseudojets, with i in branch 1 and m in branch 2
                    for (i[0] = 0; i[0]<nb1; ++i[0])
                        for (i[1] = i[0]; i[1]<nb1+nb2; ++i[1])
                            for (i[2] = i[1]; i[2]<nb1+nb2; ++i[2])
                                for (i[3] = i[2]; i[3]<nb1+nb2; ++i[3])
                                    for (i[4] = i[3]; i[4]<nb1+nb2; ++i[4])
                                        for (i[5] = std::max(i[4],nb1); i[5]<nb1+nb2; ++i[5])
                                        {
                                            int binpos = 0;
                                            //Calculate the max of all pairs of distances
                                            for (int j1=0; j1<npt; ++j1)
                                                for (int j2=j1+1; j2<npt; ++j2)
                                                    binpos = (binpos > bindist[i[j1]][i[j2]] ? binpos : bindist[i[j1]][i[j2]]);
                                            
                                            //Calculate number of permutations
                                            int perm = factorial[npt-1];
                                            
                                            int oldval=i[0], cnt=1;
                                            for (int j=1; j<npt; ++j)
                                                if (i[j]!=oldval)
                                                {
                                                    perm/=factorial[cnt-1];
                                                    
                                                    oldval=i[j];
                                                    cnt=1;
                                                } else {
                                                    cnt++;
                                                }
                                            perm/=factorial[cnt-1];
                                            
                                            //Add to histogram, including permutation factor and weight of particles transverse momenta divided by the total jet pt
                                            double efac = 1;
                                            for (int j = 0; j<npt; ++j)
                                                efac *= (b1[i[j]].pt()/jpt);
                                            res[binpos] += perm * efac;
                                        }
                }
                    break;
                case 7:
                {
                    const int npt = 7;
                    int i[npt];
                    //Loop over all pseudojets, with i in branch 1 and m in branch 2
                    for (i[0] = 0; i[0]<nb1; ++i[0])
                        for (i[1] = i[0]; i[1]<nb1+nb2; ++i[1])
                            for (i[2] = i[1]; i[2]<nb1+nb2; ++i[2])
                                for (i[3] = i[2]; i[3]<nb1+nb2; ++i[3])
                                    for (i[4] = i[3]; i[4]<nb1+nb2; ++i[4])
                                        for (i[5] = i[4]; i[5]<nb1+nb2; ++i[5])
                                            for (i[6] =  std::max(i[5],nb1); i[6]<nb1+nb2; ++i[6])
                                            {
                                                int binpos = 0;
                                                //Calculate the max of all pairs of distances
                                                for (int j1=0; j1<npt; ++j1)
                                                    for (int j2=j1+1; j2<npt; ++j2)
                                                        binpos = (binpos > bindist[i[j1]][i[j2]] ? binpos : bindist[i[j1]][i[j2]]);
                                                
                                                //Calculate number of permutations
                                                int perm = factorial[npt-1];
                                                
                                                int oldval=i[0], cnt=1;
                                                for (int j=1; j<npt; ++j)
                                                    if (i[j]!=oldval)
                                                    {
                                                        perm/=factorial[cnt-1];
                                                        
                                                        oldval=i[j];
                                                        cnt=1;
                                                    } else {
                                                        cnt++;
                                                    }
                                                perm/=factorial[cnt-1];
                                                
                                                //Add to histogram, including permutation factor and weight of particles transverse momenta divided by the total jet pt
                                                double efac = 1;
                                                for (int j = 0; j<npt; ++j)
                                                    efac *= (b1[i[j]].pt()/jpt);
                                                res[binpos] += perm * efac;
                                            }
                }
                    break;
                case 8:
                {
                    const int npt = 8;
                    int i[npt];
                    //Loop over all pseudojets, with i in branch 1 and m in branch 2
                    for (i[0] = 0; i[0]<nb1; ++i[0])
                        for (i[1] = i[0]; i[1]<nb1+nb2; ++i[1])
                            for (i[2] = i[1]; i[2]<nb1+nb2; ++i[2])
                                for (i[3] = i[2]; i[3]<nb1+nb2; ++i[3])
                                    for (i[4] = i[3]; i[4]<nb1+nb2; ++i[4])
                                        for (i[5] = i[4]; i[5]<nb1+nb2; ++i[5])
                                            for (i[6] = i[5]; i[6]<nb1+nb2; ++i[6])
                                                for (i[7] =  std::max(i[6],nb1); i[7]<nb1+nb2; ++i[7])
                                                {
                                                    int binpos = 0;
                                                    //Calculate the max of all pairs of distances
                                                    for (int j1=0; j1<npt; ++j1)
                                                        for (int j2=j1+1; j2<npt; ++j2)
                                                            binpos = (binpos > bindist[i[j1]][i[j2]] ? binpos : bindist[i[j1]][i[j2]]);
                                                    
                                                    //Calculate number of permutations
                                                    int perm = factorial[npt-1];
                                                    
                                                    int oldval=i[0], cnt=1;
                                                    for (int j=1; j<npt; ++j)
                                                        if (i[j]!=oldval)
                                                        {
                                                            perm/=factorial[cnt-1];
                                                            
                                                            oldval=i[j];
                                                            cnt=1;
                                                        } else {
                                                            cnt++;
                                                        }
                                                    perm/=factorial[cnt-1];
                                                    
                                                    //Add to histogram, including permutation factor and weight of particles transverse momenta divided by the total jet pt
                                                    double efac = 1;
                                                    for (int j = 0; j<npt; ++j)
                                                        efac *= (b1[i[j]].pt()/jpt);
                                                    res[binpos] += perm * efac;
                                                }
                }
                    break;
            }
            //Recurse on the parents
            compute(p1, jpt, n, minbin, nbins);
            compute(p2, jpt, n, minbin, nbins);
        } else {
            // Compute the correlator when all transverse momenta weights are kept on the same particle
            res[0] += std::pow((tree.pt() / jpt), n);
        }
    }
};

#endif // EEC_H
