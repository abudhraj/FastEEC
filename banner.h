// FastEEC 0.1. This file is a part of FastEEC. The 
// file generates the banner for the project developed 
// by A. Budhraja and W. Waalewijn. 
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, version 2 of the License.

#ifndef BANNER_H
#define BANNER_H

#include <iostream>

// Define the printHeader function
void printHeader() {
    std::cout << "\n";
    std::cout << "             FastEEC 0.1, developed by A. Budhraja and W. Waalewijn\n";
    std::cout << "                     https://github.com/abudhraj/FastEEC\n";
    std::cout << "\n";
    std::cout << "   Please cite arXiv:2406.XXXX if you use this package for scientific work.\n";
    std::cout << "       It uses FastJet package by M. Cacciari, G. Salam and G. Soyez.\n";
    std::cout << "                         See license for details.\n";
    std::cout << "\n";
}

#endif // BANNER_H
