//************************************************************************//
// Copyright (C) 2011-2012 Mikiya Fujii                                   // 
//                                                                        // 
// This file is part of MolDS.                                            // 
//                                                                        // 
// MolDS is free software: you can redistribute it and/or modify          // 
// it under the terms of the GNU General Public License as published by   // 
// the Free Software Foundation, either version 3 of the License, or      // 
// (at your option) any later version.                                    // 
//                                                                        // 
// MolDS is distributed in the hope that it will be useful,               // 
// but WITHOUT ANY WARRANTY; without even the implied warranty of         // 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          // 
// GNU General Public License for more details.                           // 
//                                                                        // 
// You should have received a copy of the GNU General Public License      // 
// along with MolDS.  If not, see <http://www.gnu.org/licenses/>.         // 
//************************************************************************//
#ifndef INCLUDED_UTILITIES
#define INCLUDED_UTILITIES
namespace MolDS_base{
// output welcome message
void OutputWelcomeMessage();
// output farewell message
void OutputFarewellMessage(time_t startTime, clock_t startTick, double ompStartTime, bool runingNormally);
// string for today.
std::string GetDateString();
// trim the string
std::string TrimString(const std::string str);
// number to string
// ex. Num2String(23,5) = "00023";
std::string Num2String(int number, int digit);

}
#endif

