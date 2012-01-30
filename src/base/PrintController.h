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
#ifndef INCLUDED_PRINTCONTROLLER
#define INCLUDED_PRINTCONTROLLER
namespace MolDS_base{

class PrintController{
public:
   PrintController(){
      this->printsLogs = true;
      //cout << "printController is created.\n";
   }
   virtual ~PrintController(){
      //cout << "printController is destructed.\n";
   }
   bool PrintsLogs() const{
      return this->printsLogs;
   }
   void SetPrintsLogs(bool printsLogs){
      this->printsLogs = printsLogs;
   }
private:
   bool printsLogs;
};
}
#endif
