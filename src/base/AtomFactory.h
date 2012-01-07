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
#ifndef INCLUDED_ATOMFACTORY
#define INCLUDED_ATOMFACTORY
namespace MolDS_base{

// AtomFactory is singleton
class AtomFactory: private Uncopyable{
public:
   static AtomFactory* GetInstance();
   static void DeleteInstance();
   static MolDS_base_atoms::Atom* CreateAtom(MolDS_base::AtomType atomType,
                                             double x,
                                             double y,
                                             double z,
                                             double px,
                                             double py,
                                             double pz);
   static MolDS_base_atoms::Atom* CreateAtom(MolDS_base::AtomType atomType,
                                             double x,
                                             double y,
                                             double z);
private:
   static AtomFactory* atomFactory;
   AtomFactory();
   ~AtomFactory();
};

}
#endif





