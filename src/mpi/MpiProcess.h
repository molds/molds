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
#ifndef INCLUDED_MPIPROCESS
#define INCLUDED_MPIPROCESS
namespace MolDS_mpi{

// MpiProcess is singleton
class MpiProcess: private MolDS_base::Uncopyable{
public:
   static void        CreateInstance(int argc, char *argv[]);
   static void        DeleteInstance();
   static MpiProcess* GetInstance();
   boost::mpi::communicator* GetCommunicator() const{ return this->communicator;}
private:
   static MpiProcess* mpiProcess;
   MpiProcess();
   MpiProcess(int argc, char *argv[]);
   ~MpiProcess();
   boost::mpi::environment*  environment;
   boost::mpi::communicator* communicator;
};

}
#endif

