/*

    ________/\\\\\\\\\__/\\\\\\\\\\\\\____/\\\_____________        
     _____/\\\////////__\/\\\/////////\\\_\/\\\_____________       
      ___/\\\/___________\/\\\_______\/\\\_\/\\\_____________      
       __/\\\_____________\/\\\\\\\\\\\\\/__\/\\\_____________     
        _\/\\\_____________\/\\\/////////____\/\\\_____________    
         _\//\\\____________\/\\\_____________\/\\\_____________   
          __\///\\\__________\/\\\_____________\/\\\_____________  
           ____\////\\\\\\\\\_\/\\\_____________\/\\\\\\\\\\\\\\\_ 
            _______\/////////__\///______________\///////////////__


                         C P L  -  L I B R A R Y 

           Copyright (C) 2012-2015 Edward Smith & David Trevelyan

License

    This file is part of CPL-Library.

    CPL-Library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CPL-Library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with CPL-Library.  If not, see <http://www.gnu.org/licenses/>.

Description

    "Socket" class for interfacing with CPL-Library. 

*/
#ifndef CPLSocketFOAM_H
#define CPLSocketFOAM_H

#include <vector>
#include <memory>
#include "mpi.h"
#include "fvCFD.H"
#include "error.H"
#include "cpl/cpl.h"
#include "cpl/CPL_ndArray.h"
#include "cpl/CPL_cartCreate.h"
#include "meshSearch.H"
#include "cpl/CPLSocket.h"
// Nearest integer
int nint(double d) {
    return std::floor(d+0.5);
}

class CPLSocketFOAM: public CPLSocket {   

public:
    

    // Construct from no arguments 
    CPLSocketFOAM() : CPLSocket(CPL::cfd_realm) {}
    void setTimingInfo();
    void setCartCommInfo();
    void setRealmDomainInfo();
    void setOpenFOAM(const Foam::Time &runTime, const Foam::fvMesh &mesh);

    const Foam::Time* foamRunTime;
    const Foam::fvMesh* foamMesh; 
    
    // Invoke CPL_Library and store CFD realm communicator
    void initComms(int& argc, char**& argv);
    void init();
    double CPLDensity;

	virtual ~CPLSocketFOAM() { }

};


#endif
