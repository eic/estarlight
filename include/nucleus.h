///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of starlight.
//
//    starlight is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    starlight is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with starlight. If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//
// File and Version Information:
// $Rev:: 262                         $: revision of last commit
// $Author:: jnystrand                $: author of last commit
// $Date:: 2016-06-01 14:14:20 +0100 #$: date of last commit
//
// Description:
//
//
//
///////////////////////////////////////////////////////////////////////////


#ifndef NUCLEUS_H
#define NUCLEUS_H


#include <cmath>


//This class holds the information for a target nucleus
class nucleus
{

public:
	nucleus();
	nucleus(const int    Z,
	        const int    A,
		const int     productionMode);
	~nucleus();
	
	void init();
 
	int    Z              () const { return _Z;                     }  ///< returns atomic number of nucleus
	int    A              () const { return _A;                     }  ///< returns nucleon number of nucleus
        int    productionMode () const { return _productionMode;        }

	double formFactor(const double t) const;
	// Calculates form factor for given squared 4-momentum transfer

	double dipoleFormFactor(const double t, const double t0) const;
	// Calculates dipole form factor with t0 as parameter 

	double thickness (const double b) const;
	// Calculates nuclear thickness function 

	double nuclearRadius() const { return _Radius; }
	double rho0() const { return _rho0; }
	double woodSaxonSkinDepth() const { return _woodSaxonSkinDepth;}
	/** The "wine_bottle" parameter
	 */
	double w() const {return _w;}
	
	
private:

	 
	
    double rws(const double r) const;
	int    _Z;                      ///< atomic number of nucleus
	int    _A;                      ///< nucleon number of nucleus
    int    _productionMode;

	double _r0;
	double _Radius;
	double _rho0;
	double _w;               ///< 3pF fermi model parameter
	double _woodSaxonSkinDepth;
	double _v; ///< Sum of Gaussian(SOG) model parameter 
	double *_a;  // SOG model parameter
        double *_b;  // SOG model parameter
        double *_c;  // SOG model parameter

//SOG model parameters https:doi.org/10.1016/0092-640X(87)90013-1

	double pkCarb[12] = {0.009560250143852663, 0.02161993792818375, 0.023894326875090147, 0.022871610300732596, 0.01787510497352304, 0.013266599310805718, 0.0020798983410400215, 0.0012487276970433168, 0.0001121651485588463, 0.00001856133380347859, 0.00000004241988368721027, 0.0};
        double okCarb[12] = {0.016690, 0.050325, 0.128621, 0.180515, 0.219097, 0.278416, 0.058779, 0.057817, 0.007739, 0.002001, 0.000007, 0.0};
        double qkCarb[12] = {0.0, 0.4, 1.0, 1.3, 1.7, 2.3, 2.7, 3.5, 4.3, 5.4, 6.7, 0.0};
    
	double rkHel[12] = {0.2, 0.6, 0.9, 1.4, 1.9, 2.3, 2.6, 3.1, 3.5, 4.2, 4.9, 5.2};
        double tkHel[12] = {0.034724, 0.430761, 0.203166, 0.192986, 0.083866, 0.033007, 0.014201, 0.0, 0.006860, 0.0, 0.000438, 0.0};
	double skHel[12] = {0.010229060308680777, 0.0683272860306785, 0.01954234673886537, 0.009254567928457054, 0.002338942473693466, 0.0006455191775197494, 0.00022017367746850975, 0.0, 0.000059954958779193, 0.0, 0.000001978748821697672, 0.0};
};


#endif  // NUCLEUS_H





