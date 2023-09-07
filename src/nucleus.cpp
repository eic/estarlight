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


#include <iostream>
#include <fstream>

#include "starlightconstants.h"
#include "reportingUtils.h"
#include "nucleus.h"
#include <inputParameters.h>


using namespace std;
using namespace starlightConstants;

//______________________________________________________________________________
nucleus::nucleus(const int    Z,
                 const int    A,
		 const int    productionMode)
	: _Z(Z),
	  _A(A),
	  _productionMode(productionMode)
{
  init();	
}

void nucleus::init()
{
  switch (_Z) {
	case 82:
		{
		  _Radius = 6.624;
		  _rho0 = 0.160696;
		}
		break;
	case 79:
		{
		  _Radius = 6.38;
		  _rho0 = 0.169551;
		}
		break;
	case 29:
		{
                  _Radius = 4.214;
		  _rho0 = 0.173845;
		}
		break;
	case 1: 
		{
		  //is this a proton or deuteron
		  if(_A==1){
		    _Radius = 0.7;
		    _rho0 = -1.0; //Not relevant for protons
		  }
		  // this is electron
		  else if(_A==0){
		    _Radius = 0.;
		    _rho0 = -1.0;
		  }
		else
		{
			_Radius = 2.127; // Measured rms radius for deuteron
			_rho0 = _A;
		}
	}
	break;
	case -1:
	{
		  //is this a anti-proton
		  if(_A==1){
		    _Radius = 0.7;
		    _rho0 = -1.0; //Not relevant for protons
		  }
		  // this is positron
		  else if(_A==0){
		    _Radius = 0.;
		    _rho0 = -1.0;
		  }
		  else {
		    _Radius = 2.1;
		}
	}
	break;
	case 2: 
	{
		// Helium-4 (sum of gaussian model)
		_a = rkHel;
		_c = skHel;
		_v = 0.81649;
		_rho0 = 2.0;
		_Radius = 1.676;
	}
	break;
	case 3: 
	{
		// Lithium-7 (Harmonic Oscillator Model)
		_Radius = 1.77;
		_g = 0.327;
		_rho0 = 0.1520360869; // determined by normalization
	}
	break;
	case 4: 
	{
		// Beryllium(9/10)
		if (_A == 9)
		{
			//Harmonic Oscillator Model
			_Radius = 1.791;
			_g = 0.611;
			_rho0 = 0.1467396; // determined by normalization
		}
		else
		{
			_Radius = 2.357; //measured rms radius
			_rho0 = _A;
		}
	}
	break;
	case 5: 
	{
		// Boron(10/11)
		if (_A == 10)
		{ 
			//Harmonic Oscillator Model
			_Radius = 1.71;
			_g = 0.837;
			_rho0 = 0.159172948; // determined by normalization
		}
		else
		{
			_Radius = 2.42;
			_rho0 = _A;
		}
	}
	break;
	case 6: 
	{
		// Carbon(12/13)
		if (_A == 13)
		{
			//Harmonic Oscillator Model
			_Radius = 1.635;
			_g = 1.403;
			_rho0 = 0.1719882747; // determined by normalization
		}
		else
		{
			//sum of gaussian model
			_a = qkCarb;
			_c = pkCarb;
			_v = 0.97978;
			_rho0 = 2.0;
			_Radius = 2.469;
		}
	}
	break;
	case 7: 
	{
		// Nitrogen(14)
		_Radius = 2.570;
		_woodSaxonSkinDepth = 0.5052;
		_w = -0.180;	   // 3pF fermi model wine parameter
		_rho0 = 0.1790484; // determined by normalization
	}
	break;
	case 8:
	{
		// This is for oxygen(sum of gaussian model)
		_a = xkOxyg;		
		_c = ykOxyg;
		_v = 1.061437;
		_rho0 = 2.0;
		_Radius = 2.711;

		//_Radius = 2.608;
		//_woodSaxonSkinDepth = 0.513;
		//_w= -0.051; // 3pF fermi model wine parameter
		//_rho0 = 0.16536; // determined by normalization
	}
	break;
	case 10: 
	{
		// This is for neon-20
		_Radius = 2.791;
		_woodSaxonSkinDepth = 0.698;
		_w = -0.168;	   // 3pF fermi model wime parameter
		_rho0 = 0.1798707; // determined by normalization
	}
	break;
	case 12: 
	{
		// This is for magnesium-24
		_Radius = 3.108;
		_woodSaxonSkinDepth = 0.607;
		_w = -0.163;		// 3pF fermi model wime parameter
		_rho0 = 0.16981488; // determined by normalization
	}
	break;
		// default:
		// printWarn << "density not defined for projectile with Z = " << _Z << ". using defaults." << endl;
		//       _Radius = 1.2*pow(_A, 1. / 3.);
		//_rho0 = 0.138/(1.13505-0.0004283*_A);  //This matches the radius above
		// if( _Z < 7 ){
		// This is for Gaussian form factors/densities
		//_rho0 = _A;
		//}
	}
}

//______________________________________________________________________________
nucleus::~nucleus()
{
}

//______________________________________________________________________________
double
nucleus::rws(const double r) const
{
	if (((_Z == 2) && (_A == 4)) || ((_Z == 6) && (_A == 12)) || (_Z == 8))
	{
		// Sum of Gaussian Model
		double _P = 0.0;
		for (int i = 0; i < 12; i++) {
		    _P += _c[i] * (exp(-1.0 * pow((r -_a[i]) /_v, 2)) + exp(-1.0 * pow((r +_a[i]) /_v, 2)));    
		}
		return _P;
	}
	else if (_Z < 7)
	{
		// Gaussian density distribution for light nuclei
		double norm = (3. / (2. * starlightConstants::pi)) * sqrt((3. / (2. * starlightConstants::pi)));
		norm = norm / (nuclearRadius() * nuclearRadius() * nuclearRadius());
		return norm * exp(-((3. / 2.) * r * r) / (nuclearRadius() * nuclearRadius()));
		// if( _Z <7 ){
		//  Harmonic Oscillator Model for light nuclei
		// double y=exp(-(r*r/nuclearRadius()/nuclearRadius()));
		// return y*(1.0 +g()*r*r/nuclearRadius()/nuclearRadius());
	}
	else if (_Z == 7)
	{
		// 3-parameter-Fermi Model (3pF)
		double x = exp(-(r - nuclearRadius()) / woodSaxonSkinDepth());
		return x * (1.0 + w() * r * r / nuclearRadius() / nuclearRadius()) / (1.0 + x);
		// expression adjusted to avoid problems on some machines if r is too large
	}
	else
	{

		// Fermi density distribution for heavy nuclei
		return 1.0 / (1. + exp((r - nuclearRadius()) / woodSaxonSkinDepth()));
	}
}

//______________________________________________________________________________
double
nucleus::formFactor(const double t) const
{
	// electromagnetic form factor of proton
	if ((_Z == 1) && (_A == 1))
	{
		const double rec = 1. / (1. + t / 0.71);
		return rec * rec;
	}

	if (_Z < 7)
	{
		// Gaussian form factor for light nuclei
		const double R_G = nuclearRadius();
		return exp(-R_G * R_G * t / (6. * starlightConstants::hbarc * starlightConstants::hbarc));
	}
	else
	{
		// nuclear form factor, from Klein Nystrand PRC 60 (1999) 014903, Eq. 14
		const double R = nuclearRadius();
		const double q = sqrt(t);
		const double arg1 = q * R / hbarc;
		const double arg2 = hbarc / (q * R);
		const double sph = (sin(arg1) - arg1 * cos(arg1)) * 3. * arg2 * arg2 * arg2;
		const double a0 = 0.70; // [fm]
		return sph / (1. + (a0 * a0 * t) / (hbarc * hbarc));
	}
}

//______________________________________________________________________________
double
nucleus::dipoleFormFactor(const double t, const double t0) const
{
	const double rec = 1. / (1. + t / t0);
	return rec * rec;
}

//______________________________________________________________________________
double
nucleus::thickness(const double b) const
{
	//    JS      This code calculates the nuclear thickness function as per Eq. 4 in
	//    Klein and Nystrand, PRC 60.
	//    former DOUBLE PRECISION FUNCTION T(b)

	// data for Gauss integration
	const unsigned int nmbPoints = 5;
	const double xg[nmbPoints + 1] = {0., 0.1488743390, 0.4333953941, 0.6794095683,
									  0.8650633667, 0.9739065285};
	const double ag[nmbPoints + 1] = {0., 0.2955242247, 0.2692667193, 0.2190863625,
									  0.1494513492, 0.0666713443};

	const double zMin = 0;
	const double zMax = 15;
	const double zRange = 0.5 * (zMax - zMin);
	const double zMean = 0.5 * (zMax + zMin);
	double sum = 0;
	for (unsigned int i = 1; i <= nmbPoints; ++i)
	{
		double zsp = zRange * xg[i] + zMean;
		double radius = sqrt(b * b + zsp * zsp);
		sum += ag[i] * rws(radius);
		zsp = zRange * (-xg[i]) + zMean;
		radius = sqrt(b * b + zsp * zsp);
		sum += ag[i] * rws(radius);
	}

	return 2. * zRange * sum;
}
