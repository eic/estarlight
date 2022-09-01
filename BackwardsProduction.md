# Backwards Production

## $\omega \rightarrow \pi^0+\gamma$ Decay
Default decay for $\omega$ is $\omega \rightarrow 2\pi$ with branch ratio 0.0153
$\omega \rightarrow \pi^0+\gamma$ decay (branching ratio 0.0828) can be enabled with `PROD_PID = 223022`
This neutral decay of the $\omega$ currently assumes an unpolarized $\omega$, as it is not clear whether $s$-channel helicity conservation holds in backwards production. This assumption can be changed in `src/gammaavm.cpp`


## Added Code 
Added check in `src/e_narrowResonanceCrossSection.cpp` for backwards production.

Switched to two vertices, gamma and electron, in `src/hepmc3writer.cpp` 

Added targetBeam, electronBeam and backwards_production tag for slight.in `src/inputParameters.cpp` and `include/inputParameters.h`  
  (Removed TAU+TAU decay)

Changed `src/starlightpythia.cpp` to produce eXEvent instead of UPCEvent

Backwards Production Code in `src/gammaavm.cpp`
```
if (_backwardsProduction) {
  //Fill the intial state in the lab frame
	double pz_tar_i_lab = -1.0*target_pz;  //Flip the pz of the target
	double E_tar_i_lab = sqrt(pow(target_pz,2) + pow(starlightConstants::protonMass,2.)); //Calculate proton energy
	double pz_gam_i_lab = Egam; //Photon energy
	double E_gam_i_lab = Egam;  //Photon energy

	//Calculate the final state target energy and pz
	double pz_tar_f_lab = (pz_tar_i_lab + pz_gam_i_lab) - pz;
	double E_tar_f_lab = sqrt( pow(pz_tar_f_lab,2) + pow(starlightConstants::protonMass,2.) );

	//Calculate the boost (beta,gamma) to the center of mass frame
	double beta = (pz_tar_i_lab + pz_gam_i_lab)/(E_gam_i_lab + E_tar_i_lab);
	double lorentzGamma = 1.0/sqrt( 1 - pow(beta,2.));
	  
	//Boost to the final state target to the CM frame
	double cm_frame_target_pz = lorentzGamma*( pz_tar_f_lab + beta*E_tar_f_lab );
	double cm_frame_gamma_pz  = lorentzGamma*( pz + beta*pz );
	 
	//Recalculate target energy in the CM frame (assign to VM) 
	double vm_mass = sqrt( E*E - ( px*px + py*py + pz*pz ) );
	double vm_E_cm = sqrt( vm_mass*vm_mass + px*px + py*py + cm_frame_target_pz*cm_frame_target_pz);

  //Calcuate photon energy in CM frame (assign to target)
  double cm_E = sqrt( pow(starlightConstants::protonMass,2.)  + pow(t_px,2.) + pow(t_py,2.) + pow(cm_frame_gamma_pz,2.));	  
  //Boost photon and target back to lab frame
  //target:
  t_pz = lorentzGamma*(cm_frame_gamma_pz - beta*cm_E);
  t_E = lorentzGamma*(cm_E - beta*cm_fame_gamma_pz);
  //Photon:
  pz = lorentzGamma*(cm_frame_target_pz - beta*vm_E_cm);
  E = lorentzGamma*(vm_E_cm - beta*cm_frame_target_pz);

  //Calculate the rapidity
  Y = 0.5*std::log( (E+fabs(pz))/(E-fabs(pz)) );
}  
```

## Removed STARlight Code (Conflicted with backwards production)
### DPMJET
Removed the option to include DPMJET in `CMakeLists.txt`
### BeamSystem
Changed the beam1/beam2 -> electronBeam/targetBeam. Defined in `include/beambeamsystem.h`
### UPCEvent
Removed all references to UPCEvents in:
```
include/e_starlightStandalone.h 
include/e_starlight.h 
include/eventchannel.h
include/eventfilewriter.h
include/filewriter.h
include/gammaavm.h
include/starlightpythia.h
include/pythiadecayer.h
```
### Updated kinematics:
```
src/gammaavm.cpp
```
