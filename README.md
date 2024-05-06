# eSTARlight

eSTARlight is a Monte Carlo that simulates coherent vector meson photo- and electro-production in electron-ion collisions. It can produce a variety of final states at different center of mass energies for different collision systems at arbitrary values for the photon virtuality.


## Papers and presentations
* [Exclusive vector meson production at an electron-ion collider](https://arxiv.org/abs/1803.06420): M. Lomnitz & S. Klein, Physical Review C99, 015203 (2019); also available as arXiv:1803.06420.
* [Coherent vector meson production at an EIC](https://arxiv.org/abs/1805.08586): M. Lomnitz and S. Klein, Presented at the Workshop for Deep Inelastic Scattering, Kobe Japan, 2018.

## Authors
* Michael Lomnitz<sup>1</sup>
* Spencer Klein<sup>1</sup>

<sup>1</sup>: Lawrence Berkeley National Laboratory, Relativistic Nuclear Collisions Program, Nuclear Science Division.

All rights reserved.

## Declaration
Portions of this package were originally inherited/based on the [STARlight Monte Carlo](https://starlight.hepforge.org/) generator. We would like to acknowledge the authors J. Nystrand, J. Seger and Y. Gorbunov for their contributions to STARlight. This work was funded by the U.S. Department of Energy under contract number DE-AC-76SF00098.

## Instructions for use
The following instructions illustrate the procedure to install and run eSTARlight in a *nix based environment:

* Download the code package from Github and move to the desired location, i.e.
```
~/the_path/estarlight
```
* Change to the installation directory of your choice
```
mkdir ~/my/installation/dir
cd ~/my/installation/dir
```
* Set up the installation using cmake:
```
cmake ~/the_path/estarlight
```
* Compile the code using (g)make:
```
(g)make
```
* The compilation will produce an executable in your build directory.
If you want the  code to be globally accessible on the computer, run
`CMAKE_INSTALL_PREFIX` (defaults to `/usr/local`) and then make install.  This will probably
not be possible on large installations like Perlmutter. Alternately,
you can run the executable in your build directory.
```
(g)make install
```

* Set up the desired running conditions in the input file:
```
cp ~/the_path/estarlight/slight.in .
vim slight.in
```
  * Note: The electron beam energy is set by:
```
ELECTRON_BEAM_GAMMA = (Electron Energy)/(0.000511 GeV)
```
  * For example:
```
ELECTRON_BEAM_GAMMA = 9785 # 5GeV electrons from eRHIC
ELECTRON_BEAM_GAMMA = 19569 #10GeV electrons from eRHIC
ELECTRON_BEAM_GAMMA = 35225 #18GeV electrons from eRHIC
```
  * Note: The constraints on the center-of-mass energy of the vitrual photon and ion is set by:
```
W_GP_MAX = Max value of W_gp that user can set and wants to use (GeV)
W_GP_MIN = Min value of W_gp that user can set and wants to use (GeV)
```
 * For example:
 ```
 W_GP_MAX = 30 # W_gp will have a maximum value of 30 GeV.
 W_GP_MIN = 2 # W_gp will start to sample from the minimum value of 2 GeV. W_GP_MIN should not be greater than the center-of-mass energy of the electron and ion beams.
```
  * For exclusive backward (u-channel) production, use:
```
BACKWARDS_PRODUCTION = 1
```
* Run the simulation:
```
./e_starlight > output.txt
```
  * output.txt will contain the program log and calculated cross-section for the simulation sample.
  * The event catalogue will be emptied into the file slight.out
* Interpret the result. We have provided a macro to convert the output into a ROOT TTree:
  ```
  root -b -q -l ~/the_path/estarlight/utils/eTTreeMaker.C
  ```

   * TTree is output to ntuple_slight.root

## Documentation
A more detailed version of the [README](../blob/master/README.pdf) is included as part of the software package, located in:
```
~/the_path/estarlight/Readme.pdf
~/the_path/estarlight/Readme.docx
```

Finally, the full documentation with class description and dependencies is also available
as part of the EIC-related software documentation, e.g. [here](https://eic.github.io/doxygen/d9/da4/classe__starlight.html). We have included a config file to generate a fresh version of the doxygen documentation. The following steps can be used to generate the documentation:

* Download and install [Doxygen](https://www.stack.nl/~dimitri/doxygen/manual/install.html)
* Move into the source directory:
```
cd ~/the_path/estarlight/
```
If necessary, delete previous documentation to avoid any conflicts:
```
rm ~/the_path/estarlight/doxygen/html/*
```
Generate the documentation:
```
doxygen estarlightDoxyfile.conf
```
The documentation should be generated and available in the previous location.


To view the documentation:
- Open the file ~/the_path/estarlight/doxygen/html/index.html

## Using with HepMC3 output
HepMC3 can be found at https://gitlab.cern.ch/hepmc/HepMC3. The README should be referred to for installation. For quick installation, the following can be used.
```
mkdir hepmc3
cd hepmc3
git clone https://gitlab.cern.ch/hepmc/HepMC3.git
mkdir hepmc3-build
cd hepmc3-build
cmake -DHEPMC3_ENABLE_ROOTIO=OFF -DCMAKE_INSTALL_PREFIX=../hepmc3-install -DHEPMC3_ENABLE_PYTHON=OFF -DHEPMC3_BUILD_EXAMPLES=ON -DHEPMC3_ENABLE_TEST=ON ../HepMC3
```
Running tests at this point will result in errors so first compile the test directory:
```
cd test
make
cd ..
make test
```
You should get the output: `100% tests passed, 0 tests failed out of 23`.
Now install:
```
make install
```
To compile eSTARlight with HepMC3 output enabled, use:
```
cmake /pathto/estarlight -DENABLE_HEPMC3=ON -DHepMC3_DIR=/pathto/hepmc3/hepmc3-install
```
NOTE: running `make install` on your eSTARlight build when using HepMC3 may result in errors on some systems. If this is the case, it is recommended to not install eSTARlight when using HepMC3 and simply run the code from the build directory.
