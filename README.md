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
~/the_path/eestarlight
```
* Change to the installation directory of your choice
```
mkdir ~/my/installation/dir
cd ~/my/installation/dir
```
* Set up the installation using cmake:
```
cmake ~/the_path/eestarlight
```
* Compile the code using (g)make:
```
(g)make
```
  * The compilation will produce two executables to run either STARlight or eSTARlight

* Set up the desired running conditions in the input file:
```
cp ~/the_path/eSTARlight/slight.in .
vim slight.in
```
  * Note: As of yet, the positron/electron beam must be set to beam1, it is selected by setting:
```
BEAM_1_Z = +/- 1 (as of yet both positrons and electrons are modelled identically)
BEAM_1_A = 0
```
* Run the simulation:
```
./e_starlight > output.txt
```
  * output.txt will contain the program log and calculated cross-section for the simulation sample.
  * The event catalogue will be emptied into the file slight.out
* Interpret the result. We have provided a macro to convert the output into a ROOT TTree:
  ```
  root -b -q -l ~/the_path/estarlight/utils/ConvertStarlightAsciiToTree.C
  ```

   * TTree is output to starlight.root
* Analyze your data:
   * We have included template [analysis](../blob/master/analysis)) (~/the_path/estarlight/analysis/) code that reads   slight.root and fill histograms. For more details please look at the [README](../blob/master/README.pdf).
```
cd ~/the_path/eSTARlight/analysis
sh e_run.sh ~/the_path/starlight.root
```

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
