# PENSimulation
GEANT4 Simulation for viability of PEN in LEGEND200

Simulates 128 nm photons in liquid argon in a simplified GERDA like geometry. Records the number of events that reach the region where the fiber shroud would be, and the wavelength that was detected.

## How to install

To use: Navigate to build folder and enter:

cmake -DGeant4_DIR=/opt/geant4/lib64/Geant4-10.5.0 ..

make -j8

## How to use

To run the program in interactive mode:
./PEN

If you wish to use the program with a macro file:
./PEN -m /path/to/macro

## Macro commands

### Detector geometry and materials
Additional macro commands and cases:

/PEN/det/setTargetMat [material]

Set the support plates to the chosen material. The following materials have been tested:

G4_Si
PEN


/PEN/det/setDetectorType [int 0-1]

Controls the wls properites of PEN:

  0 - WLS off
  1 - WLS on

  In testing, probably won't work

### Run commands

/run/beamOn n

Run simulation for n events
