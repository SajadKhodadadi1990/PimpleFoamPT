# pimpleFoamPT

![OpenFOAM](https://img.shields.io/badge/OpenFoam-darkblue.svg)
![OpenFOAM version](https://img.shields.io/badge/foam_extend-4.0-brightgreen.svg)

Coupled Eulerian-Lagrangian Immersed Boundary (IB) solver 

## Features
 * Coupling discrete IB method library to Eulerian-Lagrangian solver
 * Enhancing the accuracy of the default IB solver by introducing the corrected face


## Installation
### __1.  Copy files from `pimpleFoamPT` folder to  `src` and `applications/solvers/immersedBoundary/` folders of foam-extend-4.0__
   
```
git clone https://github.com/SajadKhodadadi1990/pimpleFoamPT.git
```

```
cd pimpleFoamPT/solver/
cp -r pimpleFoam  ~/foam/foam-extend-4.0/applications/solvers/immersedBoundary/
cd ../libraries/src/
cp -r immersedBoundar  ~/foam/foam-extend-4.0/src/
cp -r lagrangian  ~/foam/foam-extend-4.0/src/
```
### __2. Compile library sources__

* Immersed Boundary sources

    * immersed boundary
      ```
      cd ~/foam/foam-extend-4.0/src/immersedBoundar/immersedBoundary/
      fe40
      wclean
      wmake libso
      ```
    * Turbulence
      ```
      cd ../immersedBoundaryTurbulence/
      wclean
      wmake libso
      ```
    * Dynamic Mesh
       ```
      cd ../immersedBoundaryDynamicMesh/
      wclean
      wmake libso
       ```
    * Force
      ```
      cd ../immersedBoundaryForce/
      wclean
      wmake libso
      cd ../
      ```
* Lagrangian sources
  
  ```
  cd ../lagrangian
  ./Allwmake
  ```
### __3. Compile solver__
  ```
   cd  ~/foam/foam-extend-4.0/applications/solvers/immersedBoundary/pimpleFoam/
   wclean
   wmake
  ```

## Usage

  This solver can be utilized for various Lagrangian problems that involve complex geometries and are based on a Cartesian mesh.

## Citation
  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12520009.svg)](https://doi.org/10.5281/zenodo.12520009)

  [1]https://www.sciencedirect.com/science/article/abs/pii/S0921883124002401?via%3Dihub (https://doi.org/10.1016/j.apt.2024.104564)


  
  



       




