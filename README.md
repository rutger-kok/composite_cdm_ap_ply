# 3D Continuum Damage Mechanics VUMAT for AP-PLY Composite Materials
Explicit material subroutine (VUMAT) implementing a continuum damage mechanics (CDM) model for Adavanced Placed Ply (AP-PLY) composite materials in Abaqus (in fixed format Fortran 77).  

## Summary
This CDM model has been developed specifically to model AP-PLY laminates containing through thickness fiber connectivity. The complete numerical methodology will be documented in an as yet to be published journal paper. The CDM model implemented in this subroutine is an amalgamation of the work of several authors. Failure criteria are adopted from the work of Tan et al. [1] and Shah et al. [2]. Damage evolution is defined according to the methodlogies proposed in the work of Maimi et al. [3][4] and Shah et al. [2].  

## Usage
To run a simulation using subroutines your Abaqus installation must be linked with a Fortran compiler and compatible Visual Studio installation.  

The model requires the following material properties to be defined in the simulation input (.inp) file:  
**E<sub>11</sub>** = elastic modulus in the fiber direction  
**E<sub>22</sub>**= elastic modulus transverse direction (in-plane)  
**E<sub>33</sub>** = elastic modulus transverse direction (out-of-plane)  
**ν<sub>12</sub>** = Poisson's ratio 12 direction  
**ν<sub>13</sub>** = Poisson's ratio 13 direction  
**ν<sub>23</sub>** = Poisson's ratio 23 direction  
**G<sub>12</sub>** = shear modulus 12 direction  
**G<sub>13</sub>** = shear modulus 13 direction  
**G<sub>23</sub>** = shear modulus 23 direction  
**X<sub>T</sub>** = tensile strength fiber direction  
**X<sub>C</sub>** = compressive strength fiber direction   
**Y<sub>T</sub><sup>is</sup>** =  in-situ tensile strength transverse direction  
**Y<sub>C</sub><sup>is</sup>** = compressive strength transverse direction  
**Z<sub>T</sub><sup>is</sup>** =  in-situ tensile strength through-thickness direction  
**Z<sub>C</sub><sup>is</sup>** = compressive strength through-thickness direction  
**S<sub>L</sub><sup>is</sup>** = in-situ 12 shear strength  
**S<sub>R</sub><sup>is</sup>** = in-situ 13 shear strength  
**S<sub>T</sub><sup>is</sup>** = 23 shear strength  
**G<sub>1+</sub>** =  tensile fracture toughness fiber direction  
**G<sub>1-</sub>** =  compressive fracture toughness fiber direction  
**G<sub>2+</sub>** = tensile fracture toughness transverse direction  
**und_angle<sub>1</sub>** = in-plane angle undulation region constituent 1  
**und_angle<sub>2</sub>** = in-plane angle undulation region constituent 2  
**cpt** = cured ply thickness  
**und_length** = undulation length 

## List of Fortran source code
- **composite_cdm.for** : Main CDM model script. Imports the resin_damage and rotation_matrix files.
- **resin_damage.for** : Implements an isotropic CDM model for pure resin regions
- **rotation_matrix.for** : Populates the rotation matrix for stress/strain transformation about 2 axes.

***
Rutger Kok  
PhD Candidate  
email : rutger.kok@ed.ac.uk  

Institute for Infrastructure and Environment  
University of Edinburgh    
Thomas Bayes Road, King's Buildings, Edinburgh EH9 3FG   
United Kingdom

***
>[1] W. Tan, B. G. Falzon, L. N. S. Chiu, and M. Price  
>Predicting low velocity impact damage and Compression-After-Impact (CAI) behaviour of composite laminates  
>Composites Part A 71 (2015) 212–226.  
>http://doi.org/10.1016/j.compositesa.2015.01.025  

>[2] S. Z. H. Shah, P. S. M. Megat-Yusoff, S. Karuppanan, R.S. Choudhry, and Z. Sajid  
>Multiscale damage modelling of 3D woven composites under static and impact loads  
>Composites Part A: Applied Science and Manufacturing 151 (2021).  
>https://doi.org/10.1016/j.compositesa.2021.106659  

>[3] P. Maimi, P.P. Camanho, J.A. Mayugo, C.G. Davila  
>A continuum damage model for composite laminates: Part I – Constitutive model  
>Mechanics of Materials 39 (2007) 897–908  
>http://dx.doi.org/10.1016/j.mechmat.2007.03.005  

>[4] P. Maimi, P.P. Camanho, J.A. Mayugo, C.G. Davila  
>A continuum damage model for composite laminates: Part II – Computational implementation and validation  
>Mechanics of Materials 39 (2007) 909–919  
>http://dx.doi.org/10.1016/j.mechmat.2007.03.006  



