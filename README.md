[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.21202864-blue)](https://doi.org/10.5281/zenodo.21202864) [![License: CC BY-NC-ND 4.0](https://img.shields.io/badge/License-CC%20BY--NC--ND%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-nd/4.0/)

## UPDSF v4.4: Unified Prebiotic DNA Selection Framework Strictly Empirical Simulation of Chemical Darwinism

Environment: 🐍 Python 3.8+

## Author: Seyed Mohammad Reza Hashemi (Reza Hashemi)

## ORCID:0009-0002-0645-5180

https://doi.org/10.5281/zenodo.17273763

https://doi.org/10.5281/zenodo.18137476

https://doi.org/10.5281/zenodo.18092867

https://doi.org/10.5281/zenodo.18080826

https://doi.org/10.5281/zenodo.20988680

https://doi.org/10.5281/zenodo.20825578

https://doi.org/10.5281/zenodo.20733760

https://doi.org/10.5281/zenodo.20759622

https://doi.org/10.5281/zenodo.20771213

https://doi.org/10.5281/zenodo.18594133

https://doi.org/10.5281/zenodo.17280179

https://doi.org/10.5281/zenodo.17964430


## Overview
UPDSF v4.4 is a high-fidelity computational engine designed to model the chemical selection and evolutionary dominance of DNA nucleotides (specifically Thymine) under prebiotic conditions. Developed by Seyed Mohammad Reza Hashemi, this framework operates under the Matter World Hypothesis (MWH) and represents a case study in Intelligence-Augmented (IA) Science.

Unlike theoretical models, v4.4 is strictly empirical, utilizing kinetic parameters, activation energies ($E_a$), and half-lives derived exclusively from peer-reviewed prebiotic chemistry literature.


## The Scientific Core: 

Non-Enzymatic Transition (U → T)
The central breakthrough of UPDSF v4.4 is demonstrating the spontaneous chemical selection of Thymine over Uracil in prebiotic DNA. 

### Chemical Darwinism vs. Biological Enzymes

In modern biology, the presence of Thymine in DNA and its distinction from Uracil is maintained by complex enzymes (like Uracil-DNA Glycosylase). However, the Matter World Hypothesis (MWH)* argues that this selection occurred long before the emergence of proteins. 

This framework proves through *Chemical Darwinism that:
1.  Thermal Selection: Under high-temperature prebiotic conditions (55°C - 85°C), Uracil-bearing polymers suffer from a lower activation energy barrier ($E_a \approx 27$ kcal/mol) for hydrolysis, leading to rapid informational decay.
2.  The Methyl Advantage: The addition of the methyl group in Thymine ($E_a \approx 32$ kcal/mol) provides a 5 kcal/mol "stability premium." In a molecular mixture, this acts as a natural selection filter, where T-DNA strands survive and accumulate while U-DNA strands are chemically "extinguished."
3.  Spontaneous Code Purification: The high rate of Cytosine-to-Uracil deamination (36x faster than other mutations) creates a "noise crisis" in early genetic systems. The simulation shows that the transition to Thymine was a physicochemical necessity to distinguish permanent genetic information from ongoing chemical noise, achieving "code purity" without a single enzyme.

### Empirical Evidence in Silico

By using strictly empirical Arrhenius rates, UPDSF v4.4 shows that the "Thymine Dominance" is not a biological accident, but a deterministic outcome of prebiotic thermodynamics.

## Key Scientific Features

- Empirical Arrhenius Kinetics: Calibrated using Lindahl (1993) and Cleaves (2010) to model hydrolysis and deamination.
- Deamination Modeling: Includes the verified 36x higher deamination rate of Cytosine to Uracil (Shen et al., 1994).
- UV Photostability: Models the 3-4x higher UV resistance of Thymine compared to Uracil (Ravanat & Cadet, 1995).
- Polymer Physics: Incorporates persistence length, contour length, and conformational dynamics for long-chain stability.
- Langevin Dynamics: Simulates Brownian forces and thermal fluctuations to model real-world environmental turbulence.
- 2D Sensitivity Analysis: Optimization of Nucleotide Enrichment across a Temperature × pH gradient.

## Core Literature Calibration
| Parameter | Source | Value/Ratio |
| :--- | :--- | :--- |
| DNA Stability | Lindahl (1993) | $E_a \approx 32$ kcal/mol |
| Cytosine Deamination | Shapiro (1999) | pH-dependent curves |
| Prebiotic Clay Catalysis | Ferris (1996) | Enhanced polymerization rates |
| UV Resistance (T/U) | Ravanat (1995) | ~3.6x advantage for T |

## 
## Installation & Setup

To run the UPDSF v4.4 simulation on your local machine, follow these steps:

1. Clone the repository:
   git clone https://github.com/mrhashemi2000/UPDSF_v4.4.py.git
   cd UPDSF_v4.4
   
2. Install the dependencies:
   It is recommended to use a virtual environment. You can install all required libraries using the provided `requirements.txt` file:
   pip install -r requirements.txt
   
3. Run the simulation:
   python UPDSF_v4.4.py

## Methodology: IA-Augmented Discovery
This project utilizes a recursive collaboration between human reasoning and AI-assisted modeling. The framework evolves through continuous feedback loops, allowing for rapid interdisciplinary discovery at the intersection of astrobiology, physical chemistry, and computational physics.

## Citation
If you use this framework in your research, please cite it as:
> Reza Hashemi. (2026). mrhashemi2000/UPDSF_v4.4.py: Initial release. Zenodo. https://doi.org/10.5281/zenodo.21202864

## REFERENCES:
    - Lindahl, T. (1993). Nature, 362, 709-715. (DNA stability)
    - Shapiro, R. (1999). Chem. Rev., 99, 2501-2536. (Deamination)
    - Cleaves, H.J. (2010). Astrobiology, 10, 337-346. (Prebiotic chemistry)
    - Ravanat, J.L. (1995). J. Biol. Chem., 270, 12305-12311. (UV damage)
    - Shen, J.C. (1994). Biochemistry, 33, 10756-10764. (Cytosine deamination)
    - Ferris, J.P. (1996). Orig. Life Evol. Biosph., 26, 449-461. (Clay catalysis)


