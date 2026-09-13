# Deep-BioShield: Multiscale Computational Design for Arsenic Interception

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview
This repository contains the complete computational framework and source code for the manuscript: **"Multiscale Computational Design of an Engineered Protein (Deep-BioShield) for Arsenic Interception in Rice Rhizospheres: A Proof-of-Concept Framework"**.

The pipeline establishes a deterministic, physics-based multiscale architecture spanning generative protein design, coordinate chemisorption thermodynamics, metabolic flux, and macroscopic reactive soil transport.

## Repository Structure & Figure Mapping

The codebase is modularized to reproduce the 7 primary figures from the manuscript:

*   `00_Data_Acquisition.py`: Prerequisite script to mine the UniProt REST API for authentic ArsR, ArsC, and Acr3 sequences.
*   `01_Generative_Design.py`: Reconstructs the ESM-2 latent space (Figure 1A) and the constrained Genetic Algorithm optimization trajectory (Figure 1B).
*   `02_Biosafety_and_Folding.py`: Computes physicochemical ecotoxicity (GRAVY/Charge) and parses ESMFold pLDDT scores (Figure 2B, 2C).
*   `03_FBA_Metabolism.py`: Executes Flux Balance Analysis (FBA) using `cobrapy` to calculate genome-scale metabolic burden (Figure 3).
*   `04_Thermodynamics_Kinetics.py`: Deterministic modeling of the standard-state Langmuir chemisorption isotherm and pseudo-second-order reaction kinetics (Figure 4).
*   `05_ADRE_Transport_and_Health.py`: Numerical finite-difference solver for the 1D Advection-Dispersion-Reaction Equation (ADRE) and USEPA Lifetime Cancer Risk (LCR) (Figure 5).
*   `06_Lotka_Volterra_Dynamics.py`: Solves tri-trophic ODEs for rhizosphere systems ecology (Figure 6C).
*   `07_Genetic_Circuit.py`: Generates the pBioShield-v1 synthetic biology plasmid map (Figure 7).
*   `DeepBioShield_Evolved.pdb`: The atomic coordinates of the engineered 68-AA, 8-cysteine monomer.

## Installation & Prerequisites

1. Clone the repository:
   ```bash
   git clone [https://github.com/Imran-AMRClinical-Rep/Precision-Microbiome-MCPINN.git](https://github.com/Imran-AMRClinical-Rep/Precision-Microbiome-MCPINN.git)
   cd Precision-Microbiome-MCPINN
