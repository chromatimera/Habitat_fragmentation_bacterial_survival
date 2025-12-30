# Habitat_fragmentation_bacterial_survival

## 🌟 Overview

This repository contains the Python code used to generate the simulations and figures for:

**Verdon N, Popescu O, Titmuss S, Allen RJ (2025)**  
*Habitat fragmentation enhances microbial collective defence*  
*Journal of the Royal Society Interface* 22(223):20240611  
https://doi.org/10.1098/rsif.2024.0611

The project studies how **spatial fragmentation of microbial habitats** (e.g. droplets or isolated micro-environments) can *increase population survival* under antibiotic stress through stochastic effects and local population heterogeneity.

---

## 📌 Key Concepts

- Habitat fragmentation partitions a population into many isolated subvolumes  
- Demographic noise and uneven initial loading can increase survival probability  
- Fragmentation can enhance collective defence for enzyme-producing microbes  
- Deterministic and stochastic effects can be independently controlled and compared

---

## ⚙️ Simulation Framework

Simulations are implemented in **Python** using a modular, class-based design:

- **`strain`**  
  Implements microbial growth and antibiotic degradation dynamics

- **`Experiment_R`**  
  Stores population trajectories, antibiotic concentration, and numerical settings

- **`Droplet_R`**  
  Simulates multiple isolated subvolumes (droplets), computes survival statistics, and aggregates results

### Supported simulation modes

The framework allows independent control of:

- **Population dynamics**
  - Deterministic (ODE-based)
  - Stochastic (Gillespie birth–death)

- **Initial population loading**
  - Deterministic
  - Stochastic (Poisson partitioning)

These options yield **four simulation scenarios**, enabling clean separation of:
- stochasticity from demographic noise  
- stochasticity from spatial partitioning

---

## 🧪 Fragmentation Model

- Total volume \( V \) is divided into \( m \) isolated subvolumes  
- No exchange of microbes or antibiotic between subvolumes  
- For stochastic partitioning, cells are assigned using a Poisson distribution  
- Subvolumes are simulated independently and aggregated for population-level outcomes

---

## 📊 Output & Analysis

- Each condition is typically repeated **1000 times** to estimate survival probabilities  
- Survival is defined as at least one subvolume remaining non-empty at **t = 300 min**
- Outputs can be saved as:
  - full time series
  - endpoint statistics

All post-processing and visualization is performed in Python using **Matplotlib**, producing:
- survival curves
- phase diagrams
- heatmaps
- individual stochastic trajectories

All figures in the associated paper are generated from reproducible scripts included here.

---

## 🔁 Reproducibility

- Deterministic simulations use numerical integration of coupled ODEs  
- Stochastic simulations use a **hybrid Gillespie–ODE approach**  
- Parameters, numerical settings, and random seed handling are documented in code  

---

## 📦 Code Availability

A frozen, citable version of this codebase is archived on **Zenodo**:

**DOI:** https://doi.org/10.5281/zenodo.14514221

The code is free to reuse or adapt with appropriate citation.

---

## 📬 Contact

For questions or reuse, please don't hesitate to get in touch.




