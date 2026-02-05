# 📡 MIMO BC Weighted Sum Rate Optimization via MAC–BC Duality  
### Paper Implementation & Verification Project

This repository contains an educational implementation and verification of a published algorithm for **weighted sum rate (WSR) maximization** in multi-user MIMO broadcast channels using **MAC–BC duality transformation**.

The original non-convex Broadcast Channel (BC) optimization is transformed into an equivalent convex Multiple Access Channel (MAC) problem and solved using an iterative dual covariance update algorithm.

> ⚠️ This is a reproduction and verification project based on prior published work. All theoretical and algorithmic credit belongs to the original authors (see References section).

---

## 🎯 Purpose

- Implement MAC–BC duality based WSR optimization algorithm from literature
- Reproduce iterative dual MAC covariance update procedure
- Verify correctness using MISO special-case comparison
- Validate MIMO performance trends against textbook references
- Analyze convergence behavior across antenna/user configurations

This repo focuses on **correct implementation, validation, and analysis**, not proposing a new algorithm.

---

## 🧠 Problem Background

Weighted sum rate maximization in MIMO broadcast channels is non-convex due to inter-user interference coupling. MAC–BC duality transforms the BC problem into a convex MAC domain problem with the same sum power constraint. Optimization can then be performed over MAC covariance matrices using an iterative procedure.

---

## ⚙️ Algorithm Overview

Iterative Dual MAC Covariance Optimization:

1. Initialize user covariance matrices to zero  
2. Compute gradient matrices for each user  
3. Extract dominant eigenvalue/eigenvector  
4. Select best user update direction  
5. Solve 1-D step size optimization  
6. Update covariance matrices under sum power constraint  
7. Repeat until convergence tolerance is met  

Convergence is declared when weighted sum rate improvement falls below a threshold.

---

## ✅ Verification Methodology

### MISO Cross-Validation
- Set receive antennas **N = 1**
- Reduce system to MISO BC
- Compare results with established MISO solutions
- Mean absolute percentage error < 1%

### MIMO Validation
- Vary receive antennas (N = 1–4)
- Keep transmit antennas and users fixed
- Compare trends with textbook reference curves
- Verify diversity gain and scaling behavior

---

## 📊 Results Summary

- Weighted sum rate increases monotonically across iterations
- Convergence typically within **10–20 iterations**
- MISO verification error < 1%
- Sum rate increases with receive antennas
- Diminishing returns observed (theory-consistent)
- Stable convergence across antenna/user configurations

---

## 📈 Experiments Included

- WSR vs transmit power
- WSR vs antenna configuration
- Receive antenna scaling study
- Multi-user scenarios
- Convergence vs iterations
- MISO vs MIMO verification

---

## 🛠 Requirements

- Python with NumPy / SciPy

---
## ▶️ How to Run the Code

This project is implemented in Python with the main experiments driven through a Jupyter notebook and a supporting algorithm module.

### 📦 Files

- `[MIMO_BC.ipynb](MIMO_BC.ipynb)` — main notebook for simulations and experiments  
- `[opt_dual_mac_cov.py](opt_dual_mac_cov.py)` — implementation of the iterative dual MAC covariance optimization algorithm  

---

### ✅ Environment Setup

Create a Python environment (recommended):

```bash
python -m venv venv
source venv/bin/activate   # Linux/macOS
venv\Scripts\activate      # Windows
```

Install required packages:

```bash
pip install numpy scipy matplotlib jupyter
```

---

### ▶️ Run via Jupyter Notebook (Recommended)

Start Jupyter:

```bash
jupyter notebook
```

Open:

```
MIMO_BC.ipynb
```

Run all cells to:

- generate random channel realizations  
- execute MAC–BC duality optimization  
- compute weighted sum rate  
- run MISO verification tests  
- plot convergence and scaling results  

All experiment parameters are configurable inside the notebook:
- transmit antennas (M)
- receive antennas (N)
- users (K)
- power levels
- iteration tolerance

---

### ▶️ Run Algorithm Module Directly (Optional)

If you want to use the optimizer directly in your own script:

```python
from opt_dual_mac_cov import opt_dual_mac_cov

Q_opt = opt_dual_mac_cov(H_list, P_total, noise_var, weights, tol, max_iter)
```

Where:
- `H_list` — list of user channel matrices  
- `P_total` — sum power constraint  
- `noise_var` — noise variance  
- `weights` — user rate weights  
- `tol` — convergence tolerance  
- `max_iter` — iteration cap  

---

### 📊 Reproducing Report Results

The notebook includes experiment blocks for:

- MISO cross-check validation  
- MIMO antenna scaling study  
- WSR vs power curves  
- Convergence vs iterations plots  

Run all cells sequentially to reproduce the reported figures.

---

## 📄 Project Report

Detailed derivations, formulation, and validation methodology are documented in:

[Project Report](mini_project_report.pdf)

---

## 🙏 Credits & References

This implementation is based on MAC–BC duality and weighted sum rate optimization methods described in:

Howard Huang, Constantinos B. Papadias, and Sivarama Venkatesan,  
**MIMO Communication for Cellular Networks**, Springer, 2012.

Harish Viswanathan, Sivarama Venkatesan, and Howard Huang,  
“Downlink capacity evaluation of cellular networks with known-interference cancellation,”  
IEEE JSAC, 2003.

All algorithmic credit belongs to the original authors.  
This repository provides an independent educational implementation and verification.

---

## 📚 BibTeX Citations

```bibtex
@book{huang2012mimo,
  title={MIMO Communication for Cellular Networks},
  author={Huang, Howard and Papadias, Constantinos B. and Venkatesan, Sivarama},
  year={2012},
  publisher={Springer}
}

@article{viswanathan2003downlink,
  title={Downlink capacity evaluation of cellular networks with known-interference cancellation},
  author={Viswanathan, Harish and Venkatesan, Sivarama and Huang, Howard},
  journal={IEEE Journal on Selected Areas in Communications},
  volume={21},
  number={5},
  pages={802--811},
  year={2003}
}
```

---

## ⚖️ Academic Use Notice

This code is provided for learning and research purposes only.

If you use this repository:

- Cite the original references above  
- Do not cite this repo as the original algorithm source  
- Clearly mark usage as an implementation/reproduction

---

## 👨‍💻 Author Note

This project was completed as part of an academic study to understand and validate MAC–BC duality based optimization algorithms for multi-user MIMO systems.
