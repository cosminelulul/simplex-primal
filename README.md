# simplex-primal

A pure-Python implementation of the **Primal Simplex Method** for Linear Programming.

---

## Features:
- ✅ Exact rational arithmetic - no floating-point drift, thanks to [fractions.Fraction](https://docs.python.org/3/library/fractions.html)
- ✅ Big-M method for artificial variables ($>=$ and $=$ constraints)
- ✅ Supports MIN and MAX objectives
- ✅ Handles $<=0$ and unrestricted ($R$) variables via substitution
- ✅ Full iteration-by-iteration snapshots for teaching / debugging
- ✅ Zero third-party dependencies (stdlib only)

---

## Installation
```bash

pip install git+https://github.com/cosminelulul/simplex-primal.git
```

---
