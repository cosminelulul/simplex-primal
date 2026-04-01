# simplex-primal

A pure-Python implementation of the **Primal Simplex Method** for Linear Programming.

Available on PyPi: [simplex-primal 1.0](https://pypi.org/project/simplex-primal/)

---

## Features:
- ✅ Exact rational arithmetic - no floating-point drift, thanks to [fractions.Fraction](https://docs.python.org/3/library/fractions.html)
- ✅ Big-M method for artificial variables ($>=$ and $=$ constraints)
- ✅ Supports MIN and MAX objectives
- ✅ Handles $<=0$ and unrestricted ($R$) variables via substitution
- ✅ Full iteration-by-iteration snapshots for teaching / debugging
- ✅ Zero third-party dependencies (stdlib only)

---

## Quick Start:

```python
from simplex_primal import solve

# Maximise f = 5x1 + 4x2
# subject to:
#   6x1 + 4x2 <= 24
#    x1 + 2x2 <= 6
#   x1, x2 >= 0

result = solve(
    c=[5, 4],
    A=[[6, 4], [1, 2]],
    b=[24, 6],
    constraint_types=["<=", "<="],
    var_types=[">=0", ">=0"],
    opt="MAX",
)

print(result["status"])        # "optimal"
print(result["x_opt"])         # [Fraction(3,1), Fraction(3,2)]
print(result["f_opt"])         # Fraction(21,1)
print(result["sol_text"])      # formatted solution + verification
```
---

## Installation
```bash
pip install simplex-primal

pip install git+https://github.com/cosminelulul/simplex-primal.git
```
--- 
## License
This project is licensed under the [MIT License](LICENSE).

---


## References:
Coursework from Faculty of Applied Sciences, University POLITEHNICA of Bucharest.
1. [Curs_1_CO - PL(1) (S. Bibic).pdf](https://curs.upb.ro/2025/pluginfile.php/268349/mod_folder/content/0/Curs%201%20CO/Curs_1_CO%20-%20PL%281%29%20%5BS.%20Bibic%5D.pdf)
2. [Curs_2_CO - PL(2) (S. Bibic).pdf](https://curs.upb.ro/2025/pluginfile.php/268492/mod_folder/content/0/Curs%202%20CO/Curs_2_CO%20-%20PL%282%29%20%5BS.%20Bibic%5D.pdf)
3. [Cursul_3_CO - PL(3) (S. Bibic).pdf](https://curs.upb.ro/2025/pluginfile.php/268497/mod_folder/content/0/Curs%203%20CO/Cursul_3_CO%20-%20PL%283%29%20%5BS.%20Bibic%5D.pdf?)
4. [Seminar_1_CO - PL(1) (S. Bibic).pdf](https://curs.upb.ro/2025/pluginfile.php/268349/mod_folder/content/0/Seminar%201%20CO/Seminar_1_CO%20-%20PL%281%29%20%5BS.%20Bibic%5D.pdf)
5. [Seminar_2_CO - PL(2) (S. Bibic).pdf](https://curs.upb.ro/2025/pluginfile.php/268492/mod_folder/content/0/Seminar%202%20CO/Seminar_2_CO%20-%20PL%282%29%20%5BS.%20Bibic%5D.pdf)
