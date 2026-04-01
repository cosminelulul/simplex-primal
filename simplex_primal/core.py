"""
simplex_primal.core
Pure-algorithm implementation of the Primal Simplex Method
using the Big-M technique and exact rational arithmetic.

All computation is done with Python's ``fractions.Fraction`` so
results are exact (no floating-point drift).
"""

import math
from fractions import Fraction

# Constants

#: Numeric value used to represent M (a very large penalty coefficient).
M_VALUE: Fraction = Fraction(10**7)

#: Positive infinity sentinel.
INF = math.inf


# Fraction helpers

def _frac(x) -> Fraction:
    """Convert *x* to a :class:`~fractions.Fraction`, clamping parse errors to 0."""
    if isinstance(x, Fraction):
        return x
    try:
        return Fraction(x).limit_denominator(10**9)
    except Exception:
        return Fraction(0)


def format_fraction(v: Fraction, *, show_M: bool = True) -> str:
    if show_M and abs(v) >= M_VALUE // 2:
        return "+M" if v > 0 else "-M"
    if v.denominator == 1:
        return str(v.numerator)
    return f"{v.numerator}/{v.denominator}"


def format_fraction_plain(v: Fraction) -> str:
    """Format a fraction without M-substitution (always shows the true value)."""
    if v.denominator == 1:
        return str(int(v))
    return f"{v.numerator}/{v.denominator}"


# Step 1 - Standardisation  (LP -> Standard LP)

def standardize(
    c,
    A,
    b,
    constraint_types,
    var_types,
    opt: str = "MAX",
):

    m = len(b)
    n = len(c)

    c_std = []
    A_std = [[] for _ in range(m)]
    var_names = []
    substitutions = {}

    #  Original decision variables (with sign substitutions) 
    for j in range(n):
        vtype = var_types[j] if j < len(var_types) else ">=0"
        if vtype == ">=0":
            c_std.append(_frac(c[j]))
            for i in range(m):
                A_std[i].append(_frac(A[i][j]))
            var_names.append(f"x{j+1}")

        elif vtype == "<=0":
            # Substitute x_j = -x_j'  (x_j' >= 0)
            c_std.append(-_frac(c[j]))
            for i in range(m):
                A_std[i].append(-_frac(A[i][j]))
            var_names.append(f"x{j+1}'")
            substitutions[j] = ("nonneg_negated", len(c_std) - 1)

        elif vtype == "R":
            # Substitute x_j = x_j^+ - x_j^-  (both >= 0)
            c_std.append(_frac(c[j]))
            c_std.append(-_frac(c[j]))
            for i in range(m):
                A_std[i].append(_frac(A[i][j]))
                A_std[i].append(-_frac(A[i][j]))
            var_names.append(f"x{j+1}")
            var_names.append(f"x{j+1}'")
            substitutions[j] = ("free", len(c_std) - 2, len(c_std) - 1)

    n_original = len(c_std)

    #  Collect slack / surplus / artificial variable info per row 
    row_info = []   # (row_index, ctype, slack_number_or_None, artif_number_or_None)
    slack_count = artif_count = 0
    for i in range(m):
        ctype = constraint_types[i]
        if ctype == "<=":
            slack_count += 1
            row_info.append((i, ctype, slack_count, None))
        elif ctype == ">=":
            slack_count += 1
            artif_count += 1
            row_info.append((i, ctype, slack_count, artif_count))
        elif ctype == "=":
            artif_count += 1
            row_info.append((i, ctype, None, artif_count))

    #  Add slack / surplus columns (y) 
    slack_col_map: dict[int, int] = {}
    for i, ctype, sn, _ in row_info:
        if ctype == "<=" and sn is not None:
            c_std.append(_frac(0))
            for k in range(m):
                A_std[k].append(_frac(1) if k == i else _frac(0))
            var_names.append(f"y{sn}")
            slack_col_map[i] = len(c_std) - 1
        elif ctype == ">=" and sn is not None:
            c_std.append(_frac(0))
            for k in range(m):
                A_std[k].append(_frac(-1) if k == i else _frac(0))
            var_names.append(f"y{sn}")
            slack_col_map[i] = len(c_std) - 1

    #  Add artificial columns (z) with M penalty 
    artif_col_map: dict[int, int] = {}
    for i, ctype, _, an in row_info:
        if an is not None:
            penalty = M_VALUE if opt == "MIN" else -M_VALUE
            c_std.append(_frac(penalty))
            for k in range(m):
                A_std[k].append(_frac(1) if k == i else _frac(0))
            var_names.append(f"z{an}")
            artif_col_map[i] = len(c_std) - 1

    #  Build initial basis 
    basis_indices = []
    for i in range(m):
        ctype = constraint_types[i]
        if ctype == "<=":
            basis_indices.append(slack_col_map[i])
        else:
            basis_indices.append(artif_col_map[i])

    b_std = [_frac(x) for x in b]
    return c_std, A_std, b_std, var_names, basis_indices, n_original, substitutions



# Step 2 - Tableau calculations

def compute_tableau_row(CB, XB, A, c, m, n):

    z_obj = sum(CB[i] * XB[i] for i in range(m))
    z_j = [sum(CB[i] * A[i][j] for i in range(m)) for j in range(n)]
    delta = [c[j] - z_j[j] for j in range(n)]
    return z_obj, z_j, delta


def pivot_step(A, XB, basis, CB, c, pivot_row, pivot_col, m, n):
    P = A[pivot_row][pivot_col]
    XB[pivot_row] /= P
    for j in range(n):
        A[pivot_row][j] /= P
    for i in range(m):
        if i == pivot_row:
            continue
        factor = A[i][pivot_col]
        XB[i] -= factor * XB[pivot_row]
        for j in range(n):
            A[i][j] -= factor * A[pivot_row][j]
    basis[pivot_row] = pivot_col
    CB[pivot_row] = c[pivot_col]


# Step 3 - Solution display / verification

def format_solution(
    c_orig,
    A_orig,
    b_orig,
    x_opt,
    f_opt,
    opt: str,
    var_names_orig,
    basis_final,
    XB_final,
    var_names_std,
    A_tableau_initial,
    n_original: int,
) -> str:
    lines = []
    orig_n = len(c_orig)
    m = len(b_orig)

    lines.append("=" * 55)
    lines.append(" Optimal Solution:")
    lines.append("=" * 55)

    for j, xv in enumerate(x_opt):
        lines.append(f" x{j+1}* = {format_fraction_plain(xv)}")

    lines.append("")
    lines.append("-" * 55)
    lines.append(f" Objective function ({opt}):")
    lines.append(f" f* = {format_fraction_plain(f_opt)}")
    lines.append("")

    terms_sym = [f"{format_fraction_plain(_frac(c_orig[j]))}*x{j+1}" for j in range(orig_n)]
    lines.append(f" f = {' +'.join(terms_sym)}")

    terms_val = [
        f"{format_fraction_plain(_frac(c_orig[j]))}*{format_fraction_plain(_frac(x_opt[j]))}"
        for j in range(orig_n)
    ]
    lines.append(f" f = {' +'.join(terms_val)}")

    products = []
    total = _frac(0)
    for j in range(orig_n):
        prod = _frac(c_orig[j]) * _frac(x_opt[j])
        products.append(format_fraction_plain(prod))
        total += prod
    lines.append(f" f = {' +'.join(products)}")
    lines.append(f" f_{opt.lower()}* = {format_fraction_plain(_frac(f_opt))}")

    lines.append("")
    lines.append("-" * 55)
    lines.append(" Verification: b = S · X_B")
    lines.append("")

    # Build the basis matrix S from the initial tableau
    S = []
    for i in range(m):
        row = []
        for bi in basis_final:
            if bi < len(A_tableau_initial[i]):
                row.append(A_tableau_initial[i][bi])
            else:
                row.append(_frac(0))
        S.append(row)

    XB_f = [_frac(xv) for xv in XB_final]
    b_f  = [_frac(bv) for bv in b_orig]
    S_XB = [sum(S[i][k] * XB_f[k] for k in range(m)) for i in range(m)]

    w_b   = max(len(format_fraction_plain(v)) for v in b_f) + 2
    w_s   = max(len(format_fraction_plain(S[i][k])) for i in range(m) for k in range(m)) + 2
    w_xb  = max(len(format_fraction_plain(v)) for v in XB_f) + 2
    w_sxb = max(len(format_fraction_plain(v)) for v in S_XB) + 2

    b_col   = [format_fraction_plain(v).center(w_b)   for v in b_f]
    xb_col  = [format_fraction_plain(v).center(w_xb)  for v in XB_f]
    sxb_col = [format_fraction_plain(v).center(w_sxb) for v in S_XB]
    s_rows  = [
        "  ".join(format_fraction_plain(S[i][k]).center(w_s) for k in range(m))
        for i in range(m)
    ]

    detail_rows = []
    for i in range(m):
        parts  = " + ".join(
            f"{format_fraction_plain(S[i][k])}·{format_fraction_plain(XB_f[k])}"
            for k in range(m)
        )
        result = format_fraction_plain(S_XB[i])
        mark   = "✅" if S_XB[i] == b_f[i] else "❌"
        detail_rows.append(f"  {parts} = {result}  {mark}")

    for i in range(m):
        prefix   = "(" if i in (0, m - 1) else "|"
        suffix   = ")" if i in (0, m - 1) else "|"
        eq_sign  = "==" if i == m // 2 else "  "
        mul_sign = "·"  if i == m // 2 else " "
        lines.append(
            f"  {prefix}{b_col[i]}{suffix} {eq_sign} "
            f"{prefix}{s_rows[i]}{suffix} {mul_sign} "
            f"{prefix}{xb_col[i]}{suffix}"
        )

    lines.append("")
    lines.append("Calculation details:")
    lines.extend(detail_rows)

    all_ok = all(S_XB[i] == b_f[i] for i in range(m))
    lines.append("")
    lines.append(f"  Verification result: {'✅ CORRECT' if all_ok else '❌ ERROR'}")
    lines.append("=" * 55)

    return "\n".join(lines)


# Main solver

def solve(
    c,
    A,
    b,
    constraint_types,
    var_types,
    opt: str = "MAX",
    max_iterations: int = 200,
) -> dict:

    # Convert everything to Fraction
    c_F = [_frac(x) for x in c]
    A_F = [[_frac(x) for x in row] for row in A]
    b_F = [_frac(x) for x in b]

    (c_std, A_std, b_std,
     var_names, basis, n_original,
     substitutions) = standardize(c_F, A_F, b_F, constraint_types, var_types, opt)

    m = len(b_std)
    n = len(c_std)

    tableau      = [row[:] for row in A_std]
    tableau_init = [row[:] for row in A_std]   # kept for verification
    XB           = b_std[:]
    CB           = [c_std[basis[i]] for i in range(m)]

    iterations: list[dict] = []
    status = "optimal"

    for k in range(max_iterations + 1):
        z_obj, z_j, delta = compute_tableau_row(CB, XB, tableau, c_std, m, n)

        snap: dict = {
            "k"         : k,
            "var_names" : var_names[:],
            "basis"     : basis[:],
            "CB"        : CB[:],
            "XB"        : XB[:],
            "A"         : [row[:] for row in tableau],
            "c"         : c_std[:],
            "z_obj"     : z_obj,
            "z_j"       : z_j[:],
            "delta"     : delta[:],
            "pivot_col" : None,
            "pivot_row" : None,
        }

        # Optimality check
        if opt == "MAX":
            is_optimal = all(d <= _frac(0) for d in delta)
        else:
            is_optimal = all(d >= _frac(0) for d in delta)

        if is_optimal:
            snap["status"] = "optimal"
            iterations.append(snap)
            break

        # Entering variable (most positive / most negative Δ)
        if opt == "MAX":
            candidates = [j for j in range(n) if delta[j] > _frac(0)]
            pivot_col  = max(candidates, key=lambda j: delta[j])
        else:
            candidates = [j for j in range(n) if delta[j] < _frac(0)]
            pivot_col  = min(candidates, key=lambda j: delta[j])

        # Unboundedness check
        if all(tableau[i][pivot_col] <= _frac(0) for i in range(m)):
            snap["status"] = "unbounded"
            iterations.append(snap)
            status = "unbounded"
            break

        # Minimum-ratio test (leaving variable)
        ratios = [
            (XB[i] / tableau[i][pivot_col], i)
            for i in range(m)
            if tableau[i][pivot_col] > _frac(0)
        ]
        if not ratios:
            snap["status"] = "unbounded"
            iterations.append(snap)
            status = "unbounded"
            break

        _, pivot_row = min(ratios, key=lambda x: x[0])

        snap["pivot_col"] = pivot_col
        snap["pivot_row"] = pivot_row
        snap["status"]    = "continue"
        iterations.append(snap)

        pivot_step(tableau, XB, basis, CB, c_std, pivot_row, pivot_col, m, n)

    else:
        status = "max_iter"

    #  Recover original variable values 
    x_full = [_frac(0)] * n
    for i, bi in enumerate(basis):
        x_full[bi] = XB[i]

    x_std  = x_full[:n_original]
    x_opt  = []
    std_idx = 0
    orig_n  = len(c_F)

    for j in range(orig_n):
        vtype = var_types[j] if j < len(var_types) else ">=0"
        if vtype == ">=0":
            x_opt.append(x_std[std_idx]);       std_idx += 1
        elif vtype == "<=0":
            x_opt.append(-x_std[std_idx]);      std_idx += 1
        elif vtype == "R":
            x_opt.append(x_std[std_idx] - x_std[std_idx + 1])
            std_idx += 2

    z_final, _, _ = compute_tableau_row(CB, XB, tableau, c_std, m, n)

    sol_text = ""
    if status == "optimal":
        sol_text = format_solution(
            c_F, A_F, b_F,
            x_opt, z_final, opt,
            [f"x{j+1}" for j in range(orig_n)],
            basis, XB, var_names,
            tableau_init, n_original,
        )

    return {
        "status"     : status,
        "iterations" : iterations,
        "x_opt"      : x_opt,
        "f_opt"      : z_final,
        "sol_text"   : sol_text,
        "var_names"  : var_names,
        "n_original" : n_original,
        "orig_n"     : orig_n,
    }
