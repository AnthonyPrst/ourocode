# coding in UTF-8
# by Anthony PARISOT
"""Remplacement de handcalcs.decorator.handcalc compatible Python 3.13+.
 
Ce module implémente un décorateur @handcalc qui génère du LaTeX à partir
du code source Python via l'AST (ast module), sans dépendance à innerscope.
Il reproduit le comportement utilisé dans ourocode : rendu des assignations
d'une fonction interne avec substitution symbolique et numérique.
 
Retour du décorateur : tuple (latex_str, valeur_retournée)
  - latex_str  : chaîne LaTeX prête à injecter dans IPython.display.Latex
  - valeur     : résultat numérique de la fonction (scalaire, tuple, dict…)
"""
 
import ast
import inspect
import re
import textwrap
from math import sqrt, pi, sin, cos, tan, radians, log, log10, exp, floor, ceil, fabs
from typing import Any
 
 
# ---------------------------------------------------------------------------
# Helpers de conversion AST → LaTeX
# ---------------------------------------------------------------------------
 
_GREEK = {
    "alpha": r"\alpha", "beta": r"\beta", "gamma": r"\gamma", "delta": r"\delta",
    "epsilon": r"\epsilon", "zeta": r"\zeta", "eta": r"\eta", "theta": r"\theta",
    "lambda": r"\lambda", "mu": r"\mu", "nu": r"\nu", "xi": r"\xi",
    "pi": r"\pi", "rho": r"\rho", "sigma": r"\sigma", "tau": r"\tau",
    "phi": r"\phi", "chi": r"\chi", "psi": r"\psi", "omega": r"\omega",
    "Gamma": r"\Gamma", "Delta": r"\Delta", "Theta": r"\Theta",
    "Lambda": r"\Lambda", "Xi": r"\Xi", "Pi": r"\Pi", "Sigma": r"\Sigma",
    "Phi": r"\Phi", "Psi": r"\Psi", "Omega": r"\Omega",
}
 
_FUNC_MAP = {
    "sqrt": r"\sqrt", "sin": r"\sin", "cos": r"\cos", "tan": r"\tan",
    "log": r"\log", "log10": r"\log_{10}", "exp": r"\exp",
    "abs": r"\left|", "fabs": r"\left|",
    "floor": r"\lfloor", "ceil": r"\lceil",
    "radians": r"\frac{\pi}{180}",
    "min": r"\min", "max": r"\max",
}
 
 
def _name_to_latex(name: str) -> str:
    """Convertit un nom Python en symbole LaTeX (grecs, indices)."""
    if name in _GREEK:
        return _GREEK[name]
    # Gestion des sous-indices : taux_6_11 → taux_{6,11} ; sigma_m_d → \sigma_{m,d}
    parts = name.split("_")
    if len(parts) == 1:
        return name
    base = parts[0]
    subs = parts[1:]
    base_latex = _GREEK.get(base, base)
    sub_str = ",".join(subs)
    return f"{base_latex}_{{{sub_str}}}"
 
 
def _is_physical(value: Any) -> bool:
    """Retourne True si la valeur est un objet forallpeople.Physical."""
    return hasattr(value, "latex") and hasattr(value, "dimensions")


def _num_to_latex(value: Any, precision: int) -> str:
    """Convertit une valeur numérique en chaîne LaTeX (sans unité)."""
    try:
        if _is_physical(value):
            # str(x) -> "1.500 MPa" : on extrait la partie numérique affichée
            s = str(value)
            parts = s.rsplit(" ", 1)
            num_str = parts[0].strip()
            try:
                v = float(num_str)
                if abs(v) >= 1e6 or (abs(v) < 1e-3 and v != 0):
                    return f"{v:.{precision}e}"
                return f"{round(v, precision)}"
            except ValueError:
                return num_str
        v = value
        if isinstance(v, int):
            return str(v)
        if isinstance(v, float):
            if abs(v) >= 1e6 or (abs(v) < 1e-3 and v != 0):
                return f"{v:.{precision}e}"
            return f"{round(v, precision)}"
        return str(round(float(v), precision))
    except Exception:
        return str(value)
 
 
class _ASTNumRenderer(ast.NodeVisitor):
    """Visite l'AST d'une expression et substitue les Name par leurs valeurs numériques.

    Produit la forme numérique intermédiaire : \\min(0.77 + 3.5/25, 1).
    """

    def __init__(self, local_vars: dict, precision: int):
        self.local_vars = local_vars
        self.precision = precision

    def render(self, node: ast.AST) -> str:
        return self.visit(node)

    def visit_Name(self, node: ast.Name) -> str:
        name = node.id
        if name in self.local_vars:
            val = self.local_vars[name]
            num = _num_to_latex(val, self.precision + 1)
            unit = _get_unit_str(val)
            return f"{num}{unit}" if unit else num
        return _name_to_latex(name)

    def visit_Constant(self, node: ast.Constant) -> str:
        if isinstance(node.value, float):
            return f"{node.value}"
        return str(node.value)

    def visit_BinOp(self, node: ast.BinOp) -> str:
        sym = _ASTRenderer(self.local_vars, self.precision)
        left = self.visit(node.left)
        right = self.visit(node.right)
        op = node.op
        if isinstance(op, ast.Add):
            return f"{left} + {right}"
        if isinstance(op, ast.Sub):
            return f"{left} - {right}"
        if isinstance(op, ast.Mult):
            return f"{left} \\cdot {right}"
        if isinstance(op, ast.Div):
            return f"\\frac{{{left}}}{{{right}}}"
        if isinstance(op, ast.Pow):
            base = left if _is_simple(node.left) else f"\\left({left}\\right)"
            exp_str = right if _is_simple(node.right) else f"{{{right}}}"
            return f"{base}^{{{exp_str}}}"
        if isinstance(op, ast.FloorDiv):
            return f"\\lfloor\\frac{{{left}}}{{{right}}}\\rfloor"
        return f"{left} \\circ {right}"

    def visit_UnaryOp(self, node: ast.UnaryOp) -> str:
        operand = self.visit(node.operand)
        if isinstance(node.op, ast.USub):
            return f"-{operand}"
        return operand

    def visit_Call(self, node: ast.Call) -> str:
        func_name = ""
        if isinstance(node.func, ast.Name):
            func_name = node.func.id
        elif isinstance(node.func, ast.Attribute):
            func_name = node.func.attr
        args = [self.visit(a) for a in node.args]
        if func_name == "sqrt":
            return f"\\sqrt{{{args[0]}}}" if args else r"\sqrt{\cdot}"
        if func_name in ("sin", "cos", "tan"):
            return f"\\{func_name}\\left({args[0]}\\right)" if args else f"\\{func_name}()"
        if func_name in ("min", "max"):
            return f"\\{func_name}\\left({', '.join(args)}\\right)"
        if func_name in ("abs", "fabs"):
            return f"\\left|{args[0]}\\right|" if args else r"\left|\cdot\right|"
        if func_name == "log":
            return f"\\ln\\left({args[0]}\\right)" if args else r"\ln(\cdot)"
        if func_name == "exp":
            return f"e^{{{args[0]}}}" if args else "e^{\\cdot}"
        if func_name == "radians":
            return f"\\frac{{\\pi \\cdot {args[0]}}}{{180}}" if args else ""
        if func_name == "floor":
            return f"\\lfloor {args[0]} \\rfloor" if args else ""
        if func_name == "ceil":
            return f"\\lceil {args[0]} \\rceil" if args else ""
        return f"{_name_to_latex(func_name)}\\left({', '.join(args)}\\right)"

    def visit_Compare(self, node: ast.Compare) -> str:
        parts = [self.visit(node.left)]
        ops_map = {ast.Lt: "<", ast.LtE: "\\leq", ast.Gt: ">", ast.GtE: "\\geq",
                   ast.Eq: "=", ast.NotEq: "\\neq"}
        for op, comp in zip(node.ops, node.comparators):
            parts.append(ops_map.get(type(op), "?"))
            parts.append(self.visit(comp))
        return " ".join(parts)

    def visit_IfExp(self, node: ast.IfExp) -> str:
        body = self.visit(node.body)
        test = self.visit(node.test)
        orelse = self.visit(node.orelse)
        return f"{body} \\text{{ si }} {test} \\text{{ sinon }} {orelse}"

    def generic_visit(self, node: ast.AST) -> str:
        return "\\cdot"


class _ASTRenderer(ast.NodeVisitor):
    """Visite l'AST d'une expression et génère du LaTeX."""

    def __init__(self, local_vars: dict, precision: int):
        self.local_vars = local_vars
        self.precision = precision

    def render(self, node: ast.AST) -> str:
        return self.visit(node)

    # --- Nœuds terminaux ---

    def visit_Name(self, node: ast.Name) -> str:
        return _name_to_latex(node.id)

    def visit_Constant(self, node: ast.Constant) -> str:
        if isinstance(node.value, float):
            return f"{node.value}"
        return str(node.value)

    # --- Opérations binaires ---

    def visit_BinOp(self, node: ast.BinOp) -> str:
        left = self.visit(node.left)
        right = self.visit(node.right)
        op = node.op

        if isinstance(op, ast.Add):
            return f"{left} + {right}"
        if isinstance(op, ast.Sub):
            return f"{left} - {right}"
        if isinstance(op, ast.Mult):
            # Évite \cdot entre deux nombres
            return f"{left} \\cdot {right}"
        if isinstance(op, ast.Div):
            return f"\\frac{{{left}}}{{{right}}}"
        if isinstance(op, ast.Pow):
            # Parenthèses si la base est complexe
            base = left if _is_simple(node.left) else f"\\left({left}\\right)"
            exp_str = right if _is_simple(node.right) else f"{{{right}}}"
            return f"{base}^{{{exp_str}}}"
        if isinstance(op, ast.FloorDiv):
            return f"\\lfloor\\frac{{{left}}}{{{right}}}\\rfloor"
        if isinstance(op, ast.Mod):
            return f"{left} \\bmod {right}"
        return f"{left} \\circ {right}"

    def visit_UnaryOp(self, node: ast.UnaryOp) -> str:
        operand = self.visit(node.operand)
        if isinstance(node.op, ast.USub):
            return f"-{operand}"
        if isinstance(node.op, ast.UAdd):
            return operand
        return operand

    # --- Appels de fonctions ---

    def visit_Call(self, node: ast.Call) -> str:
        func_name = ""
        if isinstance(node.func, ast.Name):
            func_name = node.func.id
        elif isinstance(node.func, ast.Attribute):
            func_name = node.func.attr

        args_latex = [self.visit(a) for a in node.args]

        if func_name == "sqrt":
            return f"\\sqrt{{{args_latex[0]}}}" if args_latex else r"\sqrt{\cdot}"
        if func_name in ("sin", "cos", "tan"):
            return f"\\{func_name}\\left({args_latex[0]}\\right)" if args_latex else f"\\{func_name}()"
        if func_name in ("min", "max"):
            return f"\\{func_name}\\left({', '.join(args_latex)}\\right)"
        if func_name == "abs" or func_name == "fabs":
            return f"\\left|{args_latex[0]}\\right|" if args_latex else r"\left|\cdot\right|"
        if func_name == "log":
            return f"\\ln\\left({args_latex[0]}\\right)" if args_latex else r"\ln(\cdot)"
        if func_name == "log10":
            return f"\\log_{{10}}\\left({args_latex[0]}\\right)" if args_latex else r"\log_{10}(\cdot)"
        if func_name == "exp":
            return f"e^{{{args_latex[0]}}}" if args_latex else "e^{\\cdot}"
        if func_name == "radians":
            return f"\\frac{{\\pi \\cdot {args_latex[0]}}}{{180}}" if args_latex else ""
        if func_name == "floor":
            return f"\\lfloor {args_latex[0]} \\rfloor" if args_latex else ""
        if func_name == "ceil":
            return f"\\lceil {args_latex[0]} \\rceil" if args_latex else ""
        # Fonction générique
        name_latex = _name_to_latex(func_name)
        return f"{name_latex}\\left({', '.join(args_latex)}\\right)"

    # --- Comparaisons ---

    def visit_Compare(self, node: ast.Compare) -> str:
        parts = [self.visit(node.left)]
        ops_map = {
            ast.Lt: "<", ast.LtE: "\\leq", ast.Gt: ">", ast.GtE: "\\geq",
            ast.Eq: "=", ast.NotEq: "\\neq",
        }
        for op, comp in zip(node.ops, node.comparators):
            parts.append(ops_map.get(type(op), "?"))
            parts.append(self.visit(comp))
        return " ".join(parts)

    # --- Tuple/List/Dict ---

    def visit_Tuple(self, node: ast.Tuple) -> str:
        return "\\left(" + ", ".join(self.visit(e) for e in node.elts) + "\\right)"

    def visit_List(self, node: ast.List) -> str:
        return "\\left[" + ", ".join(self.visit(e) for e in node.elts) + "\\right]"

    def visit_Dict(self, node: ast.Dict) -> str:
        return ""

    def visit_Subscript(self, node: ast.Subscript) -> str:
        val = self.visit(node.value)
        sl = self.visit(node.slice)
        return f"{val}_{{{sl}}}"

    def visit_IfExp(self, node: ast.IfExp) -> str:
        body = self.visit(node.body)
        test = self.visit(node.test)
        orelse = self.visit(node.orelse)
        return f"{body} \\text{{ si }} {test} \\text{{ sinon }} {orelse}"

    def generic_visit(self, node: ast.AST) -> str:
        return "\\cdot"


def _is_simple(node: ast.AST) -> bool:
    """Retourne True si le nœud est suffisamment simple pour ne pas nécessiter de parenthèses."""
    return isinstance(node, (ast.Name, ast.Constant, ast.Attribute))


# ---------------------------------------------------------------------------
# Extraction des lignes d'assignation depuis l'AST de la fonction
# ---------------------------------------------------------------------------

def _extract_assignments(func_source: str) -> list[ast.Assign]:
    """Extrait les nœuds Assign et AugAssign du corps de la fonction."""
    tree = ast.parse(func_source)
    assignments = []
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            for stmt in node.body:
                if isinstance(stmt, ast.Assign):
                    assignments.append(stmt)
                elif isinstance(stmt, ast.AugAssign):
                    assignments.append(stmt)
                # Les Expr seuls (noms sans affectation comme `axe` ou `direction`)
                # sont ignorés — c'était des "labels" handcalcs, inutiles en LaTeX
    return assignments

# Rendu LaTeX d'une assignation
# ---------------------------------------------------------------------------

def _render_assignment(assign: ast.Assign, local_vars: dict, precision: int,
                        long: bool = False) -> str:
    """Génère une ligne LaTeX 'lhs &= rhs_sym = rhs_num = valeur_unité'.

    Si long=True, chaque étape (sym / num / résultat) est sur sa propre ligne alignée.
    """
    sym_renderer = _ASTRenderer(local_vars, precision)
    num_renderer = _ASTNumRenderer(local_vars, precision)

    # Cible (LHS)
    target = assign.targets[0] if isinstance(assign, ast.Assign) else assign.target
    lhs = sym_renderer.visit(target)

    # Expression symbolique
    rhs_sym = sym_renderer.visit(assign.value)

    # Expression numérique (substitution des variables connues)
    rhs_num = num_renderer.visit(assign.value)

    # Valeur finale calculée
    target_name = None
    if isinstance(target, ast.Name):
        target_name = target.id
    elif isinstance(target, ast.Tuple):
        target_name = None

    has_num_sub = rhs_num != rhs_sym

    if target_name and target_name in local_vars:
        val = local_vars[target_name]
        val_str = _num_to_latex(val, precision)
        unit_str = _get_unit_str(val)
        if unit_str:
            val_str = f"{val_str}{unit_str}"
        if long and has_num_sub:
            return (
                f"{lhs} &= {rhs_sym} \\\\\n"
                f"    &= {rhs_num} \\\\\n"
                f"    &= {val_str} \\\\\n"
            )
        elif long:
            return (
                f"{lhs} &= {rhs_sym} \\\\\n"
                f"    &= {val_str} \\\\\n"
            )
        elif has_num_sub:
            return f"{lhs} &= {rhs_sym} = {rhs_num} = {val_str} \\\\\n"
        else:
            return f"{lhs} &= {rhs_sym} = {val_str} \\\\\n"
    else:
        if long and has_num_sub:
            return (
                f"{lhs} &= {rhs_sym} \\\\\n"
                f"    &= {rhs_num} \\\\\n"
            )
        elif has_num_sub:
            return f"{lhs} &= {rhs_sym} = {rhs_num} \\\\\n"
        return f"{lhs} &= {rhs_sym} \\\\\n"


def _get_unit_str(value: Any) -> str:
    """Extrait l'unité d'une grandeur forallpeople en LaTeX (ex: \\mathrm{MPa})."""
    try:
        if _is_physical(value):
            s = str(value)
            parts = s.rsplit(" ", 1)
            if len(parts) == 2:
                unit = parts[1].strip()
                unit = unit.replace("**", "^").replace("*", " \\cdot ")
                return f"\\,\\mathrm{{{unit}}}"
    except Exception:
        pass
    return ""


# ---------------------------------------------------------------------------
# Décorateur principal
# ---------------------------------------------------------------------------

 
def handcalc(
    override: str = "short",
    precision: int = 2,
    jupyter_display: bool = False,
    left: str = "\\[",
    right: str = "\\]",
):
    """Remplace handcalcs.decorator.handcalc — compatible Python 3.13+.
 
    Usage identique à l'original :
 
        @handcalc(override="short", precision=2, jupyter_display=self.JUPYTER_DISPLAY,
                  left="\\\\[", right="\\\\]")
        def val():
            x = a + b
            return x
 
        latex_str, result = val()
 
    Args:
        override:        "short" (défaut) ou "long".
        precision:       Nombre de décimales pour les valeurs numériques.
        jupyter_display: Si True, affiche le LaTeX dans Jupyter via IPython.display.
        left / right:    Délimiteurs LaTeX du bloc (ex. "$$", "\\\\[").
 
    Returns:
        Décorateur qui transforme la fonction en un appelable retournant
        ``(latex_str, résultat_numérique)``.
    """
    def decorator(func):
        def wrapper(*args, **kwargs):
            # 1. Exécuter la fonction pour obtenir les valeurs numériques
            result = func(*args, **kwargs)
 
            # 2. Récupérer le code source de la fonction
            try:
                raw_source = inspect.getsource(func)
                func_source = textwrap.dedent(raw_source)
            except (OSError, TypeError):
                # Impossible de récupérer le source → retour sans LaTeX
                return ("", result)
 
            # 3. Construire le dict des variables locales disponibles
            #    = globals du scope + fermeture (closure) + arguments
            local_vars: dict = {}

            # Variables du scope module (globals) — permet la substitution numérique
            if func.__globals__:
                local_vars.update({
                    k: v for k, v in func.__globals__.items()
                    if not k.startswith("__") and not callable(v)
                })

            # Variables de fermeture (closures)
            if func.__code__.co_freevars and func.__closure__:
                for name, cell in zip(func.__code__.co_freevars, func.__closure__):
                    try:
                        local_vars[name] = cell.cell_contents
                    except ValueError:
                        pass

 
            # Arguments positionnels (cas `def comp(arg1, arg2, ...)`)
            arg_names = func.__code__.co_varnames[: func.__code__.co_argcount]
            for name, val in zip(arg_names, args):
                local_vars[name] = val
            local_vars.update(kwargs)
 
            # Résultats de la fonction — injecte les valeurs calculées
            # On récupère les noms des variables locales via co_varnames
            # En exécutant la fonction une seconde fois dans un namespace capturé
            local_capture: dict = dict(local_vars)
            try:
                _exec_and_capture(func, args, kwargs, local_capture)
            except Exception:
                pass
            local_vars.update(local_capture)
 
            # 4. Extraire les assignations de l'AST
            try:
                assignments = _extract_assignments(func_source)
            except SyntaxError:
                return ("", result)
 
            # 5. Construire le LaTeX
            is_long = override == "long"
            lines = []
            for assign in assignments:
                line = _render_assignment(assign, local_vars, precision, long=is_long)
                lines.append(line)
 
            if not lines:
                latex_body = ""
            elif override == "long" or len(lines) > 1:
                inner = "".join(lines)
                latex_body = f"\\begin{{aligned}}\n{inner}\\end{{aligned}}"
            else:
                # short : une seule ligne, format compact
                latex_body = "".join(lines).replace(" &=", " =").replace("\\\\\n", "")
 
            latex_str = f"{left}\n{latex_body}\n{right}"
 
            # 6. Affichage Jupyter si demandé
            if jupyter_display:
                try:
                    from IPython.display import display, Latex as IPyLatex
                    if display and IPyLatex:
                        display(IPyLatex(latex_str))
                except ImportError:
                    pass
 
            return (latex_str, result)
 
        return wrapper
    return decorator
 
 
# ---------------------------------------------------------------------------
# Capture des variables locales sans innerscope
# ---------------------------------------------------------------------------
 
def _exec_and_capture(func, args, kwargs, capture: dict) -> None:
    """Exécute func dans un namespace étendu et capture les variables locales.
 
    Utilise exec() sur le corps de la fonction avec les variables de fermeture
    et arguments déjà présents dans `capture`. Met à jour `capture` in-place.
    """
    try:
        raw = inspect.getsource(func)
        source = textwrap.dedent(raw)
    except (OSError, TypeError):
        return
 
    tree = ast.parse(source)
 
    # On cherche le premier FunctionDef
    func_def = None
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            func_def = node
            break
    if func_def is None:
        return
 
    # Reconstruire le corps sans le 'return' final comme module exécutable
    body_stmts = [s for s in func_def.body if not isinstance(s, ast.Return)]
    if not body_stmts:
        return
 
    # Créer un nouveau module AST avec ces instructions
    new_module = ast.Module(body=body_stmts, type_ignores=[])
    ast.fix_missing_locations(new_module)
 
    # Injecter les fonctions math disponibles + globals de la fonction (scope module)
    exec_globals = {
        "sqrt": sqrt, "pi": pi, "sin": sin, "cos": cos, "tan": tan,
        "radians": radians, "log": log, "log10": log10, "exp": exp,
        "floor": floor, "ceil": ceil, "fabs": fabs, "abs": abs,
        "min": min, "max": max, "round": round,
    }
    if func.__globals__:
        exec_globals.update(func.__globals__)
    exec_globals.update(capture)
 
    try:
        exec(compile(new_module, "<handcalc>", "exec"), exec_globals, capture)
    except Exception:
        pass