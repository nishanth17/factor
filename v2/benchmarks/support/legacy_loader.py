"""Load preserved Python 2 code on PyPy Python 3.11 for comparison only.

This follows the historical audit adapter: syntax translation, integer
division emulation, and fractions.gcd mapped to math.gcd. It does not turn
the old algorithms into a native Python 2 benchmark or repair their bugs.
Run in a dedicated process: original imports occupy legacy module names.
"""

import ast
import math
import sys
import types

from .paths import (
    REPOSITORY_ROOT,
    source_path,
)


def _python_two_division(a, b):
    return a // b if isinstance(a, int) and isinstance(b, int) else a / b


class _IntegerDivisions(ast.NodeTransformer):
    def visit_BinOp(self, node):
        self.generic_visit(node)
        if isinstance(node.op, ast.Div):
            call = ast.Call(
                ast.Name("_python_two_division", ast.Load()),
                [node.left, node.right],
                [],
            )
            return ast.copy_location(call, node)
        return node

    def visit_AugAssign(self, node):
        self.generic_visit(node)
        if isinstance(node.op, ast.Div):
            if not isinstance(node.target, ast.Name):
                raise ValueError("unsupported legacy augmented division")
            previous = ast.Name(node.target.id, ast.Load())
            call = ast.Call(
                ast.Name("_python_two_division", ast.Load()),
                [previous, node.value],
                [],
            )
            return ast.copy_location(ast.Assign([node.target], call), node)
        return node


def load_legacy(root=None):
    """Translate the frozen baseline without modifying any source file."""
    from lib2to3.refactor import RefactoringTool, get_fixers_from_package

    if root is None:
        root = REPOSITORY_ROOT / "v1"
    fixer = RefactoringTool(get_fixers_from_package("lib2to3.fixes"))
    modules = {}

    for name in (
        "constants",
        "utils",
        "primeSieve",
        "pollardRho",
        "pollardPm1",
        "ecm",
        "factor",
    ):
        path = root / f"{name}.py"
        source = source_path(path).read_text().expandtabs(8) + "\n"
        translated = str(fixer.refactor_string(source, str(path)))
        tree = ast.fix_missing_locations(
            _IntegerDivisions().visit(ast.parse(translated))
        )
        module = types.ModuleType(name)
        module._python_two_division = _python_two_division
        sys.modules[name] = module
        exec(compile(tree, str(path), "exec"), module.__dict__)
        if name == "utils":
            module.gcd = math.gcd
        modules[name] = module

    return modules
