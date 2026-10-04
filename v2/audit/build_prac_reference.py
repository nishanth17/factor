"""Extract the existing arithmetic and PRAC transformations; repair chain selection.
This generates a separate audit prototype, never edits the original repository.
"""
import ast,contextlib,io,pathlib
with contextlib.redirect_stdout(io.StringIO()): import validate as v
src=v.ROOT.joinpath('ecm.py').read_text().expandtabs(8)+'\n'
tree=ast.parse(str(v.fixer.refactor_string(src,'ecm.py')))
class IntDiv(ast.NodeTransformer):
    def visit_BinOp(self,node):
        self.generic_visit(node)
        if isinstance(node.op,ast.Div):node.op=ast.FloorDiv()
        return node
    def visit_AugAssign(self,node):
        self.generic_visit(node)
        if isinstance(node.op,ast.Div):node.op=ast.FloorDiv()
        return node
tree=IntDiv().visit(tree)
nodes={x.name:x for x in tree.body if isinstance(x,ast.FunctionDef)}
parts=['"""Audit prototype: guarded PRAC chain selection using the original point updates.\nNot a complete elliptic-curve group API; see report for exceptional-point limits.\n"""','from math import gcd','ADD_COST = 6\nDUP_COST = 5']
for name in ['point_add','point_double','scalar_multiply']:
    node=nodes[name]
    if name=='scalar_multiply':
        guard=ast.parse('if k < 0: raise ValueError("negative scalar")\nif k == 0: return (1, 0)').body
        node.body[1:1]=guard
    parts.append(ast.unparse(ast.fix_missing_locations(node)))
node=nodes['lucas_cost'];node.name='_chain_cost';node.args.args[1].arg='r';node.body=node.body[3:]
parts.append(ast.unparse(ast.fix_missing_locations(node)))
parts.append('''RATIOS = (61803398874989485, 58017872829546410, 61791440652881790, 61807966846989580)
DENOM = 10**17

def select_chain(k):
    candidates = set()
    for numerator in RATIOS:
        r = (k*numerator + DENOM//2)//DENOM
        for candidate in (r-1, r, r+1):
            if 0 < 2*candidate-k and candidate < k and gcd(k,candidate) == 1:
                candidates.add(candidate)
    if not candidates:
        return None
    return min(candidates, key=lambda r: (_chain_cost(k,r),r))
''')
node=nodes['multiply_prac']
prefix=ast.parse('''if k < 0: raise ValueError("negative scalar")
if k == 0: return (1,0)
if k == 1: return px,pz
if k == 2: return point_double(px,pz,n,a24)
r = select_chain(k)
if r is None: return scalar_multiply(k,px,pz,n,a24)
''').body
node.body=prefix+node.body[:2]+node.body[6:]
# Prove the terminal coefficient rather than silently returning at d=e>1.
node.body.insert(-2,ast.parse('assert d == e == 1').body[0])
parts.append(ast.unparse(ast.fix_missing_locations(node)))
path=pathlib.Path(__file__).with_name('prac_reference.py');path.write_text('\n\n'.join(parts)+'\n')
print(path)
