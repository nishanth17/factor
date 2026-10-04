"""Read-only audit of original Python 2 source via a Python 3 compatibility loader.
No algorithm fixes: preserve Python 2 integer division, adapt syntax and gcd API.
Run with /usr/bin/python3 (3.9, includes lib2to3).
"""
import ast, contextlib, io, json, math, pathlib, random, sys, types
from lib2to3.refactor import RefactoringTool, get_fixers_from_package
# Historical compatibility checks intentionally target the preserved baseline.
ROOT = pathlib.Path(__file__).resolve().parents[2] / 'v1'
fixer = RefactoringTool(get_fixers_from_package('lib2to3.fixes'))
def py2_div(a, b):
    return a // b if isinstance(a, int) and isinstance(b, int) else a / b
class Divisions(ast.NodeTransformer):
    def visit_BinOp(self, node):
        self.generic_visit(node)
        if isinstance(node.op, ast.Div):
            return ast.copy_location(ast.Call(func=ast.Name(id='_py2_div', ctx=ast.Load()), args=[node.left, node.right], keywords=[]), node)
        return node
    def visit_AugAssign(self, node):
        self.generic_visit(node)
        if isinstance(node.op, ast.Div):
            assert isinstance(node.target, ast.Name)
            return ast.copy_location(ast.Assign(targets=[node.target], value=ast.Call(func=ast.Name(id='_py2_div',ctx=ast.Load()), args=[ast.Name(id=node.target.id,ctx=ast.Load()),node.value],keywords=[])), node)
        return node
mods = {}
for name in ['constants', 'utils', 'primeSieve', 'pollardRho', 'pollardPm1', 'ecm', 'factor']:
    source = (ROOT / (name + '.py')).read_text().expandtabs(8) + '\n'
    translated = str(fixer.refactor_string(source, name+'.py'))
    tree = ast.fix_missing_locations(Divisions().visit(ast.parse(translated)))
    mod = types.ModuleType(name)
    mod.__dict__['_py2_div'] = py2_div
    sys.modules[name] = mod
    exec(compile(tree, str(ROOT / (name+'.py')), 'exec'), mod.__dict__)
    if name == 'utils': mod.gcd = math.gcd
    mods[name] = mod
utils, sieve, rho, pm1, ecm, factor = [mods[k] for k in ['utils','primeSieve','pollardRho','pollardPm1','ecm','factor']]
def reference_primes(lo, hi):
    a = bytearray(b'\1') * (hi+1)
    a[:2] = b'\0\0'
    for p in range(2, math.isqrt(hi)+1):
        if a[p]: a[p*p:hi+1:p] = b'\0' * ((hi-p*p)//p+1)
    return [p for p in range(max(2,lo),hi+1) if a[p]]
results = []
def record(label, func):
    try:
        with contextlib.redirect_stdout(io.StringIO()): value = func()
        results.append({'case': label, 'result': value})
    except Exception as ex:
        results.append({'case': label, 'exception': type(ex).__name__, 'message': str(ex)})

def suyama():
    n,sigma = 101,6
    u,v=(sigma*sigma-5)%n,(4*sigma)%n
    a_original=((v-u)**3*(3*u+v)//(4*u**3*v)-2)%n
    a24_original=(a_original+2)//4
    x_original=(u**3//v**3)%n
    a24_correct=((v-u)**3*(3*u+v)*pow(16*u**3*v,-1,n))%n
    x_correct=(u**3*pow(v**3,-1,n))%n
    return dict(n=n,sigma=sigma,u=u,v=v,original_A=a_original,original_a24=a24_original,correct_a24=a24_correct,original_x=x_original,correct_affine_x=x_correct)
record('ECM curve construction arithmetic', suyama)
record('direct Atkin(100)',lambda:sieve.sieve_of_atkin(100))
def atkin_dispatch():
    n=3500001
    actual=sieve.prime_sieve(n)
    expected=reference_primes(2,n)
    mismatch=next((i for i,(a,b) in enumerate(zip(actual,expected)) if a!=b),None)
    return dict(n=n,len_actual=len(actual),len_expected=len(expected),first_difference_index=mismatch,actual=actual[mismatch],expected=expected[mismatch],last_actual=actual[-1],sorted=actual==sorted(actual),values_over_bound=sum(p>n for p in actual))
record('Atkin dispatch boundary',atkin_dispatch)
for lo,hi in [(2,9),(2,10),(2,49),(10,100),(90,100),(100,101),(0,10),(1,10),(2,1000),(2,70000)]:
    record('segmented_sieve(%d,%d)'%(lo,hi),lambda lo=lo,hi=hi:dict(actual=sieve.segmented_sieve(lo,hi),expected=reference_primes(lo,hi)) if hi<1000 else dict(equal=sieve.segmented_sieve(lo,hi)==reference_primes(lo,hi)))
record('sieve endpoint conventions',lambda:{n:sieve.prime_sieve(n) for n in [2,3,59,60,61,67]})
record('factorize(0)',lambda:factor.factorize(0))
record('factorize(-15)',lambda:factor.factorize(-15))
record('factorize_bf(10**400)',lambda:factor.factorize_bf(10**400))
record('is_prime_fast(10**36+1)',lambda:utils.is_prime_fast(10**36+1))
record('pm1(9)',lambda:pm1.factorize_pm1(9))
record('scalar multiply zero',lambda:ecm.scalar_multiply(0,11,16,29,7))
record('xgcd coefficient convention',lambda:dict(xgcd_7_60=utils.xgcd(7,60),inverse60_mod7=pow(60,-1,7)))

original_rho=rho.factorize_rho
original_ecm=ecm.factorize_ecm
calls=[]
rho.factorize_rho=lambda *a,**k:-1
ecm.factorize_ecm=lambda *a,**k:(calls.append(a[0]) or -1)
record('forced rho failure below threshold',lambda:dict(n=25013*25031,factorization=factor.factorize(25013*25031),ecm_calls=calls[:]))
record('forced ECM failure above threshold',lambda:dict(n=1000000000039*1000000000061,factorization=factor.factorize(1000000000039*1000000000061)))
rho.factorize_rho=original_rho
ecm.factorize_ecm=original_ecm

# p=607, p-1=2*3*101; q=101 is the only stage-2 prime needed,
# and is omitted by the implementation's checkpoints.
old_bounds=pm1.compute_bounds
pm1.compute_bounds=lambda n:(10,200)
record('pm1 stage2 missed relation n=607*1019, B1=10 B2=200',lambda:dict(result=pm1.factorize_pm1(607*1019),individual_gcds=[(q,math.gcd(pow(pow(2,2520,607*1019),q,607*1019)-1,607*1019)) for q in reference_primes(11,200) if math.gcd(pow(pow(2,2520,607*1019),q,607*1019)-1,607*1019)>1]))
pm1.compute_bounds=old_bounds

# Deterministic reproduction of premature rho failure.
old_rand=rho.random.randint
class FixedRandom:
    def __init__(self,y,m):self.values=iter((y,m));self.calls=0
    def __call__(self,*args):self.calls+=1;return next(self.values)
def rho_failure():
    for n in [35,77,143,323,10403]:
        for y in range(1,min(n,200)):
            fixed=FixedRandom(y,n-1)
            rho.random.randint=fixed
            r=rho.factorize_rho(n)
            if r==-1:return dict(n=n,y=y,m=n-1,result=r,random_calls=fixed.calls,available_offsets=len(rho.small_primes))
record('rho exits after first failed offset',rho_failure)
rho.random.randint=old_rand
record('float prime-power undercount',lambda:[dict(B1=b,p=p,float_exponent=int(math.log(b)/math.log(p)),exact_exponent=e) for p in [2,3,5,7,11,13,17,19] for e in range(2,15) for b in [p**e] if int(math.log(b)/math.log(p))!=e][:12])
print(json.dumps({'python':sys.version,'method':'syntax translation + Python2 arithmetic emulation; utils.gcd mapped to math.gcd; not a native Python2 benchmark','results':results},indent=2))
