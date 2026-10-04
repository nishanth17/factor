import contextlib,io,json,math,signal
with contextlib.redirect_stdout(io.StringIO()):
    import validate as v

def timeout(sig,frame): raise TimeoutError('one-second budget')
signal.signal(signal.SIGALRM,timeout)
out=[]
for label,func in [
 ('segmented endpoint 61 squared',lambda:dict(has_composite_3721=3721 in v.sieve.segmented_sieve(3700,3721),primes=v.sieve.segmented_sieve(3700,3721))),
 ('empty binary search',lambda:v.utils.binary_search(3,[])),
 ('duplicate binary search',lambda:dict(gt=v.utils.binary_search(2,[1,2,2,2,3]),ge=v.utils.binary_search(2,[1,2,2,2,3],True))),
 ('is_prime_fast(29)',lambda:v.utils.is_prime_fast(29))]:
    try:out.append(dict(case=label,result=func()))
    except Exception as e:out.append(dict(case=label,error=repr(e)))
# A valid curve and point for projective equality cross-checks.
n,A=1009,6
a24=(A+2)*pow(4,-1,n)%n
P=next((x,y) for x in range(2,n) for y in range(1,n) if (y*y-x*x*x-A*x*x-x)%n==0)
def add(P,Q):
    if P is None:return Q
    if Q is None:return P
    x,y=P;u,w=Q
    if x==u and (y+w)%n==0:return None
    s=((3*x*x+2*A*x+1)*pow(2*y,-1,n) if P==Q else (w-y)*pow(u-x,-1,n))%n
    X=(s*s-A-x-u)%n
    return X,(s*(x-X)-y)%n
R=None
mismatches=[]
for k in range(1,200):
    R=add(R,P)
    X,Z=v.ecm.scalar_multiply(k,P[0],1,n,a24)
    if (Z==0 if R is None else (X-R[0]*Z)%n==0) is False:
        mismatches.append(k)
out.append(dict(case='ladder vs independent affine multiplication k=1..199',point=P,modulus=n,mismatches=mismatches))
checks=[]
for k in range(2,41):
    try:
        signal.alarm(1)
        X,Z=v.ecm.multiply_prac(k,P[0],1,n,a24)
        signal.alarm(0)
        U,V=v.ecm.scalar_multiply(k,P[0],1,n,a24)
        checks.append(dict(k=k,equal=(X*V-U*Z)%n==0,degenerate=X==0 and Z==0))
    except Exception as e:
        signal.alarm(0)
        checks.append(dict(k=k,error=repr(e)))
out.append(dict(case='PRAC vs ladder, curve modulus 1009',checks=checks))
print(json.dumps(out,indent=2))
