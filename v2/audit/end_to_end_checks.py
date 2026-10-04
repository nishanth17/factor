import contextlib,io,json,math,random,time
with contextlib.redirect_stdout(io.StringIO()): import validate as v
out={}
n=25013*25031
for seed in range(500):
    random.seed(seed)
    ans=v.factor.factorize(n)
    product=math.prod(p**e for p,e in ans) if ans!=-1 else None
    if product!=n:
        out['unmodified_algorithm_failure']={'n':n,'known_factorization':[25013,25031],'seed':seed,'returned':ans,'reconstructed':product}
        break
out['small_positive_integer_sweep']={'range':[1,20000],'mismatches':[]}
for n in range(1,20001):
    random.seed(n)
    ans=v.factor.factorize(n)
    if ans==-1 or math.prod(p**e for p,e in ans)!=n:
        out['small_positive_integer_sweep']['mismatches'].append([n,ans])

# Time equivalent point-arithmetic kernels, not full ECM, on CPython 3.9.
def early_double(px,pz,n,a24):
    u2=(px+pz)**2%n;v2=(px-pz)**2%n;t=(u2-v2)%n
    return u2*v2%n,t*(v2+a24*t)%n
def early_add(px,pz,qx,qz,rx,rz,n):
    u=(px-pz)*(qx+qz)%n;v=(px+pz)*(qx-qz)%n
    return rz*((u+v)**2%n)%n,rx*((u-v)**2%n)%n
random.seed(947)
bench=[]
for bits in [64,128,256,512,1024]:
    n=(1<<bits)-159
    vals=[random.randrange(n) for _ in range(7)]
    a=vals[:2]+[n,vals[6]]; b=vals[:6]+[n]
    assert v.ecm.point_double(*a)==early_double(*a)
    assert v.ecm.point_add(*b)==early_add(*b)
    row={'bits':bits,'iterations':10000}
    for label,fn,args in [('original_double',v.ecm.point_double,a),('early_double',early_double,a),('original_add',v.ecm.point_add,b),('early_add',early_add,b)]:
        times=[]
        for _ in range(3):
            t=time.perf_counter()
            for j in range(10000): fn(*args)
            times.append(time.perf_counter()-t)
        row[label+'_ms']=round(1000*sorted(times)[1],3)
    bench.append(row)
out['kernel_microbenchmarks']=bench
out['memory_estimates']=[dict(B2=b,estimated_prime_count=int(b/math.log(b)),estimated_list_GiB=round(36*b/math.log(b)/2**30,2)) for b in [12746592,128992510,1045563762,5706890290,20000000000]]
print(json.dumps(out,indent=2))
