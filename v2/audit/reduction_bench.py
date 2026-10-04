"""Correctness and steady-state Python modular multiply comparisons.
Montgomery encode/decode and setup are verified but excluded from timings.
"""
import json,random,statistics,sys,time
rng=random.Random(82071)
rows=[]
for bits in (64,128,192,256,512,1024,2048):
    n=rng.getrandbits(bits)|(1<<(bits-1))|1
    k=n.bit_length();R=1<<k;mask=R-1;np=(-pow(n,-1,R))&mask;mu=(1<<(2*k))//n
    def native(a,b):return a*b%n
    def barrett(a,b):
        t=a*b;q=(t*mu)>>(2*k);r=t-q*n
        return r-n if r>=n else r
    def redc(t):
        m=((t&mask)*np)&mask;u=(t+m*n)>>k
        return u-n if u>=n else u
    def mont(a,b):return redc(a*b)
    pairs=[(rng.randrange(n),rng.randrange(n)) for _ in range(1000)]
    pairs += [(a,b) for a in (0,1,n-1) for b in (0,1,n-1)]
    for a,b in pairs:
        assert barrett(a,b)==native(a,b)
        assert redc(mont(a*R%n,b*R%n))==native(a,b)
    a,b=pairs[0];am,bm=a*R%n,b*R%n
    row=dict(bits=bits,checked_products=len(pairs),iterations=20000)
    for label,func,args in [('native',native,(a,b)),('barrett',barrett,(a,b)),('montgomery',mont,(am,bm))]:
        runs=[]
        for _ in range(5):
            t=time.perf_counter()
            for j in range(20000):out=func(*args)
            runs.append(1000*(time.perf_counter()-t))
        row[label+'_ms']=round(statistics.median(runs),3)
    row['barrett_over_native']=round(row['barrett_ms']/row['native_ms'],2)
    row['montgomery_over_native']=round(row['montgomery_ms']/row['native_ms'],2)
    rows.append(row)
print(json.dumps(dict(python=sys.version,setup_and_conversion_timed=False,results=rows),indent=2))
