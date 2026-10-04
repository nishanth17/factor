import json,random,time
import prac_reference as p
# Independent affine arithmetic over prime fields.
def add(P,Q,n,A):
    if P is None:return Q
    if Q is None:return P
    x,y=P;u,w=Q
    if x==u and (y+w)%n==0:return None
    slope=((3*x*x+2*A*x+1)*pow(2*y,-1,n) if P==Q else (w-y)*pow(u-x,-1,n))%n
    X=(slope*slope-A-x-u)%n
    return X,(slope*(x-X)-y)%n
checks=0; failures=[];degenerate=[];nondegenerate_failures=[]
for n,A in [(1009,6),(1013,11),(10007,6),(10009,17)]:
    points=[]
    for x in range(2,n):
        val=(x**3+A*x*x+x)%n
        ys=[y for y in range(1,n) if y*y%n==val]
        if ys:points.append((x,ys[0]))
        if len(points)==4:break
    a24=(A+2)*pow(4,-1,n)%n
    for P in points:
        R=None
        for k in range(0,1001):
            if k:R=add(R,P,n,A)
            X,Z=p.multiply_prac(k,P[0],1,n,a24)
            checks+=1
            if X==Z==0:
                degenerate.append([n,A,P,k]);continue
            good=(Z==0) if R is None else Z!=0 and (X-R[0]*Z)%n==0
            if not good:failures.append([n,A,P,k,[X,Z],R])
# k=9 chosen-r evidence and exact integer invariant for many scalar values.
print(json.dumps(dict(checks=checks,nondegenerate_mismatch_count=len(failures),first_mismatches=failures[:8],degenerate_count=len(degenerate),first_degenerate=degenerate[:8],k9_selected_r=p.select_chain(9),k9_cost=p._chain_cost(9,p.select_chain(9))),indent=2))
