#!/usr/bin/env python3
"""Exact checks for the conditional mathematics in inertia.html."""
from pathlib import Path
import json, math

def outer(v): return [[v[i]*v[j] for j in range(3)] for i in range(3)]
def add(a,b): return [[a[i][j]+b[i][j] for j in range(3)] for i in range(3)]
def scale(c,a): return [[c*a[i][j] for j in range(3)] for i in range(3)]
def sub(a,b): return add(a,scale(-1,b))
def norm(v): return math.sqrt(sum(x*x for x in v))
def unit(v):
    z=norm(v); return [x/z for x in v]
def eye(): return [[1.0 if i==j else 0.0 for j in range(3)] for i in range(3)]
def op(records):
    total=sum(a for _,a in records); assert total>0
    M=[[0.0]*3 for _ in range(3)]
    for v,a in records:
        n=unit(v); M=add(M,scale(a/total,sub(eye(),outer(n))))
    return M
def err(a,b): return math.sqrt(sum((a[i][j]-b[i][j])**2 for i in range(3) for j in range(3)))
def trace(a): return sum(a[i][i] for i in range(3))

axes=[([1,0,0],1),([-1,0,0],1),([0,1,0],1),([0,-1,0],1),([0,0,1],1),([0,0,-1],1)]
iso=scale(2/3,eye()); Miso=op(axes)
base=[([1,0,0],17),([0,1,0],31),([0,0,1],43),([1,2,3],59)]
split=[]
for v,a in base:
    q,r=divmod(a,7); split.extend((v,q+(k<r)) for k in range(7))
Mb,Ms=op(base),op(split)
parity=err(op([([1,2,3],1)]),op([([-1,-2,-3],1)]))

# Affine covariant response family K_k/m = I+k[(3/2)M-I].  For every
# admissible M, lambda(M) is in [0,1], so positivity for all accounts is
# equivalent to -2 <= k <= 1.  Isotropy cannot distinguish any member.
def response_eigs(lambdas,k): return [1+k*(1.5*x-1) for x in lambdas]
kappas=(-2,-1,0,0.5,1)
kappa_countermodels={str(k):response_eigs([0,1,1],k) for k in kappas}
isotropic_kappa={str(k):response_eigs([2/3,2/3,2/3],k) for k in kappas}
payload={
 "status":"conditional mathematics; no physical inertia measurement",
 "isotropic_operator":Miso,"isotropic_target":iso,"isotropic_error":err(Miso,iso),
 "trace":trace(Mb),"receipt_partition_error":err(Mb,Ms),"parity_error":parity,
 "kappa_positive_interval":[-2,1],
 "single_axis_kappa_countermodels":kappa_countermodels,
 "isotropic_kappa_indistinguishability":isotropic_kappa,
 "adopted_selector":"NOIR, ratified 2026-10-04, forces kappa=1 because the zero-comparison direction must have zero response",
 "estimator_guard":"detector receipt measure is not automatically the constitutive custody measure",
 "synthetic_warning":"locked source/hardware rotation is structurally aliased; see exhibit design section"
}
Path(__file__).with_name('inertia_operator_results.json').write_text(json.dumps(payload,indent=2)+'\n')
assert payload['isotropic_error']<1e-15
assert abs(payload['trace']-2)<1e-14
assert payload['receipt_partition_error']<1e-15
assert payload['parity_error']<1e-15
assert all(v==[1.0,1.0,1.0] for v in isotropic_kappa.values())
assert all(min(v)>=0 for v in kappa_countermodels.values())
assert kappa_countermodels['1'][0]==0
print(json.dumps(payload,indent=2));print('PASS: isotropy, trace, receipt refinement, parity, kappa countermodels, conditional NOIR selector')
