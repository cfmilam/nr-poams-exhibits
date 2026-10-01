#!/usr/bin/env python3
"""Checks for the conditional POAMS/CCK scalar-pair Bell theorem."""
from itertools import product
from math import pi, sin, cos, sqrt
import numpy as np

tol=2e-12; checks=0
I=np.eye(2,dtype=complex); X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]],complex); Z=np.array([[1,0],[0,-1]],complex)
sig=(X,Y,Z); up=np.array([1,0],complex); dn=np.array([0,1],complex)
psi=(np.kron(up,dn)-np.kron(dn,up))/sqrt(2); rho=np.outer(psi,psi.conj())
def close(a,b,eps=tol):
 global checks
 assert abs(a-b)<eps,(a,b); checks+=1
close(np.vdot(psi,psi),1)
swap=np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]],complex)
assert np.linalg.norm(swap@psi+psi)<tol; checks+=1
J=[]
for s in sig:
 j=(np.kron(s,I)+np.kron(I,s))/2; J.append(j)
 assert np.linalg.norm(j@psi)<tol; checks+=1
J2=sum(j@j for j in J); evals=np.linalg.eigvalsh(J2)
assert sum(abs(x)<tol for x in evals)==1; checks+=1
close(evals[-1],2)
def vec(t): return np.array([sin(t),0,cos(t)])
def obs(v): return sum(v[i]*sig[i] for i in range(3))
def proj(v,s): return (I+s*obs(v))/2
for alpha,beta in product((0,.23,.71,1.4),(0,.4,1.1,2.2)):
 a,b=vec(alpha),vec(beta); probs={}
 for s,t in product((-1,1),repeat=2):
  p=float(np.real(np.trace(rho@np.kron(proj(a,s),proj(b,t)))))
  close(p,(1-s*t*np.dot(a,b))/4); probs[s,t]=p
 close(sum(probs.values()),1)
 for s in (-1,1): close(sum(probs[s,t] for t in (-1,1)),.5)
 for t in (-1,1): close(sum(probs[s,t] for s in (-1,1)),.5)
 close(sum(s*t*probs[s,t] for s,t in probs),-np.dot(a,b))
A0=vec(0); A1=vec(pi/2); B0=vec(pi/4); B1=vec(-pi/4)
E=lambda a,b,r=rho: float(np.real(np.trace(r@np.kron(obs(a),obs(b)))))
close(abs(E(A0,B0)+E(A0,B1)+E(A1,B0)-E(A1,B1)),2*sqrt(2),1e-10)
for q in (0,.2,.7,1):
 rw=q*rho+(1-q)*np.eye(4)/4
 for s in sig:
  close(np.trace(rw@np.kron(s,I)),0); close(np.trace(rw@np.kron(I,s)),0)
 close(-q,float(np.real(np.trace(rw@np.kron(Z,Z)))))
rho_cl=(np.outer(np.kron(up,dn),np.kron(up,dn).conj())+np.outer(np.kron(dn,up),np.kron(dn,up).conj()))/2
close(np.trace(rho_cl@np.kron(Z,Z)),-1); close(np.trace(rho_cl@np.kron(X,X)),0)

# Operational scalar receipt: Pi_s=(I-XX-YY-ZZ)/4 and J^2=2(I-Pi_s).
I4=np.eye(4,dtype=complex)
Pi_s=(I4-np.kron(X,X)-np.kron(Y,Y)-np.kron(Z,Z))/4
close(np.linalg.norm(Pi_s@Pi_s-Pi_s),0)
close(np.linalg.norm(Pi_s-rho),0)
J2_receipt=1.5*I4+0.5*(np.kron(X,X)+np.kron(Y,Y)+np.kron(Z,Z))
close(np.linalg.norm(J2_receipt-2*(I4-Pi_s)),0)
for q in (0,.2,.7,1):
 rw=q*rho+(1-q)*I4/4
 csum=sum(np.trace(rw@np.kron(s,s)) for s in sig)
 fidelity=float(np.real((1-csum)/4))
 close(fidelity,float(np.real(np.trace(rw@Pi_s))))
 close(fidelity,(1+3*q)/4)
 close(fidelity,1-float(np.real(np.trace(rw@J2_receipt)))/2)
print(f"PASS: {checks} scalar-pair, joint-law, CHSH, and source-countermodel checks")
