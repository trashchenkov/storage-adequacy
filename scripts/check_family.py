#!/usr/bin/env python3
"""Exact physical checks of the explanatory example and general family.

These are rational-arithmetic feasibility checks, not additional MILP solves or
a replacement for the universal proofs in the manuscript.
"""

if not __debug__:
    raise SystemExit("Run without -O or PYTHONOPTIMIZE: assertions validate the results.")

import argparse
import json
from fractions import Fraction as F
from pathlib import Path

checks=[]

def physical(x, eta, energy, pc, pd, c, d, r):
    states=[F(0) for _ in eta]
    shortage=[]
    for t,xt in enumerate(x):
        s=xt+r
        assert sum(ci[t] for ci in c)<=max(s,0)
        assert sum(di[t] for di in d)<=max(-s,0)
        for i in range(len(eta)):
            assert 0<=c[i][t]<=pc[i] and 0<=d[i][t]<=pd[i]
            states[i]+=eta[i]*c[i][t]-d[i][t]
            assert 0<=states[i]<=energy[i]
        shortage.append(max(-s,0)-sum(di[t] for di in d))
    assert states==[0 for _ in eta]
    return sum(shortage),sum(z>0 for z in shortage)

def tutorial():
    args=([F(2),F(-1),F(-1)],[F(1),F(1)],[F(1),F(1)],[F(1),F(1)],[F(1),F(1,2)])
    old=physical(*args,[[F(1),F(0),F(0)],[F(1,2),F(0),F(0)]],[[F(0),F(1),F(0)],[F(0),F(0),F(1,2)]],F(0))
    new=physical(*args,[[F(1),F(0),F(0)],[F(1),F(0),F(0)]],[[F(0),F(1,2),F(1,2)],[F(0),F(1,2),F(1,2)]],F(0))
    assert old==(F(1,2),1) and new==(0,0)
    checks.append(dict(kind='tutorial',ENS_before=str(old[0]),ENS_after=str(new[0]),passed=True))

def family(a,b,p,offset,N,slack_fraction):
    k=b/a; q=k-1; K=N*k+1+b; L=k*N*p; Q=L*(N+1)/b
    emin=(N+1)*p+offset
    slack=slack_fraction*N*p*q
    r=max(F(0),(N*p*q-slack)/K)
    x=[Q]+[-(L+p)]*N+[-(L+p+offset)]
    E=[N*p,L*(N+1)]; pc=[N*p/a,Q]; pd=[p,L]
    cA=N*(p-r)/a; cB=Q+r-cA; y=(K-1)*r
    c=[[cA]+[F(0)]*(N+1),[cB]+[F(0)]*(N+1)]
    d=[[F(0)]+[p-r]*N+[F(0)],[F(0)]+[L]*N+[y]]
    ens,hrs=physical(x,[a,b],E,pc,pd,c,d,r)
    assert hrs==1 and ens==emin+N*p*q-K*r and ens<=emin+slack
    # Check the separate schedules against the same full physical model.
    energy_only=physical(x,[a,b],E,pc,pd,[[F(0)]*(N+2),[Q]+[F(0)]*(N+1)],[[F(0)]*(N+2),[F(0)]+[L]*(N+1)],F(0))
    hour_only=physical(x,[a,b],E,pc,pd,[[N*p/a]+[F(0)]*(N+1),[Q-N*p/a]+[F(0)]*(N+1)],[[F(0)]+[p]*N+[F(0)],[F(0)]+[L]*N+[F(0)]],F(0))
    assert energy_only[0]==emin and hour_only[1]==1
    delta=p*(1-a/b); bound=delta*N/(N+a)
    assert 0<=r<=bound<delta
    if slack==0:
        assert r/delta==N*k/K and r/bound==(N*k+b)/K
    checks.append(dict(kind='general_family',a=str(a),b=str(b),p=str(p),offset=str(offset),N=N,slack_fraction=str(slack_fraction),r=str(r),B=str(bound),ENS=str(ens),LOLH=hrs,passed=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    tutorial()
    for a,b in [(F(4,5),F(19,20)),(F(1,2),F(9,10)),(F(9,10),F(1)),(F(1,100),F(1)),(F(99,100),F(1))]:
        for N in [1,4,16]:
            for p in [F(1),F(3,2)]:
                for offset in [F(1,100),F(1),F(3)]:
                    for slack_fraction in [F(0),F(1,2),F(1),F(2)]:
                        family(a,b,p,offset,N,slack_fraction)
    result=dict(all_passed=True,arithmetic='fractions.Fraction',family_cases=len(checks)-1,tutorial_cases=1,checks=checks)
    path=args.output
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='checks'}))


if __name__ == "__main__":
    main()
