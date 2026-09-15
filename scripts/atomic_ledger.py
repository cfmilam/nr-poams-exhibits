#!/usr/bin/env python3
"""Reproduce the isotope-count table and circular two-end energy/action ledger.

Only standard-library Python is required. Run from any directory. Data outputs
are placed beside this script. The isotope list is inherited exhibit data;
its proposed physical partition is not inferred by these arithmetic checks.
"""
import ast
import csv
import json
import re
from decimal import Decimal as D, localcontext
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'poams-periodic-table.html'
INPUTS = json.loads((HERE / 'coulomb_calibrated_inputs.json').read_text())
PI = D('3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628')
CHECKS = []


def check(label, condition):
    if not condition:
        raise AssertionError(label)
    CHECKS.append(label)


def close(label, actual, expected):
    scale = max(abs(actual), abs(expected), D('1e-100'))
    check(label, abs(actual - expected) / scale < D('1e-65'))


def data():
    raw = re.search(r'const RAW_ELEMENTS = (\[.*?\n        \]);', SOURCE.read_text(), re.S)
    if not raw:
        raise ValueError('Cannot locate inherited element data')
    records = ast.literal_eval(raw[1])
    check('118 distinct sequential element labels', [r[0] for r in records] == list(range(1, 119)))
    rows = []
    for Z, symbol, name, _kind, _mass, A, isotopes in records:
        available = list(dict.fromkeys([A, *(int(re.search(r'-(\d+)', s)[1]) for s in isotopes)]))
        check(f'{symbol}: valid isotope labels', all(a >= Z and isinstance(a, int) for a in available))
        rows.append(dict(Z=Z, symbol=symbol, name=name, A=A, depth=A-Z, pattern=Z,
                         ratio=str(D(A-Z)/Z), depth_fraction=str(D(A-Z)/A),
                         available_isotopes=available))
    check('carbon isotope example', [a-6 for a in rows[5]['available_isotopes']] == [6,7,8])
    check('hydrogen zero surplus is a count result', rows[0]['depth'] == 0)
    return rows


def circular(model, n):
    value = lambda key: D(INPUTS['inputs'][key]['value'])
    h, c, alpha = value('h'), value('c'), value('alpha')
    m1 = value('m_light')
    m2 = value('m_heavy') if model == 'pair' else m1
    M = m1+m2
    mu = m1*m2/M
    hbar = h/(2*PI)
    L = n*hbar
    K = alpha*hbar*c
    # Solve dE/dr=0 at fixed angular momentum, then calculate the endpoints
    # directly. Binding is computed from the sum, not stipulated as T.
    r = L*L/(mu*K)
    omega = L/(mu*r*r)
    r1, r2 = m2*r/M, m1*r/M
    T1, T2 = m1*(omega*r1)**2/2, m2*(omega*r2)**2/2
    L1, L2 = m1*r1*r1*omega, m2*r2*r2*omega
    U = -K/r
    B = -(T1+T2+U)
    f = omega/(2*PI)
    nu = B/h
    S1, S2 = 2*PI*L1, 2*PI*L2
    label = f'{model}, n={n}: '
    close(label+'stationarity', L*L/(mu*r**3), K/r**2)
    check(label+'positive radial curvature', 3*L*L/(mu*r**4)-2*K/r**3 > 0)
    close(label+'barycentre', m1*r1, m2*r2)
    close(label+'relative radius sum', r1+r2, r)
    close(label+'angular-momentum sum', L1+L2, L)
    close(label+'virial and one joint energy', 2*(T1+T2), -U)
    close(label+'binding expression', B, mu*c*c*alpha*alpha/(2*n*n))
    close(label+'rate distinction', nu/f, D(n)/2)
    close(label+'action sum', S1+S2, n*h)
    close(label+'energy share', T1/(T1+T2), m2/M)
    close(label+'angular-momentum share', L1/L, m2/M)
    close(label+'radius/reference ratio', r/(hbar/(mu*c)), n*n/alpha)
    close(label+'L omega is twice kinetic energy', L*omega, 2*(T1+T2))
    if model == 'equal':
        close(label+'equal endpoint split', T1, T2)
    return dict(model=model,n=n,mu=mu,K=K,radius=r,omega=omega,
                binding=B,bindingEV=B/value('e_SI'),bindingRate=nu,f=f,
                r1=r1,r2=r2,T1=T1,T2=T2,L=L,L1=L1,L2=L2,S1=S1,S2=S2,
                U=U,share1=m2/M,share2=m1/M)


def serialize(value):
    if isinstance(value, D):
        return format(value, '.30g')
    raise TypeError(type(value).__name__)


def main():
    with localcontext() as ctx:
        ctx.prec = 80
        rows = data()
        cases = [circular(model,n) for model in ('pair','equal') for n in range(1,7)]
        # Scaling checks compare separately solved cases, not a display label.
        for start in (0,6):
            base = cases[start]
            for q in cases[start+1:start+6]:
                n = q['n']
                close(f"{q['model']} n={n}: radius progression", q['radius']/base['radius'], D(n*n))
                close(f"{q['model']} n={n}: rate progression", base['f']/q['f'], D(n**3))
        result = dict(
            title='Atomic count and calibrated two-ended ledger', date='2026-09-15',
            inputs=INPUTS['inputs'],
            used_inputs=['alpha','c','h','e_SI','m_light','m_heavy'],
            source=INPUTS['source'],
            precision='80 working decimal digits, 30 printed; input uncertainty limits physical accuracy',
            assumptions=[
                'Counts are inherited neutral-isotope assignments, not measured regional angular momenta.',
                'Integer mass-number labels are not endpoint mass calibrations.',
                'E=L^2/(2 mu r^2)-K/r, K=alpha hbar c, L=n hbar; nonrelativistic circular scale model.',
                'The pair endpoints are not identified as a many-entry atom\'s core and envelope.',
                'Intrinsic spin, many-entry state coupling and region sizes are not calculated.',
                'The action here is angular action per circuit, not the full time-integrated Lagrangian action.'
            ],
            representatives=rows, pair_cases=cases,
            checks_passed=len(CHECKS), checks=CHECKS)
        (HERE/'atomic_ledger_results.json').write_text(json.dumps(result,default=serialize,indent=2)+'\n')
        with (HERE/'atomic_count_partitions.csv').open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=['Z','symbol','name','A','depth','pattern','ratio','depth_fraction'])
            writer.writeheader()
            writer.writerows({key:row[key] for key in writer.fieldnames} for row in rows)
        print(json.dumps(dict(checks_passed=len(CHECKS),representatives=len(rows),cases=len(cases),
                              fundamental=cases[0]),default=serialize,indent=2))


if __name__ == '__main__':
    main()
