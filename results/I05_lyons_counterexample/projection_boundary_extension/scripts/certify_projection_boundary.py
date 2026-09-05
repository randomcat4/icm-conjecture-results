#!/usr/bin/env python3
"""Uniform projection-boundary certificate, exact rational arithmetic only.

Run alongside certify.py. This author self-check does not claim independent
review. No floating-point determinant, eigensolver, epsilon scan or old JSON
certificate is used. Polynomial coefficient identities certify an entire
interval 0 < epsilon <= 2**(-N), not just one very small epsilon.
"""
from __future__ import annotations
import argparse
import json
from fractions import Fraction as F
from itertools import combinations
from math import comb
from pathlib import Path
import certify as c

Matrix = list[list[F]]

def transpose(a: Matrix) -> Matrix:
    return [list(row) for row in zip(*a)]

def madd(a: Matrix, b: Matrix, sign: int = 1) -> Matrix:
    return [[x+sign*y for x,y in zip(r,s)] for r,s in zip(a,b)]

def mscale(a: Matrix, v: F) -> Matrix:
    return [[v*x for x in row] for row in a]

def zeros(n: int, m: int) -> Matrix:
    return [[F(0) for _ in range(m)] for _ in range(n)]

def sub(a: Matrix, ids: list[int]) -> Matrix:
    return [[a[i][j] for j in ids] for i in ids]

def row_bound(a: Matrix) -> F:
    """|x^T A x| <= max absolute row sum * ||C||_F^2 for C=sum x_i E_i."""
    return max((sum(map(abs,row),F(0)) for row in a), default=F(0))

def poly_weight(k: int, j: int, derivative: int, length: int) -> list[F]:
    # epsilon^(k-j) (1-2 epsilon)^(j-derivative)
    out = [F(0)]*length
    for d in range(j-derivative+1):
        out[k-j+d] = F(comb(j-derivative,d)*(-2)**d)
    return out

def show(x: F) -> str:
    return str(x)

def pack_interval(x: c.Interval, places: int = 30) -> list[str]:
    scale=10**places
    return [str(F(c.floor(x[0]*scale),scale)),str(F(c.ceil(x[1]*scale),scale))]

def check_nongeneric_stratum() -> dict:
    """One exact 4D test: connected but not full spark; non-isotropic inward N."""
    V=[[F(1),F(0)],[F(1),F(0)],[F(1),F(1)],[F(1),F(-1)]]
    gi,_=c.inverse_det(c.matmul(transpose(V),V))
    P=c.matmul(c.matmul(V,gi),transpose(V));Q=madd(c.eye(4),P,-1)
    C0=zeros(4,4);C0[0][0]=F(1)
    E0=madd(c.matmul(c.matmul(P,C0),Q),c.matmul(c.matmul(Q,C0),P))
    N=madd(madd(c.eye(4),mscale(P,F(-2))),mscale(E0,F(1,3)))
    c.need(c.matmul(c.matmul(P,N),P)==mscale(P,F(-1)),'inward P block')
    c.need(c.matmul(c.matmul(Q,N),Q)==Q,'inward Q block')
    c.need(any(x for row in E0 for x in row),'test lacks transverse inward part')
    pairs=[(i,j) for i in range(4) for j in range(i,4)]
    dim=len(pairs);basis=[];db=[];eb=[]
    for i,j in pairs:
        C=zeros(4,4);C[i][j]=C[j][i]=F(1);basis.append(C)
        D=madd(c.matmul(c.matmul(P,C),P),c.matmul(c.matmul(Q,C),Q))
        E=madd(C,D,-1)
        db.append([D[u][v] for u,v in pairs]);eb.append([E[u][v] for u,v in pairs])
    J=zeros(dim,dim);G=zeros(dim,dim);T=zeros(dim,dim);orders=[]
    zero_bases=[]
    for mask in range(16):
        M=c.event_matrix(P,mask);sign=(-1)**(4-mask.bit_count())
        poly=[sign*x for x in c.determinant_polynomial(M,N)]
        w=next(i for i,x in enumerate(poly) if x);a=poly[w]
        ids=[i for i in range(4) if mask>>i&1]
        rank=0
        if ids:rank=1
        if any(P[i][i]*P[j][j]-P[i][j]**2>0 for i,j in combinations(ids,2)):rank=2
        c.need(w==2+len(ids)-2*rank and a>0,'nongeneric event order')
        g,H=c.leibniz_jet2(M,pairs);g=[sign*x for x in g];H=mscale(H,F(sign))
        J=madd(J,mscale(H,F(w)))
        if w==1:
            gd=[sum(x*y for x,y in zip(g,v)) for v in db]
            for i in range(dim):
                for j in range(dim):G[i][j]+=gd[i]*gd[j]/a
            c.need(all(sum(x*y for x,y in zip(g,v))==0 for v in eb),'transverse first defect derivative')
        if w==0:
            ge=[sum(x*y for x,y in zip(g,v)) for v in eb]
            for i in range(dim):
                for j in range(dim):T[i][j]+=ge[i]*ge[j]/a
        if len(ids)==2 and w:zero_bases.append(mask)
        orders.append(w)
    def bil(A: Matrix,x: list[F],y: list[F]) -> F:
        return sum((x[i]*A[i][j]*y[j] for i in range(dim) for j in range(dim)),F(0))
    for i in range(dim):
        for j in range(dim):
            c.need(bil(J,db[i],eb[j])==0,'mixed D,E logarithmic coefficient')
            c.need(bil(J,eb[i],eb[j])==-T[i][j],'transverse log coefficient is not -T')
    # Exact LDL of G+T on all ten real directions. G is applied to D, T to E.
    A=madd(G,T);L=c.eye(dim);pivots=[]
    for i in range(dim):
        pivot=A[i][i]-sum((L[i][k]**2*pivots[k] for k in range(i)),F(0))
        c.need(pivot>0,'G+T failed exact positive definiteness')
        pivots.append(pivot)
        for j in range(i+1,dim):
            L[j][i]=(A[j][i]-sum((L[j][k]*L[i][k]*pivots[k] for k in range(i)),F(0)))/pivot
    c.need(bool(zero_bases),'test unexpectedly full spark')
    return {'V_integer':[[int(x) for x in row] for row in V],
            'N_definition':'I-2P+(P C0 Q+Q C0 P)/3; C0=diag(1,0,0,0)',
            'zero_rank_two_base_masks':zero_bases,'event_vanishing_orders':orders,
            'all_mixed_and_transverse_logarithmic_entries':True,
            'G_plus_T_exact_positive_LDL_pivots':[show(x) for x in pivots]}


def main(output: Path) -> None:
    n,r=5,2
    V=[[F(row[0]),F(row[1])] for row in c.K_INT]
    gram=c.matmul(transpose(V),V)
    gi, gd=c.inverse_det(gram)
    P=c.matmul(c.matmul(V,gi),transpose(V)); Q=madd(c.eye(n),P,-1)
    B0=[[F(v,50) for v in row] for row in c.B_INT]
    B=madd(c.matmul(c.matmul(P,B0),P),c.matmul(c.matmul(Q,B0),Q))
    c.need(P==transpose(P) and c.matmul(P,P)==P,'not an orthogonal projection')
    c.need(c.trace(P)==r and c.trace(Q)==n-r,'incorrect ranks')
    c.need(B==mscale(transpose(B),F(-1)),'direction not skew')
    c.need(c.matmul(P,B)==c.matmul(B,P),'direction does not commute with P')
    pairs=[(i,j) for i in range(n) for j in range(i,n)]
    lookup={ij:a for a,ij in enumerate(pairs)};dim=len(pairs)
    basis=[]
    for i,j in pairs:
        e=zeros(n,n);e[i][j]=1;e[j][i]=1;basis.append(e)
    fullmask=2**n-1

    # Boundary principal-minor values and derivatives. There are only
    # sum_{k=0}^5 binom(5,k)*k! = 326 permutation terms in the value pass.
    boundary=[]
    for mask in range(2**n):
        ids=[i for i in range(n) if mask>>i&1];k=len(ids)
        pu,bu=sub(P,ids),sub(B,ids)
        val=F(1) if k==0 else c.determinant_polynomial(pu,zeros(k,k))[0]
        g=[F(0)]*dim;H=zeros(dim,dim)
        if k:
            local_pairs=[(i,j) for i in range(k) for j in range(i,k)]
            lg,lh=c.leibniz_jet2(pu,local_pairs)
            mapping=[lookup[(ids[i],ids[j])] for i,j in local_pairs]
            for a,ga in enumerate(mapping):
                g[ga]=lg[a]
                for b,gb in enumerate(mapping):H[ga][gb]=lh[a][b]
        dh=F(0)
        if k>=2:dh=-2*c.determinant_polynomial(pu,bu)[2]
        boundary.append((val,g,H,dh))

    # Inclusion polynomials; then exact Mobius inversion to complete events.
    inc=[]
    for mask in range(2**n):
        k=mask.bit_count()
        ps=[F(0)]*(n+1);gs=zeros(n,dim);hs=[zeros(dim,dim) for _ in range(n-1)]
        im=[F(0)]*(n-1)
        u=mask
        while True:
            j=u.bit_count();val,g,H,dh=boundary[u]
            for d,w in enumerate(poly_weight(k,j,0,n+1)):
                ps[d]+=w*val
            if j>=1:
                for d,w in enumerate(poly_weight(k,j,1,n)):
                    if w:
                        for a in range(dim):gs[d][a]+=w*g[a]
            if j>=2:
                for d,w in enumerate(poly_weight(k,j,2,n-1)):
                    if w:
                        for a in range(dim):
                            for b in range(dim):hs[d][a][b]+=w*H[a][b]
                        im[d]+=w*dh
            if u==0:break
            u=(u-1)&mask
        inc.append((ps,gs,hs,im))
    events=[]
    for mask in range(2**n):
        ps=[F(0)]*(n+1);gs=zeros(n,dim);hs=[zeros(dim,dim) for _ in range(n-1)]
        im=[F(0)]*(n-1)
        rest=fullmask^mask;u=rest
        while True:
            sign=(-1)**u.bit_count();data=inc[mask|u]
            for d in range(n+1):ps[d]+=sign*data[0][d]
            for d in range(n):
                for a in range(dim):gs[d][a]+=sign*data[1][d][a]
            for d in range(n-1):
                for a in range(dim):
                    for b in range(dim):hs[d][a][b]+=sign*data[2][d][a][b]
                im[d]+=sign*data[3][d]
            if u==0:break
            u=(u-1)&rest
        w=abs(mask.bit_count()-r)
        c.need(all(x==0 for x in ps[:w]) and ps[w]>0,'full-spark leading term failed')
        c.need(all(x==0 for row in gs[:max(0,w-1)] for x in row),'first-jet vanishing order')
        c.need(all(x==0 for m in hs[:max(0,w-2)] for row in m for x in row),'second-jet order')
        a=ps[w];f=ps[w:]
        R=sum(map(abs,f[1:]),F(0))
        expected=boundary[mask][0] if mask.bit_count()<=r else None
        if mask.bit_count()>r:
            ids=[i for i in range(n) if not mask>>i&1]
            expected=F(1) if not ids else c.determinant_polynomial(sub(Q,ids),zeros(len(ids),len(ids)))[0]
        c.need(a==expected,'leading coefficient disagrees with projection minor')
        events.append(dict(mask=mask,w=w,a=a,f=f,R=R,p=ps,g=gs,H=hs,im=im))

    for d in range(n+1):c.need(sum(e['p'][d] for e in events)==(1 if d==0 else 0),'normalization')
    for d in range(n):
        for a in range(dim):c.need(sum(e['g'][d][a] for e in events)==0,'gradient normalization')
    for d in range(n-1):
        c.need(sum(e['im'][d] for e in events)==0,'imaginary second normalization')
        for a in range(dim):
            for b in range(dim):c.need(sum(e['H'][d][a][b] for e in events)==0,'Hessian normalization')

    # Second algorithm: evaluate every event polynomial/jet at epsilon=1/4
    # directly with Gaussian elimination and the inverse trace formulas.
    check_eps=F(1,4)
    Kcheck=madd(mscale(c.eye(n),check_eps),mscale(P,1-2*check_eps))
    def evaluate(poly: list[F]) -> F:
        value=F(0)
        for coefficient in reversed(poly):value=value*check_eps+coefficient
        return value
    for e in events:
        inv,det=c.inverse_det(c.event_matrix(Kcheck,e['mask']))
        prob=(-1)**(n-e['mask'].bit_count())*det
        c.need(prob==evaluate(e['p']),'direct event determinant mismatch')
        ts=[c.matmul(inv,E) for E in basis]
        scores=[c.trace(t) for t in ts]
        for a in range(dim):
            c.need(prob*scores[a]==evaluate([m[a] for m in e['g']]),'direct event gradient mismatch')
            for b in range(dim):
                target=prob*(scores[a]*scores[b]-c.trace_product(ts[a],ts[b]))
                c.need(target==evaluate([m[a][b] for m in e['H']]),'direct event second-jet mismatch')
        rb=c.matmul(inv,B)
        c.need(prob*c.trace_product(rb,rb)==evaluate(e['im']),'direct imaginary second-jet mismatch')

    c.need(sum(e['a'] for e in events if e['w']==1)==n,'one-defect mass is not n')
    imaginary_fourth_leading=3*sum(e['im'][0]**2/e['a'] for e in events if e['w']==2)
    c.need(imaginary_fourth_leading>0,'nonnegative imaginary fourth leading coefficient')
    log2=c.atanh_log(F(2))
    for e in events:e['loga']=c.log_interval(e['a'],log2)
    limit=c.ZERO
    logcoef=sum(e['w']*e['im'][0] for e in events)
    c.need(logcoef==0,'nonzero logarithmic divergence in commuting imaginary direction')
    for e in events:limit=c.add(limit,c.scale(e['loga'],-e['im'][0]))
    c.need(limit[0]>F(1,30) and limit[1]<F(1,29),'positive limit certificate failed')

    # All 225 entries of the logarithmic quadratic form are checked.
    J=zeros(dim,dim)
    for e in events:J=madd(J,mscale(e['H'][0],F(e['w'])))
    pa=[c.trace_product(P,E) for E in basis]
    qa=[c.trace_product(Q,E) for E in basis]
    for a in range(dim):
        for b in range(dim):
            target=2*(pa[a]*qa[b]+pa[b]*qa[a])-4*c.trace(c.matmul(c.matmul(c.matmul(P,basis[a]),Q),basis[b]))
            c.need(J[a][b]==target,'logarithmic real Hessian form mismatch')

    # Bounds sufficient uniformly for the WHOLE interval, without evaluating
    # probabilities at epsilon=2^(-N).
    safe=F(1,4)
    for e in events:
        if e['R']:safe=min(safe,e['a']/(2*e['R']))
    minedge=min(abs(P[i][j]) for i in range(n) for j in range(i))
    c.need(minedge>0,'dense-Jacobian coercivity bound unavailable')
    g_lower=1/(4*n*n+9*n*n*(n-1)/(minedge*minedge))
    R1=F(0);Aerr=F(0);Berr=F(0);Ai=F(0);Bi=F(0)
    Q0=[[c.ZERO for _ in range(dim)] for _ in range(dim)]
    for e in events:
        abslp=max(abs(e['loga'][0]),abs(e['loga'][1]))
        if e['w']==1:
            cg=sum((abs(x) for row in e['g'][1:] for x in row),F(0))
            R1+=2*cg*cg/e['a']
        hh=[row_bound(m) for m in e['H']]
        aa=sum(hh[1:],F(0));bb=sum(hh,F(0))
        Aerr+=e['w']*aa
        Berr+=aa*abslp+2*bb*e['R']/e['a']
        ii=sum(map(abs,e['im'][1:]),F(0));jj=sum(map(abs,e['im']),F(0))
        Ai+=e['w']*ii;Bi+=ii*abslp+2*jj*e['R']/e['a']
        for a in range(dim):
            for b in range(dim):Q0[a][b]=c.add(Q0[a][b],c.scale(e['loga'],-e['H'][0][a][b]))
    q0=max(sum((max(abs(iv[0]),abs(iv[1])) for iv in row),F(0)) for row in Q0)
    N=2
    while True:
        eps=F(1,2**N)
        delta=eps*(Aerr*N+Berr+R1)
        imag_error=eps*(Ai*N+Bi)
        if (eps<=safe and delta<1 and F(N)>q0+delta+1
                and eps*(6*N+q0+delta+1)<g_lower/4 and imag_error<F(1,120)):
            break
        N*=2
        c.need(N<=16384,'bounds unexpectedly too large')
    # log(2) in (1/2,1) implies N/2 < log(1/eps0) < N.
    c.need(log2[0]>F(1,2) and log2[1]<1,'binary logarithm bound')

    # A small exact check of the additional mixed-fourth identity.
    K1=[[F(1,3),F(1,10)],[F(1,10),F(2,5)]]
    K2=[[F(1,2),F(1,7)],[F(1,7),F(3,5)]]
    X=[[F(1),F(2)],[F(-1),F(1)]]
    Y=[[F(2),F(0)],[F(1),F(-1)]]
    def block_events(k: Matrix) -> list[tuple[F,Matrix]]:
        out=[]
        for mask in range(4):
            inv,dd=c.inverse_det(c.event_matrix(k,mask))
            pp=(-1)**(2-mask.bit_count())*dd
            c.need(pp>0,'mixed fourth check: nonpositive block probability')
            out.append((pp,inv))
        for i in range(2):
            for j in range(2):
                for k2 in range(2):
                    for ell in range(2):
                        lhs=sum(pp*R[i][j]*R[k2][ell] for pp,R in out)
                        rhs=sum(pp*R[i][ell]*R[k2][j] for pp,R in out)
                        c.need(lhs==rhs,'Fisher tensor not fully symmetric')
        return out
    mixed_product=F(0);mixed_square=F(0)
    for pp,R in block_events(K1):
        for qq,T in block_events(K2):
            def bil(a: Matrix,b: Matrix) -> F:
                return c.trace(c.matmul(c.matmul(c.matmul(R,a),T),transpose(b)))
            mixed_product+=pp*qq*bil(X,X)*bil(Y,Y)
            mixed_square+=pp*qq*bil(X,Y)**2
    c.need(mixed_product==mixed_square and mixed_square>0,'mixed quartic identity failed')

    payload={
      'status':'EXACT_AUTHOR_SELF_CHECK_PASSED; FRESH_INDEPENDENT_REVIEW_PENDING',
      'base_commit':c.BASE,'previous_mechanism_commit':'a296f7cc555312e1dbd5a28af08836c4756659c7',
      'construction':'V=old integer K[:,0:2]; P=V(V^T V)^(-1)V^T; B=P B_old P+(I-P)B_old(I-P)',
      'V_integer':[[int(x) for x in row] for row in V],
      'P_rational':[[show(x) for x in row] for row in P],
      'B_parallel_rational':[[show(x) for x in row] for row in B],
      'exact_checks':{'projection':True,'rank':r,'skew':True,'commuting':True,
        'full_spark_and_all_leading_coefficients_positive':True,'all_polynomial_normalizations':True,
        'all_225_logarithmic_Hessian_entries':True,'direct_event_inverse_trace_cross_check_at_epsilon_1_over_4':True,'imaginary_logarithmic_coefficient':show(logcoef)},
      'imaginary_Hessian_limit':pack_interval(limit),
      'limit_epsilon_squared_times_imaginary_H4':show(-imaginary_fourth_leading),
      'uniform_parameter_range':{'lower':'0 (strict)','upper':f'2^(-{N})','binary_exponent_N':N},
      'uniform_conclusions':{'real':'D2 H(K_epsilon)[C,C] < -||C||_F^2 for every nonzero real symmetric C',
        'imaginary':'D2 H(K_epsilon)[i B_parallel,i B_parallel] > 1/40',
        'domain':'0 < epsilon <= 2^(-N); K_epsilon=epsilon I+(1-2 epsilon)P'},
      'bound_constants':{name:show(val) for name,val in dict(safe_epsilon=safe,min_nonzero_offdiag=minedge,
         real_Fisher_coercivity_g=g_lower,Q0_operator_upper=q0,first_jet_error_R1=R1,
         real_log_error_A=Aerr,real_log_error_B=Berr,imaginary_error_A=Ai,imaginary_error_B=Bi,
         endpoint_error_delta=delta,endpoint_imaginary_error=imag_error).items()},
      'events':[{'mask':e['mask'],'vanishing_order':e['w'],'leading_probability':show(e['a']),
         'boundary_imaginary_second_derivative':show(e['im'][0])} for e in events],
      'nongeneric_projection_exact_test':check_nongeneric_stratum(),
      'mixed_fourth_exact_test':{'E_qXX_qYY':show(mixed_product),'E_qXY_squared':show(mixed_square),
        'D4_H_cross_X_X_Y_Y':show(-12*mixed_square)},
      'scope':'All-parameter certificate for this rational family; the general projection theorem has an analytic proof. No global real-concavity theorem or real counterexample is asserted.'}
    output.write_text(json.dumps(payload,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'N':N,'limit_interval':pack_interval(limit),'q0_upper_ceiling':c.ceil(q0),
                      'g_lower':show(g_lower),'status':'PASS'},ensure_ascii=False,indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path(__file__).with_name('projection_boundary_certificate.json'))
    main(parser.parse_args().output)
