#!/usr/bin/env python3
# Author: Carles Marin (Claude, Anthropic, as AI research assistant).
# ghu_preflight.py - an EXACT decoupling/feasibility oracle for orbifold fermion sectors.
#
# Given a library of candidate additions, each carrying a vector of "hard" anomaly
# coefficients A in Q**a and a "soft" localized functional Theta in Q**t (e.g. the fixed-
# point tadpole), and a base block, it decides -- by exact linear algebra + a Farkas /
# cone test, NO search -- whether Theta is an INDEPENDENT OBSTRUCTION or is always
# absorbable by anomaly-neutral additions:
#
#   * rank[A;Theta] - rank A  <  t   =>  a CONSERVED INVARIANT exists (phi . Theta).
#         If phi(Theta_base) != 0 for the base  =>  NO-GO certificate (skip any search).
#   * cone test: is Theta(anomaly-neutral, x>=0) = R**t ?  (all signed axes reachable)
#         yes => EXISTENCE guaranteed (an unbounded search WILL succeed);
#         no  => pointed cone: obstruction for residuals outside it (with a witness).
#
# Ideal uses: (1) PREFLIGHT before an expensive integer search -- decision-complete on
# both sides; (2) LANDSCAPE CLASSIFIER -- sweep models to locate genuinely obstructed
# ones (where anomaly-free =/=> tadpole-free, a real no-go with a low-energy signature);
# (3) CONSTRUCTIVE completions from the Farkas cone certificates.
#
# The physics model (which reps/parities/chiralities, and how A, Theta are computed) is
# a pluggable callback `library(model)` returning rows_A (list of Q**a) and rows_T
# (list of Q**t), plus base_A, base_T. Run under SageMath.  This file ships the SU(4)
# on T**2/Z2 model of the (3,60) paper as the default.

from sage.all import *  # run: sage ghu_preflight.py

def su4_t2z2_library():
    """Default model: SU(4) on T**2/Z2, additions of dim<=35, all colors/parities/chirs.
    Returns (rows_A [Q**8], rows_T [Q**2], base_A [Q**8], base_T [Q**2])."""
    A3 = WeylCharacterRing("A3", style="coroots")
    def su2_decompose(mb):
        rem=dict(mb); out=[]
        while any(v>0 for v in rem.values()):
            top=max(k for k,v in rem.items() if v>0); m=rem[top]
            for m2 in range(top,-top-1,-2): assert rem.get(m2,0)>=m; rem[m2]-=m
            out.extend([top+1]*m)
        return out
    def ew_modes(labels):
        a,b,c=labels; boxes=a+2*b+3*c; buckets={}
        for w,mult in A3(*labels).weight_multiplicities().items():
            n1,n2,n3,n4=[ZZ(w[i]+QQ(boxes)/4) for i in range(4)]
            key=((-1)**n4,(-1)**(n3+n4),n1+n2-2*n3,n1+n2+n3-3*n4)
            buckets.setdefault(key,{})[n1-n2]=buckets.setdefault(key,{}).get(n1-n2,0)+mult
        out=[]
        for (p0,p2,q8,q15),bym in buckets.items():
            for d in su2_decompose(bym): out.append((d,p0,p2,q8,q15))
        return out
    hyp=lambda q8,q15: QQ(-7*q8+4*q15)/18
    T2=lambda d: QQ(d*(d**2-1))/12
    H1,H2=(1,-1,0,0),(1,1,-1,-1)
    tr4=lambda L,H: sum(m*sum(w[i]*H[i] for i in range(4))**4 for w,m in A3(*L).weight_multiplicities().items())
    def A4(L):
        r1,r2=tr4(L,H1),tr4(L,H2); return ZZ((r1-4*((r2-2*r1)/8))/2)
    def zero_anom(L,col,eta,chi6):
        _,cdim,cubic=col; out=[QQ(0)]*5
        for d,p0,p2,q8,q15 in ew_modes(L):
            if (p0,p2)==eta: s=QQ(chi6)
            elif (p0,p2)==(-eta[0],-eta[1]): s=QQ(-chi6)
            else: continue
            y=hyp(q8,q15)
            out[0]+=s*cubic*d; out[1]+=s*QQ(d)/2*y if cdim==3 else 0
            out[2]+=s*cdim*T2(d)*y; out[3]+=s*cdim*d*y**3; out[4]+=s*cdim*d*y
        return out
    def tad(L,col,eta):
        _,cdim,cubic=col; a,b,c=L; boxes=a+2*b+3*c; t01=t23=ZZ(0)
        for w,m in A3(*L).weight_multiplicities().items():
            n1,n2,n3,n4=[ZZ(w[i]+QQ(boxes)/4) for i in range(4)]
            p0,p2=(-1)**n4,(-1)**(n3+n4); q15=n1+n2+n3-3*n4; qP=n1+n2-n3-n4
            t01+=m*eta[0]*p0*q15*cdim; t23+=m*eta[1]*p2*qP*cdim
        return (QQ(t01),QQ(t23))
    def vec8(L,col,eta,chi6):
        _,cdim,cubic=col; ac=ZZ(0) if col[0]=="1" else ZZ(1)
        return [QQ(chi6*ac*A3(*L).degree()),QQ(chi6*cdim*A4(L)),QQ(chi6*cdim*A3(*L).degree())]+[QQ(x) for x in zero_anom(L,col,eta,chi6)]
    COLORS=[("1",1,0),("3",3,1),("3bar",3,-1)]; ETAS=[(1,1),(1,-1),(-1,1),(-1,-1)]
    rows_A=[]; rows_T=[]
    for a in range(6):
      for b in range(5):
        for c in range(6):
          L=(a,b,c)
          if L==(0,0,0) or A3(*L).degree()>35: continue
          for col in COLORS:
            for eta in ETAS:
              for chi6 in [-1,1]:
                v=vec8(L,col,eta,chi6); th=tad(L,col,eta)
                if all(x==0 for x in v) and th==(0,0): continue
                rows_A.append(v); rows_T.append(list(th))
    base_A=vec8((0,2,1),("3",3,1),(1,1),1); base_T=list(tad((0,2,1),("3",3,1),(1,1)))
    return rows_A, rows_T, base_A, base_T

def preflight(rows_A, rows_T, base_A, base_T):
    """The exact oracle. Prints a verdict with certificates. Returns a dict."""
    A=Matrix(QQ,rows_A).transpose(); T=Matrix(QQ,rows_T).transpose()
    a,t,N=A.nrows(),T.nrows(),A.ncols()
    rA=A.rank(); rAT=A.stack(T).rank(); free=rAT-rA
    print(f"library size {N} | anomalies {a} | soft functionals {t}")
    print(f"rank A = {rA} | rank[A;Theta] = {rAT} | free tadpole rank = {free} (=t means unobstructed)")
    out={"free":free,"t":t,"verdict":None}
    if free<t:                                     # conserved invariant(s)
        M=A.stack(T); inv=[]
        for u in M.left_kernel().basis():
            ap,tp=u[:a],u[a:]
            if any(x!=0 for x in tp):
                val=sum(tp[i]*base_T[i] for i in range(t))+sum(ap[i]*base_A[i] for i in range(a))
                inv.append((list(tp),list(ap),val))
        blocking=[iv for iv in inv if iv[2]!=0]
        if blocking:
            print("VERDICT: NO-GO. conserved invariant phi.Theta blocks joint cancellation:")
            for tp,ap,val in blocking: print(f"   phi_Theta={tp}  value on base={val} (!=0)")
            out["verdict"]="NO-GO"; out["invariants"]=blocking
        else:
            print("VERDICT: constrained but SATISFIABLE (invariant vanishes on base); tadpole on a sublattice.")
            out["verdict"]="CONSTRAINED"; out["invariants"]=inv
        return out
    # full rank: cone (Farkas) reachability of every signed axis by anomaly-neutral x>=0
    reach={}
    for axis in range(t):
        for sign in (+1,-1):
            lp=MixedIntegerLinearProgram(maximization=False,solver="GLPK")
            x=lp.new_variable(nonnegative=True,real=True)
            for k in range(a): lp.add_constraint(sum(A[k,j]*x[j] for j in range(N))==0)
            for o in range(t):
                if o!=axis: lp.add_constraint(sum(T[o,j]*x[j] for j in range(N))==0)
            lp.add_constraint(sign*sum(T[axis,j]*x[j] for j in range(N))>=1)
            try: lp.solve(); reach[(axis,sign)]=True
            except Exception: reach[(axis,sign)]=False
    full=all(reach.values())
    print("cone reachability (anomaly-neutral, x>=0):", {f"{'+-'[s<0]}e{ax}":ok for (ax,s),ok in reach.items()})
    if full:
        print("VERDICT: UNOBSTRUCTED. Theta(anomaly-neutral cone)=R**t; every anomaly-free")
        print("         completion is tadpole-compatible. An unbounded search WILL succeed.")
        out["verdict"]="UNOBSTRUCTED"
    else:
        print("VERDICT: CONE-OBSTRUCTED. pointed cone; residuals outside it are unreachable.")
        out["verdict"]="CONE-OBSTRUCTED"
    out["reach"]=reach
    return out

if __name__=="__main__":
    print("== GHU tadpole/anomaly preflight oracle ==  (default model: SU(4) on T^2/Z2, the (3,60) paper)")
    rows_A, rows_T, base_A, base_T = su4_t2z2_library()
    preflight(rows_A, rows_T, base_A, base_T)
