# Author: Carles Marin (Claude, Anthropic, as AI research assistant).
# STRUCTURAL question: is joint anomaly+tadpole cancellation for (3,60) EXISTENCE
# (tadpole freely adjustable among anomaly-neutral additions) or OBSTRUCTION (a
# conserved invariant)? Decided by  rank(Theta | ker A) = rank[A;Theta] - rank A
# over the full library of SU(4) additions. Exact rational linear algebra.

A3 = WeylCharacterRing("A3", style="coroots")
def su2_decompose(mb):
    rem=dict(mb); out=[]
    while any(v>0 for v in rem.values()):
        top=max(k for k,v in rem.items() if v>0); m=rem[top]
        for m2 in range(top,-top-1,-2):
            assert rem.get(m2,0)>=m; rem[m2]-=m
        out.extend([top+1]*m)
    return out
def ew_modes(labels):
    a,b,c=labels; boxes=a+2*b+3*c; buckets={}
    for w,mult in A3(*labels).weight_multiplicities().items():
        n1,n2,n3,n4=[ZZ(w[i]+QQ(boxes)/4) for i in range(4)]
        key=((-1)^n4,(-1)^(n3+n4),n1+n2-2*n3,n1+n2+n3-3*n4)
        buckets.setdefault(key,{})[n1-n2]=buckets.setdefault(key,{}).get(n1-n2,0)+mult
    out=[]
    for (p0,p2,q8,q15),bym in buckets.items():
        for d in su2_decompose(bym): out.append((d,p0,p2,q8,q15))
    return out
def hyp(q8,q15): return QQ(-7*q8+4*q15)/18
def T2(d): return QQ(d*(d^2-1))/12
H1=(1,-1,0,0); H2=(1,1,-1,-1)
def trace4(labels,H): return sum(mult*sum(w[i]*H[i] for i in range(4))^4 for w,mult in A3(*labels).weight_multiplicities().items())
def A4(labels):
    r1,r2=trace4(labels,H1),trace4(labels,H2); B=(r2-2*r1)/8
    return ZZ((r1-4*B)/2)
COLORS=[("1",1,0),("3",3,1),("3bar",3,-1)]; ETAS=[(1,1),(1,-1),(-1,1),(-1,-1)]
def zero_anom(labels,color,eta,chi6):
    cname,cdim,cubic=color; out=[QQ(0)]*5
    for d,p0,p2,q8,q15 in ew_modes(labels):
        if (p0,p2)==eta: sign=QQ(chi6)
        elif (p0,p2)==(-eta[0],-eta[1]): sign=QQ(-chi6)
        else: continue
        y=hyp(q8,q15)
        out[0]+=sign*cubic*d; out[1]+=sign*QQ(d)/2*y if cdim==3 else 0
        out[2]+=sign*cdim*T2(d)*y; out[3]+=sign*cdim*d*y^3; out[4]+=sign*cdim*d*y
    return out
def tadpole(labels,color,eta):
    cname,cdim,cubic=color; a,b,c=labels; boxes=a+2*b+3*c; t01=ZZ(0); t23=ZZ(0)
    for w,mult in A3(*labels).weight_multiplicities().items():
        n1,n2,n3,n4=[ZZ(w[i]+QQ(boxes)/4) for i in range(4)]
        p0=(-1)^n4; p2=(-1)^(n3+n4); q15=n1+n2+n3-3*n4; qP=n1+n2-n3-n4
        t01 += mult*eta[0]*p0*q15*cdim; t23 += mult*eta[1]*p2*qP*cdim
    return (ZZ(t01),ZZ(t23))
def vec8(labels,color,eta,chi6):
    cname,cdim,cubic=color; a4=A4(labels); acolor=ZZ(0) if cname=="1" else ZZ(1)
    z=zero_anom(labels,color,eta,chi6)
    return [QQ(chi6*acolor*A3(*labels).degree()),QQ(chi6*cdim*a4),QQ(chi6*cdim*A3(*labels).degree())]+[QQ(x) for x in z]

# build library
rows_anom=[]; rows_tad=[]
for a in range(6):
  for b in range(5):
    for c in range(6):
      labels=(a,b,c);
      if labels==(0,0,0) or A3(*labels).degree()>35: continue
      for color in COLORS:
        for eta in ETAS:
          for chi6 in [-1,1]:
            v=vec8(labels,color,eta,chi6)
            if all(x==0 for x in v) and tadpole(labels,color,eta)==(0,0): continue
            rows_anom.append(v)
            th=tadpole(labels,color,eta); rows_tad.append([QQ(th[0]),QQ(th[1])])
N=len(rows_anom); print("library species:",N)
A=Matrix(QQ,rows_anom).transpose()      # 8 x N
T=Matrix(QQ,rows_tad).transpose()        # 2 x N
rA=A.rank(); rAT=A.stack(T).rank()
free=rAT-rA
print("rank A            :",rA)
print("rank [A;Theta]    :",rAT)
print("rank(Theta|ker A) :",free,"  (2=EXISTENCE-flavored, 1=OBSTRUCTION invariant, 0=tadpole determined)")

# base (3,60) vectors
base_v=vec8((0,2,1),("3",3,1),(1,1),1)
base_th=tadpole((0,2,1),("3",3,1),(1,1))
print("\nbase anomalies:",base_v)
print("base tadpole  :",base_th)

if free<2:
    # find conserved invariant: covector (a|c) on (anomaly|tadpole) with a*A + c*T = 0 on every species,
    # i.e. left null vector of stacked [A;T] with nonzero tadpole part c. Then phi=c*Theta is invariant.
    M=A.stack(T)                          # 10 x N
    LK=M.left_kernel().basis()            # covectors killing every species
    print("\nleft-kernel (conserved covectors on [anom(8)|tad(2)]):")
    for u in LK:
        anom_part=u[:8]; tad_part=u[8:]
        if any(x!=0 for x in tad_part):
            print("  tad-covector c =",list(tad_part)," anom-part =",list(anom_part))
            phi_base = tad_part[0]*base_th[0]+tad_part[1]*base_th[1] + sum(anom_part[i]*base_v[i] for i in range(8))
            print("    invariant value on base (must be 0 for a joint completion to exist):",phi_base)
else:
    print("\n=> Theta is full-rank on ker A: over Q, ANY anomaly-free content is tadpole-completable.")
    print("   Only obstruction possible = the non-negative-integer cone. The 45-witness realizes it for (3,60).")
