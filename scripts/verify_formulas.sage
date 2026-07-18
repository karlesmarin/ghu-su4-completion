# Author: Carles Marin (Claude, Anthropic, as AI research assistant).
# FULL FORMULA AUDIT for the GHU (3,60) paper. Recomputes every number the paper
# reports, in BOTH the /9 (script) and /18 (physical) hypercharge normalizations,
# and re-verifies the 45-multiplet joint witness cancellation. Prints, asserts nothing
# hardcoded, so the truth is on screen.

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
    assert sum(x[0] for x in out)==A3(*labels).degree()
    return out

DIV=[9]  # set at call time
def hyp(q8,q15,div): return QQ(-7*q8+4*q15)/div
def T2(d): return QQ(d*(d^2-1))/12
H1=(1,-1,0,0); H2=(1,1,-1,-1)
def trace4(labels,H): return sum(mult*sum(w[i]*H[i] for i in range(4))^4 for w,mult in A3(*labels).weight_multiplicities().items())
def A4(labels):
    r1,r2=trace4(labels,H1),trace4(labels,H2); B=(r2-2*r1)/8
    return ZZ((r1-4*B)/2)
COLORS=[("1",1,0),("3",3,1),("3bar",3,-1)]

def zero_anom(labels,color,eta,chi6,div):
    cname,cdim,cubic=color; out=[QQ(0)]*5
    for d,p0,p2,q8,q15 in ew_modes(labels):
        if (p0,p2)==eta: sign=QQ(chi6)
        elif (p0,p2)==(-eta[0],-eta[1]): sign=QQ(-chi6)
        else: continue
        y=hyp(q8,q15,div)
        out[0]+=sign*cubic*d; out[1]+=sign*QQ(d)/2*y if cdim==3 else 0
        out[2]+=sign*cdim*T2(d)*y; out[3]+=sign*cdim*d*y^3; out[4]+=sign*cdim*d*y
    return out

def p0_loc(labels,color,eta,chi6):
    cname,cdim,cubic=color; a,b,c=labels; boxes=a+2*b+3*c; c3=s3=QQ(0)
    for w,mult in A3(*labels).weight_multiplicities().items():
        n1,n2,n3,n4=[ZZ(w[i]+QQ(boxes)/4) for i in range(4)]
        p0=(-1)^n4; h3=n1+n2-2*n3; s=QQ(chi6*eta[0]*p0)*mult/4
        c3+=cubic*s; s3+=cdim*s*h3^3/(-6)
    return [c3,s3]

def tadpole(labels,color,eta):
    cname,cdim,cubic=color; a,b,c=labels; boxes=a+2*b+3*c; t01=ZZ(0); t23=ZZ(0)
    for w,mult in A3(*labels).weight_multiplicities().items():
        n1,n2,n3,n4=[ZZ(w[i]+QQ(boxes)/4) for i in range(4)]
        p0=(-1)^n4; p2=(-1)^(n3+n4); q15=n1+n2+n3-3*n4; qP=n1+n2-n3-n4
        t01 += mult*eta[0]*p0*q15*cdim; t23 += mult*eta[1]*p2*qP*cdim
    return (ZZ(t01),ZZ(t23))

def vec(labels,color,eta,chi6,div):
    cname,cdim,cubic=color; a4=A4(labels); acolor=ZZ(0) if cname=="1" else ZZ(1)
    z=zero_anom(labels,color,eta,chi6,div)
    return [ZZ(chi6*acolor*A3(*labels).degree()),ZZ(chi6*cdim*a4),ZZ(chi6*cdim*A3(*labels).degree())]+z

print("="*70)
print("1) SM QUARK CHARGES under physical Y=(-7q8+4q15)/18")
for name,(q8,q15) in [("Q_L",(-1,-1)),("u_R",(0,3)),("d_R",(-2,-5))]:
    print("   %-4s (q8,q15)=(%2d,%2d)  Y=%s" % (name,q8,q15,hyp(q8,q15,18)))

print("="*70)
print("2) BASE (3,60) 8 anomaly coefficients, both normalizations")
labels=(0,2,1); col=("3",3,1); eta=(1,1); chi=1
b9 =vec(labels,col,eta,chi,9)
b18=vec(labels,col,eta,chi,18)
names=["SU3C^4","SU4^4","grav6","SU3^3","SU3^2Y","SU2^2Y","Y^3","grav^2Y"]
for k in range(8):
    print("   %-8s  /9=%-10s   /18(physical)=%-10s" % (names[k],b9[k],b18[k]))
print("   local cubics [SU3C^3,SU3_123^3] =", p0_loc(labels,col,eta,chi))
print("   base tadpole Theta =", tadpole(labels,col,eta))

print("="*70)
print("3) 45-MULTIPLET WITNESS cancellation, in BOTH normalizations")
W=[(8,(0,0,1),"1",(-1,-1),1),(1,(0,0,1),"3",(-1,1),1),(4,(0,0,1),"3",(-1,-1),-1),
   (5,(0,0,1),"3bar",(1,1),1),(3,(0,0,1),"3bar",(-1,1),1),(3,(0,0,3),"1",(1,-1),-1),
   (3,(0,0,3),"1",(-1,1),-1),(2,(0,0,3),"3",(-1,-1),1),(1,(0,0,3),"3bar",(1,-1),-1),
   (1,(0,0,3),"3bar",(-1,1),-1),(2,(0,0,4),"1",(-1,1),1),(2,(0,1,0),"1",(1,1),-1),
   (2,(0,1,1),"1",(1,1),-1),(2,(0,1,1),"1",(1,1),1),(2,(0,1,1),"3bar",(-1,-1),-1),
   (2,(0,2,0),"3",(-1,1),-1),(2,(1,0,1),"1",(1,1),1)]
cmap={"1":("1",1,0),"3":("3",3,1),"3bar":("3bar",3,-1)}
print("   total multiplets:", sum(n for n,*_ in W))
for div in (9,18):
    base=vec(labels,col,eta,chi,div); tot=[QQ(0)]*8
    for n,lab,cn,e,ch in W:
        v=vec(lab,cmap[cn],e,ch,div)
        for k in range(8): tot[k]+=n*QQ(v[k])
    res=[tot[k]+QQ(base[k]) for k in range(8)]
    print("   div=%2d : 8 anomalies Sum+base = %s  -> %s" % (div,res,"ALL ZERO" if all(x==0 for x in res) else "NONZERO!"))
# tadpole + local (Y-independent, compute once)
base_th=tadpole(labels,col,eta); base_lc=p0_loc(labels,col,eta,chi)
tth=[0,0]; tlc=[QQ(0),QQ(0)]
for n,lab,cn,e,ch in W:
    th=tadpole(lab,cmap[cn],e); lc=p0_loc(lab,cmap[cn],e,ch)
    tth[0]+=n*int(th[0]); tth[1]+=n*int(th[1])
    for m in range(2): tlc[m]+=n*QQ(lc[m])
print("   tadpole Sum+base = (%s,%s) -> %s" % (tth[0]+base_th[0],tth[1]+base_th[1],
      "ZERO" if tth[0]+base_th[0]==0 and tth[1]+base_th[1]==0 else "NONZERO!"))
res_lc=[tlc[m]+QQ(base_lc[m]) for m in range(2)]
print("   local cubic residue =", res_lc, "-> INTEGER" if all(x in ZZ for x in res_lc) else "-> NOT integer!")

print("="*70)
print("4) 6D anomaly factorization Q-matrix rank (gauge, Y-independent)")
Qm=matrix(ZZ,[[30,213],[213,864]])
print("   Q =", Qm.list(), " rank =", Qm.rank(), " det =", Qm.det())
