# Author: Carles Marin (Claude as AI assistant).
# Method pivot: OR-Tools CP-SAT (exact, integer-native, multithread) for the
# bounded (3,60) completion — the verdict PPL's exact B&B could not deliver in
# 3000 s.  Branch (a), physical per Codex's GS gate: impose SU(3)_C^4, integrality
# of the two f=0,1 non-Abelian cubics.  Same problem as r60_completion_ppl_v2.sage.
# A feasible witness is re-verified exactly in Sage rationals.

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

def hyp(q8,q15): return QQ(-7*q8+4*q15)/18   # physical hypercharge (paper convention)
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
        p0=(-1)^n4; p2=(-1)^(n3+n4)
        q15=n1+n2+n3-3*n4; qP=n1+n2-n3-n4
        t01 += mult*eta[0]*p0*q15*cdim
        t23 += mult*eta[1]*p2*qP*cdim
    return (ZZ(t01),ZZ(t23))

def vec(labels,color,eta,chi6):
    cname,cdim,cubic=color; a4=A4(labels); acolor=ZZ(0) if cname=="1" else ZZ(1)
    z=zero_anom(labels,color,eta,chi6)
    return [ZZ(chi6*acolor*A3(*labels).degree()),ZZ(chi6*cdim*a4),ZZ(chi6*cdim*A3(*labels).degree())]+z

base=vec((0,2,1),("3",3,1),(1,1),1); base_local=p0_loc((0,2,1),("3",3,1),(1,1),1)
base_theta=tadpole((0,2,1),("3",3,1),(1,1))
assert base==[ZZ(60),ZZ(-579),ZZ(180),QQ(-3),QQ(17)/6,QQ(23)/2,QQ(2293)/18,QQ(17)]
assert base_local==[QQ(-3)/2,QQ(-3)/4]

lib=[]
for a in range(6):
  for b in range(5):
    for c in range(6):
      labels=(a,b,c); dim=A3(*labels).degree()
      if labels==(0,0,0) or dim>35: continue
      for color in COLORS:
        for eta in ETAS:
          for chi6 in [-1,1]:
            v=vec(labels,color,eta,chi6); loc=p0_loc(labels,color,eta,chi6)
            if any(x!=0 for x in v): lib.append((labels,color,eta,chi6,v,loc))
uniq={}
for it in lib:
    key=tuple(it[4])+tuple(it[5])+tuple(tadpole(it[0],it[1],it[2]))
    if key not in uniq: uniq[key]=it
lib=list(uniq.values()); N=len(lib)
print("library=%d" % N)

# scale 8 equalities to integers
den=lcm([QQ(x).denominator() for x in base]+[QQ(it[4][k]).denominator() for it in lib for k in range(8)])
tgt=[ZZ(den*QQ(x)) for x in base]
A=[[ZZ(den*QQ(it[4][k])) for it in lib] for k in range(8)]
# local: 4*local integer
dl=4
locb=[ZZ(dl*QQ(x)) for x in base_local]
Aloc=[[ZZ(dl*QQ(it[5][m])) for it in lib] for m in range(2)]


# Independent verification of the 45-multiplet JOINT witness (anomaly+tadpole+integrality)
W=[(8,(0,0,1),"1",(-1,-1),1),(1,(0,0,1),"3",(-1,1),1),(4,(0,0,1),"3",(-1,-1),-1),
   (5,(0,0,1),"3bar",(1,1),1),(3,(0,0,1),"3bar",(-1,1),1),(3,(0,0,3),"1",(1,-1),-1),
   (3,(0,0,3),"1",(-1,1),-1),(2,(0,0,3),"3",(-1,-1),1),(1,(0,0,3),"3bar",(1,-1),-1),
   (1,(0,0,3),"3bar",(-1,1),-1),(2,(0,0,4),"1",(-1,1),1),(2,(0,1,0),"1",(1,1),-1),
   (2,(0,1,1),"1",(1,1),-1),(2,(0,1,1),"1",(1,1),1),(2,(0,1,1),"3bar",(-1,-1),-1),
   (2,(0,2,0),"3",(-1,1),-1),(2,(1,0,1),"1",(1,1),1)]
cmap={"1":("1",1,0),"3":("3",3,1),"3bar":("3bar",3,-1)}
totv=[QQ(0)]*8; totth=[0,0]; totloc=[QQ(0),QQ(0)]
for n,lab,cn,eta,chi in W:
    col=cmap[cn]; v=vec(lab,col,eta,chi); th=tadpole(lab,col,eta); lc=p0_loc(lab,col,eta,chi)
    for k in range(8): totv[k]+=n*QQ(v[k])
    totth[0]+=n*int(th[0]); totth[1]+=n*int(th[1])
    for m in range(2): totloc[m]+=n*QQ(lc[m])
print("total multiplets:", sum(n for n,*_ in W))
okA=all(totv[k]+QQ(base[k])==0 for k in range(8))
okT=(totth[0]+base_theta[0]==0 and totth[1]+base_theta[1]==0)
resloc=[totloc[m]+QQ(base_local[m]) for m in range(2)]
okL=all(x in ZZ for x in resloc)
print("(1) 8 anomalies cancel (Sum+base=0)? %s   -> %s" % (okA,[totv[k]+QQ(base[k]) for k in range(8)]))
print("(2) TADPOLE cancels (Sum theta + base_theta=0)? %s  -> Sum theta=%s, base_theta=%s" % (okT,tuple(totth),base_theta))
print("(3) local cubic residue INTEGER? %s  -> residue=%s" % (okL,resloc))
print("")
if okA and okT and okL:
    print("*** JOINT WITNESS CONFIRMED (independent verification) ***")
    print("=> OPTION 2 CLOSED: (3,60) HAS a completion cancelling 8 anomalies + tadpole")
    print("   + integer local cubics. 45 added multiplets. Next physical filter = N7.")
else:
    print("!!! verification FAILED on one constraint -> the 'witness' does not hold. Do NOT claim.")
