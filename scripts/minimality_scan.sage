# Author: Carles Marin (Claude, Anthropic, as AI research assistant).
# Minimality check for the (3,60) choice. Over all SU(4) irreps R up to dimension 63,
# with the color-triplet assignment and eta=(+,+), chi6=+ (so LH=(p0,p2)=(+,+),
# RH=(-,-), the convention of r60_sm_embedding_test.sage), test whether the T^2/Z2
# zero modes of (3,R) admit a full SM quark generation Q_L(2,1/6), u_R(1,2/3),
# d_R(1,-1/3) under a SINGLE unbroken hypercharge Y=a*q8+b*q15. Report the smallest
# such R. Only if 60 is the minimum may the paper claim (3,60) is smallest.

A3 = WeylCharacterRing("A3", style="coroots")
def su2_decompose(mb):
    rem=dict(mb); out=[]
    while any(v>0 for v in rem.values()):
        top=max(k for k,v in rem.items() if v>0); m=rem[top]
        for m2 in range(top,-top-1,-2):
            assert rem.get(m2,0)>=m; rem[m2]-=m
        out.extend([top+1]*m)
    return out

def modes_of(labels,parity):
    boxes=sum((i+1)*labels[i] for i in range(3)); buckets={}
    for w,mult in A3(*labels).weight_multiplicities().items():
        n1,n2,n3,n4=[ZZ(w[i]+QQ(boxes)/4) for i in range(4)]
        p0=(-1)^n4; p2=(-1)^(n3+n4)
        if (p0,p2)!=parity: continue
        q8=n1+n2-2*n3; q15=n1+n2+n3-3*n4
        buckets.setdefault((q8,q15),{})[n1-n2]=buckets.setdefault((q8,q15),{}).get(n1-n2,0)+mult
    out=[]
    for (q8,q15),bym in buckets.items():
        for d in su2_decompose(bym): out.append((d,q8,q15))
    return out

def admits_SM(labels):
    LH=modes_of(labels,(1,1)); RH=modes_of(labels,(-1,-1))
    LHd=[m for m in LH if m[0]==2]; RHs=[m for m in RH if m[0]==1]
    for (dq,q8Q,q15Q) in LHd:
        for (du,q8u,q15u) in RHs:
            M=Matrix(QQ,[[q8Q,q15Q],[q8u,q15u]])
            if M.det()==0: continue
            a,b=M.solve_right(vector(QQ,[QQ(1)/6,QQ(2)/3]))
            for (dd,q8d,q15d) in RHs:
                if a*q8d+b*q15d==QQ(-1)/3 and (dd,q8d,q15d)!=(du,q8u,q15u):
                    return (a,b)
    return None

hits=[]
for a in range(7):
  for b in range(6):
    for c in range(7):
      labels=(a,b,c)
      if labels==(0,0,0): continue
      dim=A3(*labels).degree()
      if dim>63: continue
      Y=admits_SM(labels)
      if Y is not None:
        hits.append((dim,labels,Y))
hits.sort()
print("SU(4) irreps (dim<=63) whose (3,R) zero modes admit a full SM quark generation")
print("under a single Y=a*q8+b*q15  [eta=(+,+), chi6=+ convention]:")
for dim,labels,Y in hits:
    tag=" <-- (3,60), OUR CHOICE" if labels==(0,2,1) else ""
    print("   dim=%-3d  Dynkin=%s  Y=(%s)q8+(%s)q15%s" % (dim,labels,Y[0],Y[1],tag))
if hits:
    print("\nSMALLEST admitting rep: dim=%d  Dynkin=%s" % (hits[0][0],hits[0][1]))
    print("Is (3,60) the minimum? %s" % ("YES" if hits[0][1]==(0,2,1) else "NO -- do NOT claim minimality"))
else:
    print("\nNO rep up to dim 63 admits the SM generation in this convention.")
