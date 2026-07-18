# Author: Carles Marin (Claude as AI assistant).
# STRUCTURAL compatibility test (critical path, ahead of N7): does the (3,60) admit
# the SM quark content Q(3,2,Y=1/6)_L, u_R(3,1,2/3), d_R(3,1,-1/3) under SOME global
# unbroken hypercharge Y = a*q8 + b*q15 and a consistent chirality/parity split?
# The single (3,60) 6D field (eta=(+,+), chi6=+) yields zero modes: (p0,p2)=(+,+)
# with one 4D chirality, (p0,p2)=(-,-) with the opposite.  We enumerate both,
# then scan for (a,b) making Q,u_R,d_R appear with SM charges.  If NO global Y
# works, that is a STRUCTURAL NO-GO for (3,60) as the quark embedding -- stronger
# and earlier than the brane matching (Codex's point).

A3 = WeylCharacterRing("A3", style="coroots")
def su2_decompose(mb):
    rem=dict(mb); out=[]
    while any(v>0 for v in rem.values()):
        top=max(k for k,v in rem.items() if v>0); m=rem[top]
        for m2 in range(top,-top-1,-2):
            assert rem.get(m2,0)>=m; rem[m2]-=m
        out.extend([top+1]*m)
    return out

labels=(0,2,1); boxes=sum((i+1)*labels[i] for i in range(3))
# collect modes by (p0,p2); LH = (+,+) [chi +1], RH = (-,-) [chi -1]
def modes_of(parity):
    buckets={}
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

LH=modes_of((1,1)); RH=modes_of((-1,-1))
print("(3,60) zero modes:")
print("  LH (d,q8,q15):", sorted(LH))
print("  RH (d,q8,q15):", sorted(RH))
LH_doublets=[m for m in LH if m[0]==2]; RH_singlets=[m for m in RH if m[0]==1]
LH_singlets=[m for m in LH if m[0]==1]; RH_doublets=[m for m in RH if m[0]==2]
print("  LH doublets (Q candidates): %d ; RH singlets (u_R,d_R candidates): %d" %
      (len(LH_doublets), len(RH_singlets)))
print("  LH singlets: %d ; RH doublets: %d" % (len(LH_singlets), len(RH_doublets)))

if not RH_singlets:
    print("\nNO RH singlets in the (3,60) -> cannot provide u_R,d_R as SU(2)_L singlets.")
    print("STRUCTURAL NO-GO candidate: check the mirror chirality assignment before concluding.")

# scan global Y = a*q8 + b*q15 : need a LH doublet with Y=1/6, and two RH singlets
# with Y=2/3 (u_R) and Y=-1/3 (d_R).  Solve (a,b) from Q & u_R, then test d_R.
found=[]
for (dq,q8Q,q15Q) in LH_doublets:
    for (du,q8u,q15u) in RH_singlets:
        M=Matrix(QQ,[[q8Q,q15Q],[q8u,q15u]]); rhs=vector(QQ,[QQ(1)/6,QQ(2)/3])
        if M.det()==0: continue
        a,b=M.solve_right(rhs)
        # does some RH singlet give d_R with Y=-1/3 under this (a,b)?
        for (dd,q8d,q15d) in RH_singlets:
            if a*q8d+b*q15d==QQ(-1)/3 and (dd,q8d,q15d)!=(du,q8u,q15u):
                found.append((a,b,(q8Q,q15Q),(q8u,q15u),(q8d,q15d)))
if found:
    print("\nCOMPATIBLE: exists global Y=a*q8+b*q15 giving SM Q,u_R,d_R:")
    for a,b,Q,u,d in found[:3]:
        print("   Y=%s*q8+%s*q15 ; Q at q8,q15=%s ; u_R=%s ; d_R=%s" % (a,b,Q,u,d))
    print("   => (3,60) CAN embed SM quarks; N7 proceeds with this Y.")
else:
    print("\nNO global Y=a*q8+b*q15 yields Q(2,1/6)+u_R(1,2/3)+d_R(1,-1/3) from the (3,60).")
    print("=> STRUCTURAL NO-GO (candidate): the (3,60) cannot host the SM quark content")
    print("   under any single unbroken hypercharge -- earlier & stronger than brane N7.")
    print("   CAVEAT: this uses Y in the (q8,q15) plane; if AHMN's Y needs a 3rd Cartan")
    print("   direction (qP) or a different chirality map, hand that to refine before final.")
