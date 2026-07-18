# Author: Carles Marin (Claude, Anthropic, as AI research assistant).
# Publication-grade figures for the GHU (3,60) completion paper. Three figures:
#  (A) anomaly cancellation base->completion->0 (8 coefficients, normalized);
#  (B) localized-tadpole journey on the (Theta01,Theta23) plane;
#  (C) the 45-multiplet spectrum, grouped by color representation.
# Colorblind-safe Okabe-Ito palette; vector PDF for the manuscript + 400 dpi PNG.
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from pathlib import Path

OUT = Path(r'E:/proyectos/Curiosity/research/smeft_formalization/paper')
plt.rcParams.update({
    'font.family': 'serif', 'font.size': 11, 'mathtext.fontset': 'cm',
    'axes.edgecolor': '#333333', 'axes.linewidth': 0.8, 'figure.dpi': 120,
})
OK = {'blue':'#0072B2','orange':'#E69F00','pink':'#CC79A7','green':'#009E73',
      'red':'#D55E00','grey':'#7A7A7A','sky':'#56B4E9','yellow':'#F0E442'}

# ---------------------------------------------------------------- Figure A
# 8 bulk/zero-mode anomaly coefficients of the base (3,60) block; the
# completion contributes exactly the negative of each, so the total is 0.
names = [r'$SU(3)_C^4$', r'$SU(4)^4$', r'$\mathrm{grav}_6$', r'$SU3^3$',
         r'$SU3^2Y$', r'$SU2_L^2Y$', r'$Y^3$', r'$\mathrm{grav}^2Y$']
blab = ['$+60$', '$-579$', '$+180$', '$-3$', r'$+\frac{17}{6}$', r'$+\frac{23}{2}$',
        r'$+\frac{2293}{18}$', '$+17$']           # exact base coefficients (physical Y)
fig, ax = plt.subplots(figsize=(8.6, 4.6))
y = list(range(len(names)))
# mirror design: base contribution to the right (+1), completion to the left (-1);
# equal magnitude => their sum at the centre is exactly zero.
ax.barh(y, [1]*len(y), color=OK['blue'], edgecolor='black', lw=.5, height=.60,
        label='base $(3,\\mathbf{60})$ coefficient', zorder=3)
ax.barh(y, [-1]*len(y), color=OK['orange'], edgecolor='black', lw=.5, height=.60,
        label='45-multiplet completion (cancels it)', zorder=3, alpha=.96)
for i, lb in enumerate(blab):
    ax.text(1.04, i, lb, ha='left', va='center', fontsize=9, color=OK['blue'])
    ax.text(-1.04, i, r'$-(\,$'+lb+r'$\,)$', ha='right', va='center', fontsize=8, color=OK['orange'])
ax.axvline(0, color='black', lw=1.4, zorder=5)
ax.scatter([0]*len(y), y, s=70, color=OK['green'], zorder=6, marker='o',
           edgecolor='black', lw=.7)
ax.set_yticks(y); ax.set_yticklabels(names, fontsize=12)
ax.set_xlim(-1.9, 1.9); ax.set_xticks([-1, 0, 1])
ax.set_xticklabels(['completion', r'sum $=\mathbf{0}$', 'base'], fontsize=10)
ax.invert_yaxis()
ax.set_title(r'All eight anomaly coefficients cancel exactly (green $=0$)', fontsize=12)
ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.20), ncol=2, frameon=True, fontsize=9)
for s in ('top', 'right', 'left'): ax.spines[s].set_visible(False)
ax.tick_params(axis='y', length=0)
fig.tight_layout(); fig.savefig(OUT/'fig_anomaly_cancel.pdf')
fig.savefig(OUT/'fig_anomaly_cancel.png', dpi=400); plt.close(fig)

# ---------------------------------------------------------------- Figure B
# The localized-tadpole plane. base=(126,-36); the 8-anomaly witness of 27
# multiplets reaches (126,-20) and STILL fails; the 45-multiplet completion
# takes the total to the origin.
fig, ax = plt.subplots(figsize=(6.6, 5.6))
ax.axhspan(-4, 4, color=OK['green'], alpha=.06)
ax.axvspan(-4, 4, color=OK['green'], alpha=.06)
for gx in range(-20, 141, 20): ax.axvline(gx, color='#e8e8e8', lw=.6, zorder=0)
for gy in range(-40, 41, 20): ax.axhline(gy, color='#e8e8e8', lw=.6, zorder=0)
ax.axhline(0, color='#555', lw=.9); ax.axvline(0, color='#555', lw=.9)
# journey arrow base -> origin (the 45-multiplet completion)
ax.add_patch(FancyArrowPatch((126, -36), (2, -1), arrowstyle='-|>', mutation_scale=18,
             lw=2.0, color=OK['green'], zorder=4, shrinkA=8, shrinkB=8))
ax.text(66, -12, '45-multiplet\ncompletion', color=OK['green'], fontsize=9.5,
        rotation=14, ha='center', va='bottom', fontweight='bold')
pts = [((126, -36), 'base $(3,\\mathbf{60})$', OK['blue'], 'FAIL tadpole'),
       ((126, -20), '27-multiplet\n(8 anomalies only)', OK['red'], 'FAIL tadpole'),
       ((0, 0),   'joint target', OK['green'], 'PASS')]
for (x, y0), lab, col, tag in pts:
    ax.scatter(x, y0, s=150, color=col, edgecolor='black', lw=.9, zorder=6)
    dy = 6 if y0 <= 0 else -10
    ax.annotate(f'{lab}\n{tag} $({x},{y0})$', (x, y0), (x, y0+dy),
                fontsize=8.5, ha='center', va='bottom' if dy > 0 else 'top',
                color=col, fontweight='bold')
ax.set_xlim(-24, 150); ax.set_ylim(-50, 34)
ax.set_xlabel(r'$\Theta_{01}$  (tadpole at $f=0,1$)')
ax.set_ylabel(r'$\Theta_{23}$  (tadpole at $f=2,3$)')
ax.set_title('The localized tadpole can be cancelled jointly', fontsize=11.5)
for s in ('top', 'right'): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(OUT/'fig_tadpole_journey.pdf')
fig.savefig(OUT/'fig_tadpole_journey.png', dpi=400); plt.close(fig)

# ---------------------------------------------------------------- Figure C
# The 45-multiplet spectrum, grouped by color representation (readable rows).
witness = [
 (8,'(0,0,1)','1','--','+'),(1,'(0,0,1)','3','-+','+'),(4,'(0,0,1)','3','--','-'),
 (5,'(0,0,1)','3bar','++','+'),(3,'(0,0,1)','3bar','-+','+'),(3,'(0,0,3)','1','+-','-'),
 (3,'(0,0,3)','1','-+','-'),(2,'(0,0,3)','3','--','+'),(1,'(0,0,3)','3bar','+-','-'),
 (1,'(0,0,3)','3bar','-+','-'),(2,'(0,0,4)','1','-+','+'),(2,'(0,1,0)','1','++','-'),
 (2,'(0,1,1)','1','++','-'),(2,'(0,1,1)','1','++','+'),(2,'(0,1,1)','3bar','--','-'),
 (2,'(0,2,0)','3','-+','-'),(2,'(1,0,1)','1','++','+')]
order = {'1':0,'3':1,'3bar':2}
cmap = {'1':OK['blue'],'3':OK['orange'],'3bar':OK['pink']}
cname = {'1':r'color singlet $\mathbf{1}$','3':r'triplet $\mathbf{3}$','3bar':r'antitriplet $\bar{\mathbf{3}}$'}
ws = sorted(witness, key=lambda t: (order[t[2]], -t[0]))
fig, ax = plt.subplots(figsize=(8.4, 6.0))
ypos = list(range(len(ws)))[::-1]
totals = {'1':0,'3':0,'3bar':0}
for yv, (n, dyn, col, eta, chi) in zip(ypos, ws):
    totals[col] += n
    ax.barh(yv, n, color=cmap[col], edgecolor='black', lw=.5, height=.66, zorder=3)
    ax.text(n+.12, yv, f'{n}', va='center', ha='left', fontsize=8.5, color='#333')
labels = [rf'$\mathbf{{{dyn}}}$  $\eta={eta}$, $\chi_6={chi}$' for _, dyn, _, eta, chi in ws]
ax.set_yticks(ypos); ax.set_yticklabels(labels, fontsize=8)
ax.set_xlabel('multiplicity in the completion')
ax.set_xlim(0, 9)
handles = [plt.Rectangle((0,0),1,1,color=cmap[c],ec='black',lw=.4) for c in ('1','3','3bar')]
leg = [f'{cname[c]}  (total {totals[c]})' for c in ('1','3','3bar')]
ax.legend(handles, leg, loc='lower right', frameon=True, fontsize=9)
ax.set_title(r'The 45-multiplet completion by color rep  (total $=45$)', fontsize=11.5)
for s in ('top', 'right'): ax.spines[s].set_visible(False)
ax.grid(axis='x', alpha=.2)
fig.tight_layout(); fig.savefig(OUT/'fig_spectrum45.pdf')
fig.savefig(OUT/'fig_spectrum45.png', dpi=400); plt.close(fig)
print('WROTE 3 figures to', OUT)
print('color totals:', totals, 'sum', sum(totals.values()))
