# Author: Carles Marin (Claude, Anthropic, as AI research assistant).
# Figure for Theorem 1 (Existence): the anomaly-neutral tadpole cone fills R^2.
# Data = the four Farkas witnesses (each an anomaly-neutral, non-negative addition):
#   +Theta01=(288,0) [2 multiplets], -Theta01=(-288,0) [2], +Theta23=(0,48) [4], -Theta23=(0,-16) [4].
# These +-axis generators positively span the plane, so 0 is interior and any residual
# tadpole -- in particular -Theta_base=(-126,36) -- is cancellable. EN + ES.
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from pathlib import Path
OUT=Path(r'E:/proyectos/Curiosity/research/smeft_formalization/paper')
plt.rcParams.update({'font.family':'serif','font.size':11,'mathtext.fontset':'cm','axes.linewidth':0.8})
OK={'blue':'#0072B2','orange':'#E69F00','pink':'#CC79A7','green':'#009E73','red':'#D55E00','grey':'#666'}

# witnesses: (vec, n_multiplets, label)
W=[((288,0),2,r'$+\Theta_{01}$'),((-288,0),2,r'$-\Theta_{01}$'),
   ((0,48),4,r'$+\Theta_{23}$'),((0,-16),4,r'$-\Theta_{23}$')]

def make(es=False):
    fig,ax=plt.subplots(figsize=(7.0,5.4))
    ax.set_xlim(-360,360); ax.set_ylim(-90,90)
    # whole plane reachable
    ax.axhspan(-90,90,color=OK['green'],alpha=.07,zorder=0)
    ax.axhline(0,color='#888',lw=.8); ax.axvline(0,color='#888',lw=.8)
    for gx in range(-300,301,100): ax.axvline(gx,color='#ececec',lw=.6,zorder=0)
    for gy in range(-80,81,20): ax.axhline(gy,color='#ececec',lw=.6,zorder=0)
    # the four Farkas generators (directions)
    lab_mult = "mult." if not es else "mult."
    for (vx,vy),n,nm in W:
        ax.add_patch(FancyArrowPatch((0,0),(vx,vy),arrowstyle='-|>',mutation_scale=15,
                     lw=2.2,color=OK['blue'],zorder=5,shrinkA=0,shrinkB=0))
        dx = 18 if vx>=0 else -18; dy = 10 if vy>=0 else -14
        ha='left' if vx>0 else ('right' if vx<0 else 'center')
        ax.annotate(f'{nm}\n{n} mult.',(vx,vy),(vx+dx,vy+dy),fontsize=8.5,color=OK['blue'],
                    ha=ha,va='center',fontweight='bold')
    # base residual target that must be reached, -Theta_base
    tx,ty=-126,36
    ax.scatter([tx],[ty],marker='*',s=340,color=OK['red'],edgecolor='black',lw=.7,zorder=7)
    tlab = (r'$-\Theta_{\rm base}=(-126,36)$'+'\n'+('alcanzable' if es else 'reachable'))
    ax.annotate(tlab,(tx,ty),(tx-14,ty-20),fontsize=9,color=OK['red'],fontweight='bold',ha='right',va='top')
    # a couple of interior combos (non-neg combos of generators) as small dots -> the cone is 2D-dense
    combos=[(288-288,48-16),(288,48),(-288,48),(-288,-16),(288,-16),(-126,36),(150,-40),(-200,60)]
    ax.scatter([c[0] for c in combos],[c[1] for c in combos],s=16,color=OK['green'],
               edgecolor='none',alpha=.8,zorder=4)
    ax.scatter([0],[0],s=45,color='black',zorder=6)
    ax.set_xlabel(r'$\Theta_{01}$'); ax.set_ylabel(r'$\Theta_{23}$')
    if es:
        ax.set_title('El cono del tadpole anomalía-neutro llena el plano  '+r'($=\mathbb{R}^2$)',fontsize=11.5)
        note='Generadores de Farkas (± ejes) por adiciones anomalía-neutras no-negativas;\n0 en el interior '+r'$\Rightarrow$ cualquier residuo, en particular $-\Theta_{\rm base}$, es cancelable.'
    else:
        ax.set_title('The anomaly-neutral tadpole cone fills the plane  '+r'($=\mathbb{R}^2$)',fontsize=11.5)
        note='Farkas generators ('+r'$\pm$'+' axes) via anomaly-neutral non-negative additions;\n0 interior '+r'$\Rightarrow$ any residual, in particular $-\Theta_{\rm base}$, is cancellable.'
    ax.text(0.5,-0.19,note,transform=ax.transAxes,ha='center',va='top',fontsize=8.5,color='#444')
    for s in ('top','right'): ax.spines[s].set_visible(False)
    fig.tight_layout()
    tag='_es' if es else ''
    fig.savefig(OUT/f'fig_cone{tag}.pdf',bbox_inches='tight')
    fig.savefig(OUT/f'fig_cone{tag}.png',dpi=400,bbox_inches='tight'); plt.close(fig)

make(False); make(True)
print("wrote fig_cone.pdf/png and fig_cone_es.pdf/png")
