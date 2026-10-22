"""Plot audited fixed-cycle R-factor comparisons with standard matplotlib."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',required=True);p.add_argument('--output',required=True)
    args=p.parse_args();data=json.loads(Path(args.input).read_text())
    native=[next(r for r in data['rows'] if r['label']==label) for label in ['PhaseBuild 0.6','PhaseBuild 0.7']]
    control=data['buccaneer_pipeline'];out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(11.2,5.4))
    for axis,key,title in zip(axes,['r_work','r_free'],['Rwork (%)','Rfree (%)']):
        for x,rows,label,color in [([10,20],native,'PhaseBuild: original-map rebuilding','#1764a0'),
                                   ([10*s['round'] for s in control],control,'Buccaneer: alternating updated maps','#bd5413')]:
            y=[100*r[key] for r in rows]
            axis.plot(x,y,'o-',label=label,color=color,lw=2,ms=6)
            for xx,yy in zip(x,y):
                offset=8 if color=='#1764a0' else -17
                left=color=='#bd5413' and xx==30
                axis.annotate('%.2f'%yy,(xx,yy),xytext=(-5 if left else 0,offset),textcoords='offset points',
                              ha='right' if left else 'center',color=color,fontsize=9)
        axis.set_title(title);axis.set_xlabel('Cumulative REFMAC cycles')
        axis.set_xticks([10,20,30,40,50]);axis.set_xlim(7,53);axis.grid(alpha=.2)
        axis.set_ylim((19,39) if key=='r_work' else (24,42))
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.91),ncol=2,frameon=False,fontsize=10)
    fig.suptitle('INS development-set comparison',fontsize=15,y=.98)
    fig.text(.08,.06,'Fixed last cycle; identical native observations, free flags and refinement settings.\n'
             'Different building/map-update workflows and budgets; Rfree is not an independent blind-test metric.',fontsize=9)
    fig.subplots_adjust(top=.78,bottom=.23,wspace=.25)
    fig.savefig(out/'r_factors.png',dpi=180)
    fig.savefig(out/'r_factors.svg')
    plt.close(fig)


if __name__=='__main__':
    main()
