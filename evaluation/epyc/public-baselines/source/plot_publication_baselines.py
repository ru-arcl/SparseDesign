#!/usr/bin/env python3
"""Generate publication figures and a table from verified native baseline data."""
import argparse
import json
from pathlib import Path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--analysis',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--timeout',type=float,default=120)
    args = ap.parse_args()
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.transforms import ScaledTranslation
    summary = json.loads((args.analysis/'summary.json').read_text())
    records = json.loads((args.analysis/'results.json').read_text())
    if not summary['completion'].get('repeated',{}).get('complete'):
        raise ValueError('Refusing publication figures for an incomplete repeated grid')
    args.output.mkdir(parents=True,exist_ok=True)
    tools = [('sparsedesign','SparseDesign','#171717','o'),('lineardesign','LinearDesign','#3274a1','s'),
             ('linearcdsfold','LinearCDSfold','#e1812c','D'),('derna','DERNA','#3a923a','v')]
    plt.rcParams.update({'font.size':8,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    fig, axes = plt.subplots(2,2,figsize=(7.2,5.3),sharex=True,layout='constrained')
    for col, lam in enumerate((0,4)):
        timeout_conditions = {(r['tool'],r['aa_length']) for r in records
                              if r['phase'] in ('feasibility','repeated')
                              and r['lambda_value']==lam and r['status']=='timeout'}
        timeout_tools = {length:[t for t,_,_,_ in tools if (t,length) in timeout_conditions]
                         for _,length in timeout_conditions}
        near_cap_lengths = {s['aa_length'] for s in summary['groups']
                            if s['phase']=='repeated' and s['lambda_value']==lam
                            and s['successful']==5
                            and s.get('wall_seconds',{}).get('median',0)>=.8*args.timeout}
        for tool,label,color,marker in tools:
            groups = sorted((s for s in summary['groups'] if s['phase']=='repeated' and s['tool']==tool and s['lambda_value']==lam and s['successful']==5),key=lambda s:s['aa_length'])
            for row,metric,factor in [(0,'wall_seconds',1),(1,'max_rss_kib',1/2**20)]:
                good = [s for s in groups if metric in s]
                x = [3*s['aa_length'] for s in good]
                y = [factor*s[metric]['median'] for s in good]
                low = [factor*(s[metric]['median']-s[metric]['minimum']) for s in good]
                high = [factor*(s[metric]['maximum']-s[metric]['median']) for s in good]
                axes[row,col].errorbar(x,y,yerr=[low,high],label=label,color=color,marker=marker,markersize=4,linewidth=1,capsize=2)
            for length in sorted(n for t,n in timeout_conditions if t==tool):
                colocated = timeout_tools[length]
                offset = 10*(colocated.index(tool)-(len(colocated)-1)/2)
                if len(colocated)==1 and length in near_cap_lengths:
                    offset = -7
                # Shift only glyphs in typographic points. The stored data
                # coordinates retain the actual RNA length and wall cap.
                transform = axes[0,col].transData + ScaledTranslation(offset/72,0,fig.dpi_scale_trans)
                axes[0,col].scatter([3*length],[args.timeout],marker='^',facecolors='none',
                                    edgecolors=color,s=40,transform=transform)
        axes[0,col].set_title(r'$\lambda='+str(lam)+'$ (native objectives)',loc='left')
        axes[1,col].set_xlabel('RNA length (nt)')
        for row in range(2):
            axes[row,col].set_xscale('log');axes[row,col].set_yscale('log');axes[row,col].grid(axis='y',alpha=.18)
    axes[0,0].set_ylabel('Complete native CLI time (s)')
    axes[1,0].set_ylabel('Maximum resident set (GiB)')
    attempted_lengths = [3*r['aa_length'] for r in records if r['status'] in ('ok','timeout')]
    axes[0,0].set_xlim(min(attempted_lengths)*.85,max(attempted_lengths)*1.18)
    axes[0,0].legend(frameon=False,fontsize=7)
    fig.savefig(args.output/'native-cli-time-memory.pdf')
    fig.savefig(args.output/'native-cli-time-memory.png',dpi=180)
    plt.close(fig)
    tex = [r'\begin{longtable}{llrrrrr}',r'\caption{Native complete-CLI measurements, conditional on successful validated fresh repetitions. The column $r$ reports the successful count out of five attempted repetitions. Times are medians (minimum--maximum); RSS is the median peak resident set. Native models differ as specified in the protocol; these columns are not matched-model speed ratios.}\\',r'\toprule',r'Tool & Accession & AA & $\lambda$ & $r$ & Time (s) & RSS (GiB) \\',r'\midrule',r'\endfirsthead',r'\toprule Tool & Accession & AA & $\lambda$ & $r$ & Time (s) & RSS (GiB) \\ \midrule',r'\endhead']
    labels = {t:l for t,l,_,_ in tools}
    for s in sorted((s for s in summary['groups'] if s['phase']=='repeated'),key=lambda s:(s['aa_length'],s['lambda_value'],s['tool'])):
        time = s.get('wall_seconds',{})
        times = f"{time['median']:.2f} ({time['minimum']:.2f}--{time['maximum']:.2f})" if time else '--'
        rss = s.get('max_rss_kib',{}).get('median')
        memory = f'{rss/2**20:.3f}' if rss is not None else '--'
        tex.append(f"{labels[s['tool']]} & {s['accession']} & {s['aa_length']} & {s['lambda_value']:g} & {s['successful']}/5 & {times} & {memory}" + r' \\')
    tex += [r'\bottomrule',r'\end{longtable}','']
    (args.output/'native-cli-table.tex').write_text('\n'.join(tex))
    (args.output/'figure-caption.txt').write_text(f'Complete native CLI time and peak RSS on the prespecified final model-organism panel. Points show medians and whiskers show observed minima and maxima only for conditions with five successful validated repetitions after an excluded feasibility attempt. Hollow upward triangles mark {args.timeout:g}-second timeouts in feasibility or fresh repetitions. Overlapping timeout glyphs are shifted horizontally in display coordinates for visibility; their underlying RNA lengths and wall caps are unchanged. Conditions with fewer than five successes remain in the outcome table and are excluded from ordinary plotted points. Extension stops after the first feasibility timeout or allocation failure per tool and weight; unattempted lengths are not plotted. The 32-GiB cap applies to virtual address space and must not be interpreted as measured RSS. Native LinearDesign and LinearCDSfold impose additional structural restrictions; all lambda-4 objectives retain their native frequency arithmetic.\n')


if __name__ == '__main__':
    main()
