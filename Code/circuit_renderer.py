"""Matplotlib renderer for joint measurements and decoded Pauli feedback.

    fig = draw_measurement_circuit(decoder)
    fig.savefig('measurement.svg', bbox_inches='tight')

Uses only decoder data. No circuit mutation or extra quantum optimization.
Connected X/X or Z/Z boxes denote ONE joint measurement, not two singles.
Only the boxes' two wires participate; crossing intervening wires does not.
Dashed boxes show classically controlled OUTPUT corrections.
"""
import textwrap


def _expr(expression):
    return str(expression).replace(' xor ', ' ⊕ ')


def _tex(expression):
    terms = (['1'] if expression.constant else [])
    terms += [f'b_{{{i}}}' for i in expression.indices]
    return r'\oplus '.join(terms) or '0'


def _schedule(circuit, row_of):
    """ASAP on each wire; reserve vertical gate spans to keep drawing clear."""
    next_time = [0]*circuit.num_qubits
    occupied = {}
    result = []
    measurement_index = 0
    for operation in circuit.operations:
        rows = [row_of[q] for q in operation.targets]
        lo,hi = min(rows),max(rows)
        t = max(next_time[q] for q in operation.targets)
        while any(not (hi<a or lo>b) for a,b in occupied.get(t,())):
            t += 1
        occupied.setdefault(t,[]).append((lo,hi))
        for q in operation.targets:
            next_time[q] = t+1
        bit = None
        if operation.name in ('M','MX','MXX','MZZ'):
            bit = measurement_index
            measurement_index += 1
        result.append((operation,t,bit))
    return result


def draw_measurement_circuit(decoder, *, title=None, wire_order=None,
                             show_checks=True, show_record_legend=False):
    """Return a matplotlib Figure of the operational measurement circuit.

    Physical wire IDs p0,... are distinct from logical input/output labels D0,...
    wire_order optionally specifies a permutation of all physical wire IDs.
    Feedback expressions longer than three terms get names defined below.
    For large circuits SVG is recommended; the canvas expands with circuit size.
    Circuit order is preserved on every wire; independent gates may be aligned.
    Result indices b_i retain the ORIGINAL decoder record ordering.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch, Rectangle
    from matplotlib.lines import Line2D
    c = decoder.compilation.circuit
    qcount = c.num_qubits
    order = tuple(range(qcount)) if wire_order is None else tuple(wire_order)
    if len(order)!=qcount or set(order)!=set(range(qcount)):
        raise ValueError('wire_order must be a permutation of physical wire IDs')
    row_of = {q:r for r,q in enumerate(order)}
    scheduled = _schedule(c,row_of)
    width = max((t for _,t,_ in scheduled),default=0)+1
    # One grid column is 1.7 drawing units; output corrections are wider.
    x_of = lambda t: 1.7*t
    correction_x = 1.7*width+.8
    end_x = correction_x+6.0
    inputs = {q:j for j,q in enumerate(c.inputs)}
    outputs = {q:j for j,q in enumerate(c.outputs)}
    first,last = {},{}
    for op,t,_ in scheduled:
        for q in op.targets:
            first.setdefault(q,x_of(t))
            last[q]=x_of(t)
    footer = [f'Measured observable: {decoder.certificate["signed_pauli"]}',
              f'Outcome bit: r = {_expr(decoder.outcome)}    Eigenvalue: (-1)^r',
              'Connected X/X or Z/Z boxes = one joint measurement; crossing wires are not targets.',
              'Dashed boxes = classical feedback. Apply it OR track the Pauli frame, not both.']
    if show_checks:
        footer += [f'Noiseless check: {_expr(e)} = 0' for e in decoder.checks]
    correction_labels = {}
    for j,(x,z) in enumerate(zip(decoder.x_corrections,decoder.z_corrections)):
        for letter,e in (('X',x),('Z',z)):
            if e.mask:
                if len(e.indices)+e.constant>3:
                    symbol = ('f' if letter=='X' else 'g')+'_{'+str(j)+'}'
                    correction_labels[j,letter] = f'${letter}^{{{symbol}}}$'
                    footer.append(f'{"f" if letter=="X" else "g"}{j} = {_expr(e)}')
                else:
                    correction_labels[j,letter] = f'${letter}^{{{_tex(e)}}}$'
    if show_record_legend:
        footer += [f'b{e["index"]}: {e["name"]} on physical wires {e["physical_wires"]}'
                   for e in decoder.measurement_operations]
    footer = [part for line in footer for part in textwrap.wrap(line,width=120)]
    plot_height = max(3.5,.60*qcount+1.6)
    footer_height = .23*len(footer)+.35
    fig = plt.figure(figsize=(max(12,(end_x+5)*.43),plot_height+footer_height))
    gs=fig.add_gridspec(2,1,height_ratios=[plot_height,footer_height],hspace=.04)
    ax=fig.add_subplot(gs[0])
    notes=fig.add_subplot(gs[1]); notes.axis('off')
    notes.text(0,1,'\n'.join(footer),va='top',ha='left',fontsize=9,
               linespacing=1.45,transform=notes.transAxes)
    for row,q in enumerate(order):
        y=-row
        start=-1 if q in inputs else first.get(q,0)-.3
        end=end_x if q in outputs else last.get(q,0)+.35
        ax.plot([start,end],[y,y],color='#64748b',lw=1,zorder=1)
        left=f'p{q}  '+(f'D{inputs[q]} in' if q in inputs else 'ancilla')
        ax.text(-1.3,y,left,ha='right',va='center',fontsize=9)
        if q in outputs:
            ax.annotate('',xy=(end_x+.2,y),xytext=(end_x-.3,y),
                        arrowprops={'arrowstyle':'->','color':'#334155'})
            ax.text(end_x+.35,y,f'D{outputs[q]} out',ha='left',va='center',fontsize=9)
        else:
            ax.text(last.get(q,0)+.6,y,'discard',ha='left',va='center',fontsize=7,color='#64748b')

    def box(x,y,text,*,color='#fff3c4',wide=.70,dashed=False):
        patch=FancyBboxPatch((x-wide/2,y-.23),wide,.46,
                            boxstyle='round,pad=0.03,rounding_size=0.04',
                            facecolor=color,edgecolor='#334155',lw=1.1,
                            linestyle='--' if dashed else '-',zorder=4)
        ax.add_patch(patch)
        ax.text(x,y,text,ha='center',va='center',fontsize=9,zorder=5)

    for op,t,bit in scheduled:
        x=x_of(t); ys=[-row_of[q] for q in op.targets]
        if op.name in ('MXX','MZZ','CX'):
            # White under-stroke distinguishes mere wire crossings.
            ax.plot([x,x],[min(ys),max(ys)],color='white',lw=4,zorder=2)
            ax.plot([x,x],[min(ys),max(ys)],color='#2563a6',lw=1.5,zorder=3)
        if op.name in ('MXX','MZZ'):
            letter='X' if op.name=='MXX' else 'Z'
            for y in ys:
                box(x,y,letter,color='#deedff')
            ax.text(x,max(ys)+.31,f'$b_{{{bit}}}=M_{{{letter}{letter}}}$',
                    ha='center',va='bottom',fontsize=8,color='#17436d',zorder=5)
        elif op.name in ('M','MX'):
            letter='Z' if op.name=='M' else 'X'
            box(x,ys[0],f'$M_{letter}$',color='#deedff')
            ax.text(x,ys[0]+.31,f'$b_{{{bit}}}$',ha='center',va='bottom',fontsize=8)
        elif op.name in ('R','RX'):
            box(x,ys[0],'$|0\\rangle$' if op.name=='R' else '$|+\\rangle$',color='#edf1f5')
        elif op.name=='CX':
            ax.scatter([x],[ys[0]],s=22,c='#334155',zorder=5)
            box(x,ys[1],'$\\oplus$',color='white')
        else:
            box(x,ys[0],r'$S^\dagger$' if op.name=='S_DAG' else op.name)
    for j,q in enumerate(c.outputs):
        for offset,letter in ((0,'X'),(2.9,'Z')):
            if (j,letter) in correction_labels:
                box(correction_x+offset,-row_of[q],correction_labels[j,letter],
                    color='#e2f4e7',wide=2.5,dashed=True)
    ax.set_title(title or 'Pauli measurement circuit with classical feedback',fontsize=14,pad=22)
    ax.set_xlim(-4.5,end_x+2.4)
    ax.set_ylim(-qcount+.25,.95)
    ax.axis('off')
    fig.subplots_adjust(left=.04,right=.98,top=.94,bottom=.025)
    return fig
