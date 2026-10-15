"""Benchmark random Pauli measurements and plot depth/ancilla counts."""
import argparse
import json
from pathlib import Path
from pauli_compiler.analysis.resources import benchmark_resources


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--min-weight', type=int, default=4)
    parser.add_argument('--max-weight', type=int, default=10)
    parser.add_argument('--samples', type=int, default=10)
    parser.add_argument('--seed', type=int, default=20261015)
    parser.add_argument('--output', type=Path, default=Path('outputs/benchmark'))
    args = parser.parse_args()
    if args.min_weight < 1 or args.max_weight < args.min_weight:
        parser.error('Require 1 <= min-weight <= max-weight')
    rows, summary = benchmark_resources(weights=range(args.min_weight,args.max_weight+1),
                                         samples_per_weight=args.samples, seed=args.seed)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output/'circuits.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    weights = [row['weight'] for row in summary]
    fig, (depth, ancillas) = plt.subplots(1,2,figsize=(12,4.5))
    for key,label in [('quantum_depth','All operations'),('measurement_depth','Measurements only'),('depth_with_feedback','Aggregated X/Z feedback')]:
        line, = depth.plot(weights,[r[key]['mean'] for r in summary],'o-',label=label)
        depth.fill_between(weights,[r[key]['min'] for r in summary],
                           [r[key]['max'] for r in summary],color=line.get_color(),alpha=.15)
    depth.set_ylabel('Critical-path depth (stated unit costs)')
    depth.set_title('Fixed wire order; mean and min–max')
    depth.legend()
    ancillas.plot(weights,[r['ancilla_qubits']['mean'] for r in summary],'o-')
    ancillas.set_ylabel('Ancilla qubits')
    ancillas.set_title('Ancilla count without qubit reuse')
    for ax in (depth,ancillas):
        ax.set_xlabel('Pauli weight'); ax.set_xticks(weights); ax.grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(args.output/'resources.svg',bbox_inches='tight')
    plt.close(fig)
    print(f'Analyzed {len(rows)} circuits. Results: {args.output.resolve()}')


if __name__ == '__main__':
    main()
