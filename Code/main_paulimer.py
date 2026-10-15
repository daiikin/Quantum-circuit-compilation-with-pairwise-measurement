"""ZX -> MC flow -> pair measurements -> Paulimer outcome/feedback extraction.

Examples:
    python main_paulimer.py --weight 9 --seed 54
    python main_paulimer.py --pauli=-XYZXYZ
"""
import argparse
import json
from pathlib import Path

from zx_r4_r5 import generate_pauli_measurement, compile_measurement
from measurement_decoder import derive_decoder


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--pauli', help='Signed Pauli; for a minus sign use --pauli=-XYZ')
    group.add_argument('--weight', type=int, default=None)
    parser.add_argument('--seed', type=int, default=54)
    parser.add_argument('--output', type=Path, default=Path('paulimer_output'))
    parser.add_argument('--no-plots', action='store_true')
    args = parser.parse_args()
    measurement = args.pauli or generate_pauli_measurement(4 if args.weight is None else args.weight, seed=args.seed)
    compilation = compile_measurement(measurement)
    decoder = derive_decoder(compilation)
    folder = args.output
    folder.mkdir(parents=True, exist_ok=True)
    (folder/'measurement.stim').write_text(str(decoder.corrected_stim()), encoding='utf-8')
    (folder/'measurement.json').write_text(json.dumps(decoder.metadata(), indent=2), encoding='utf-8')
    (folder/'formulas.txt').write_text(decoder.formulas(), encoding='utf-8')
    print(decoder.formulas())
    print('\nBackend:', decoder.certificate['backend'])
    print('Output physical wires:', compilation.circuit.outputs)
    if not args.no_plots:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import pyzx as zx
        from flow_layout import draw_flow
        from circuit_renderer import draw_measurement_circuit
        m = compilation.measurement
        def save(graph, title, filename):
            fig = zx.draw_matplotlib(graph, labels=True)
            fig.suptitle(title)
            fig.savefig(folder/filename, bbox_inches='tight')
            plt.close(fig)
        save(m.original, 'Initial ZX measurement branch', '01_initial.svg')
        if m.rewritten is not None:
            save(m.rewritten, 'Optional direct R4/R5 preview', '02_direct_rewrite.svg')
        save(compilation.prepared.graph, 'MC-flow preparation', '03_prepared.svg')
        fig = draw_flow(compilation.reduced.graph, compilation.reduced.paths,
                        labels=False, title='Reduced ZX diagram with MC-flow layout')
        fig.savefig(folder/'04_reduced.svg', bbox_inches='tight')
        plt.close(fig)
        fig = draw_measurement_circuit(decoder, title='Paulimer-derived outcome and feedback',
                                       show_checks=True)
        fig.savefig(folder/'05_circuit.svg', bbox_inches='tight')
        plt.close(fig)
    print('\nSaved output to', folder.resolve())


if __name__ == '__main__':
    main()
