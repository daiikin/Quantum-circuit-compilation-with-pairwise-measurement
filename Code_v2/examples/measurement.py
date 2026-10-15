"""Optional command-line interface to the same run_pipeline function."""
import argparse
from pauli_compiler import run_pipeline


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument('--pauli')
    choice.add_argument('--weight', type=int)
    parser.add_argument('--seed', type=int)
    parser.add_argument('--total-qubits', type=int)
    parser.add_argument('--sign', type=int, choices=(-1,1))
    parser.add_argument('--outcome', type=int, choices=(-1,1), default=1)
    parser.add_argument('--no-simplify', action='store_true')
    parser.add_argument('--check-each', action='store_true')
    parser.add_argument('--no-resources', action='store_true')
    parser.add_argument('--no-plots', action='store_true')
    parser.add_argument('--show', action='store_true')
    parser.add_argument('--labels', action='store_true')
    parser.add_argument('--hide-checks', action='store_true')
    parser.add_argument('--record-legend', action='store_true')
    parser.add_argument('--output', default='outputs/cli')
    parser.add_argument('--format', choices=('svg','png'), default='svg')
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()
    if not args.show:
        import matplotlib
        matplotlib.use('Agg')
    result = run_pipeline(
        args.pauli, weight=args.weight, seed=args.seed, total_qubits=args.total_qubits,
        sign=args.sign, outcome=args.outcome, simplify_measurements=not args.no_simplify,
        check_each=args.check_each, calculate_resources=not args.no_resources,
        render=not args.no_plots, show=args.show, labels=args.labels,
        show_checks=not args.hide_checks, show_record_legend=args.record_legend,
        output_dir=args.output, image_format=args.format, overwrite=args.overwrite)
    print(result.decoder.formulas())
    print('\nSimplification:', result.statistics)
    print('Resources:', result.resources)
    print('Saved files:', result.files)
    result.close_figures()


if __name__ == '__main__':
    main()
