"""Edit these arguments and run: python main.py"""
from pauli_compiler import run_pipeline


def main():
    result = run_pipeline(
        weight=9,
        seed=54,
        # For a specific Pauli, replace weight and seed with: pauli='-XYZXYZ'
        simplify_measurements=True,
        check_each=False,
        calculate_resources=True,
        render=True,
        show=False,
        output_dir='outputs/demo',
        overwrite=True,
    )
    print(result.decoder.formulas())
    print('\nSimplification:', result.statistics)
    print('Resources:', result.resources)
    print('Saved files:', result.files)
    result.close_figures()


if __name__ == '__main__':
    main()
