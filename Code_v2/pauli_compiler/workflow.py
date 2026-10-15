"""One-call workflow from a signed Pauli or random seed/weight to operational output.

No hidden global settings: all rendering, validation and export options are
explicit keyword arguments. See docs/PIPELINE_API.md for examples and semantics.
"""
from dataclasses import dataclass, field
from pathlib import Path

from .pipeline import Compilation, compile_measurement
from .zx.spiders import generate_pauli_measurement, pauli_measurement


@dataclass
class PipelineResult:
    """All stages of one run; use close_figures() when finished with plots."""
    compilation: Compilation
    decoder: object
    corrected_circuit: object
    resources: dict | None
    statistics: dict
    parameters: dict
    figures: dict = field(default_factory=dict)
    files: dict = field(default_factory=dict)

    @property
    def measurement(self):
        return self.compilation.measurement

    @property
    def raw_circuit(self):
        return self.compilation.raw_circuit

    @property
    def simplified_circuit(self):
        return self.compilation.circuit

    def close_figures(self):
        """Release only this run's Matplotlib figures; retain saved-file paths."""
        if self.figures:
            import matplotlib.pyplot as plt
            for figure in self.figures.values():
                plt.close(figure)
            self.figures.clear()


def run_pipeline(
    pauli=None, *, weight=None, seed=None, total_qubits=None, sign=None,
    outcome=1, simplify_measurements=True, check_each=False,
    calculate_resources=True, render=True, show=False, labels=False,
    show_checks=True, show_record_legend=False, output_dir=None,
    image_format='svg', overwrite=False, timing=None,
):
    """Run construction -> MC flow -> reduction -> extraction -> Paulimer.

    Supply EITHER a signed string pauli, OR weight with optional seed.
    Random generation accepts total_qubits >= weight and sign=+1/-1.
    outcome=+1/-1 chooses the initial ZX branch, never a physical outcome.
    Redundant measurements are removed before Paulimer by default.

    render=True creates and returns figures; show=True also calls plt.show().
    output_dir=None writes no files. Otherwise export circuits, formulas,
    metadata, resources (if requested), and figures (if rendered).
    A nonempty output directory requires overwrite=True. Only fixed names
    produced by this workflow may be overwritten/removed; other files survive.
    image_format is 'svg' or 'png'. The caller owns the returned figures.
    """
    if (pauli is None) == (weight is None):
        raise ValueError('Supply exactly one of pauli or weight')
    if outcome not in (-1, 1) or isinstance(outcome, bool):
        raise ValueError('outcome must be +1 or -1')
    if seed is not None and (not isinstance(seed, int) or isinstance(seed, bool)):
        raise ValueError('seed must be an integer or None')
    if sign is not None and (sign not in (-1, 1) or isinstance(sign, bool)):
        raise ValueError('sign must be +1 or -1')
    if image_format not in ('svg', 'png'):
        raise ValueError("image_format must be 'svg' or 'png'")
    if show and not render:
        raise ValueError('show=True requires render=True')
    if total_qubits is not None and (not isinstance(total_qubits, int) or isinstance(total_qubits, bool)):
        raise ValueError('total_qubits must be an integer or None')
    if pauli is not None:
        if seed is not None or total_qubits is not None or sign is not None:
            raise ValueError('seed, total_qubits and sign apply only to random weight generation; '
                             'encode signs/identities in the explicit Pauli string')
        measurement = pauli_measurement(pauli, outcome=outcome)
    else:
        measurement = generate_pauli_measurement(weight, seed=seed,
                                                total_qubits=total_qubits,
                                                sign=1 if sign is None else sign,
                                                outcome=outcome)
    destination = None if output_dir is None else Path(output_dir).expanduser().resolve()
    if destination is not None:
        if destination.exists() and not destination.is_dir():
            raise ValueError('output_dir must be a directory')
        if destination.exists() and any(destination.iterdir()) and not overwrite:
            raise FileExistsError('Output directory is not empty; choose a new directory or overwrite=True')

    compilation = compile_measurement(measurement, check_each=check_each,
                                      simplify_measurements=simplify_measurements)
    from .decoding.paulimer_backend import derive_decoder
    decoder = derive_decoder(compilation)
    resources = None
    if calculate_resources:
        from .analysis.resources import measure_circuit_resources
        resources = measure_circuit_resources(decoder, timing=timing)
    raw = compilation.raw_circuit
    optimized = compilation.circuit
    count = lambda c: sum(op.name in ('M', 'MX', 'MY', 'MXX', 'MYY', 'MZZ') for op in c.operations)
    statistics = {'raw_measurements': count(raw), 'simplified_measurements': count(optimized),
                  'removed_measurements': count(raw)-count(optimized),
                  'raw_operations': len(raw.operations), 'simplified_operations': len(optimized.operations)}
    parameters = {'input_pauli': pauli, 'weight': weight, 'seed': seed,
                  'total_qubits': total_qubits, 'sign': sign, 'outcome': outcome,
                  'simplify_measurements': simplify_measurements, 'check_each': check_each,
                  'calculate_resources': calculate_resources, 'render': render, 'show': show,
                  'labels': labels, 'show_checks': show_checks,
                  'show_record_legend': show_record_legend, 'image_format': image_format,
                  'output_dir': str(destination) if destination is not None else None,
                  'overwrite': overwrite,
                  'timing_model': resources['timing_model'] if resources is not None else None}
    result = PipelineResult(compilation, decoder, decoder.corrected_stim(),
                            resources, statistics, parameters)
    try:
        if render:
            from .visualization.stages import render_stages
            result.figures = render_stages(decoder, labels=labels, show_checks=show_checks,
                                           show_record_legend=show_record_legend)
        if destination is not None:
            from .export import export_result
            result.files = export_result(result, destination, image_format=image_format,
                                         overwrite=overwrite)
        if show:
            import matplotlib.pyplot as plt
            plt.show()
    except Exception:
        result.close_figures()
        raise
    return result
