"""File export for a completed workflow; no simulation or compilation logic."""
import json
from pathlib import Path

_FIGURE_NAMES = {'initial':'01_initial', 'direct_rewrite':'02_direct_rewrite',
                 'prepared':'03_prepared', 'reduced':'04_reduced', 'circuit':'05_circuit'}
_TEXT_NAMES = ('raw_extracted.stim', 'before_feedback.stim', 'measurement.stim',
               'measurement.json', 'formulas.txt', 'resources.json', 'run.json')


def export_result(result, output_dir, *, image_format='svg', overwrite=False):
    """Write named outputs and return {artifact: absolute Path}.

    When overwriting, clear only this workflow's known artifact names so an
    optional preview or resources file from an older run cannot appear current.
    Unrelated files in the output directory are preserved.
    """
    if image_format not in ('svg', 'png'):
        raise ValueError('Unsupported image format')
    folder = Path(output_dir).expanduser().resolve()
    if folder.exists() and any(folder.iterdir()) and not overwrite:
        raise FileExistsError('Output directory is not empty')
    folder.mkdir(parents=True, exist_ok=True)
    if overwrite:
        names = list(_TEXT_NAMES) + [stem+'.'+ext for stem in _FIGURE_NAMES.values()
                                    for ext in ('svg','png')]
        for name in names:
            path = folder/name
            if path.is_file() or path.is_symlink():
                path.unlink()
    files = {}
    def write(key, name, text):
        path = folder/name
        path.write_text(text, encoding='utf-8')
        files[key] = path

    write('raw_circuit', 'raw_extracted.stim', str(result.raw_circuit.to_stim()))
    write('before_feedback', 'before_feedback.stim', str(result.simplified_circuit.to_stim()))
    write('corrected_circuit', 'measurement.stim', str(result.corrected_circuit))
    write('decoder', 'measurement.json', json.dumps(result.decoder.metadata(), indent=2))
    write('formulas', 'formulas.txt', result.decoder.formulas())
    if result.resources is not None:
        write('resources', 'resources.json', json.dumps(result.resources, indent=2))
    for key, figure in result.figures.items():
        path = folder/(_FIGURE_NAMES[key]+'.'+image_format)
        figure.savefig(path, bbox_inches='tight')
        files[key] = path
    manifest = {'parameters':result.parameters, 'statistics':result.statistics,
                'signed_pauli':result.decoder.certificate['signed_pauli'],
                'rewrite_trace':result.compilation.reduced.trace,
                'artifacts':{key:path.name for key,path in files.items()},
                'note':'Record indices refer to before_feedback.stim and measurement.stim, '
                       'not raw_extracted.stim when measurements were removed.'}
    write('run', 'run.json', json.dumps(manifest, indent=2))
    return files
