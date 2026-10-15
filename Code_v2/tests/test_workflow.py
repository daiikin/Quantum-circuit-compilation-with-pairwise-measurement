"""Public workflow integration and file/export parameter tests."""
import json
from pathlib import Path
import tempfile
import unittest
import matplotlib
matplotlib.use('Agg')
from pauli_compiler import run_pipeline
from tests.validate_decoder import check_shot


class WorkflowTests(unittest.TestCase):
    def test_random_seed_and_explicit_pauli_agree(self):
        generated = run_pipeline(weight=6, seed=54, total_qubits=8, sign=-1,
                                 outcome=-1, render=False)
        word = '-' + generated.measurement.pauli
        explicit = run_pipeline(word, outcome=-1, render=False)
        self.assertEqual(generated.corrected_circuit, explicit.corrected_circuit)
        self.assertEqual(len(generated.measurement.support),6)
        self.assertEqual(len(generated.measurement.pauli),8)
        self.assertEqual(generated.files,{})
        self.assertEqual(generated.figures,{})
        check_shot(generated.decoder,123,bell=True)

    def test_raw_and_simplified_paths_both_work(self):
        raw = run_pipeline('XXXX', simplify_measurements=False, render=False,
                           calculate_resources=False)
        simplified = run_pipeline('XXXX', render=False)
        self.assertEqual(raw.statistics['removed_measurements'],0)
        self.assertGreater(simplified.statistics['removed_measurements'],0)
        self.assertEqual(raw.raw_circuit, simplified.raw_circuit)
        self.assertIsNone(raw.resources)
        self.assertEqual(len(simplified.decoder.measurement_operations),
                         simplified.statistics['simplified_measurements'])
        for decoder in (raw.decoder,simplified.decoder):
            check_shot(decoder,84,bell=False,feedback=False)

    def test_render_export_and_overwrite_optional_artifacts(self):
        import matplotlib.pyplot as plt
        with tempfile.TemporaryDirectory() as directory:
            first = run_pipeline('XXXX',output_dir=directory)
            self.assertIn('direct_rewrite',first.figures)
            self.assertTrue(all(path.exists() for path in first.files.values()))
            self.assertEqual(first.files['corrected_circuit'].read_text(),str(first.corrected_circuit))
            manifest = json.loads(first.files['run'].read_text())
            self.assertEqual(manifest['statistics'],first.statistics)
            first.close_figures()
            self.assertEqual(plt.get_fignums(),[])
            sentinel = Path(directory)/'user_notes.txt';sentinel.write_text('preserve me')
            with self.assertRaises(FileExistsError):
                run_pipeline('Z',render=False,output_dir=directory)
            second = run_pipeline('-Y',render=False,calculate_resources=False,
                                  output_dir=directory,overwrite=True)
            self.assertEqual(second.figures,{})
            self.assertFalse((Path(directory)/'02_direct_rewrite.svg').exists())
            self.assertFalse((Path(directory)/'resources.json').exists())
            self.assertTrue(sentinel.exists())
            self.assertEqual(plt.get_fignums(),[])

    def test_high_weight_render_has_no_direct_preview(self):
        result = run_pipeline(weight=6,seed=54,render=True,show_checks=False)
        self.assertIsNone(result.measurement.rewritten)
        self.assertNotIn('direct_rewrite',result.figures)
        self.assertEqual(set(result.figures),{'initial','prepared','reduced','circuit'})
        result.close_figures()

    def test_reject_ambiguous_or_invalid_inputs(self):
        for kwargs in [dict(),dict(pauli='X',weight=1),dict(pauli='X',seed=4),
                       dict(pauli='X',sign=-1),dict(pauli='X',total_qubits=1),
                       dict(weight=0),dict(weight=3,total_qubits=2),
                       dict(weight=2,seed=True),dict(weight=2,sign=True),
                       dict(weight=2,outcome=0),dict(pauli='X',image_format='bad'),
                       dict(pauli='X',show=True,render=False)]:
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                run_pipeline(**kwargs)


if __name__ == '__main__':
    unittest.main()
