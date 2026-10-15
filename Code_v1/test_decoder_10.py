"""Ten circuit instances; symbolic all-record certification + sampled Stim checks."""
import json
from pathlib import Path
import random
import stim
from zx_r4_r5 import compile_measurement
from pauli_measurements import generate_pauli_measurement
from measurement_decoder import derive_decoder


def check_shot(decoder, seed, *, bell=True, feedback=True):
    m,c = decoder.compilation.measurement,decoder.compilation.circuit
    n,q = len(m.pauli),c.num_qubits
    actual = stim.TableauSimulator(seed=seed)
    ideal = stim.TableauSimulator(seed=seed)
    actual.set_num_qubits(q+n if bell else q)
    ideal.set_num_qubits(2*n if bell else n)
    if bell:
        for j,w in enumerate(c.inputs):
            actual.h(q+j)
            actual.cnot(q+j,w)
            ideal.h(n+j)
            ideal.cnot(n+j,j)
    else:
        # Independent random stabilizer inputs, including entanglement.
        rng = random.Random(seed)
        for _ in range(3*n):
            j = rng.randrange(n)
            gate = rng.choice(['H','S','CX']) if n>1 else rng.choice(['H','S'])
            targets = [j]
            if gate=='CX':
                targets += [rng.choice([i for i in range(n) if i!=j])]
            a,b = stim.Circuit(),stim.Circuit()
            a.append(gate,[c.inputs[t] for t in targets])
            b.append(gate,targets)
            actual.do_circuit(a)
            ideal.do_circuit(b)
    actual.do_circuit(decoder.corrected_stim() if feedback else c.to_stim())
    record = actual.current_measurement_record()
    bit = decoder.outcome_bit(record)
    if not feedback:
        frame = decoder.output_pauli_frame(record)
        correction = stim.Circuit()
        for j,w in enumerate(c.outputs):
            letter = {0:'I',1:'X',2:'Y',3:'Z'}[frame[j]]
            if letter!='I':
                correction.append(letter,[w])
        actual.do_circuit(correction)
    observable = stim.PauliString(m.pauli + ('I'*n if bell else ''))
    observable.sign = m.sign
    ideal.postselect_observable(observable,desired_value=bool(bit))
    mapping = tuple(c.outputs) + (tuple(q+j for j in range(n)) if bell else ())
    for generator in ideal.canonical_stabilizers():
        lifted = stim.PauliString(q+n if bell else q)
        lifted.sign = generator.sign
        for j,w in enumerate(mapping):
            lifted[w] = generator[j]
        assert actual.peek_observable_expectation(lifted)==1, 'Output state mismatch'
    return bit


def main():
    folder = Path('decoder_examples')
    folder.mkdir(exist_ok=True)
    results = []
    for i,weight in enumerate((1,2,3,4,5,6,7,9,13,20)):
        # Include an explicit XXXX example whose bit numbering is documented.
        if weight==4:
            from pauli_measurements import pauli_measurement
            m=pauli_measurement('XXXX')
        else:
            m=generate_pauli_measurement(weight,seed=20261009+i,
                                        total_qubits=weight+i%3,
                                        sign=-1 if i%3==0 else 1,
                                        outcome=-1 if i%2 else 1)
        d=derive_decoder(compile_measurement(m))
        outcomes=[]
        for shot in range(32):
            outcomes.append(check_shot(d,1000*i+shot,bell=True,feedback=shot%2==0))
        for shot in range(8):
            check_shot(d,50000+1000*i+shot,bell=False,feedback=shot%2==0)
        assert set(outcomes)=={0,1}, 'Both outcomes should be exercised'
        row={'weight':weight,'pauli':('-' if m.sign<0 else '+')+m.pauli,
             'measurements':len(d.measurement_operations),'bell_shots':32,'other_input_shots':8,
             'passed':True,'outcome_formula':str(d.outcome),'certificate':d.certificate}
        results.append(row)
        stem=f'{i+1:02d}_weight_{weight}'
        (folder/(stem+'.stim')).write_text(str(d.corrected_stim()))
        (folder/(stem+'.json')).write_text(json.dumps(d.metadata(),indent=2))
        (folder/(stem+'_formulas.txt')).write_text(d.formulas())
        print(json.dumps(row),flush=True)
    (folder/'validation.json').write_text(json.dumps(results,indent=2))
    print('10 symbolic instrument certificates and 400 sampled state checks passed.')


if __name__=='__main__':
    main()
