"""Small value type for displaying and evaluating XOR expressions.
Bit 0 of mask is the constant; bit i+1 is physical record bit b_i."""
from dataclasses import dataclass

@dataclass(frozen=True)
class Affine:
    mask: int = 0

    def __xor__(self, other):
        return Affine(self.mask ^ other.mask)

    @property
    def constant(self):
        return self.mask & 1

    @property
    def indices(self):
        return tuple(i for i in range(self.mask.bit_length()-1)
                     if (self.mask >> (i+1)) & 1)

    def evaluate(self, record):
        value = self.constant
        for i in self.indices:
            value ^= int(record[i])
        return value

    def __str__(self):
        return ' xor '.join((['1'] if self.constant else []) +
                            [f'b{i}' for i in self.indices]) or '0'

    def as_dict(self):
        return {'constant': self.constant, 'record_indices': self.indices,
                'formula': str(self)}
