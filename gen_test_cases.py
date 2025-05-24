import os
import json
import random

BASES = {
    "dna": "ACGT",
    "rna": "ACGU",
    "protein": "ACDEFGHIKLMNPQRSTVWY"
}

def random_sequence(seq_type, length):
    return ''.join(random.choice(BASES[seq_type]) for _ in range(length))

if __name__ == "__main__":
    cases = []
    for i in range(1000):
        seq_type = random.choice(list(BASES.keys()))
        seq1 = random_sequence(seq_type, random.randint(1, 30))
        seq2 = random_sequence(seq_type, random.randint(1, 30))
        params = {
            "match_score": random.randint(1, 5),
            "mismatch_penalty": -random.randint(1, 3),
            "open_penalty": -random.randint(3, 5),
            "extend_penalty": -random.randint(1, 3),
        }
        test_case = {
            "seq1": seq1,
            "seq2": seq2,
            "seq_type": seq_type,
            "params": params
        }
        cases.append(test_case)
    with open("test_cases.json", "w") as f:
        json.dump(cases, f, indent=2)

