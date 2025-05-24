import json
from Bio import pairwise2
from main import needleman_wunsch


test_cases = []
with open('test_cases.json', 'r') as f:
    test_cases = json.load(f)

for tc in test_cases:
    expected = pairwise2.align.globalms(
        tc['seq1'],
        tc['seq2'],
        tc['params']['match_score'],
        tc['params']['mismatch_penalty'],
        tc['params']['open_penalty'],
        tc['params']['extend_penalty'],
    )[0]
    actual1, actual2, actual_score = needleman_wunsch(
        seq1=tc['seq1'],
        seq2=tc['seq2'],
        match_score=tc['params']['match_score'],
        mismatch_penalty=tc['params']['mismatch_penalty'],
        gap_penalty=tc['params']['open_penalty'],
        #TODO: tc['params']['extend_penalty'],
    )

    assert actual1 == expected.seqA
    assert actual2 == expected.seqB
    assert actual_score == expected.score


