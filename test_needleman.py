import json
from Bio import pairwise2
from main import needleman_wunsch


test_cases = []
with open('test_cases.json', 'r') as f:
    test_cases = json.load(f)

full = 0
only_score = 0

for tc in test_cases:


    actual1, actual2, actual_score = needleman_wunsch(
        seq1=tc['seq1'],
        seq2=tc['seq2'],
        match_score=tc['params']['match_score'],
        mismatch_penalty=tc['params']['mismatch_penalty'],
        open_penalty=tc['params']['open_penalty'],
        extend_penalty = tc['params']['extend_penalty'],
    )
    is_full = 0
    is_only_score = 0
    for expected in pairwise2.align.globalms(
        tc['seq1'],
        tc['seq2'],
        tc['params']['match_score'],
        tc['params']['mismatch_penalty'],
        tc['params']['open_penalty'],
        tc['params']['extend_penalty'],
    ):
        if (actual1 == expected.seqA and actual2 == expected.seqB) and actual_score == expected.score:
            # full += 1
            is_full = 1
            continue
        if actual_score == expected.score:
            # only_score += 1
            is_only_score = 1
            continue
        
        


        print(f"Actual: {actual1=}, {expected.seqA=}, {actual_score=}")
        print(f"Actual: {actual2=}, {expected.seqB=}, {actual_score=}")
        print(f'Given: {tc["seq1"]=}, {tc["seq2"]=}')
        assert False, "Scores do not match"

    if is_full:
        is_only_score = 0
    full += is_full
    only_score += is_only_score


print(f"Full matches: {full}, Only score matches: {only_score}")
