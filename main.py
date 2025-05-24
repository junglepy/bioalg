# чтение FASTA файлов в main

import argparse
import io

import numpy as np
from Bio import SeqIO


def get_seq_from_fasta(path):
    with open(path, "r") as f:
        return str(next(SeqIO.parse(f, "fasta")).seq)

def needleman_wunsch(seq1, seq2, match_score, mismatch_penalty, gap_penalty):
    len1, len2 = len(seq1), len(seq2)
    matrix = np.zeros((len2 + 1, len1 + 1), dtype=int)
    for i in range(len2 + 1):
        matrix[i][0] = i * gap_penalty
    for i in range(len1 + 1):
        matrix[0][i] = i * gap_penalty

    for i in range(1, len2 + 1):
        for j in range(1, len1 + 1):
            if seq1[j - 1] == seq2[i - 1]:
                score_diag = matrix[i - 1][j - 1] + match_score
            else:
                score_diag = matrix[i - 1][j - 1] + mismatch_penalty

            score_up = matrix[i - 1][j] + gap_penalty
            score_left = matrix[i][j - 1] + gap_penalty

            matrix[i][j] = max(score_diag, score_up, score_left)
    align1, align2, score = needleman_perform_traceback(
        seq1, seq2, matrix, match_score, mismatch_penalty, gap_penalty
    )
    return align1, align2, score

def needleman_perform_traceback(
    seq1: str,
    seq2: str,
    score_matrix: np.ndarray,
    match_score: int,
    mismatch_penalty: int,
    gap_penalty: int,
) -> tuple[str, str, int]:

    i, j = len(seq1), len(seq2)
    aligned_seq1 = []
    aligned_seq2 = []
    final_score = score_matrix[i][j]

    while i > 0 or j > 0:
        symb_seq1 = seq1[i - 1] if i > 0 else "-"
        symb_seq2 = seq2[j - 1] if j > 0 else "-"
        char_score = match_score if symb_seq1 == symb_seq2 else mismatch_penalty

        if i > 0 and j > 0 and score_matrix[i][j] == score_matrix[i - 1][j - 1] + char_score:
            i -= 1
            j -= 1
            aligned_seq1.append(symb_seq1)
            aligned_seq2.append(symb_seq2)
        elif i > 0 and score_matrix[i][j] == score_matrix[i - 1][j] + gap_penalty:
            i -= 1
            aligned_seq1.append(symb_seq1)
            aligned_seq2.append("-")
        else:
            j -= 1
            aligned_seq1.append("-")
            aligned_seq2.append(symb_seq2)

    return "".join(aligned_seq1[::-1]), "".join(aligned_seq2[::-1]), final_score


def print_alignment(align1, align2, score):
    match_line = []
    for a, b in zip(align1, align2):
        if a == b:
            match_line.append('|')
        elif a == '-' or b == '-':
            match_line.append(' ')
        else:
            match_line.append('*')
    print(align1)
    print(''.join(match_line))
    print(align2)
    print(f"Alignment score: {score}")

def smith_waterman(seq1, seq2, match_score, mismatch_penalty, open_penalty, extend_penalty):
    len1, len2 = len(seq1), len(seq2)
    
    M = np.zeros((len2 + 1, len1 + 1), dtype=float)
    Ix = np.full((len2 + 1, len1 + 1), float('-inf'), dtype=float)
    Iy = np.full((len2 + 1, len1 + 1), float('-inf'), dtype=float)
    
    Ix[0][0] = Iy[0][0] = 0
    
    for i in range(1, len2 + 1):
        for j in range(1, len1 + 1):
            if seq1[j-1] == seq2[i-1]:
                score = match_score
            else:
                score = mismatch_penalty
            
            M[i][j] = max(
                M[i-1][j-1] + score,
                Ix[i-1][j-1] + score,
                Iy[i-1][j-1] + score,
                0
            )
            
            Ix[i][j] = max(
                M[i-1][j] + open_penalty,
                Ix[i-1][j] + extend_penalty
            )
            
            Iy[i][j] = max(
                M[i][j-1] + open_penalty,
                Iy[i][j-1] + extend_penalty
            )
    
    align1, align2, score = waterman_perform_traceback_affine(
        seq1, seq2, M, Ix, Iy, match_score, mismatch_penalty, open_penalty, extend_penalty
    )
    return align1, align2, int(score)

def waterman_perform_traceback_affine(
    seq1: str,
    seq2: str,
    M: np.ndarray,
    Ix: np.ndarray,
    Iy: np.ndarray,
    match_score: int,
    mismatch_penalty: int,
    open_penalty: int,
    extend_penalty: int,
) -> tuple[str, str, float]:
    
    score_matrix = np.maximum(np.maximum(M, Ix), Iy)
    max_i, max_j = np.unravel_index(np.argmax(score_matrix), score_matrix.shape)
    max_score = score_matrix[max_i, max_j]
    
    if max_score == M[max_i, max_j]:
        state = 'M'
    elif max_score == Ix[max_i, max_j]:
        state = 'Ix'
    else:
        state = 'Iy'
    
    i, j = max_i, max_j
    local_align1 = []
    local_align2 = []
    end_i, end_j = i, j
    
    while i > 0 and j > 0:
        current_max = max(M[i][j], Ix[i][j], Iy[i][j])
        
        if current_max <= 0:
            break
            
        if state == 'M':
            symb_seq1 = seq1[j - 1]
            symb_seq2 = seq2[i - 1]
            local_align1.append(symb_seq1)
            local_align2.append(symb_seq2)
            
            score = match_score if symb_seq1 == symb_seq2 else mismatch_penalty
            
            prev_M = M[i-1][j-1] + score if i > 0 and j > 0 else float('-inf')
            prev_Ix = Ix[i-1][j-1] + score if i > 0 and j > 0 else float('-inf')
            prev_Iy = Iy[i-1][j-1] + score if i > 0 and j > 0 else float('-inf')
            
            if abs(M[i][j] - prev_M) < 1e-9:
                state = 'M'
            elif abs(M[i][j] - prev_Ix) < 1e-9:
                state = 'Ix'
            elif abs(M[i][j] - prev_Iy) < 1e-9:
                state = 'Iy'
            else:
                state = 'M'
            
            i -= 1
            j -= 1
            
        elif state == 'Ix':
            local_align1.append('-')
            local_align2.append(seq2[i-1])
            
            prev_M = M[i-1][j] + open_penalty if i > 0 else float('-inf')
            prev_Ix = Ix[i-1][j] + extend_penalty if i > 0 else float('-inf')
            
            if abs(Ix[i][j] - prev_M) < 1e-9:
                state = 'M'
            else:
                state = 'Ix'
            
            i -= 1
            
        elif state == 'Iy':
            local_align1.append(seq1[j-1])
            local_align2.append('-')
            
            prev_M = M[i][j-1] + open_penalty if j > 0 else float('-inf')
            prev_Iy = Iy[i][j-1] + extend_penalty if j > 0 else float('-inf')
            
            if abs(Iy[i][j] - prev_M) < 1e-9:
                state = 'M'
            else:
                state = 'Iy'
            
            j -= 1
    
    start_i, start_j = i, j
    
    local_align1 = local_align1[::-1]
    local_align2 = local_align2[::-1]
    
    seq1_aligned_start = start_j
    seq1_aligned_end = end_j
    seq2_aligned_start = start_i
    seq2_aligned_end = end_i
    
    full_align1_parts = []
    full_align2_parts = []
    
    seq1_prefix = seq1[:seq1_aligned_start]
    seq2_prefix = seq2[:seq2_aligned_start]
    
    seq1_suffix = seq1[seq1_aligned_end:]
    seq2_suffix = seq2[seq2_aligned_end:]
    
    max_prefix_len = max(len(seq1_prefix), len(seq2_prefix))
    if max_prefix_len > 0:
        seq1_prefix_padded = '-' * (max_prefix_len - len(seq1_prefix)) + seq1_prefix
        seq2_prefix_padded = '-' * (max_prefix_len - len(seq2_prefix)) + seq2_prefix
        full_align1_parts.append(seq1_prefix_padded)
        full_align2_parts.append(seq2_prefix_padded)
    
    full_align1_parts.append(''.join(local_align1))
    full_align2_parts.append(''.join(local_align2))
    
    max_suffix_len = max(len(seq1_suffix), len(seq2_suffix))
    if max_suffix_len > 0:
        seq1_suffix_padded = seq1_suffix + '-' * (max_suffix_len - len(seq1_suffix))
        seq2_suffix_padded = seq2_suffix + '-' * (max_suffix_len - len(seq2_suffix))
        full_align1_parts.append(seq1_suffix_padded)
        full_align2_parts.append(seq2_suffix_padded)
    
    full_align1 = ''.join(full_align1_parts)
    full_align2 = ''.join(full_align2_parts)
    
    return full_align1, full_align2, max_score

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sequence alignment tool")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--fasta", nargs=2, metavar=("FILE1", "FILE2"), help="Input FASTA files"
    )
    group.add_argument(
        "--seq", nargs=2, metavar=("SEQ1", "SEQ2"), help="Input raw sequenses"
    )
    parser.add_argument(
        "--method",
        choices=["nw", "sw"],
        required=True,
        help="Alignment method: nw (Needleman-Wunsch) or sw (Smith-Waterman)",
    )
    parser.add_argument(
        "--match", type=int, metavar=("match_score"), help="Input match score"
    )
    parser.add_argument(
        "--mismatch",
        type=int,
        metavar=("mismatch_penalty"),
        help="Input mismatch penalty",
    )
    parser.add_argument(
        "--gap", type=int, metavar=("gap_penalty"), help="Input gap penalty"
    )
    parser.add_argument("--output", type=str, help="Output file to save result")

    args = parser.parse_args()


    match_score = 1
    mismatch_penalty = -1
    gap_penalty = -2

    if args.fasta:
        seq1 = get_seq_from_fasta(args.fasta[0])
        seq2 = get_seq_from_fasta(args.fasta[1])
    elif args.seq:
        seq1, seq2 = args.seq

    if args.match:
        match_score = args.match

    if args.mismatch:
        mismatch_penalty = args.mismatch

    if args.gap:
        gap_penalty = args.gap


    if args.method == "nw":
        align1, align2, score = needleman_wunsch(seq1, seq2, match_score, mismatch_penalty, gap_penalty)
    elif args.method == "sw":
        align1, align2, score = smith_waterman(seq1, seq2, match_score, mismatch_penalty, gap_penalty)

    if args.output:
        with open(args.output, "w") as f:
            f.write(f"{align1}\n{align2}\nScore: {score}\n")
        print_alignment(align1, align2, score)
    else:
        print_alignment(align1, align2, score)

