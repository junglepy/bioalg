# чтение FASTA файлов в main

import argparse
import io

import numpy as np
from Bio import SeqIO


def get_seq_from_fasta(path):
    with open(path, "r") as f:
        return str(next(SeqIO.parse(f, "fasta")).seq)

def needleman_wunsch(seq1, seq2, match_score, mismatch_penalty, gap_penalty, gap_extend_penalty):
    len1, len2 = len(seq1), len(seq2)

    # Основная матрица очков
    score = np.zeros((len2 + 1, len1 + 1), dtype=int)
    # Матрицы для хранения информации о продолжении гэпа
    gapA = np.zeros((len2 + 1, len1 + 1), dtype=bool)
    gapB = np.zeros((len2 + 1, len1 + 1), dtype=bool)

    # Инициализация первой строки и столбца
    for i in range(1, len2 + 1):
        score[i][0] = gap_penalty + (i - 1) * gap_extend_penalty
        gapA[i][0] = True  # гэп в первой последовательности

    for j in range(1, len1 + 1):
        score[0][j] = gap_penalty + (j - 1) * gap_extend_penalty
        gapB[0][j] = True  # гэп во второй последовательности

    # Заполнение матрицы
    for i in range(1, len2 + 1):
        for j in range(1, len1 + 1):
            if seq1[j - 1] == seq2[i - 1]:
                match = match_score
            else:
                match = mismatch_penalty

            score_diag = score[i - 1][j - 1] + match

            # Гэп в seq1 (вертикальный)
            if gapA[i - 1][j]:
                score_up = score[i - 1][j] + gap_extend_penalty
            else:
                score_up = score[i - 1][j] + gap_penalty
            # Гэп в seq2 (горизонтальный)
            if gapB[i][j - 1]:
                score_left = score[i][j - 1] + gap_extend_penalty
            else:
                score_left = score[i][j - 1] + gap_penalty

            best = max(score_diag, score_up, score_left)
            score[i][j] = best
            gapA[i][j] = best == score_up
            gapB[i][j] = best == score_left

    align1, align2, final_score = needleman_perform_traceback_affine(
        seq1, seq2, score, match_score, mismatch_penalty, gap_penalty, gap_extend_penalty
    )
    return align1, align2, final_score

def smith_waterman(seq1, seq2, match_score, mismatch_penalty, gap_penalty):
    len1, len2 = len(seq1), len(seq2)
    matrix = np.zeros((len2 + 1, len1 + 1), dtype=int)
    
    for i in range(1, len2 + 1):
        for j in range(1, len1 + 1):
            if seq1[j-1] == seq2[i-1]:
                diag = matrix[i-1][j-1] + match_score
            else:
                diag = matrix[i-1][j-1] + mismatch_penalty
            
            up = matrix[i-1][j] + gap_penalty
            left = matrix[i][j-1] + gap_penalty
            
            matrix[i][j] = max(diag, up, left, 0)
    align1, align2, score = waterman_perform_traceback(
        seq1, seq2, matrix, match_score, mismatch_penalty, gap_penalty
    )
    return align1, align2, score


def waterman_perform_traceback(
    seq1: str,
    seq2: str,
    score_matrix: np.ndarray,
    match_score: int,
    mismatch_penalty: int,
    gap_penalty: int,
) -> tuple[str, str, int]:

    i, j = np.unravel_index(np.argmax(score_matrix), score_matrix.shape)
    max_score = score_matrix[i][j]
    aligned_seq1 = []
    aligned_seq2 = []

    while score_matrix[i][j] != 0:
        curr_score = score_matrix[i][j]
        symb_seq1 = seq1[i - 1]
        symb_seq2 = seq2[j - 1]
        char_score = match_score if symb_seq1 == symb_seq2 else mismatch_penalty
        if curr_score == score_matrix[i - 1][j - 1] + char_score:
            i -= 1
            j -= 1
            aligned_seq1.append(symb_seq1)
            aligned_seq2.append(symb_seq2)
        elif curr_score == score_matrix[i - 1][j] + gap_penalty:
            i -= 1
            aligned_seq1.append(symb_seq1)
            aligned_seq2.append("-")
        elif curr_score == score_matrix[i][j - 1] + gap_penalty:
            j -= 1
            aligned_seq1.append("-")
            aligned_seq2.append(symb_seq2)
        else:
            i -= 1
            j -= 1
            aligned_seq1.append(symb_seq1)
            aligned_seq2.append(symb_seq2)

    return "".join(aligned_seq1[::-1]), "".join(aligned_seq2[::-1]), max_score


def needleman_perform_traceback_affine(seq1, seq2, score_matrix, match_score, mismatch_penalty, gap_open_penalty, gap_extend_penalty):
    i, j = len(seq2), len(seq1)
    aligned_seq1 = []
    aligned_seq2 = []

    while i > 0 or j > 0:
        current_score = score_matrix[i][j]
        if i > 0 and j > 0:
            char1 = seq1[j - 1]
            char2 = seq2[i - 1]
            match = match_score if char1 == char2 else mismatch_penalty
            if current_score == score_matrix[i - 1][j - 1] + match:
                aligned_seq1.append(char1)
                aligned_seq2.append(char2)
                i -= 1
                j -= 1
                continue

        if i > 0:
            # Гэп в seq1
            gap = gap_extend_penalty if i >= 2 and score_matrix[i][j] == score_matrix[i - 1][j] + gap_extend_penalty else gap_open_penalty
            if current_score == score_matrix[i - 1][j] + gap:
                aligned_seq1.append("-")
                aligned_seq2.append(seq2[i - 1])
                i -= 1
                continue

        if j > 0:
            # Гэп в seq2
            gap = gap_extend_penalty if j >= 2 and score_matrix[i][j] == score_matrix[i][j - 1] + gap_extend_penalty else gap_open_penalty
            if current_score == score_matrix[i][j - 1] + gap:
                aligned_seq1.append(seq1[j - 1])
                aligned_seq2.append("-")
                j -= 1
                continue

        break

    return "".join(aligned_seq1[::-1]), "".join(aligned_seq2[::-1]), score_matrix[len(seq2)][len(seq1)]


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

    # дефолтные значения
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

    # Запуск нужного алгоритма
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

