#чтение FASTA файлов в main

import argparse
from Bio import SeqIO
import io
import numpy as np


def get_seq_from_fasta(path):
    with open(path, "r") as f:
        return str(next(SeqIO.parse(f, "fasta")).seq)


def smith_waterman(seq1, seq2, match_score, mismatch_penalty, gap_penalty):
    return ("---", "---", 0)

#заполнение матрицы Нидлман
def needleman_wunsch(seq1, seq2, match_score, mismatch_penalty, gap_penalty):
    len1, len2 = len(seq1), len(seq2)
    #создаем матрицу
    matrix = np.zeros((len2 + 1, len1 + 1), dtype =int)
    #инициализируем 1 строку и 1 столбец
    for i in range(len2+1):
        matrix[i][0] = i * gap_penalty
    for i in range(len1+1):
        matrix[0][i] = i * gap_penalty

    #заполнение матрицы
    for i in range(1, len2 + 1):
        for j in range(1, len1 + 1): 
            if seq1[j-1] == seq2[i-1]:
                score_diag = matrix[i-1][j-1] + match_score
            else:
                score_diag = matrix[i-1][j-1] + mismatch_penalty
            
            score_up = matrix[i-1][j] + gap_penalty
            score_left = matrix[i][j-1] + gap_penalty

            matrix[i][j] = max(score_diag, score_up, score_left)
    return matrix   
    

def waterman_perform_traceback(self, seq1: str, seq2: str, score_matrix: list[list[int]], start_pos: tuple[int, int]) -> tuple[str, str]:
    return

def needleman_perform_traceback(self, seq1: str, seq2: str, score_matrix: list[list[int]], start_pos: tuple[int, int]) -> tuple[str, str]:
    return




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sequence alignment tool")
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--fasta", nargs=2, metavar=("FILE1", "FILE2"), help="Input FASTA files")
    group.add_argument("--seq", nargs=2, metavar=("SEQ1", "SEQ2"), help="Input raw sequenses")
    parser.add_argument("--method", choices=["nw", "sw"], required=True, help="Alignment method: nw (Needleman-Wunsch) or sw (Smith-Waterman)")
    parser.add_argument("--match", type=int, metavar=("match_score"), help="Input match score")
    parser.add_argument("--mismatch", type=int, metavar=("mismatch_penalty"), help="Input mismatch penalty")
    parser.add_argument("--gap", type=int, metavar=("gap_penalty"), help="Input gap penalty")
    parser.add_argument("--output", type=str, help="Output file to save result")

    args = parser.parse_args()

    #дефолтные значения
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
        result = needleman_wunsch(seq1, seq2, match_score, mismatch_penalty, gap_penalty)
    elif args.method == "sw":
        result = smith_waterman(seq1, seq2, match_score, mismatch_penalty, gap_penalty)

    # Вывод
    if args.output:
        with open(args.output, "w") as f:
            f.write(result)
        print(result) 
    else:
        print(result)    