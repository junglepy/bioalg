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


