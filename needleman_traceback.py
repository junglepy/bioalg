def needleman_perform_traceback(
    seq1: str,
    seq2: str,
    score_matrix: list[list[int]],
    start_pos: tuple[int, int],
) -> tuple[str, str]:

    i, j = start_pos
    aligned_seq1 = []
    aligned_seq2 = []

    while i > 0 or j > 0:
        symb_seq1 = seq1[i - 1] if i > 0 else "-"
        symb_seq2 = seq2[j - 1] if j > 0 else "-"
        
        if i > 0 and j > 0 and score_matrix[i][j] == score_matrix[i - 1][j - 1] + (
            match_score if symb_seq1 == symb_seq2 else mismatch_penalty
        ):
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

    return "".join(aligned_seq1[::-1]), "".join(aligned_seq2[::-1])


