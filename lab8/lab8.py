import numpy as np
import matplotlib.pyplot as plt

MATCH_SCORE = 1
MISMATCH_PENALTY = -1
GAP_PENALTY = -2


def needleman_wunsch(s1: str, s2: str, match: int = 1, mismatch: int = -1, gap: int = -2):
    n = len(s1)
    m = len(s2)

    score_matrix = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(n + 1):
        score_matrix[i][0] = i * gap
    for j in range(m + 1):
        score_matrix[0][j] = j * gap

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if s1[i - 1] == s2[j - 1]:
                diag = score_matrix[i - 1][j - 1] + match
            else:
                diag = score_matrix[i - 1][j - 1] + mismatch

            up = score_matrix[i - 1][j] + gap      
            left = score_matrix[i][j - 1] + gap    

            score_matrix[i][j] = max(diag, up, left)

    aligned_s1 = []
    aligned_s2 = []
    traceback_path = [(n, m)]

    i, j = n, m
    while i > 0 or j > 0:
        current = score_matrix[i][j]

        if i > 0 and j > 0:
            step_score = match if s1[i - 1] == s2[j - 1] else mismatch
            if current == score_matrix[i - 1][j - 1] + step_score:
                aligned_s1.append(s1[i - 1])
                aligned_s2.append(s2[j - 1])
                i -= 1
                j -= 1
                traceback_path.append((i, j))
                continue

        if i > 0 and current == score_matrix[i - 1][j] + gap:
            aligned_s1.append(s1[i - 1])
            aligned_s2.append("-")
            i -= 1
            traceback_path.append((i, j))
            continue

        if j > 0 and current == score_matrix[i][j - 1] + gap:
            aligned_s1.append("-")
            aligned_s2.append(s2[j - 1])
            j -= 1
            traceback_path.append((i, j))
            continue

    aligned_s1 = "".join(reversed(aligned_s1))
    aligned_s2 = "".join(reversed(aligned_s2))
    traceback_path.reverse()
    final_score = score_matrix[n][m]

    return score_matrix, aligned_s1, aligned_s2, final_score, traceback_path


def print_scoring_matrix(score_matrix: list, s1: str, s2: str):
    header = "      -  " + "  ".join(f"{ch:>2}" for ch in s2)
    print(header)
    print("-" * len(header))

    row_labels = ["-"] + list(s1)
    for r_idx, row in enumerate(score_matrix):
        row_str = "  ".join(f"{val:>2}" for val in row)
        print(f"{row_labels[r_idx]:>2} | {row_str}")


def plot_alignment_matrix(score_matrix: list, s1: str, s2: str, path: list):
    mat = np.array(score_matrix)
    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(mat, cmap="viridis")

    ax.set_xticks(range(len(s2) + 1))
    ax.set_yticks(range(len(s1) + 1))
    ax.set_xticklabels(["-"] + list(s2))
    ax.set_yticklabels(["-"] + list(s1))
    ax.xaxis.tick_top()

    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            color = "black" if (i, j) in path else "white"
            weight = "bold" if (i, j) in path else "normal"
            ax.text(j, i, str(mat[i, j]), ha="center", va="center", color=color, fontsize=8, fontweight=weight)

    path_y = [p[0] for p in path]
    path_x = [p[1] for p in path]
    ax.plot(path_x, path_y, color="#ff3b30", linewidth=2.5, marker="o", markersize=6, alpha=0.75, label="Traceback Optim")

    plt.colorbar(im, ax=ax, label="Scor aliniere")
    ax.set_title("Matricea Needleman-Wunsch & Traceback (Lab 8)", pad=25, fontweight="bold")
    ax.legend(loc="lower left")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    s1 = "ACCGTGAAGCCAATAC"
    s2 = "AGCGTGCAGCCAATAC"

    matrix, aln1, aln2, score, path = needleman_wunsch(
        s1, s2,
        match=MATCH_SCORE,
        mismatch=MISMATCH_PENALTY,
        gap=GAP_PENALTY
    )

    match_line = []
    matches = 0
    mismatches = 0
    gaps = 0

    for c1, c2 in zip(aln1, aln2):
        if c1 == c2:
            match_line.append("|")
            matches += 1
        elif c1 == "-" or c2 == "-":
            match_line.append(" ")
            gaps += 1
        else:
            match_line.append(".")
            mismatches += 1

    similarity = (matches / len(aln1)) * 100

    print("=== MATRICEA DE SCORURI NEEDLEMAN-WUNSCH ===\n")
    print_scoring_matrix(matrix, s1, s2)

    print("\n=== REZULTAT ALINIERE GLOBALĂ ===")
    print(f"S1:   {aln1}")
    print(f"      {''.join(match_line)}")
    print(f"S2:   {aln2}")
    print("-" * 40)
    print(f"Scor final:     {score}")
    print(f"Matches:        {matches} / {len(aln1)}")
    print(f"Mismatches:     {mismatches}")
    print(f"Gaps:           {gaps}")
    print(f"Similaritate:   {similarity:.2f}%")

    plot_alignment_matrix(matrix, s1, s2, path)