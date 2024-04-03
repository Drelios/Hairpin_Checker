import sys

def nussinov_jacobson_dna(sequence):
    n = len(sequence)
    matrix = [[0] * n for _ in range(n)]

    for length in range(1, n):
        for i in range(n - length):
            j = i + length

            # Case 1: Pairing between i and j
            if (sequence[i] == 'A' and sequence[j] == 'T') or (sequence[i] == 'T' and sequence[j] == 'A') or \
               (sequence[i] == 'G' and sequence[j] == 'C') or (sequence[i] == 'C' and sequence[j] == 'G'):
                matrix[i][j] = matrix[i + 1][j - 1] + 1

            # Case 2: Pairing between i and k, and between k+1 and j (where i <= k < j)
            for k in range(i, j):
                matrix[i][j] = max(matrix[i][j], matrix[i][k] + matrix[k + 1][j])

    return matrix, sequence

def nussinov_jacobson_rna(sequence):
    n = len(sequence)
    matrix = [[0] * n for _ in range(n)]

    for length in range(1, n):
        for i in range(n - length):
            j = i + length

            # Case 1: Pairing between i and j
            if (sequence[i] == 'A' and sequence[j] == 'U') or (sequence[i] == 'U' and sequence[j] == 'A') or \
               (sequence[i] == 'G' and sequence[j] == 'C') or (sequence[i] == 'C' and sequence[j] == 'G'):
                matrix[i][j] = matrix[i + 1][j - 1] + 1

            # Case 2: Pairing between i and k, and between k+1 and j (where i <= k < j)
            for k in range(i, j):
                matrix[i][j] = max(matrix[i][j], matrix[i][k] + matrix[k + 1][j])

    return matrix, sequence

def get_dot_bracket_notation_dna(matrix, sequence):
    n = len(sequence)
    dot_bracket = ['.'] * n

    def traceback(i, j):
        if i >= j:
            return

        if matrix[i][j] == matrix[i + 1][j - 1] + 1 and \
                ((sequence[i] == 'A' and sequence[j] == 'T') or (sequence[i] == 'T' and sequence[j] == 'A') or
                 (sequence[i] == 'G' and sequence[j] == 'C') or (sequence[i] == 'C' and sequence[j] == 'G')):
            dot_bracket[i] = '('
            dot_bracket[j] = ')'
            traceback(i + 1, j - 1)
        else:
            for k in range(i, j):
                if matrix[i][j] == matrix[i][k] + matrix[k + 1][j]:
                    traceback(i, k)
                    traceback(k + 1, j)
                    break

    traceback(0, n - 1)
    return ''.join(dot_bracket)

def get_dot_bracket_notation_rna(matrix, sequence):
    n = len(sequence)
    dot_bracket = ['.'] * n

    def traceback(i, j):
        if i >= j:
            return

        if matrix[i][j] == matrix[i + 1][j - 1] + 1 and \
                ((sequence[i] == 'A' and sequence[j] == 'U') or (sequence[i] == 'U' and sequence[j] == 'A') or
                 (sequence[i] == 'G' and sequence[j] == 'C') or (sequence[i] == 'C' and sequence[j] == 'G')):
            dot_bracket[i] = '('
            dot_bracket[j] = ')'
            traceback(i + 1, j - 1)
        else:
            for k in range(i, j):
                if matrix[i][j] == matrix[i][k] + matrix[k + 1][j]:
                    traceback(i, k)
                    traceback(k + 1, j)
                    break

    traceback(0, n - 1)
    return ''.join(dot_bracket)

def check_dot_bracket(string, threshold, gap_authorized):
    count = 0
    gap_count = 0
    for char in string:
        if char == '(':
            count += 1
            if count > threshold and gap_count <= gap_authorized:
                return "True"
            gap_count = 0
        elif char == ')':
            count = 0
            gap_count = 0
        else:
            if count > 0 and count <= threshold:
                gap_count += 1
                if gap_count > gap_authorized:
                    count = 0
                    gap_count = 0
            else:
                count = 0
                gap_count = 0
    return "False"


threshold=5                             #Number of consecutives paired bases
gap_authorized=2                        #Number of non-paired bases between 2 paires bases

type=str(input("Is your sequence DNA or RNA?: ")).upper()
sequence=str(input("Enter your sequence: ")).upper()


if type == 'DNA':
    if all(base in 'ATGC' for base in sequence):
        algo = nussinov_jacobson_dna
        dot_bracket = get_dot_bracket_notation_dna
    else:
        print("Error in the provided sequence")
        sys.exit()
elif type == 'RNA':
    if all(base in 'AUGC' for base in sequence):
        algo = nussinov_jacobson_rna
        dot_bracket = get_dot_bracket_notation_rna
    else:
        print("Error in the provided sequence")
        sys.exit()
else:
    print("Error in the sequence type")
    sys.exit()


result_matrix, result_sequence = algo(sequence)
print(algo(sequence))
hairpin_structure = dot_bracket(result_matrix, result_sequence)
result_hairpin = check_dot_bracket(hairpin_structure, threshold, gap_authorized)


print (f"The most likely hairpin scheme is {hairpin_structure}. \nHairpin structure: {result_hairpin}")
