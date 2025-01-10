from datasketch import MinHash, LeanMinHash


def create_minhash(sequence, shingle=5, permutations=128):
    m = MinHash(num_perm=permutations)
    shingles = [sequence[i:i+shingle].encode('utf8') for i in range(len(sequence)-(shingle-1))]
    m.update_batch(shingles)
    return LeanMinHash(m)


# def one_hot_encode(sequence, amino_acids):
#     encoding = np.zeros((len(sequence), len(amino_acids)), dtype=int)  # Initialize a zero matrix
#     for i, amino_acid in enumerate(sequence):
#         if amino_acid in amino_acids:
#             encoding[i, amino_acids.index(amino_acid)] = 1  # Set 1 for the matching amino acid
#     return encoding


# def create_minhash(sequence, amino_acids):
#     encoded_sequence = one_hot_encode(sequence, amino_acids)
#     m = MinHash(num_perm=NUM_PERM)
#     for i in range(encoded_sequence.shape[0]):
#         for j in np.where(encoded_sequence[i] == 1)[0]:
#             # Construct a "shingle" (a substring or feature) from the one-hot encoding
#             shingle = f'{i}_{j}'.encode('utf8')
#             m.update(shingle)
#     return LeanMinHash(m)