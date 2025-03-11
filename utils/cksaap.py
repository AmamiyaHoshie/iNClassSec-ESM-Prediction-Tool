# utils/cksaap.py
import re

def calculate_cksaap(sequence, gap=3, aa_order='ACDEFGHIKLMNPQRSTVWY'):
    """
    计算 CKSAAP 特征向量，gap 默认3
    """
    sequence = re.sub('-', '', sequence)
    sequence_length = len(sequence)
    if sequence_length < gap + 2:
        return None

    aa_pairs = [aa1 + aa2 for aa1 in aa_order for aa2 in aa_order]
    cksaap_vector = []

    for g in range(gap + 1):
        pair_count = {pair: 0 for pair in aa_pairs}
        sum_pairs = 0

        for i in range(sequence_length - g - 1):
            if sequence[i] in aa_order and sequence[i + g + 1] in aa_order:
                pair = sequence[i] + sequence[i + g + 1]
                pair_count[pair] += 1
                sum_pairs += 1

        cksaap_vector.extend([pair_count[pair] / sum_pairs if sum_pairs > 0 else 0 for pair in aa_pairs])
    if all(value == 0 for value in cksaap_vector):
        return None
    return cksaap_vector
