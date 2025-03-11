# utils/aac.py
import re
from collections import Counter

def calculate_aac(sequence, aa_order='ACDEFGHIKLMNPQRSTVWY'):
    """
    计算氨基酸组成（AAC）特征向量
    """
    sequence = re.sub('-', '', sequence)
    sequence_length = len(sequence)
    if sequence_length == 0:
        return None
    aa_count = Counter(sequence)
    aac_vector = [aa_count[aa] / sequence_length for aa in aa_order]
    # 如果全部为0，则返回 None
    if all(value == 0 for value in aac_vector):
        return None
    return aac_vector
