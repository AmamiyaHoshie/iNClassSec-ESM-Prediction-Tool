# utils/gaac.py
import re
from collections import Counter

def calculate_gaac(sequence):
    """
    计算 GAAC 特征：基于氨基酸分组的组成比例
    """
    group = {
        'alphatic': 'GAVLMI',
        'aromatic': 'FYW',
        'positivecharge': 'KRH',
        'negativecharge': 'DE',
        'uncharge': 'STCPNQ'
    }
    sequence = re.sub('-', '', sequence.upper())
    sequence_length = len(sequence)
    if sequence_length == 0:
        return None
    count = Counter(sequence)
    gaac_vector = []
    for key in group:
        group_count = sum([count[aa] for aa in group[key]])
        gaac_vector.append(group_count / sequence_length)
    return gaac_vector
