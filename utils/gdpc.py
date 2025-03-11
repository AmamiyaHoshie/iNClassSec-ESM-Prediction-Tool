# utils/gdpc.py
import re

def calculate_gdpc(sequence):
    """
    计算 GDPC 特征：基于氨基酸分组的二肽组成比例
    """
    group = {
        'alphaticr': 'GAVLMI',
        'aromatic': 'FYW',
        'postivecharger': 'KRH',
        'negativecharger': 'DE',
        'uncharger': 'STCPNQ'
    }
    group_keys = list(group.keys())
    dipeptides = [g1 + '.' + g2 for g1 in group_keys for g2 in group_keys]
    index = {}
    for key in group_keys:
        for aa in group[key]:
            index[aa] = key
    sequence = re.sub('-', '', sequence.upper())
    sequence_length = len(sequence)
    if sequence_length < 2:
        return None
    myDict = {dipeptide: 0 for dipeptide in dipeptides}
    total = 0
    for j in range(sequence_length - 1):
        aa1 = sequence[j]
        aa2 = sequence[j + 1]
        if aa1 in index and aa2 in index:
            dipeptide = index[aa1] + '.' + index[aa2]
            myDict[dipeptide] += 1
            total += 1
    if total == 0:
        gdpc_vector = [0] * len(dipeptides)
    else:
        gdpc_vector = [myDict[dipeptide] / total for dipeptide in dipeptides]
    return gdpc_vector
