# utils/dpc.py
import re

def calculate_dpc(sequence, aa_order='ACDEFGHIKLMNPQRSTVWY'):
    """
    计算 DPC（二肽组成）特征向量
    """
    sequence = re.sub('-', '', sequence)
    if len(sequence) < 2:
        return None
    dipeptides = [aa1 + aa2 for aa1 in aa_order for aa2 in aa_order]
    dipeptide_count = {dipeptide: 0 for dipeptide in dipeptides}
    for i in range(len(sequence) - 1):
        dipeptide = sequence[i] + sequence[i + 1]
        if dipeptide in dipeptide_count:
            dipeptide_count[dipeptide] += 1
    total_dipeptides = sum(dipeptide_count.values())
    dpc_vector = [dipeptide_count[dipeptide] / total_dipeptides if total_dipeptides > 0 else 0 for dipeptide in dipeptides]
    if all(value == 0 for value in dpc_vector):
        return None
    return dpc_vector
