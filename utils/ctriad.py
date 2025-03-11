# utils/ctriad.py
import re

def calculate_ctriad(sequence):
    """
    计算 CTriad 特征：基于氨基酸三联体统计特征
    """
    AAGroup = {
        'g1': 'AGV',
        'g2': 'ILFP',
        'g3': 'YMTS',
        'g4': 'HNQW',
        'g5': 'RK',
        'g6': 'DE',
        'g7': 'C'
    }
    myGroups = sorted(AAGroup.keys())  # ['g1', 'g2', 'g3', 'g4', 'g5', 'g6', 'g7']
    AADict = {}
    for g in myGroups:
        for aa in AAGroup[g]:
            AADict[aa] = g
    features = [f1 + '.' + f2 + '.' + f3 for f1 in myGroups for f2 in myGroups for f3 in myGroups]
    sequence = re.sub('-', '', sequence.upper())
    sequence_length = len(sequence)
    if sequence_length < 3:
        return None
    myDict = {f: 0 for f in features}
    for i in range(sequence_length - 2):
        aa1 = sequence[i]
        aa2 = sequence[i + 1]
        aa3 = sequence[i + 2]
        if aa1 not in AADict or aa2 not in AADict or aa3 not in AADict:
            continue
        key = AADict[aa1] + '.' + AADict[aa2] + '.' + AADict[aa3]
        myDict[key] += 1
    total = sum(myDict.values())
    if total == 0:
        return [0] * len(features)
    ctriad_vector = [myDict[f] / total for f in features]
    return ctriad_vector
