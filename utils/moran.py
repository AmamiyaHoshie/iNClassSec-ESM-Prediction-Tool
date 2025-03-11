# utils/moran.py
import re
import numpy as np

def read_fasta(file_path):
    """
    读取 FASTA 文件，返回序列字符串（适用于单序列读取）
    """
    sequence = ''
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if not line.startswith('>'):
                sequence += line
    return sequence

def load_AAidx(file_path):
    """
    读取 AAindex 文件，返回属性列表、AAidx 数组以及氨基酸索引映射
    """
    AA = 'ARNDCQEGHILKMFPSTWYV'
    with open(file_path) as f:
        records = f.readlines()
    print(f'Total AAindex read: {len(records) - 1}')
    myDict = {}
    props = []
    for i in records[1:]:
        array = i.rstrip().split('\t')
        if 'NA' in array:
            continue
        myDict[array[0]] = array[1:]
        props.append(array[0])
    print(f'Total valid AAindex: {len(props)}')
    AAidx = []
    for i in props:
        AAidx.append([float(x) for x in myDict[i]])
    AAidx = np.array(AAidx)
    # 标准化处理
    propMean = np.mean(AAidx, axis=1)
    propStd = np.std(AAidx, axis=1)
    AAidx = (AAidx - propMean[:, None]) / propStd[:, None]
    index = {}
    for i in range(len(AA)):
        index[AA[i]] = i
    return props, AAidx, index

def calculate_moran(sequence, nlag, props, AAidx, index):
    """
    计算 Moran 特征，nlag 为滞后数
    """
    sequence = re.sub('-', '', sequence.upper())
    N = len(sequence)
    if N < nlag + 1:
        return None
    code = []
    for prop in range(len(props)):
        xmean = np.mean([AAidx[prop][index.get(aa, 0)] for aa in sequence])
        fenmu = np.sum([(AAidx[prop][index.get(aa, 0)] - xmean) ** 2 for aa in sequence]) / N
        for n in range(1, nlag + 1):
            if N > nlag:
                fenzi = np.sum([(AAidx[prop][index.get(sequence[j], 0)] - xmean) *
                                 (AAidx[prop][index.get(sequence[j + n], 0)] - xmean)
                                 for j in range(N - n)]) / (N - n)
                rn = fenzi / fenmu if fenmu != 0 else 0
            else:
                rn = 0
            code.append(rn)
    return code
