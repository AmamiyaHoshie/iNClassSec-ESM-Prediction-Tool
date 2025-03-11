# utils/file_io.py
import os
import re
import numpy as np

def read_fasta(file_path):
    ids = []
    sequences = []
    seq = ""
    seq_id = ""
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq:
                    sequences.append(seq)
                    seq = ""
                seq_id = line[1:].split()[0]  # 取 header 第一个 token 作为 id
                ids.append(seq_id)
            else:
                seq += line
        if seq:
            sequences.append(seq)
    return ids, sequences

def read_pssm_file(file_path):
    """读取单个 PSSM 文件，并返回矩阵（仅提取第3到第22列）"""
    PSSM_matrix = []
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) > 42:
                # 提取第3到第22列（注意：索引从0开始）
                PSSM_matrix.append(parts[2:22])
    return np.array(PSSM_matrix)
