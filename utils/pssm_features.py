# utils/pssm_features.py
import numpy as np

def normalizePSSM(PSSM_matrix):
    PSSM = PSSM_matrix[:, :20].astype(float)
    mean_matrix = np.mean(PSSM, axis=1, keepdims=True)
    std_matrix = np.std(PSSM, axis=1, keepdims=True)
    std_matrix[std_matrix == 0] = 1
    PSSM_norm = (PSSM - mean_matrix) / std_matrix
    return PSSM_norm

def pse_pssm(PSSM_matrix, alpha=1):
    PSSM_norm = normalizePSSM(PSSM_matrix)
    L = PSSM_norm.shape[0]  # 序列长度
    avg_pssm = np.mean(PSSM_norm, axis=0)  # 计算均值
    diff_pssm = np.zeros(20)
    for i in range(L - alpha):
        diff_pssm += (PSSM_norm[i] - PSSM_norm[i + alpha]) ** 2
    diff_pssm /= (L - alpha)
    pse_pssm_vector = np.hstack((avg_pssm, diff_pssm))
    return pse_pssm_vector
