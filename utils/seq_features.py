from utils.aac import calculate_aac
from utils.cksaap import calculate_cksaap
from utils.ctdc import calculate_ctdc
from utils.ctdd import calculate_ctdd
from utils.ctdt import calculate_ctdt
from utils.ctriad import calculate_ctriad
from utils.dpc import calculate_dpc
from utils.gaac import calculate_gaac
from utils.gdpc import calculate_gdpc
from utils.moran import calculate_moran


def calculate_all_seq_features(sequence, moran_params):
    """
    对单个序列计算所有预定义的序列特征，并将各个特征向量拼接成一个长向量。
    :param sequence: 字符串格式的蛋白质序列
    :param moran_params: 一个字典，包含计算 Moran 特征需要的参数，比如 { 'nlag': 2, 'props': props, 'AAidx': AAidx, 'index': index }
    :return: 拼接后的特征向量（列表或一维数组）
    """
    feature_list = []
    # 每个元组包括特征名称和对应的计算函数
    feature_funcs = [
        ('AAC', calculate_aac),
        ('CKSAAP', calculate_cksaap),
        ('CTDC', calculate_ctdc),
        ('CTDD', calculate_ctdd),
        ('CTDT', calculate_ctdt),
        ('CTriad', calculate_ctriad),
        ('DPC', calculate_dpc),
        ('GAAC', calculate_gaac),
        ('GDPC', calculate_gdpc),
        ('Moran', lambda seq: calculate_moran(seq, **moran_params))
    ]

    for name, func in feature_funcs:
        feat = func(sequence)
        if feat is None:
            # 如果某个特征计算返回 None，则可选择填充全零（这里简化处理：跳过该特征）
            print(f"Warning: {name} 特征计算失败，跳过。")
            continue
        feature_list.extend(feat)
    return feature_list