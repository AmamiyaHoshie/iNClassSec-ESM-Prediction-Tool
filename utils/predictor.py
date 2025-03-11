# utils/predictor.py
import os
import joblib
import torch
import torch.nn as nn
import numpy as np


# 定义 DNN 模型，注意 input_dim 需与 ESM3 嵌入维度匹配
class DeepFeedForwardClassifier(nn.Module):
    def __init__(self, input_dim, hidden_dim=128, dropout=0.5):
        super(DeepFeedForwardClassifier, self).__init__()
        self.layer1 = nn.Linear(input_dim, hidden_dim)
        self.bn1 = nn.BatchNorm1d(hidden_dim)
        self.relu1 = nn.LeakyReLU()
        self.dropout1 = nn.Dropout(dropout)

        self.layer2 = nn.Linear(hidden_dim, hidden_dim)
        self.bn2 = nn.BatchNorm1d(hidden_dim)
        self.relu2 = nn.LeakyReLU()
        self.dropout2 = nn.Dropout(dropout)

        self.layer3 = nn.Linear(hidden_dim, hidden_dim)
        self.bn3 = nn.BatchNorm1d(hidden_dim)
        self.relu3 = nn.LeakyReLU()
        self.dropout3 = nn.Dropout(dropout)

        self.fc = nn.Linear(hidden_dim, 1)

    def forward(self, x):
        out = self.layer1(x)
        out = self.bn1(out)
        out = self.relu1(out)
        out = self.dropout1(out)

        out = self.layer2(out)
        out = self.bn2(out)
        out = self.relu2(out)
        out = self.dropout2(out)

        out = self.layer3(out)
        out = self.bn3(out)
        out = self.relu3(out)
        out = self.dropout3(out)

        out = self.fc(out)
        return out.squeeze(1)


def load_models(input_dim=1536):
    """
    加载 XGBoost 管道、元模型和 DNN 模型。
    参数：
      input_dim: DNN 模型输入维度，需与训练时使用的 ESM3 嵌入维度一致
    """
    xgb_pipeline = joblib.load(os.path.join("models", "xgb_pipeline.pkl"))
    meta_model = joblib.load(os.path.join("models", "meta_model.pkl"))
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    dnn_model = DeepFeedForwardClassifier(input_dim=input_dim, hidden_dim=128, dropout=0.5).to(device)
    dnn_model.load_state_dict(torch.load(os.path.join("models", "best_model.pth"), map_location=device))
    dnn_model.eval()
    return xgb_pipeline, meta_model, dnn_model, device


def predict_features(xgb_pipeline, meta_model, dnn_model, device, xgb_features, esm3_embeddings):
    """
    参数：
      xgb_features: numpy 数组，来自最终拼接后的特征（例如 pse-pssm 与其他序列特征的拼接），用于 XGBoost 模型
      esm3_embeddings: numpy 数组，来自 ESM3 npy 文件，用于 DNN 模型
    两者行数需与序列数一致。

    预测流程：
      1. 分别用 XGBoost 管道和 DNN 模型得到预测概率；
      2. 将两个模型的预测结果堆叠，再通过元模型进行融合，返回最终预测概率。
    """
    # 第一阶段：XGBoost 预测
    test_preds_xgb = xgb_pipeline.predict_proba(xgb_features)[:, 1]

    # 第二阶段：DNN 模型预测（使用 ESM3 嵌入）
    from torch.utils.data import DataLoader, TensorDataset
    tensor_embeddings = torch.tensor(esm3_embeddings, dtype=torch.float32)
    dataset = TensorDataset(tensor_embeddings)
    loader = DataLoader(dataset, batch_size=64, shuffle=False)

    dnn_preds = []
    with torch.no_grad():
        for batch in loader:
            emb = batch[0].to(device)
            outputs = dnn_model(emb)
            preds = torch.sigmoid(outputs)
            dnn_preds.extend(preds.cpu().numpy())
    dnn_preds = np.array(dnn_preds)

    # 堆叠两个模型的预测结果，再通过元模型融合（stacking）
    X_stack = np.vstack((test_preds_xgb, dnn_preds)).T
    meta_probs = meta_model.predict_proba(X_stack)[:, 1]
    return meta_probs
