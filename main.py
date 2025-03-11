# main.py
import os
import sys
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import numpy as np
import pandas as pd

# Import utility functions
from utils.file_io import read_fasta, read_pssm_file
from utils.pssm_features import pse_pssm
from utils.seq_features import calculate_all_seq_features
from utils.moran import load_AAidx
from utils.predictor import load_models, predict_features


class PredictorApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("iNClassSec-ESM Prediction Tool")
        self.geometry("600x300")
        self.resizable(True, True)

        # Variables to store paths
        self.fasta_file = tk.StringVar()
        self.pssm_folder = tk.StringVar()
        self.esm3_file = tk.StringVar()

        # Prediction results variables
        self.fasta_ids = []
        self.prediction_df = None
        self.positive_ids = []
        self.negative_ids = []

        self.create_input_frame()
        self.create_result_frame()

        self.input_frame.pack(fill="both", expand=True)
        self.result_frame.pack_forget()

    def create_input_frame(self):
        self.input_frame = tk.Frame(self)

        title = tk.Label(self.input_frame, text="Select Samples for Prediction", font=("Arial", 16))
        title.pack(pady=10)

        # FASTA File selection
        frame1 = tk.Frame(self.input_frame)
        frame1.pack(pady=5, padx=10, fill="x", anchor="center")
        tk.Label(frame1, text="FASTA File:", width=20).pack(side="left")
        tk.Entry(frame1, textvariable=self.fasta_file, width=50).pack(side="left", padx=5)
        tk.Button(frame1, text="Browse", command=self.browse_fasta).pack(side="left")

        # PSSM Folder selection
        frame2 = tk.Frame(self.input_frame)
        frame2.pack(pady=5, padx=10, fill="x", anchor="center")
        tk.Label(frame2, text="PSSM Folder:", width=20).pack(side="left")
        tk.Entry(frame2, textvariable=self.pssm_folder, width=50).pack(side="left", padx=5)
        tk.Button(frame2, text="Browse", command=self.browse_pssm_folder).pack(side="left")

        # ESM3 Embedding File selection
        frame3 = tk.Frame(self.input_frame)
        frame3.pack(pady=5, padx=10, fill="x", anchor="center")
        tk.Label(frame3, text="ESM3 Embedding File:", width=20).pack(side="left")
        tk.Entry(frame3, textvariable=self.esm3_file, width=50).pack(side="left", padx=5)
        tk.Button(frame3, text="Browse", command=self.browse_esm3_file).pack(side="left")

        # Predict button
        self.predict_button = tk.Button(self.input_frame, text="Predict", command=self.run_prediction, width=20,
                                        height=2)
        self.predict_button.pack(pady=10)

        # Progress bar placed under the Predict button
        self.progress_bar = ttk.Progressbar(self.input_frame, orient="horizontal", length=400, mode="determinate")
        self.progress_bar.pack(pady=5)
        self.progress_bar["value"] = 0

    def create_result_frame(self):
        self.result_frame = tk.Frame(self)

        title = tk.Label(self.result_frame, text="Prediction Results", font=("Arial", 16))
        title.pack(pady=10)

        # Frame for two listboxes (with scrollbars)
        lists_frame = tk.Frame(self.result_frame)
        lists_frame.pack(fill="both", expand=True, padx=10, pady=10)

        # Positive samples listbox
        pos_frame = tk.Frame(lists_frame)
        pos_frame.pack(side="left", fill="both", expand=True, padx=5)
        tk.Label(pos_frame, text="Positive Samples").pack()
        self.pos_listbox = tk.Listbox(pos_frame)
        self.pos_listbox.pack(side="left", fill="both", expand=True)
        pos_scroll = tk.Scrollbar(pos_frame, orient="vertical", command=self.pos_listbox.yview)
        pos_scroll.pack(side="right", fill="y")
        self.pos_listbox.config(yscrollcommand=pos_scroll.set)

        # Negative samples listbox
        neg_frame = tk.Frame(lists_frame)
        neg_frame.pack(side="left", fill="both", expand=True, padx=5)
        tk.Label(neg_frame, text="Negative Samples").pack()
        self.neg_listbox = tk.Listbox(neg_frame)
        self.neg_listbox.pack(side="left", fill="both", expand=True)
        neg_scroll = tk.Scrollbar(neg_frame, orient="vertical", command=self.neg_listbox.yview)
        neg_scroll.pack(side="right", fill="y")
        self.neg_listbox.config(yscrollcommand=neg_scroll.set)

        # Save button
        tk.Button(self.result_frame, text="Save Results", command=self.save_results, width=20).pack(pady=10)

    def browse_fasta(self):
        path = filedialog.askopenfilename(title="Select FASTA File",
                                          filetypes=[("FASTA files", "*.fasta;*.fa;*.txt"), ("All files", "*.*")])
        if path:
            self.fasta_file.set(path)

    def browse_pssm_folder(self):
        folder = filedialog.askdirectory(title="Select PSSM Folder")
        if folder:
            self.pssm_folder.set(folder)

    def browse_esm3_file(self):
        path = filedialog.askopenfilename(title="Select ESM3 Embedding File",
                                          filetypes=[("Numpy files", "*.npy"), ("All files", "*.*")])
        if path:
            self.esm3_file.set(path)

    def run_prediction(self):
        # Validate inputs
        if not self.fasta_file.get():
            messagebox.showerror("Error", "No FASTA file selected!")
            return
        if not self.pssm_folder.get():
            messagebox.showerror("Error", "No PSSM folder selected!")
            return
        if not self.esm3_file.get():
            messagebox.showerror("Error", "No ESM3 embedding file selected!")
            return

        try:
            # Reset progress bar and set max to 100%
            self.progress_bar["value"] = 0
            self.update_idletasks()

            # 1. Read FASTA file; expecting (ids, sequences)
            fasta_ids, sequences = read_fasta(self.fasta_file.get())
            if not sequences:
                messagebox.showerror("Error", "No sequences found in FASTA file!")
                return
            n_seq = len(sequences)
            print(f"Read {n_seq} sequences.")

            # We'll update progress bar per sequence processed
            progress_increment = 100 / n_seq

            # 2. Load AAidx data for Moran feature
            aaidx_file = os.path.join("resources", "AAidx.txt")
            if not os.path.exists(aaidx_file):
                messagebox.showerror("Error", "AAidx.txt not found!")
                return
            props, AAidx, index = load_AAidx(aaidx_file)
            moran_params = {'nlag': 2, 'props': props, 'AAidx': AAidx, 'index': index}

            # 3. For each sequence, calculate sequence features and PSSM feature, then combine them
            final_features_list = []
            for i, seq in enumerate(sequences, start=1):
                print(f"Processing sequence {i}, length: {len(seq)}")
                seq_features = calculate_all_seq_features(seq, moran_params)
                if seq_features is None:
                    print(f"Warning: Sequence {i} feature calculation failed.")
                    # 若失败，可填充全0向量，长度需与预期保持一致
                    seq_features = [0] * (self.get_expected_seq_feature_dim())  # 请实现 get_expected_seq_feature_dim()
                # 读取对应 PSSM 文件（文件名为纯数字）
                pssm_file = os.path.join(self.pssm_folder.get(), f"{i}")
                if not os.path.exists(pssm_file):
                    messagebox.showerror("Error", f"PSSM file not found: {pssm_file}")
                    return
                print(f"Processing PSSM file: {pssm_file}")
                pssm_matrix = read_pssm_file(pssm_file)
                pssm_feature = pse_pssm(pssm_matrix)
                if pssm_feature is None:
                    print(f"Warning: PSSM feature calculation failed for file: {pssm_file}")
                    # 填充全0向量（长度按预期）
                    pssm_feature = [0] * self.get_expected_pssm_feature_dim()  # 请实现 get_expected_pssm_feature_dim()
                # Concatenate sequence features and PSSM feature for this sequence
                combined_feature = np.hstack((seq_features, pssm_feature))
                final_features_list.append(combined_feature)

                # Update progress bar
                self.progress_bar["value"] += progress_increment
                self.update_idletasks()

            final_features = np.array(final_features_list)
            print(f"Final features shape (XGB input): {final_features.shape}")

            # 4. Load ESM3 embeddings (for DNN input)
            esm3_embeddings = np.load(self.esm3_file.get(), allow_pickle=True)
            if esm3_embeddings.shape[0] != n_seq:
                messagebox.showerror("Error", "Number of ESM3 embeddings does not match number of sequences!")
                return
            print(f"ESM3 embeddings shape: {esm3_embeddings.shape}")

            # 5. Load models; set DNN input_dim to match ESM3 embedding dimension
            esm3_dim = esm3_embeddings.shape[1]
            xgb_pipeline, meta_model, dnn_model, device = load_models(input_dim=esm3_dim)
            print("Models loaded.")

            # 6. Predict using both models and stacking via meta_model
            meta_probs = predict_features(xgb_pipeline, meta_model, dnn_model, device, final_features, esm3_embeddings)
            print("Prediction complete.")

            # 7. Apply threshold to get prediction labels
            threshold = 0.15
            preds = (meta_probs >= threshold).astype(int)

            # Build positive and negative ID lists using FASTA ids
            pos_ids = [fid for fid, p in zip(fasta_ids, preds) if p == 1]
            neg_ids = [fid for fid, p in zip(fasta_ids, preds) if p == 0]
            self.positive_ids = pos_ids
            self.negative_ids = neg_ids
            self.prediction_df = pd.DataFrame({
                "Sequence_ID": fasta_ids,
                "Prediction_Probability": meta_probs,
                "Prediction_Label": preds
            })

            # Reset progress bar to 100%
            self.progress_bar["value"] = 100
            self.update_idletasks()

            self.show_results()

        except Exception as e:
            messagebox.showerror("Error", f"An error occurred:\n{e}")
            print(f"Error: {e}")

    def show_results(self):
        # Clear listboxes
        self.pos_listbox.delete(0, tk.END)
        self.neg_listbox.delete(0, tk.END)

        for sid in self.positive_ids:
            self.pos_listbox.insert(tk.END, str(sid))
        for sid in self.negative_ids:
            self.neg_listbox.insert(tk.END, str(sid))

        # Resize window for result view (self-adaptive)
        self.geometry("800x600")
        self.input_frame.pack_forget()
        self.result_frame.pack(fill="both", expand=True)

    def create_result_frame(self):
        self.result_frame = tk.Frame(self)

        title = tk.Label(self.result_frame, text="Prediction Results", font=("Arial", 16))
        title.pack(pady=10)

        # Frame for two listboxes (with scrollbars)
        lists_frame = tk.Frame(self.result_frame)
        lists_frame.pack(fill="both", expand=True, padx=10, pady=10)

        # Positive samples listbox
        pos_frame = tk.Frame(lists_frame)
        pos_frame.pack(side="left", fill="both", expand=True, padx=5)
        tk.Label(pos_frame, text="Positive Samples").pack()
        self.pos_listbox = tk.Listbox(pos_frame)
        self.pos_listbox.pack(side="left", fill="both", expand=True)
        pos_scroll = tk.Scrollbar(pos_frame, orient="vertical", command=self.pos_listbox.yview)
        pos_scroll.pack(side="right", fill="y")
        self.pos_listbox.config(yscrollcommand=pos_scroll.set)

        # Negative samples listbox
        neg_frame = tk.Frame(lists_frame)
        neg_frame.pack(side="left", fill="both", expand=True, padx=5)
        tk.Label(neg_frame, text="Negative Samples").pack()
        self.neg_listbox = tk.Listbox(neg_frame)
        self.neg_listbox.pack(side="left", fill="both", expand=True)
        neg_scroll = tk.Scrollbar(neg_frame, orient="vertical", command=self.neg_listbox.yview)
        neg_scroll.pack(side="right", fill="y")
        self.neg_listbox.config(yscrollcommand=neg_scroll.set)

        # Save button
        tk.Button(self.result_frame, text="Save Results", command=self.save_results, width=20).pack(pady=10)

    def browse_fasta(self):
        path = filedialog.askopenfilename(title="Select FASTA File",
                                          filetypes=[("FASTA files", "*.fasta;*.fa;*.txt"), ("All files", "*.*")])
        if path:
            self.fasta_file.set(path)

    def browse_pssm_folder(self):
        folder = filedialog.askdirectory(title="Select PSSM Folder")
        if folder:
            self.pssm_folder.set(folder)

    def browse_esm3_file(self):
        path = filedialog.askopenfilename(title="Select ESM3 Embedding File",
                                          filetypes=[("Numpy files", "*.npy"), ("All files", "*.*")])
        if path:
            self.esm3_file.set(path)

    def save_results(self):
        save_path = filedialog.asksaveasfilename(title="Save Prediction Results", defaultextension=".csv",
                                                 filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
        if save_path:
            self.prediction_df.to_csv(save_path, index=False)
            messagebox.showinfo("Saved", f"Results saved to {save_path}")
            self.destroy()  # Exit application


if __name__ == "__main__":
    app = PredictorApp()
    app.mainloop()
