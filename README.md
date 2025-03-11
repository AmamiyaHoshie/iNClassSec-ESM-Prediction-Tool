# iNClassSec-ESM Prediction Tool

## Introduction

The **iNClassSec-ESM Prediction Tool** is a companion software derived from the research article *“iNClassSec-ESM: Discovering Potential Non-Classical Secreted Proteins through a Novel Protein Language Model”*. The tool aims to provide a user-friendly and highly accessible platform for predicting Non-Classical Secreted Proteins (NCSPs).

Based on the methods described in the paper, this tool integrates traditional handcrafted sequence features with hidden-layer embeddings derived from the ESM3 protein language model, enabling users to rapidly and effectively predict potential NCSPs.

**Original Research Project:**  
[github.com/AmamiyaHoshie/iNClassSec-ESM](https://github.com/AmamiyaHoshie/iNClassSec-ESM)

## Requirements

- **Operating System:** Windows 10 or later, Linux, macOS
- **Python Version:** Python 3.10
- **Dependencies:** See `requirements.txt` provided in the project repository

## Installation

Follow these steps to install and run the **iNClassSec-ESM Prediction Tool** from source:

### Step 1: Download Source Code

```bash
git clone https://github.com/AmamiyaHoshie/iNClassSec-ESM-Prediction-Tool.git
cd iNClassSec-ESM-Prediction-Tool
```

### Step 2: Create and Activate a Virtual Environment (Recommended)

- **Using Conda:**
  ```bash
  conda create -n inclasssec python=3.10
  conda activate inclasssec
  ```

- **Using Python venv:**
  ```bash
  python -m venv env
  # Activate the environment:
  source env/bin/activate        # Linux/macOS
  .\env\Scripts\activate         # Windows
  ```

### Step 3: Install Dependencies

```bash
pip install -r requirements.txt
```

### Step 4: Launch the Application

```bash
python main.py
```

## Usage

The tool interface consists of two main pages: **Input Page** and **Results Page**.

### Input Page

1. **FASTA File:**  
   - Click the "Browse" button to select a FASTA file containing multiple protein sequences.
   - Ensure the file is in standard FASTA format, where each header line starts with `>` and includes a unique identifier.

2. **PSSM Folder:**  
   - Click "Browse" to choose the folder containing the PSSM matrix files.
   - Each file in the folder should be named using a pure number that corresponds to the sequence order in the FASTA file.

3. **ESM3 Embedding File:**  
   - Click "Browse" to select the `.npy` file containing ESM3 embedding features.
   - Ensure that the `.npy` file is generated using the provided `ESM3_embedding_extraction.py` script.

4. **Start Prediction:**  
   - After selecting the files, click the "Predict" button.
   - A progress bar will appear below the button, updating in real time as each sequence is processed.

### Results Page

After the prediction process is complete, the application will automatically switch to the Results Page:

1. **Prediction Results Display:**  
   - The left panel displays sequence IDs predicted as positive (NCSPs).
   - The right panel displays sequence IDs predicted as negative.

2. **Save Results:**  
   - Click the "Save Results" button to open a dialog for saving the prediction results as a CSV file.
   - Once saved, the application will automatically exit.

## FAQ

**Q1: I get an error related to the ESM3 embedding file. What should I do?**  
A: Ensure that the `.npy` file is generated using the provided `ESM3_embedding_extraction.py` script without any modifications.

**Q2: My PSSM files are not being recognized.**  
A: Verify that the PSSM folder contains files named solely with pure numbers that correspond exactly to the sequence order in the FASTA file.

**Q3: The prediction process is slow or the progress bar does not update.**  
A: The processing time depends on the number of sequences and the computational complexity of feature extraction and model prediction. For more details, check the command-line output during execution.
