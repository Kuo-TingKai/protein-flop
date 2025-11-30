# Protein Flop: Topological Data Analysis for Protein Conformation Transitions

這是一個使用拓撲數據分析 (Topological Data Analysis, TDA) 來分析蛋白質構象轉變的 Python 專案。專案目標是將代數幾何中的極小模型綱領 (Minimal Model Program, MMP) 概念與拓撲數據分析結合，用於研究蛋白質折疊過程中的拓撲特徵變化。

## 專案概述

本專案實現了以下核心功能：

1. **數據處理**：處理分子動力學 (MD) 模擬數據，提取蛋白質構象座標
2. **拓撲數據分析**：使用持續同調 (Persistent Homology) 和 Vietoris-Rips 過濾來計算拓撲特徵
3. **構象比較**：比較折疊態和展開態的拓撲差異
4. **可視化**：生成持續圖 (Persistence Diagrams) 和 Betti 曲線

## 理論背景

根據 `protein-flop.md` 中的討論，本專案旨在：

- 使用 β-髮夾肽 (Beta-Hairpin Peptide) 作為計算範例
- 分析折疊態 (Folded State) 和展開中間態 (Unfolded State) 的拓撲差異
- 觀察 Betti 數 (β₀, β₁, β₂) 的變化，特別是轉角區域 (Turn Region) 的拓撲特徵
- 為未來建立典範因子 K_X 與分子自由能 G 之間的定量橋樑提供初步證據

## 安裝

### 系統需求

- Python 3.8 或更高版本
- pip 或 conda

### 安裝步驟

1. 克隆或下載此專案

2. 安裝依賴套件：

```bash
pip install -r requirements.txt
```

主要依賴套件包括：
- `numpy`: 數值計算
- `scipy`: 科學計算
- `matplotlib`: 繪圖
- `ripser`: 持續同調計算（主要 TDA 工具）
- `mdtraj`: 處理 MD 軌跡數據
- `biopython`: 生物資訊學工具
- `pandas`: 數據處理
- `scikit-learn`: 機器學習工具

## 使用方法

### 基本使用

執行主程式進行分析：

```bash
python main.py
```

使用預設參數（1000 個樣本，6 個殘基）生成合成數據並進行分析。

### 自訂參數

```bash
python main.py --n-samples 2000 --n-residues 8 --max-dim 2 --output-dir results
```

參數說明：
- `--n-samples`: 生成的構象數量（預設：1000）
- `--n-residues`: 轉角區域的殘基數量（預設：6）
- `--max-dim`: 最大同調維度（預設：2）
- `--output-dir`: 輸出目錄（預設：output）
- `--no-show`: 不顯示互動式圖表
- `--use-rmsd`: 使用 RMSD 距離矩陣（預設：True）

### 範例腳本

執行範例腳本查看詳細使用範例：

```bash
python example_usage.py
```

此腳本包含三個範例：
1. 單一構象狀態的基本分析
2. 折疊態與展開態的比較
3. 使用 RMSD 距離矩陣的自訂分析

## 專案結構

```
protein-flop/
├── src/
│   ├── __init__.py
│   ├── data_processing.py      # 數據處理模組
│   ├── tda_computation.py       # TDA 計算模組
│   ├── visualization.py         # 可視化模組
│   └── analysis.py              # 主要分析類別
├── main.py                      # 主程式
├── example_usage.py             # 使用範例
├── requirements.txt             # 依賴套件
├── README.md                    # 本文件
└── protein-flop.md            # 理論背景文件
```

## 核心功能

### 1. 數據處理 (`data_processing.py`)

- `extract_ca_coordinates()`: 提取 C-alpha 原子座標
- `extract_turn_region_coordinates()`: 提取轉角區域座標
- `compute_rmsd_matrix()`: 計算 RMSD 距離矩陣
- `generate_synthetic_conformation_data()`: 生成合成構象數據（用於測試）

### 2. TDA 計算 (`tda_computation.py`)

- `compute_persistent_homology()`: 計算持續同調
- `compute_betti_numbers()`: 計算 Betti 數
- `compute_persistence_statistics()`: 計算持續性統計
- `compare_persistence_diagrams()`: 比較兩個持續圖

### 3. 可視化 (`visualization.py`)

- `plot_persistence_diagram()`: 繪製持續圖
- `plot_betti_curve()`: 繪製 Betti 曲線
- `compare_persistence_diagrams_plot()`: 並排比較兩個持續圖
- `plot_betti_comparison()`: 比較兩個狀態的 Betti 曲線

### 4. 分析類別 (`analysis.py`)

`ProteinConformationAnalyzer` 類別提供完整的分析流程：

```python
from src.analysis import ProteinConformationAnalyzer

analyzer = ProteinConformationAnalyzer(max_dim=2)
comparison_results = analyzer.compare_states(
    folded_coordinates,
    unfolded_coordinates,
    labels=('Folded', 'Unfolded')
)
analyzer.visualize_comparison(comparison_results)
analyzer.print_comparison_summary(comparison_results)
```

## 輸出結果

分析完成後，會在輸出目錄中生成以下圖表：

1. **持續圖 (Persistence Diagrams)**: 
   - `persistence_diagram_dim0.png`: 連通分量 (β₀)
   - `persistence_diagram_dim1.png`: 環/洞 (β₁)
   - `persistence_diagram_dim2.png`: 空腔 (β₂)

2. **Betti 曲線**:
   - `betti_curve_dim1.png`: 維度 1 的 Betti 數隨過濾值變化
   - `betti_curve_dim2.png`: 維度 2 的 Betti 數隨過濾值變化

## 使用真實 MD 數據

要使用真實的 MD 模擬數據，請修改數據載入部分：

```python
from src.data_processing import load_trajectory_from_file, extract_turn_region_coordinates
import mdtraj as md

# 載入 MD 軌跡
trajectory = load_trajectory_from_file('trajectory.xtc', top_file='topology.pdb')

# 提取轉角區域座標（例如殘基 4-10）
turn_coords = extract_turn_region_coordinates(trajectory, turn_start=4, turn_end=10)

# 進行分析
analyzer = ProteinConformationAnalyzer(max_dim=2)
results = analyzer.analyze_conformation(turn_coords, use_rmsd=True)
```

## 預期結果

根據理論分析，我們預期觀察到：

1. **持續圖差異**: 折疊態和展開態的持續圖應該不同，證明兩種構象在局部拓撲上不是同胚的。

2. **β₁ 特徵的持久性差異**: 折疊態中持久性最強的 β₁ 特徵（環/洞）的 L 值應該比展開態大，代表折疊態的拓撲更穩定或更緊密。

3. **β₂ 特徵的存在**: 折疊態中可能出現穩定的 β₂ 特徵（空腔），代表轉角區域形成緊密的封閉結構，而在展開態中 β₂ 特徵可能消失。

## 未來發展

- [ ] 整合真實的 MD 模擬數據
- [ ] 實現更複雜的距離度量（如 Gromov-Wasserstein 距離）
- [ ] 開發與 MMP 理論的連接框架
- [ ] 添加機器學習模型來預測構象轉變
- [ ] 支援更多蛋白質結構類型

## 參考文獻

本專案基於以下理論背景：
- 拓撲數據分析 (Topological Data Analysis)
- 持續同調 (Persistent Homology)
- 極小模型綱領 (Minimal Model Program, MMP)
- 蛋白質折疊機制

詳細理論討論請參見 `protein-flop.md`。

## 授權

本專案僅供學術研究使用。

## 聯絡資訊

如有問題或建議，請開啟 Issue 或提交 Pull Request。

