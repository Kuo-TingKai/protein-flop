# Protein Flop: Topological Data Analysis for Protein Conformation Transitions

這是一個使用拓撲數據分析 (Topological Data Analysis, TDA) 來分析蛋白質構象轉變的 Python 專案。專案目標是將代數幾何中的極小模型綱領 (Minimal Model Program, MMP) 概念與拓撲數據分析結合，用於研究蛋白質折疊過程中的拓撲特徵變化。

## 專案概述

本專案實現了以下核心功能：

1. **數據處理**：處理分子動力學 (MD) 模擬數據，提取蛋白質構象座標
2. **拓撲數據分析**：使用持續同調 (Persistent Homology) 和 Vietoris-Rips 過濾來計算拓撲特徵
3. **Flip/Flop 計算**：實現阿蒂亞 Flop 和蛋白質構象的 Flop 類比操作
4. **構象比較**：比較折疊態和展開態的拓撲差異
5. **可視化**：生成持續圖 (Persistence Diagrams) 和 Betti 曲線

## 理論背景

根據 `protein-flop.md` 中的討論，本專案旨在：

- 實現阿蒂亞 Flop (Atiyah Flop) 的計算範例，展示代數幾何操作如何影響拓撲結構
- 使用 β-髮夾肽 (Beta-Hairpin Peptide) 作為計算範例
- 分析折疊態 (Folded State) 和展開中間態 (Unfolded State) 的拓撲差異
- 觀察 Betti 數 (β₀, β₁, β₂) 的變化，特別是轉角區域 (Turn Region) 的拓撲特徵
- 建立 TDA 與 MMP 之間的理論橋樑
- 為未來建立典範因子 K_X 與分子自由能 G 之間的定量橋樑提供初步證據

## 實驗結果

本專案已實現並驗證了以下計算範例：

### 阿蒂亞 Flop 範例
- 成功實現了阿蒂亞 Flop 的點雲採樣（Flop 前後各 300 個點）
- 觀察到 Flop 後 $H_1$ 特徵的最大持久性從 0.0644 增加到 0.1672
- 證明了 TDA 能夠捕捉代數幾何操作（如 Flop）的拓撲效果

### 蛋白質 Flop 類比範例
- 對 β-髮夾肽的折疊態和展開態各生成了 300 個構象
- 折疊態顯示了更多的 $\beta_2$ 特徵（280 vs 205），表明轉角區域形成緊密結構
- 使用 RMSD 距離矩陣進行 TDA 分析，觀察到明顯的拓撲差異

詳細結果請參見 `protein-flop.tex` 和 `output/` 目錄中的可視化圖表。

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

#### 基本使用範例

執行基本範例腳本：

```bash
python example_usage.py
```

此腳本包含三個範例：
1. 單一構象狀態的基本分析
2. 折疊態與展開態的比較
3. 使用 RMSD 距離矩陣的自訂分析

#### Flip/Flop 計算範例

執行 Flip/Flop 計算範例（**推薦**）：

```bash
cd examples
python flop_example.py
```

此腳本包含兩個主要範例：
1. **阿蒂亞 Flop**：數學基礎範例，展示代數幾何中的 Flop 操作
2. **蛋白質 Flop 類比**：將 Flop 操作應用到 β-髮夾肽的折疊/展開轉變

範例會生成持續圖比較，保存在 `output/` 目錄中。

## 專案結構

```
protein-flop/
├── src/
│   ├── __init__.py
│   ├── data_processing.py      # 數據處理模組
│   ├── tda_computation.py       # TDA 計算模組
│   ├── visualization.py         # 可視化模組
│   ├── analysis.py              # 主要分析類別
│   └── flop_computation.py      # Flip/Flop 計算模組
├── examples/
│   ├── flop_example.py          # Flip/Flop 計算範例（推薦）
│   └── README.md                # 範例說明
├── output/                      # 輸出目錄（可視化圖表）
├── main.py                      # 主程式
├── example_usage.py             # 基本使用範例
├── requirements.txt             # 依賴套件
├── README.md                    # 本文件
├── protein-flop.md            # 理論背景文件
└── protein-flop.tex            # 學術論文（含實驗結果）
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

### 5. Flip/Flop 計算 (`flop_computation.py`)

實現阿蒂亞 Flop 和蛋白質 Flop 類比：

```python
from src.flop_computation import (
    sample_atiyah_flop_before,
    sample_atiyah_flop_after,
    protein_flop_analogy_folded,
    protein_flop_analogy_unfolded
)

# 阿蒂亞 Flop
points_before = sample_atiyah_flop_before(n_samples=300)
points_after = sample_atiyah_flop_after(n_samples=300)

# 蛋白質 Flop 類比
folded_coords = protein_flop_analogy_folded(n_samples=300)
unfolded_coords = protein_flop_analogy_unfolded(n_samples=300)
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

3. **Flip/Flop 範例結果**:
   - `atiyah_flop_dim1.png`: 阿蒂亞 Flop 維度 1 持續圖比較
   - `atiyah_flop_dim2.png`: 阿蒂亞 Flop 維度 2 持續圖比較
   - `protein_flop_dim1.png`: 蛋白質 Flop 類比維度 1 持續圖比較
   - `protein_flop_dim2.png`: 蛋白質 Flop 類比維度 2 持續圖比較

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

## 實驗結果摘要

根據實際計算實驗，我們觀察到：

1. **阿蒂亞 Flop 結果**:
   - Flop 後 $H_1$ 特徵的最大持久性從 0.0644 增加到 0.1672
   - Flop 後 $H_1$ 特徵數量從 111 增加到 143
   - 證明了 TDA 能夠捕捉代數幾何操作（如 Flop）的拓撲效果

2. **蛋白質 Flop 類比結果**:
   - 折疊態有更多的 $\beta_2$ 特徵（280 vs 205），表明轉角區域形成緊密結構
   - 折疊態和展開態在拓撲上確實不同，證明兩種構象在局部拓撲上不是同胚的
   - 展開態的 $\beta_1$ 最大持久性更高（0.0248 vs 0.0177），可能反映部分折疊中間態的複雜拓撲

3. **持續圖差異**: 所有範例都顯示了明顯的持續圖差異，證明了 TDA 方法的有效性

詳細結果請參見 `protein-flop.tex` 中的實驗結果部分。

## 未來發展

- [x] 實現阿蒂亞 Flop 的計算範例
- [x] 實現蛋白質 Flop 類比
- [x] 驗證 TDA 能夠捕捉 Flop 操作的拓撲效果
- [ ] 整合真實的 MD 模擬數據
- [ ] 實現更複雜的距離度量（如 Gromov-Wasserstein 距離）
- [ ] 建立典範因子 $K_X$ 與自由能 $G$ 之間的定量對應關係
- [ ] 添加機器學習模型來預測構象轉變
- [ ] 支援更多蛋白質結構類型
- [ ] 擴展到更複雜的蛋白質系統（如打結結構）

## 相關文件

- **理論背景**: `protein-flop.md` - 詳細的理論討論和對話記錄
- **學術論文**: `protein-flop.tex` / `protein-flop.pdf` - 完整的學術論文，包含實驗結果
- **範例說明**: `examples/README.md` - Flip/Flop 範例的詳細說明

## 參考文獻

本專案基於以下理論背景：
- 拓撲數據分析 (Topological Data Analysis)
- 持續同調 (Persistent Homology)
- 極小模型綱領 (Minimal Model Program, MMP)
- 阿蒂亞 Flop (Atiyah Flop)
- 蛋白質折疊機制

詳細理論討論和參考文獻請參見 `protein-flop.md` 和 `protein-flop.tex`。

## 授權

本專案僅供學術研究使用。

## 聯絡資訊

如有問題或建議，請開啟 Issue 或提交 Pull Request。

