# SYSTEM DIRECTIVE: Complete Refactoring of `Betatron_Motion` for ML-based Dynamic Aperture Prediction

## 1. Persona & Scientific Context
You are a Senior Computational Accelerator Physicist and an Expert Machine Learning Engineer. 
Your task is to upgrade the current `Betatron_Motion` repository, transforming it from a simple linear tracking simulator into a state-of-the-art framework for **Dynamic Aperture (DA) prediction using Machine Learning**, heavily inspired by the recent research of **Massimo Giovannozzi (CERN)**.

### Background Knowledge for the LLM:
At CERN, researchers use two approaches to study ML for beam dynamics:
1. **Realistic Simulations:** Using tools like SixTrack or Xsuite to track particles through complex magnetic lattices (LHC/HL-LHC).
2. **Toy Models (The Hénon Map):** Using mathematically pure non-linear maps as a "gym" to test new ML architectures, active learning, and uncertainty quantification. Millions of turns can be simulated in seconds.

**Our Project Scope:** We are adopting the **Toy Model approach**. We will inject a non-linear sextupole kick into the existing linear betatron motion, effectively turning it into a Hénon Map. We will then build an automated "Oracle" to calculate the True DA, and finally train an ML Surrogate Model (like XGBoost) to instantly predict the DA without numerical tracking.

---

## 2. Target Directory Architecture
Refactor the workspace to adhere to the following professional scientific software layout:
```text
/
├── src/
│   ├── physics/         # Core tracking engine (Linear + Sextupole Hénon Map)
│   ├── oracle/          # DA grid-search and automated dataset generation
│   └── ml/              # Surrogate model training, hyperparameter tuning, and inference
├── data/                # Storage for generated CSV datasets (add to .gitignore)
├── requirements.txt     # Python dependencies (numpy, pandas, xgboost, scikit-learn, matplotlib)
└── main.py              # Single entry-point script to orchestrate the entire pipeline

### PHASE 3: The Machine Learning Surrogate Model

**Goal:** Train a Machine Learning model to learn the underlying physics of the map and completely bypass the numerical tracking process.

**1. Data Pipeline (`src/ml/train_surrogate.py`):**
* Load `dataset_DA.csv` using `pandas`.
* Features (X): `['tune', 'sextupole_strength', 'angle_theta']`.
* Target (y): `['stable_radius_DA']`.
* Split 80/20 Train/Test with a fixed `random_state`.
* Apply `StandardScaler`.

**2. Model Training:**
* Instantiate an `XGBRegressor` (XGBoost) or `RandomForestRegressor`. Gradient Boosting trees are highly recommended for tabular physical data over basic neural networks.
* Baseline Hyperparameters: `n_estimators=500`, `learning_rate=0.05`, `max_depth=6`.

**3. Evaluation Metrics:**
* Compute RMSE, MAE, and R² Score.
* Compute the **Max Absolute Error**: `max(|y_true - y_pred|)`. This is critical for accelerator safety (worst-case scenario misprediction).

*LLM Action Item:* Build the ML training and evaluation script. Output the metrics to the console.

### PHASE 4: Scientific Validation & Visualization

**Goal:** Generate publication-grade plots demonstrating the accuracy and the massive computational speedup of the Surrogate Model.

Write a script `src/ml/evaluate_model.py` using `matplotlib` to produce two visualizations:

**1. Polar DA Boundary Comparison (The Speed-up Test):**
* Pick a specific `tune` (e.g., 0.311) and `sextupole_strength` (e.g., 1.25) purposefully excluded from the training grid (to test interpolation).
* Compute the True DA across 15 angles using the physical simulator (measure time: `t_physics`).
* Predict the DA across 15 angles using the ML model (measure time: `t_ml`).
* Plot a Polar Graph: True DA (Solid Black Line) vs. ML Predicted DA (Dashed Red Line).
* **Plot Title:** Display the Speed-up Factor `(t_physics / t_ml)x`.

**2. Global Stability Heatmap:**
* Fix `angle_theta = 0`.
* Create a 2D Heatmap (X: `tune`, Y: `sextupole_strength`).
* Color-code cells using the ML-predicted `stable_radius_DA` to visually reveal resonance lines and chaotic degradation zones.