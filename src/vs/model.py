from __future__ import annotations

from pathlib import Path
from typing import Any

import joblib
import numpy as np
from sklearn.ensemble import HistGradientBoostingRegressor


def make_regressor(cfg: dict[str, Any]) -> HistGradientBoostingRegressor:
    m = cfg["model"]
    return HistGradientBoostingRegressor(
        max_iter=m["max_iter"],
        learning_rate=m["learning_rate"],
        max_depth=m["max_depth"],
        min_samples_leaf=m["min_samples_leaf"],
        l2_regularization=m["l2_regularization"],
        early_stopping=m["early_stopping"],
        validation_fraction=m["validation_fraction"],
        n_iter_no_change=m["n_iter_no_change"],
        random_state=m["random_state"],
    )


def train_model(X: np.ndarray, y: np.ndarray, cfg: dict[str, Any]) -> HistGradientBoostingRegressor:
    model = make_regressor(cfg)
    model.fit(X, y)
    return model


def save_artifact(path: Path, model, meta: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump({"model": model, "meta": meta}, path)


def load_artifact(path: Path):
    return joblib.load(path)
