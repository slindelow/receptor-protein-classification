from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Iterable

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_INDEX = {a: i for i, a in enumerate(AA)}
AA_SET = set(AA)

PROTEIN_TYPES = ("aac", "dpc", "aac_dpc", "esm2")
ESM2_MODEL_NAME = "esm2_t6_8M_UR50D"
ESM2_EMBED_DIM = 320  # esm2_t6_8M_UR50D
ESM2_MAX_LEN = 1022  # model context; longer sequences truncated (N-term)


def smiles_to_ecfp(smiles: str, radius: int = 2, n_bits: int = 2048) -> np.ndarray | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    arr = np.zeros((n_bits,), dtype=np.float32)
    AllChem.DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def sequence_to_aac(seq: str) -> np.ndarray:
    vec = np.zeros(20, dtype=np.float32)
    if not seq:
        return vec
    n = 0
    for ch in seq.upper():
        idx = AA_INDEX.get(ch)
        if idx is not None:
            vec[idx] += 1.0
            n += 1
    if n:
        vec /= float(n)
    return vec


def sequence_to_dpc(seq: str) -> np.ndarray:
    """Dipeptide composition: 20x20 normalized counts of consecutive AA pairs."""
    vec = np.zeros(400, dtype=np.float32)
    if not seq:
        return vec
    s = "".join(ch for ch in seq.upper() if ch in AA_SET)
    if len(s) < 2:
        return vec
    n = 0
    for i in range(len(s) - 1):
        a = AA_INDEX[s[i]]
        b = AA_INDEX[s[i + 1]]
        vec[a * 20 + b] += 1.0
        n += 1
    if n:
        vec /= float(n)
    return vec


def sequence_to_aac_dpc(seq: str) -> np.ndarray:
    return np.concatenate([sequence_to_aac(seq), sequence_to_dpc(seq)])


def protein_feat_dim(protein_type: str) -> int:
    t = (protein_type or "aac").lower()
    if t == "aac":
        return 20
    if t == "dpc":
        return 400
    if t == "aac_dpc":
        return 420
    if t == "esm2":
        return ESM2_EMBED_DIM
    raise ValueError(f"Unknown protein_features.type={protein_type!r}; expected one of {PROTEIN_TYPES}")


def seq_hash(seq: str) -> str:
    return hashlib.sha256(seq.encode("utf-8")).hexdigest()


class Esm2Embedder:
    """Frozen ESM-2 mean-pooled embeddings, disk-cached by sequence hash."""

    def __init__(self, cache_dir: Path, model_name: str = ESM2_MODEL_NAME):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.model_name = model_name
        self._model = None
        self._alphabet = None
        self._batch_converter = None
        self._device = "cpu"

    def _ensure_model(self) -> None:
        if self._model is not None:
            return
        import torch
        import esm

        self._model, self._alphabet = esm.pretrained.load_model_and_alphabet(self.model_name)
        self._model.eval()
        self._model = self._model.to(self._device)
        self._batch_converter = self._alphabet.get_batch_converter()

    def _cache_path(self, h: str) -> Path:
        return self.cache_dir / f"{h}.npy"

    def embed_one(self, seq: str) -> np.ndarray:
        s = (seq or "").upper()
        # keep standard AAs only for ESM alphabet compatibility
        s = "".join(ch if ch in AA_SET else "X" for ch in s)
        if not s:
            return np.zeros(ESM2_EMBED_DIM, dtype=np.float32)
        if len(s) > ESM2_MAX_LEN:
            s = s[:ESM2_MAX_LEN]
        h = seq_hash(s)
        path = self._cache_path(h)
        if path.exists():
            return np.load(path).astype(np.float32, copy=False)
        self._ensure_model()
        import torch

        _, _, tokens = self._batch_converter([("prot", s)])
        tokens = tokens.to(self._device)
        with torch.no_grad():
            out = self._model(tokens, repr_layers=[6], return_contacts=False)
        # mean pool residue tokens (exclude BOS/EOS)
        token_reps = out["representations"][6]  # (1, L, C)
        # tokens: 0=BOS, 1..len=residues, len+1=EOS
        reps = token_reps[0, 1 : len(s) + 1].mean(dim=0).cpu().numpy().astype(np.float32)
        if reps.shape[0] != ESM2_EMBED_DIM:
            raise RuntimeError(f"Unexpected ESM-2 dim {reps.shape[0]}, expected {ESM2_EMBED_DIM}")
        np.save(path, reps)
        return reps

    def embed_many(self, sequences: Iterable[str]) -> dict[str, np.ndarray]:
        """Embed unique sequences; return map seq -> vector (original seq key)."""
        uniq = {}
        for seq in sequences:
            key = seq or ""
            if key not in uniq:
                uniq[key] = self.embed_one(key)
        return uniq


def try_init_esm2(cache_dir: Path) -> tuple[Esm2Embedder | None, str | None]:
    """Return (embedder, None) on success, or (None, cut_reason) on failure."""
    try:
        import torch  # noqa: F401
        import esm  # noqa: F401
    except Exception as e:
        return None, f"import failed: {type(e).__name__}: {e}"
    try:
        emb = Esm2Embedder(cache_dir=cache_dir)
        # smoke: tiny sequence to force weight download
        _ = emb.embed_one("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGL")
        return emb, None
    except Exception as e:
        return None, f"init/download failed: {type(e).__name__}: {e}"


def sequence_to_protein_features(
    seq: str,
    protein_type: str,
    esm_embedder: Esm2Embedder | None = None,
    esm_cache: dict[str, np.ndarray] | None = None,
) -> np.ndarray:
    t = (protein_type or "aac").lower()
    if t == "aac":
        return sequence_to_aac(seq)
    if t == "dpc":
        return sequence_to_dpc(seq)
    if t == "aac_dpc":
        return sequence_to_aac_dpc(seq)
    if t == "esm2":
        if esm_cache is not None and (seq or "") in esm_cache:
            return esm_cache[seq or ""]
        if esm_embedder is None:
            raise RuntimeError("esm2 protein features requested but no Esm2Embedder provided")
        return esm_embedder.embed_one(seq or "")
    raise ValueError(f"Unknown protein_features.type={protein_type!r}")


def murcko_scaffold(smiles: str) -> str:
    from rdkit.Chem.Scaffolds import MurckoScaffold

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ""
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
    except Exception:
        return ""


def build_pair_matrix(
    smiles_list: Iterable[str],
    seq_list: Iterable[str] | None,
    radius: int,
    n_bits: int,
    include_protein: bool,
    protein_type: str = "aac",
    esm_embedder: Esm2Embedder | None = None,
    esm_cache: dict[str, np.ndarray] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (X, valid_mask). Invalid SMILES rows are False in mask."""
    smiles_list = list(smiles_list)
    n = len(smiles_list)
    pdim = protein_feat_dim(protein_type) if include_protein else 0
    feat_dim = n_bits + pdim
    X = np.zeros((n, feat_dim), dtype=np.float32)
    mask = np.ones(n, dtype=bool)
    seqs = list(seq_list) if seq_list is not None else [None] * n
    for i, (smi, seq) in enumerate(zip(smiles_list, seqs)):
        fp = smiles_to_ecfp(smi, radius=radius, n_bits=n_bits)
        if fp is None:
            mask[i] = False
            continue
        if include_protein:
            prot = sequence_to_protein_features(
                seq or "", protein_type, esm_embedder=esm_embedder, esm_cache=esm_cache
            )
            X[i] = np.concatenate([fp, prot])
        else:
            X[i] = fp
    return X, mask
