"""Gene-family-aware cross-validation.

Tubulin paralogs (or all sarcomere / all MMR genes) must never appear in both
train and test - otherwise the model memorises families instead of learning the
mechanism -> penetrance relationship. Every fold therefore keeps whole gene
families on one side of the split.
"""

from __future__ import annotations

from typing import Iterator, List, Sequence, Tuple

import numpy as np


class GeneFamilyKFold:
    """K-fold splitter that keeps each gene family intact within a fold.

    Families are greedily assigned to the currently-smallest fold (a balanced
    bin-packing), which keeps fold sizes even even when families differ wildly in
    size. Deterministic given ``shuffle`` / ``random_state``.
    """

    def __init__(self, n_splits: int = 5, shuffle: bool = True, random_state: int = 0):
        if n_splits < 2:
            raise ValueError("n_splits must be >= 2")
        self.n_splits = n_splits
        self.shuffle = shuffle
        self.random_state = random_state

    def split(
        self, X: Sequence, y: Sequence, groups: Sequence
    ) -> Iterator[Tuple[np.ndarray, np.ndarray]]:
        groups = np.asarray(groups)
        n = len(groups)
        unique, counts = np.unique(groups, return_counts=True)
        order = np.argsort(-counts)  # largest families first
        unique = unique[order]
        if self.shuffle:
            rng = np.random.default_rng(self.random_state)
            # Shuffle within equal-count blocks to break ties reproducibly.
            perm = rng.permutation(len(unique))
            unique = unique[perm]
            counts = counts[order][perm]
        else:
            counts = counts[order]

        fold_of_family = {}
        fold_sizes = np.zeros(self.n_splits, dtype=int)
        for fam, cnt in zip(unique, counts):
            target = int(np.argmin(fold_sizes))
            fold_of_family[fam] = target
            fold_sizes[target] += cnt

        fold_assignment = np.array([fold_of_family[g] for g in groups])
        indices = np.arange(n)
        for k in range(self.n_splits):
            test_mask = fold_assignment == k
            if not test_mask.any():
                continue
            yield indices[~test_mask], indices[test_mask]

    def get_n_splits(self, X=None, y=None, groups=None) -> int:
        return self.n_splits


def family_aware_split(
    groups: Sequence, test_families: Sequence[str]
) -> Tuple[np.ndarray, np.ndarray]:
    """Single train/test split holding out the named families entirely."""

    groups = np.asarray(groups)
    test_families = set(test_families)
    test_mask = np.array([g in test_families for g in groups])
    indices = np.arange(len(groups))
    return indices[~test_mask], indices[test_mask]
