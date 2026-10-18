"""
pdxe_utils.py — shared helpers for the PDXE growth-regression project.

WHY THIS FILE EXISTS
--------------------
In the notebook, `safe_slope`, the BAR formulas and the feature extractors were
each defined inside several different cells. That has three consequences:

  1. Two cells can drift into computing the same quantity differently. That is
     exactly what happened with BAR: one cell capped it at zero, another used a
     trapezoid rule, and a third used the authors' column.
  2. A kernel restart loses every definition, so nothing can be re-run partway.
  3. Nothing can be unit-tested, because there is no importable surface.

Importing from one module fixes all three. Every notebook cell should call these
rather than redefine them.

USAGE
-----
    import pdxe_utils as U
    U.run_self_test()          # verifies the maths on synthetic trajectories

CONVENTIONS
-----------
  * Volume is in mm3. All likelihood-scale work happens on log volume.
  * `day` is days from treatment start (T0), taken from the raw file as-is.
  * An "arm" is one (model_id, treatment) pair. Under the 1x1x1 design that is
    usually one animal; a handful of control arms carry three, and those are
    collapsed on the log scale by `collapse_duplicate_days`.
"""

from __future__ import annotations

import hashlib
import json
import pickle
import re
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Constants — every magic number in the project lives here, once.
# ---------------------------------------------------------------------------

SEED = 20260905

VOLUME_FLOOR = 1.0
VOL_FLOOR = VOLUME_FLOOR

"""mm3. Volumes at or below zero are complete responses. They are FLOORED and
FLAGGED, never dropped: deleting a zero mid-series fabricates a nadir and then
fabricates regrowth after it. 2,265 observations across 135 arms are affected."""

MIN_BAR_DAY = 10.0
MIN_DAY_FOR_BEST_AVERAGE_RESPONSE = MIN_BAR_DAY

"""Gao's minimum evaluation day for best-average-response. Day 0 must not be
eligible to be the argmin: percent-change is identically zero there, which caps
the statistic at zero and makes it degenerate for more than half of all arms."""

HORIZON_DAYS = 28

"""Default horizon for translating a growth-rate difference into a fold
difference in volume, which is the units a reader can act on."""

EARLY_WINDOWS = (7, 10, 14, 21, 28)

CONTROL_PATTERN = r"untreat|vehicle|control|^none$|^no.?treat|mock|dmso"

# Durability-phenotype thresholds (the KT0 definitions).
REG_FRAC = 0.10          # nadir must be >=10% below baseline
MIN_POST_NADIR = 2       # >=2 observations after the nadir
REGROWTH_THRESH = 1.20   # >=20% above nadir


# ---------------------------------------------------------------------------
# Normalisation
# ---------------------------------------------------------------------------

def norm_trt(s) -> str:
    """Normalise a treatment label.

    MUST be applied identically to the observation table and to the metrics
    sheet, or the join silently drops arms. In the raw file `figitumumab"`
    carries a stray double quote, and combination labels have inconsistent
    internal whitespace.
    """
    s = str(s).strip().replace('"', "").replace("'", "")
    return re.sub(r"\s+", " ", s)

normalise_treatment = norm_trt

def is_control(series: pd.Series) -> pd.Series:
    """Boolean mask for control/untreated arms."""
    return series.astype(str).str.contains(CONTROL_PATTERN, case=False, regex=True)


# ---------------------------------------------------------------------------
# Trajectory primitives
# ---------------------------------------------------------------------------

def safe_slope(x, y) -> float:
    """Least-squares slope of y on x. NaN when fewer than two distinct x."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 2 or np.unique(x).size < 2:
        return np.nan
    xm, ym = x.mean(), y.mean()
    denom = ((x - xm) ** 2).sum()
    if denom == 0:
        return np.nan
    return float(((x - xm) * (y - ym)).sum() / denom)


def bar_time_weighted(t, v, split_day=np.inf, min_day=MIN_BAR_DAY) -> float:
    """Best average response over observations with day <= split_day.

    Percent change from baseline is averaged over ELAPSED TIME by the trapezoid
    rule, not over measurement index, so irregular sampling does not silently
    reweight the statistic. The minimum is taken only over timepoints at or
    beyond `min_day`.

    Returns NaN when the window contains no eligible timepoint.

    Validated against Gao's own BestAvgResponse column at Spearman +0.985 on
    4,512 matched arms, so this is a faithful reimplementation.
    """
    t = np.asarray(t, dtype=float)
    v = np.asarray(v, dtype=float)
    keep = t <= split_day
    t, v = t[keep], v[keep]
    if len(t) < 2 or not np.isfinite(v[0]) or v[0] <= 0:
        return np.nan
    resp = 100.0 * (v / v[0] - 1.0)
    dt = np.diff(t)
    if np.any(dt <= 0):
        return np.nan
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (resp[:-1] + resp[1:]) * dt)])
    elig = t >= min_day
    if not elig.any():
        return np.nan
    # An eligible timepoint can coincide with t[0] when an arm's first
    # measurement falls on or after the day-10 eligibility threshold, which
    # makes the denominator zero. Those terms are dropped rather than allowed
    # to produce inf or nan and silently poison the minimum.
    denominator = t[elig] - t[0]
    valid = denominator > 0
    if not valid.any():
        return np.nan
    return float(np.min(cum[elig][valid] / denominator[valid]))

best_average_response = bar_time_weighted

def t2x_within(t, v, split_day=np.inf):
    """Day at which volume first reaches twice baseline.

    Linear interpolation on log volume between the bracketing measurements.
    Restricted to day <= split_day.

    Returns (time, reached). When 2x is never reached inside the window the
    return is (last observed day, False) — a RIGHT-CENSORED value. The caller
    must treat it as censored, not as an observed event time; treating it as
    observed is the informative-censoring error in plan section 3.5.
    """
    t = np.asarray(t, dtype=float)
    v = np.asarray(v, dtype=float)
    keep = t <= split_day
    t, v = t[keep], v[keep]
    if len(t) < 2 or v[0] <= 0:
        return np.nan, False
    y = np.log(v)
    target = np.log(2.0 * v[0])
    hit = np.where(y >= target)[0]
    if len(hit) == 0:
        return float(t[-1]), False
    j = int(hit[0])
    if j == 0:
        return float(t[0]), True
    y0, y1, t0, t1 = y[j - 1], y[j], t[j - 1], t[j]
    if y1 == y0:
        return float(t1), True
    return float(t0 + (target - y0) * (t1 - t0) / (y1 - y0)), True


def logv_at_horizon(t, v, horizon):
    """Log volume at a fixed horizon.

    Interpolates on log volume only when the horizon is BRACKETED by observed
    measurements. Never extrapolates: an arm that ended before the horizon
    returns (nan, False) and is right-censored, not imputed.

    Returns (value, was_exact_observation).
    """
    t = np.asarray(t, dtype=float)
    v = np.asarray(v, dtype=float)
    order = np.argsort(t)
    t, v = t[order], v[order]
    if len(t) == 0 or horizon < t.min() or horizon > t.max():
        return np.nan, False
    y = np.log(np.clip(v, VOL_FLOOR, None))
    exact = np.where(np.isclose(t, horizon))[0]
    if len(exact):
        return float(y[exact[0]]), True
    return float(np.interp(horizon, t, y)), False


# ---------------------------------------------------------------------------
# Per-arm feature extraction
# ---------------------------------------------------------------------------

def window_features(t, v, split_day, prefix="") -> dict:
    """All early-window predictors for one arm, computed from day <= split_day.

    Includes LEVEL features (logV_baseline, logV_last, logV_min). Those matter:
    the forecast outcome is an absolute log volume, so a predictor set built
    only from scale-free quantities is being asked to guess the tumour's size,
    and the comparison then measures who knows the size rather than who knows
    where the trajectory is going.
    """
    t = np.asarray(t, dtype=float)
    v = np.asarray(v, dtype=float)
    keep = t <= split_day
    t, v = t[keep], v[keep]
    p = prefix
    if len(t) < 2:
        return {}

    y = np.log(np.clip(v, VOL_FLOOR, None))
    out = {
        f"{p}n": len(t),
        f"{p}last_day": float(t[-1]),
        # level
        f"{p}logV_baseline": float(y[0]),
        f"{p}logV_last": float(y[-1]),
        f"{p}logV_min": float(y.min()),
        # depth
        f"{p}BAR": bar_time_weighted(t, v, split_day=split_day),
        f"{p}nadir_depth_log": float(y[0] - y.min()),
        f"{p}nadir_ratio": float(v.min() / v[0]),
        # shape
        f"{p}slope": safe_slope(t, y),
    }

    n_init = min(3, len(t))
    out[f"{p}initial_slope"] = safe_slope(t[:n_init], y[:n_init])

    dt, dy = np.diff(t), np.diff(y)
    ok = dt > 0
    out[f"{p}minimum_slope"] = float((dy[ok] / dt[ok]).min()) if ok.any() else np.nan

    nadir_idx = int(np.argmin(y))
    out[f"{p}time_to_nadir"] = float(t[nadir_idx] - t[0])
    if nadir_idx < len(t) - 1 and np.ptp(t[nadir_idx:]) > 0:
        out[f"{p}post_nadir_slope"] = safe_slope(t[nadir_idx:], y[nadir_idx:])
    else:
        out[f"{p}post_nadir_slope"] = np.nan

    out[f"{p}curvature"] = (
        float(np.polyfit(t, y, 2)[0]) if len(t) >= 3 and np.ptp(t) > 0 else np.nan
    )
    a, b = out[f"{p}post_nadir_slope"], out[f"{p}initial_slope"]
    out[f"{p}slope_change"] = a - b if np.isfinite(a) and np.isfinite(b) else np.nan

    tt, reached = t2x_within(t, v, split_day=split_day)
    out[f"{p}t2x"] = tt
    out[f"{p}t2x_reached"] = bool(reached)
    return out


def arm_features(g: pd.DataFrame) -> dict | None:
    """Full-curve features for one arm. `g` holds that arm's observations.

    Returns None for arms with fewer than three usable measurements.
    """
    g = g.sort_values("day").drop_duplicates("day")
    if len(g) < 3:
        return None

    t = g["day"].to_numpy(float)
    v_raw = g["volume"].to_numpy(float)
    ok = np.isfinite(t) & np.isfinite(v_raw)
    t, v_raw = t[ok], v_raw[ok]
    if len(v_raw) < 3:
        return None

    n_floored = int((v_raw < VOL_FLOOR).sum())
    v = np.clip(v_raw, VOL_FLOOR, None)
    baseline = v[0]
    if baseline <= 0:
        return None

    y = np.log(v)
    nadir_idx = int(np.argmin(v))
    nadir, terminal = v[nadir_idx], v[-1]

    post = np.arange(len(v)) > nadir_idx
    if post.sum() >= 2:
        g_late = safe_slope(t[post], y[post])
        yy = y[post]
        pred = np.polyval(np.polyfit(t[post], yy, 1), t[post])
        ss_res = float(((yy - pred) ** 2).sum())
        ss_tot = float(((yy - yy.mean()) ** 2).sum())
        post_r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    else:
        g_late, post_r2 = np.nan, np.nan
    post_n = int(post.sum())

    t2x, t2x_reached = t2x_within(t, v)

    row = {
        "model_id": g["model_id"].iloc[0],
        "indication": g["indication"].iloc[0],
        "treatment": g["treatment"].iloc[0],
        "arm_id": g["arm_id"].iloc[0],
        "n_obs": len(t),
        "duration": float(t[-1] - t[0]),
        "n_floored": n_floored,
        "any_floored": n_floored > 0,
        "baseline": float(baseline),
        "logV_baseline": float(y[0]),
        "terminal_volume": float(terminal),
        "terminal_ratio": float(terminal / baseline),
        "BAR_self": bar_time_weighted(t, v),
        "nadir_ratio": float(nadir / baseline),
        "nadir_day": float(t[nadir_idx]),
        "time_to_nadir": float(t[nadir_idx] - t[0]),
        "nadir_depth_log": float(y[0] - y[nadir_idx]),
        "regrowth_ratio": float(terminal / nadir) if nadir > 0 else np.nan,
        "g_late": g_late,
        "post_n": post_n,
        "post_r2": post_r2,
        "overall_slope": safe_slope(t, y),
        "initial_slope": safe_slope(t[: min(3, len(t))], y[: min(3, len(t))]),
        "t2x_self": t2x,
        "t2x_reached": t2x_reached,
    }

    row["regressed_10"] = row["nadir_ratio"] <= (1 - REG_FRAC)
    row["interior_nadir"] = 0 < nadir_idx < len(v) - 1
    row["two_post"] = post_n >= MIN_POST_NADIR
    row["regrew_120"] = (
        np.isfinite(row["regrowth_ratio"]) and row["regrowth_ratio"] >= REGROWTH_THRESH
    )
    row["durable_estimable"] = (
        row["regressed_10"] and row["interior_nadir"]
        and row["two_post"] and row["regrew_120"]
    )
    """`durable_estimable` selects on the OUTCOME. It is fine for describing how
    often the phenotype exists, and wrong as a primary analysis population —
    conditioning on an interior nadir conditions on a descendant of the
    trajectory. Primary outcomes go through the censoring likelihood instead."""

    for T in EARLY_WINDOWS:
        row.update(window_features(t, v, split_day=T, prefix=f"w{T}_"))
    return row


def collapse_duplicate_days(obs: pd.DataFrame, keys=("arm_id", "day")) -> pd.DataFrame:
    """Collapse multi-animal rows on the LOG scale (geometric mean).

    A handful of control arms carry three animals with no animal identifier, so
    they cannot be split into separate series. The geometric mean is the
    arithmetic mean on log volume, which is the scale the likelihood lives on.
    Adds `n_animals` and `multi_animal_arm` so the choice can be varied in the
    robustness table.
    """
    obs = obs.copy()
    obs["log_volume"] = np.log(obs["volume"].clip(lower=VOL_FLOOR))
    dup_arms = set(obs.loc[obs.duplicated(subset=list(keys), keep=False), "arm_id"])
    group_cols = [c for c in ["model_id", "indication", "treatment", "arm_id", "day"]
                  if c in obs.columns]
    out = (obs.groupby(group_cols, dropna=False)
              .agg(log_volume=("log_volume", "mean"),
                   n_animals=("volume", "size"),
                   vol_raw_min=("volume", "min"))
              .reset_index())
    out["volume"] = np.exp(out["log_volume"])
    out["multi_animal_arm"] = out["arm_id"].isin(dup_arms)
    out["floored"] = out["vol_raw_min"] < VOL_FLOOR
    return out.sort_values(["arm_id", "day"]).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Checkpointing — so a kernel death costs one cell, not two hours
# ---------------------------------------------------------------------------

class Checkpoint:
    """Tiny disk cache keyed by name.

    Every expensive cell should wrap its result:

        samples = ckpt.get_or_compute("m2_samples", lambda: run_fit(...))

    Re-running the cell after a crash loads from disk. Call `ckpt.clear(name)`
    to force a genuine recomputation. For the published pipeline, clear the
    whole directory first so nothing is served from cache — plan section 3.11
    forbids conditional "exists, skipping" behaviour in the final run.
    """

    def __init__(self, root="pdxe_phase0/ckpt"):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)

    def _p(self, name):
        return self.root / f"{name}.pkl"

    def exists(self, name):
        return self._p(name).exists()

    def save(self, name, obj):
        with open(self._p(name), "wb") as f:
            pickle.dump(obj, f, protocol=4)
        mb = self._p(name).stat().st_size / 1e6
        print(f"[ckpt] saved {name} ({mb:.1f} MB)")
        return obj

    def load(self, name):
        with open(self._p(name), "rb") as f:
            obj = pickle.load(f)
        print(f"[ckpt] loaded {name}")
        return obj

    def get_or_compute(self, name, fn, force=False):
        if self.exists(name) and not force:
            return self.load(name)
        return self.save(name, fn())

    def clear(self, name=None):
        targets = [self._p(name)] if name else list(self.root.glob("*.pkl"))
        for p in targets:
            if p.exists():
                p.unlink()
                print(f"[ckpt] cleared {p.name}")


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

def sha256_file(path, chunk=1 << 20) -> str:
    """Streaming SHA256 so a 117 MB workbook never lands in memory."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def content_hash(items, n=12) -> str:
    """Short fingerprint of a sorted list — used for whitelist and CV freezes."""
    joined = "\n".join(sorted(str(x) for x in items))
    return hashlib.sha1(joined.encode()).hexdigest()[:n]


def require_columns(df: pd.DataFrame, cols, where="") -> None:
    """Fail loudly on missing columns.

    Replaces the `[c for c in cols if c in df.columns]` idiom, which silently
    converts a bug into a smaller analysis. That idiom already fired once here:
    a frozen "predictor table" ended up containing one predictor because every
    K0 column had failed to merge and nothing complained.
    """
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise KeyError(f"missing columns {missing}{' in ' + where if where else ''}. "
                       f"Available: {sorted(df.columns)[:40]}")


def grouped_folds(df, group_col="model_id", k=5, seed=SEED) -> pd.Series:
    """Assign CV folds by GROUP, never by row, and assert no leakage.

    One PDX model appears in up to 21 arms, so splitting by arm would put the
    same tumour in train and test.
    """
    rng = np.random.default_rng(seed)
    groups = np.sort(df[group_col].astype(str).unique())
    assign = {g: int(i % k) for i, g in enumerate(rng.permutation(groups))}
    folds = df[group_col].astype(str).map(assign)
    for f in range(k):
        tr = set(df.loc[folds != f, group_col].astype(str))
        te = set(df.loc[folds == f, group_col].astype(str))
        assert tr.isdisjoint(te), f"group leakage in fold {f}"
    return folds


def audit_posterior_memory(samples_by_name: dict, max_trailing=10_000) -> float:
    """Assert no per-observation array survived into a posterior.

    `mu` is (draws x n_obs). At 1,200 draws and ~50k observations that is a
    quarter of a gigabyte per model in float32, doubled by the host copy. That
    is what killed the kernel. Nothing of length ~n_obs belongs in the samples.
    """
    total = 0.0
    for name, s in samples_by_name.items():
        big = {k: tuple(v.shape) for k, v in s.items()
               if getattr(v, "ndim", 1) > 1 and v.shape[-1] > max_trailing}
        if big:
            raise AssertionError(
                f"{name} stores per-observation arrays: {big}. "
                "Remove the matching numpyro.deterministic() call."
            )
        mb = sum(v.nbytes for v in s.values()) / 1e6
        total += mb
        print(f"  {name:<14} {len(s):>2} params  {mb:8.1f} MB")
    print(f"  {'TOTAL':<14} {total:>21.1f} MB")
    return total


# ---------------------------------------------------------------------------
# Self-test — synthetic trajectories with known answers
# ---------------------------------------------------------------------------

def run_self_test(verbose=True) -> bool:
    """Verify the maths on trajectories whose correct answers are known.

    Run this once after importing. It takes under a second and it is the only
    thing standing between a silent formula error and a wrong headline number.
    """
    ok = True

    def check(label, cond):
        nonlocal ok
        ok &= bool(cond)
        if verbose:
            print(f"  [{'PASS' if cond else 'FAIL'}] {label}")

    # --- safe_slope on an exact exponential -------------------------------
    t = np.arange(0, 29, 3.5)
    lam = 0.05
    v = 100 * np.exp(lam * t)
    check("safe_slope recovers exponential rate",
          abs(safe_slope(t, np.log(v)) - lam) < 1e-9)
    check("safe_slope returns NaN on a single point",
          np.isnan(safe_slope([1.0], [2.0])))

    # --- BAR must NOT be capped at zero -----------------------------------
    bar_grow = bar_time_weighted(t, v)
    check("BAR is positive for a purely growing tumour", bar_grow > 0)
    # The old formula returned exactly 0 here because day 0 was eligible.
    check("BAR is not identically zero for a grower", not np.isclose(bar_grow, 0.0))

    v_shrink = 100 * np.exp(-0.05 * t)
    check("BAR is negative for a shrinking tumour",
          bar_time_weighted(t, v_shrink) < 0)
    check("BAR is NaN when the window ends before min_day",
          np.isnan(bar_time_weighted([0., 3., 6.], [100., 90., 80.], split_day=6)))

    # --- time to double ----------------------------------------------------
    tt, reached = t2x_within(np.array([0., 10., 20.]), np.array([100., 100., 400.]))
    check("t2x interpolates on the log scale (expect ~15)",
          reached and abs(tt - 15.0) < 0.5)
    tt2, reached2 = t2x_within(np.array([0., 10.]), np.array([100., 110.]))
    check("t2x reports censored when 2x is never reached",
          (not reached2) and tt2 == 10.0)

    # --- horizon: bracketed only, never extrapolated -----------------------
    yh, exact = logv_at_horizon([0., 14., 28.], [100., 200., 400.], 28)
    check("horizon uses the exact observation when present",
          exact and abs(yh - np.log(400)) < 1e-9)
    yh2, exact2 = logv_at_horizon([0., 14.], [100., 200.], 28)
    check("horizon refuses to extrapolate", np.isnan(yh2) and not exact2)

    # --- arm_features on a regress-then-regrow curve ------------------------
    t3 = np.array([0., 7., 14., 21., 28., 35., 42.])
    v3 = np.array([200., 120., 80., 100., 150., 230., 350.])
    g = pd.DataFrame({"day": t3, "volume": v3,
                      "model_id": "X-TEST", "indication": "TEST",
                      "treatment": "drugA", "arm_id": "X-TEST || drugA"})
    f = arm_features(g)
    check("nadir located at day 14", f["nadir_day"] == 14.0)
    check("interior nadir detected", f["interior_nadir"])
    check("regression flagged (nadir 0.4x baseline)", f["regressed_10"])
    check("regrowth flagged", f["regrew_120"])
    check("durability phenotype estimable", f["durable_estimable"])
    check("g_late is positive after the nadir", f["g_late"] > 0)
    check("post-nadir fit uses 4 points", f["post_n"] == 4)
    check("early-window features present", "w14_logV_last" in f)
    check("14d window sees only days <= 14", f["w14_last_day"] == 14.0)
    check("14d nadir_ratio is 80/200", abs(f["w14_nadir_ratio"] - 0.4) < 1e-9)

    # --- flooring, not dropping -------------------------------------------
    v4 = np.array([200., 50., 0., 0., 60., 180.])
    t4 = np.array([0., 7., 14., 21., 28., 35.])
    g4 = g.iloc[:6].copy()
    g4["day"], g4["volume"] = t4, v4
    f4 = arm_features(g4)
    check("zero volumes floored, not dropped", f4["n_obs"] == 6)
    check("floored observations counted", f4["n_floored"] == 2)

    # --- collapse ----------------------------------------------------------
    dup = pd.DataFrame({
        "model_id": ["X-025"] * 4, "indication": ["GC"] * 4,
        "treatment": ["untreated"] * 4, "arm_id": ["X-025 || untreated"] * 4,
        "day": [0, 0, 0, 5], "volume": [100.0, 200.0, 400.0, 300.0]})
    col = collapse_duplicate_days(dup)
    check("duplicate days collapsed to one row per day", len(col) == 2)
    check("collapse is the geometric mean (100,200,400 -> 200)",
          abs(col.loc[col.day == 0, "volume"].iloc[0] - 200.0) < 1e-6)
    check("multi-animal arm flagged", col["multi_animal_arm"].all())

    # --- grouped folds -----------------------------------------------------
    df = pd.DataFrame({"model_id": np.repeat([f"M{i}" for i in range(20)], 5)})
    folds = grouped_folds(df, k=5)
    check("grouped folds cover every row", folds.notna().all())
    check("grouped folds keep a model in one fold",
          df.assign(f=folds).groupby("model_id").f.nunique().eq(1).all())

    # --- require_columns ---------------------------------------------------
    try:
        require_columns(pd.DataFrame({"a": [1]}), ["a", "b"])
        check("require_columns raises on missing", False)
    except KeyError:
        check("require_columns raises on missing", True)

    print(f"\nSELF-TEST {'PASSED' if ok else 'FAILED'}")
    return ok


if __name__ == "__main__":
    run_self_test()
