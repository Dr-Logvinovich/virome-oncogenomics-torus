#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
HUMAN T2T-CHM13v2.0 24-CHROMOSOME TORUS ATLAS  v1.5.3  —  human_t2t_atlas.py
================================================================================
Autonomous, Zero-RNG Mathematical Engine for Complete Whole-Chromosome Architecture
Across All 24 Human T2T-CHM13v2.0 Chromosomes (NC_060925.1 – NC_060948.1).

Assembly Scope Note:
  T2T-CHM13v2.0 (GCF_009914755.1) combines the 46,XX CHM13 cell line (chr1–chr22,
  chrX) with T2T-chrY from HG002 (46,XY); acrocentric rDNA morph arrays include
  consensus-modeled copies in RefSeq v2.0.

Implemented Theorems (Zero-RNG, Chunked Arithmetic; D15 <50 MB Aux RAM, ~0.4-0.7 GB Peak RSS):
  [THM-D1]  Exact Klein V_4 decomposition of Chargaff PR2 & positional irreversibility.
  [THM-D2]  Exact within-sector heterogeneity null E_het^(d) (W_-+ = W_-- = 0).
  [THM-D3]  Two-sided 16x16 dinucleotide screen (df=3969) & exact KL chain rules.
  [THM-D4]  Blind 3-bp frame-phase lock detector (Omega_phase = Omega_marg + Omega_cond).
  [THM-D5]  Sector scanner with intrinsic right-singular-vector SVD polarity walks.
  [THM-D6]  Exact zero-padded FFT GC autocorrelation (max_lag=600; 2500 on satellites).
  [THM-D7]  Data-driven heteroscedastic Kittler-Illingworth BIC + entropy-inversion partition.
  [THM-D8]  Multi-scale WLS Hurst with Delete-1 Macro-Block Jackknife 95% CI.
  [THM-D9]  Exact U-statistic debiasing of Klein V_4 energies & fixed-point traces.
  [THM-D10] 3.12-Gb 8-strata entropy phase diagram across H_3mer bins.
  [THM-D11] Orthogonal screen/CpG-free anatomy & exact Klein-orbit Top-15 6-mers.
  [THM-D12] 3.12-Gb Whole & Pure-Euchromatin 10-lag spectrum + WW/SS nucleosome pitch.
  [THM-D13] V_4-closed 496-D screen envelope Pi_env & 3600-D strictly closed complement.
  [THM-D14] Exact 4-quadrant noise floor, Lugannani-Rice Bowker chi^2, stride sweep,
            Lemma D14.2 pooling floor, & euchromatin LOO Newey-West HAC 1-Mb inference.
  [THM-D15] Single-pass memory-bounded compositional current null E_loc (whole & euch).
  [THM-D16] Quadratic+3-bp pre-whitened Hann local-band ENBW-adj. F_{2,2M_eff} & WW-SS.
  [THM-D17] Level-2 path signature Levy area with O(n log n) FFT circular-shift null.
  [THM-D18] Exact delayed-pruning LRV-calibrated PELT L2 + cross-channel concordance.
  [THM-D19] Exact mod-p tandem-repeat Klein aliasing (W_-+ = W_-- = 0 for d == 0 mod p).

Authors: Victor Logvinovich (ORCID 0009-0007-2040-6955)
         Jean-Claude Perez  (ORCID 0000-0001-6446-2042)
================================================================================
"""
import argparse
import datetime
import gzip
import hashlib
import json
import math
import os
import platform
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAVE_MPL = True
except Exception:
    HAVE_MPL = False

VERSION = "2.9.1"
SCRIPT_DIR = Path(__file__).resolve().parent
TINY = 1e-300
UA = "human-t2t-atlas/1.0 (academic zero-RNG whole-chromosome audit)"


# =====================================================================
# 0. Safe formatting & atomic I/O
# =====================================================================
def fnum(x, nd=4, sign=False):
    if x is None:
        return "n/a"
    try:
        v = float(x)
    except (TypeError, ValueError):
        return "n/a"
    if not math.isfinite(v):
        return "n/a"
    return f"{v:+.{nd}f}" if sign else f"{v:.{nd}f}"


def _finite(x):
    try:
        return x is not None and math.isfinite(float(x))
    except (TypeError, ValueError):
        return False


def _fval(x, default=np.nan):
    return float(x) if _finite(x) else default


def _get_cmap(name):
    try:
        return matplotlib.colormaps[name]
    except Exception:
        return plt.get_cmap(name)


def _clean(o):
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, (np.floating, float)):
        v = float(o)
        return v if math.isfinite(v) else None
    if isinstance(o, np.ndarray):
        return _clean(o.tolist())
    if isinstance(o, (np.bool_, bool)):
        return bool(o)
    return o


def atomic_write(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)


def sha256_file(path, buf_size=4 * 1024 * 1024):
    p = Path(path)
    if not p.exists():
        return None
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        while True:
            chunk = fh.read(buf_size)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def collect_provenance():
    sp_ver = None
    try:
        import scipy
        sp_ver = getattr(scipy, "__version__", None)
    except Exception:
        pass
    ref_file = ref_sha = None
    for _rel in (
        "human_t2t_euch_ref_v2_9.npz",
        "human_t2t_euch_ref_v2_8.npz",
        "human_t2t_euch_ref_v2_6.npz",
        "atlas_b3_cache.npz",
    ):
        _cand = (SCRIPT_DIR / _rel).resolve()
        if _cand.exists():
            ref_file = _cand.name
            ref_sha = sha256_file(_cand)
            break
    return {
        "version": VERSION,
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "python_version": sys.version.split()[0],
        "numpy_version": getattr(np, "__version__", None),
        "scipy_version": sp_ver,
        "matplotlib_version": getattr(matplotlib, "__version__", None) if HAVE_MPL else None,
        "platform": platform.platform(),
        "script_sha256": sha256_file(__file__),
        "human_ref_npz_file": ref_file,
        "human_ref_npz_sha256": ref_sha,
        "zero_rng": True,
        "partition_method": "log10_inb_raw_kittler_bic_v1.5.3",
        "pelt_penalty": "MBIC_3ln_n_LRV",
    }


_T_CRIT_975 = {
    3: 3.182446, 4: 2.776445, 5: 2.570582, 6: 2.446912, 7: 2.364624,
    8: 2.306004, 9: 2.262157, 10: 2.228139, 11: 2.200985, 12: 2.178813,
    13: 2.160369, 14: 2.144787, 15: 2.131450, 16: 2.119905, 17: 2.109816,
    18: 2.100922, 19: 2.093024, 20: 2.085963,
}


def t_crit_95(df):
    df_i = int(df)
    if df_i in _T_CRIT_975:
        return _T_CRIT_975[df_i]
    try:
        from scipy.stats import t as _sp_t
        return float(_sp_t.ppf(0.975, df_i))
    except Exception:
        z = 1.959963984540054
        d = max(float(df_i), 1.0)
        return z + (z**3 + z) / (4.0 * d) + (5.0 * z**5 + 16.0 * z**3 + 3.0 * z) / (96.0 * d * d)


# =====================================================================
# 1. Discrete Torus T_4^3 = (Z/4Z)^3 Geometry & Klein V_4 Operators
# =====================================================================
BASES = "ACGT"
TRIPLETS = [a + b + c for a in BASES for b in BASES for c in BASES]
TRIPLET_IDX = {c: i for i, c in enumerate(TRIPLETS)}
NC = 64

NT_LUT = np.full(256, 255, dtype=np.uint8)
for _i, _b in enumerate(BASES):
    NT_LUT[ord(_b)] = _i
    NT_LUT[ord(_b.lower())] = _i

_COMP = str.maketrans("ACGT", "TGCA")
MIRROR = np.array([TRIPLET_IDX[c.translate(_COMP)[::-1]] for c in TRIPLETS], dtype=np.int64)

C_N1 = np.arange(NC) // 16
C_N2 = (np.arange(NC) // 4) % 4
C_N3 = np.arange(NC) % 4
SUF16 = 4 * C_N2 + C_N3          # suffix dinucleotide of left triplet
OUT16 = 4 * C_N1 + C_N2          # prefix dinucleotide of right triplet
RHO_BOX = np.array([4 * (3 - (s % 4)) + (3 - (s // 4)) for s in range(16)], dtype=np.int64)


# =====================================================================
# 2. Complete 24-Chromosome Human T2T-CHM13v2.0 Catalog (RefSeq NC_ accessions)
# =====================================================================
_HUMAN24_SPEC = [
    ("chr1",  "NC_060925.1", "metacentric",    "248.4 Mb; largest chr, 1q12 HSat2/3"),
    ("chr2",  "NC_060926.1", "submetacentric", "242.7 Mb; ancestral fusion 2q13"),
    ("chr3",  "NC_060927.1", "metacentric",    "201.1 Mb; metacentric autosome"),
    ("chr4",  "NC_060928.1", "submetacentric", "193.6 Mb; AT-rich euchromatin, compact cen"),
    ("chr5",  "NC_060929.1", "submetacentric", "182.0 Mb; low gene density autosome"),
    ("chr6",  "NC_060930.1", "submetacentric", "172.1 Mb; MHC locus 6p21"),
    ("chr7",  "NC_060931.1", "submetacentric", "160.6 Mb; high segdup density"),
    ("chr8",  "NC_060932.1", "submetacentric", "146.3 Mb; first T2T-resolved autosome"),
    ("chr9",  "NC_060933.1", "submetacentric", "150.6 Mb; massive 9q12 HSat3 block"),
    ("chr10", "NC_060934.1", "submetacentric", "134.8 Mb; submetacentric autosome"),
    ("chr11", "NC_060935.1", "submetacentric", "135.1 Mb; olfactory receptor clusters"),
    ("chr12", "NC_060936.1", "submetacentric", "133.3 Mb; HOXC cluster"),
    ("chr13", "NC_060937.1", "acrocentric",    "113.6 Mb; acrocentric rDNA + AT-rich 13q"),
    ("chr14", "NC_060938.1", "acrocentric",    "101.2 Mb; acrocentric rDNA short arm"),
    ("chr15", "NC_060939.1", "acrocentric",    "99.8 Mb; acrocentric rDNA short arm"),
    ("chr16", "NC_060940.1", "metacentric",    "96.3 Mb; 16q11.2 HSat2 pericentromeric"),
    ("chr17", "NC_060941.1", "submetacentric", "84.3 Mb; gene-dense GC-rich autosome"),
    ("chr18", "NC_060942.1", "submetacentric", "80.5 Mb; gene-poor AT-rich autosome"),
    ("chr19", "NC_060943.1", "metacentric",    "61.7 Mb; highest gene density & GC"),
    ("chr20", "NC_060944.1", "metacentric",    "66.2 Mb; compact metacentric autosome"),
    ("chr21", "NC_060945.1", "acrocentric",    "45.1 Mb; smallest autosome, rDNA + alpha-sat"),
    ("chr22", "NC_060946.1", "acrocentric",    "51.3 Mb; acrocentric rDNA + GC-rich 22q"),
    ("chrX",  "NC_060947.1", "sex",            "154.3 Mb; first T2T chr, LINE-1 rich"),
    ("chrY",  "NC_060948.1", "sex",            "62.5 Mb; palindromes + 35-Mb DYZ1/2 Yq12"),
]
HUMAN_T2T_24 = {
    f"Human_T2T_{c}": {
        "acc": acc, "chr": c, "label": c, "arm_class": arm,
        "default_sector_bp": 100000, "topology": "linear", "domain": "mammal",
        "organism": "Homo sapiens T2T-CHM13v2.0", "notes": note,
    }
    for c, acc, arm, note in _HUMAN24_SPEC
}
PILOT_10_CHROMOSOMES = HUMAN_T2T_24  # backward-compat alias for 24-chr catalog

# =====================================================================
# 3. Chunked Streaming NCBI RefSeq Downloader & FASTA Reader
# =====================================================================
def _ncbi_esummary_nuccore(acc, timeout=60):
    eutils = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    url = eutils + "esummary.fcgi?db=nuccore&retmode=json&id=" + urllib.parse.quote(acc)
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        js = json.loads(r.read().decode("utf-8", "ignore"))
    res = js.get("result", {})
    uids = res.get("uids", [])
    if not uids:
        raise ValueError(f"No nuccore UID found for {acc}")
    rec = res[str(uids[0])]
    return int(rec.get("slen", 0)), str(rec.get("title", acc))


def download_chromosome_fasta(key, meta, cache_dir, chunk_bp=5_000_000, timeout=180):
    """Chunked, timeout-proof NCBI E-utilities downloader saving <key>_<acc>.fna.gz."""
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    acc = meta["acc"]
    dest = cache_dir / f"{key}_{acc}.fna.gz"
    if dest.exists() and dest.stat().st_size > 150:
        return dest

    print(f"  [NCBI] Querying metadata for {key} ({acc})...")
    slen, title = _ncbi_esummary_nuccore(acc, timeout=min(timeout, 60))
    if slen <= 0:
        raise RuntimeError(f"Could not determine sequence length for {acc}")
    n_chunks = int(math.ceil(slen / float(chunk_bp)))
    print(f"  [NCBI] {title} | length={slen:,} bp ({n_chunks} chunk(s))")

    eutils = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    tmp = dest.with_suffix(".part")
    written_bp = 0
    try:
        with gzip.open(tmp, "wt", encoding="ascii", compresslevel=5) as gz:
            gz.write(f">{acc} {title}\n")
            for ci in range(n_chunks):
                s_start = ci * chunk_bp + 1
                s_stop = min(slen, (ci + 1) * chunk_bp)
                url = (f"{eutils}?db=nuccore&id={urllib.parse.quote(acc)}"
                       f"&rettype=fasta&retmode=text&seq_start={s_start}&seq_stop={s_stop}")
                chunk_seq = None
                for attempt in range(1, 5):
                    try:
                        req = urllib.request.Request(url, headers={"User-Agent": UA})
                        with urllib.request.urlopen(req, timeout=timeout) as r:
                            raw_txt = r.read().decode("ascii", "ignore")
                        lines = [ln.strip() for ln in raw_txt.splitlines()
                                 if ln.strip() and not ln.startswith(">")]
                        chunk_seq = "".join(lines)
                        expected = s_stop - s_start + 1
                        if len(chunk_seq) != expected:
                            raise ValueError(f"len mismatch {len(chunk_seq)} != {expected}")
                        break
                    except Exception as e:
                        print(f"    [retry {attempt}/4] {acc} chunk {ci+1}/{n_chunks}: {e}")
                        time.sleep(1.5 * attempt)
                if chunk_seq is None:
                    raise RuntimeError(f"Failed to download {acc} chunk {ci+1}/{n_chunks}")
                for off in range(0, len(chunk_seq), 80):
                    gz.write(chunk_seq[off:off + 80] + "\n")
                written_bp += len(chunk_seq)
                if n_chunks > 1:
                    print(f"    chunk {ci+1}/{n_chunks}: {written_bp:,}/{slen:,} bp "
                          f"({100.0*written_bp/slen:.1f}%)")
                time.sleep(0.35)
        if written_bp != slen:
            raise RuntimeError(f"Total bp mismatch {written_bp} != {slen}")
        os.replace(tmp, dest)
        print(f"  [NCBI] Saved {dest.name} ({dest.stat().st_size/1e6:.2f} MB)")
        return dest
    except Exception:
        if tmp.exists():
            try:
                tmp.unlink()
            except Exception:
                pass
        raise


def load_chromosome_nt(path, buf_bytes=8 * 1024 * 1024):
    """Fast 8-MB block-streaming FASTA reader (~5x faster, low RAM, stops at 2nd record)."""
    path = Path(path)
    opener = gzip.open if str(path).endswith(".gz") else open
    hdr = ""
    chunks = []
    strip_tbl = bytes.maketrans(b"", b"")
    ws_del = b"\r\n \t"
    started = False
    with opener(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                hdr = line[1:].decode("utf-8", "ignore").strip()
                started = True
                break
            s = line.strip()
            if s:
                chunks.append(NT_LUT[np.frombuffer(s, dtype=np.uint8)])
                started = True
                break
        if not started:
            return np.empty(0, dtype=np.uint8), hdr
        tail_nl = True
        while True:
            blk = fh.read(buf_bytes)
            if not blk:
                break
            cut = -1
            if tail_nl and blk.startswith(b">"):
                cut = 0
            else:
                pos = blk.find(b"\n>")
                if pos >= 0:
                    cut = pos
            if cut >= 0:
                blk = blk[:cut]
            cleaned = blk.translate(strip_tbl, ws_del)
            if cleaned:
                chunks.append(NT_LUT[np.frombuffer(cleaned, dtype=np.uint8)])
            if cut >= 0:
                break
            tail_nl = blk.endswith(b"\n")
    if not chunks:
        return np.empty(0, dtype=np.uint8), hdr
    return np.concatenate(chunks), hdr


def reverse_complement_nt(nt):
    """Exact reverse complement on uint8 array (preserving 255 for N)."""
    rc = nt[::-1].copy()
    valid = rc < 4
    rc[valid] = 3 - rc[valid]
    return rc


# =====================================================================
# 4. Exact Chunked Triplet-Pair Counting (Zero Cross-Gap Contamination)
# =====================================================================
def count_triplet_pairs_lag(nt, d=3, chunk_size=4_000_000, phase=None):
    # Counts exact 64x64 triplet pairs (c_i, c_{i+d}) across nt without crossing N-gaps.
    # Fast gap-free path (~3x faster on 100%-valid T2T chromosomes).
    L = int(nt.size)
    B = np.zeros((NC, NC), dtype=np.float64)
    if L <= d + 2:
        return B
    max_i = L - d - 2
    step = max(3, (int(chunk_size) // 3) * 3)
    for s0 in range(0, max_i, step):
        s1 = min(max_i, s0 + step)
        sub_span = nt[s0:s1 + d + 2]
        has_gap = bool(sub_span.max() >= 4) if sub_span.size else False
        if not has_gap:
            if phase is not None:
                sl = slice(s0 + phase, s1, 3)
                sl1 = slice(s0 + phase + 1, s1 + 1, 3)
                sl2 = slice(s0 + phase + 2, s1 + 2, 3)
                rl = slice(s0 + d + phase, s1 + d, 3)
                rl1 = slice(s0 + d + phase + 1, s1 + d + 1, 3)
                rl2 = slice(s0 + d + phase + 2, s1 + d + 2, 3)
                u = (nt[sl].astype(np.int32) << 4) | (nt[sl1].astype(np.int32) << 2) | nt[sl2].astype(np.int32)
                v = (nt[rl].astype(np.int32) << 4) | (nt[rl1].astype(np.int32) << 2) | nt[rl2].astype(np.int32)
            else:
                u = ((nt[s0:s1].astype(np.int32) << 4)
                     | (nt[s0 + 1:s1 + 1].astype(np.int32) << 2)
                     | nt[s0 + 2:s1 + 2].astype(np.int32))
                v = ((nt[s0 + d:s1 + d].astype(np.int32) << 4)
                     | (nt[s0 + d + 1:s1 + d + 1].astype(np.int32) << 2)
                     | nt[s0 + d + 2:s1 + d + 2].astype(np.int32))
            if u.size:
                B += np.bincount((u << 6) | v, minlength=NC * NC).reshape(NC, NC)
            continue
        a0 = nt[s0:s1]
        a1 = nt[s0 + 1:s1 + 1]
        a2 = nt[s0 + 2:s1 + 2]
        b0 = nt[s0 + d:s1 + d]
        b1 = nt[s0 + d + 1:s1 + d + 1]
        b2 = nt[s0 + d + 2:s1 + d + 2]
        ok = (a0 < 4) & (a1 < 4) & (a2 < 4) & (b0 < 4) & (b1 < 4) & (b2 < 4)
        if phase is not None:
            idx_p = np.arange(phase, s1 - s0, 3)
            ok_p = ok[idx_p]
            if not np.any(ok_p):
                continue
            ii = idx_p[ok_p]
            u = (a0[ii].astype(np.int32) << 4) | (a1[ii].astype(np.int32) << 2) | a2[ii].astype(np.int32)
            v = (b0[ii].astype(np.int32) << 4) | (b1[ii].astype(np.int32) << 2) | b2[ii].astype(np.int32)
        else:
            if not np.any(ok):
                continue
            u = (a0[ok].astype(np.int32) << 4) | (a1[ok].astype(np.int32) << 2) | a2[ok].astype(np.int32)
            v = (b0[ok].astype(np.int32) << 4) | (b1[ok].astype(np.int32) << 2) | b2[ok].astype(np.int32)
        B += np.bincount((u << 6) | v, minlength=NC * NC).reshape(NC, NC)
    return B

def count_triplet_pairs_disjoint(nt, d=3, stride=None, chunk_size=4_000_000):
    """[THM-D14] Strictly non-overlapping (disjoint) triplet pairs at given stride (default d+3)."""
    L = int(nt.size)
    B = np.zeros((NC, NC), dtype=np.float64)
    st = int(d) + 3 if stride is None else max(1, int(stride))
    if L <= d + 2:
        return B
    max_i = L - d - 2
    step = max(st, (int(chunk_size) // st) * st)
    for s0 in range(0, max_i, step):
        s1 = min(max_i, s0 + step)
        idx = np.arange(s0, s1, st, dtype=np.int64)
        if idx.size == 0:
            continue
        a0, a1, a2 = nt[idx], nt[idx + 1], nt[idx + 2]
        b0, b1, b2 = nt[idx + d], nt[idx + d + 1], nt[idx + d + 2]
        ok = (a0 < 4) & (a1 < 4) & (a2 < 4) & (b0 < 4) & (b1 < 4) & (b2 < 4)
        if not np.any(ok):
            continue
        u = (a0[ok].astype(np.int32) << 4) | (a1[ok].astype(np.int32) << 2) | a2[ok].astype(np.int32)
        v = (b0[ok].astype(np.int32) << 4) | (b1[ok].astype(np.int32) << 2) | b2[ok].astype(np.int32)
        B += np.bincount((u << 6) | v, minlength=NC * NC).reshape(NC, NC)
    return B


def klein_noise_floor_from_null(E, N):
    """[THM-D14] Exact finite-sample multinomial noise floor E[||W_eps||^2] - ||W_eps(P)||^2."""
    N = float(N)
    se = float(np.sum(E))
    if N <= 0 or se <= 0:
        return {"W_pp": 0.0, "W_pm": 0.0, "W_mp": 0.0, "W_mm": 0.0}
    P = np.asarray(E, dtype=np.float64) / se
    A = float(np.sum(P * P))
    Tp = float(np.sum(P * P.T))
    Rp = float(np.sum(P * P.T[np.ix_(MIRROR, MIRROR)]))
    Sp = float(np.sum(P * P[np.ix_(MIRROR, MIRROR)]))
    d = float(np.trace(P))
    l = float(np.sum(P[np.arange(NC), MIRROR]))
    return {
        "W_pp": (1.0 - A + (d - Tp) + (l - Rp) - Sp) / (4.0 * N),
        "W_pm": (1.0 - A + (d - Tp) - (l - Rp) + Sp) / (4.0 * N),
        "W_mp": (1.0 - A - (d - Tp) + (l - Rp) + Sp) / (4.0 * N),
        "W_mm": (1.0 - A - (d - Tp) - (l - Rp) - Sp) / (4.0 * N),
    }


def klein_energy_sd_from_null(P_null, N):
    # [THM-V18] Exact finite-N standard deviation of degenerate U-statistic energies
    # u_pm, u_mp, u_mm under H_0: P in W_++ (exact factor 1/2 per orthogonal V_4 quadrant):
    #   Var(u_pm) = (p2 + d2 - l2) / (2 N (N - 1))
    #   Var(u_mp) = (p2 - d2 + l2) / (2 N (N - 1))
    #   Var(u_mm) = (p2 - d2 - l2) / (2 N (N - 1))
    N = float(N)
    if N <= 1.0:
        return {"W_pm": 0.0, "W_mp": 0.0, "W_mm": 0.0, "var_W_pm": 0.0, "var_W_mp": 0.0, "var_W_mm": 0.0}
    P = np.asarray(P_null, dtype=np.float64)
    s = float(np.sum(P))
    if s <= 0.0:
        return {"W_pm": 0.0, "W_mp": 0.0, "W_mm": 0.0, "var_W_pm": 0.0, "var_W_mp": 0.0, "var_W_mm": 0.0}
    P = P / s
    P_pp = 0.25 * (P + P.T + P.T[np.ix_(MIRROR, MIRROR)] + P[np.ix_(MIRROR, MIRROR)])
    p2 = float(np.sum(P_pp * P_pp))
    d2 = float(np.sum(np.diag(P_pp) ** 2))
    l2 = float(np.sum(P_pp[np.arange(NC), MIRROR] ** 2))
    denom = 2.0 * N * (N - 1.0)
    v_pm = max(p2 + d2 - l2, 0.0) / denom
    v_mp = max(p2 - d2 + l2, 0.0) / denom
    v_mm = max(p2 - d2 - l2, 0.0) / denom
    return {
        "W_pm": math.sqrt(v_pm),
        "W_mp": math.sqrt(v_mp),
        "W_mm": math.sqrt(v_mm),
        "var_W_pm": v_pm,
        "var_W_mp": v_mp,
        "var_W_mm": v_mm,
    }


def klein_design_effect_audit(nt, d=3):
    # [THM-V18] Zero-RNG serial-correlation / overlap design-effect audit via
    # strictly non-overlapping stride=6 pairs and contiguous Delete-1 Block Jackknife.
    L = int(nt.size)
    B_all = count_triplet_pairs_lag(nt, d=d)
    kv_all = klein_v4_pr2_audit(B_all)
    if kv_all is None:
        return None
    B_disj = count_triplet_pairs_disjoint(nt, d=d, stride=int(d) + 3)
    kv_disj = klein_v4_pr2_audit(B_disj)
    K = min(10, max(5, L // 500))
    se_jk_mp = se_jk_mm = z_jk_mp = z_jk_mm = deff_mp = deff_mm = None
    n_blk_ok = 0
    if L >= K * 60:
        edges = np.linspace(0, L, K + 1, dtype=int)
        blk_mats = [count_triplet_pairs_lag(nt[edges[b]:edges[b + 1]], d=d) for b in range(K)]
        B_sum = sum(blk_mats)
        u_mp_loo, u_mm_loo = [], []
        for b in range(K):
            kv_l = klein_v4_pr2_audit(B_sum - blk_mats[b])
            if kv_l is not None:
                u_mp_loo.append(kv_l["norm2_W_mp_unb"])
                u_mm_loo.append(kv_l["norm2_W_mm_unb"])
        if len(u_mp_loo) >= 4:
            n_blk_ok = len(u_mp_loo)
            am = np.asarray(u_mp_loo, dtype=np.float64)
            ac = np.asarray(u_mm_loo, dtype=np.float64)
            fac = float(n_blk_ok - 1) / float(n_blk_ok)
            se_jk_mp = math.sqrt(max(fac * float(np.sum((am - float(np.mean(am))) ** 2)), 0.0))
            se_jk_mm = math.sqrt(max(fac * float(np.sum((ac - float(np.mean(ac))) ** 2)), 0.0))
            if se_jk_mp > 1e-36:
                z_jk_mp = kv_all["norm2_W_mp_unb"] / se_jk_mp
            if se_jk_mm > 1e-36:
                z_jk_mm = kv_all["norm2_W_mm_unb"] / se_jk_mm
            sd_mp = kv_all.get("sd_iid_W_mp") or 0.0
            sd_mm = kv_all.get("sd_iid_W_mm") or 0.0
            if sd_mp > 1e-36:
                deff_mp = (se_jk_mp / sd_mp) ** 2
            if sd_mm > 1e-36:
                deff_mm = (se_jk_mm / sd_mm) ** 2
    return {
        "N_disj6": kv_disj["N"] if kv_disj else 0,
        "u_mp_disj6": kv_disj["norm2_W_mp_unb"] if kv_disj else None,
        "z_iid_mp_disj6": kv_disj.get("z_iid_W_mp") if kv_disj else None,
        "u_mm_disj6": kv_disj["norm2_W_mm_unb"] if kv_disj else None,
        "z_iid_mm_disj6": kv_disj.get("z_iid_W_mm") if kv_disj else None,
        "n_blocks_jk": n_blk_ok,
        "se_jk_W_mp": se_jk_mp,
        "z_jk_W_mp": z_jk_mp,
        "deff_W_mp": deff_mp,
        "se_jk_W_mm": se_jk_mm,
        "z_jk_W_mm": z_jk_mm,
        "deff_W_mm": deff_mm,
    }


def chi2_tail_log10p(chi2_val, df):
    """Exact Lugannani-Rice saddlepoint log10(P(chi^2_df > chi2_val)) (<=1e-6 dex vs SciPy, never underflows)."""
    if chi2_val <= 0.0 or df <= 0:
        return 0.0
    k = float(df)
    x = float(chi2_val)
    if abs(x - k) < 1e-5 * math.sqrt(2.0 * k):
        p0 = 0.5 - 1.0 / (3.0 * math.sqrt(math.pi * k))
        return math.log10(min(max(p0, 1e-300), 1.0))
    sign = 1.0 if x > k else -1.0
    arg = max(x - k - k * math.log(x / k), 0.0)
    w = sign * math.sqrt(arg)
    u = (x - k) / math.sqrt(2.0 * k)
    if w >= 8.0:
        corr = (1.0 / u) - (1.0 / (w ** 3)) + (3.0 / (w ** 5))
        if corr <= 0.0:
            corr = 1.0 / max(u, 1e-30)
        ln_p = -0.5 * w * w - 0.5 * math.log(2.0 * math.pi) + math.log(corr)
        return ln_p / math.log(10.0)
    phi_w = math.exp(-0.5 * w * w) / math.sqrt(2.0 * math.pi)
    tail = 0.5 * math.erfc(w / math.sqrt(2.0)) + phi_w * (1.0 / u - 1.0 / w)
    return math.log10(min(max(tail, 1e-300), 1.0))


def bowker_klein_chi2_test(B_disj):
    """[THM-D14] Descriptive Bowker symmetry chi^2 (df=2016) on disjoint pairs
    (Rao-Scott clustering adjustment not applied; see 1-Mb cluster inference for p-values)."""
    B = np.asarray(B_disj, dtype=np.float64)
    N = float(B.sum())
    if N <= 0:
        return None
    iu, iv = np.triu_indices(NC, k=1)
    a_t, b_t = B[iu, iv], B[iv, iu]
    s_t = a_t + b_t
    m_t = s_t > 0
    chi2_T = float(np.sum(((a_t[m_t] - b_t[m_t]) ** 2) / s_t[m_t]))
    df_T = int(np.count_nonzero(m_t))

    B_R = B.T[np.ix_(MIRROR, MIRROR)]
    pal_mask = np.zeros((NC, NC), dtype=bool)
    pal_mask[np.arange(NC), MIRROR] = True
    s_r = B + B_R
    m_r = (~pal_mask) & (s_r > 0)
    chi2_R = 0.5 * float(np.sum(((B[m_r] - B_R[m_r]) ** 2) / s_r[m_r]))
    df_R = int(np.count_nonzero(m_r) // 2)
    lp_T = chi2_tail_log10p(chi2_T, max(df_T, 1))
    lp_R = chi2_tail_log10p(chi2_R, max(df_R, 1))
    return {
        "N_disjoint": int(N),
        "bowker_T_chi2": chi2_T,
        "bowker_T_df": df_T,
        "bowker_T_chi2_per_N": chi2_T / max(N, 1.0),
        "bowker_T_log10p_unadj": lp_T,
        "bowker_T_log10p": lp_T,
        "bowker_R_chi2": chi2_R,
        "bowker_R_df": df_R,
        "bowker_R_chi2_per_N": chi2_R / max(N, 1.0),
        "bowker_R_log10p_unadj": lp_R,
        "bowker_R_log10p": lp_R,
        "rao_scott_note": "Descriptive only: disjoint pairs remove overlap but not Mb-scale clustering.",
    }


def bowker_stride_sweep(nt, d=3, strides=(6, 12, 24, 48, 96, 192)):
    """[THM-D14] Disjoint stride sweep diagnosing Rao-Scott design effect via chi^2_s / N_s."""
    out = []
    for st in strides:
        Bd = count_triplet_pairs_disjoint(nt, d=d, stride=st)
        bk = bowker_klein_chi2_test(Bd)
        if bk is not None:
            out.append({
                "stride": int(st), "N": bk["N_disjoint"],
                "chi2_T": bk["bowker_T_chi2"], "chi2_T_per_N": bk["bowker_T_chi2_per_N"],
                "chi2_R": bk["bowker_R_chi2"], "chi2_R_per_N": bk["bowker_R_chi2_per_N"],
            })
    return out


def newey_west_t_stats(scores):
    """Exact Bartlett-kernel Newey-West HAC SE & t-stat (L = floor(K^{1/3}), zero RNG)."""
    x = np.asarray(scores, dtype=np.float64)
    x = x[np.isfinite(x)]
    K = int(x.size)
    if K < 5:
        return None
    mu = float(np.mean(x))
    xc = x - mu
    g0 = float(np.dot(xc, xc)) / float(K)
    var_iid = float(np.dot(xc, xc)) / float(K - 1)
    se_iid = math.sqrt(max(var_iid, 0.0) / float(K))
    L_nw = min(max(1, int(math.floor(K ** (1.0 / 3.0)))), K - 2)
    lrv = g0
    for lag in range(1, L_nw + 1):
        w = 1.0 - float(lag) / float(L_nw + 1)
        gl = float(np.dot(xc[:-lag], xc[lag:])) / float(K)
        lrv += 2.0 * w * gl
    lrv_adj = max(lrv * (float(K) / float(K - 1)), g0 / float(K))
    se_hac = math.sqrt(max(lrv_adj, 0.0) / float(K))
    return {
        "L_nw": int(L_nw), "se_iid": se_iid, "se_hac": se_hac,
        "t_iid": (mu / se_iid) if se_iid > 1e-30 else None,
        "t_hac": (mu / se_hac) if se_hac > 1e-30 else None,
    }


def cluster_1mb_arrow_test(block_mats):
    """[THM-D14 v1.5] Zero-RNG cluster inference across ~1-Mb blocks using
    leave-one-block-out W_-+ axis v_{-b} & Bartlett Newey-West HAC SE."""
    valid_W, valid_N = [], []
    for Bm in block_mats:
        Nb = float(np.sum(Bm))
        if Nb < 1000.0:
            continue
        kv = klein_v4_pr2_audit(Bm)
        if kv is not None:
            valid_W.append(kv["_W_mp"] * Nb)
            valid_N.append(Nb)
    K = len(valid_W)
    if K < 5:
        return None
    W_tot = sum(valid_W)
    scores = np.empty(K, dtype=np.float64)
    for i in range(K):
        W_loo = W_tot - valid_W[i]
        n_loo = float(np.linalg.norm(W_loo, ord="fro"))
        scores[i] = 0.0 if n_loo <= 1e-30 else float(np.sum((valid_W[i] / valid_N[i]) * (W_loo / n_loo)))
    mean_a = float(np.mean(scores))
    nw = newey_west_t_stats(scores)
    n_pos = int(np.count_nonzero(scores > 0.0))
    k_lo = min(n_pos, K - n_pos)
    if K <= 1000:
        num = sum(math.comb(K, j) for j in range(k_lo + 1))
        p_sign = min(1.0, float(2 * num) / float(1 << K))
    else:
        z_s = (abs(n_pos - 0.5 * K) - 0.5) / math.sqrt(0.25 * K)
        p_sign = min(1.0, math.erfc(max(z_s, 0.0) / math.sqrt(2.0)))
    return {
        "n_blocks": K, "n_pos": n_pos, "frac_pos": float(n_pos) / float(K),
        "mean_proj": mean_a, "se_proj_iid": nw["se_iid"], "se_proj_hac": nw["se_hac"],
        "se_proj": nw["se_hac"], "newey_west_L": nw["L_nw"],
        "cluster_t_iid": nw["t_iid"], "cluster_t_hac": nw["t_hac"],
        "cluster_t": nw["t_hac"], "sign_p_twosided": p_sign,
    }


def compositional_current_multiscale(nt, windows=(200, 1000, 10000, 100000), d=3, trip_budget=4_000_000):
    """[THM-D15] Exact single-pass, zero-RNG, strictly memory-bounded (<50 MB RAM at ANY M)
    nested multiscale compositional current null E_loc = E_sym + J."""
    ws = tuple(sorted(int(m) for m in windows if int(m) > int(d)))
    if not ws:
        return {}
    M0 = ws[0]
    assert all(M % M0 == 0 for M in ws), f"Windows {ws} must be multiples of {M0}"
    ratios = [M // M0 for M in ws]
    max_r = max(ratios)
    super_bp = max_r * M0
    if super_bp <= trip_budget:
        chunk_base_wins = max(1, trip_budget // super_bp) * max_r
    else:
        chunk_base_wins = max(1, trip_budget // M0)

    E_syms = {M: np.zeros((NC, NC), dtype=np.float64) for M in ws}
    Js = {M: np.zeros((NC, NC), dtype=np.float64) for M in ws}
    masses = {M: 0.0 for M in ws}
    w_syms = {M: float(M - d) / float(M * (M - 1)) for M in ws}

    L = int(nt.size)
    if L <= d + 2:
        return {M: (E_syms[M], 0.0) for M in ws}

    valid = nt < 4
    padded = np.concatenate([[False], valid, [False]])
    starts = np.flatnonzero(~padded[:-1] & padded[1:])
    ends = np.flatnonzero(padded[:-1] & ~padded[1:])

    for s0, s1 in zip(starts, ends):
        n_trip = int(s1 - s0) - 2
        nb0_tot = n_trip // M0
        if nb0_tot <= 0:
            continue
        p_prev = {M: None for M in ws}
        rem_buf = {M: np.zeros(NC, dtype=np.float64) for M in ws}
        rem_cnt = {M: 0 for M in ws}
        for M, r in zip(ws, ratios):
            nb_m = nb0_tot // r
            if nb_m > 0:
                masses[M] += float(nb_m * M - d)
        for w0 in range(0, nb0_tot, chunk_base_wins):
            nb0 = min(chunk_base_wins, nb0_tot - w0)
            base0 = int(s0) + w0 * M0
            n_u = nb0 * M0
            u = ((nt[base0:base0 + n_u].astype(np.int32) << 4)
                 | (nt[base0 + 1:base0 + n_u + 1].astype(np.int32) << 2)
                 | nt[base0 + 2:base0 + n_u + 2].astype(np.int32))
            offsets = (np.arange(nb0, dtype=np.int32) * NC)[:, None]
            C0 = np.bincount((u.reshape(nb0, M0) + offsets).ravel(),
                             minlength=nb0 * NC).reshape(nb0, NC).astype(np.float64)
            del u, offsets
            for M, r in zip(ws, ratios):
                if r == 1:
                    Cm = C0
                else:
                    rows_out = []
                    pos_c = 0
                    if rem_cnt[M] > 0:
                        need = r - rem_cnt[M]
                        take = min(need, nb0)
                        rem_buf[M] += C0[:take].sum(axis=0)
                        rem_cnt[M] += take
                        pos_c = take
                        if rem_cnt[M] == r:
                            rows_out.append(rem_buf[M].copy())
                            rem_buf[M].fill(0.0)
                            rem_cnt[M] = 0
                    avail = nb0 - pos_c
                    full = avail // r
                    if full > 0:
                        Cm_mid = C0[pos_c:pos_c + full * r].reshape(full, r, NC).sum(axis=1)
                        pos_c += full * r
                    else:
                        Cm_mid = None
                    if pos_c < nb0:
                        rem_buf[M] += C0[pos_c:].sum(axis=0)
                        rem_cnt[M] += (nb0 - pos_c)
                    if rows_out and Cm_mid is not None:
                        Cm = np.vstack([rows_out[0][None, :], Cm_mid])
                    elif rows_out:
                        Cm = rows_out[0][None, :]
                    elif Cm_mid is not None:
                        Cm = Cm_mid
                    else:
                        continue
                E_syms[M] += w_syms[M] * (Cm.T @ Cm - np.diag(Cm.sum(axis=0)))
                Pm = Cm / float(M)
                if p_prev[M] is not None:
                    Js[M] += float(d) * np.outer(p_prev[M], Pm[0])
                if Pm.shape[0] >= 2:
                    Js[M] += float(d) * (Pm[:-1].T @ Pm[1:])
                p_prev[M] = Pm[-1].copy()
    return {M: (E_syms[M] + Js[M], masses[M]) for M in ws}


def compositional_current_null(nt, window_m=1000, d=3, trip_budget=4_000_000, max_win_chunk=None):
    """[THM-D15] Single-scale wrapper backed by memory-bounded compositional_current_multiscale."""
    tb = int(trip_budget)
    if max_win_chunk is not None:
        tb = min(tb, max(int(window_m), int(max_win_chunk) * int(window_m)))
    res = compositional_current_multiscale(nt, windows=(int(window_m),), d=d, trip_budget=tb)
    return res[int(window_m)]


def fisher_g_pvalue(powers):
    """[THM-D16] Exact closed-form Fisher g-test p-value (Decimal precision, zero RNG)."""
    from decimal import Decimal, getcontext
    p = np.asarray(powers, dtype=np.float64)
    p = p[np.isfinite(p) & (p >= 0.0)]
    n = int(p.size)
    s = float(p.sum())
    if n < 3 or s <= 0.0:
        return None, None
    g = float(np.max(p) / s)
    m_max = min(int(math.floor(1.0 / max(g, 1e-12))), n)
    getcontext().prec = 60
    gD = Decimal(str(g))
    total = Decimal(0)
    comb = Decimal(1)
    for k in range(1, m_max + 1):
        comb = comb * Decimal(n - k + 1) / Decimal(k)
        term = comb * ((Decimal(1) - Decimal(k) * gD) ** (n - 1))
        total = (total + term) if (k % 2 == 1) else (total - term)
    return g, max(0.0, min(1.0, float(total)))


def bartlett_acf_se(acf, n_eff, lag_k=10):
    """[THM-D16] Bartlett asymptotic SE for sample autocorrelation r_k."""
    acf = np.asarray(acf, dtype=np.float64)
    if n_eff <= 10 or acf.size < lag_k:
        return None
    r = np.concatenate([[1.0], acf])
    K = min(acf.size, 60)
    rk = float(r[lag_k])
    s = 0.0
    for v in range(-K, K + 1):
        rv = float(r[abs(v)])
        r_vk = float(r[abs(v + lag_k)]) if abs(v + lag_k) <= acf.size else 0.0
        r_vm = float(r[abs(v - lag_k)]) if abs(v - lag_k) <= acf.size else 0.0
        s += rv * rv + r_vk * r_vm - 4.0 * rk * rv * r_vm + 2.0 * (rk * rk) * (rv * rv)
    return math.sqrt(max(s, 0.0) / float(n_eff))


def _detrend_and_remove_3bp(lags, y):
    """Exact 5-param OLS removal of quadratic baseline (kills curved isochore ACF decay)
    plus 3-bp harmonic (cos/sin 2pi*d/3) before Hann-tapered spectral estimation."""
    t_norm = (lags - float(np.mean(lags))) / max(float(np.std(lags)), 1e-12)
    X = np.column_stack([
        np.ones_like(lags),
        t_norm,
        t_norm * t_norm,
        np.cos(2.0 * math.pi * lags / 3.0),
        np.sin(2.0 * math.pi * lags / 3.0),
    ])
    beta, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ beta


def _hann_taper(y):
    """Unit-RMS Hann window (ENBW = 1.5 bins) suppressing smooth baseline sidelobes."""
    w = np.hanning(int(y.size))
    return y * (w / math.sqrt(max(float(np.mean(w * w)), 1e-30)))


def notched_fisher_and_band_peak_p(acf, p_peak_nt, spec_lim=200, m_eff=14):
    """[THM-D16 v1.5.1] Quadratic+3-bp OLS pre-whitened, Hann-tapered spectral test:
    calibrates peak power against same-window (4..200) local in-band PSD using the
    ENBW-adjusted (Welch-Satterthwaite, Hann ENBW=1.5) F_{2, 2*M_eff} approximation
    (1 + r/M_eff)^(-M_eff) with conservative grid Sidak factor m_eff=14 (>=3x the
    ~3-4 Rayleigh Fourier bins in [9.5, 11.7] nt). Notched Fisher g is exploratory
    under Hann tapering (adjacent bins correlated)."""
    acf = np.asarray(acf, dtype=np.float64)
    lim = min(int(spec_lim), acf.size)
    if lim < 25:
        return None, None, None, None, None
    lags = np.arange(4, lim + 1, dtype=np.float64)
    yw = _hann_taper(_detrend_and_remove_3bp(lags, acf[3:lim]))
    n_y = yw.size
    freqs = np.fft.rfftfreq(n_y, d=1.0)
    fft_p = (np.abs(np.fft.rfft(yw)) ** 2) / float(n_y)
    pwr_clean = float(np.abs(np.sum(yw * np.exp(-2j * math.pi * lags / float(p_peak_nt)))) ** 2) / float(n_y)

    keep = np.ones(freqs.size, dtype=bool)
    keep[:min(4, freqs.size)] = False
    for h_per in (3.0, 6.0, 9.0, 12.0, 15.0, 18.0):
        keep &= (np.abs(freqs - (1.0 / h_per)) > 0.015)
    fg_notch, fp_notch = fisher_g_pvalue(fft_p[keep]) if int(np.count_nonzero(keep)) >= 5 else (None, None)

    f_pk = 1.0 / max(float(p_peak_nt), 1.0)
    band_loc = (
        (freqs >= 0.05) & (freqs <= 0.22)
        & (np.abs(freqs - f_pk) > 0.018)
        & (np.abs(freqs - 1.0 / 6.0) > 0.010)
        & (np.abs(freqs - 1.0 / 9.0) > 0.010)
        & (np.abs(freqs - 1.0 / 12.0) > 0.010)
        & (np.abs(freqs - 1.0 / 15.0) > 0.008)
        & (np.abs(freqs - 1.0 / 18.0) > 0.008)
    )
    m_raw = int(np.count_nonzero(band_loc))
    if m_raw < 5:
        band_loc = keep & (np.abs(freqs - f_pk) > 0.018)
        m_raw = int(np.count_nonzero(band_loc))
    psd_noise = float(np.mean(fft_p[band_loc])) if m_raw >= 3 else None

    p_raw = p_adj = None
    if psd_noise is not None and psd_noise > 1e-30:
        M_eff = max(3.0, float(m_raw) / 1.5)
        ratio = pwr_clean / psd_noise
        ln_p = -M_eff * math.log1p(ratio / M_eff)
        p_raw = math.exp(max(ln_p, -700.0))
        if p_raw < 1e-8:
            p_adj = min(1.0, float(m_eff) * p_raw)
        else:
            p_adj = min(1.0, max(0.0, -math.expm1(float(m_eff) * math.log1p(-min(p_raw, 1.0 - 1e-15)))))
    return fg_notch, fp_notch, psd_noise, p_raw, p_adj


def otsu_threshold_1d(vals, n_bins=None):
    """Exact bin-free O(n log n) zero-RNG Otsu inter-class variance threshold (no histogram binning)."""
    a = np.asarray(vals, dtype=np.float64)
    a = np.sort(a[np.isfinite(a)])
    n = int(a.size)
    if n < 10 or float(a[-1] - a[0]) < 1e-9:
        return None
    jump = a[1:] > (a[:-1] + 1e-15)
    if not np.any(jump):
        return None
    csum = np.cumsum(a)
    tot = float(csum[-1])
    k = np.arange(1, n, dtype=np.float64)
    w0 = k / float(n)
    w1 = 1.0 - w0
    mu0 = csum[:-1] / k
    mu1 = (tot - csum[:-1]) / (float(n) - k)
    var_b = np.where(jump, w0 * w1 * ((mu0 - mu1) ** 2), -1.0)
    best_i = int(np.argmax(var_b))
    return float(0.5 * (a[best_i] + a[best_i + 1]))


def bic_2vs1_gaussian(vals, thr):
    """Exact 2-component hard-split vs 1-component Gaussian BIC difference (BIC1 - BIC2 > 0 favors 2)."""
    a = np.asarray(vals, dtype=np.float64)
    a = a[np.isfinite(a)]
    n = int(a.size)
    if n < 10:
        return -1e9
    v1 = max(float(np.var(a, ddof=0)), 1e-18)
    ll1 = -0.5 * n * (math.log(2.0 * math.pi * v1) + 1.0)
    bic1 = -2.0 * ll1 + 2.0 * math.log(float(n))
    c0, c1 = a[a <= thr], a[a > thr]
    if c0.size < 3 or c1.size < 3:
        return -1e9
    v0_raw = float(np.var(c0, ddof=0))
    v2_raw = float(np.var(c1, ddof=0))
    if v0_raw < 1e-10 or v2_raw < 1e-10:
        return -1e9
    v0 = max(v0_raw, 1e-18)
    v2 = max(v2_raw, 1e-18)
    ll2 = (
        c0.size * math.log(c0.size / float(n)) - 0.5 * c0.size * (math.log(2.0 * math.pi * v0) + 1.0)
        + c1.size * math.log(c1.size / float(n)) - 0.5 * c1.size * (math.log(2.0 * math.pi * v2) + 1.0)
    )
    bic2 = -2.0 * ll2 + 5.0 * math.log(float(n))
    return float(bic1 - bic2)


def _kittler_bic_constrained(feat, h3, min_k=3, min_dbic=10.0):
    """Exact O(n log n) zero-RNG heteroscedastic 2-component Gaussian profile-likelihood
    threshold (Kittler & Illingworth 1986, Pattern Recognition) constrained to the
    entropy-inversion half-space (mean & median H_3 strictly lower in high-I class),
    gated by Kass-Raftery (1995) very-strong-evidence delta_BIC > 10.0."""
    order = np.argsort(feat, kind="mergesort")
    a = feat[order]
    h = h3[order]
    n = int(a.size)
    if n < 2 * min_k + 4 or float(a[-1] - a[0]) < 1e-9:
        return None, -1e9, 0, 0, False
    v1 = max(float(np.var(a, ddof=0)), 1e-18)
    ll1 = -0.5 * n * (math.log(2.0 * math.pi * v1) + 1.0)
    bic1 = -2.0 * ll1 + 2.0 * math.log(float(n))
    cs1 = np.cumsum(a)
    cs2 = np.cumsum(a * a)
    tot1, tot2 = float(cs1[-1]), float(cs2[-1])
    ch1 = np.cumsum(h)
    tot_h = float(ch1[-1])
    k = np.arange(min_k, n - min_k + 1, dtype=np.int64)
    k = k[a[k] > (a[k - 1] + 1e-15)]
    if k.size == 0:
        return None, -1e9, 0, 0, False
    n0 = k.astype(np.float64)
    n1 = float(n) - n0
    s1_0 = cs1[k - 1]
    s2_0 = cs2[k - 1]
    s1_1 = tot1 - s1_0
    s2_1 = tot2 - s2_0
    v0 = s2_0 / n0 - (s1_0 / n0) ** 2
    v2 = s2_1 / n1 - (s1_1 / n1) ** 2
    v_floor = max(1e-10, 1e-4 * v1)
    mean_h_lo = ch1[k - 1] / n0
    mean_h_hi = (tot_h - ch1[k - 1]) / n1
    ok = (v0 >= v_floor) & (v2 >= v_floor) & (mean_h_hi < mean_h_lo)
    if not np.any(ok):
        ok_u = (v0 >= v_floor) & (v2 >= v_floor)
        if not np.any(ok_u):
            return None, -1e9, 0, 0, False
        ll2_u = (
            n0[ok_u] * np.log(n0[ok_u] / float(n)) - 0.5 * n0[ok_u] * (np.log(2.0 * math.pi * v0[ok_u]) + 1.0)
            + n1[ok_u] * np.log(n1[ok_u] / float(n)) - 0.5 * n1[ok_u] * (np.log(2.0 * math.pi * v2[ok_u]) + 1.0)
        )
        db_u = float(bic1 - (-2.0 * float(np.max(ll2_u)) + 5.0 * math.log(float(n))))
        return None, db_u, 0, 0, False
    kf, n0f, n1f, v0f, v2f = k[ok], n0[ok], n1[ok], v0[ok], v2[ok]
    ll2 = (
        n0f * np.log(n0f / float(n)) - 0.5 * n0f * (np.log(2.0 * math.pi * v0f) + 1.0)
        + n1f * np.log(n1f / float(n)) - 0.5 * n1f * (np.log(2.0 * math.pi * v2f) + 1.0)
    )
    for ri in np.argsort(ll2)[::-1]:
        ki = int(kf[ri])
        db = float(bic1 - (-2.0 * float(ll2[ri]) + 5.0 * math.log(float(n))))
        if db <= min_dbic:
            return None, db, int(n0f[ri]), int(n1f[ri]), False
        if float(np.median(h[ki:])) < float(np.median(h[:ki])):
            return float(0.5 * (a[ki - 1] + a[ki])), db, int(n0f[ri]), int(n1f[ri]), True
    return None, -1e9, 0, 0, False


def partition_euch_sat_parameter_free(H3, INB, INB_RAW=None):
    """[THM-D7 v1.5.3] Data-driven Euchromatin vs Satellite/Tandem partition (no hand-tuned per-chr cutoffs):
      1) Exact Shannon half-entropy low-complexity guard H_3 < 3.0 bits (= 0.5 * log2(64),
         i.e. <8 effective 3-mers) flags ultra-low-complexity tracts where bridge MI saturates.
      2) On H_3 >= 3.0 sectors, exact heteroscedastic Kittler-Illingworth (1986) 2-component
         Gaussian BIC argmax on log10(I_nonbr_raw) constrained to the entropy-inversion
         half-space (mean & median H_3 lower in high-I class, delta_BIC > 10.0), resolving
         both compact ~1% centromeres (chr2, chr3) and large 10-62% blocks (chr1, chr9, chrY)
         without equal-variance Otsu imbalance bias.
    """
    H3 = np.asarray(H3, dtype=np.float64)
    INB = np.asarray(INB, dtype=np.float64)
    use_raw = INB_RAW is not None
    SRC = np.asarray(INB_RAW, dtype=np.float64) if use_raw else INB
    is_valid = np.isfinite(H3) & np.isfinite(SRC)
    is_micro_sat = is_valid & (H3 < 3.0)
    rem = is_valid & (~is_micro_sat)
    n_rem = int(np.count_nonzero(rem))
    sat_high_mem = np.zeros_like(is_valid, dtype=bool)
    inb_thr = None
    delta_bic = None
    inv_ok = False
    frac_clipped = 0.0
    n_hi = n_lo = 0
    if n_rem >= 12:
        inb_r = SRC[rem]
        h3_r = H3[rem]
        if use_raw:
            floor_v = 1e-12
            frac_clipped = float(np.mean(inb_r <= floor_v))
        else:
            pos_inb = inb_r[inb_r > 0.0]
            floor_v = max(1e-4, 0.05 * float(np.median(pos_inb))) if pos_inb.size else 1e-4
            frac_clipped = float(np.mean(inb_r <= floor_v))
        feat = np.log10(np.maximum(inb_r, floor_v))
        thr_f, delta_bic, n_lo, n_hi, inv_ok = _kittler_bic_constrained(feat, h3_r, min_k=3, min_dbic=10.0)
        if thr_f is not None and inv_ok:
            hi = feat > thr_f
            lo = ~hi
            if (not use_raw) and frac_clipped > 0.15 and float(np.mean(inb_r[lo] <= floor_v)) > 0.50:
                delta_bic = -1e9
                inv_ok = False
            else:
                inb_thr = float(10.0 ** thr_f)
                idx_rem = np.flatnonzero(rem)
                sat_high_mem[idx_rem[hi]] = True
    is_sat = is_micro_sat | sat_high_mem
    is_euch = is_valid & (~is_sat)
    return is_euch, is_sat, {
        "partition_method": "log10_inb_raw_kittler_bic_v1.5.3" if use_raw else "log10_inb_kittler_bic_v1.5.3",
        "half_entropy_guard_bits": 3.0,
        "n_micro_sat_h3_lt_3": int(np.count_nonzero(is_micro_sat)),
        "inb_threshold_otsu": inb_thr,
        "bic_2vs1_delta": delta_bic,
        "entropy_inversion_ok": inv_ok,
        "frac_clipped": frac_clipped,
        "n_hi": n_hi,
        "n_lo": n_lo,
        "used_raw_inb": use_raw,
    }


def levy_area_signed_test(a, b):
    """[THM-D17 & D17.2] Exact Level-2 path signature Levy area, Stratonovich shuffle check,
    and EXACT O(n log n) FFT exhaustive circular-shift null preserving each channel's Hurst ACF.
    Note: FFT circular-correlation equals exact shifted Levy area STRICTLY because increments
    ac, bc are mean-centered (dX = dY = 0 => shuffle s2_xy + s2_yx = 0 => A = s2_xy)."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    n = int(a.size)
    if n < 6:
        return None
    ac = a - float(np.mean(a))
    bc = b - float(np.mean(b))
    ca_prev = np.concatenate([[0.0], np.cumsum(ac[:-1])])
    cb_prev = np.concatenate([[0.0], np.cumsum(bc[:-1])])
    s2_xy = float(np.sum((ca_prev + 0.5 * ac) * bc))
    s2_yx = float(np.sum((cb_prev + 0.5 * bc) * ac))
    dX, dY = float(np.sum(ac)), float(np.sum(bc))
    shuffle_err = abs((s2_xy + s2_yx) - dX * dY)
    A_levy = 0.5 * (s2_xy - s2_yx)

    w_a = ca_prev + 0.5 * ac
    A_shifts = np.fft.irfft(np.fft.rfft(w_a) * np.conj(np.fft.rfft(bc)), n=n)
    sd_shift = float(np.std(A_shifts[1:], ddof=1)) if n > 2 else 0.0
    z_shift = (A_levy / sd_shift) if sd_shift > 1e-30 else None
    cnt_ge = int(np.count_nonzero(np.abs(A_shifts[1:]) >= abs(A_levy) - 1e-12))
    p_shift = float(1 + cnt_ge) / float(n)

    cov = np.cov(ac, bc, ddof=1)
    det_sig = max(float(cov[0, 0] * cov[1, 1] - cov[0, 1] * cov[1, 0]), 0.0)
    var_A = (n * (n - 1) / 4.0) * det_sig
    z_iid = (A_levy / math.sqrt(var_A)) if var_A > 1e-30 else None
    return {
        "n_steps": n,
        "levy_area": float(A_levy),
        "levy_p_shift": p_shift,
        "levy_z_shift": float(z_shift) if z_shift is not None else None,
        "abs_levy_z_shift": abs(float(z_shift)) if z_shift is not None else None,
        "levy_z_iid_unadj": float(z_iid) if z_iid is not None else None,
        "abs_levy_z": abs(float(z_shift)) if z_shift is not None else None,
        "shuffle_err": float(shuffle_err),
    }


def _top_svd_intrinsic(X):
    """Exact deterministic top singular pair (u1, var_exp, anchor_idx) using partial SciPy eigh
    if available, with 100% INTRINSIC sign anchored on max-|v1| in torus-orbit space R^{4096}.
    Canonical BLAS-independent row-major tie-break (u <= v across symmetric/antisymmetric
    (u,v)/(v,u) pairs via 1e-12 relative tolerance); strictly invariant to sector permutation."""
    try:
        from scipy.linalg import eigh as sp_eigh
        have_sp = True
    except Exception:
        have_sp = False
    n_r, n_c = X.shape
    tot_ev = max(float(np.sum(X * X)), TINY)
    if n_r <= n_c:
        G = X @ X.T
        if have_sp and n_r > 32:
            ev, evec = sp_eigh(G, subset_by_index=(n_r - 1, n_r - 1))
            top_ev = max(float(ev[-1]), 0.0)
            u_unit = evec[:, -1]
        else:
            ev, evec = np.linalg.eigh(G)
            top_ev = max(float(ev[-1]), 0.0)
            u_unit = evec[:, -1]
        s1 = math.sqrt(top_ev)
        u1 = u_unit * s1
        v1 = (X.T @ u_unit) / max(s1, TINY)
    else:
        C = X.T @ X
        if have_sp and n_c > 32:
            ev, evec = sp_eigh(C, subset_by_index=(n_c - 1, n_c - 1))
            top_ev = max(float(ev[-1]), 0.0)
            v1 = evec[:, -1]
        else:
            ev, evec = np.linalg.eigh(C)
            top_ev = max(float(ev[-1]), 0.0)
            v1 = evec[:, -1]
        u1 = X @ v1
    abs_v = np.abs(v1)
    max_abs = float(abs_v.max())
    k_max = int(np.flatnonzero(abs_v >= max_abs * (1.0 - 1e-12))[0]) if max_abs > 0.0 else 0
    if float(v1[k_max]) < 0.0:
        u1 = -u1
    return u1, float(top_ev / tot_ev), k_max


def count_triplets_1d(nt, chunk_size=4_000_000):
    L = int(nt.size)
    c = np.zeros(NC, dtype=np.float64)
    if L < 3:
        return c
    max_i = L - 2
    for s0 in range(0, max_i, chunk_size):
        s1 = min(max_i, s0 + chunk_size)
        sub = nt[s0:s1 + 2]
        if sub.size and int(sub.max()) < 4:
            u = ((nt[s0:s1].astype(np.int32) << 4)
                 | (nt[s0 + 1:s1 + 1].astype(np.int32) << 2)
                 | nt[s0 + 2:s1 + 2].astype(np.int32))
            c += np.bincount(u, minlength=NC)
            continue
        if sub.size and int(sub.min()) >= 4:
            continue
        a0, a1, a2 = nt[s0:s1], nt[s0 + 1:s1 + 1], nt[s0 + 2:s1 + 2]
        ok = (a0 < 4) & (a1 < 4) & (a2 < 4)
        if np.any(ok):
            u = (a0[ok].astype(np.int32) << 4) | (a1[ok].astype(np.int32) << 2) | a2[ok].astype(np.int32)
            c += np.bincount(u, minlength=NC)
    return c


def sector_het_expected_counts(nt, sector_bp, d=3):
    L = int(nt.size)
    n_sec = int(math.ceil(L / float(sector_bp)))
    C = np.zeros((n_sec, NC), dtype=np.float64)
    w_vec = np.zeros(n_sec, dtype=np.float64)
    B_within = np.zeros((NC, NC), dtype=np.float64)
    for w in range(n_sec):
        seg = nt[w * sector_bp:min(L, (w + 1) * sector_bp)]
        if seg.size <= d + 2 or int(seg.min()) >= 4:
            continue
        cw = count_triplets_1d(seg, chunk_size=sector_bp)
        Mw = float(cw.sum())
        if Mw < 2.0:
            continue
        Bw = count_triplet_pairs_lag(seg, d=d, chunk_size=sector_bp)
        Nw = float(Bw.sum())
        if Nw <= 0.0:
            continue
        C[w] = cw
        w_vec[w] = Nw / (Mw * (Mw - 1.0))
        B_within += Bw
    E_het = np.maximum(C.T @ (w_vec[:, None] * C) - np.diag(C.T @ w_vec), 0.0)
    return B_within, E_het

# =====================================================================
# 5. Information-Theoretic & Klein V_4 Primitives (Theorems D1, D3, D4)
# =====================================================================
def _mi_tab(obs, den, N):
    m = obs > 0
    r = np.where(m, obs / np.maximum(den, TINY), 1.0)
    return float(np.sum(np.where(m, (obs / N) * np.log(np.maximum(r, TINY)), 0.0)))


def _kl_joint(P, Q):
    m = P > 0
    return float(np.sum(P[m] * np.log(P[m] / np.maximum(Q[m], TINY))))


def klein_v4_pr2_audit(B):
    """[THM-D1] Exact Klein V_4 decomposition of a 64x64 pair table B into
    static symmetric (W_++), static PR2-defect (W_+-), PR2-respecting time arrow (W_-+),
    and chiral PR2-breaking flow (W_--)."""
    N = float(np.sum(B))
    if N <= 0:
        return None
    M = np.asarray(B, dtype=np.float64) / N
    TP = M.T
    RP = M.T[np.ix_(MIRROR, MIRROR)]
    SP = M[np.ix_(MIRROR, MIRROR)]
    W_pp = 0.25 * (M + TP + RP + SP)
    W_pm = 0.25 * (M + TP - RP - SP)
    W_mp = 0.25 * (M - TP + RP - SP)
    W_mm = 0.25 * (M - TP - RP + SP)
    n_tot = float(np.sum(M * M))
    n_pp = float(np.sum(W_pp * W_pp))
    n_pm = float(np.sum(W_pm * W_pm))
    n_mp = float(np.sum(W_mp * W_mp))
    n_mm = float(np.sum(W_mm * W_mm))
    e_pr2 = 0.25 * float(np.sum((M - RP) ** 2))
    e_time = 0.25 * float(np.sum((M - TP) ** 2))
    a_pr2 = 0.5 * float(np.sum(np.abs(M - RP)))
    a_time = 0.5 * float(np.sum(np.abs(M - TP)))
    pr2_sum = n_pm + n_mm
    time_sum = n_mp + n_mm
    p_diag = float(np.trace(M))
    p_pal = float(np.sum(M[np.arange(NC), MIRROR]))
    if N > 1.0:
        fI = (N * n_tot - 1.0) / (N - 1.0)
        fT = (N * float(np.sum(M * TP)) - p_diag) / (N - 1.0)
        fR = (N * float(np.sum(M * RP)) - p_pal) / (N - 1.0)
        fS = (N * float(np.sum(M * SP))) / (N - 1.0)
        u_pp = 0.25 * (fI + fT + fR + fS)
        u_pm = 0.25 * (fI + fT - fR - fS)
        u_mp = 0.25 * (fI - fT + fR - fS)
        u_mm = 0.25 * (fI - fT - fR + fS)
    else:
        u_pp, u_pm, u_mp, u_mm = n_pp, n_pm, n_mp, n_mm
    u_time = max(u_mp, 0.0) + max(u_mm, 0.0)
    sd_null = klein_energy_sd_from_null(W_pp, N)
    sd_pm, sd_mp, sd_mm = sd_null["W_pm"], sd_null["W_mp"], sd_null["W_mm"]
    z_pm = (float(u_pm) / sd_pm) if sd_pm > 1e-36 else None
    z_mp = (float(u_mp) / sd_mp) if sd_mp > 1e-36 else None
    z_mm = (float(u_mm) / sd_mm) if sd_mm > 1e-36 else None
    u_mm_upper95 = max(float(u_mm), 0.0) + 1.6448536269514722 * sd_mm
    u_time_u95 = max(float(u_mp), 0.0) + u_mm_upper95
    pr2_lower95 = (max(float(u_mp), 0.0) / u_time_u95) if u_time_u95 > 1e-30 else None
    return {
        "N": int(N),
        "norm2_total": n_tot,
        "norm2_W_pp": n_pp,
        "norm2_W_pm": n_pm,
        "norm2_W_mp": n_mp,
        "norm2_W_mm": n_mm,
        "norm2_W_pp_unb": float(u_pp),
        "norm2_W_pm_unb": float(u_pm),
        "norm2_W_mp_unb": float(u_mp),
        "norm2_W_mm_unb": float(u_mm),
        "sd_iid_W_pm": float(sd_pm),
        "sd_iid_W_mp": float(sd_mp),
        "sd_iid_W_mm": float(sd_mm),
        "z_iid_W_pm": float(z_pm) if z_pm is not None else None,
        "z_iid_W_mp": float(z_mp) if z_mp is not None else None,
        "z_iid_W_mm": float(z_mm) if z_mm is not None else None,
        "W_mm_unb_censored": bool(u_mm <= 0.0),
        "W_mm_95pct_upper": float(u_mm_upper95),
        "time_PR2_compliant_share_unb_lower95": float(pr2_lower95) if pr2_lower95 is not None else None,
        "p_diag_homodimer": p_diag,
        "p_pal_6bp": p_pal,
        "PR2_defect_F2": pr2_sum,
        "time_arrow_F2": time_sum,
        "PR2_dyn_share": (n_mm / pr2_sum) if pr2_sum > 1e-30 else None,
        "time_PR2_compliant_share": (n_mp / time_sum) if time_sum > 1e-30 else None,
        "time_PR2_compliant_share_unb": (max(u_mp, 0.0) / u_time) if u_time > 1e-30 else None,
        "A_PR2_L1": a_pr2,
        "A_time_L1": a_time,
        "pythagorean_err": abs(n_tot - (n_pp + n_pm + n_mp + n_mm)),
        "PR2_split_err": abs(e_pr2 - pr2_sum),
        "time_split_err": abs(e_time - time_sum),
        "_W_pm": W_pm,
        "_W_mp": W_mp,
        "_W_mm": W_mm,
    }


def strip_private(d):
    return {k: v for k, v in d.items() if not str(k).startswith("_")}


def two_sided_screen_full(B):
    """[THM-D3] Four-term 16x16 dinucleotide screen split on the full 64x64 torus
    (df = 225 + 720 + 720 + 2304 = 3969 = 63^2)."""
    N = float(np.sum(B))
    if N <= 0:
        return None
    xs = np.arange(NC)
    nx, ny = B.sum(axis=1), B.sum(axis=0)
    Bso = np.bincount((SUF16[:, None] * 16 + OUT16[None, :]).ravel(),
                      weights=B.ravel(), minlength=256).reshape(16, 16)
    Bxo = np.bincount((xs[:, None] * 16 + OUT16[None, :]).ravel(),
                      weights=B.ravel(), minlength=NC * 16).reshape(NC, 16)
    Bsy = np.bincount((SUF16[:, None] * NC + xs[None, :]).ravel(),
                      weights=B.ravel(), minlength=16 * NC).reshape(16, NC)
    ns, no = Bso.sum(axis=1), Bso.sum(axis=0)
    I_cc = _mi_tab(B, np.maximum(nx[:, None] * ny[None, :], TINY) / N, N)
    I_br = _mi_tab(Bso, np.maximum(ns[:, None] * no[None, :], TINY) / N, N)
    D_lead = _mi_tab(Bxo, np.maximum(nx[:, None] * Bso[SUF16, :], TINY) / np.maximum(ns[SUF16][:, None], TINY), N)
    D_wob = _mi_tab(Bsy, np.maximum(Bso[:, OUT16] * ny[None, :], TINY) / np.maximum(no[OUT16][None, :], TINY), N)
    D_13 = _mi_tab(B, np.maximum(Bxo[:, OUT16] * Bsy[SUF16, :], TINY) / np.maximum(Bso[np.ix_(SUF16, OUT16)], TINY), N)

    df_br = max(0, int(np.count_nonzero(ns > 0)) - 1) * max(0, int(np.count_nonzero(no > 0)) - 1)
    df_lead = sum(
        max(0, int(np.count_nonzero(Bxo[SUF16 == s_, :].sum(axis=1) > 0)) - 1)
        * max(0, int(np.count_nonzero(Bso[s_, :] > 0)) - 1)
        for s_ in range(16) if ns[s_] > 0
    )
    df_wob = sum(
        max(0, int(np.count_nonzero(Bso[:, o_] > 0)) - 1)
        * max(0, int(np.count_nonzero(ny[OUT16 == o_] > 0)) - 1)
        for o_ in range(16) if no[o_] > 0
    )
    df_13 = sum(
        max(0, int(np.count_nonzero(Bxo[SUF16 == s_, o_] > 0)) - 1)
        * max(0, int(np.count_nonzero(Bsy[s_, OUT16 == o_] > 0)) - 1)
        for s_ in range(16) for o_ in range(16) if Bso[s_, o_] > 0
    )
    df_nonbr = int(df_lead + df_wob + df_13)
    I_nonbr = D_lead + D_wob + D_13
    return {
        "N": int(N),
        "I_CC": I_cc,
        "I_bridge": I_br,
        "D_5prime_lead": D_lead,
        "D_3prime_wobble": D_wob,
        "D_1_3": D_13,
        "I_nonbridge": I_nonbr,
        "I_nonbridge_MM": I_nonbr - df_nonbr / (2.0 * N),
        "lead_wobble_PR2_asym": abs(D_lead - D_wob),
        "four_term_err": abs(I_cc - (I_br + D_lead + D_wob + D_13)),
        "df_eff": {"bridge": df_br, "lead": df_lead, "wobble": df_wob, "d13": df_13,
                   "nonbridge": df_nonbr, "sum": df_br + df_nonbr},
    }


def screen_kl_chain_rule(B, E, df_eff=None):
    """[THM-D3] Exact forward & reverse 3-term screen-stratified KL chain rules:
        D_KL(P || Q) = K_bridge + K_lead|SO + K_wob|XO
                     = K_bridge + K_wob|SO  + K_lead|SY."""
    N1, N2 = float(np.sum(B)), float(np.sum(E))
    if N1 <= 0 or N2 <= 0:
        return None
    P = np.asarray(B, dtype=np.float64) / N1
    Q = np.maximum(np.asarray(E, dtype=np.float64) / N2, TINY)
    Q = Q / Q.sum()
    xs = np.arange(NC)

    def _projs(M):
        Mso = np.bincount((SUF16[:, None] * 16 + OUT16[None, :]).ravel(),
                          weights=M.ravel(), minlength=256).reshape(16, 16)
        Mxo = np.bincount((xs[:, None] * 16 + OUT16[None, :]).ravel(),
                          weights=M.ravel(), minlength=NC * 16).reshape(NC, 16)
        Msy = np.bincount((SUF16[:, None] * NC + xs[None, :]).ravel(),
                          weights=M.ravel(), minlength=16 * NC).reshape(16, NC)
        return Mso, Mxo, Msy

    Pso, Pxo, Psy = _projs(P)
    Qso, Qxo, Qsy = _projs(Q)
    kl_tot = _kl_joint(P, Q)
    K_br = _kl_joint(Pso, np.maximum(Qso, TINY))

    m_xo = Pxo > 0
    den_l = np.where(m_xo, Qxo * Pso[SUF16, :] / np.maximum(Qso[SUF16, :], TINY), TINY)
    K_lead_SO = float(np.sum(Pxo[m_xo] * np.log(Pxo[m_xo] / np.maximum(den_l[m_xo], TINY))))

    m_xy = P > 0
    den_wx = np.where(m_xy, Q * Pxo[:, OUT16] / np.maximum(Qxo[:, OUT16], TINY), TINY)
    K_wob_XO = float(np.sum(P[m_xy] * np.log(P[m_xy] / np.maximum(den_wx[m_xy], TINY))))

    m_sy = Psy > 0
    den_w = np.where(m_sy, Qsy * Pso[:, OUT16] / np.maximum(Qso[:, OUT16], TINY), TINY)
    K_wob_SO = float(np.sum(Psy[m_sy] * np.log(Psy[m_sy] / np.maximum(den_w[m_sy], TINY))))

    den_ly = np.where(m_xy, Q * Psy[SUF16, :] / np.maximum(Qsy[SUF16, :], TINY), TINY)
    K_lead_SY = float(np.sum(P[m_xy] * np.log(P[m_xy] / np.maximum(den_ly[m_xy], TINY))))

    K_nonbr = kl_tot - K_br
    out = {
        "KL_total": kl_tot,
        "K_bridge": K_br,
        "K_lead_given_SO": K_lead_SO,
        "K_wob_given_XO": K_wob_XO,
        "K_wob_given_SO": K_wob_SO,
        "K_lead_given_SY": K_lead_SY,
        "K_nonbridge": K_nonbr,
        "nonbridge_share": (K_nonbr / kl_tot) if kl_tot > 1e-15 else None,
        "fwd_err": abs(kl_tot - (K_br + K_lead_SO + K_wob_XO)),
        "rev_err": abs(kl_tot - (K_br + K_wob_SO + K_lead_SY)),
        "all_nonnegative": bool(min(K_br, K_lead_SO, K_wob_XO, K_wob_SO, K_lead_SY) > -1e-12),
    }
    if df_eff is not None:
        df_nb = int(df_eff["nonbridge"])
        df_lsy = int(df_eff["lead"] + df_eff["d13"])
        out["K_nonbridge_MM"] = K_nonbr - df_nb / (2.0 * N1)
        out["K_lead_given_SY_MM"] = K_lead_SY - df_lsy / (2.0 * N1)
    return out


def frame_phase_lock_sector(seg, d=3):
    """[THM-D4] Exact 3-frame Jensen-Shannon divergence Omega_phase = Omega_marg + Omega_cond
    on one chromosome sector (blind executable/3n-periodicity detector).
    Aligns all 3 frames to an identical pair span (divisible by 3) so boundary counts match."""
    n_avail = int(seg.size) - d - 2
    if n_avail < 300:
        return None
    n_aligned = (n_avail // 3) * 3
    seg_al = seg[:n_aligned + d + 2]
    Bf = [count_triplet_pairs_lag(seg_al, d=d, chunk_size=seg_al.size + 6, phase=ph) for ph in (0, 1, 2)]
    ns = [float(b.sum()) for b in Bf]
    if min(ns) < 100.0:
        return None
    Pf = [b / n for b, n in zip(Bf, ns)]
    Pbar = (Pf[0] + Pf[1] + Pf[2]) / 3.0
    om_phase = sum(_kl_joint(P, Pbar) for P in Pf) / 3.0
    pX = [P.sum(axis=1) for P in Pf]
    pXbar = (pX[0] + pX[1] + pX[2]) / 3.0
    om_marg = sum(_kl_joint(px, pXbar) for px in pX) / 3.0
    om_cond = 0.0
    for P, px in zip(Pf, pX):
        for x in np.flatnonzero(px > 0):
            om_cond += px[x] * _kl_joint(P[x] / px[x], Pbar[x] / max(pXbar[x], TINY))
    om_cond /= 3.0
    Nmean = float(np.mean(ns))
    Bu = Bf[0] + Bf[1] + Bf[2]
    qu = int(np.count_nonzero(Bu.sum(axis=1) > 0))
    mu = int(np.count_nonzero(Bu > 0))
    df_cond = 2 * max(mu - qu, 0)
    bias_cond = df_cond / (6.0 * max(Nmean, 1.0))
    return {
        "Omega_phase": float(om_phase),
        "Omega_marg": float(om_marg),
        "Omega_cond": float(om_cond),
        "Omega_cond_excess": float(om_cond - bias_cond),
        "bias_cond": float(bias_cond),
        "chain_err": abs(om_phase - (om_marg + om_cond)),
    }


# =====================================================================
# 6. Exact FFT Nucleotide Autocorrelation (1..200 bp) & Helical Pitch
# =====================================================================
def exact_gc_autocorr_fft(nt, max_lag=600, chunk_size=2_000_000):
    """[THM-D6/D16] Exact zero-padded FFT GC/AT autocorrelation over d=1..max_lag bp
    across valid ACGT positions (never crosses N-gaps; O(L log chunk_size)).
    Note: acf_values stores the first min(max_lag, 250) lags for compact JSON;
    higher-order satellite HOR harmonics (342/513/2052 bp) are extracted before truncation."""
    valid = nt < 4
    n_valid = int(np.count_nonzero(valid))
    if n_valid <= max_lag + 10:
        max_lag = min(max_lag, n_valid - 15)
    if max_lag < 20:
        return None
    is_gc = (nt == 1) | (nt == 2)
    mu = float(np.count_nonzero(is_gc)) / float(n_valid)
    var = mu * (1.0 - mu)
    if var < 1e-12:
        return None

    L = int(nt.size)
    num = np.zeros(max_lag + 1, dtype=np.float64)
    den = np.zeros(max_lag + 1, dtype=np.float64)

    step = int(chunk_size)
    pos = 0
    while pos < L:
        end = min(L, pos + step)
        ext_end = min(L, end + max_lag)
        seg = nt[pos:ext_end]
        if seg.size == 0 or int(seg.min()) >= 4:
            pos = end
            continue
        m_ext = (seg < 4).astype(np.float64)
        x_ext = np.where(seg < 4, ((seg == 1) | (seg == 2)).astype(np.float64) - mu, 0.0)
        n_left = end - pos
        m_left = m_ext.copy()
        x_left = x_ext.copy()
        m_left[n_left:] = 0.0
        x_left[n_left:] = 0.0
        nfft = 1 << int(math.ceil(math.log2(seg.size + max_lag + 2)))
        fx_l = np.fft.rfft(x_left, nfft)
        fx_e = np.fft.rfft(x_ext, nfft)
        corr_x = np.fft.irfft(np.conj(fx_l) * fx_e, nfft)[:max_lag + 1]
        fm_l = np.fft.rfft(m_left, nfft)
        fm_e = np.fft.rfft(m_ext, nfft)
        corr_m = np.fft.irfft(np.conj(fm_l) * fm_e, nfft)[:max_lag + 1]
        num += corr_x
        den += corr_m
        pos = end

    acf = np.zeros(max_lag, dtype=np.float64)
    ok = den[1:] > 100.0
    acf[ok] = (num[1:][ok] / den[1:][ok]) / var

    spec_lim = min(max_lag, 200)
    lags = np.arange(4, spec_lim + 1, dtype=np.float64)
    yd = _hann_taper(_detrend_and_remove_3bp(lags, acf[3:spec_lim]))
    periods = np.round(np.arange(7.0, 14.05, 0.1), 2)
    clean = (np.abs(periods - 9.0) >= 0.35) & (np.abs(periods - 12.0) >= 0.35)
    pc = periods[clean]
    pwrs = np.array([
        float(np.abs(np.sum(yd * np.exp(-2j * math.pi * lags / P))) ** 2) / yd.size
        for P in pc
    ])
    med = float(np.median(pwrs))
    m_hlx = (pc >= 9.5) & (pc <= 11.7)
    m_euk = (pc >= 9.5) & (pc <= 11.0)
    m_prk = (pc >= 10.9) & (pc <= 11.6)
    i_hlx = int(np.flatnonzero(m_hlx)[np.argmax(pwrs[m_hlx])])
    i_euk = int(np.flatnonzero(m_euk)[np.argmax(pwrs[m_euk])])
    i_prk = int(np.flatnonzero(m_prk)[np.argmax(pwrs[m_prk])])
    fg_n, fp_n, psd_n, p_pk_raw, p_pk_adj = notched_fisher_and_band_peak_p(acf, float(pc[i_hlx]), spec_lim=spec_lim)
    b_se10 = bartlett_acf_se(acf[:spec_lim], n_eff=float(den[10]), lag_k=10)

    beat_3bp = float(acf[2] - 0.5 * (acf[1] + acf[3])) if max_lag >= 4 else None
    sat171_ex = float(acf[170] - 0.5 * (acf[167] + acf[173])) if max_lag >= 174 else None
    sat342_ex = float(acf[341] - 0.5 * (acf[338] + acf[344])) if max_lag >= 345 else None
    sat513_ex = float(acf[512] - 0.5 * (acf[509] + acf[515])) if max_lag >= 516 else None
    sat2052_ex = float(acf[2051] - 0.5 * (acf[2048] + acf[2054])) if max_lag >= 2056 else None

    return {
        "max_lag_bp": int(max_lag),
        "gc_global": mu,
        "helix_peak_period_nt": float(pc[i_hlx]),
        "helix_peak_snr": float(pwrs[i_hlx] / max(med, 1e-30)),
        "euk_band_period_nt": float(pc[i_euk]),
        "euk_band_snr": float(pwrs[i_euk] / max(med, 1e-30)),
        "prok_band_period_nt": float(pc[i_prk]),
        "prok_band_snr": float(pwrs[i_prk] / max(med, 1e-30)),
        "fisher_g": fg_n,
        "fisher_p": fp_n,
        "fisher_p_notched": fp_n,
        "peak_p_raw": p_pk_raw,
        "peak_p_adj": p_pk_adj,
        "band_noise_psd": psd_n,
        "bartlett_se_d10": b_se10,
        "beat_3bp": beat_3bp,
        "acf_d10": float(acf[9]),
        "acf_d11": float(acf[10]),
        "acf_d147": float(acf[146]) if max_lag >= 147 else None,
        "acf_d171": float(acf[170]) if max_lag >= 171 else None,
        "acf_d185": float(acf[184]) if max_lag >= 185 else None,
        "acf_d200": float(acf[199]) if max_lag >= 200 else None,
        "sat171_peak_excess": sat171_ex,
        "sat342_hor2_excess": sat342_ex,
        "sat513_hor3_excess": sat513_ex,
        "sat2052_hor12_excess": sat2052_ex,
        "acf_values_truncated_at": min(max_lag, 250),
        "acf_values": acf[:min(max_lag, 250)].tolist(),
    }



# =====================================================================
# 6b. Exact WW (AA/TT/TA) Nucleosome Dinucleotide FFT & THM-D11 Arrow Audit
# =====================================================================
_DS_PROJ = np.zeros((NC, NC), dtype=np.float64)
_DO_PROJ = np.zeros((NC, NC), dtype=np.float64)
for _s in range(16):
    _ix = np.flatnonzero(SUF16 == _s); _DS_PROJ[np.ix_(_ix, _ix)] = 1.0 / _ix.size
    _iy = np.flatnonzero(OUT16 == _s); _DO_PROJ[np.ix_(_iy, _iy)] = 1.0 / _iy.size

_HAS_CG_TRIP = np.array(["CG" in c for c in TRIPLETS], dtype=bool)
_JUNC_CG_MASK = (C_N3[:, None] == 1) & (C_N1[None, :] == 2)
_ANY_CG_PAIR_MASK = _HAS_CG_TRIP[:, None] | _HAS_CG_TRIP[None, :] | _JUNC_CG_MASK

ENTROPY_STRATA_EDGES = [0.0, 4.50, 5.00, 5.30, 5.50, 5.65, 5.75, 5.85, 6.01]
ENTROPY_STRATA_LABELS = ["<4.50", "4.50-5.00", "5.00-5.30", "5.30-5.50",
                         "5.50-5.65", "5.65-5.75", "5.75-5.85", ">=5.85"]


_DC2_PROJ = _DS_PROJ @ _DO_PROJ


def pi_env_project(M):
    """[THM-D13] Exact V_4-closed orthogonal screen envelope projector:
        Pi_env(M) = D_S M D_O + D_O M D_S - D_C2 M D_C2,
    commuting with ALL of V_4 = {I, T, R, S} (rank 496; complement 3600;
    quadrant dims (132, 118, 130, 116) + (924, 906, 894, 876))."""
    return _DS_PROJ @ M @ _DO_PROJ + _DO_PROJ @ M @ _DS_PROJ - _DC2_PROJ @ M @ _DC2_PROJ


def _build_wmp_env_diag_tables():
    """[THM-V16] Exact 64x64 self-energy tables D_env[u,v] = ||P_-+ Pi_env(E_{u,v})||_F^2
    and D_nonenv[u,v] = ||P_-+ (I - Pi_env)(E_{u,v})||_F^2 (traces = 130 and 894)."""
    u = np.arange(NC)[:, None]
    v = np.arange(NC)[None, :]
    ru = MIRROR[u]
    rv = MIRROR[v]
    m_I = (_DS_PROJ[u, u] * _DO_PROJ[v, v]
           + _DO_PROJ[u, u] * _DS_PROJ[v, v]
           - _DC2_PROJ[u, u] * _DC2_PROJ[v, v])
    m_T = 2.0 * _DS_PROJ[v, u] * _DO_PROJ[v, u] - (_DC2_PROJ[v, u] ** 2)
    m_R = (_DS_PROJ[rv, u] * _DO_PROJ[v, ru]
           + _DO_PROJ[rv, u] * _DS_PROJ[v, ru]
           - _DC2_PROJ[rv, u] * _DC2_PROJ[v, ru])
    m_S = (_DS_PROJ[ru, u] * _DO_PROJ[v, rv]
           + _DO_PROJ[ru, u] * _DS_PROJ[v, rv]
           - _DC2_PROJ[ru, u] * _DC2_PROJ[v, rv])
    d_env = 0.25 * (m_I - m_T + m_R - m_S)
    d_wmp = 0.25 * (1.0 - (u == v).astype(np.float64) + (v == ru).astype(np.float64))
    return d_env, d_wmp - d_env


_D_ENV_WMP_DIAG, _D_NONENV_WMP_DIAG = _build_wmp_env_diag_tables()


def exact_ww_autocorr_fft(nt, max_lag=600, chunk_size=2_000_000):
    """[THM-D12/D16] Exact zero-RNG FFT autocorrelation of WW={AA,AT,TA,TT} and
    SS={CC,CG,GC,GG} dinucleotides plus WW-SS cross-correlation (antiphase ~180 deg)
    and WW-SS composite contrast, with 3-bp OLS pre-whitening, notched Fisher g,
    and targeted empirical-PSD nucleosome-band p_adj (Sidak m_eff=14)."""
    L = int(nt.size)
    if L <= max_lag + 12:
        max_lag = min(max_lag, L - 15)
    if max_lag < 20:
        return None
    n_val = n_ww = n_ss = 0
    for s0 in range(0, L - 1, chunk_size):
        s1 = min(L - 1, s0 + chunk_size)
        a, b = nt[s0:s1], nt[s0 + 1:s1 + 1]
        if int(a.min()) >= 4 and int(b.min()) >= 4:
            continue
        ok = (a < 4) & (b < 4)
        nv = int(np.count_nonzero(ok))
        if nv:
            n_val += nv
            ww = ((a == 0) | (a == 3)) & ((b == 0) | (b == 3))
            ss = ((a == 1) | (a == 2)) & ((b == 1) | (b == 2))
            n_ww += int(np.count_nonzero(ww & ok))
            n_ss += int(np.count_nonzero(ss & ok))
    if n_val <= max_lag + 10:
        return None
    mu_w = float(n_ww) / float(n_val)
    mu_s = float(n_ss) / float(n_val)
    var_w = mu_w * (1.0 - mu_w)
    var_s = mu_s * (1.0 - mu_s)
    if var_w < 1e-12 or var_s < 1e-12:
        return None

    num_w = np.zeros(max_lag + 1, dtype=np.float64)
    num_s = np.zeros(max_lag + 1, dtype=np.float64)
    num_ws = np.zeros(max_lag + 1, dtype=np.float64)
    den = np.zeros(max_lag + 1, dtype=np.float64)
    pos, L1 = 0, L - 1
    while pos < L1:
        end = min(L1, pos + chunk_size)
        ext_end = min(L1, end + max_lag)
        a, b = nt[pos:ext_end], nt[pos + 1:ext_end + 1]
        if a.size == 0 or (int(a.min()) >= 4 and int(b.min()) >= 4):
            pos = end; continue
        ok = (a < 4) & (b < 4)
        if not np.any(ok):
            pos = end; continue
        ww = ((a == 0) | (a == 3)) & ((b == 0) | (b == 3))
        ss = ((a == 1) | (a == 2)) & ((b == 1) | (b == 2))
        m_ext = ok.astype(np.float64)
        xw_ext = np.where(ok, ww.astype(np.float64) - mu_w, 0.0)
        xs_ext = np.where(ok, ss.astype(np.float64) - mu_s, 0.0)
        n_left = end - pos
        m_left = m_ext.copy(); m_left[n_left:] = 0.0
        xw_left = xw_ext.copy(); xw_left[n_left:] = 0.0
        xs_left = xs_ext.copy(); xs_left[n_left:] = 0.0
        nfft = 1 << int(math.ceil(math.log2(a.size + max_lag + 2)))
        fw_l = np.fft.rfft(xw_left, nfft); fw_e = np.fft.rfft(xw_ext, nfft)
        fs_l = np.fft.rfft(xs_left, nfft); fs_e = np.fft.rfft(xs_ext, nfft)
        fm_l = np.fft.rfft(m_left, nfft); fm_e = np.fft.rfft(m_ext, nfft)
        num_w += np.fft.irfft(np.conj(fw_l) * fw_e, nfft)[:max_lag + 1]
        num_s += np.fft.irfft(np.conj(fs_l) * fs_e, nfft)[:max_lag + 1]
        num_ws += np.fft.irfft(np.conj(fw_l) * fs_e, nfft)[:max_lag + 1]
        den += np.fft.irfft(np.conj(fm_l) * fm_e, nfft)[:max_lag + 1]
        pos = end

    acf_w = np.zeros(max_lag, dtype=np.float64)
    acf_s = np.zeros(max_lag, dtype=np.float64)
    ccf_ws = np.zeros(max_lag, dtype=np.float64)
    ok_d = den[1:] > 100.0
    acf_w[ok_d] = (num_w[1:][ok_d] / den[1:][ok_d]) / var_w
    acf_s[ok_d] = (num_s[1:][ok_d] / den[1:][ok_d]) / var_s
    ccf_ws[ok_d] = (num_ws[1:][ok_d] / den[1:][ok_d]) / math.sqrt(var_w * var_s)
    acf_ctr = 0.25 * (acf_w + acf_s - 2.0 * ccf_ws)

    spec_lim = min(max_lag, 200)
    lags = np.arange(4, spec_lim + 1, dtype=np.float64)
    yd_w = _hann_taper(_detrend_and_remove_3bp(lags, acf_w[3:spec_lim]))
    yd_s = _hann_taper(_detrend_and_remove_3bp(lags, acf_s[3:spec_lim]))
    yd_ws = _hann_taper(_detrend_and_remove_3bp(lags, ccf_ws[3:spec_lim]))
    yd_c = _hann_taper(_detrend_and_remove_3bp(lags, acf_ctr[3:spec_lim]))

    periods = np.round(np.arange(7.0, 14.05, 0.1), 2)
    clean = (np.abs(periods - 9.0) >= 0.35) & (np.abs(periods - 12.0) >= 0.35)
    pc = periods[clean]
    m_hlx = (pc >= 9.5) & (pc <= 11.7)
    m_euk = (pc >= 9.5) & (pc <= 11.0)

    def _scan_band(yd_vec):
        pw = np.array([float(np.abs(np.sum(yd_vec * np.exp(-2j * math.pi * lags / P))) ** 2) / yd_vec.size for P in pc])
        med_ = float(np.median(pw))
        ih = int(np.flatnonzero(m_hlx)[np.argmax(pw[m_hlx])])
        ie = int(np.flatnonzero(m_euk)[np.argmax(pw[m_euk])])
        return float(pc[ih]), float(pw[ih] / max(med_, 1e-30)), float(pc[ie]), float(pw[ie] / max(med_, 1e-30))

    p_hlx, snr_hlx, p_nuc, snr_nuc = _scan_band(yd_w)
    _, _, ss_nuc_nt, ss_nuc_snr = _scan_band(yd_s)
    _, _, ctr_nuc_nt, ctr_nuc_snr = _scan_band(yd_c)

    zw_pk = complex(np.sum(yd_w * np.exp(-2j * math.pi * lags / p_nuc)))
    zws_pk = complex(np.sum(yd_ws * np.exp(-2j * math.pi * lags / p_nuc)))
    ph_diff = abs(math.degrees(math.atan2(zws_pk.imag, zws_pk.real) - math.atan2(zw_pk.imag, zw_pk.real)))
    if ph_diff > 180.0:
        ph_diff = 360.0 - ph_diff

    fg_w, fp_w, psd_w, p_w_raw, p_w_adj = notched_fisher_and_band_peak_p(acf_w, p_nuc, spec_lim=spec_lim)
    _, _, _, p_c_raw, p_c_adj = notched_fisher_and_band_peak_p(acf_ctr, ctr_nuc_nt, spec_lim=spec_lim)
    b_se10 = bartlett_acf_se(acf_w[:spec_lim], n_eff=float(den[10]), lag_k=10)

    return {
        "ww_global": mu_w,
        "ss_global": mu_s,
        "max_lag_bp": int(max_lag),
        "ww_helix_peak_nt": p_hlx,
        "ww_helix_snr": snr_hlx,
        "ww_nuc_band_nt": p_nuc,
        "ww_nuc_band_snr": snr_nuc,
        "ss_nuc_band_nt": ss_nuc_nt,
        "ss_nuc_band_snr": ss_nuc_snr,
        "ww_minus_ss_nuc_nt": ctr_nuc_nt,
        "ww_minus_ss_nuc_snr": ctr_nuc_snr,
        "ww_ss_phase_diff_deg": float(ph_diff),
        "ww_fisher_g": fg_w,
        "ww_fisher_p": fp_w,
        "ww_fisher_p_notched": fp_w,
        "ww_peak_p_raw": p_w_raw,
        "ww_peak_p_adj": p_w_adj,
        "ww_minus_ss_peak_p_adj": p_c_adj,
        "ww_band_noise_psd": psd_w,
        "ww_bartlett_se_d10": b_se10,
        "ww_acf_d10": float(acf_w[9]),
        "ww_acf_d11": float(acf_w[10]),
        "ss_acf_d10": float(acf_s[9]),
        "ss_acf_d11": float(acf_s[10]),
        "ww_acf_d147": float(acf_w[146]) if max_lag >= 147 else None,
        "ww_acf_d185": float(acf_w[184]) if max_lag >= 185 else None,
        "ww_acf_d200": float(acf_w[199]) if max_lag >= 200 else None,
        "ww_acf_d342": float(acf_w[341]) if max_lag >= 342 else None,
        "ww_acf_d513": float(acf_w[512]) if max_lag >= 513 else None,
    }


def thm_D11_arrow_screen_and_cpg_audit(B):
    """[THM-D11 & THM-D13] Orthogonal decomposition of W_-+ via V_4-closed 496-D screen envelope
    Pi_env (zero quadrant leakage) alongside legacy 256-D one-sided screen, CpG-free core & Top-15 6-mers."""
    kv = klein_v4_pr2_audit(B)
    if kv is None:
        return None
    Wmp = kv["_W_mp"]
    n_tot = float(np.sum(Wmp * Wmp))

    Wmp_env = pi_env_project(Wmp)
    Wmp_nonenv = Wmp - Wmp_env
    n_env = float(np.sum(Wmp_env * Wmp_env))
    n_nonenv = float(np.sum(Wmp_nonenv * Wmp_nonenv))
    N_b = float(kv["N"])
    if N_b > 1.0:
        B_arr = np.asarray(B, dtype=np.float64)
        n_env_unb = (N_b * N_b * n_env - float(np.sum(B_arr * _D_ENV_WMP_DIAG))) / (N_b * (N_b - 1.0))
        n_nonenv_unb = (N_b * N_b * n_nonenv - float(np.sum(B_arr * _D_NONENV_WMP_DIAG))) / (N_b * (N_b - 1.0))
    else:
        n_env_unb, n_nonenv_unb = n_env, n_nonenv
    u_sum_env = max(n_env_unb, 0.0) + max(n_nonenv_unb, 0.0)
    pure_nonenv_unb = (max(n_nonenv_unb, 0.0) / u_sum_env) if u_sum_env > 1e-30 else None
    leak_T = float(np.max(np.abs(Wmp_env + Wmp_env.T)))
    leak_R = float(np.max(np.abs(Wmp_env - Wmp_env.T[np.ix_(MIRROR, MIRROR)])))

    Wmp_br = _DS_PROJ @ Wmp @ _DO_PROJ
    Wmp_nb = Wmp - Wmp_br
    n_br = float(np.sum(Wmp_br * Wmp_br))
    n_nb = float(np.sum(Wmp_nb * Wmp_nb))

    n_cg = float(np.sum(Wmp[_ANY_CG_PAIR_MASK] ** 2))
    n_cgfree = float(np.sum(Wmp[~_ANY_CG_PAIR_MASK] ** 2))
    n_nb_cgfree = float(np.sum(Wmp_nonenv[~_ANY_CG_PAIR_MASK] ** 2))

    flat_idx = np.argsort(Wmp.ravel())[::-1]
    top_all, top_cgfree, seen = [], [], set()
    for fi in flat_idx:
        u, v = divmod(int(fi), NC)
        val = float(Wmp[u, v])
        if val <= 0:
            break
        ru, rv = int(MIRROR[v]), int(MIRROR[u])
        canon = min((u, v), (ru, rv))
        if canon in seen:
            continue
        seen.add(canon)
        orb_mult = 2.0 if (u == ru and v == rv) else 4.0
        is_cg = bool(_ANY_CG_PAIR_MASK[u, v])
        rec = {
            "left": TRIPLETS[u], "right": TRIPLETS[v],
            "motif_6mer": TRIPLETS[u] + "-" + TRIPLETS[v],
            "mirror_6mer": TRIPLETS[ru] + "-" + TRIPLETS[rv],
            "W_mp_flux": val,
            "flux_share_pct": 100.0 * (orb_mult * val * val) / max(n_tot, TINY),
            "has_CpG": is_cg,
        }
        if len(top_all) < 15:
            top_all.append(rec)
        if (not is_cg) and len(top_cgfree) < 12:
            top_cgfree.append(rec)
        if len(top_all) >= 15 and len(top_cgfree) >= 12:
            break

    return {
        "norm2_W_mp": n_tot,
        "norm2_W_mp_env_496": n_env,
        "norm2_W_mp_nonenv_3600": n_nonenv,
        "norm2_W_mp_env_496_unb": float(n_env_unb),
        "norm2_W_mp_nonenv_3600_unb": float(n_nonenv_unb),
        "pure_nonenv_arrow_share": (n_nonenv / n_tot) if n_tot > 1e-30 else None,
        "pure_nonenv_arrow_share_unb": pure_nonenv_unb,
        "env_pythagorean_err": abs(n_tot - (n_env + n_nonenv)),
        "env_quadrant_leak_err": max(leak_T, leak_R),
        "norm2_W_mp_bridge": n_br,
        "norm2_W_mp_nonbridge": n_nb,
        "nonbridge_arrow_share": (n_nb / n_tot) if n_tot > 1e-30 else None,
        "pythagorean_err": abs(n_tot - (n_br + n_nb)),
        "cpg_pair_arrow_share": (n_cg / n_tot) if n_tot > 1e-30 else None,
        "cpg_free_arrow_share": (n_cgfree / n_tot) if n_tot > 1e-30 else None,
        "nonbridge_cpg_free_energy": n_nb_cgfree,
        "top_motifs_all": top_all,
        "top_motifs_cpg_free": top_cgfree,
    }


# =====================================================================
# 7. Sector "Hard Drive" Scanner, Toric SVD Walks & Macro Autocorrelation
# =====================================================================
def _autocorr_1d(arr, lags):
    arr = np.asarray(arr, dtype=np.float64)
    ok = np.isfinite(arr)
    if int(ok.sum()) < 10 or float(np.mean(~ok)) > 0.25:
        return {str(k): None for k in lags}
    x = np.arange(arr.size, dtype=np.float64)
    y = np.interp(x, x[ok], arr[ok])
    y = y - float(np.mean(y))
    denom = float(np.dot(y, y))
    if denom < 1e-30:
        return {str(k): None for k in lags}
    out = {}
    for k in lags:
        if k >= y.size - 3:
            out[str(k)] = None
        else:
            out[str(k)] = float(np.dot(y[:-k], y[k:]) / denom)
    return out


def scan_chromosome_sectors(nt, sector_bp, topology="linear"):
    """[THM-D5] Scans chromosome sector-by-sector as a physical hard drive.
    Returns per-sector tracks, exact SVD-signed toric walks S_+- & S_--,
    blind oriC/ter extrema, and macro sector-to-sector autocorrelation."""
    L = int(nt.size)
    n_sec = int(L // sector_bp)
    if n_sec < 6:
        sector_bp = max(2000, L // 10)
        n_sec = int(L // sector_bp)

    tracks = {
        "pos_mb": [], "valid_frac": [], "gc": [], "gc_skew": [], "at_skew": [],
        "H_3mer_bits": [],
        "W_pp": [], "W_pm": [], "W_mp": [], "W_mm": [],
        "W_mp_unb": [], "W_mm_unb": [], "p_pal_6bp": [],
        "W_mp_floor_sec": [], "W_mp_snr_sec": [],
        "PR2_defect_F2": [], "time_arrow_F2": [], "PR2_dyn_share": [],
        "A_PR2_L1": [], "A_time_L1": [],
        "I_bridge": [], "D_5lead": [], "D_3wob": [], "D_1_3": [], "I_nonbr_MM": [],
        "Omega_phase": [], "Omega_cond": [], "Omega_cond_excess": [],
    }
    Wpm_rows, Wmm_rows, valid_sec_idx = [], [], []
    B_strata = np.zeros((len(ENTROPY_STRATA_LABELS), NC, NC), dtype=np.float64)
    strata_counts = np.zeros(len(ENTROPY_STRATA_LABELS), dtype=np.int64)
    strata_neg_mp = np.zeros(len(ENTROPY_STRATA_LABELS), dtype=np.int64)
    block_bp = 1_000_000 if L >= 5_000_000 else max(sector_bp, L // 10)
    secs_per_blk = max(1, int(round(block_bp / float(sector_bp))))
    block_mats_1mb, B_blk_cur, blk_sec_cnt = [], np.zeros((NC, NC), dtype=np.float64), 0

    for w in range(n_sec):
        s0 = w * sector_bp
        s1 = s0 + sector_bp
        seg = nt[s0:s1]
        pos_mb = (s0 + 0.5 * sector_bp) / 1e6
        tracks["pos_mb"].append(pos_mb)
        v_mask = seg < 4
        v_frac = float(np.mean(v_mask))
        tracks["valid_frac"].append(v_frac)
        if v_frac < 0.5:
            for k in tracks:
                if k not in ("pos_mb", "valid_frac"):
                    tracks[k].append(float("nan"))
            continue
        cnt4 = np.bincount(seg[v_mask], minlength=4).astype(np.float64)
        n_acgt = float(cnt4.sum())
        gc = (cnt4[1] + cnt4[2]) / max(n_acgt, 1.0)
        gc_sk = (cnt4[2] - cnt4[1]) / max(cnt4[2] + cnt4[1], 1.0)
        at_sk = (cnt4[0] - cnt4[3]) / max(cnt4[0] + cnt4[3], 1.0)
        c64 = count_triplets_1d(seg, chunk_size=sector_bp)
        p64 = c64 / max(float(c64.sum()), 1.0)
        m64 = p64 > 0
        H3 = float(-np.sum(p64[m64] * np.log2(p64[m64])))

        Bw = count_triplet_pairs_lag(seg, d=3, chunk_size=sector_bp)
        kv = klein_v4_pr2_audit(Bw)
        sc = two_sided_screen_full(Bw)
        om = frame_phase_lock_sector(seg, d=3)

        tracks["gc"].append(gc)
        tracks["gc_skew"].append(gc_sk)
        tracks["at_skew"].append(at_sk)
        tracks["H_3mer_bits"].append(H3)
        tracks["W_pp"].append(kv["norm2_W_pp"])
        tracks["W_pm"].append(kv["norm2_W_pm"])
        tracks["W_mp"].append(kv["norm2_W_mp"])
        tracks["W_mm"].append(kv["norm2_W_mm"])
        tracks["W_mp_unb"].append(kv["norm2_W_mp_unb"])
        tracks["W_mm_unb"].append(kv["norm2_W_mm_unb"])
        tracks["p_pal_6bp"].append(kv["p_pal_6bp"])
        Mw_s = float(c64.sum())
        if Mw_s > 1.0:
            Ew_sec = np.maximum(np.outer(c64, c64) - np.diag(c64), 0.0)
            fl_w = klein_noise_floor_from_null(Ew_sec, kv["N"])["W_mp"]
        else:
            fl_w = (1.0 - kv["p_diag_homodimer"] + kv["p_pal_6bp"]) / (4.0 * max(float(kv["N"]), 1.0))
        tracks["W_mp_floor_sec"].append(fl_w)
        tracks["W_mp_snr_sec"].append(kv["norm2_W_mp_unb"] / max(fl_w, 1e-30))
        b_idx = int(np.searchsorted(ENTROPY_STRATA_EDGES[1:], H3, side="right"))
        b_idx = min(max(b_idx, 0), len(ENTROPY_STRATA_LABELS) - 1)
        B_strata[b_idx] += Bw
        strata_counts[b_idx] += 1
        if kv["norm2_W_mp_unb"] < 0.0:
            strata_neg_mp[b_idx] += 1
        B_blk_cur += Bw
        blk_sec_cnt += 1
        if blk_sec_cnt >= secs_per_blk:
            block_mats_1mb.append(B_blk_cur)
            B_blk_cur = np.zeros((NC, NC), dtype=np.float64)
            blk_sec_cnt = 0
        tracks["PR2_defect_F2"].append(kv["PR2_defect_F2"])
        tracks["time_arrow_F2"].append(kv["time_arrow_F2"])
        tracks["PR2_dyn_share"].append(kv["PR2_dyn_share"] if kv["PR2_dyn_share"] is not None else float("nan"))
        tracks["A_PR2_L1"].append(kv["A_PR2_L1"])
        tracks["A_time_L1"].append(kv["A_time_L1"])
        tracks["I_bridge"].append(sc["I_bridge"])
        tracks["D_5lead"].append(sc["D_5prime_lead"])
        tracks["D_3wob"].append(sc["D_3prime_wobble"])
        tracks["D_1_3"].append(sc["D_1_3"])
        tracks["I_nonbr_MM"].append(sc["I_nonbridge_MM"])
        tracks["Omega_phase"].append(om["Omega_phase"] if om else float("nan"))
        tracks["Omega_cond"].append(om["Omega_cond"] if om else float("nan"))
        tracks["Omega_cond_excess"].append(om["Omega_cond_excess"] if om else float("nan"))

        Wpm_rows.append(kv["_W_pm"].ravel())
        Wmm_rows.append(kv["_W_mm"].ravel())
        valid_sec_idx.append(w)

    if blk_sec_cnt > 0 and float(B_blk_cur.sum()) > 0.0:
        block_mats_1mb.append(B_blk_cur)
    cluster_1mb = cluster_1mb_arrow_test(block_mats_1mb)
    if cluster_1mb is not None:
        cluster_1mb["block_bp"] = int(secs_per_blk * sector_bp)
    del block_mats_1mb, B_blk_cur
    _wmp_u_v = np.asarray(tracks["W_mp_unb"], dtype=np.float64)
    _wmp_u_ok = _wmp_u_v[np.isfinite(_wmp_u_v)]
    frac_neg_u_mp = float(np.mean(_wmp_u_ok < 0.0)) if _wmp_u_ok.size else None

    gc_sk_arr = np.asarray(tracks["gc_skew"], dtype=np.float64)
    gc_sk_clean = np.where(np.isfinite(gc_sk_arr), gc_sk_arr, 0.0)
    if topology == "circular":
        gc_sk_walk = np.cumsum(gc_sk_clean - float(np.mean(gc_sk_clean)))
    else:
        gc_sk_walk = np.cumsum(gc_sk_clean)

    def _svd_signed_walk(rows_list, anchor_vec):
        if len(rows_list) < 5:
            z = np.zeros(n_sec, dtype=np.float64)
            return z, z, z, 0.0, None, None
        X = np.vstack(rows_list)
        u1, var_exp, _ = _top_svd_intrinsic(X)
        proj = np.full(n_sec, float("nan"), dtype=np.float64)
        proj[valid_sec_idx] = u1
        sd = float(np.nanstd(proj))
        proj_z = (proj / sd) if sd > 1e-18 else np.zeros(n_sec)
        fill = np.where(np.isfinite(proj_z), proj_z, 0.0)
        walk_ctr = np.cumsum(fill - float(np.mean(fill)))
        walk = walk_ctr if topology == "circular" else np.cumsum(fill)
        anc = anchor_vec[valid_sec_idx]
        ok_a = np.isfinite(anc) & np.isfinite(u1)
        r_gc = (float(np.corrcoef(u1[ok_a], anc[ok_a])[0, 1])
                if int(ok_a.sum()) >= 5 and np.std(anc[ok_a]) > 1e-15 and np.std(u1[ok_a]) > 1e-15
                else None)
        abs_r = abs(r_gc) if r_gc is not None else None
        return proj_z, walk, walk_ctr, var_exp, r_gc, abs_r

    u_pm, walk_pm, walk_pm_ctr, var_pm, r_pm_gc, abs_r_pm_gc = _svd_signed_walk(Wpm_rows, gc_sk_arr)
    u_mm, walk_mm, walk_mm_ctr, var_mm, r_mm_gc, abs_r_mm_gc = _svd_signed_walk(Wmm_rows, gc_sk_arr)
    levy_all = levy_area_signed_test(u_pm, u_mm)
    tracks["u_W_pm_signed"] = u_pm.tolist()
    tracks["u_W_mm_signed"] = u_mm.tolist()
    tracks["walk_gc_skew"] = gc_sk_walk.tolist()
    tracks["walk_W_pm"] = walk_pm.tolist()
    tracks["walk_W_mm"] = walk_mm.tolist()

    pos_arr = np.asarray(tracks["pos_mb"], dtype=np.float64)

    def _walk_extrema(w_arr):
        imin = int(np.argmin(w_arr))
        imax = int(np.argmax(w_arr))
        return {
            "min_mb": float(pos_arr[imin]),
            "max_mb": float(pos_arr[imax]),
            "span": float(w_arr[imax] - w_arr[imin]),
        }

    macro_lags = [1, 2, 5, 10, 20, 50]
    macro_acf = {
        ch: _autocorr_1d(tracks[ch], macro_lags)
        for ch in ("gc", "W_pp", "W_pm", "W_mp", "W_mm", "I_nonbr_MM", "Omega_cond", "H_3mer_bits")
    }

    summary = {}
    for ch in ("gc", "H_3mer_bits", "W_pp", "W_pm", "W_mp", "W_mm",
               "PR2_dyn_share", "I_nonbr_MM", "Omega_cond", "Omega_cond_excess"):
        a = np.asarray(tracks[ch], dtype=np.float64)
        a = a[np.isfinite(a)]
        if a.size:
            summary[ch] = {
                "median": float(np.median(a)),
                "p10": float(np.percentile(a, 10)),
                "p90": float(np.percentile(a, 90)),
                "max": float(np.max(a)),
            }

    return {
        "sector_bp": int(sector_bp),
        "n_sectors": int(n_sec),
        "n_valid_sectors": len(valid_sec_idx),
        "svd_W_pm_var1": float(var_pm),
        "svd_W_mm_var1": float(var_mm),
        "corr_u_W_pm_vs_gc_skew": r_pm_gc,
        "abs_corr_u_W_pm_vs_gc_skew": abs_r_pm_gc,
        "corr_u_W_mm_vs_gc_skew": r_mm_gc,
        "abs_corr_u_W_mm_vs_gc_skew": abs_r_mm_gc,
        "extrema_W_pm_walk_centered": _walk_extrema(walk_pm_ctr),
        "extrema_W_mm_walk_centered": _walk_extrema(walk_mm_ctr),
        "levy_D17": levy_all,
        "corr_u_W_pm_vs_W_mm": (
            float(np.corrcoef(u_pm[valid_sec_idx], u_mm[valid_sec_idx])[0, 1])
            if len(valid_sec_idx) >= 5 and np.std(u_pm[valid_sec_idx]) > 1e-12 and np.std(u_mm[valid_sec_idx]) > 1e-12
            else None
        ),
        "extrema_gc_skew_walk": _walk_extrema(gc_sk_walk),
        "extrema_W_pm_walk": _walk_extrema(walk_pm),
        "extrema_W_mm_walk": _walk_extrema(walk_mm),
        "sector_summary": summary,
        "macro_autocorr": macro_acf,
        "tracks": tracks,
        "_Wpm_rows": Wpm_rows,
        "_Wmm_rows": Wmm_rows,
        "_valid_sec_idx": valid_sec_idx,
        "_B_strata": B_strata,
        "entropy_strata_sector_counts": strata_counts.tolist(),
        "entropy_strata_neg_mp_counts": strata_neg_mp.tolist(),
        "frac_negative_u_mp": frac_neg_u_mp,
        "cluster_1mb_D14": cluster_1mb,
    }


# =====================================================================
# 8. Single-Chromosome Complete Driver
# =====================================================================

# =====================================================================
# 7b. THEOREM D7 (Blind Euchromatin vs Satellite Partition) &
#     THEOREM D8 (Multi-Scale Fractal / Hurst Disk Scaling)
# =====================================================================
EUCH_ENTROPY_THRESHOLD_BITS = 5.65


def _svd_unit_normalized_walk(rows_list, valid_idx, n_sec, anchor_vec, topology="linear"):
    # Unit-normalizes each sector Klein matrix before SVD with INTRINSIC v1 sign (zero GC circularity).
    if len(rows_list) < 5:
        z = np.zeros(n_sec, dtype=np.float64)
        return z, z, 0.0, None
    X_raw = np.vstack(rows_list)
    norms = np.linalg.norm(X_raw, axis=1, keepdims=True)
    X = X_raw / np.maximum(norms, 1e-30)
    u1, var_exp, _ = _top_svd_intrinsic(X)
    proj = np.full(n_sec, float("nan"), dtype=np.float64)
    proj[valid_idx] = u1
    sd = float(np.nanstd(proj))
    proj_z = (proj / sd) if sd > 1e-18 else np.zeros(n_sec)
    fill = np.where(np.isfinite(proj_z), proj_z, 0.0)
    if topology == "circular":
        walk = np.cumsum(fill - float(np.mean(fill)))
    else:
        walk = np.cumsum(fill)
    anc = anchor_vec[valid_idx]
    ok_a = np.isfinite(anc) & np.isfinite(u1)
    r_gc = (float(np.corrcoef(u1[ok_a], anc[ok_a])[0, 1])
            if int(ok_a.sum()) >= 5 and np.std(anc[ok_a]) > 1e-15 and np.std(u1[ok_a]) > 1e-15
            else None)
    return proj_z, walk, var_exp, r_gc


def thm_D7_euchromatin_vs_satellite_audit(nt, disk_scan, sector_bp, topology="linear",
                                          h_thresh=None, micro_lags=None):
    """[THM-D7 v1.5] Unified parameter-free Euchromatin vs Satellite/Tandem partition.
    If h_thresh is None or <= 0 (default): uses log10(I_nonbr_MM) Otsu-BIC + entropy inversion
    + exact Shannon half-entropy guard (H_3 < 3.0 bits). Computes euchromatin-only 1-Mb HAC
    cluster test, euchromatin-only D15 (M=1000) residual, and optional 10-lag euchromatin counts."""
    tr = disk_scan["tracks"]
    H3 = np.asarray(tr["H_3mer_bits"], dtype=np.float64)
    INB = np.asarray(tr["I_nonbr_MM"], dtype=np.float64)
    INB_RAW = (np.asarray(tr["D_5lead"], dtype=np.float64)
               + np.asarray(tr["D_3wob"], dtype=np.float64)
               + np.asarray(tr["D_1_3"], dtype=np.float64)) if "D_5lead" in tr else None
    n_sec = int(H3.size)
    is_valid = np.isfinite(H3) & np.isfinite(INB)
    use_legacy = (h_thresh is not None) and (float(h_thresh) > 0.0)
    if use_legacy:
        is_euch = is_valid & (H3 >= float(h_thresh))
        is_sat = is_valid & (H3 < float(h_thresh))
        pf_info = {
            "partition_method": "legacy_fixed_h_thresh",
            "half_entropy_guard_bits": float(h_thresh),
            "n_micro_sat_h3_lt_3": int(np.count_nonzero(is_valid & (H3 < 3.0))),
            "inb_threshold_otsu": None,
            "bic_2vs1_delta": None,
            "entropy_inversion_ok": None,
            "frac_clipped": 0.0,
        }
        h_report = float(h_thresh)
    else:
        is_euch, is_sat, pf_info = partition_euch_sat_parameter_free(H3, INB, INB_RAW=INB_RAW)
        h_report = float(np.nanmin(H3[is_euch])) if np.any(is_euch) else 3.0

    n_val = int(np.count_nonzero(is_valid))
    n_euch = int(np.count_nonzero(is_euch))
    n_sat = int(np.count_nonzero(is_sat))

    def _audit_partition(mask_sec, label):
        idx_sec = np.flatnonzero(mask_sec)
        if idx_sec.size < 3:
            return {"label": label, "n_sectors": int(idx_sec.size),
                    "frac_of_valid": float(idx_sec.size) / max(1, n_val),
                    "skipped": True}
        L = int(nt.size)
        nt_part = np.full(L, 255, dtype=np.uint8)
        for w in idx_sec:
            s0 = int(w) * sector_bp
            s1 = min(L, s0 + sector_bp)
            nt_part[s0:s1] = nt[s0:s1]
        B3 = count_triplet_pairs_lag(nt_part, d=3)
        kv3_full = klein_v4_pr2_audit(B3)
        kv3 = strip_private(kv3_full)
        sc3 = two_sided_screen_full(B3)
        Bw3, Eh3 = sector_het_expected_counts(nt_part, sector_bp=sector_bp, d=3)
        skl3 = screen_kl_chain_rule(Bw3, Eh3, df_eff=sc3["df_eff"])
        B90 = count_triplet_pairs_lag(nt_part, d=90)
        sc90 = two_sided_screen_full(B90)
        Bw90, Eh90 = sector_het_expected_counts(nt_part, sector_bp=sector_bp, d=90)
        skl90 = screen_kl_chain_rule(Bw90, Eh90, df_eff=sc90["df_eff"])
        fft_lag = 2500 if label == "satellite_firmware" else 600
        fft_p = exact_gc_autocorr_fft(nt_part, max_lag=fft_lag)
        fft_compact = {k: v for k, v in fft_p.items() if k != "acf_values"} if fft_p is not None else None
        ww_p = exact_ww_autocorr_fft(nt_part, max_lag=600) if label == "euchromatin" else None
        d11_p = thm_D11_arrow_screen_and_cpg_audit(B3)

        cl_euch = None
        d15_euch = None
        b_lags_part = None
        if label == "euchromatin":
            block_bp = 1_000_000 if L >= 5_000_000 else max(sector_bp, L // 10)
            secs_per_blk = max(1, int(round(block_bp / float(sector_bp))))
            blk_mats, B_cur, cnt_b = [], np.zeros((NC, NC), dtype=np.float64), 0
            for w in idx_sec:
                s0 = int(w) * sector_bp
                s1 = min(L, s0 + sector_bp)
                B_cur += count_triplet_pairs_lag(nt[s0:s1], d=3, chunk_size=sector_bp)
                cnt_b += 1
                if cnt_b >= secs_per_blk:
                    blk_mats.append(B_cur)
                    B_cur = np.zeros((NC, NC), dtype=np.float64)
                    cnt_b = 0
            if cnt_b > 0 and float(B_cur.sum()) > 0.0:
                blk_mats.append(B_cur)
            cl_euch = cluster_1mb_arrow_test(blk_mats)
            if cl_euch is not None:
                cl_euch["block_bp"] = int(secs_per_blk * sector_bp)
            del blk_mats, B_cur

            n_part_val = int(idx_sec.size * sector_bp)
            if n_part_val >= 2020 and kv3_full is not None:
                E_loc_e, m_ex_e = compositional_current_null(nt_part, window_m=1000, d=3, trip_budget=4_000_000)
                if float(E_loc_e.sum()) > 0.0:
                    kv_le = klein_v4_pr2_audit(E_loc_e)
                    W_om, W_cm = kv3_full["_W_mp"], kv_le["_W_mp"]
                    n_om = float(np.sum(W_om * W_om))
                    n_cm = float(np.sum(W_cm * W_cm))
                    n_rm = float(np.sum((W_om - W_cm) ** 2))
                    cos_e = (float(np.sum(W_om * W_cm)) / math.sqrt(n_om * n_cm)) if (n_om > 1e-30 and n_cm > 1e-30) else None
                    cov_e = float(m_ex_e) / max(float(B3.sum()), 1.0)
                    d15_euch = {
                        "window_m": 1000,
                        "coverage_frac": cov_e,
                        "coverage_ok": bool(cov_e >= 0.95),
                        "E_loc_W_mp": n_cm,
                        "W_mp_current_ratio": n_cm / max(n_om, TINY),
                        "W_mp_residual_share": n_rm / max(n_om, TINY),
                        "cos_W_mp_obs_vs_cur": cos_e,
                    }

            if micro_lags:
                b_lags_part = {}
                for lag_d in micro_lags:
                    ld = int(lag_d)
                    if ld == 3:
                        b_lags_part[ld] = B3
                    elif ld == 90:
                        b_lags_part[ld] = B90
                    else:
                        b_lags_part[ld] = count_triplet_pairs_lag(nt_part, d=ld)

        del nt_part
        return {
            "label": label,
            "n_sectors": int(idx_sec.size),
            "frac_of_valid": float(idx_sec.size) / max(1, n_val),
            "mb_total": float(idx_sec.size * sector_bp) / 1e6,
            "skipped": False,
            "klein_d3": kv3,
            "screen_d3": sc3,
            "screen_kl_d3": skl3,
            "screen_kl_d90": skl90,
            "fft_pitch": fft_compact,
            "ww_nucleosome_fft": ww_p,
            "arrow_D11": d11_p,
            "cluster_1mb_euch_D14": cl_euch,
            "comp_current_euch_M1000": d15_euch,
            "_B3_part": B3,
            "_B_lags_part": b_lags_part,
        }

    euch_res = _audit_partition(is_euch, "euchromatin")
    sat_res = _audit_partition(is_sat, "satellite_firmware")

    Wpm_e, Wmm_e, idx_e = [], [], []
    if "_Wpm_rows" in disk_scan and "_valid_sec_idx" in disk_scan:
        for row_pm, row_mm, w in zip(disk_scan["_Wpm_rows"], disk_scan["_Wmm_rows"], disk_scan["_valid_sec_idx"]):
            if is_euch[w]:
                Wpm_e.append(row_pm); Wmm_e.append(row_mm); idx_e.append(int(w))
    else:
        L = int(nt.size)
        for w in np.flatnonzero(is_euch):
            s0 = int(w) * sector_bp
            Bw = count_triplet_pairs_lag(nt[s0:min(L, s0 + sector_bp)], d=3, chunk_size=sector_bp)
            kv = klein_v4_pr2_audit(Bw)
            if kv is not None:
                Wpm_e.append(kv["_W_pm"].ravel()); Wmm_e.append(kv["_W_mm"].ravel()); idx_e.append(int(w))

    gc_sk = np.asarray(tr["gc_skew"], dtype=np.float64)
    gc_sk_e = np.where(is_euch & np.isfinite(gc_sk), gc_sk, 0.0)
    walk_gc_e = (
        np.cumsum(gc_sk_e - (float(np.mean(gc_sk_e[is_euch])) if n_euch else 0.0))
        if topology == "circular" else np.cumsum(gc_sk_e)
    )

    u_pm_e, walk_pm_e, var_pm_e, r_pm_gc_e = _svd_unit_normalized_walk(
        Wpm_e, idx_e, n_sec, gc_sk, topology=topology)
    u_mm_e, walk_mm_e, var_mm_e, r_mm_gc_e = _svd_unit_normalized_walk(
        Wmm_e, idx_e, n_sec, gc_sk, topology=topology)
    r_pm_mm_e = (
        float(np.corrcoef(u_pm_e[idx_e], u_mm_e[idx_e])[0, 1])
        if len(idx_e) >= 5 and np.std(u_pm_e[idx_e]) > 1e-12 and np.std(u_mm_e[idx_e]) > 1e-12
        else None
    )

    sat_blocks = []
    pos = np.asarray(tr["pos_mb"], dtype=np.float64)
    half_mb = 0.5 * sector_bp / 1e6
    in_blk = False
    b_start = 0
    for w in range(n_sec):
        if is_sat[w] and not in_blk:
            in_blk = True
            b_start = w
        elif not is_sat[w] and in_blk:
            if w - b_start >= 2:
                sat_blocks.append([float(pos[b_start] - half_mb), float(pos[w - 1] + half_mb)])
            in_blk = False
    if in_blk and n_sec - b_start >= 2:
        sat_blocks.append([float(pos[b_start] - half_mb), float(pos[n_sec - 1] + half_mb)])

    otsu_h3 = otsu_threshold_1d(H3)
    wmp_u_arr = np.asarray(tr.get("W_mp_unb", tr["W_mp"]), dtype=np.float64)
    wmm_u_arr = np.asarray(tr.get("W_mm_unb", tr["W_mm"]), dtype=np.float64)
    thresh_sweep = {}
    for th_v in (5.50, 5.65, 5.80):
        me = is_valid & (H3 >= th_v)
        ne = int(np.count_nonzero(me))
        sm = float(np.nansum(np.maximum(wmp_u_arr[me], 0.0))) if ne else 0.0
        sc = float(np.nansum(np.maximum(wmm_u_arr[me], 0.0))) if ne else 0.0
        thresh_sweep[f"{th_v:.2f}"] = {
            "euch_frac": float(ne) / max(1, n_val),
            "n_euch": ne,
            "pr2_ok_unb_sector_share": (sm / (sm + sc)) if (sm + sc) > 1e-30 else None,
            "frac_negative_u_mp": float(np.mean(wmp_u_arr[me] < 0.0)) if ne else None,
        }
    levy_euch = levy_area_signed_test(u_pm_e[idx_e], u_mm_e[idx_e]) if len(idx_e) >= 6 else None

    return {
        "partition_method": pf_info["partition_method"],
        "inb_threshold_otsu": pf_info["inb_threshold_otsu"],
        "bic_2vs1_delta": pf_info["bic_2vs1_delta"],
        "entropy_inversion_ok": pf_info["entropy_inversion_ok"],
        "n_micro_sat_h3_lt_3": pf_info["n_micro_sat_h3_lt_3"],
        "h_threshold_bits": h_report,
        "min_euch_h3_bits": h_report,
        "frac_clipped": pf_info.get("frac_clipped", 0.0),
        "otsu_h3_bits": otsu_h3,
        "h_threshold_sensitivity": thresh_sweep,
        "n_valid_sectors": n_val,
        "n_euch_sectors": n_euch,
        "n_sat_sectors": n_sat,
        "euch_fraction": float(n_euch) / max(1, n_val),
        "sat_fraction": float(n_sat) / max(1, n_val),
        "satellite_blocks_mb": sat_blocks,
        "euchromatin": euch_res,
        "satellite_firmware": sat_res,
        "euch_svd_walk": {
            "svd1_W_pm_unit": float(var_pm_e),
            "svd1_W_mm_unit": float(var_mm_e),
            "corr_W_pm_vs_gc_skew": r_pm_gc_e,
            "abs_corr_W_pm_vs_gc_skew": abs(r_pm_gc_e) if r_pm_gc_e is not None else None,
            "corr_W_mm_vs_gc_skew": r_mm_gc_e,
            "abs_corr_W_mm_vs_gc_skew": abs(r_mm_gc_e) if r_mm_gc_e is not None else None,
            "corr_W_pm_vs_W_mm": r_pm_mm_e,
            "levy_D17_euch": levy_euch,
            "walk_gc_euch": walk_gc_e.tolist(),
            "walk_W_pm_euch": walk_pm_e.tolist(),
            "walk_W_mm_euch": walk_mm_e.tolist(),
        },
    }


def _hurst_wls_on_scales(arr, scales, min_blks=8):
    ln10_sq = math.log(10.0) ** 2
    log_m, log_v, n_blks, used, pts = [], [], [], [], []
    for m in scales:
        n_blk = int(arr.size) // int(m)
        if n_blk < min_blks:
            continue
        mat = arr[:n_blk * m].reshape(n_blk, m)
        valid_blk = np.mean(np.isfinite(mat), axis=1) >= 0.75
        nb_ok = int(np.count_nonzero(valid_blk))
        if nb_ok < min_blks:
            continue
        means = np.nanmean(mat[valid_blk], axis=1)
        v = float(np.var(means, ddof=1))
        if v > 1e-30:
            log_m.append(math.log10(float(m)))
            log_v.append(math.log10(v))
            n_blks.append(nb_ok)
            used.append(int(m))
            pts.append({"m_sectors": int(m), "n_blocks": nb_ok, "var": v})
    if len(log_m) < 3:
        return None
    lm = np.asarray(log_m, dtype=np.float64)
    lv = np.asarray(log_v, dtype=np.float64)
    slope, ic = np.polyfit(lm, lv, 1)
    hurst_ols = 1.0 + 0.5 * float(slope)
    pred = slope * lm + ic
    ss_res = float(np.sum((lv - pred) ** 2))
    ss_tot = float(np.sum((lv - np.mean(lv)) ** 2))
    r2 = 1.0 - ss_res / max(ss_tot, 1e-30)
    df = np.maximum(np.asarray(n_blks, dtype=np.float64) - 1.0, 1.0)
    w = (ln10_sq * df) / 2.0
    Sw = float(np.sum(w))
    Sx = float(np.sum(w * lm))
    Sy = float(np.sum(w * lv))
    Sxx = float(np.sum(w * lm * lm))
    Sxy = float(np.sum(w * lm * lv))
    det = Sw * Sxx - Sx * Sx
    h_wls = se_naive = None
    if det > 1e-30:
        sl_w = (Sw * Sxy - Sx * Sy) / det
        se_naive = 0.5 * math.sqrt(max(Sw / det, 0.0))
        h_wls = 1.0 + 0.5 * sl_w
    return {
        "hurst_H": hurst_ols, "hurst_H_wls": h_wls, "se_naive": se_naive,
        "r2": r2, "used_scales": used, "pts": pts,
    }


def thm_D8_multiscale_hurst_audit(disk_scan, sector_bp, euch_mask=None):
    """[THM-D8/D9 v1.5] Multi-scale WLS Hurst across m in {1,2,4,8,16,32,64} sectors
    with fixed-scale-grid Delete-1 Macro-Block Jackknife SE & Student-t 95% CI (zero RNG)."""
    tr = disk_scan["tracks"]
    scales = [1, 2, 4, 8, 16, 32, 64]
    out = {"sector_bp": int(sector_bp), "channels": {}}
    for ch in ("gc", "W_pp", "W_mp", "W_mp_unb", "W_mm", "I_nonbr_MM", "Omega_cond_excess"):
        arr = np.asarray(tr[ch], dtype=np.float64)
        if euch_mask is not None:
            arr = np.where(euch_mask, arr, np.nan)
        full = _hurst_wls_on_scales(arr, scales, min_blks=8)
        if full is None:
            out["channels"][ch] = {
                "hurst_H": None, "hurst_H_wls": None,
                "hurst_se_jackknife": None, "hurst_se_wls_naive": None,
                "jackknife_n_blocks": 0, "n_scales_used": 0, "used_scales": [],
                "hurst_ci95_lo": None, "hurst_ci95_hi": None,
                "above_white_noise_95": None, "loglog_r2": None, "points": [],
            }
            continue
        for pt in full["pts"]:
            pt["scale_kb"] = int(pt["m_sectors"] * sector_bp // 1000)
        h_wls = full["hurst_H_wls"]
        se_naive = full["se_naive"]
        used = full["used_scales"]
        se_jk = None
        b_eff = 0
        vpos = np.flatnonzero(np.isfinite(arr))
        if h_wls is not None and vpos.size >= 40:
            B = min(10, max(5, int(vpos.size) // max(used[-1], 16)))
            edges = np.linspace(0, int(vpos.size), B + 1, dtype=int)
            h_loo = []
            for b in range(B):
                a_del = arr.copy()
                a_del[vpos[edges[b]:edges[b + 1]]] = np.nan
                rb = _hurst_wls_on_scales(a_del, used, min_blks=6)
                if rb is not None and rb["hurst_H_wls"] is not None and len(rb["used_scales"]) == len(used):
                    h_loo.append(rb["hurst_H_wls"])
            if len(h_loo) >= 4:
                hl = np.asarray(h_loo, dtype=np.float64)
                b_eff = int(hl.size)
                se_jk = math.sqrt(((b_eff - 1.0) / float(b_eff)) * float(np.sum((hl - float(np.mean(hl))) ** 2)))
        if se_jk is not None and se_naive is not None:
            se_use = max(se_jk, se_naive)
        else:
            se_use = se_jk if (se_jk is not None and se_jk > 1e-12) else se_naive
        t_mult = t_crit_95(b_eff - 1) if (se_jk is not None and b_eff >= 4) else 1.959963984540054
        ci_lo = (h_wls - t_mult * se_use) if (h_wls is not None and se_use is not None) else None
        ci_hi = (h_wls + t_mult * se_use) if (h_wls is not None and se_use is not None) else None
        above_wn = bool(ci_lo > 0.50) if ci_lo is not None else None
        out["channels"][ch] = {
            "hurst_H": full["hurst_H"], "hurst_H_wls": h_wls,
            "hurst_se_jackknife": se_jk, "hurst_se_wls_naive": se_naive,
            "se_floor_applied": True,
            "jackknife_n_blocks": b_eff, "n_scales_used": len(used), "used_scales": used,
            "hurst_ci95_lo": ci_lo, "hurst_ci95_hi": ci_hi,
            "above_white_noise_95": above_wn,
            "loglog_r2": full["r2"], "points": full["pts"],
        }
    return out


MICRO_LAGS_DEFAULT = [3, 6, 9, 10, 11, 15, 30, 90, 147, 171]


def analyze_chromosome(key, meta, path, args):
    t0 = time.time()
    print(f"\n[{key}] {meta['label']} ({meta['acc']}) | file={Path(path).name}")
    fasta_sha = sha256_file(path)
    nt, hdr = load_chromosome_nt(path)
    L = int(nt.size)
    n_valid = int(np.count_nonzero(nt < 4))
    if n_valid < 5000:
        print(f"  [skip] too few valid bases ({n_valid})")
        return None
    cnt4 = np.bincount(nt[nt < 4], minlength=4).astype(np.float64)
    gc_global = float((cnt4[1] + cnt4[2]) / n_valid)
    gc_skew_global = float((cnt4[2] - cnt4[1]) / max(cnt4[2] + cnt4[1], 1.0))
    at_skew_global = float((cnt4[0] - cnt4[3]) / max(cnt4[0] + cnt4[3], 1.0))
    sector_bp = args.sector_bp if args.sector_bp > 0 else int(meta["default_sector_bp"])
    print(f"  L={L:,} bp | valid={n_valid:,} ({100.0*n_valid/L:.2f}%) | "
          f"GC={100*gc_global:.2f}% | sector={sector_bp//1000} kb | sha256={(fasta_sha or '')[:12]}")

    lags = [int(x) for x in args.micro_lags.split(",") if x.strip()]
    if 3 not in lags:
        lags = [3] + lags
    lag_rows = []
    B3_global = None
    B_lags = {}
    for d in lags:
        B = count_triplet_pairs_lag(nt, d=d)
        B_lags[int(d)] = B
        if d == 3:
            B3_global = B
        kv = strip_private(klein_v4_pr2_audit(B))
        sc = two_sided_screen_full(B)
        row = {"lag_bp": int(d), "klein": kv, "screen": sc}
        if d in (3, 30, 90, 147, 171):
            Bw, Eh = sector_het_expected_counts(nt, sector_bp=sector_bp, d=d)
            kv_eh = strip_private(klein_v4_pr2_audit(Eh))
            skl = screen_kl_chain_rule(Bw, Eh, df_eff=sc["df_eff"])
            sc_eh = two_sided_screen_full(Eh)
            row["het_null"] = {
                "E_het_W_mp": kv_eh["norm2_W_mp"],
                "E_het_W_mm": kv_eh["norm2_W_mm"],
                "E_het_A_time": kv_eh["A_time_L1"],
                "E_het_p_pal": kv_eh["p_pal_6bp"],
                "E_het_p_diag": kv_eh["p_diag_homodimer"],
                "p_pal_enrichment": kv["p_pal_6bp"] / max(kv_eh["p_pal_6bp"], TINY),
                "p_diag_enrichment": kv["p_diag_homodimer"] / max(kv_eh["p_diag_homodimer"], TINY),
                "I_nonbridge_het": sc_eh["I_nonbridge"],
                "I_nonbridge_intra_MM": sc["I_nonbridge_MM"] - sc_eh["I_nonbridge"],
                "noise_floor_D14": klein_noise_floor_from_null(Eh, kv["N"]),
            }
            row["screen_kl_vs_het"] = skl
        lag_rows.append(row)

    r3 = next(row for row in lag_rows if row["lag_bp"] == 3)
    kv3, sc3 = r3["klein"], r3["screen"]
    skl3 = r3.get("screen_kl_vs_het") or {}
    print(f"  [D1/PR2 d=3] W_++={kv3['norm2_W_pp']:.2e} W_+-={kv3['norm2_W_pm']:.2e} "
          f"W_-+={kv3['norm2_W_mp']:.2e} W_--={kv3['norm2_W_mm']:.2e} | "
          f"PR2_dyn={fnum(100*(kv3['PR2_dyn_share'] or 0),2)}% | "
          f"Time_PR2ok={fnum(100*(kv3['time_PR2_compliant_share'] or 0),2)}%")
    print(f"  [D2/D3  d=3] E_het(W_-++W_--)={r3['het_null']['E_het_W_mp']+r3['het_null']['E_het_W_mm']:.2e} | "
          f"I_br={sc3['I_bridge']:.5f} I_nonbr_MM={sc3['I_nonbridge_MM']:.5f} | "
          f"K_nonbr_MM={fnum(skl3.get('K_nonbridge_MM'),5)} | "
          f"|lead-wob|={sc3['lead_wobble_PR2_asym']:.2e}")

    fft_ac = exact_gc_autocorr_fft(nt, max_lag=600)
    bowker_D14 = bowker_klein_chi2_test(count_triplet_pairs_disjoint(nt, d=3))
    if bowker_D14 is not None:
        bowker_D14["stride_sweep"] = bowker_stride_sweep(nt, d=3)
    kv3_full = klein_v4_pr2_audit(B3_global)
    W_obs_mp, W_obs_mm = kv3_full["_W_mp"], kv3_full["_W_mm"]
    n_obs_mp = float(np.sum(W_obs_mp * W_obs_mp))
    n_obs_mm = float(np.sum(W_obs_mm * W_obs_mm))
    n_b3_tot = max(float(B3_global.sum()), 1.0)
    d15_ws = tuple(wm for wm in (200, 1000, 10000, 100000) if n_valid >= 2 * wm + 20)
    d15_raw = compositional_current_multiscale(nt, windows=d15_ws, d=3, trip_budget=4_000_000) if d15_ws else {}
    d15_sweep = {}
    for wm, (E_loc, m_ex) in d15_raw.items():
        if float(E_loc.sum()) <= 0.0:
            continue
        kv_loc = klein_v4_pr2_audit(E_loc)
        W_cur_mp, W_cur_mm = kv_loc["_W_mp"], kv_loc["_W_mm"]
        n_cur_mp = float(np.sum(W_cur_mp * W_cur_mp))
        n_cur_mm = float(np.sum(W_cur_mm * W_cur_mm))
        n_res_mp = float(np.sum((W_obs_mp - W_cur_mp) ** 2))
        n_res_mm = float(np.sum((W_obs_mm - W_cur_mm) ** 2))
        cos_mp = (float(np.sum(W_obs_mp * W_cur_mp)) / math.sqrt(n_obs_mp * n_cur_mp)) if (n_obs_mp > 1e-30 and n_cur_mp > 1e-30) else None
        cos_mm = (float(np.sum(W_obs_mm * W_cur_mm)) / math.sqrt(n_obs_mm * n_cur_mm)) if (n_obs_mm > 1e-30 and n_cur_mm > 1e-30) else None
        cov = float(m_ex) / n_b3_tot
        d15_sweep[f"M_{wm}"] = {
            "window_m": int(wm),
            "coverage_frac": cov,
            "coverage_ok": bool(cov >= 0.95),
            "E_loc_W_mp": n_cur_mp,
            "E_loc_W_mm": n_cur_mm,
            "W_mp_current_ratio": n_cur_mp / max(n_obs_mp, TINY),
            "W_mm_current_ratio": n_cur_mm / max(n_obs_mm, TINY),
            "W_mp_residual_share": n_res_mp / max(n_obs_mp, TINY),
            "W_mm_residual_share": n_res_mm / max(n_obs_mm, TINY),
            "cos_W_mp_obs_vs_cur": cos_mp,
            "cos_W_mm_obs_vs_cur": cos_mm,
        }
    d15_1k = d15_sweep.get("M_1000") or {}
    if bowker_D14:
        print(f"  [D14/D15]    Bowker(unadj) T log10p={fnum(bowker_D14['bowker_T_log10p_unadj'],1)} "
              f"R log10p={fnum(bowker_D14['bowker_R_log10p_unadj'],1)} | "
              f"D15(M=1k) W_res={fnum(100.0*_fval(d15_1k.get('W_mp_residual_share'), np.nan),2)}% "
              f"(cos={fnum(d15_1k.get('cos_W_mp_obs_vs_cur'),3,True)}, cov={fnum(100.0*_fval(d15_1k.get('coverage_frac'), np.nan),1)}%)")
    if fft_ac:
        print(f"  [D6/FFT]     Pitch={fft_ac['helix_peak_period_nt']:.1f} nt "
              f"(SNR={fft_ac['helix_peak_snr']:.2f}; Euk={fft_ac['euk_band_snr']:.2f} vs "
              f"Prok={fft_ac['prok_band_snr']:.2f}) | beat_3bp={fnum(fft_ac['beat_3bp'],5,True)} | "
              f"sat171={fnum(fft_ac['sat171_peak_excess'],5,True)}")

    disk = scan_chromosome_sectors(nt, sector_bp=sector_bp, topology=meta["topology"])
    ex_gc = disk["extrema_gc_skew_walk"]
    ex_pm = disk["extrema_W_pm_walk"]
    ex_mm = disk["extrema_W_mm_walk"]
    print(f"  [D5/DISK]    sectors={disk['n_valid_sectors']}/{disk['n_sectors']} | "
          f"SVD1(W_+-)={100*disk['svd_W_pm_var1']:.1f}% (r_GCsk={fnum(disk['corr_u_W_pm_vs_gc_skew'],3,True)}) | "
          f"SVD1(W_--)={100*disk['svd_W_mm_var1']:.1f}% (r_pm_mm={fnum(disk['corr_u_W_pm_vs_W_mm'],3,True)})")
    print(f"  [D5/WALKS]   GCsk min/max={ex_gc['min_mb']:.2f}/{ex_gc['max_mb']:.2f} Mb | "
          f"W_+- min/max={ex_pm['min_mb']:.2f}/{ex_pm['max_mb']:.2f} Mb | "
          f"W_-- min/max={ex_mm['min_mb']:.2f}/{ex_mm['max_mb']:.2f} Mb")

    # 4. Theorem D7: Unified Parameter-Free Partition + Euchromatin 10-lag, 1-Mb HAC & D15(M=1000)
    h_arg = getattr(args, "euch_h_thresh", 0.0)
    h_opt = float(h_arg) if (h_arg is not None and float(h_arg) > 0.0) else None
    d7 = thm_D7_euchromatin_vs_satellite_audit(
        nt, disk, sector_bp=sector_bp, topology=meta["topology"],
        h_thresh=h_opt, micro_lags=lags)
    e_fp = (d7["euchromatin"].get("fft_pitch") or {})
    s_fp = (d7["satellite_firmware"].get("fft_pitch") or {})
    ew = d7["euch_svd_walk"]
    print(f"  [D7/PART]    method={d7['partition_method']} | I_nonbr*={fnum(d7.get('inb_threshold_otsu'),4)} "
          f"(dBIC={fnum(d7.get('bic_2vs1_delta'),1)}) | Euch={100*d7['euch_fraction']:.1f}% ({d7['n_euch_sectors']} sec) | "
          f"Sat={100*d7['sat_fraction']:.1f}% ({d7['n_sat_sectors']} sec, {len(d7['satellite_blocks_mb'])} block(s))")
    if not d7["euchromatin"].get("skipped"):
        ek = d7["euchromatin"]["klein_d3"]
        esk = d7["euchromatin"]["screen_kl_d3"]
        print(f"  [D7/EUCH]    W_-+={ek['norm2_W_mp']:.2e} Time_PR2ok={fnum(100*(ek.get('time_PR2_compliant_share_unb') or 0),2)}% | "
              f"K_nonbr_MM={fnum(esk.get('K_nonbridge_MM'),5)} | "
              f"Pitch={fnum(e_fp.get('helix_peak_period_nt'),1)} nt (SNR={fnum(e_fp.get('helix_peak_snr'),2)}) | "
              f"unitSVD1(W_+-)={100*ew['svd1_W_pm_unit']:.1f}% (r_GCsk={fnum(ew['corr_W_pm_vs_gc_skew'],3,True)})")
    if not d7["satellite_firmware"].get("skipped"):
        sk_ = d7["satellite_firmware"]["klein_d3"]
        ssk = d7["satellite_firmware"]["screen_kl_d3"]
        print(f"  [D7/SAT]     W_-+={sk_['norm2_W_mp']:.2e} W_--={sk_['norm2_W_mm']:.2e} | "
              f"K_nonbr_MM={fnum(ssk.get('K_nonbridge_MM'),5)} | "
              f"Pitch={fnum(s_fp.get('helix_peak_period_nt'),1)} nt | sat171={fnum(s_fp.get('sat171_peak_excess'),5,True)}")

    # 5. Theorem D8: Multi-Scale WLS Hurst with Delete-1 Macro-Block Jackknife (Unified THM-D7 mask)
    H3_arr = np.asarray(disk["tracks"]["H_3mer_bits"], dtype=np.float64)
    INB_arr = np.asarray(disk["tracks"]["I_nonbr_MM"], dtype=np.float64)
    INB_RAW_arr = (np.asarray(disk["tracks"]["D_5lead"], dtype=np.float64)
                   + np.asarray(disk["tracks"]["D_3wob"], dtype=np.float64)
                   + np.asarray(disk["tracks"]["D_1_3"], dtype=np.float64)) if "D_5lead" in disk["tracks"] else None
    if h_opt is not None:
        euch_mask = np.isfinite(H3_arr) & np.isfinite(INB_arr) & (H3_arr >= h_opt)
    else:
        euch_mask, _, _ = partition_euch_sat_parameter_free(H3_arr, INB_arr, INB_RAW=INB_RAW_arr)
    d8_all = thm_D8_multiscale_hurst_audit(disk, sector_bp=sector_bp, euch_mask=None)
    d8_euch = thm_D8_multiscale_hurst_audit(disk, sector_bp=sector_bp, euch_mask=euch_mask)
    c_e = d8_euch["channels"]
    cu_e = c_e["W_mp_unb"]
    print(f"  [D8/HURST]   Euch WLS(JK): GC={fnum(c_e['gc']['hurst_H_wls'],3)} | "
          f"W_++={fnum(c_e['W_pp']['hurst_H_wls'],3)} | "
          f"W_-+_unb={fnum(cu_e['hurst_H_wls'],3)} [{fnum(cu_e['hurst_ci95_lo'],2)},{fnum(cu_e['hurst_ci95_hi'],2)}] "
          f"(B={cu_e['jackknife_n_blocks']}, scales={cu_e['n_scales_used']}) | "
          f"NonBr={fnum(c_e['I_nonbr_MM']['hurst_H_wls'],3)}")

    # 6. Theorem D18: Exact Delayed-Pruning PELT L2 Segmentation (MBIC = 3*ln n)
    pos_arr = np.asarray(disk["tracks"]["pos_mb"], dtype=np.float64)
    p_h3 = exact_pelt_l2(H3_arr, min_seg=3)
    wmp_u = np.asarray(disk["tracks"].get("W_mp_unb", disk["tracks"]["W_mp"]), dtype=np.float64)
    lwmp = np.where(np.isfinite(wmp_u) & (wmp_u > 0.0), np.log10(np.maximum(wmp_u, 1e-10)), np.nan)
    p_wmp = exact_pelt_l2(lwmp, min_seg=3)
    if p_h3 is not None:
        p_h3["changepoints_mb"] = [float(pos_arr[min(c, pos_arr.size - 1)]) for c in p_h3["changepoints_sector"]]
    if p_wmp is not None:
        p_wmp["changepoints_mb"] = [float(pos_arr[min(c, pos_arr.size - 1)]) for c in p_wmp["changepoints_sector"]]
    cc_18 = pelt_cross_channel_concordance(
        p_h3["changepoints_sector"] if p_h3 else [],
        p_wmp["changepoints_sector"] if p_wmp else [],
        int(H3_arr.size), tol=2,
    )
    d18_pelt = {"H_3mer": p_h3, "log10_W_mp_unb": p_wmp, "cross_channel_concordance": cc_18}

    e_d11 = (d7["euchromatin"].get("arrow_D11") or {})
    e_ww = (d7["euchromatin"].get("ww_nucleosome_fft") or {})
    cl_euch = d7["euchromatin"].get("cluster_1mb_euch_D14") or {}
    d15_eu = d7["euchromatin"].get("comp_current_euch_M1000") or {}
    lv_e = (d7.get("euch_svd_walk") or {}).get("levy_D17_euch") or disk.get("levy_D17") or {}
    if e_d11:
        print(f"  [D13-D18]    NonEnv3600={fnum(100.0*_fval(e_d11.get('pure_nonenv_arrow_share'), np.nan),1)}% | "
              f"Euch-1Mb HAC t={fnum(cl_euch.get('cluster_t_hac'),2)} (sign_p={_fsci(cl_euch.get('sign_p_twosided'))}) | "
              f"Euch-D15 Res={fnum(100.0*_fval(d15_eu.get('W_mp_residual_share'), np.nan),2)}% (cos={fnum(d15_eu.get('cos_W_mp_obs_vs_cur'),2,True)}) | "
              f"WW={fnum(e_ww.get('ww_nuc_band_nt'),1)}nt (Δφ={fnum(e_ww.get('ww_ss_phase_diff_deg'),1)}°, p_F={_fsci(e_ww.get('ww_peak_p_adj'))}) | "
              f"PELT={p_h3['n_segments'] if p_h3 else 'n/a'}/{p_wmp['n_segments'] if p_wmp else 'n/a'} | Lévy p={fnum(lv_e.get('levy_p_shift'),3)}")

    b_strata = disk.pop("_B_strata", None)
    disk.pop("_Wpm_rows", None); disk.pop("_Wmm_rows", None); disk.pop("_valid_sec_idx", None)
    b3_euch = d7["euchromatin"].pop("_B3_part", None)
    b_euch_lags = d7["euchromatin"].pop("_B_lags_part", None)
    b3_sat = d7["satellite_firmware"].pop("_B3_part", None)
    d7["satellite_firmware"].pop("_B_lags_part", None)
    P3 = B3_global / max(float(B3_global.sum()), TINY)
    elapsed = time.time() - t0
    print(f"  [DONE] {key} in {elapsed:.1f} s")

    return {
        "key": key,
        "meta": meta,
        "header": hdr,
        "fasta_file": Path(path).name,
        "fasta_sha256": fasta_sha,
        "length_bp": L,
        "valid_bp": n_valid,
        "valid_frac": n_valid / float(L),
        "gc_global": gc_global,
        "gc_skew_global": gc_skew_global,
        "at_skew_global": at_skew_global,
        "sector_bp": int(sector_bp),
        "micro_lags": lag_rows,
        "fft_autocorr": fft_ac,
        "disk_scan": disk,
        "theorem_D7": d7,
        "theorem_D8": {"all_sectors": d8_all, "euchromatin_only": d8_euch},
        "theorem_D14_bowker": bowker_D14,
        "theorem_D14_cluster_1mb": disk.get("cluster_1mb_D14"),
        "theorem_D14_cluster_1mb_euch": cl_euch or None,
        "theorem_D15_comp_current": d15_sweep,
        "theorem_D15_comp_current_euch": d15_eu or None,
        "theorem_D18_pelt": d18_pelt,
        "runtime_seconds": elapsed,
        "_P3_global": P3,
        "_B3_all": B3_global,
        "_B3_euch": b3_euch,
        "_B3_sat": b3_sat,
        "_B_strata": b_strata,
        "_B_lags": B_lags,
        "_B_euch_lags": b_euch_lags,
    }


# =====================================================================
# 9. Deterministic Zero-RNG Self-Tests (Theorems D1-D18)
# =====================================================================
def self_test(verbose=False):
    ok = {}

    # 1. Klein V_4 exact reverse-complement equivariance: B^(d)(revcomp(s)) == R(B^(d)(s))
    L = 12000
    t = np.arange(L, dtype=np.int64)
    nt = ((t * 7 + (t // 3) * 5 + (t * t) // 11) % 4).astype(np.uint8)
    nt[5000:5020] = 255
    rc = reverse_complement_nt(nt)
    max_rc_err = 0.0
    for d in (3, 7, 11, 30):
        B = count_triplet_pairs_lag(nt, d=d, chunk_size=2400)
        Brc = count_triplet_pairs_lag(rc, d=d, chunk_size=2400)
        RB = B.T[np.ix_(MIRROR, MIRROR)]
        max_rc_err = max(max_rc_err, float(np.abs(Brc - RB).max()))
    assert max_rc_err == 0.0, f"revcomp != R(B): {max_rc_err}"
    ok["thmD1_revcomp_equals_R_max_err"] = max_rc_err

    # 2. Exact PR2-symmetric fixture with NON-ZERO W_-+ time arrow & D_5'lead == D_3'wobble
    nt_pr2 = np.concatenate([nt, np.full(50, 255, dtype=np.uint8), rc])
    B_pr2 = count_triplet_pairs_lag(nt_pr2, d=3, chunk_size=3000)
    kv_pr2 = klein_v4_pr2_audit(B_pr2)
    sc_pr2 = two_sided_screen_full(B_pr2)
    assert kv_pr2["pythagorean_err"] < 1e-15
    assert kv_pr2["PR2_split_err"] < 1e-15
    assert kv_pr2["time_split_err"] < 1e-15
    assert kv_pr2["norm2_W_pm"] < 1e-30 and kv_pr2["norm2_W_mm"] < 1e-30
    assert kv_pr2["norm2_W_mp"] > 1e-5
    assert sc_pr2["four_term_err"] < 1e-14
    assert sc_pr2["lead_wobble_PR2_asym"] < 1e-14
    ok["thmD1_PR2_exact_Wpm_Wmm_zero"] = max(kv_pr2["norm2_W_pm"], kv_pr2["norm2_W_mm"])
    ok["thmD1_PR2_nonzero_Wmp_time_arrow"] = kv_pr2["norm2_W_mp"]
    ok["thmD3_PR2_lead_equals_wobble_err"] = sc_pr2["lead_wobble_PR2_asym"]

    # 3. Theorem D2: Sector-heterogeneity null E_het^(d) strictly annihilates W_-+ and W_--
    Bw, Eh = sector_het_expected_counts(nt, sector_bp=2000, d=3)
    kv_eh = klein_v4_pr2_audit(Eh)
    assert kv_eh["norm2_W_mp"] < 1e-30 and kv_eh["norm2_W_mm"] < 1e-30
    assert kv_eh["A_time_L1"] < 1e-15
    skl = screen_kl_chain_rule(Bw, Eh)
    assert skl["all_nonnegative"] and skl["fwd_err"] < 1e-14 and skl["rev_err"] < 1e-14
    ok["thmD2_Ehet_directed_annihilation"] = kv_eh["norm2_W_mp"] + kv_eh["norm2_W_mm"]
    ok["thmD3_screen_kl_chain_max_err"] = max(skl["fwd_err"], skl["rev_err"])

    # 4. Theorem D4: Frame-phase lock chain rule (both positive coding fixture & exact null)
    n_cod = 2000
    nt_cod = np.empty(n_cod * 3, dtype=np.uint8)
    c_prev = 0
    for i in range(n_cod):
        b1 = (c_prev + (i % 2)) % 4
        b2 = (i * 3 + b1) % 4
        b3 = (i * 5 + 2 * b1 + 1) % 4
        nt_cod[3 * i] = b1
        nt_cod[3 * i + 1] = b2
        nt_cod[3 * i + 2] = b3
        c_prev = (b1 + b2) % 4
    om = frame_phase_lock_sector(nt_cod, d=3)
    assert om["chain_err"] < 1e-14 and om["Omega_cond"] > 0.05, om
    nt_null = np.concatenate([np.tile(np.array([0, 1, 2, 3], dtype=np.uint8), 1500),
                              np.array([0, 1, 2, 3, 0], dtype=np.uint8)])
    om_null = frame_phase_lock_sector(nt_null, d=3)
    assert om_null["Omega_phase"] < 1e-14 and om_null["Omega_cond"] < 1e-14, om_null
    ok["thmD4_frame_lock_chain_err"] = om["chain_err"]
    ok["thmD4_synthetic_coding_Omega_cond"] = om["Omega_cond"]
    ok["thmD4_phase_invariant_null_Omega_phase"] = om_null["Omega_phase"]

    # 5. Theorem D5: Synthetic circular chromosome with exact oriC at 0.015 Mb and ter at 0.045 Mb
    motif_fwd = np.array([2, 2, 2, 0, 0, 2, 2, 1, 0, 2, 0, 3, 2, 2, 0], dtype=np.uint8)
    motif_rev = reverse_complement_nt(motif_fwd)
    arm_a = np.tile(motif_fwd, 2000)
    arm_b = np.tile(motif_rev, 2000)
    syn_chr = np.concatenate([arm_b[:15000], arm_a, arm_b[15000:]])
    dsk = scan_chromosome_sectors(syn_chr, sector_bp=3000, topology="circular")
    assert dsk["svd_W_pm_var1"] > 0.99 and dsk["svd_W_mm_var1"] > 0.99
    pm_pair = sorted([dsk["extrema_W_pm_walk"]["min_mb"], dsk["extrema_W_pm_walk"]["max_mb"]])
    mm_pair = sorted([dsk["extrema_W_mm_walk"]["min_mb"], dsk["extrema_W_mm_walk"]["max_mb"]])
    assert abs(pm_pair[0] - 0.015) <= 0.004 and abs(pm_pair[1] - 0.045) <= 0.004, pm_pair
    assert abs(mm_pair[0] - 0.015) <= 0.004 and abs(mm_pair[1] - 0.045) <= 0.004, mm_pair
    ok["thmD5_synthetic_oriC_ter_Wpm_min_max_mb"] = [
        dsk["extrema_W_pm_walk"]["min_mb"], dsk["extrema_W_pm_walk"]["max_mb"]]
    ok["thmD5_synthetic_oriC_ter_Wmm_min_max_mb"] = [
        dsk["extrema_W_mm_walk"]["min_mb"], dsk["extrema_W_mm_walk"]["max_mb"]]

    # 6. Theorem D6: Exact FFT autocorrelation vs direct autocorrelation + 10.5-bp peak
    idx_f = np.arange(30000, dtype=np.float64)
    syn_105 = np.where(np.cos(2.0 * math.pi * idx_f / 10.5) > 0, 2, 0).astype(np.uint8)
    syn_105[12000:12030] = 255
    fft_r = exact_gc_autocorr_fft(syn_105, max_lag=180, chunk_size=8000)
    v_f = syn_105 < 4
    mu_f = float(np.mean((syn_105[v_f] == 1) | (syn_105[v_f] == 2)))
    var_f = mu_f * (1.0 - mu_f)
    x_f = np.where(v_f, ((syn_105 == 1) | (syn_105 == 2)).astype(np.float64) - mu_f, 0.0)
    max_fft_err = 0.0
    for dd in (1, 10, 21, 105, 171):
        ok_d = v_f[:-dd] & v_f[dd:]
        direct = float(np.mean(x_f[:-dd][ok_d] * x_f[dd:][ok_d]) / var_f)
        max_fft_err = max(max_fft_err, abs(fft_r["acf_values"][dd - 1] - direct))
    assert max_fft_err < 1e-12, max_fft_err
    assert 10.4 <= fft_r["helix_peak_period_nt"] <= 10.6
    assert fft_r["euk_band_snr"] > 30.0
    ok["thmD6_fft_vs_direct_max_err"] = max_fft_err
    ok["thmD6_synthetic_10_5bp_peak"] = fft_r["helix_peak_period_nt"]

    # 7. Theorem D7 & D8: Synthetic hybrid chromosome (euchromatin + low-entropy satellite block)
    sat_unit = np.array([2, 2, 0, 0, 3], dtype=np.uint8)
    euch_x = np.arange(30000, dtype=np.uint64) + np.uint64(1)
    euch_z = (euch_x ^ (euch_x >> np.uint64(30))) * np.uint64(0xbf58476d1ce4e5b9)
    euch_z = (euch_z ^ (euch_z >> np.uint64(27))) * np.uint64(0x94d049bb133111eb)
    euch_seq = ((euch_z ^ (euch_z >> np.uint64(31))) & np.uint64(3)).astype(np.uint8)
    sat_seq = np.tile(sat_unit, 3000)
    hyb_chr = np.concatenate([euch_seq[:15000], sat_seq, euch_seq[15000:]])
    dsk_h = scan_chromosome_sectors(hyb_chr, sector_bp=3000, topology="linear")
    d7_h = thm_D7_euchromatin_vs_satellite_audit(hyb_chr, dsk_h, sector_bp=3000, topology="linear")
    assert d7_h["n_sat_sectors"] == 5 and d7_h["n_euch_sectors"] == 10, (d7_h["n_sat_sectors"], d7_h["n_euch_sectors"])
    assert len(d7_h["satellite_blocks_mb"]) == 1, d7_h["satellite_blocks_mb"]
    d8_h = thm_D8_multiscale_hurst_audit(dsk_h, sector_bp=3000)
    ok["thmD7_synthetic_euch_sat_sectors"] = [d7_h["n_euch_sectors"], d7_h["n_sat_sectors"]]
    ok["thmD7_synthetic_sat_block_mb"] = d7_h["satellite_blocks_mb"][0]
    ok["thmD8_multiscale_channels_ok"] = (len(d8_h["channels"]) == 7)

    assert abs(kv_pr2["norm2_W_pm_unb"]) < 1e-4 and abs(kv_pr2["norm2_W_mm_unb"]) < 1e-4
    M_t11 = B_pr2 / B_pr2.sum()
    R_of_Pi = (_DS_PROJ @ M_t11 @ _DO_PROJ).T[np.ix_(MIRROR, MIRROR)]
    Pi_of_R = _DS_PROJ @ (M_t11.T[np.ix_(MIRROR, MIRROR)]) @ _DO_PROJ
    comm_err = float(np.max(np.abs(R_of_Pi - Pi_of_R)))
    assert comm_err < 1e-15, comm_err
    d11_t = thm_D11_arrow_screen_and_cpg_audit(B_pr2)
    assert d11_t["pythagorean_err"] < 1e-15 and d11_t["nonbridge_arrow_share"] > 0.5
    ww_syn = np.where(np.cos(2.0 * math.pi * idx_f / 10.3) > 0, 0, 2).astype(np.uint8)
    ww_r = exact_ww_autocorr_fft(ww_syn, max_lag=160, chunk_size=8000)
    assert 10.1 <= ww_r["ww_nuc_band_nt"] <= 10.5 and ww_r["ww_nuc_band_snr"] > 20.0, ww_r
    ok["thmD9_unbiased_PR2_zero_ok"] = True
    ok["thmD11_commutator_R_Pi_err"] = comm_err
    ok["thmD11_Wmp_screen_pythag_err"] = d11_t["pythagorean_err"]
    ok["thmD12_WW_synthetic_10_3bp_peak"] = ww_r["ww_nuc_band_nt"]

    # 8. Theorem D13: V_4-closed screen envelope Pi_env (rank 496; [Pi_env, T]=[Pi_env, R]=[Pi_env, S]=0)
    M_d13 = np.arange(NC * NC, dtype=np.float64).reshape(NC, NC) / (NC * NC)
    M_d13 = np.sin(M_d13 * 13.7) + 0.3 * np.cos(M_d13.T * 7.3)
    pe = pi_env_project(M_d13)
    pe2 = pi_env_project(pe)
    cT = float(np.max(np.abs(pi_env_project(M_d13.T) - pe.T)))
    cR = float(np.max(np.abs(pi_env_project(M_d13.T[np.ix_(MIRROR, MIRROR)]) - pe.T[np.ix_(MIRROR, MIRROR)])))
    cS = float(np.max(np.abs(pi_env_project(M_d13[np.ix_(MIRROR, MIRROR)]) - pe[np.ix_(MIRROR, MIRROR)])))
    assert float(np.max(np.abs(pe2 - pe))) < 1e-14 and max(cT, cR, cS) < 1e-14
    tr_I = tr_T = tr_R = tr_S = 0.0
    for u_ in range(NC):
        for v_ in range(NC):
            tr_I += _DS_PROJ[u_, u_] * _DO_PROJ[v_, v_] + _DO_PROJ[u_, u_] * _DS_PROJ[v_, v_] - _DC2_PROJ[u_, u_] * _DC2_PROJ[v_, v_]
            tr_T += _DS_PROJ[v_, u_] * _DO_PROJ[v_, u_] + _DO_PROJ[v_, u_] * _DS_PROJ[v_, u_] - _DC2_PROJ[v_, u_] * _DC2_PROJ[v_, u_]
            ru_, rv_ = int(MIRROR[v_]), int(MIRROR[u_])
            tr_R += _DS_PROJ[ru_, u_] * _DO_PROJ[v_, rv_] + _DO_PROJ[ru_, u_] * _DS_PROJ[v_, rv_] - _DC2_PROJ[ru_, u_] * _DC2_PROJ[v_, rv_]
            su_, sv_ = int(MIRROR[u_]), int(MIRROR[v_])
            tr_S += _DS_PROJ[su_, u_] * _DO_PROJ[v_, sv_] + _DO_PROJ[su_, u_] * _DS_PROJ[v_, sv_] - _DC2_PROJ[su_, u_] * _DC2_PROJ[v_, sv_]
    env_dims = [int(round(0.25 * (tr_I + tr_T + tr_R + tr_S))),
                int(round(0.25 * (tr_I + tr_T - tr_R - tr_S))),
                int(round(0.25 * (tr_I - tr_T + tr_R - tr_S))),
                int(round(0.25 * (tr_I - tr_T - tr_R + tr_S)))]
    assert env_dims == [132, 118, 130, 116] and d11_t["env_pythagorean_err"] < 1e-15 and d11_t["env_quadrant_leak_err"] < 1e-15
    ok["thmD13_Pi_env_V4_comm_err"] = max(cT, cR, cS)
    ok["thmD13_Pi_env_quadrant_dims"] = env_dims

    # 9. Theorems D14-D17: Bowker symmetry, E_loc exact mass, WW-SS antiphase & Fisher g, Levy shuffle
    bk_sym = bowker_klein_chi2_test(B_pr2 + B_pr2.T)
    assert bk_sym["bowker_T_chi2"] == 0.0 and bk_sym["bowker_R_chi2"] == 0.0
    E_l1, m_l1 = compositional_current_null(nt, window_m=200, d=3, max_win_chunk=7)
    E_l2, m_l2 = compositional_current_null(nt, window_m=200, d=3, max_win_chunk=1000)
    assert abs(float(E_l1.sum()) - m_l1) < 1e-9 and float(np.max(np.abs(E_l1 - E_l2))) < 1e-11
    fg_ex, fp_ex = fisher_g_pvalue([0.5] + [0.5 / 9.0] * 9)
    assert abs(fp_ex - 10.0 * (0.5 ** 9)) < 1e-14 and ww_r["ww_ss_phase_diff_deg"] > 170.0
    lv_t = levy_area_signed_test(np.sin(np.linspace(0, 8, 50)), np.cos(np.linspace(0, 8, 50)))
    assert lv_t["shuffle_err"] < 1e-13 and lv_t["abs_levy_z"] > 1.0
    ok["thmD14_bowker_sym_zero"] = 0.0
    ok["thmD15_Eloc_mass_err"] = abs(float(E_l1.sum()) - m_l1)
    ok["thmD16_WW_SS_antiphase_deg"] = ww_r["ww_ss_phase_diff_deg"]
    ok["thmD16_fisher_g_exact_err"] = abs(fp_ex - 10.0 * (0.5 ** 9))
    ok["thmD17_levy_shuffle_err"] = lv_t["shuffle_err"]

    # 10. Patch 5B Zero-RNG tests (i)-(vi): D15 nested budget, D14 LR & pooling floor,
    #     SVD permutation invariance, D17.2 FFT circular-shift Levy, D16 3-bp negative vs positive
    ms_a = compositional_current_multiscale(nt, windows=(200, 1000), d=3, trip_budget=600)
    ms_b = compositional_current_multiscale(nt, windows=(200, 1000), d=3, trip_budget=50000)
    e1k_s, m1k_s = compositional_current_null(nt, window_m=1000, d=3, trip_budget=50000)
    d15_nest_err = max(
        float(np.max(np.abs(ms_a[200][0] - ms_b[200][0]))),
        float(np.max(np.abs(ms_a[1000][0] - ms_b[1000][0]))),
        float(np.max(np.abs(ms_a[1000][0] - e1k_s))),
    )
    assert d15_nest_err < 1e-10 and abs(ms_a[1000][1] - m1k_s) < 1e-9, d15_nest_err
    cov_1k = m1k_s / max(float(count_triplet_pairs_lag(nt, d=3).sum()), 1.0)
    assert 0.80 <= cov_1k <= 1.0 + 1e-12, cov_1k

    lr_50k = chi2_tail_log10p(50000.0, 2016)
    lr_k = chi2_tail_log10p(2016.0, 2016)
    assert math.isfinite(lr_50k) and -9025.0 < lr_50k < -9010.0, lr_50k
    lr_max_err = abs(lr_k - math.log10(0.5 - 1.0 / (3.0 * math.sqrt(math.pi * 2016.0))))
    try:
        from scipy.stats import chi2 as _sp_chi2
        for _xv in (1800.0, 2016.0, 2200.0, 2500.0, 3500.0):
            lr_max_err = max(lr_max_err, abs(chi2_tail_log10p(_xv, 2016) - (_sp_chi2.logsf(_xv, 2016) / math.log(10.0))))
        assert lr_max_err < 1e-4, lr_max_err
    except ImportError:
        pass

    fl_unif = klein_noise_floor_from_null(np.ones((NC, NC), dtype=np.float64), 100000.0)
    fl_err_rel = abs(fl_unif["W_mp"] - 1.0 / (4.0 * 100000.0)) / (1.0 / (4.0 * 100000.0))
    assert fl_err_rel < 0.01, fl_err_rel
    cl_t = cluster_1mb_arrow_test([B_pr2 * (1.0 + 0.05 * j_) for j_ in range(8)])
    assert cl_t is not None and cl_t["n_pos"] == 8 and cl_t["sign_p_twosided"] < 0.01, cl_t

    _x_half = np.sin(np.arange(1, 30 * 50 + 1, dtype=np.float64).reshape(30, 50) * 0.73)
    _X_tie = np.hstack([_x_half, _x_half])
    _u_orig, _, _k1 = _top_svd_intrinsic(_X_tie)
    _perm = np.arange(30)[::-1]
    _u_perm, _, _k2 = _top_svd_intrinsic(_X_tie[_perm])
    svd_perm_err = float(np.max(np.abs(_u_orig[_perm] - _u_perm)))
    assert _k1 == _k2 and svd_perm_err < 1e-11, (svd_perm_err, _k1, _k2)

    _t5 = np.linspace(0.0, 10.0, 64, endpoint=False)
    _a5 = np.sin(_t5) + 0.3 * np.cos(2.0 * _t5)
    _b5 = np.cos(_t5) - 0.2 * np.sin(3.0 * _t5)
    _ac5, _bc5 = _a5 - float(np.mean(_a5)), _b5 - float(np.mean(_b5))
    _ca5 = np.concatenate([[0.0], np.cumsum(_ac5[:-1])])
    _A_fft5 = np.fft.irfft(np.fft.rfft(_ca5 + 0.5 * _ac5) * np.conj(np.fft.rfft(_bc5)), n=64)
    _A_brk5 = np.array([
        0.5 * float(np.sum(_ca5 * np.roll(_bc5, k_) - np.concatenate([[0.0], np.cumsum(np.roll(_bc5, k_)[:-1])]) * _ac5))
        for k_ in range(64)
    ])
    levy_fft_err = float(np.max(np.abs(_A_fft5 - _A_brk5)))
    assert levy_fft_err < 1e-11, levy_fft_err

    _l600 = np.arange(1, 601, dtype=np.float64)
    _noise600 = 2e-4 * np.sin(_l600 * _l600 * 0.37)
    _acf_neg = 0.04 * ((_l600 / 4.0) ** (-0.25)) + 0.01 * np.cos(2.0 * math.pi * _l600 / 3.0) + _noise600
    _acf_pos = _acf_neg + 4e-4 * np.cos(2.0 * math.pi * _l600 / 10.3) * np.exp(-_l600 / 150.0)
    _, _, _, _, _p_neg = notched_fisher_and_band_peak_p(_acf_neg, 10.3, spec_lim=200)
    _, _, _, _, _p_pos = notched_fisher_and_band_peak_p(_acf_pos, 10.3, spec_lim=200)
    assert _p_neg > 0.1 and _p_pos < 1e-6, (_p_neg, _p_pos)

    # Newey-West HAC check: positively autocorrelated scores -> se_hac > se_iid
    _nw_ac = newey_west_t_stats(np.sin(np.linspace(0.0, 4.0 * math.pi, 40)) + 2.0)
    assert _nw_ac is not None and _nw_ac["se_hac"] > _nw_ac["se_iid"], _nw_ac

    # Parameter-free THM-D7 checks: 3-mode (euch+alpha-sat+HSat) & pure bimodal Alu euchromatin (0 FP)
    _h_3m = np.concatenate([np.full(200, 5.82), np.full(20, 5.66), np.full(20, 4.80)])
    _i_3m = np.concatenate([0.02 + 0.01 * np.sin(np.arange(200)), np.full(20, 0.45), np.full(20, 0.90)])
    _eu3, _sa3, _in3 = partition_euch_sat_parameter_free(_h_3m, _i_3m)
    assert int(_eu3.sum()) == 200 and int(_sa3.sum()) == 40 and _in3["entropy_inversion_ok"], _in3
    _h_alu = np.concatenate([np.full(120, 5.78), np.full(80, 5.86)])
    _i_alu = np.concatenate([0.025 + 0.005 * np.sin(np.arange(120)), 0.075 + 0.01 * np.sin(np.arange(80))])
    _eua, _saa, _ina = partition_euch_sat_parameter_free(_h_alu, _i_alu)
    assert int(_saa.sum()) == 0 and (not _ina["entropy_inversion_ok"]), _ina

    # THM-D8 Jackknife check on deterministic white noise vs persistent track
    _ix_jk = np.arange(500, dtype=np.uint64) + np.uint64(1)
    _zz_jk = (_ix_jk ^ (_ix_jk >> np.uint64(30))) * np.uint64(0xbf58476d1ce4e5b9)
    _zz_jk = (_zz_jk ^ (_zz_jk >> np.uint64(27))) * np.uint64(0x94d049bb133111eb)
    _wn_jk = ((_zz_jk ^ (_zz_jk >> np.uint64(31))) & np.uint64(0xFFFFFF)).astype(np.float64) / float(1 << 24)
    _ps_jk = np.convolve(_wn_jk, np.ones(25) / 25.0, mode="same") + 0.1 * np.sin(np.arange(500) / 40.0)
    _d8_syn = thm_D8_multiscale_hurst_audit(
        {"tracks": {"gc": _ps_jk, "W_pp": _ps_jk, "W_mp": _ps_jk, "W_mp_unb": _wn_jk,
                    "W_mm": _wn_jk, "I_nonbr_MM": _ps_jk, "Omega_cond_excess": _wn_jk}},
        sector_bp=100000,
    )
    assert _d8_syn["channels"]["W_mp_unb"]["above_white_noise_95"] is False, _d8_syn["channels"]["W_mp_unb"]
    assert _d8_syn["channels"]["W_mp"]["above_white_noise_95"] is True, _d8_syn["channels"]["W_mp"]

    ok["thmD15_nested_budget_max_err"] = d15_nest_err
    ok["thmD14_lugannani_rice_max_err"] = lr_max_err
    ok["thmD14_lemma_pooling_floor_rel_err"] = fl_err_rel
    ok["thmD14_newey_west_hac_gt_iid"] = True
    ok["thmD5_svd_perm_invariance_err"] = svd_perm_err
    ok["thmD17_fft_vs_roll_max_err"] = levy_fft_err
    ok["thmD16_curved_3bp_neg_vs_pos_p_adj"] = [_p_neg, _p_pos]
    ok["thmD7_param_free_3mode_and_alu_null_ok"] = True
    ok["thmD8_jackknife_wn_vs_pers_ok"] = True

    # 11. Theorem D18: Exact delayed-pruning PELT L2 == Exact O(n^2) DP (30 adversarial tests + step + flat null)
    _y_step = np.concatenate([np.full(300, 5.82), np.full(60, 4.70), np.full(240, 5.82)]) + 0.02 * np.sin(np.arange(600) * 0.7)
    _r_step = exact_pelt_l2(_y_step, min_seg=3)
    assert _r_step["changepoints_sector"] == [300, 360] and _r_step["r2_explained"] > 0.99, _r_step
    _ix_fl = np.arange(300, dtype=np.uint64) + np.uint64(1)
    _zz_fl = (_ix_fl ^ (_ix_fl >> np.uint64(30))) * np.uint64(0xbf58476d1ce4e5b9)
    _zz_fl = (_zz_fl ^ (_zz_fl >> np.uint64(27))) * np.uint64(0x94d049bb133111eb)
    _u_fl = ((_zz_fl ^ (_zz_fl >> np.uint64(31))) & np.uint64(0xFFFFFF)).astype(np.float64) / float(1 << 24)
    _r_flat = exact_pelt_l2(_u_fl, min_seg=3)
    assert _r_flat["n_segments"] == 1, _r_flat
    _pelt_dp_max_err = 0.0
    for _tr in range(30):
        _nn = 12 + (_tr % 17)
        _ms = (1, 2, 3, 5)[_tr % 4]
        _is = np.arange(_nn, dtype=np.uint64) + np.uint64(_tr * 97 + 1)
        _hz = (_is ^ (_is >> np.uint64(30))) * np.uint64(0xbf58476d1ce4e5b9)
        _hz = (_hz ^ (_hz >> np.uint64(27))) * np.uint64(0x94d049bb133111eb)
        _vv = np.array([-10.0, 0.0, 10.0])[((_hz ^ (_hz >> np.uint64(31))) % np.uint64(3)).astype(np.int64)]
        _pt = 8.0 + float(_tr % 7) * 5.0
        _rp = exact_pelt_l2(_vv, beta=_pt, min_seg=_ms)
        _f_dp, _ = _exact_dp_l2_z((_vv - float(np.mean(_vv))) / _rp["sigma_eff"], _pt, min_seg=_ms)
        _pelt_dp_max_err = max(_pelt_dp_max_err, abs(_rp["opt_cost"] - _f_dp))
    assert _pelt_dp_max_err < 1e-9, _pelt_dp_max_err
    ok["thmD18_pelt_exact_cps"] = _r_step["changepoints_sector"]
    ok["thmD18_pelt_flat_null_k1"] = (_r_flat["n_segments"] == 1)
    ok["thmD18_pelt_vs_exact_dp_30trials_max_err"] = _pelt_dp_max_err

    # 12. v1.5.2 Zero-RNG Audit Tests: THM-D19 mod-5 aliasing, D7 clipping-atom null,
    #     D18 persistent-null + concordance, D14 zero-cluster null, D15 M=10k, legacy post-audit
    _hs5 = np.tile(np.array([2, 2, 0, 0, 3], dtype=np.uint8), 4000)
    _k10 = klein_v4_pr2_audit(count_triplet_pairs_lag(_hs5, d=10))
    _k11 = klein_v4_pr2_audit(count_triplet_pairs_lag(_hs5, d=11))
    assert _k10["norm2_W_mp"] < 1e-30 and _k10["norm2_W_mm"] < 1e-30
    assert _k11["norm2_W_mp"] > 1e-4
    assert abs(_k10["norm2_W_pp"] + _k10["norm2_W_pm"] - _k10["norm2_total"]) < 1e-14
    ok["thmD19_mod5_HSat_aliasing"] = [_k10["norm2_W_mp"], _k11["norm2_W_mp"]]

    _inb_atom = np.concatenate([np.full(200, -0.002), 0.001 + 0.002 * np.abs(np.sin(np.arange(300)))])
    _h3_atom = np.concatenate([np.full(200, 5.85), np.full(300, 5.80)])
    _, _sa_atom, _inf_atom = partition_euch_sat_parameter_free(_h3_atom, _inb_atom)
    assert int(_sa_atom.sum()) == 0 and _inf_atom["frac_clipped"] >= 0.35, _inf_atom
    ok["thmD7_clipping_atom_null_rejected"] = True

    _cc_t = pelt_cross_channel_concordance([300, 360], [301, 360], 600, tol=2)
    assert _cc_t["n_matched"] == 2 and _cc_t["jaccard"] == 1.0 and _cc_t["hypergeom_p"] < 1e-3, _cc_t
    _r_ps_raw = exact_pelt_l2(_ps_jk, min_seg=3, lrv_inflate=False)
    _r_ps_lrv = exact_pelt_l2(_ps_jk, min_seg=3, lrv_inflate=True)
    assert _r_ps_lrv["n_segments"] < _r_ps_raw["n_segments"], (_r_ps_lrv["n_segments"], _r_ps_raw["n_segments"])
    ok["thmD18_lrv_reduces_persistent_splits"] = [_r_ps_raw["n_segments"], _r_ps_lrv["n_segments"]]
    ok["thmD18_concordance_exact_ok"] = True

    _B_sym0 = B_pr2 + B_pr2.T
    _cl_zero = cluster_1mb_arrow_test([_B_sym0 for _ in range(8)])
    assert _cl_zero is not None and abs(_cl_zero["mean_proj"]) < 1e-15 and _cl_zero["n_pos"] == 0, _cl_zero
    ok["thmD14_cluster_1mb_zero_null_ok"] = True

    _nt25 = np.tile(nt[:5000], 5)
    _ms10_a = compositional_current_multiscale(_nt25, windows=(1000, 10000), d=3, trip_budget=5000)
    _ms10_b = compositional_current_multiscale(_nt25, windows=(1000, 10000), d=3, trip_budget=50000)
    _e10_s, _ = compositional_current_null(_nt25, window_m=10000, d=3, trip_budget=50000)
    _err10 = max(float(np.max(np.abs(_ms10_a[10000][0] - _ms10_b[10000][0]))),
                 float(np.max(np.abs(_ms10_a[10000][0] - _e10_s))))
    assert _err10 < 1e-9, _err10
    ok["thmD15_multiscale_10k_max_err"] = _err10

    _leg_pay = {"results": [{"sector_bp": 3000, "disk_scan": dsk_h, "theorem_D7": {"partition_method": "legacy_fixed_h_thresh", "h_threshold_bits": 5.65}}]}
    apply_v15_post_audit(_leg_pay)
    ok["thmD7_legacy_post_audit_ok"] = True

    # Deterministic chr2 1%-centromere imbalance test (2000 wide euch + 20 compact cen)
    _ix2 = np.arange(2000, dtype=np.float64)
    _i2_e = 0.036 * np.exp(0.35 * np.sin(_ix2 * 0.13) + 0.15 * np.cos(_ix2 * 0.41))
    _h2_e = 5.76 + 0.15 * (np.log10(_i2_e) - math.log10(0.036)) + 0.01 * np.sin(_ix2 * 0.7)
    _i2_c = 0.35 + 0.15 * np.abs(np.sin(np.arange(20) * 0.5))
    _h2_c = 5.64 + 0.02 * np.cos(np.arange(20) * 0.9)
    _eu2, _sa2, _in2 = partition_euch_sat_parameter_free(
        np.concatenate([_h2_e, _h2_c]), np.concatenate([_i2_e, _i2_c]) - 0.01872,
        INB_RAW=np.concatenate([_i2_e, _i2_c]),
    )
    assert int(_sa2.sum()) == 20 and int(_eu2.sum()) == 2000 and _in2["entropy_inversion_ok"], _in2
    ok["thmD7_chr2_1pct_centromere_imbalance_ok"] = True

    ok.update(self_test_virome_onco(verbose=False))
    ok["pilot_catalog_count"] = len(PILOT_10_CHROMOSOMES)
    assert len(HUMAN_T2T_24) == 24
    ok["all_passed"] = True
    if verbose:
        print(json.dumps(_clean(ok), indent=2))
    return ok


def _exact_dp_l2_z(z, pen, min_seg=3):
    """Reference O(n^2) exact dynamic programming L2 solver for verifying PELT optimality."""
    n = int(z.size)
    ps1 = np.concatenate([[0.0], np.cumsum(z)])
    ps2 = np.concatenate([[0.0], np.cumsum(z * z)])
    F = np.full(n + 1, np.inf, dtype=np.float64)
    F[0] = -pen
    cp_prev = np.zeros(n + 1, dtype=np.int32)
    for t in range(min_seg, n + 1):
        cand = np.arange(0, t - min_seg + 1, dtype=np.int32)
        lens = (t - cand).astype(np.float64)
        s1 = ps1[t] - ps1[cand]
        s2 = ps2[t] - ps2[cand]
        sse = np.maximum(s2 - (s1 * s1) / lens, 0.0)
        vals = F[cand] + sse + pen
        bi = int(np.argmin(vals))
        F[t] = float(vals[bi])
        cp_prev[t] = int(cand[bi])
    cps = []
    cur = n
    while cur > 0:
        p = int(cp_prev[cur])
        if p > 0:
            cps.append(p)
        cur = p
    cps.reverse()
    return float(F[n]), cps


def _nw_lrv_ratio(v):
    x = np.asarray(v, dtype=np.float64)
    K = int(x.size)
    if K < 8:
        return 1.0
    xc = x - float(np.mean(x))
    g0 = float(np.dot(xc, xc)) / float(K)
    if g0 <= 1e-30:
        return 1.0
    L = min(max(1, int(math.floor(K ** (1.0 / 3.0)))), K - 2)
    lrv = g0
    for lag in range(1, L + 1):
        w = 1.0 - float(lag) / float(L + 1)
        lrv += 2.0 * w * float(np.dot(xc[:-lag], xc[lag:])) / float(K)
    return max(1.0, lrv / g0)


def pelt_cross_channel_concordance(cps_a, cps_b, n_sec, tol=2):
    """[THM-D18 v1.5.2] Exact tolerance-matched (+/- tol sectors) Jaccard & hypergeometric
    p-value quantifying boundary agreement between two independent PELT channels."""
    a = sorted(set(int(c) for c in (cps_a or []) if 0 < int(c) < n_sec))
    b = sorted(set(int(c) for c in (cps_b or []) if 0 < int(c) < n_sec))
    na, nb = len(a), len(b)
    if n_sec <= 2 or na == 0 or nb == 0:
        return {"n_a": na, "n_b": nb, "tol_sectors": int(tol), "n_matched": 0, "jaccard": 0.0, "hypergeom_p": 1.0}
    dil_b = np.zeros(n_sec, dtype=bool)
    for c in b:
        dil_b[max(1, c - tol):min(n_sec, c + tol + 1)] = True
    K_dil = int(np.count_nonzero(dil_b[1:]))
    M = n_sec - 1
    k_match = sum(1 for c in a if dil_b[c])
    jacc = float(k_match) / float(max(na + nb - k_match, 1))
    k_hi = min(na, K_dil)
    denom = math.comb(M, na)
    num = sum(math.comb(K_dil, j) * math.comb(M - K_dil, na - j)
              for j in range(k_match, k_hi + 1) if (na - j) <= (M - K_dil))
    p_hyp = min(1.0, num / denom) if denom > 0 else 1.0
    return {
        "n_a": na, "n_b": nb, "tol_sectors": int(tol),
        "n_matched": int(k_match), "jaccard": jacc, "hypergeom_p": p_hyp,
    }


def exact_pelt_l2(y, beta=None, min_seg=3, lrv_inflate=True):
    """[THM-D18 v1.5.2] Exact Zero-RNG PELT L2 solver with valid delayed pruning at
    u = t - min_seg + 1 (Killick, Fearnhead & Eckley 2012, JASA Thm 3.1 with K=0 for
    additive SSE; provably identical to exact O(n^2) DP for any min_seg >= 1),
    Bartlett LRV persistence inflation sigma_eff = sigma_mad * sqrt(LRV/gamma_0),
    and Zhang-Siegmund (2007, Biometrics) MBIC penalty beta = 3 * ln(n)."""
    a = np.asarray(y, dtype=np.float64)
    ok = np.isfinite(a)
    if int(np.count_nonzero(ok)) < 2 * min_seg + 2:
        return None
    x = np.arange(a.size, dtype=np.float64)
    v = np.interp(x, x[ok], a[ok])
    n = int(v.size)
    dv = np.diff(v)
    mad = float(np.median(np.abs(dv - np.median(dv))))
    sig_mad = mad / (0.67448975 * math.sqrt(2.0))
    sig_fl = float(np.std(v)) * 0.05
    sig = max(sig_mad, sig_fl, 1e-12)
    lrv_r = _nw_lrv_ratio(v) if lrv_inflate else 1.0
    sig_eff = sig * math.sqrt(lrv_r)
    z = (v - float(np.mean(v))) / sig_eff
    pen = float(beta) if beta is not None else 3.0 * math.log(float(n))

    ps1 = np.concatenate([[0.0], np.cumsum(z)])
    ps2 = np.concatenate([[0.0], np.cumsum(z * z)])
    F = np.full(n + 1, np.inf, dtype=np.float64)
    F[0] = -pen
    cp_prev = np.zeros(n + 1, dtype=np.int32)
    R = np.array([0], dtype=np.int32)

    for t in range(min_seg, n + 1):
        cand = R[t - R >= min_seg]
        if cand.size > 0:
            lens = (t - cand).astype(np.float64)
            s1 = ps1[t] - ps1[cand]
            s2 = ps2[t] - ps2[cand]
            sse = np.maximum(s2 - (s1 * s1) / lens, 0.0)
            vals = F[cand] + sse + pen
            best_i = int(np.argmin(vals))
            F[t] = float(vals[best_i])
            cp_prev[t] = int(cand[best_i])
        u = t - min_seg + 1
        if u >= min_seg and math.isfinite(F[u]):
            pc = R[u - R >= min_seg]
            if pc.size > 0:
                l_u = (u - pc).astype(np.float64)
                s1_u = ps1[u] - ps1[pc]
                s2_u = ps2[u] - ps2[pc]
                sse_u = np.maximum(s2_u - (s1_u * s1_u) / l_u, 0.0)
                kept = pc[(F[pc] + sse_u) <= F[u] + 1e-12]
                recent = R[u - R < min_seg]
                R = np.concatenate([kept, recent, np.array([t], dtype=np.int32)])
            else:
                R = np.append(R, np.int32(t))
        else:
            R = np.append(R, np.int32(t))

    cps = []
    cur = n
    while cur > 0:
        p = int(cp_prev[cur])
        if p > 0:
            cps.append(p)
        cur = p
    cps.reverse()
    bounds = [0] + cps + [n]
    sse_seg = sum(float(np.sum((v[b0:b1] - np.mean(v[b0:b1])) ** 2)) for b0, b1 in zip(bounds[:-1], bounds[1:]))
    sse_tot = float(np.sum((v - np.mean(v)) ** 2))
    return {
        "n_sectors": n,
        "n_interpolated_nans": int(np.count_nonzero(~ok)),
        "n_segments": len(cps) + 1,
        "changepoints_sector": cps,
        "sigma_mad": float(sig),
        "sigma_floor_active": bool(sig_mad < sig_fl),
        "lrv_ratio": float(lrv_r),
        "sigma_eff": float(sig_eff),
        "lrv_calibrated": bool(lrv_inflate),
        "penalty_beta": float(pen),
        "penalty_rule": "MBIC_3ln_n_LRV" if (beta is None and lrv_inflate) else ("MBIC_3ln_n" if beta is None else "custom"),
        "opt_cost": float(F[n]),
        "r2_explained": 1.0 - sse_seg / max(sse_tot, 1e-30),
    }


def apply_v15_post_audit(payload):
    """[v1.5.0] Unified Zero-RNG post-audit (no hardcoded 5.70/0.08 override):
      - Ensures THM-D8 Jackknife Hurst uses the exact unified THM-D7 euchromatin mask.
      - Recomputes THM-D18 exact delayed-pruning PELT L2 (MBIC = 3*ln n)."""
    for r in payload.get("results", []):
        disk = r.get("disk_scan") or {}
        tr = disk.get("tracks") or {}
        if not tr or "H_3mer_bits" not in tr:
            continue
        H3 = np.asarray(tr["H_3mer_bits"], dtype=np.float64)
        INB = np.asarray(tr["I_nonbr_MM"], dtype=np.float64)
        INB_RAW = (np.asarray(tr["D_5lead"], dtype=np.float64)
                   + np.asarray(tr["D_3wob"], dtype=np.float64)
                   + np.asarray(tr["D_1_3"], dtype=np.float64)) if "D_5lead" in tr else None
        pos = np.asarray(tr["pos_mb"], dtype=np.float64)
        sec_bp = int(r.get("sector_bp", disk.get("sector_bp", 100000)))
        n_sec = int(H3.size)
        d7 = r.setdefault("theorem_D7", {})
        if d7.get("partition_method") == "legacy_fixed_h_thresh" and _finite(d7.get("h_threshold_bits")):
            is_euch = np.isfinite(H3) & np.isfinite(INB) & (H3 >= float(d7["h_threshold_bits"]))
        else:
            is_euch, is_sat, pf_info = partition_euch_sat_parameter_free(H3, INB, INB_RAW=INB_RAW)
            if "partition_method" not in d7:
                d7.update(pf_info)
                n_val_i = max(1, int(np.count_nonzero(np.isfinite(H3) & np.isfinite(INB))))
                d7["euch_fraction"] = float(np.count_nonzero(is_euch)) / float(n_val_i)
                d7["sat_fraction"] = 1.0 - d7["euch_fraction"]
                sat_blks = []
                half_mb = 0.5 * sec_bp / 1e6
                in_b, b0 = False, 0
                for w in range(n_sec):
                    if is_sat[w] and not in_b:
                        in_b, b0 = True, w
                    elif (not is_sat[w]) and in_b:
                        if w - b0 >= 2:
                            sat_blks.append([float(pos[b0] - half_mb), float(pos[w - 1] + half_mb)])
                        in_b = False
                if in_b and n_sec - b0 >= 2:
                    sat_blks.append([float(pos[b0] - half_mb), float(pos[n_sec - 1] + half_mb)])
                d7["satellite_blocks_mb"] = sat_blks

        if "frac_clipped" not in d7:
            pos_i = INB[np.isfinite(INB) & (INB > 0.0)]
            fl_i = max(1e-4, 0.05 * float(np.median(pos_i))) if pos_i.size else 1e-4
            d7["frac_clipped"] = float(np.mean(INB[np.isfinite(INB)] <= fl_i)) if np.any(np.isfinite(INB)) else 0.0
        if "min_euch_h3_bits" not in d7 and _finite(d7.get("h_threshold_bits")):
            d7["min_euch_h3_bits"] = float(d7["h_threshold_bits"])

        d8 = r.setdefault("theorem_D8", {})
        eu_ch = ((d8.get("euchromatin_only") or {}).get("channels") or {}).get("W_mp_unb") or {}
        if not eu_ch.get("se_floor_applied"):
            d8["euchromatin_only"] = thm_D8_multiscale_hurst_audit(disk, sector_bp=sec_bp, euch_mask=is_euch)

        p18 = r.get("theorem_D18_pelt") or {}
        if (p18.get("H_3mer") or {}).get("lrv_calibrated") is not True:
            p_h3 = exact_pelt_l2(H3, min_seg=3, lrv_inflate=True)
            wmp_u = np.asarray(tr.get("W_mp_unb", tr["W_mp"]), dtype=np.float64)
            lwmp = np.where(np.isfinite(wmp_u) & (wmp_u > 0.0), np.log10(np.maximum(wmp_u, 1e-10)), np.nan)
            p_wmp = exact_pelt_l2(lwmp, min_seg=3, lrv_inflate=True)
            if p_h3 is not None:
                p_h3["changepoints_mb"] = [float(pos[min(c, n_sec - 1)]) for c in p_h3["changepoints_sector"]]
            if p_wmp is not None:
                p_wmp["changepoints_mb"] = [float(pos[min(c, n_sec - 1)]) for c in p_wmp["changepoints_sector"]]
            cc_18 = pelt_cross_channel_concordance(
                p_h3["changepoints_sector"] if p_h3 else [],
                p_wmp["changepoints_sector"] if p_wmp else [],
                n_sec, tol=2,
            )
            r["theorem_D18_pelt"] = {"H_3mer": p_h3, "log10_W_mp_unb": p_wmp, "cross_channel_concordance": cc_18}
    return payload


apply_sprint6_post_audit = apply_v15_post_audit


def run_euch_2d_multilag_control(payload, chr_dir, chr_list_str, lags=None, h2d_thresh=0.0, inb2d_max=0.0):
    """Optional fast per-chromosome pure-euchromatin 10-lag control on cached chromosomes
    (uses parameter-free THM-D7 mask when h2d_thresh <= 0)."""
    if not chr_list_str:
        return None
    if lags is None:
        lags = MICRO_LAGS_DEFAULT
    chr_dir = Path(chr_dir)
    by_lbl = {r["meta"]["label"]: r for r in payload.get("results", [])}
    by_key = {r["key"]: r for r in payload.get("results", [])}
    sel = [s.strip() for s in chr_list_str.split(",") if s.strip()]
    use_pf = (h2d_thresh <= 0.0)
    out = {
        "partition_method": "log10_inb_raw_kittler_bic_v1.5.3" if use_pf else "custom_2d",
        "h2d_thresh": None if use_pf else float(h2d_thresh),
        "inb2d_max": None if use_pf else float(inb2d_max),
        "chromosomes": {},
    }
    for raw in sel:
        r = by_lbl.get(raw) or by_key.get(raw) or by_key.get(f"Human_T2T_{raw}")
        if r is None:
            continue
        k, meta = r["key"], r["meta"]
        fp = chr_dir / f"{k}_{meta['acc']}.fna.gz"
        if not fp.exists() and chr_dir.exists():
            hits = list(chr_dir.glob(f"*{meta['acc']}*.fna*"))
            if hits:
                fp = hits[0]
        if not fp.exists():
            print(f"  [EUCH-LAG] Skipping {meta['label']}: FASTA not found in {chr_dir}")
            continue
        tr = r["disk_scan"]["tracks"]
        H3 = np.asarray(tr["H_3mer_bits"], dtype=np.float64)
        INB = np.asarray(tr["I_nonbr_MM"], dtype=np.float64)
        INB_RAW = (np.asarray(tr["D_5lead"], dtype=np.float64)
                   + np.asarray(tr["D_3wob"], dtype=np.float64)
                   + np.asarray(tr["D_1_3"], dtype=np.float64)) if "D_5lead" in tr else None
        sec_bp = int(r.get("sector_bp", 100000))
        if use_pf:
            mask2d, _, _ = partition_euch_sat_parameter_free(H3, INB, INB_RAW=INB_RAW)
        else:
            mask2d = np.isfinite(H3) & np.isfinite(INB) & (H3 >= float(h2d_thresh)) & (INB <= float(inb2d_max))
        print(f"  [EUCH-LAG] Loading {meta['label']} ({fp.name}) | pure euchromatin = {int(np.count_nonzero(mask2d))}/{H3.size} sectors...")
        nt, _ = load_chromosome_nt(fp)
        L = int(nt.size)
        nt_pure = np.full(L, 255, dtype=np.uint8)
        for w in np.flatnonzero(mask2d):
            s0 = int(w) * sec_bp
            s1 = min(L, s0 + sec_bp)
            nt_pure[s0:s1] = nt[s0:s1]
        del nt
        rows = []
        for d in lags:
            B = count_triplet_pairs_lag(nt_pure, d=int(d))
            kv = strip_private(klein_v4_pr2_audit(B))
            sc = two_sided_screen_full(B)
            a11 = thm_D11_arrow_screen_and_cpg_audit(B)
            top1 = (a11.get("top_motifs_all") or [None])[0]
            rows.append({
                "lag_bp": int(d),
                "N": kv["N"],
                "norm2_W_mp_unb": kv["norm2_W_mp_unb"],
                "norm2_W_mm_unb": kv["norm2_W_mm_unb"],
                "time_PR2_compliant_share_unb": kv["time_PR2_compliant_share_unb"],
                "I_nonbridge_MM": sc["I_nonbridge_MM"],
                "pure_nonenv_arrow_share": a11["pure_nonenv_arrow_share"],
                "cpg_free_arrow_share": a11["cpg_free_arrow_share"],
                "top1_motif": top1["motif_6mer"] if top1 else None,
                "top1_share_pct": top1["flux_share_pct"] if top1 else None,
            })
        del nt_pure
        out["chromosomes"][meta["label"]] = {
            "n_pure_euch_sectors": int(np.count_nonzero(mask2d)),
            "n_total_sectors": int(H3.size),
            "lags": rows,
        }
    return out if out["chromosomes"] else None


def _fsci(x):
    return "n/a" if not _finite(x) else f"{float(x):.2e}"


def grand_genome_audit(results, b3_store):
    B_all = np.zeros((NC, NC), dtype=np.float64)
    B_euch = np.zeros((NC, NC), dtype=np.float64)
    B_sat = np.zeros((NC, NC), dtype=np.float64)
    n_strata = len(ENTROPY_STRATA_LABELS)
    B_strata_tot = np.zeros((n_strata, NC, NC), dtype=np.float64)
    strata_sec_tot = np.zeros(n_strata, dtype=np.int64)
    strata_neg_tot = np.zeros(n_strata, dtype=np.int64)
    B_lags_tot = {}
    B_euch_lags_tot = {}

    tot_bp = tot_valid = n_euch_sec = n_sat_sec = 0
    ww_pitches, ww_snrs, ww_ctr_snrs, ww_dphi = [], [], [], []
    h_raw_list, h_unb_list, h_wls_list, pal_enr_list, diag_enr_list = [], [], [], [], []
    above_wn_cnt = 0

    for r in results:
        k = r["key"]
        tot_bp += int(r["length_bp"]); tot_valid += int(r["valid_bp"])
        d7 = r.get("theorem_D7") or {}
        n_euch_sec += int(d7.get("n_euch_sectors", 0)); n_sat_sec += int(d7.get("n_sat_sectors", 0))
        if f"{k}__all" in b3_store: B_all += b3_store[f"{k}__all"]
        if f"{k}__euch" in b3_store: B_euch += b3_store[f"{k}__euch"]
        if f"{k}__sat" in b3_store: B_sat += b3_store[f"{k}__sat"]
        if f"{k}__strata" in b3_store:
            st_arr = np.asarray(b3_store[f"{k}__strata"], dtype=np.float64)
            if st_arr.shape == B_strata_tot.shape:
                B_strata_tot += st_arr
        sc_cnt = (r.get("disk_scan") or {}).get("entropy_strata_sector_counts")
        if sc_cnt and len(sc_cnt) == n_strata:
            strata_sec_tot += np.asarray(sc_cnt, dtype=np.int64)
        sn_cnt = (r.get("disk_scan") or {}).get("entropy_strata_neg_mp_counts")
        if sn_cnt and len(sn_cnt) == n_strata:
            strata_neg_tot += np.asarray(sn_cnt, dtype=np.int64)
        pref = f"{k}__lag_"
        epref = f"{k}__euch_lag_"
        for bk, bmat in b3_store.items():
            if bk.startswith(epref):
                try:
                    lag_d = int(bk[len(epref):])
                except ValueError:
                    continue
                if lag_d not in B_euch_lags_tot:
                    B_euch_lags_tot[lag_d] = np.zeros((NC, NC), dtype=np.float64)
                B_euch_lags_tot[lag_d] += np.asarray(bmat, dtype=np.float64)
            elif bk.startswith(pref):
                try:
                    lag_d = int(bk[len(pref):])
                except ValueError:
                    continue
                if lag_d not in B_lags_tot:
                    B_lags_tot[lag_d] = np.zeros((NC, NC), dtype=np.float64)
                B_lags_tot[lag_d] += np.asarray(bmat, dtype=np.float64)

        eu_ww = (d7.get("euchromatin") or {}).get("ww_nucleosome_fft") or {}
        if _finite(eu_ww.get("ww_nuc_band_nt")): ww_pitches.append(float(eu_ww["ww_nuc_band_nt"]))
        if _finite(eu_ww.get("ww_nuc_band_snr")): ww_snrs.append(float(eu_ww["ww_nuc_band_snr"]))
        if _finite(eu_ww.get("ww_minus_ss_nuc_snr")): ww_ctr_snrs.append(float(eu_ww["ww_minus_ss_nuc_snr"]))
        if _finite(eu_ww.get("ww_ss_phase_diff_deg")): ww_dphi.append(float(eu_ww["ww_ss_phase_diff_deg"]))
        d8_ch = ((r.get("theorem_D8") or {}).get("euchromatin_only") or {}).get("channels") or {}
        if _finite((d8_ch.get("W_mp") or {}).get("hurst_H")):
            h_raw_list.append(float(d8_ch["W_mp"]["hurst_H"]))
        if _finite((d8_ch.get("W_mp_unb") or {}).get("hurst_H")):
            h_unb_list.append(float(d8_ch["W_mp_unb"]["hurst_H"]))
        if _finite((d8_ch.get("W_mp_unb") or {}).get("hurst_H_wls")):
            h_wls_list.append(float(d8_ch["W_mp_unb"]["hurst_H_wls"]))
        if (d8_ch.get("W_mp_unb") or {}).get("above_white_noise_95") is True:
            above_wn_cnt += 1
        ml = r.get("micro_lags") or []
        if ml:
            hn0 = ml[0].get("het_null") or {}
            if _finite(hn0.get("p_pal_enrichment")): pal_enr_list.append(float(hn0["p_pal_enrichment"]))
            if _finite(hn0.get("p_diag_enrichment")): diag_enr_list.append(float(hn0["p_diag_enrichment"]))

    def _pack(B):
        if float(B.sum()) <= 0:
            return None
        return {
            "klein_d3": strip_private(klein_v4_pr2_audit(B)),
            "screen_d3": two_sided_screen_full(B),
            "arrow_D11": thm_D11_arrow_screen_and_cpg_audit(B),
        }

    entropy_phase_D10 = []
    for bi, lbl in enumerate(ENTROPY_STRATA_LABELS):
        Bb = B_strata_tot[bi]
        Nb = float(Bb.sum())
        if Nb <= 0:
            entropy_phase_D10.append({
                "bin_idx": bi, "label": lbl, "n_sectors": int(strata_sec_tot[bi]),
                "N": 0, "skipped": True,
            })
        else:
            n_sb = int(strata_sec_tot[bi])
            entropy_phase_D10.append({
                "bin_idx": bi, "label": lbl, "n_sectors": n_sb,
                "n_neg_u_mp": int(strata_neg_tot[bi]),
                "frac_negative_u_mp": (float(strata_neg_tot[bi]) / float(n_sb)) if n_sb > 0 else None,
                "N": int(Nb), "skipped": False,
                "klein_d3": strip_private(klein_v4_pr2_audit(Bb)),
                "screen_d3": two_sided_screen_full(Bb),
                "arrow_D11": thm_D11_arrow_screen_and_cpg_audit(Bb),
            })

    def _pack_multilag(b_map):
        rows = []
        for lag_d in sorted(b_map.keys()):
            Bl = b_map[lag_d]
            if float(Bl.sum()) <= 0:
                continue
            rows.append({
                "lag_bp": int(lag_d),
                "klein": strip_private(klein_v4_pr2_audit(Bl)),
                "screen": two_sided_screen_full(Bl),
                "arrow_D11": thm_D11_arrow_screen_and_cpg_audit(Bl),
            })
        return rows

    multilag_D12 = _pack_multilag(B_lags_tot)
    euch_multilag_D12 = _pack_multilag(B_euch_lags_tot)

    def _med_minmax(vals):
        if not vals:
            return None
        a = np.asarray(vals, dtype=np.float64)
        return {"n": int(a.size), "median": float(np.median(a)),
                "min": float(np.min(a)), "max": float(np.max(a))}

    return {
        "n_chromosomes": len(results), "total_bp": int(tot_bp), "total_gb": float(tot_bp) / 1e9,
        "total_valid_bp": int(tot_valid), "n_euch_sectors": int(n_euch_sec), "n_sat_sectors": int(n_sat_sec),
        "euch_sector_share": float(n_euch_sec) / max(1, n_euch_sec + n_sat_sec),
        "whole_genome": _pack(B_all), "euchromatin_grand": _pack(B_euch), "satellite_grand": _pack(B_sat),
        "pooling_floor_D14_2": {
            "whole_genome_1_over_4N": (1.0 / (4.0 * float(B_all.sum()))) if float(B_all.sum()) > 0 else None,
            "euchromatin_1_over_4N": (1.0 / (4.0 * float(B_euch.sum()))) if float(B_euch.sum()) > 0 else None,
        },
        "entropy_phase_diagram_D10": entropy_phase_D10,
        "grand_multilag_spectrum_D12": multilag_D12,
        "grand_euch_multilag_spectrum_D12": euch_multilag_D12,
        "grand_D9_D12_summary": {
            "ww_nuc_pitch_nt": _med_minmax(ww_pitches),
            "ww_nuc_snr": _med_minmax(ww_snrs),
            "ww_minus_ss_nuc_snr": _med_minmax(ww_ctr_snrs),
            "ww_ss_phase_diff_deg": _med_minmax(ww_dphi),
            "hurst_W_mp_raw": _med_minmax(h_raw_list),
            "hurst_W_mp_unb": _med_minmax(h_unb_list),
            "hurst_W_mp_unb_wls": _med_minmax(h_wls_list),
            "hurst_W_mp_above_wn_count": int(above_wn_cnt),
            "p_pal_enrichment_d3": _med_minmax(pal_enr_list),
            "p_diag_enrichment_d3": _med_minmax(diag_enr_list),
        },
    }


def cross_chromosome_matrix_human24(results, b3_store):
    keys = [r["key"] for r in results if f"{r['key']}__all" in b3_store]
    S = len(keys)
    if S < 2:
        return None
    by_k = {r["key"]: r for r in results}
    labels = [by_k[k]["meta"]["label"] for k in keys]
    arms = [by_k[k]["meta"].get("arm_class", "autosome") for k in keys]
    DS, DO = np.zeros((NC, NC)), np.zeros((NC, NC))
    for s in range(16):
        ix = np.flatnonzero(SUF16 == s); DS[np.ix_(ix, ix)] = 1.0 / ix.size
        iy = np.flatnonzero(OUT16 == s); DO[np.ix_(iy, iy)] = 1.0 / iy.size

    def _normP(B):
        s = float(B.sum())
        return (B / s) if s > 0 else np.zeros((NC, NC), dtype=np.float64)

    P_all = [_normP(b3_store[f"{k}__all"]) for k in keys]
    P_euch = [_normP(b3_store.get(f"{k}__euch", b3_store[f"{k}__all"])) for k in keys]

    def _cos_mat(mats):
        U = [M / max(float(np.linalg.norm(M, ord="fro")), TINY) for M in mats]
        D = np.zeros((S, S))
        for i in range(S):
            for j in range(i + 1, S):
                D[i, j] = D[j, i] = 1.0 - float(np.sum(U[i] * U[j]))
        return D

    D_stat_all = _cos_mat([0.5 * (P + P.T) for P in P_all])
    D_dir_all = _cos_mat([0.5 * (P - P.T) for P in P_all])
    D_dir_euch = _cos_mat([0.5 * (P - P.T) for P in P_euch])
    D_nonbr_euch = _cos_mat([P - DS @ P @ DO for P in P_euch])

    def _nn(D):
        out = []
        for i in range(S):
            row = D[i].copy(); row[i] = np.inf; j = int(np.argmin(row))
            out.append({"chr": labels[i], "nearest": labels[j], "dist": float(row[j])})
        return out

    return {
        "keys": keys, "labels": labels, "arm_classes": arms,
        "dist_static_all": D_stat_all.tolist(), "dist_directed_all": D_dir_all.tolist(),
        "dist_directed_euch": D_dir_euch.tolist(), "dist_nonbridge_euch": D_nonbr_euch.tolist(),
        "nn_directed_all": _nn(D_dir_all), "nn_directed_euch": _nn(D_dir_euch),
    }


def write_human24_markdown(payload, path):
    R = payload.get("results", [])
    gg = payload.get("grand_genome") or {}
    prov = payload.get("provenance") or {}
    L = [
        f"# Complete 24-Chromosome Human T2T-CHM13v2.0 Torus Disk Atlas (v{VERSION})", "",
        "100% Zero-RNG · Gapless Telomere-to-Telomere Audit (`NC_060925.1`–`NC_060948.1`) · Theorems D1–D19.", "",
        "> **Assembly Scope Note:** T2T-CHM13v2.0 combines the 46,XX CHM13 cell line (chr1–chr22, chrX) with T2T-chrY from HG002 (46,XY); acrocentric rDNA morph arrays include consensus-modeled copies in RefSeq v2.0.", "",
    ]
    if prov:
        L += [
            f"> **Provenance & Reproducibility:** Python `{prov.get('python_version')}` · NumPy `{prov.get('numpy_version')}` · SciPy `{prov.get('scipy_version')}` · `script_sha256={(prov.get('script_sha256') or 'n/a')[:16]}` · Partition=`{prov.get('partition_method')}` · PELT=`{prov.get('pelt_penalty')}` · UTC=`{prov.get('timestamp_utc')}`.", "",
        ]

    def _fpct(val, nd=1):
        return f"{100.0 * float(val):.{nd}f}%" if _finite(val) else "n/a"

    if gg.get("whole_genome"):
        L += [f"## 0. Exact Grand-Pooled Human Genome Summary ({gg['n_chromosomes']} Chromosomes, {gg['total_gb']:.3f} Gb)", "",
              "| Partition | Triplet Pairs N | ‖W_++‖² | ‖W_+-‖² | ‖W_-+‖²_unb | ‖W_--‖²_unb | PR2-ok (unb) % | Non-Env 3600-D % | NonBr 256-D % | CpG-free % | I_nonbr_MM | 6-bp Pal % |",
              "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
        for nm, blk in (("Whole Human Genome", gg["whole_genome"]),
                        ("Euchromatin (THM-D7 KI-BIC)", gg["euchromatin_grand"]),
                        ("Satellite/Tandem (THM-D7 KI-BIC)", gg["satellite_grand"])):
            if not blk:
                continue
            kv, sc, a11 = blk["klein_d3"], blk["screen_d3"], (blk.get("arrow_D11") or {})
            pr2_u = kv.get("time_PR2_compliant_share_unb")
            if pr2_u is None:
                pr2_u = kv.get("time_PR2_compliant_share")
            L.append(f"| **{nm}** | {kv['N']:,} | {kv['norm2_W_pp']:.2e} | {kv['norm2_W_pm']:.2e} | "
                     f"{_fsci(kv.get('norm2_W_mp_unb', kv['norm2_W_mp']))} | {_fsci(kv.get('norm2_W_mm_unb', kv['norm2_W_mm']))} | "
                     f"{fnum(100.0*_fval(pr2_u, 0.0),2)}% | {fnum(100.0*_fval(a11.get('pure_nonenv_arrow_share'), np.nan),1)}% | "
                     f"{fnum(100.0*_fval(a11.get('nonbridge_arrow_share'), np.nan),1)}% | "
                     f"{fnum(100.0*_fval(a11.get('cpg_free_arrow_share'), np.nan),1)}% | {sc['I_nonbridge_MM']:.5f} | "
                     f"{fnum(100.0*_fval(kv.get('p_pal_6bp'), np.nan),3)}% |")
        L.append("")

    L += ["## 1. Whole-Chromosome Klein V_4 & PR2 Decomposition (d = 3 bp, All 3 Phases)", "",
          "| Chr | Arm Class | Length (Mb) | GC % | ‖W_++‖² | ‖W_+-‖² | ‖W_-+‖²_unb | ‖W_--‖²_unb | Whole PR2-ok % | Euch PR2-ok (unb) % | 6-bp Pal Enrich |",
          "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for r in R:
        kv = r["micro_lags"][0]["klein"]
        hn0 = r["micro_lags"][0].get("het_null") or {}
        eu = (r.get("theorem_D7") or {}).get("euchromatin") or {}
        ekv = (eu.get("klein_d3") or {}) if not eu.get("skipped") else {}
        ep = ekv.get("time_PR2_compliant_share_unb")
        if ep is None:
            ep = ekv.get("time_PR2_compliant_share")
        L.append(f"| {r['meta']['label']} | {r['meta'].get('arm_class','')} | {r['length_bp']/1e6:.2f} | "
                 f"{100*r['gc_global']:.2f}% | {kv['norm2_W_pp']:.2e} | {kv['norm2_W_pm']:.2e} | "
                 f"{_fsci(kv.get('norm2_W_mp_unb', kv['norm2_W_mp']))} | {_fsci(kv.get('norm2_W_mm_unb', kv['norm2_W_mm']))} | "
                 f"{fnum(100*(kv['time_PR2_compliant_share'] or 0),2)}% | {fnum(100*ep,2) if ep is not None else 'n/a'}% | "
                 f"{fnum(hn0.get('p_pal_enrichment'),3)}× |")

    L += ["", "## 2. Theorem D7 — Data-Driven Euchromatin vs. Satellite/Tandem Partition (Constrained Kittler–Illingworth BIC on `log10(I_nonbr_raw)` + Entropy Inversion + `H_3 < 3.0` Guard)", "",
          "| Chr | KI-BIC I_nonbr* | ΔBIC | Euch % | Sat % | Sat Blocks | Euch ‖W_-+‖² | Sat ‖W_--‖² | Euch K_nonbr^MM (d=3/90) | Sat 171-bp | Euch unitSVD1 W_+- |",
          "|---|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|"]
    for r in R:
        d7 = r.get("theorem_D7") or {}
        eu, sa, ew = d7.get("euchromatin") or {}, d7.get("satellite_firmware") or {}, d7.get("euch_svd_walk") or {}
        ek, sk_, sfp = eu.get("klein_d3") or {}, sa.get("klein_d3") or {}, sa.get("fft_pitch") or {}
        k3e = (eu.get("screen_kl_d3") or {}).get("K_nonbridge_MM")
        k90e = (eu.get("screen_kl_d90") or {}).get("K_nonbridge_MM")
        L.append(f"| {r['meta']['label']} | {fnum(d7.get('inb_threshold_otsu'),4)} | {fnum(d7.get('bic_2vs1_delta'),1)} | "
                 f"{100*d7.get('euch_fraction',0):.1f}% | {100*d7.get('sat_fraction',0):.1f}% | "
                 f"{len(d7.get('satellite_blocks_mb') or [])} | {_fsci(ek.get('norm2_W_mp'))} | {_fsci(sk_.get('norm2_W_mm'))} | "
                 f"{fnum(k3e,5)} / {fnum(k90e,5)} | {fnum(sfp.get('sat171_peak_excess'),4,True)} | "
                 f"{fnum(100*(ew.get('svd1_W_pm_unit') or 0),1)}% |")

    L += ["", "## 3. Theorems D8 & D9 — Multi-Scale Fractal / Hurst Scaling H (Unified THM-D7 Euchromatin; Delete-1 Macro-Block Jackknife 95% CI; i.i.d. = 0.50)", "",
          "| Chr | H_wls(GC) | H_wls(W_++) | H_wls(W_-+_unb) [JK 95% CI] | Scales | Above WN? | H_wls(NonBr) | H_wls(Omega_cond) | W_++ (Δ=10) | W_-+ (Δ=10) |",
          "|---|---:|---:|---|---:|:---:|---:|---:|---:|---:|"]
    for r in R:
        ch = ((r.get("theorem_D8") or {}).get("euchromatin_only") or {}).get("channels") or {}
        cu = ch.get("W_mp_unb") or {}
        ma = r["disk_scan"]["macro_autocorr"]
        hw_s = (f"{fnum(cu.get('hurst_H_wls'),3)} [{fnum(cu.get('hurst_ci95_lo'),2)},{fnum(cu.get('hurst_ci95_hi'),2)}]"
                if _finite(cu.get("hurst_H_wls")) else "n/a")
        wn_s = "Yes" if cu.get("above_white_noise_95") else ("No" if cu.get("above_white_noise_95") is False else "n/a")
        L.append(f"| {r['meta']['label']} | {fnum((ch.get('gc') or {}).get('hurst_H_wls'),3)} | "
                 f"{fnum((ch.get('W_pp') or {}).get('hurst_H_wls'),3)} | {hw_s} | {cu.get('n_scales_used','n/a')} | {wn_s} | "
                 f"{fnum((ch.get('I_nonbr_MM') or {}).get('hurst_H_wls'),3)} | {fnum((ch.get('Omega_cond_excess') or {}).get('hurst_H_wls'),3)} | "
                 f"{fnum(ma['W_pp']['10'],3,True)} | {fnum(ma['W_mp']['10'],3,True)} |")

    xm = payload.get("cross_chromosome") or {}
    if xm.get("nn_directed_all") and xm.get("nn_directed_euch"):
        L += ["", "## 4. Alignment-Free Directed Nearest Neighbors (Full Disk vs. Euchromatin-Clean W_-+)", "",
              "| Chr | Nearest (Full Disk) | Dist | Nearest (Euchromatin-Clean) | Dist |", "|---|---|---:|---|---:|"]
        for a, b in zip(xm["nn_directed_all"], xm["nn_directed_euch"]):
            L.append(f"| {a['chr']} | {a['nearest']} | {a['dist']:.4f} | {b['nearest']} | {b['dist']:.4f} |")

    epd = gg.get("entropy_phase_diagram_D10") or []
    active_epd = [row for row in epd if not row.get("skipped") and row.get("N", 0) > 0]
    if active_epd:
        L += ["", "## 5. Theorem D10 — 3.12-Gb 8-Strata Entropy Phase Diagram (H_3mer Resolution)", "",
              "| H_3mer Stratum (bits) | Sectors | Neg W_-+_unb % | Pairs N | ‖W_++‖² | ‖W_+-‖² | ‖W_-+‖²_unb | ‖W_--‖²_unb | PR2-ok (unb) % | I_nonbr_MM | Non-Env 3600-D % | CpG-free % |",
              "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
        for row in active_epd:
            kv, sc, a11 = row["klein_d3"], row["screen_d3"], (row.get("arrow_D11") or {})
            pr2_u = kv.get("time_PR2_compliant_share_unb")
            if pr2_u is None:
                pr2_u = kv.get("time_PR2_compliant_share")
            ne_sh = a11.get("pure_nonenv_arrow_share", a11.get("nonbridge_arrow_share"))
            L.append(f"| `{row['label']}` | {row['n_sectors']:,} | {fnum(100.0*_fval(row.get('frac_negative_u_mp'), np.nan),1)}% | "
                     f"{row['N']:,} | {kv['norm2_W_pp']:.2e} | {kv['norm2_W_pm']:.2e} | "
                     f"{_fsci(kv.get('norm2_W_mp_unb'))} | {_fsci(kv.get('norm2_W_mm_unb'))} | "
                     f"{fnum(100.0*_fval(pr2_u, 0.0),2)}% | {sc['I_nonbridge_MM']:.5f} | "
                     f"{fnum(100.0*_fval(ne_sh, np.nan),1)}% | {fnum(100.0*_fval(a11.get('cpg_free_arrow_share'), np.nan),1)}% |")

    L += ["", "## 6. Theorems D11 & D13 — V_4-Closed 3600-D Screen Complement & CpG-Free Anatomy of Human W_-+", "",
          "| Chr | Euch ‖W_-+‖² | Euch Non-Env 3600-D % | Euch NonBr 256-D % | Euch CpG-free % | Euch Top-1 Motif (Share %) | Sat Non-Env 3600-D % | Sat CpG-free % |",
          "|---|---:|---:|---:|---:|---|---:|---:|"]
    for r in R:
        d7 = r.get("theorem_D7") or {}
        ea = (d7.get("euchromatin") or {}).get("arrow_D11") or {}
        sa = (d7.get("satellite_firmware") or {}).get("arrow_D11") or {}
        top1 = (ea.get("top_motifs_all") or [None])[0]
        t1_s = f"`{top1['motif_6mer']}` ({top1['flux_share_pct']:.2f}%)" if top1 else "n/a"
        L.append(f"| {r['meta']['label']} | {_fsci(ea.get('norm2_W_mp'))} | "
                 f"{_fpct(ea.get('pure_nonenv_arrow_share'))} | {_fpct(ea.get('nonbridge_arrow_share'))} | "
                 f"{_fpct(ea.get('cpg_free_arrow_share'))} | {t1_s} | "
                 f"{_fpct(sa.get('pure_nonenv_arrow_share'))} | {_fpct(sa.get('cpg_free_arrow_share'))} |")

    grand_a11 = ((gg.get("euchromatin_grand") or {}).get("arrow_D11")
                 or (gg.get("whole_genome") or {}).get("arrow_D11") or {})
    if grand_a11.get("top_motifs_all"):
        L += ["", "### 6B. Grand-Pooled Euchromatin Top-15 Canonical Directed 6-mer Flux Motifs (Exact Klein Orbit Share)", "",
              "| Rank | Motif 6-mer | Mirror 6-mer | W_-+ Flux | Orbit Share % | Contains CpG? |",
              "|---:|---|---|---:|---:|:---:|"]
        for rank, m in enumerate(grand_a11["top_motifs_all"], 1):
            L.append(f"| {rank} | `{m['motif_6mer']}` | `{m['mirror_6mer']}` | {m['W_mp_flux']:.3e} | "
                     f"{m['flux_share_pct']:.2f}% | {'Yes' if m['has_CpG'] else 'No (CpG-free)'} |")
    if grand_a11.get("top_motifs_cpg_free"):
        L += ["", "### 6C. Grand-Pooled Euchromatin Top-12 Strictly CpG-Free Directed 6-mer Flux Motifs", "",
              "| Rank | Motif 6-mer | Mirror 6-mer | W_-+ Flux | Orbit Share % |",
              "|---:|---|---|---:|---:|"]
        for rank, m in enumerate(grand_a11["top_motifs_cpg_free"], 1):
            L.append(f"| {rank} | `{m['motif_6mer']}` | `{m['mirror_6mer']}` | {m['W_mp_flux']:.3e} | {m['flux_share_pct']:.2f}% |")

    L += ["", "## 7. Theorem D12 — 3.12-Gb Grand Multi-Lag Toric Spectrum (Whole vs. 24-Chr Pure Euchromatin) & WW–SS Nucleosome Pitch"]
    ml_d12 = gg.get("grand_multilag_spectrum_D12") or []
    if ml_d12:
        L += ["", "### 7A. Whole-Genome Grand-Pooled 10-Lag Spectrum (d = 3 to 171 bp)", "",
              "> **Theorem D19 (Exact mod-5 Tandem-Repeat Aliasing of HSat2/3 `(GGAAT)_n`):** For any ideal tandem repeat of period `p`, `c_{i+d} = c_i` iff `d ≡ 0 (mod p)`, confining `B^(d)` to the T-symmetric diagonal (`u == v`) where antisymmetric `W_-+ = W_-- = 0`. On the real whole-genome pool, lags `d ≡ 0 (mod 5)` (`10, 15, 30, 90`) strongly suppress HSat2/3 directed flux, whereas `d ≢ 0 (mod 5)` (`147 ≡ 2`, `171 ≡ 1`) restore it. Table 7B (pure euchromatin) removes this aliasing completely.", "",
              "| Lag d (bp) | d mod 5 | ‖W_++‖² | ‖W_+-‖² | ‖W_-+‖²_unb | ‖W_--‖²_unb | PR2-ok (unb) % | I_bridge | I_nonbr_MM | Non-Env 3600-D % | CpG-free % |",
              "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
        for row in ml_d12:
            kv, sc, a11 = row["klein"], row["screen"], (row.get("arrow_D11") or {})
            pr2_u = kv.get("time_PR2_compliant_share_unb")
            if pr2_u is None:
                pr2_u = kv.get("time_PR2_compliant_share")
            L.append(f"| {row['lag_bp']} | {int(row['lag_bp']) % 5} | {kv['norm2_W_pp']:.2e} | {kv['norm2_W_pm']:.2e} | "
                     f"{_fsci(kv.get('norm2_W_mp_unb'))} | {_fsci(kv.get('norm2_W_mm_unb'))} | "
                     f"{fnum(100.0*_fval(pr2_u, 0.0),2)}% | {sc['I_bridge']:.5f} | {sc['I_nonbridge_MM']:.5f} | "
                     f"{_fpct(a11.get('pure_nonenv_arrow_share'))} | {_fpct(a11.get('cpg_free_arrow_share'))} |")

    eml_d12 = gg.get("grand_euch_multilag_spectrum_D12") or []
    if eml_d12:
        L += ["", "### 7B. Complete 24-Chromosome Pure-Euchromatin Grand 10-Lag Spectrum (THM-D7 Satellite-Stripped)", "",
              "| Lag d (bp) | Pure Pairs N | ‖W_-+‖²_unb | ‖W_--‖²_unb | PR2-ok (unb) % | I_bridge | I_nonbr_MM | Non-Env 3600-D % | CpG-free % | Top-1 Motif (Share %) |",
              "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|"]
        for row in eml_d12:
            kv, sc, a11 = row["klein"], row["screen"], (row.get("arrow_D11") or {})
            pr2_u = kv.get("time_PR2_compliant_share_unb")
            if pr2_u is None:
                pr2_u = kv.get("time_PR2_compliant_share")
            top1 = (a11.get("top_motifs_all") or [None])[0]
            t1 = f"`{top1['motif_6mer']}` ({top1['flux_share_pct']:.2f}%)" if top1 else "n/a"
            L.append(f"| {row['lag_bp']} | {kv['N']:,} | {_fsci(kv.get('norm2_W_mp_unb'))} | "
                     f"{_fsci(kv.get('norm2_W_mm_unb'))} | {fnum(100.0*_fval(pr2_u, 0.0),2)}% | "
                     f"{sc['I_bridge']:.5f} | {sc['I_nonbridge_MM']:.5f} | "
                     f"{_fpct(a11.get('pure_nonenv_arrow_share'))} | {_fpct(a11.get('cpg_free_arrow_share'))} | {t1} |")

    L += ["", "### 7C. Per-Chromosome Euchromatin WW={AA,AT,TA,TT} & SS={CC,CG,GC,GG} Nucleosome Pitch Audit (Hann + Local-Band ENBW-Adj. F_{2,2M_eff})", "",
          "| Chr | WW % | WW Nuc (nt) | WW SNR | SS Nuc (nt) | WW−SS SNR | WW–SS Δφ (°) | WW p_adj (F) | WW−SS p_adj (F) | WW ACF(185/200) |",
          "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|"]
    for r in R:
        ww = ((r.get("theorem_D7") or {}).get("euchromatin") or {}).get("ww_nucleosome_fft") or {}
        L.append(f"| {r['meta']['label']} | {fnum(100.0*_fval(ww.get('ww_global'), np.nan),2)}% | "
                 f"{fnum(ww.get('ww_nuc_band_nt'),1)} | {fnum(ww.get('ww_nuc_band_snr'),2)} | "
                 f"{fnum(ww.get('ss_nuc_band_nt'),1)} | {fnum(ww.get('ww_minus_ss_nuc_snr'),2)} | "
                 f"{fnum(ww.get('ww_ss_phase_diff_deg'),1)} | {_fsci(ww.get('ww_peak_p_adj'))} | "
                 f"{_fsci(ww.get('ww_minus_ss_peak_p_adj'))} | {fnum(ww.get('ww_acf_d185'),4,True)} / {fnum(ww.get('ww_acf_d200'),4,True)} |")

    ec_ml = payload.get("euchromatin_multilag_control_D12") or {}
    if ec_ml.get("chromosomes"):
        L += ["", "### 7D. Optional Per-Chromosome Pure-Euchromatin Multi-Lag Check", "",
              "| Chr | Lag d (bp) | Pure Pairs N | ‖W_-+‖²_unb | ‖W_--‖²_unb | PR2-ok (unb) % | I_nonbr_MM | Non-Env 3600-D % | Top-1 Motif (Share %) |",
              "|---|---:|---:|---:|---:|---:|---:|---:|---|"]
        for clbl, cdat in ec_ml["chromosomes"].items():
            for row in cdat["lags"]:
                t1 = f"`{row['top1_motif']}` ({row['top1_share_pct']:.2f}%)" if row.get("top1_motif") else "n/a"
                L.append(f"| {clbl} | {row['lag_bp']} | {row['N']:,} | {_fsci(row['norm2_W_mp_unb'])} | "
                         f"{_fsci(row['norm2_W_mm_unb'])} | {fnum(100.0*_fval(row['time_PR2_compliant_share_unb'], np.nan),2)}% | "
                         f"{fnum(row['I_nonbridge_MM'],5)} | {_fpct(row['pure_nonenv_arrow_share'])} | {t1} |")

    L += ["", "## 8. Theorems D13–D17 — V_4-Closed Envelope (496⊕3600), Euchromatin 1-Mb HAC Inference, Euchromatin D15 Residual, WW–SS Antiphase & Circular-Shift Lévy Area", "",
          "| Chr | Euch Non-Env 3600-D % | Euch 1-Mb HAC t (sign p unadj.) | Bowker T unadj log10p | Euch D15(M=1k) Res % (cos) | WW–SS Δφ (°) / p_adj (F) | JK WLS H(W_-+_unb) [95% CI] | Lévy p_shift |",
          "|---|---:|---|---:|---|---|---|---:|"]
    for r in R:
        d7 = r.get("theorem_D7") or {}
        ea = (d7.get("euchromatin") or {}).get("arrow_D11") or {}
        ww = (d7.get("euchromatin") or {}).get("ww_nucleosome_fft") or {}
        bk = r.get("theorem_D14_bowker") or {}
        cl = r.get("theorem_D14_cluster_1mb_euch") or (d7.get("euchromatin") or {}).get("cluster_1mb_euch_D14") or r.get("theorem_D14_cluster_1mb") or {}
        d15 = r.get("theorem_D15_comp_current_euch") or (d7.get("euchromatin") or {}).get("comp_current_euch_M1000") or (r.get("theorem_D15_comp_current") or {}).get("M_1000") or {}
        ch_u = (((r.get("theorem_D8") or {}).get("euchromatin_only") or {}).get("channels") or {}).get("W_mp_unb") or {}
        lv = (d7.get("euch_svd_walk") or {}).get("levy_D17_euch") or (r.get("disk_scan") or {}).get("levy_D17") or {}
        cl_s = f"{fnum(cl.get('cluster_t_hac', cl.get('cluster_t')),2)} ({_fsci(cl.get('sign_p_twosided'))})" if cl else "n/a"
        d15_s = (f"{fnum(100.0*_fval(d15.get('W_mp_residual_share'), np.nan),2)}% ({fnum(d15.get('cos_W_mp_obs_vs_cur'),2,True)})"
                 if d15 else "n/a")
        ww_s = f"{fnum(ww.get('ww_ss_phase_diff_deg'),1)}° / {_fsci(ww.get('ww_peak_p_adj'))}" if ww else "n/a"
        hw_s = (f"{fnum(ch_u.get('hurst_H_wls'),3)} [{fnum(ch_u.get('hurst_ci95_lo'),2)},{fnum(ch_u.get('hurst_ci95_hi'),2)}]"
                if _finite(ch_u.get("hurst_H_wls")) else "n/a")
        L.append(f"| {r['meta']['label']} | "
                 f"{fnum(100.0*_fval(ea.get('pure_nonenv_arrow_share'), np.nan),1)}% | {cl_s} | "
                 f"{fnum(bk.get('bowker_T_log10p_unadj', bk.get('bowker_T_log10p')),1)} | {d15_s} | {ww_s} | "
                 f"{hw_s} | {fnum(lv.get('levy_p_shift'),3)} |")

    L += ["", "## 9. Theorem D18 & Cryptographic Provenance — LRV-Calibrated PELT L2 (`MBIC = 3 ln n`), Cross-Channel Boundary Concordance & FASTA SHA-256", "",
          "| Chr | PELT-LRV H_3 Segs (R²) | PELT-LRV logW_-+ Segs (R²) | Concordance (Jaccard / p) | THM-D7 I_nonbr* (clip %) | THM-D7 ΔBIC | FASTA SHA-256 |",
          "|---|---|---|---|---:|---:|---|"]
    for r in R:
        d7 = r.get("theorem_D7") or {}
        p18 = r.get("theorem_D18_pelt") or {}
        ph3 = p18.get("H_3mer") or {}
        pw2 = p18.get("log10_W_mp_unb") or {}
        cc = p18.get("cross_channel_concordance") or {}
        ph_s = f"{ph3['n_segments']} ({ph3['r2_explained']:.3f})" if ph3 else "n/a"
        pw_s = f"{pw2['n_segments']} ({pw2['r2_explained']:.3f})" if pw2 else "n/a"
        cc_s = f"{cc['n_matched']}/{cc['n_a']} ({cc['jaccard']:.2f}, p={_fsci(cc.get('hypergeom_p'))})" if cc else "n/a"
        cl_pct = f"{100.0*_fval(d7.get('frac_clipped'), 0.0):.1f}%"
        sha_s = f"`{r['fasta_sha256'][:16]}`" if r.get("fasta_sha256") else "n/a"
        L.append(f"| {r['meta']['label']} | {ph_s} | {pw_s} | {cc_s} | "
                 f"{fnum(d7.get('inb_threshold_otsu'),4)} ({cl_pct}) | {fnum(d7.get('bic_2vs1_delta'),1)} | {sha_s} |")

    L.append("")
    atomic_write(path, "\n".join(L) + "\n")


def make_human24_figures(payload, out_dir, dpi=300, skip_diskmaps=False):
    if not HAVE_MPL or not payload.get("results"):
        return []
    R = payload["results"]
    fd = Path(out_dir) / "figures"; dm_dir = fd / "diskmaps"; dm_dir.mkdir(parents=True, exist_ok=True)
    saved, lbls, x = [], [r["meta"]["label"] for r in R], np.arange(len(R))

    # FIG 1: 4-panel 24-Chromosome Executive Overview
    fig, axes = plt.subplots(2, 2, figsize=(18, 11))
    fig.suptitle(f"Complete Human T2T-CHM13v2.0 Torus Disk Atlas v{VERSION} ({len(R)} Chromosomes, Zero-RNG)",
                 fontsize=15, fontweight="bold", y=0.98)
    ax = axes[0, 0]
    p_all = [100.0 * _fval(r["micro_lags"][0]["klein"].get("time_PR2_compliant_share_unb", r["micro_lags"][0]["klein"].get("time_PR2_compliant_share")), 0.0) for r in R]
    p_eu = [100.0 * _fval(((r.get("theorem_D7") or {}).get("euchromatin", {}).get("klein_d3") or {}).get("time_PR2_compliant_share_unb", ((r.get("theorem_D7") or {}).get("euchromatin", {}).get("klein_d3") or {}).get("time_PR2_compliant_share")), np.nan) for r in R]
    ax.bar(x - 0.19, p_all, 0.38, color="#ff7f0e", alpha=0.85, label="Whole Chr (incl. Satellites)")
    ax.bar(x + 0.19, p_eu, 0.38, color="#1f77b4", label="Euchromatin Only (THM-D7 KI-BIC)")
    ax.axhline(99.0, color="green", ls="--", lw=0.9, label="99% PR2-Compliant Flux Regime")
    ax.set_ylim(50, 101.5); ax.set_xticks(x); ax.set_xticklabels(lbls, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("PR2-Compliant Directed Flux (%)"); ax.set_title("A. PR2-Compliant Directed Flux Share: Whole Chr vs. Euchromatin")
    ax.legend(fontsize=8, loc="lower left", framealpha=0.88); ax.grid(True, axis="y", ls=":", alpha=0.35)

    ax = axes[0, 1]
    k_e = [_fval(((r.get("theorem_D7") or {}).get("euchromatin", {}).get("screen_kl_d3") or {}).get("K_nonbridge_MM"), np.nan) for r in R]
    k_s = [np.nan if ((r.get("theorem_D7") or {}).get("satellite_firmware") or {}).get("skipped", True)
           else _fval(((r.get("theorem_D7") or {}).get("satellite_firmware", {}).get("screen_kl_d3") or {}).get("K_nonbridge_MM"), np.nan)
           for r in R]
    ax.bar(x - 0.19, k_e, 0.38, color="#17becf", label="Euchromatin K_nonbr^MM (d=3)")
    ax.bar(x + 0.19, k_s, 0.38, color="#d62728", label="Satellite/Tandem K_nonbr^MM (d=3)")
    ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels(lbls, rotation=45, ha="right", fontsize=8)
    valid_pos = [v for v in (k_e + k_s) if _finite(v) and v > 0]
    y_bot = min(valid_pos) * 0.55 if valid_pos else 1e-3
    y_top = max(valid_pos) * 4.5 if valid_pos else 5.0
    ax.set_ylim(y_bot, y_top)
    for xi_, ks_ in zip(x, k_s):
        if not _finite(ks_):
            ax.text(xi_ + 0.19, y_bot * 1.18, "n/a", ha="center", va="bottom", fontsize=6.5, color="#777777", rotation=90)
    ax.set_ylabel("nats / triplet (log)"); ax.set_title("B. Non-Bridge Screen Memory: Euchromatin vs. Satellite/Tandem Compartment")
    ax.legend(fontsize=8, loc="upper center", ncol=2, framealpha=0.92); ax.grid(True, axis="y", ls=":", alpha=0.35)

    ax = axes[1, 0]
    sat_pct = [100.0 * ((r.get("theorem_D7") or {}).get("sat_fraction") or 0.0) for r in R]
    sat171 = [(((r.get("theorem_D7") or {}).get("satellite_firmware", {}).get("fft_pitch") or {}).get("sat171_peak_excess") or 0.0) for r in R]
    ax.bar(x, sat_pct, 0.65, color="#c44e52", alpha=0.78)
    ax.set_xticks(x); ax.set_xticklabels(lbls, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Satellite Disk Share (%)", color="#c44e52")
    ax2 = ax.twinx(); ax2.plot(x, sat171, color="#1a1a1a", marker="o", ms=4.5, lw=1.4)
    ax2.set_ylabel("171-bp Peak Excess", color="#1a1a1a")
    ax.set_title("C. Satellite Firmware Share & 171-bp Alpha-Satellite Resonance"); ax.grid(True, axis="y", ls=":", alpha=0.35)

    ax = axes[1, 1]
    def _hv(r, ch, prefer_2d=False):
        d8_ = r.get("theorem_D8") or {}
        blk_ = d8_.get("euchromatin_only") or {}
        c_ = (blk_.get("channels") or {}).get(ch) or {}
        v = c_.get("hurst_H_wls")
        if v is None:
            v = c_.get("hurst_H")
        return _fval(v, np.nan)
    h_all_vals = []
    for ch, mk, col, nm in (("gc", "s", "#555555", "H_wls(GC isochores)"), ("W_pp", "^", "#1f77b4", "H_wls(W_++ static)"),
                            ("W_mp_unb", "o", "#d62728", "H_wls(W_-+_unb arrow)"), ("I_nonbr_MM", "d", "#17becf", "H_wls(I_nonbr^MM)")):
        ys_ = [_hv(r, ch) for r in R]
        h_all_vals.extend([v for v in ys_ if _finite(v)])
        ax.plot(x, ys_, marker=mk, ms=4.5, lw=1.6, color=col, label=nm)
    ax.axhline(0.50, color="k", ls="--", lw=1.0, label="White Noise (H=0.50)")
    ylo_d = min(0.42, min(h_all_vals) - 0.04) if h_all_vals else 0.42
    yhi_d = max(1.05, max(h_all_vals) + 0.04) if h_all_vals else 1.05
    ax.set_ylim(ylo_d, yhi_d); ax.set_xticks(x); ax.set_xticklabels(lbls, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("WLS Hurst Exponent H (100 kb -> 6.4 Mb)"); ax.set_title("D. Theorem D8: Unified Euchromatin Megabase Persistence (WLS)")
    ax.legend(fontsize=8, loc="lower right", ncol=2, framealpha=0.88); ax.grid(True, ls=":", alpha=0.35)
    fig.tight_layout(rect=[0, 0, 1, 0.955])
    p1 = fd / "fig1_human24_overview.png"; fig.savefig(p1, dpi=dpi, bbox_inches="tight"); plt.close(fig); saved.append(str(p1))
    print(f"  [FIG] {p1}")

    # FIG 2: Complete 24-Chromosome Unbiased Karyotype Ideogram (SNR>=3 Noise-Floor Masked)
    max_mb = max(r["length_bp"] / 1e6 for r in R)
    fig, ax = plt.subplots(figsize=(17, max(8.5, 0.48 * len(R) + 2.2)))
    fig.suptitle(f"Human T2T-CHM13v2.0 Unbiased Toric Ideogram — All {len(R)} Chromosomes (Zero-RNG, SNR>=3 Masked)",
                 fontsize=14, fontweight="bold", y=0.985)
    all_surv = []
    for r_ in R:
        tr_ = r_["disk_scan"]["tracks"]
        wmp_u_ = np.asarray(tr_.get("W_mp_unb", tr_["W_mp"]), dtype=np.float64)
        snr_s_ = np.asarray(tr_.get("W_mp_snr_sec", np.full_like(wmp_u_, 10.0)), dtype=np.float64)
        ok_s_ = np.isfinite(wmp_u_) & (wmp_u_ > 0.0) & np.isfinite(snr_s_) & (snr_s_ >= 3.0)
        if np.any(ok_s_):
            all_surv.append(np.log10(wmp_u_[ok_s_]))
    if all_surv:
        cat_s = np.concatenate(all_surv)
        vmin, vmax = float(np.percentile(cat_s, 2)), float(np.percentile(cat_s, 98))
        if vmax - vmin < 0.25:
            vmin, vmax = vmin - 0.25, vmax + 0.25
    else:
        vmin, vmax = -5.4, -2.6
    cmap = _get_cmap("viridis").copy()
    cmap.set_bad("#dcdcdc")
    for idx, r in enumerate(R):
        y0 = len(R) - 1 - idx
        tr = r["disk_scan"]["tracks"]
        pos = np.asarray(tr["pos_mb"], dtype=np.float64)
        wmp_u = np.asarray(tr.get("W_mp_unb", tr["W_mp"]), dtype=np.float64)
        snr_s = np.asarray(tr.get("W_mp_snr_sec", np.full_like(wmp_u, 10.0)), dtype=np.float64)
        val = np.where((wmp_u > 0.0) & (snr_s >= 3.0), np.log10(np.maximum(wmp_u, 1e-8)), np.nan)
        if pos.size >= 2:
            x_edges = np.concatenate([[0.0], pos + 0.5 * float(pos[1] - pos[0])])
            ax.pcolormesh(x_edges, [y0 - 0.34, y0 + 0.34], np.ma.masked_invalid(val[None, :]),
                          cmap=cmap, vmin=vmin, vmax=vmax, shading="flat")
            for b0, b1 in ((r.get("theorem_D7") or {}).get("satellite_blocks_mb") or []):
                ax.plot([b0, b1], [y0 + 0.37, y0 + 0.37], color="#ff0055", lw=2.4, solid_capstyle="butt")
                ax.plot([b0, b1], [y0 - 0.37, y0 - 0.37], color="#ff0055", lw=2.4, solid_capstyle="butt")
            L_mb = r["length_bp"] / 1e6
            ax.plot([0, L_mb, L_mb, 0, 0], [y0 - 0.34, y0 - 0.34, y0 + 0.34, y0 + 0.34, y0 - 0.34], color="#222222", lw=0.6)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax)); sm.set_array(np.array([]))
    fig.colorbar(sm, ax=ax, pad=0.015, aspect=35).set_label("log10( ‖W_-+‖²_unb ) [Unbiased PR2-Compliant Flux; Gray = SNR < 3]", fontsize=10)
    ax.set_yticks(range(len(R))); ax.set_yticklabels(lbls[::-1], fontsize=9, fontweight="bold")
    ax.set_xlim(-0.01 * max_mb, max_mb * 1.02); ax.set_ylim(-0.8, len(R) - 0.2)
    ax.set_xlabel("Physical Coordinate (Mb)  |  Magenta brackets = Data-Driven THM-D7 Satellite/Centromere/rDNA Blocks (Constrained Kittler–Illingworth BIC on log10(I_nonbr_raw))", fontsize=10)
    ax.grid(True, axis="x", ls=":", alpha=0.35)
    fig.tight_layout(rect=[0, 0, 1, 0.965])
    p2 = fd / "fig2_human24_karyotype_ideogram.png"; fig.savefig(p2, dpi=dpi, bbox_inches="tight"); plt.close(fig); saved.append(str(p2))
    print(f"  [FIG] {p2}")

    # FIG 3: 24x24 Alignment-Free Similarity Heatmaps
    xm = payload.get("cross_chromosome")
    if xm and len(xm.get("labels", [])) >= 2:
        fig, axes = plt.subplots(1, 2, figsize=(16.5, 7.2))
        fig.suptitle("Alignment-Free 24×24 Directed Toric Distance (Cosine Distance in W_-+ ⊕ W_--)", fontsize=13.5, fontweight="bold", y=0.98)
        xl = xm["labels"]
        for ax_i, key_m, ttl in ((axes[0], "dist_directed_all", "A. Full Disk (Satellites Included)"),
                                 (axes[1], "dist_directed_euch", "B. Satellite-Stripped Euchromatin Only (THM-D7 KI-BIC)")):
            M = np.asarray(xm[key_m], dtype=np.float64)
            im = ax_i.imshow(M, cmap="magma_r", vmin=0.0, vmax=max(0.15, float(np.percentile(M, 98))))
            ax_i.set_xticks(range(len(xl))); ax_i.set_xticklabels(xl, rotation=90, fontsize=7.5)
            ax_i.set_yticks(range(len(xl))); ax_i.set_yticklabels(xl, fontsize=7.5)
            ax_i.set_title(ttl, fontsize=11); fig.colorbar(im, ax=ax_i, fraction=0.046, pad=0.04)
        fig.tight_layout(rect=[0, 0, 1, 0.94])
        p3 = fd / "fig3_human24_similarity_heatmaps.png"; fig.savefig(p3, dpi=dpi, bbox_inches="tight"); plt.close(fig); saved.append(str(p3))
        print(f"  [FIG] {p3}")

    # FIG 4: Theorem D10 — 3.12-Gb 8-Strata Entropy Phase Diagram
    gg = payload.get("grand_genome") or {}
    epd = [row for row in (gg.get("entropy_phase_diagram_D10") or [])
           if not row.get("skipped") and row.get("N", 0) > 0]
    if len(epd) >= 2:
        fig, axes = plt.subplots(1, 2, figsize=(17, 6.2))
        fig.suptitle(f"Theorem D10: 8-Strata Entropy Phase Diagram Across {len(R)} Human T2T Chromosomes (Zero-RNG)",
                     fontsize=14, fontweight="bold", y=0.98)
        ex = np.arange(len(epd))
        elbls = [f"{row['label']}\n(n={row['n_sectors']:,})" for row in epd]
        wmp_u = [max(float((row["klein_d3"].get("norm2_W_mp_unb") or 0.0)), 1e-11) for row in epd]
        wmm_u = [max(float((row["klein_d3"].get("norm2_W_mm_unb") or 0.0)), 1e-11) for row in epd]
        pr2_u = [100.0 * float((row["klein_d3"].get("time_PR2_compliant_share_unb")
                                if row["klein_d3"].get("time_PR2_compliant_share_unb") is not None
                                else (row["klein_d3"].get("time_PR2_compliant_share") or 0.0))) for row in epd]
        inb = [max(float(row["screen_d3"]["I_nonbridge_MM"]), 1e-5) for row in epd]
        nb_sh = [100.0 * _fval(((row.get("arrow_D11") or {}).get("pure_nonenv_arrow_share", (row.get("arrow_D11") or {}).get("nonbridge_arrow_share"))), 0.0) for row in epd]
        cg_sh = [100.0 * _fval(((row.get("arrow_D11") or {}).get("cpg_free_arrow_share")), 0.0) for row in epd]

        split_x = None
        for i_e, row in enumerate(epd):
            if row["label"].startswith("5.65"):
                split_x = i_e - 0.5
                break

        ax = axes[0]
        ax.plot(ex, wmp_u, marker="o", lw=2.0, color="#1f77b4", label="‖W_-+‖²_unb (PR2-ok arrow)")
        ax.plot(ex, wmm_u, marker="s", lw=2.0, color="#ff7f0e", label="‖W_--‖²_unb (chiral flow)")
        ax.set_yscale("log"); ax.set_xticks(ex); ax.set_xticklabels(elbls, fontsize=8)
        ax.set_ylabel("Unbiased Klein Energy (F², log)"); ax.set_xlabel("Sector 3-mer Entropy H_3mer (bits)")
        if split_x is not None:
            ax.axvline(split_x, color="crimson", ls="--", lw=1.4, label="H = 5.65 bits reference")
        ax2 = ax.twinx()
        ax2.plot(ex, pr2_u, marker="^", lw=2.0, ls="-.", color="#2ca02c", label="PR2-ok Share (unb %)")
        ax2.set_ylim(0, 105); ax2.set_ylabel("PR2-Compliant Time-Arrow (%)", color="#2ca02c")
        ax.set_title("A. Klein V_4 Energies & PR2-Compliant Share Across H_3mer Strata")
        h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="lower left", framealpha=0.9); ax.grid(True, ls=":", alpha=0.35)

        ax = axes[1]
        ax.bar(ex - 0.18, nb_sh, 0.36, color="#17becf", label="THM-D13 Non-Env 3600-D (%)")
        ax.bar(ex + 0.18, cg_sh, 0.36, color="#9467bd", label="CpG-Free Arrow (%)")
        ax.set_ylim(0, 124); ax.set_xticks(ex); ax.set_xticklabels(elbls, fontsize=8)
        ax.set_ylabel("Share of ‖W_-+‖² (%)"); ax.set_xlabel("Sector 3-mer Entropy H_3mer (bits)")
        if split_x is not None:
            ax.axvline(split_x, color="crimson", ls="--", lw=1.4)
        ax3 = ax.twinx()
        ax3.plot(ex, inb, marker="D", lw=2.0, color="#d62728", label="I_nonbr^MM (nats, log)")
        ax3.set_yscale("log"); ax3.set_ylim(min(inb) * 0.6, max(inb) * 4.0)
        ax3.set_ylabel("I_nonbridge^MM (nats/triplet)", color="#d62728")
        ax.set_title("B. Non-Bridge Screen Memory & CpG-Free Core Across Entropy Strata")
        h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax3.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper center", ncol=3, framealpha=0.92); ax.grid(True, axis="y", ls=":", alpha=0.35)

        fig.tight_layout(rect=[0, 0, 1, 0.94])
        p4 = fd / "fig4_human24_D10_entropy_phase_diagram.png"
        fig.savefig(p4, dpi=dpi, bbox_inches="tight"); plt.close(fig); saved.append(str(p4))
        print(f"  [FIG] {p4}")

    # FIG 5: Theorem D11 — Orthogonal Non-Bridge & CpG-Free Arrow Anatomy + Top-15 6-mers
    grand_a11 = ((gg.get("euchromatin_grand") or {}).get("arrow_D11")
                 or (gg.get("whole_genome") or {}).get("arrow_D11") or {})
    if grand_a11.get("top_motifs_all"):
        fig, axes = plt.subplots(1, 2, figsize=(18, 6.8), gridspec_kw={"width_ratios": [1.15, 1.0]})
        fig.suptitle(f"Theorem D11: Exact Orthogonal Non-Bridge & CpG-Free Anatomy of Human W_-+ ({len(R)} Chr, Zero-RNG)",
                     fontsize=14, fontweight="bold", y=0.98)
        ax = axes[0]
        nb_chr = [100.0 * _fval((((r.get("theorem_D7") or {}).get("euchromatin") or {}).get("arrow_D11") or {}).get("pure_nonenv_arrow_share", (((r.get("theorem_D7") or {}).get("euchromatin") or {}).get("arrow_D11") or {}).get("nonbridge_arrow_share")), np.nan) for r in R]
        cg_chr = [100.0 * _fval((((r.get("theorem_D7") or {}).get("euchromatin") or {}).get("arrow_D11") or {}).get("cpg_free_arrow_share"), np.nan) for r in R]
        ax.bar(x - 0.19, nb_chr, 0.38, color="#17becf", label="THM-D13 Pure Non-Env 3600-D Share (‖W_-+^3600‖² / ‖W_-+‖²)")
        ax.bar(x + 0.19, cg_chr, 0.38, color="#9467bd", label="CpG-Free Arrow Share (strictly no CpG)")
        ax.set_ylim(0, 118); ax.set_xticks(x); ax.set_xticklabels(lbls, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("Share of Euchromatin ‖W_-+‖² (%)")
        ax.set_title("A. Non-Envelope 3600-D & CpG-Free Share of Euchromatin W_-+")
        ax.legend(fontsize=8.5, loc="upper center", ncol=2, framealpha=0.92); ax.grid(True, axis="y", ls=":", alpha=0.35)

        ax = axes[1]
        motifs = grand_a11["top_motifs_all"][:15][::-1]
        y_m = np.arange(len(motifs))
        shares = [float(m["flux_share_pct"]) for m in motifs]
        cols = ["#d62728" if m["has_CpG"] else "#1f77b4" for m in motifs]
        mlbls = [f"{m['motif_6mer']} / {m['mirror_6mer']}" if m["motif_6mer"] != m["mirror_6mer"] else m["motif_6mer"] for m in motifs]
        bars = ax.barh(y_m, shares, color=cols, alpha=0.88)
        ax.set_yticks(y_m); ax.set_yticklabels(mlbls, fontsize=8, family="monospace")
        ax.set_xlabel("Exact Klein-Orbit Share of Grand Euchromatin ‖W_-+‖² (%)")
        ax.set_title("B. Top-15 Directed 6-mer Flux Orbits (Red = CpG, Blue = CpG-Free)")
        xmax = max(shares) * 1.18 if shares else 1.0
        ax.set_xlim(0, xmax)
        for b_, sh_ in zip(bars, shares):
            ax.text(sh_ + xmax * 0.012, b_.get_y() + b_.get_height() / 2.0, f"{sh_:.2f}%", va="center", fontsize=7.5)
        ax.grid(True, axis="x", ls=":", alpha=0.35)

        fig.tight_layout(rect=[0, 0, 1, 0.94])
        p5 = fd / "fig5_human24_D11_arrow_anatomy_motifs.png"
        fig.savefig(p5, dpi=dpi, bbox_inches="tight"); plt.close(fig); saved.append(str(p5))
        print(f"  [FIG] {p5}")

    # FIG 6: Theorems D8, D9, D12 & D16 — 24-Chr Multi-Lag (Whole vs Pure Euchromatin), WW-SS & Jackknife Hurst
    ml_d12 = gg.get("grand_multilag_spectrum_D12") or []
    eml_d12 = gg.get("grand_euch_multilag_spectrum_D12") or []
    fig, axes = plt.subplots(1, 3, figsize=(19.5, 6.0))
    fig.suptitle(f"Theorems D8, D9, D12 & D16: Multi-Lag Decay (Whole vs. Pure Euchromatin), WW−SS Antiphase & Jackknife Hurst ({len(R)} Chr)",
                 fontsize=13.5, fontweight="bold", y=0.985)

    ax = axes[0]
    if ml_d12:
        lx = np.arange(len(ml_d12))
        llbls = [str(row["lag_bp"]) for row in ml_d12]
        wmp_l = [max(float((row["klein"].get("norm2_W_mp_unb") or 0.0)), 1e-10) for row in ml_d12]
        inb_l = [max(float(row["screen"]["I_nonbridge_MM"]), 1e-6) for row in ml_d12]
        ax.plot(lx, wmp_l, marker="o", lw=1.6, ls="--", color="#ff7f0e", alpha=0.85, label="Whole ‖W_-+‖²_unb (mod-5 HSat)")
        if eml_d12 and len(eml_d12) == len(ml_d12):
            wmp_eu_l = [max(float((row["klein"].get("norm2_W_mp_unb") or 0.0)), 1e-10) for row in eml_d12]
            nb_l = [100.0 * _fval(((row.get("arrow_D11") or {}).get("pure_nonenv_arrow_share")), 0.0) for row in eml_d12]
            ax.plot(lx, wmp_eu_l, marker="o", lw=2.2, color="#d62728", label="24-Chr Pure Euch ‖W_-+‖²_unb")
        else:
            nb_l = [100.0 * _fval(((row.get("arrow_D11") or {}).get("pure_nonenv_arrow_share")), 0.0) for row in ml_d12]
        ax.plot(lx, inb_l, marker="s", lw=1.8, color="#17becf", label="Whole I_nonbr^MM (left, log)")
        ax.set_yscale("log"); ax.set_xticks(lx); ax.set_xticklabels(llbls, fontsize=8.5)
        ax.set_xlabel("Lag d (bp)  [d ≡ 0 mod 5 nulls HSat2/3 (GGAAT)_n]"); ax.set_ylabel("Flux / Mutual Info (log)")
        for special_lag, scol in ((147, "#2ca02c"), (171, "#ff7f0e")):
            for ii, row in enumerate(ml_d12):
                if row["lag_bp"] == special_lag:
                    ax.axvline(ii, color=scol, ls="--", lw=1.0, alpha=0.7)
        ax2 = ax.twinx()
        ax2.plot(lx, nb_l, marker="^", lw=1.5, ls=":", color="#1f77b4", label="Non-Env 3600-D % (right)")
        ax2.set_ylim(0, 110); ax2.set_ylabel("Non-Env 3600-D Arrow Share (%)", color="#1f77b4")
        h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, fontsize=6.8, loc="lower left", framealpha=0.92)
    else:
        ax.text(0.5, 0.5, "No multi-lag cache available", ha="center", va="center", transform=ax.transAxes)
    ax.set_title("A. 3.12-Gb Multi-Lag Spectrum (Whole vs. Pure Euchromatin)"); ax.grid(True, ls=":", alpha=0.35)

    ax = axes[1]
    ctr_snr = [_fval((((r.get("theorem_D7") or {}).get("euchromatin") or {}).get("ww_nucleosome_fft") or {}).get("ww_minus_ss_nuc_snr",
               (((r.get("theorem_D7") or {}).get("euchromatin") or {}).get("ww_nucleosome_fft") or {}).get("ww_nuc_band_snr")), np.nan) for r in R]
    ww_snr = [_fval((((r.get("theorem_D7") or {}).get("euchromatin") or {}).get("ww_nucleosome_fft") or {}).get("ww_nuc_band_snr"), np.nan) for r in R]
    ph_deg = [_fval((((r.get("theorem_D7") or {}).get("euchromatin") or {}).get("ww_nucleosome_fft") or {}).get("ww_ss_phase_diff_deg"), np.nan) for r in R]
    ax.bar(x - 0.18, np.nan_to_num(ctr_snr, nan=0.0), 0.36, color="#2ca02c", alpha=0.85, label="WW−SS Contrast SNR")
    ax.bar(x + 0.18, np.nan_to_num(ww_snr, nan=0.0), 0.36, color="#98df8a", alpha=0.85, label="WW-Only SNR")
    max_snr = max([v for v in (ctr_snr + ww_snr) if _finite(v)], default=5.0)
    ax.set_ylim(0.0, max_snr * 1.35)
    ax.set_xticks(x); ax.set_xticklabels(lbls, rotation=45, ha="right", fontsize=7.5)
    ax.set_ylabel("Nucleosome-Band SNR", color="#2ca02c")
    ax2 = ax.twinx()
    ax2.plot(x, ph_deg, marker="o", ms=4.5, lw=1.5, color="#1a1a1a", label="WW–SS Δφ (°)")
    ax2.axhspan(150.0, 180.0, color="#ff7f0e", alpha=0.14, label="Antiphase Band (150–180°)")
    ax2.axhline(180.0, color="#d62728", ls="--", lw=0.9, alpha=0.75)
    ax2.set_ylim(0.0, 238.0); ax2.set_ylabel("WW–SS Relative Phase Δφ (°)", color="#1a1a1a")
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=6.6, loc="upper center", ncol=2, framealpha=0.92)
    ax.set_title("B. Euchromatin WW−SS Antiphase Contrast & Phase Lock"); ax.grid(True, axis="y", ls=":", alpha=0.35)

    ax = axes[2]
    def _h_eu(r_, ch_, field_="hurst_H_wls"):
        blk_ = (r_.get("theorem_D8") or {}).get("euchromatin_only") or {}
        c_ = (blk_.get("channels") or {}).get(ch_) or {}
        v_ = c_.get(field_)
        if v_ is None and field_ == "hurst_H_wls":
            v_ = c_.get("hurst_H")
        return _fval(v_, np.nan)
    h_gc = [_h_eu(r, "gc", "hurst_H_wls") for r in R]
    h_u = [_h_eu(r, "W_mp_unb", "hurst_H_wls") for r in R]
    h_lo = [_h_eu(r, "W_mp_unb", "hurst_ci95_lo") for r in R]
    h_hi = [_h_eu(r, "W_mp_unb", "hurst_ci95_hi") for r in R]
    yerr_lo = [max(0.0, hv_ - lo_) if (_finite(hv_) and _finite(lo_)) else 0.0 for hv_, lo_ in zip(h_u, h_lo)]
    yerr_hi = [max(0.0, hi_ - hv_) if (_finite(hv_) and _finite(hi_)) else 0.0 for hv_, hi_ in zip(h_u, h_hi)]
    ax.plot(x, h_gc, marker="s", ms=4.0, lw=1.3, color="#666666", alpha=0.85, label="H_wls(GC, THM-D7 Euch)")
    ax.errorbar(x, h_u, yerr=[yerr_lo, yerr_hi], fmt="-o", ms=4.2, lw=1.6, capsize=2.2, color="#1f77b4", label="H_wls(W_-+_unb) ±95% JK CI")
    ax.axhline(0.50, color="k", ls=":", lw=1.0, label="White Noise (H=0.50)")
    ylo, yhi = 0.22, 1.22
    ax.set_ylim(ylo, yhi); ax.set_xticks(x); ax.set_xticklabels(lbls, rotation=45, ha="right", fontsize=7.5)
    for xi_, hv_, lo_, hi_, lb_ in zip(x, h_u, h_lo, h_hi, lbls):
        if _finite(hv_) and hv_ < ylo:
            ax.plot(xi_, ylo + 0.025, marker="v", ms=6.5, color="#d62728", zorder=6)
            ax.annotate(f"{lb_}: H={hv_:.2f} [{lo_:.2f},{hi_:.2f}]", xy=(xi_, ylo + 0.025),
                        xytext=(max(0.5, xi_ - 10.5), ylo + 0.055), fontsize=6.8, color="#b2182b",
                        bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="#d62728", alpha=0.92),
                        arrowprops=dict(arrowstyle="->", color="#d62728", lw=0.9))
    ax.set_ylabel("WLS Hurst Exponent H"); ax.set_title("C. THM-D8/D9: Unified Euchromatin WLS Hurst (±95% Jackknife CI)")
    ax.legend(fontsize=7.0, loc="lower left", bbox_to_anchor=(0.05, 0.02), framealpha=0.92); ax.grid(True, ls=":", alpha=0.35)

    fig.tight_layout(rect=[0, 0, 1, 0.935])
    p6 = fd / "fig6_human24_D12_multilag_ww_hurst.png"
    fig.savefig(p6, dpi=dpi, bbox_inches="tight"); plt.close(fig); saved.append(str(p6))
    print(f"  [FIG] {p6}")

    # FIG 7+: All individual 5-track hard-drive maps
    if skip_diskmaps:
        return saved
    for r in R:
        k, tr = r["key"], r["disk_scan"]["tracks"]
        pos = np.asarray(tr["pos_mb"], dtype=np.float64)
        d7 = r.get("theorem_D7") or {}
        has_ew = bool(d7.get("euch_svd_walk") and d7.get("sat_fraction", 0.0) > 0.005)
        fig, axes = plt.subplots(5 if has_ew else 4, 1, figsize=(15, 11.2 if has_ew else 9.2), sharex=True)
        fig.suptitle(f"Human T2T-CHM13v2.0 Disk Map: {r['meta']['label']} ({r['meta']['acc']}, "
                     f"{r['length_bp']/1e6:.2f} Mb, {r['sector_bp']//1000}-kb sectors)",
                     fontsize=13, fontweight="bold", y=0.985)
        def _shade(a_):
            for b0, b1 in (d7.get("satellite_blocks_mb") or []): a_.axvspan(b0, b1, color="crimson", alpha=0.10, lw=0)
        def _nw(v):
            a = np.asarray(v, dtype=np.float64); return a / max(float(np.nanmax(np.abs(a))), 1e-12)
        ax = axes[0]; _shade(ax)
        ax.plot(pos, _nw(tr["walk_gc_skew"]), color="#7f7f7f", lw=1.3, ls="--", label="GC-skew walk")
        ax.plot(pos, _nw(tr["walk_W_pm"]), color="#1f77b4", lw=1.8, label="W_+- Toric walk (1024-D)")
        ax.plot(pos, _nw(tr["walk_W_mm"]), color="#ff7f0e", lw=1.6, label="W_-- Toric walk (992-D)")
        ax.set_ylabel("Norm. Walk"); ax.legend(fontsize=8, loc="upper right", ncol=3, framealpha=0.88); ax.grid(True, ls=":", alpha=0.35)
        ax.set_title("1. Full-Disk Strand Polarity (pink = THM-D7 Satellite/Tandem)", fontsize=10)
        ro = 0
        if has_ew:
            ro = 1; ew = d7["euch_svd_walk"]; axe = axes[1]; _shade(axe)
            axe.plot(pos, _nw(ew["walk_gc_euch"]), color="#7f7f7f", lw=1.3, ls="--", label="Euch GC-skew")
            axe.plot(pos, _nw(ew["walk_W_pm_euch"]), color="#1f77b4", lw=1.8, label="Euch W_+- unit-SVD")
            axe.plot(pos, _nw(ew["walk_W_mm_euch"]), color="#ff7f0e", lw=1.6, label="Euch W_-- unit-SVD")
            axe.set_ylabel("Euch Walk"); axe.legend(fontsize=8, loc="upper right", ncol=3, framealpha=0.88); axe.grid(True, ls=":", alpha=0.35)
            axe.set_title(f"1B. Satellite-Stripped Euchromatin-Only Toric Walk (THM-D7 KI-BIC; r_GCsk={fnum(ew['corr_W_pm_vs_gc_skew'],3,True)})", fontsize=10)
        at = axes[1 + ro:]
        ax = at[0]; _shade(ax)
        ax.plot(pos, tr["W_mp"], color="#d62728", lw=1.1, label="‖W_-+‖² (PR2-ok)"); ax.plot(pos, tr["W_mm"], color="#ff7f0e", lw=0.95, alpha=0.85, label="‖W_--‖² (chiral)")
        ax.set_yscale("log"); ax.set_ylabel("Flux (F²)"); ax.legend(fontsize=8, loc="upper right", framealpha=0.88); ax.grid(True, ls=":", alpha=0.35)
        ax2 = ax.twinx(); ax2.plot(pos, tr["H_3mer_bits"], color="#2ca02c", lw=0.85, alpha=0.65)
        ax2.set_ylabel("H_3mer (bits)", color="#2ca02c")
        ax.set_title("2. Directed Within-Sector Syntax (W_-+, W_--) & Entropy Drops", fontsize=10)
        ax = at[1]; _shade(ax)
        ax.plot(pos, tr["Omega_cond_excess"], color="#9467bd", lw=1.1, label="Ω_cond − bias (3-bp phase lock)")
        ax.axhline(0.0, color="k", ls="--", lw=0.8); ax.set_ylabel("nats"); ax.legend(fontsize=8, loc="upper right", framealpha=0.88); ax.grid(True, ls=":", alpha=0.35)
        ax.set_title("3. Blind 3-bp Frame-Phase Lock (Coding / 3n-Satellite Detector)", fontsize=10)
        ax = at[2]; _shade(ax)
        use_raw_dm = ("D_5lead" in tr) and ("raw" in str(d7.get("partition_method", "")))
        inb_dm = (np.asarray(tr["D_5lead"], dtype=np.float64) + np.asarray(tr["D_3wob"], dtype=np.float64)
                  + np.asarray(tr["D_1_3"], dtype=np.float64)) if use_raw_dm else tr["I_nonbr_MM"]
        ax.plot(pos, inb_dm, color="#17becf", lw=1.1,
                label="I_nonbridge_raw (screen memory)" if use_raw_dm else "I_nonbridge^MM (non-bridge screen memory)")
        if _finite(d7.get("inb_threshold_otsu")):
            thr_dm = float(d7["inb_threshold_otsu"])
            ax.axhline(thr_dm, color="crimson", ls=":", lw=1.1, alpha=0.85, label=f"KI-BIC thr={thr_dm:.3f}")
        ax.set_ylabel("nats/triplet", color="#17becf"); ax.legend(fontsize=8, loc="upper right", framealpha=0.88); ax.grid(True, ls=":", alpha=0.35)
        ax3 = ax.twinx(); ax3.plot(pos, tr["gc"], color="#555555", lw=0.85, alpha=0.6); ax3.set_ylabel("GC fraction", color="#555555")
        ax.set_xlabel("Chromosome coordinate (Mb)"); ax.set_title("4. Non-Bridge 16×16 Screen Memory vs GC Isochores", fontsize=10)
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        pk = dm_dir / f"fig_diskmap_{k}.png"; fig.savefig(pk, dpi=dpi, bbox_inches="tight"); plt.close(fig); saved.append(str(pk))
    print(f"  [FIG] Saved {len(R)} individual 5-track disk maps in {dm_dir}")
    return saved




# =====================================================================
# 10. PAPER 5 (v2.0.0): ZERO-RNG VIROME & ONCOGENOMICS TORIC ENGINE
#     [THM-V1] Viral Toric Impedance & 24-Chr Alignment Spectrum
#     [THM-V2] Deterministic Proviral Insertion & Chiral Flare Operator
#     [THM-V3] Exact Analytic Mutagenesis Jacobian J_S (7 COSMIC Processes)
#     [THM-V4] Fine-Scale Human Cancer Driver & Fragile-Site Toric Audit
# =====================================================================
_VIRAL_33_SPEC = [
    # --- Integrating Retroviruses (ssRNA-RT -> dsDNA provirus) ---
    ("HIV1_HXB2",     "NC_001802.1", "retrovirus_int",  "integrating", "Human immunodeficiency virus 1 (HXB2)"),
    ("HIV2_BEN",      "NC_001722.1", "retrovirus_int",  "integrating", "Human immunodeficiency virus 2"),
    ("HTLV1",         "NC_001436.1", "retrovirus_int",  "onco_int",    "Human T-lymphotropic virus 1 (ATLL leukemia)"),
    ("HTLV2",         "NC_001488.1", "retrovirus_int",  "integrating", "Human T-lymphotropic virus 2"),
    ("SIV",           "NC_001549.1", "retrovirus_int",  "integrating", "Simian immunodeficiency virus"),
    ("MoMLV",         "NC_001501.1", "retrovirus_int",  "onco_int",    "Moloney murine leukemia virus"),
    # --- Papillomaviruses & Polyomaviruses (DNA oncoviruses + benign controls) ---
    ("HPV16",         "NC_001526.4", "dna_onco_small",  "onco_int",    "HPV-16 high-risk cervical/HNSCC oncovirus"),
    ("HPV18",         "NC_001357.1", "dna_onco_small",  "onco_int",    "HPV-18 high-risk oncovirus"),
    ("HPV6b",         "NC_001355.1", "dna_benign_ctrl", "non_onco",    "HPV-6b low-risk benign papillomavirus control"),
    ("HPV11",         "NC_001525.1", "dna_benign_ctrl", "non_onco",    "HPV-11 low-risk benign papillomavirus control"),
    ("MCPyV",         "NC_010277.2", "dna_onco_small",  "onco_int",    "Merkel cell polyomavirus (integrating MCC)"),
    ("SV40",          "NC_001669.1", "dna_onco_small",  "onco_model",  "Simian virus 40 (classic p53/RB1 model)"),
    ("JCV",           "NC_001699.1", "dna_benign_ctrl", "non_onco",    "JC polyomavirus"),
    ("BKV",           "NC_001538.1", "dna_benign_ctrl", "non_onco",    "BK polyomavirus"),
    # --- Hepatitis (HCC liver cancer: integrating HBV vs non-integrating RNA HCV) ---
    ("HBV",           "NC_003977.2", "hepatitis_hcc",   "onco_int",    "Hepatitis B virus (dsDNA-RT, TERT/MYC integrator)"),
    ("HCV",           "NC_004102.1", "hepatitis_hcc",   "onco_rna",    "Hepatitis C virus (+ssRNA, non-integrating HCC)"),
    # --- Herpesviruses (Large dsDNA oncoviruses + non-onco controls) ---
    ("EBV_HHV4",      "NC_007605.1", "herpes_large",    "onco_episome","Epstein-Barr virus (Burkitt/NPC, ~172 kb)"),
    ("KSHV_HHV8",     "NC_009333.1", "herpes_large",    "onco_episome","Kaposi sarcoma herpesvirus (~138 kb)"),
    ("HSV1",          "NC_001806.2", "herpes_large",    "non_onco",    "Herpes simplex virus 1 (~152 kb, non-onco)"),
    ("HCMV_HHV5",     "NC_006273.2", "herpes_large",    "non_onco",    "Human cytomegalovirus (~236 kb)"),
    # --- Gene-Therapy / Transformation DNA Viruses ---
    ("AAV2",          "NC_001401.2", "vector_dna",      "site_int",    "Adeno-associated virus 2 (AAVS1 chr19 integrator)"),
    ("HAdV5",         "NC_001405.1", "vector_dna",      "onco_model",  "Human adenovirus 5 (~36 kb)"),
    # --- Acute Non-Integrating RNA Viruses (Controls) ---
    ("SARS_CoV_2",    "NC_045512.2", "rna_acute_ctrl",  "non_int_rna", "SARS-CoV-2 Wuhan-Hu-1 (~29.9 kb)"),
    ("SARS_CoV_1",    "NC_004718.3", "rna_acute_ctrl",  "non_int_rna", "SARS-CoV-1 Tor2 (~29.7 kb)"),
    ("MERS_CoV",      "NC_019843.3", "rna_acute_ctrl",  "non_int_rna", "MERS-CoV (~30.1 kb)"),
    ("Ebola_Zaire",   "NC_002549.1", "rna_acute_ctrl",  "non_int_rna", "Zaire ebolavirus (-ssRNA, ~19 kb)"),
    ("Measles",       "NC_001498.1", "rna_acute_ctrl",  "non_int_rna", "Measles morbillivirus (-ssRNA, ~15.9 kb)"),
    ("Poliovirus",    "NC_002058.3", "rna_acute_ctrl",  "non_int_rna", "Poliovirus (+ssRNA, ~7.4 kb)"),
    ("Zika",          "NC_012532.1", "rna_acute_ctrl",  "non_int_rna", "Zika flavivirus (+ssRNA, ~10.8 kb)"),
    ("Dengue2",       "NC_001474.2", "rna_acute_ctrl",  "non_int_rna", "Dengue virus 2 (+ssRNA, ~10.7 kb)"),
    # --- Bacteriophage Outgroup Controls ---
    ("Phage_Lambda",  "NC_001416.1", "phage_outgroup",  "phage_int",   "Enterobacteria phage lambda (~48.5 kb)"),
    ("Phage_T4",      "NC_000866.4", "phage_outgroup",  "phage_lytic", "Escherichia phage T4 (~169 kb)"),
    ("Phage_PhiX174", "NC_001422.1", "phage_outgroup",  "phage_ssdna", "Escherichia phage phiX174 (~5.4 kb)"),
]
VIRAL_CATALOG_33 = {
    name: {"acc": acc, "label": name, "vclass": vcls, "mode": mode, "notes": notes}
    for name, acc, vcls, mode, notes in _VIRAL_33_SPEC
}

_ONCO_28_SPEC = [
    # --- Core Tumor Suppressors ---
    ("TP53",   "NG_017013.2", "chr17", "tumor_suppressor", "Guardian of genome; 86% mutations in exons 5-8"),
    ("TP63",   "NG_007550.1", "chr3",  "tumor_suppressor", "p53 family epithelial regulator"),
    ("TP73",   "NG_017035.2", "chr1",  "tumor_suppressor", "p53 family suppressor"),
    ("BRCA1",  "NG_005905.2", "chr17", "tumor_suppressor", "HR repair; high Alu density / Alu deletions"),
    ("BRCA2",  "NG_012772.3", "chr13", "tumor_suppressor", "HR repair; AT-rich chr13"),
    ("PTEN",   "NG_007466.2", "chr10", "tumor_suppressor", "PI3K/AKT phosphatase suppressor"),
    ("RB1",    "NG_009009.1", "chr13", "tumor_suppressor", "Retinoblastoma cell-cycle gatekeeper"),
    ("APC",    "NG_008481.4", "chr5",  "tumor_suppressor", "Colorectal Wnt gatekeeper"),
    ("VHL",    "NG_008212.3", "chr3",  "tumor_suppressor", "Renal carcinoma; Alu-mediated deletion hotspot"),
    ("MLH1",   "NG_007109.2", "chr3",  "dna_repair",       "Mismatch repair / Lynch MSI locus"),
    ("CDKN2A", "NG_007485.1", "chr9",  "tumor_suppressor", "p16/ARF melanoma & pan-cancer locus"),
    ("ATM",    "NG_009830.1", "chr11", "dna_repair",       "DSB checkpoint kinase"),
    # --- Major Oncogenes & Viral Integration / APOBEC / ecDNA Hotspots ---
    ("MYC",    "NG_007161.2", "chr8",  "oncogene_viral",   "Master oncogene; HPV-16/18, HBV, EBV hotspot"),
    ("TERT",   "NG_009265.1", "chr5",  "oncogene_viral",   "Telomerase; #1 HBV liver integration hotspot"),
    ("KRAS",   "NG_007524.2", "chr12", "oncogene",         "Pancreatic/lung/colorectal GTPase driver"),
    ("NRAS",   "NG_007572.1", "chr1",  "oncogene",         "Melanoma/AML RAS driver"),
    ("HRAS",   "NG_007666.1", "chr11", "oncogene",         "Bladder/HNSCC RAS driver"),
    ("PIK3CA", "NG_012113.2", "chr3",  "oncogene_apobec",  "Canonical APOBEC3A/B kataegis hotspot (E542K/E545K)"),
    ("EGFR",   "NG_007726.3", "chr7",  "oncogene_ecdna",   "Lung/GBM driver & classic ecDNA amplicon"),
    ("BRAF",   "NG_007873.3", "chr7",  "oncogene_uv",      "Melanoma V600E driver"),
    ("ERBB2",  "NG_007503.1", "chr17", "oncogene",         "HER2 breast cancer amplicon on chr17"),
    ("BCL2",   "NG_009361.1", "chr18", "oncogene",         "Follicular lymphoma t(14;18) driver"),
    ("ALK",    "NG_009445.1", "chr2",  "oncogene",         "Lung/neuroblastoma rearrangement locus"),
    ("FHIT",   "NG_007551.2",     "chr3", "fragile_viral", "FRA3B #1 common fragile site & HPV integration site"),
    # --- Non-Cancer Housekeeping / Structural Controls ---
    ("GAPDH",  "NG_007073.2", "chr12", "control_hk",       "Housekeeping glycolytic control"),
    ("ACTB",   "NG_007992.1", "chr7",  "control_hk",       "Housekeeping beta-actin control"),
    ("HBB",    "NG_000007.3", "chr11", "control_nononco",  "Beta-globin locus control"),
    ("APOE",   "NG_007084.2", "chr19", "control_gc19",     "GC-rich chr19 non-cancer control"),
]
ONCO_LOCI_28 = {
    name: {"acc": acc, "label": name, "chr": ch, "role": role, "notes": notes}
    for name, acc, ch, role, notes in _ONCO_28_SPEC
}


def _build_cosmic_rate_tables():
    """[THM-V3] Exact deterministic (64, 4) trinucleotide context -> target base rate tables
    for 7 canonical COSMIC mutational signatures (A=0, C=1, G=2, T=3).
    Every row sums to <= 1.0, guaranteeing strict nonnegativity of B + eps*J for eps <= 1/6."""
    def _empty():
        return np.zeros((64, 4), dtype=np.float64)

    def _ctx(a, b, c):
        return (a << 4) | (b << 2) | c

    tables = {}

    # 1. SBS1 (Symmetric CpG deamination: NCG -> NTG and CGN -> CAN)
    t = _empty()
    for a in range(4):
        t[_ctx(a, 1, 2), 3] = 1.0
    for c in range(4):
        t[_ctx(1, 2, c), 0] = 1.0
    tables["SBS1_CpG_sym"] = {"rates": t, "strand_mode": "symmetric", "desc": "Clock-like CpG C>T / G>A"}

    # 2. SBS2 APOBEC3A/B kataegis (+ strand only: TCA/TCT -> TTA/TTT)
    t = _empty()
    for c in (0, 3):
        t[_ctx(3, 1, c), 3] = 1.0
    tables["SBS2_APOBEC_plus"] = {"rates": t, "strand_mode": "strand_coord", "desc": "APOBEC3 kataegis TCW>TTW (+ strand)"}

    # 3. SBS13 APOBEC3A/B kataegis (+ strand only: TCA/TCT -> TGA/TGT)
    t = _empty()
    for c in (0, 3):
        t[_ctx(3, 1, c), 2] = 1.0
    tables["SBS13_APOBEC_plus"] = {"rates": t, "strand_mode": "strand_coord", "desc": "APOBEC3 kataegis TCW>TGW (+ strand)"}

    # 4. SBS2 APOBEC symmetric control (TCW>TTW and WGA>WAA equally)
    t = _empty()
    for c in (0, 3):
        t[_ctx(3, 1, c), 3] = 0.5
    for a in (0, 3):
        t[_ctx(a, 2, 0), 0] = 0.5
    tables["SBS2_APOBEC_sym"] = {"rates": t, "strand_mode": "symmetric", "desc": "APOBEC3 symmetric duplex control"}

    # 5. SBS7a UV dipyrimidine (+ strand only: TC/CC -> TT/CT)
    t = _empty()
    for a in (1, 3):
        for c in range(4):
            t[_ctx(a, 1, c), 3] = 1.0
    tables["SBS7a_UV_plus"] = {"rates": t, "strand_mode": "strand_coord", "desc": "UV melanoma YC>YT (+ strand)"}

    # 6. SBS4 Tobacco benzo[a]pyrene (+ strand only: G -> T)
    t = _empty()
    for a in range(4):
        for c in range(4):
            t[_ctx(a, 2, c), 3] = 1.0
    tables["SBS4_Tobacco_plus"] = {"rates": t, "strand_mode": "strand_coord", "desc": "Tobacco lung G>T (+ strand)"}

    # 7. SBS4 Tobacco symmetric control (G>T and C>A equally)
    t = _empty()
    for a in range(4):
        for c in range(4):
            t[_ctx(a, 2, c), 3] = 0.5
            t[_ctx(a, 1, c), 0] = 0.5
    tables["SBS4_Tobacco_sym"] = {"rates": t, "strand_mode": "symmetric", "desc": "Tobacco G>T / C>A symmetric control"}

    return tables


COSMIC_SIG_TABLES = _build_cosmic_rate_tables()


def _cos_mat_frob(A, B):
    na = float(np.linalg.norm(A, ord="fro"))
    nb = float(np.linalg.norm(B, ord="fro"))
    if na <= 1e-30 or nb <= 1e-30:
        return None
    return float(np.sum(A * B) / (na * nb))


def _cos_wmp_unbiased_vs_ref(B_obj, W_ref, z_gate=3.0):
    # [THM-V18] De-attenuated cosine cos_unb = <W_-+^obj, w_hat_ref> / sqrt(u_-+(obj))
    # with strict z_iid(W_-+) >= z_gate stability gate and analytical delta-method SE/95% CI.
    kv = klein_v4_pr2_audit(B_obj)
    if kv is None or W_ref is None:
        return None
    n_ref = float(np.linalg.norm(W_ref, ord="fro"))
    if n_ref <= 1e-30:
        return None
    w_hat = np.asarray(W_ref, dtype=np.float64) / n_ref
    W_obj = kv["_W_mp"]
    n_obj = float(np.linalg.norm(W_obj, ord="fro"))
    cos_raw = float(np.sum(W_obj * w_hat) / n_obj) if n_obj > 1e-30 else None
    u_mp = float(kv["norm2_W_mp_unb"])
    z_mp = kv.get("z_iid_W_mp")
    sd_u = float(kv.get("sd_iid_W_mp") or 0.0)
    deatt = math.sqrt(kv["norm2_W_mp"] / u_mp) if u_mp > 1e-36 else None
    if u_mp <= 1e-36 or z_mp is None or float(z_mp) < float(z_gate):
        return {
            "cos_raw": cos_raw,
            "cos_unb": None,
            "se_cos_unb": None,
            "ci95_lo": None,
            "ci95_hi": None,
            "deattenuation_factor": deatt,
            "z_iid_W_mp": float(z_mp) if z_mp is not None else None,
            "z_gate": float(z_gate),
            "status": "unresolved",
        }
    N = float(kv["N"])
    P = np.asarray(B_obj, dtype=np.float64) / N
    num = float(np.sum(W_obj * w_hat))
    cos_unb = num / math.sqrt(u_mp)
    var_num = max(float(np.sum(P * (w_hat * w_hat))) - num * num, 0.0) / max(N, 1.0)
    var_cos = (var_num / u_mp) + (cos_unb * cos_unb * (sd_u * sd_u)) / (4.0 * u_mp * u_mp)
    se_cos = math.sqrt(max(var_cos, 0.0))
    return {
        "cos_raw": cos_raw,
        "cos_unb": float(cos_unb),
        "se_cos_unb": float(se_cos),
        "ci95_lo": float(cos_unb - 1.959963984540054 * se_cos),
        "ci95_hi": float(cos_unb + 1.959963984540054 * se_cos),
        "deattenuation_factor": float(deatt) if deatt is not None else None,
        "z_iid_W_mp": float(z_mp),
        "z_gate": float(z_gate),
        "status": "resolved",
    }


def benjamini_hochberg_fdr(p_values):
    # [THM-V18] Exact deterministic step-up Benjamini-Hochberg adjusted q-values.
    m = len(p_values)
    if m == 0:
        return []
    indexed = [(i, float(p) if _finite(p) else 1.0) for i, p in enumerate(p_values)]
    indexed.sort(key=lambda t: t[1])
    q_sorted = [0.0] * m
    running = 1.0
    for rank_idx in range(m - 1, -1, -1):
        rank = rank_idx + 1
        val = indexed[rank_idx][1] * float(m) / float(rank)
        running = min(running, val)
        q_sorted[rank_idx] = min(1.0, max(0.0, running))
    out = [1.0] * m
    for (orig_i, _), qv in zip(indexed, q_sorted):
        out[orig_i] = qv
    return out




_P_U = np.arange(NC * NC, dtype=np.int32) >> 6
_P_V = np.arange(NC * NC, dtype=np.int32) & 63
_PAIR_T = (_P_V << 6) | _P_U
_PAIR_R = (MIRROR[_P_V].astype(np.int32) << 6) | MIRROR[_P_U].astype(np.int32)
_PAIR_S = (MIRROR[_P_U].astype(np.int32) << 6) | MIRROR[_P_V].astype(np.int32)


def _resolve_refseqgene_acc(gene_symbol, expected_chr=None, timeout=60):
    """Deterministically resolves the exact human RefSeqGene NG_ accession via NCBI E-utilities."""
    eutils = "https:/" + "/eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    term = f'"{gene_symbol}"[Gene Name] AND "Homo sapiens"[Organism] AND RefSeqGene[Keyword]'
    url = eutils + "esearch.fcgi?db=nuccore&retmode=json&retmax=15&term=" + urllib.parse.quote(term)
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        js = json.loads(r.read().decode("utf-8", "ignore"))
    uids = (js.get("esearchresult") or {}).get("idlist") or []
    if not uids:
        raise RuntimeError(f"NCBI esearch returned no RefSeqGene UIDs for {gene_symbol}")
    s_url = eutils + "esummary.fcgi?db=nuccore&retmode=json&id=" + ",".join(uids[:10])
    req2 = urllib.request.Request(s_url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req2, timeout=timeout) as r2:
        sjs = json.loads(r2.read().decode("utf-8", "ignore"))
    res = sjs.get("result") or {}
    chr_num = str(expected_chr).replace("chr", "") if expected_chr else ""
    best_acc = None
    for uid in uids[:10]:
        rec = res.get(str(uid)) or {}
        acc = str(rec.get("accessionversion") or "")
        title = str(rec.get("title") or "")
        if not acc.startswith("NG_"):
            continue
        if f"({gene_symbol})" in title or f"({gene_symbol}@)" in title:
            if (not chr_num) or (f"chromosome {chr_num}" in title.lower()):
                return acc
            if best_acc is None:
                best_acc = acc
    if best_acc:
        return best_acc
    raise RuntimeError(f"Could not verify NG_ RefSeqGene accession for {gene_symbol} ({expected_chr})")


def exact_signature_jacobian_diag(nt, rate_tbl, site_block=40000):
    """[THM-V3.2] Exact O(L) Zero-RNG per-site diagonal self-term sum_t ||Pi_eps(Delta_t)||^2
    for U-statistic debiasing of the mutagenesis Jacobian J = sum_t Delta_t."""
    L = int(nt.size)
    zero_d = {"W_pp": 0.0, "W_pm": 0.0, "W_mp": 0.0, "W_mm": 0.0}
    if L < 8:
        return zero_d
    # Fast exact vectorized path when each trinucleotide context maps to at most 1 target base (all 7 COSMIC tables)
    row_nz = np.count_nonzero(rate_tbl > 0.0, axis=1)
    assert int(row_nz.max()) <= 1, "rate_tbl must have <= 1 target base per trinucleotide context"
    target_b = np.argmax(rate_tbl, axis=1).astype(np.int32)
    target_w = rate_tbl[np.arange(64), target_b]

    max_i = L - 6  # valid interior 6-mer windows i in [1, L-7]
    win_ok = np.zeros(L, dtype=bool)
    win_ok[1:max_i] = (
        (nt[0:max_i - 1] < 4) & (nt[1:max_i] < 4) & (nt[2:max_i + 1] < 4)
        & (nt[3:max_i + 2] < 4) & (nt[4:max_i + 3] < 4) & (nt[5:max_i + 4] < 4)
        & (nt[6:max_i + 5] < 4) & (nt[7:max_i + 6] < 4)
    )
    pair_all = np.zeros(L, dtype=np.int32)
    v_idx = np.flatnonzero(win_ok)
    if v_idx.size == 0:
        return zero_d
    pair_all[v_idx] = (
        (nt[v_idx].astype(np.int32) << 10)
        | (nt[v_idx + 1].astype(np.int32) << 8)
        | (nt[v_idx + 2].astype(np.int32) << 6)
        | (nt[v_idx + 3].astype(np.int32) << 4)
        | (nt[v_idx + 4].astype(np.int32) << 2)
        | nt[v_idx + 5].astype(np.int32)
    )

    # All candidate mutable sites t in [1, L-2]
    site_ok = (nt[0:L - 2] < 4) & (nt[1:L - 1] < 4) & (nt[2:L] < 4)
    t_cand = np.flatnonzero(site_ok) + 1
    if t_cand.size == 0:
        return zero_d
    ctx_t = (nt[t_cand - 1].astype(np.int32) << 4) | (nt[t_cand].astype(np.int32) << 2) | nt[t_cand + 1].astype(np.int32)
    w_t = target_w[ctx_t]
    elig = w_t > 0.0
    if not np.any(elig):
        return zero_d
    t_el = t_cand[elig]
    w_el = w_t[elig]
    bn_el = target_b[ctx_t[elig]]

    q_I = q_T = q_R = q_S = 0.0
    shifts = np.array([10, 8, 6, 4, 2, 0], dtype=np.int32)
    clear_masks = ~(3 << shifts)
    ks = np.arange(6, dtype=np.int64)[None, :]  # (1, 6)

    for b0 in range(0, t_el.size, int(site_block)):
        b1 = min(t_el.size, b0 + int(site_block))
        tb = t_el[b0:b1]
        wb = w_el[b0:b1]
        bnb = bn_el[b0:b1]
        ik = tb[:, None] - ks  # (B, 6): window indices containing site tb at position k
        in_rng = (ik >= 1) & (ik < max_i)
        ik_safe = np.where(in_rng, ik, 0)
        m6 = in_rng & win_ok[ik_safe]
        if not np.any(m6):
            continue
        p_old = pair_all[ik_safe]
        p_new = (p_old & clear_masks[None, :]) | (bnb[:, None] << shifts[None, :])
        cI = np.concatenate([p_old, p_new], axis=1)  # (B, 12)
        cT = _PAIR_T[cI]
        cR = _PAIR_R[cI]
        cS = _PAIR_S[cI]
        m_f = m6.astype(np.float64)
        s12 = np.concatenate([-m_f, m_f], axis=1)  # (B, 12)
        S_mat = s12[:, :, None] * s12[:, None, :]  # (B, 12, 12)
        w2 = wb * wb
        q_I += float(np.sum(w2 * np.sum(S_mat * (cI[:, :, None] == cI[:, None, :]), axis=(1, 2))))
        q_T += float(np.sum(w2 * np.sum(S_mat * (cI[:, :, None] == cT[:, None, :]), axis=(1, 2))))
        q_R += float(np.sum(w2 * np.sum(S_mat * (cI[:, :, None] == cR[:, None, :]), axis=(1, 2))))
        q_S += float(np.sum(w2 * np.sum(S_mat * (cI[:, :, None] == cS[:, None, :]), axis=(1, 2))))

    return {
        "W_pp": 0.25 * (q_I + q_T + q_R + q_S),
        "W_pm": 0.25 * (q_I + q_T - q_R - q_S),
        "W_mp": 0.25 * (q_I - q_T + q_R - q_S),
        "W_mm": 0.25 * (q_I - q_T - q_R + q_S),
    }


def self_test_virome_onco_v25(verbose=False):
    """[THM-V3.2 & V8.2] 2 additional exact Zero-RNG self-tests for v2.5.0 (65/65 total)."""
    ok = {}
    L = 160
    _ix = np.arange(L, dtype=np.uint64) + np.uint64(777)
    _zz = (_ix ^ (_ix >> np.uint64(30))) * np.uint64(0xbf58476d1ce4e5b9)
    _zz = (_zz ^ (_zz >> np.uint64(27))) * np.uint64(0x94d049bb133111eb)
    s_tiny = ((_zz ^ (_zz >> np.uint64(31))) & np.uint64(3)).astype(np.uint8)
    s_tiny[slice(70, 74)] = 255  # include an N-gap to test exact boundary masking

    # 1. Verify vectorized exact_signature_jacobian_diag against exhaustive site-by-site brute force
    rate_ap = COSMIC_SIG_TABLES["SBS2_APOBEC_plus"]["rates"]
    d_vec = exact_signature_jacobian_diag(s_tiny, rate_ap, site_block=17)
    B0_v, J_v, _ = exact_signature_jacobian(s_tiny, rate_ap)
    d_brute = {"W_pp": 0.0, "W_pm": 0.0, "W_mp": 0.0, "W_mm": 0.0}
    max_i = L - 6
    B0_in = np.zeros((NC, NC), dtype=np.float64)
    for ii in range(1, max_i):
        if int(s_tiny[ii - 1:ii + 7].max()) < 4:
            u = (int(s_tiny[ii]) << 4) | (int(s_tiny[ii + 1]) << 2) | int(s_tiny[ii + 2])
            v = (int(s_tiny[ii + 3]) << 4) | (int(s_tiny[ii + 4]) << 2) | int(s_tiny[ii + 5])
            B0_in[u, v] += 1.0
    assert float(np.max(np.abs(B0_v - B0_in))) == 0.0
    for t_pos in range(1, L - 1):
        if int(s_tiny[t_pos - 1:t_pos + 2].max()) >= 4:
            continue
        ctx = (int(s_tiny[t_pos - 1]) << 4) | (int(s_tiny[t_pos]) << 2) | int(s_tiny[t_pos + 1])
        for b_new in range(4):
            w = float(rate_ap[ctx, b_new])
            if w <= 0.0:
                continue
            s_mut = s_tiny.copy()
            s_mut[t_pos] = b_new
            B_mut = np.zeros((NC, NC), dtype=np.float64)
            for ii in range(1, max_i):
                if int(s_tiny[ii - 1:ii + 7].max()) < 4:
                    u = (int(s_mut[ii]) << 4) | (int(s_mut[ii + 1]) << 2) | int(s_mut[ii + 2])
                    v = (int(s_mut[ii + 3]) << 4) | (int(s_mut[ii + 4]) << 2) | int(s_mut[ii + 5])
                    B_mut[u, v] += 1.0
            D = w * (B_mut - B0_in)
            DT = D.T
            DR = D.T[np.ix_(MIRROR, MIRROR)]
            DS = D[np.ix_(MIRROR, MIRROR)]
            d_brute["W_pp"] += float(np.sum((0.25 * (D + DT + DR + DS)) ** 2))
            d_brute["W_pm"] += float(np.sum((0.25 * (D + DT - DR - DS)) ** 2))
            d_brute["W_mp"] += float(np.sum((0.25 * (D - DT + DR - DS)) ** 2))
            d_brute["W_mm"] += float(np.sum((0.25 * (D - DT - DR + DS)) ** 2))
    diag_err = max(abs(d_vec[k] - d_brute[k]) for k in d_brute)
    assert diag_err < 1e-12, (diag_err, d_vec, d_brute)
    ok["thmV3_2_jacobian_diag_vs_brute_err"] = diag_err

    # 2. Verify U-statistic chi_mm_unb removes O(1/N) length bias (8x 2-kb blocks vs pooled 16-kb)
    _ix2 = np.arange(16000, dtype=np.uint64) + np.uint64(2026)
    _z2 = (_ix2 ^ (_ix2 >> np.uint64(30))) * np.uint64(0xbf58476d1ce4e5b9)
    _z2 = (_z2 ^ (_z2 >> np.uint64(27))) * np.uint64(0x94d049bb133111eb)
    _s16 = ((_z2 ^ (_z2 >> np.uint64(31))) & np.uint64(3)).astype(np.uint8)
    _blks = [_s16[slice(k_ * 2000, (k_ + 1) * 2000)] for k_ in range(8)]
    _gap10 = np.full(10, 255, dtype=np.uint8)
    _spooled = np.concatenate([x_ for blk_ in _blks for x_ in (blk_, _gap10)])
    _tbl_ap = {"SBS2_APOBEC_plus": COSMIC_SIG_TABLES["SBS2_APOBEC_plus"]}
    _mb = [analyze_mutagenesis_signatures(blk_, _tbl_ap)["SBS2_APOBEC_plus"] for blk_ in _blks]
    _mp = analyze_mutagenesis_signatures(_spooled, _tbl_ap)["SBS2_APOBEC_plus"]
    _raw_b = float(np.mean([d_["chi_mm_x1e4"] for d_ in _mb]))
    _unb_b = float(np.mean([d_["chi_mm_unb_x1e4"] for d_ in _mb]))
    rel_diff_raw = abs(_raw_b - _mp["chi_mm_x1e4"]) / max(abs(_mp["chi_mm_x1e4"]), 1e-12)
    rel_diff_unb = abs(_unb_b - _mp["chi_mm_unb_x1e4"]) / max(abs(_mp["chi_mm_unb_x1e4"]), 1e-12)
    assert rel_diff_raw > 0.50 and rel_diff_unb < 0.25 and rel_diff_unb < 0.25 * rel_diff_raw, (
        rel_diff_unb, rel_diff_raw, _unb_b, _mp["chi_mm_unb_x1e4"], _raw_b, _mp["chi_mm_x1e4"]
    )
    ok["thmV3_2_unbiased_chi_length_invariance_ok"] = True
    if verbose:
        print(json.dumps(_clean(ok), indent=2))
    return ok


def exact_signature_jacobian(nt, rate_tbl, chunk_size=2_000_000):
    """[THM-V3] Exact O(L) Zero-RNG first-order mutagenesis Jacobian J in R^{64x64}
    for lag d=3 bp over all interior valid 8-mers nt[i-1 : i+7] (1 <= i <= L-7).
    Satisfies sum(J) == 0.0 identically and equals the exact exhaustive sum of
    single-base substitutions weighted by rate_tbl[ctx(t), b_new]."""
    L = int(nt.size)
    B0 = np.zeros((NC, NC), dtype=np.float64)
    J_flat = np.zeros(NC * NC, dtype=np.float64)
    if L < 8:
        return B0, J_flat.reshape(NC, NC), 0.0
    max_i = L - 6  # i in [1, L-7], uses nt[i-1 : i+7]
    mut_mass = 0.0
    for s0 in range(1, max_i, int(chunk_size)):
        s1 = min(max_i, s0 + int(chunk_size))
        span = nt[s0 - 1:s1 + 6]
        if span.size < 8:
            continue
        ok = (
            (nt[s0 - 1:s1 - 1] < 4) & (nt[s0:s1] < 4) & (nt[s0 + 1:s1 + 1] < 4)
            & (nt[s0 + 2:s1 + 2] < 4) & (nt[s0 + 3:s1 + 3] < 4)
            & (nt[s0 + 4:s1 + 4] < 4) & (nt[s0 + 5:s1 + 5] < 4)
            & (nt[s0 + 6:s1 + 6] < 4)
        )
        if not np.any(ok):
            continue
        idx = np.flatnonzero(ok) + s0
        b_arr = [nt[idx + p].astype(np.int32) for p in range(-1, 7)]  # b_arr[k] = nt[idx + k - 1]
        pair_old = (
            (b_arr[1] << 10) | (b_arr[2] << 8) | (b_arr[3] << 6)
            | (b_arr[4] << 4) | (b_arr[5] << 2) | b_arr[6]
        )
        B0 += np.bincount(pair_old, minlength=NC * NC).reshape(NC, NC)
        for p in range(6):
            ctx = (b_arr[p] << 4) | (b_arr[p + 1] << 2) | b_arr[p + 2]
            shift = 2 * (5 - p)
            clear_mask = ~(3 << shift)
            for b_new in range(4):
                w = rate_tbl[ctx, b_new]
                m = w > 0.0
                if not np.any(m):
                    continue
                wm = w[m]
                po = pair_old[m]
                pn = (po & clear_mask) | (b_new << shift)
                J_flat += np.bincount(pn, weights=wm, minlength=NC * NC)
                J_flat -= np.bincount(po, weights=wm, minlength=NC * NC)
                mut_mass += float(wm.sum())
    return B0, J_flat.reshape(NC, NC), mut_mass / 6.0


def analyze_mutagenesis_signatures(nt, sig_tables=None, eps_eval=0.05, compute_denovo=True):
    """[THM-V3 & THM-V3.3] Evaluates exact Zero-RNG Klein V_4 linear response,
    U-statistic unbiased susceptibility chi_--^unb, and exact de-novo vs substrate-coupled
    chiral decomposition J = J_denovo + J_sub across all canonical COSMIC signatures."""
    if sig_tables is None:
        sig_tables = COSMIC_SIG_TABLES
    rc = reverse_complement_nt(nt) if compute_denovo else None
    out = {}
    B0_ref = None
    kv0 = None
    for sname, smeta in sig_tables.items():
        B0, J, exp_muts = exact_signature_jacobian(nt, smeta["rates"])
        N = float(B0.sum())
        if N <= 10.0:
            out[sname] = {"skipped": True}
            continue
        if kv0 is None:
            B0_ref = B0
            kv0 = klein_v4_pr2_audit(B0_ref)
        dM = J / N
        TP = dM.T
        RP = dM.T[np.ix_(MIRROR, MIRROR)]
        SP = dM[np.ix_(MIRROR, MIRROR)]
        J_pp = 0.25 * (dM + TP + RP + SP)
        J_pm = 0.25 * (dM + TP - RP - SP)
        J_mp = 0.25 * (dM - TP + RP - SP)
        J_mm = 0.25 * (dM - TP - RP + SP)
        n_pp = float(np.sum(J_pp * J_pp))
        n_pm = float(np.sum(J_pm * J_pm))
        n_mp = float(np.sum(J_mp * J_mp))
        n_mm = float(np.sum(J_mm * J_mm))
        dir_sum = n_mp + n_mm
        pr2_sum = n_pm + n_mm
        diag_k = exact_signature_jacobian_diag(nt, smeta["rates"])
        denom_u = max(N * (N - 1.0), 1.0)
        u_mp = (N * N * n_mp - diag_k["W_mp"]) / denom_u
        u_mm = (N * N * n_mm - diag_k["W_mm"]) / denom_u
        u_dir = max(u_mp, 0.0) + max(u_mm, 0.0)

        # THM-V3.3: Exact de novo vs substrate-coupled chiral split
        chi_dn_raw = chi_dn_unb = chi_sub_raw = dn_share = None
        if compute_denovo and rc is not None:
            _, J_rc, _ = exact_signature_jacobian(rc, smeta["rates"])
            J_dn = 0.5 * (J + J_rc) / N
            J_sb = 0.5 * (J - J_rc) / N
            J_dn_mm = 0.25 * (J_dn - J_dn.T - J_dn.T[np.ix_(MIRROR, MIRROR)] + J_dn[np.ix_(MIRROR, MIRROR)])
            J_sb_mm = 0.25 * (J_sb - J_sb.T - J_sb.T[np.ix_(MIRROR, MIRROR)] + J_sb[np.ix_(MIRROR, MIRROR)])
            n_dn_mm = float(np.sum(J_dn_mm * J_dn_mm))
            n_sb_mm = float(np.sum(J_sb_mm * J_sb_mm))
            chi_dn_raw = 1e4 * n_dn_mm
            chi_sub_raw = 1e4 * n_sb_mm
            if smeta["strand_mode"] == "symmetric":
                u_dn_mm = n_dn_mm  # identically 0.0 to machine precision
            else:
                diag_rc = exact_signature_jacobian_diag(rc, smeta["rates"])
                u_dn_mm = (N * N * n_dn_mm - 0.25 * (diag_k["W_mm"] + diag_rc["W_mm"])) / max(N * (N - 0.5), 1.0)
            chi_dn_unb = 1e4 * u_dn_mm
            dn_share = (n_dn_mm / (n_dn_mm + n_sb_mm)) if (n_dn_mm + n_sb_mm) > 1e-36 else 0.0

        # THM-V18: Exact algebraic baseline 1/2 decomposition for T-odd A = J_-+ + J_--:
        #   ||J_--||^2 / ||A||^2 = 1/2 * (1 - <A, R A> / ||A||^2)
        A_odd = 0.5 * (dM - TP)
        RA_odd = A_odd.T[np.ix_(MIRROR, MIRROR)]
        r_ov_cos = (float(np.sum(A_odd * RA_odd)) / dir_sum) if dir_sum > 1e-36 else 0.0
        chiral_excess_half = -0.5 * r_ov_cos

        B_eps = np.maximum(B0 + float(eps_eval) * J, 0.0)
        kv_eps = klein_v4_pr2_audit(B_eps)
        floor_4n = 1.0 / (4.0 * N)
        delta_wmm_raw = float(kv_eps["norm2_W_mm"] - kv0["norm2_W_mm"])
        chi_eff = u_mm if u_mm > 1e-36 else n_mm
        eps_crit = (1.0 / math.sqrt(4.0 * N * chi_eff)) if chi_eff > 1e-36 else None
        eps_phys = bool(eps_crit is not None and eps_crit <= 1.0)
        mut_e_fl = _fval(kv_eps.get("norm2_W_mm_unb"), 0.0) / max(floor_4n, TINY)
        base_e_fl = _fval(kv0.get("norm2_W_mm_unb"), 0.0) / max(floor_4n, TINY)
        out[sname] = {
            "strand_mode": smeta["strand_mode"],
            "desc": smeta["desc"],
            "N_interior": int(N),
            "eligible_site_frac": float(exp_muts) / N,
            "mass_err": abs(float(J.sum())),
            "J_norm2_W_pp": n_pp,
            "J_norm2_W_pm": n_pm,
            "J_norm2_W_mp": n_mp,
            "J_norm2_W_mm": n_mm,
            "chi_mm_x1e4": 1e4 * n_mm,
            "chi_mm_unb_x1e4": 1e4 * u_mm,
            "chi_mp_unb_x1e4": 1e4 * u_mp,
            "chi_mm_denovo_x1e4": chi_dn_raw,
            "chi_mm_denovo_unb_x1e4": chi_dn_unb,
            "chi_mm_sub_x1e4": chi_sub_raw,
            "denovo_chiral_share": dn_share,
            "J_chiral_share_of_directed_unb": (max(u_mm, 0.0) / u_dir) if u_dir > 1e-36 else 0.0,
            "J_chiral_share_of_directed": (n_mm / dir_sum) if dir_sum > 1e-36 else 0.0,
            "J_R_overlap_cos": float(r_ov_cos),
            "J_chiral_excess_over_half": float(chiral_excess_half),
            "J_PR2_dyn_share": (n_mm / pr2_sum) if pr2_sum > 1e-36 else 0.0,
            "cos_Wmp_native_vs_J": _cos_mat_frob(kv0["_W_mp"], J_mp),
            "eps_eval": float(eps_eval),
            "eps_crit_4N": eps_crit,
            "eps_crit_physical_transition": eps_phys,
            "eps_crit_status": "transition_in_unit_dose" if eps_phys else "no_physical_transition",
            "delta_W_mm_raw_injected": delta_wmm_raw,
            "delta_W_mm_over_human_euch": delta_wmm_raw / 5.86e-11,
            "base_PR2_ok_unb": kv0.get("time_PR2_compliant_share_unb"),
            "mut_PR2_ok_unb": kv_eps.get("time_PR2_compliant_share_unb"),
            "delta_PR2_ok_pct": 100.0 * (
                _fval(kv_eps.get("time_PR2_compliant_share_unb"), 0.0)
                - _fval(kv0.get("time_PR2_compliant_share_unb"), 0.0)
            ),
            "mut_W_mm_energy_over_floor": mut_e_fl,
            "base_W_mm_energy_over_floor": base_e_fl,
            "mut_W_mm_snr_over_4N": mut_e_fl,
            "base_W_mm_snr_over_4N": base_e_fl,
            "delta_W_mp_rel_pct": 100.0 * (
                (_fval(kv_eps.get("norm2_W_mp_unb"), 0.0) - _fval(kv0.get("norm2_W_mp_unb"), 0.0))
                / max(abs(_fval(kv0.get("norm2_W_mp_unb"), 0.0)), 1e-30)
            ),
        }
    return out


def deterministic_proviral_insertion_audit(nt_host, nt_virus):
    """[THM-V2] Exact Zero-RNG deterministic insertion of a real viral genome into a real
    human genomic locus in +, -, and balanced head-to-head (+-) orientations, including
    both host-normalized ΔW and exact dilution-free count-space ΔB operators."""
    Lh, Lv = int(nt_host.size), int(nt_virus.size)
    if Lh < 2000 or Lv < 500:
        return None
    mid = Lh // 2
    rc_v = reverse_complement_nt(nt_virus)
    seq_p = np.concatenate([nt_host[:mid], nt_virus, nt_host[mid:]])
    seq_m = np.concatenate([nt_host[:mid], rc_v, nt_host[mid:]])
    seq_hh = np.concatenate([nt_host[:mid], nt_virus, rc_v, nt_host[mid:]])

    Bh = count_triplet_pairs_lag(nt_host, d=3)
    Bv = count_triplet_pairs_lag(nt_virus, d=3)
    Bp = count_triplet_pairs_lag(seq_p, d=3)
    Bm = count_triplet_pairs_lag(seq_m, d=3)
    Bhh = count_triplet_pairs_lag(seq_hh, d=3)

    kv_h = klein_v4_pr2_audit(Bh)
    kv_v = klein_v4_pr2_audit(Bv)
    kv_p = klein_v4_pr2_audit(Bp)
    kv_m = klein_v4_pr2_audit(Bm)
    kv_hh = klein_v4_pr2_audit(Bhh)

    dWmm_p = kv_p["_W_mm"] - kv_h["_W_mm"]
    dWmm_m = kv_m["_W_mm"] - kv_h["_W_mm"]
    dWmp_p = kv_p["_W_mp"] - kv_h["_W_mp"]
    dWmp_m = kv_m["_W_mp"] - kv_h["_W_mp"]

    # Dilution-free count-space viral injection: ΔB_eps = N_ins * W_eps(ins) - N_h * W_eps(host)
    Bmm_h, Bmm_v, Bmm_p = kv_h["N"] * kv_h["_W_mm"], kv_v["N"] * kv_v["_W_mm"], kv_p["N"] * kv_p["_W_mm"]
    Bmp_h, Bmp_v, Bmp_p = kv_h["N"] * kv_h["_W_mp"], kv_v["N"] * kv_v["_W_mp"], kv_p["N"] * kv_p["_W_mp"]
    dBmm_p = Bmm_p - Bmm_h
    dBmm_m = kv_m["N"] * kv_m["_W_mm"] - Bmm_h
    dBmp_p = Bmp_p - Bmp_h
    dBmp_m = kv_m["N"] * kv_m["_W_mp"] - Bmp_h
    dB_junc_p = (Bp - Bh) - Bv
    junc_share = float(np.sum(dB_junc_p * dB_junc_p)) / max(float(np.sum(Bv * Bv)), 1.0)

    def _concat_3term(Bh_q, Bv_q, Bp_q):
        e_h = float(np.sum(Bh_q * Bh_q))
        e_v = float(np.sum(Bv_q * Bv_q))
        e_hv = 2.0 * float(np.sum(Bh_q * Bv_q))
        Bj_q = Bp_q - Bh_q - Bv_q
        e_j = float(np.sum(Bj_q * Bj_q)) + 2.0 * float(np.sum((Bh_q + Bv_q) * Bj_q))
        e_tot = float(np.sum(Bp_q * Bp_q))
        return {
            "host_E": e_h, "virus_E": e_v, "interference_2_hv": e_hv,
            "junction_E": e_j, "total_ins_E": e_tot,
            "decomp_err": abs(e_tot - (e_h + e_v + e_hv + e_j)),
        }

    fl_p = 1.0 / (4.0 * max(float(kv_p["N"]), 1.0))
    fl_hh = 1.0 / (4.0 * max(float(kv_hh["N"]), 1.0))
    p_ef = _fval(kv_p.get("norm2_W_mm_unb"), 0.0) / fl_p
    hh_ef = _fval(kv_hh.get("norm2_W_mm_unb"), 0.0) / fl_hh
    return {
        "model_note": "in-silico concatenation model (R-involution equivariance per THM-V1)",
        "host_bp": Lh,
        "virus_bp": Lv,
        "virus_fraction": float(Lv) / float(Lh + Lv),
        "host_PR2_ok_unb": kv_h.get("time_PR2_compliant_share_unb"),
        "plus_PR2_ok_unb": kv_p.get("time_PR2_compliant_share_unb"),
        "minus_PR2_ok_unb": kv_m.get("time_PR2_compliant_share_unb"),
        "head2head_PR2_ok_unb": kv_hh.get("time_PR2_compliant_share_unb"),
        "plus_W_mm_energy_over_floor": p_ef,
        "head2head_W_mm_energy_over_floor": hh_ef,
        "plus_W_mm_z_iid": kv_p.get("z_iid_W_mm"),
        "head2head_W_mm_z_iid": kv_hh.get("z_iid_W_mm"),
        "plus_W_mm_snr_over_4N": p_ef,
        "head2head_W_mm_snr_over_4N": hh_ef,
        "chiral_cancellation_ratio_hh": (
            max(_fval(kv_p.get("norm2_W_mm_unb"), 0.0), 1e-30)
            / max(abs(_fval(kv_hh.get("norm2_W_mm_unb"), 0.0)), 1e-30)
        ),
        "cos_dWmm_plus_vs_minus": _cos_mat_frob(dWmm_p, dWmm_m),
        "cos_dWmp_plus_vs_minus": _cos_mat_frob(dWmp_p, dWmp_m),
        "cos_dBmm_plus_vs_minus": _cos_mat_frob(dBmm_p, dBmm_m),
        "cos_dBmp_plus_vs_minus": _cos_mat_frob(dBmp_p, dBmp_m),
        "junction_energy_share": junc_share,
        "concat_energy_3term_W_mm": _concat_3term(Bmm_h, Bmm_v, Bmm_p),
        "concat_energy_3term_W_mp": _concat_3term(Bmp_h, Bmp_v, Bmp_p),
    }


def _load_human_reference_arrows(atlas_npz_path=None, require_all_24=True):
    """Loads exact 24-chromosome & grand-pooled human euchromatin W_-+ tensors
    STRICTLY from SCRIPT_DIR (or explicit --atlas-npz), enforcing 24/24 completeness."""
    candidates = []
    explicit = bool(atlas_npz_path)
    if explicit:
        candidates.append(Path(atlas_npz_path))
    for rel in (
        "human_t2t_euch_ref_v2_8.npz",
        "human_t2t_euch_ref_v2_6.npz",
        "atlas_b3_cache.npz",
        "results_human_t2t_atlas_v1_5_3/atlas_b3_cache.npz",
    ):
        candidates.append((SCRIPT_DIR / rel).resolve())
    for cand in candidates:
        if cand.exists():
            npz = np.load(cand)
            chr_wmp = {}
            B_grand = np.zeros((NC, NC), dtype=np.float64)
            D_grand = np.zeros((4, 4), dtype=np.float64)
            T3_grand = np.zeros(NC, dtype=np.float64)
            has_all_lags = True
            for k in HUMAN_T2T_24:
                ek = f"{k}__euch"
                if ek in npz.files:
                    Bm = np.array(npz[ek], dtype=np.float64)
                    B_grand += Bm
                    kv = klein_v4_pr2_audit(Bm)
                    if kv is not None:
                        chr_wmp[HUMAN_T2T_24[k]["label"]] = kv["_W_mp"]
                dk = f"{k}__euch_dinuc"
                if dk in npz.files:
                    D_grand += np.array(npz[dk], dtype=np.float64)
                tk = f"{k}__euch_trinuc"
                if tk in npz.files:
                    T3_grand += np.array(npz[tk], dtype=np.float64)
                if f"{k}__euch_lag_171" not in npz.files:
                    has_all_lags = False
            npz.close()
            if require_all_24 and (not explicit):
                if len(chr_wmp) < 24 or float(D_grand.sum()) <= 0.0 or float(T3_grand.sum()) <= 0.0 or (not has_all_lags):
                    continue
            if float(B_grand.sum()) > 0.0:
                kv_g = klein_v4_pr2_audit(B_grand)
                P_m1_g = exact_markov1_from_dinuc_counts(D_grand) if float(D_grand.sum()) > 0.0 else exact_markov1_from_pair_table(B_grand)
                kv_m1_g = klein_v4_pr2_audit(P_m1_g)
                P_m2_g = exact_markov2_from_trinuc_counts(T3_grand) if float(T3_grand.sum()) > 0.0 else exact_markov2_from_pair_table(B_grand)
                kv_m2_g = klein_v4_pr2_audit(P_m2_g)
                a11_g = thm_D11_arrow_screen_and_cpg_audit(B_grand)
                if require_all_24 and (not explicit):
                    assert len(chr_wmp) == 24, f"Expected 24 human chr in {cand}, got {len(chr_wmp)}"
                return {
                    "source": str(cand),
                    "source_sha256": sha256_file(cand),
                    "n_chr": len(chr_wmp),
                    "grand_B": B_grand,
                    "grand_W_mp": kv_g["_W_mp"],
                    "grand_W_mp_m1": kv_m1_g["_W_mp"],
                    "grand_W_mp_perp_m1": kv_g["_W_mp"] - kv_m1_g["_W_mp"],
                    "grand_W_mp_m2": kv_m2_g["_W_mp"],
                    "grand_W_mp_perp_m2": kv_g["_W_mp"] - kv_m2_g["_W_mp"],
                    "grand_PR2_ok_unb": kv_g.get("time_PR2_compliant_share_unb"),
                    "grand_nonenv_3600_unb": a11_g.get("pure_nonenv_arrow_share_unb") if a11_g else None,
                    "grand_nonenv_3600_raw": a11_g.get("pure_nonenv_arrow_share") if a11_g else None,
                    "chr_W_mp": chr_wmp,
                }
    return None


def build_local_human_reference(chr_dir, ref_npz_path, auto_download=True, sector_bp=100000):
    """[STANDALONE v2.6.1] Downloads any missing Human T2T-CHM13v2.0 chromosomes into SCRIPT_DIR/chr_cache,
    verifies RefSeq accession headers, applies exact THM-D7 Kittler-Illingworth BIC euchromatin masking,
    and saves a reusable local human_t2t_euch_ref_v2_6.npz (24 chr, 10 lags, dinuc & trinuc) with checkpointing."""
    import gc
    chr_dir = Path(chr_dir)
    chr_dir.mkdir(parents=True, exist_ok=True)
    ref_npz_path = Path(ref_npz_path)
    store = {}
    if ref_npz_path.exists():
        try:
            npz = np.load(ref_npz_path)
            for fn in npz.files:
                store[fn] = np.array(npz[fn])
            npz.close()
        except Exception:
            store = {}

    def _chr_complete(ck):
        return (
            f"{ck}__euch" in store
            and f"{ck}__euch_dinuc" in store
            and f"{ck}__euch_trinuc" in store
            and all(f"{ck}__euch_lag_{int(ld)}" in store for ld in V21_LAGS if int(ld) != 3)
        )

    missing = [k for k in HUMAN_T2T_24 if not _chr_complete(k)]
    if not missing:
        return ref_npz_path

    print(f"  [REF BUILD] Building standalone 24-chr Human T2T euchromatin reference in {ref_npz_path.name} ({len(missing)} chr to process)...")
    for idx_k, k in enumerate(HUMAN_T2T_24, 1):
        if _chr_complete(k):
            continue
        meta = HUMAN_T2T_24[k]
        fp = chr_dir / f"{k}_{meta['acc']}.fna.gz"
        if not fp.exists() and chr_dir.exists():
            hits = list(chr_dir.glob(f"*{meta['acc']}*.fna*"))
            if hits:
                fp = hits[0]
        if not fp.exists() or fp.stat().st_size <= 150:
            if not auto_download:
                print(f"  [REF BUILD skip] {k} not in {chr_dir} (pass --auto-download)")
                continue
            fp = download_chromosome_fasta(k, meta, chr_dir)
        print(f"  [REF BUILD {idx_k:2d}/24] Masking THM-D7 euchromatin on {meta['label']} ({fp.name})...")
        nt, hdr = load_chromosome_nt(fp)
        assert meta["acc"] in hdr, f"T2T chr header mismatch for {k}: expected {meta['acc']}, got {hdr!r}"
        L = int(nt.size)
        n_sec = L // int(sector_bp)
        H3 = np.full(n_sec, np.nan, dtype=np.float64)
        INB = np.full(n_sec, np.nan, dtype=np.float64)
        INB_RAW = np.full(n_sec, np.nan, dtype=np.float64)
        for w in range(n_sec):
            seg = nt[w * sector_bp:(w + 1) * sector_bp]
            if float(np.mean(seg < 4)) < 0.5:
                continue
            c64 = count_triplets_1d(seg, chunk_size=sector_bp)
            p64 = c64 / max(float(c64.sum()), 1.0)
            m64 = p64 > 0
            H3[w] = float(-np.sum(p64[m64] * np.log2(p64[m64])))
            sc = two_sided_screen_full(count_triplet_pairs_lag(seg, d=3, chunk_size=sector_bp))
            if sc is not None:
                INB[w] = sc["I_nonbridge_MM"]
                INB_RAW[w] = sc["I_nonbridge"]
        is_euch, _, _ = partition_euch_sat_parameter_free(H3, INB, INB_RAW=INB_RAW)
        nt_eu = np.full(L, 255, dtype=np.uint8)
        for w in np.flatnonzero(is_euch):
            s0 = int(w) * sector_bp
            nt_eu[s0:min(L, s0 + sector_bp)] = nt[s0:min(L, s0 + sector_bp)]
        del nt
        store[f"{k}__euch"] = count_triplet_pairs_lag(nt_eu, d=3)
        store[f"{k}__euch_dinuc"] = count_dinuc_4x4(nt_eu)
        store[f"{k}__euch_trinuc"] = count_triplets_1d(nt_eu)
        for lag_d in V21_LAGS:
            if int(lag_d) != 3:
                store[f"{k}__euch_lag_{int(lag_d)}"] = count_triplet_pairs_lag(nt_eu, d=int(lag_d))
        del nt_eu
        store["__atlas_version__"] = np.array(VERSION)
        store["__script_sha256__"] = np.array(sha256_file(__file__) or "")
        tmp_npz = ref_npz_path.with_name(ref_npz_path.stem + ".tmp.npz")
        np.savez_compressed(tmp_npz, **store)
        os.replace(tmp_npz, ref_npz_path)
        gc.collect()
    return ref_npz_path


def _spaced_palindrome_spectrum(nt, lags=(3, 6, 9, 15, 30)):
    """[THM-V4] Computes Theorem-5 palindromic anti-diagonal share of W_-+ across micro-lags."""
    out = {}
    pal_i = np.arange(NC)
    pal_j = MIRROR
    for d in lags:
        B = count_triplet_pairs_lag(nt, d=int(d))
        kv = klein_v4_pr2_audit(B)
        if kv is None:
            continue
        Wmp = kv["_W_mp"]
        tot = float(np.sum(Wmp * Wmp))
        pal_e = float(np.sum(Wmp[pal_i, pal_j] ** 2))
        out[str(d)] = {
            "W_mp_unb": kv.get("norm2_W_mp_unb"),
            "p_pal_frac": kv.get("p_pal_6bp"),
            "pal_Wmp_energy_share": (pal_e / tot) if tot > 1e-30 else None,
        }
    return out


def _download_with_retry_v2(key, meta, cache_dir):
    """Rate-limit-safe wrapper around download_chromosome_fasta for short viral/RefSeqGene files."""
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    acc = meta["acc"]
    dest = cache_dir / f"{key}_{acc}.fna.gz"
    if dest.exists() and dest.stat().st_size > 150:
        return dest
    for attempt in range(1, 5):
        try:
            time.sleep(0.35 * attempt)
            return download_chromosome_fasta(key, meta, cache_dir)
        except Exception as e:
            print(f"    [NCBI retry {attempt}/4] {key} ({acc}): {e}")
            time.sleep(2.0 * attempt)
    raise RuntimeError(f"Failed NCBI fetch for {key} ({acc}) after 4 retries")


def self_test_virome_onco(verbose=False):
    """[THM-V1..V4] 6 exact deterministic Zero-RNG self-tests for Paper 5."""
    ok = {}
    L = 6000
    t = np.arange(L, dtype=np.int64)
    nt = ((t * 7 + (t // 3) * 5 + (t * t) // 11) % 4).astype(np.uint8)
    rc = reverse_complement_nt(nt)
    nt_pr2 = np.concatenate([nt, np.full(40, 255, dtype=np.uint8), rc])

    # 1. THM-V1: Reverse-complement preserves W_-+ identically and negates W_-- identically
    kv_a = klein_v4_pr2_audit(count_triplet_pairs_lag(nt, d=3))
    kv_b = klein_v4_pr2_audit(count_triplet_pairs_lag(rc, d=3))
    err_mp = float(np.max(np.abs(kv_b["_W_mp"] - kv_a["_W_mp"])))
    err_mm = float(np.max(np.abs(kv_b["_W_mm"] + kv_a["_W_mm"])))
    assert err_mp == 0.0 and err_mm == 0.0, (err_mp, err_mm)
    ok["thmV1_viral_revcomp_Wmp_inv_Wmm_flip_err"] = max(err_mp, err_mm)

    # 2. THM-V2: Proviral insertion into PR2-exact host -> chiral flare + head-to-head cancellation
    v_syn = np.tile(np.array([2, 2, 0, 1, 3, 2, 0, 0, 3, 1, 2], dtype=np.uint8), 180)
    ins = deterministic_proviral_insertion_audit(nt_pr2, v_syn)
    assert ins["plus_W_mm_snr_over_4N"] > 15.0, ins
    assert ins["head2head_W_mm_snr_over_4N"] < 1.0, ins
    assert ins["cos_dWmm_plus_vs_minus"] < -0.999, ins
    assert ins["cos_dWmp_plus_vs_minus"] > 0.999, ins
    assert ins["chiral_cancellation_ratio_hh"] > 20.0, ins
    ok["thmV2_provirus_chiral_flare_and_hh_cancel_ok"] = True
    ok.update(self_test_virome_onco_v21(verbose=False))

    # 3. THM-V3a: Symmetric SBS1 (CpG) & SBS2_sym preserve PR2 identically on PR2-exact duplex
    mut_pr2 = analyze_mutagenesis_signatures(nt_pr2)
    sbs1 = mut_pr2["SBS1_CpG_sym"]
    ap_sym = mut_pr2["SBS2_APOBEC_sym"]
    assert sbs1["mass_err"] == 0.0 and sbs1["J_norm2_W_mm"] < 1e-34 and sbs1["J_norm2_W_pm"] < 1e-34, sbs1
    assert ap_sym["mass_err"] == 0.0 and ap_sym["J_norm2_W_mm"] < 1e-34 and ap_sym["J_norm2_W_pm"] < 1e-34, ap_sym
    ok["thmV3_sbs1_and_apobec_sym_zero_chiral_err"] = max(sbs1["J_norm2_W_mm"], ap_sym["J_norm2_W_mm"])

    # 4. THM-V3b: Strand-coordinated APOBEC kataegis strictly breaks PR2 and injects chiral W_--
    ap_plus = mut_pr2["SBS2_APOBEC_plus"]
    assert ap_plus["J_norm2_W_mm"] > 1e-4 and ap_plus["J_chiral_share_of_directed"] > 0.2, ap_plus
    assert ap_plus["delta_W_mm_raw_injected"] > 1e-6, ap_plus
    quad_err = abs(ap_plus["delta_W_mm_raw_injected"] - (0.05 ** 2) * ap_plus["J_norm2_W_mm"])
    assert quad_err < 1e-18, quad_err
    assert ap_plus["eps_crit_4N"] is not None and 0.05 < ap_plus["eps_crit_4N"] < 0.25, ap_plus
    ok["thmV3_apobec_kataegis_chiral_injection_ok"] = True
    ok["thmV3_quadratic_chiral_response_err"] = quad_err

    # 5. THM-V3c: Exact algebraic identity of vectorized Jacobian vs exhaustive site-by-site single-hit sum
    s_tiny = nt[:140].copy()
    rate_ap = COSMIC_SIG_TABLES["SBS2_APOBEC_plus"]["rates"]
    B0_v, J_v, _ = exact_signature_jacobian(s_tiny, rate_ap)
    J_brute = np.zeros((NC, NC), dtype=np.float64)
    B0_in = count_triplet_pairs_lag(s_tiny[1:-1], d=3)
    for pos_t in range(1, s_tiny.size - 1):
        ctx = (int(s_tiny[pos_t - 1]) << 4) | (int(s_tiny[pos_t]) << 2) | int(s_tiny[pos_t + 1])
        for b_new in range(4):
            w = float(rate_ap[ctx, b_new])
            if w > 0.0:
                s_mut = s_tiny.copy()
                s_mut[pos_t] = b_new
                J_brute += w * (count_triplet_pairs_lag(s_mut[1:-1], d=3) - B0_in)
    jac_err = float(np.max(np.abs(J_v - J_brute)))
    assert jac_err < 1e-12, jac_err
    ok["thmV3_jacobian_vs_exhaustive_single_hit_err"] = jac_err

    # 6. THM-V4: Catalog integrity (33 viruses + 28 human onco/control loci)
    assert len(VIRAL_CATALOG_33) == 33 and len(ONCO_LOCI_28) == 28
    assert len({m["acc"] for m in VIRAL_CATALOG_33.values()}) == 33
    assert len({m["acc"] for m in ONCO_LOCI_28.values()}) == 28
    ok["thmV4_catalog_33vir_28onco_ok"] = True
    if verbose:
        print(json.dumps(_clean(ok), indent=2))
    return ok


def _layout_scatter_labels(xs, ys, labels, xlim, ylim, log_y=False, min_dx=0.075, min_dy=0.042, iters=120):
    """Deterministic Zero-RNG label-repulsion layout in normalized [0,1]^2 axes coordinates."""
    n = len(xs)
    if n == 0:
        return []
    x0, x1 = float(xlim[0]), float(xlim[1])
    y0, y1 = (math.log10(ylim[0]), math.log10(ylim[1])) if log_y else (float(ylim[0]), float(ylim[1]))
    px = np.array([(float(x) - x0) / max(x1 - x0, 1e-12) for x in xs], dtype=np.float64)
    py = np.array([((math.log10(float(y)) if log_y else float(y)) - y0) / max(y1 - y0, 1e-12) for y in ys], dtype=np.float64)
    lx = px + np.where(px > 0.80, -0.028, 0.014)
    ly = py + 0.010
    widths = np.array([max(min_dx, 0.0085 * len(str(st)) + 0.016) for st in labels], dtype=np.float64)
    for _ in range(iters):
        for i in range(n):
            for j in range(i + 1, n):
                dx = lx[j] - lx[i]
                dy = ly[j] - ly[i]
                req_x = 0.5 * (widths[i] + widths[j])
                req_y = min_dy
                if abs(dx) < req_x and abs(dy) < req_y:
                    ov_y = req_y - abs(dy)
                    ov_x = req_x - abs(dx)
                    if ov_y <= ov_x * 1.45:
                        sgn = 1.0 if (dy > 0 or (dy == 0 and py[j] >= py[i])) else -1.0
                        sh = 0.52 * ov_y * sgn
                        ly[i] -= sh
                        ly[j] += sh
                    else:
                        sgn = 1.0 if (dx > 0 or (dx == 0 and px[j] >= px[i])) else -1.0
                        sh = 0.35 * ov_x * sgn
                        lx[i] -= sh
                        lx[j] += sh
        lx = 0.96 * lx + 0.04 * (px + np.where(px > 0.80, -0.028, 0.014))
        ly = 0.96 * ly + 0.04 * (py + 0.008)
        lx = np.clip(lx, 0.015, 0.985 - widths * 0.65)
        ly = np.clip(ly, 0.025, 0.965)
    out = []
    for i in range(n):
        rx = x0 + lx[i] * (x1 - x0)
        ry = (10.0 ** (y0 + ly[i] * (y1 - y0))) if log_y else (y0 + ly[i] * (y1 - y0))
        out.append((rx, ry, math.hypot(lx[i] - px[i], ly[i] - py[i]) > 0.030))
    return out


def make_virome_onco_figures(payload, out_dir, dpi=300):
    if not HAVE_MPL:
        return []
    fd = Path(out_dir) / "figures"
    fd.mkdir(parents=True, exist_ok=True)
    saved = []
    vr = payload.get("viruses") or []
    oc = payload.get("onco_loci") or []

    # FIG V1: 33-Virus Toric Landscape (Impedance vs Human Euchromatin, PR2 Compliance & 3600-D Syntax)
    if vr:
        fig, axes = plt.subplots(1, 2, figsize=(18.5, 8.2))
        fig.suptitle(f"Theorem V1: 33-Virus Toric Impedance & Chargaff PR2 Landscape vs. Human T2T Euchromatin (v{VERSION}, Zero-RNG)",
                     fontsize=13.5, fontweight="bold", y=0.98)
        cls_col = {
            "retrovirus_int": "#d62728", "dna_onco_small": "#ff7f0e", "dna_benign_ctrl": "#2ca02c",
            "hepatitis_hcc": "#9467bd", "herpes_large": "#1f77b4", "vector_dna": "#17becf",
            "rna_acute_ctrl": "#7f7f7f", "phage_outgroup": "#8c564b",
        }
        vr_sorted = sorted(vr, key=lambda r: _fval(r.get("cos_Wmp_vs_human_grand"), -2.0))
        ys = np.arange(len(vr_sorted))
        lbls = [f"{r['label']} ({r['vclass'][:6]})" for r in vr_sorted]
        cos_h = [_fval(r.get("cos_Wmp_vs_human_grand"), 0.0) for r in vr_sorted]
        cols = [cls_col.get(r["vclass"], "#333333") for r in vr_sorted]
        ax = axes[0]
        ax.barh(ys, cos_h, color=cols, alpha=0.88)
        ax.axvline(0.0, color="k", lw=0.9)
        ax.set_yticks(ys); ax.set_yticklabels(lbls, fontsize=7.5)
        ax.set_xlabel("Toric Arrow Alignment cos( W_-+^virus , W_-+^human_euch )")
        ax.set_title("A. Viral-to-Human Directional Arrow Impedance (1024-D W_-+)")
        ax.grid(True, axis="x", ls=":", alpha=0.35)

        ax = axes[1]
        pts_v1 = []
        for r in vr:
            pr2 = 100.0 * _fval(r["klein_d3"].get("time_PR2_compliant_share_unb", r["klein_d3"].get("time_PR2_compliant_share")), np.nan)
            a11_r = r.get("arrow_D11") or {}
            ne = 100.0 * _fval(a11_r.get("pure_nonenv_arrow_share_unb", a11_r.get("pure_nonenv_arrow_share")), np.nan)
            if _finite(pr2) and _finite(ne):
                ax.scatter(pr2, ne, s=55, color=cls_col.get(r["vclass"], "#333"), alpha=0.88, edgecolors="k", lw=0.4, zorder=4)
                pts_v1.append((pr2, ne, r["label"]))
        v1_xlim, v1_ylim = (17.0, 114.0), (-3.5, 63.5)
        ax.set_xlim(*v1_xlim); ax.set_ylim(*v1_ylim)
        pos_v1 = _layout_scatter_labels([p[0] for p in pts_v1], [p[1] for p in pts_v1], [p[2] for p in pts_v1],
                                        v1_xlim, v1_ylim, log_y=False, min_dx=0.078, min_dy=0.042, iters=220)
        for (px_, py_, lb_), (lx_, ly_, draw_ln) in zip(pts_v1, pos_v1):
            if draw_ln:
                ax.plot([px_, lx_], [py_, ly_], color="#888888", lw=0.55, alpha=0.65, zorder=2)
            ax.text(lx_, ly_, lb_, fontsize=6.5, va="center", zorder=5)
        hb = payload.get("human_reference_baselines") or {}
        h_pr2 = 100.0 * _fval(hb.get("grand_PR2_ok_unb"), 0.999995)
        h_ne = 100.0 * _fval(hb.get("grand_nonenv_3600_unb", hb.get("grand_nonenv_3600_raw")), 0.370)
        ax.axvline(h_pr2, color="#1f77b4", ls="--", lw=1.1, label=f"Human Euchromatin PR2 ({h_pr2:.4f}%)")
        ax.axhline(h_ne, color="#17becf", ls=":", lw=1.1, label=f"Human Euchromatin NonEnv 3600-D unb ({h_ne:.1f}%)")
        ax.set_xlabel("Unbiased PR2-Compliant Directed Flux (%)")
        ax.set_ylabel("THM-V16 Unbiased 3600-D Non-Envelope Share (894/1024 in W_-+, %)")
        ax.set_title("B. Viral Chiral Parity Break vs. Unbiased Higher-Order 3600-D Syntax")
        ax.legend(fontsize=8, loc="upper left", framealpha=0.92); ax.grid(True, ls=":", alpha=0.35)

        fig.tight_layout(rect=[0, 0, 1, 0.945])
        p1 = fd / "fig_V1_virome_landscape.png"
        fig.savefig(p1, dpi=dpi, bbox_inches="tight"); plt.close(fig); saved.append(str(p1))
        print(f"  [FIG] {p1}")

    # FIG V2: 28 Human Cancer Loci — Mutagenesis Chiral Injection & Native Arrow Anatomy
    if oc:
        fig, axes = plt.subplots(1, 2, figsize=(18.5, 7.6))
        fig.suptitle(f"Theorems V3 & V4: Exact Zero-RNG COSMIC Signature Chiral Injection Across {len(oc)} Human Cancer & Control Loci (v{VERSION})",
                     fontsize=13.5, fontweight="bold", y=0.98)
        sigs = list(COSMIC_SIG_TABLES.keys())
        mat_ch = np.zeros((len(oc), len(sigs)), dtype=np.float64)
        for i, r in enumerate(oc):
            mj = r.get("mutagenesis_V3") or {}
            for j, sn in enumerate(sigs):
                mat_ch[i, j] = 100.0 * _fval((mj.get(sn) or {}).get("J_chiral_share_of_directed_unb"), 0.0)
        ax = axes[0]
        im = ax.imshow(mat_ch, aspect="auto", cmap="inferno", vmin=0.0, vmax=max(50.0, float(np.max(mat_ch))))
        ax.set_xticks(range(len(sigs))); ax.set_xticklabels(sigs, rotation=35, ha="right", fontsize=8)
        ax.set_yticks(range(len(oc))); ax.set_yticklabels([f"{r['label']} ({r['chr']})" for r in oc], fontsize=7.5)
        ax.set_title("A. Unbiased Jacobian Chiral Share (Zero 1/N Length Bias, %)")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)

        ax = axes[1]
        x_o = np.arange(len(oc))
        ap_dn = [max(_fval(((r.get("mutagenesis_V3") or {}).get("SBS2_APOBEC_plus") or {}).get("chi_mm_denovo_unb_x1e4"), 0.0), 1e-4) for r in oc]
        s1_chi = [max(_fval(((r.get("mutagenesis_V3") or {}).get("SBS1_CpG_sym") or {}).get("chi_mm_unb_x1e4"), 0.0), 1e-4) for r in oc]
        ax.bar(x_o - 0.19, ap_dn, 0.38, color="#d62728", alpha=0.88, label="APOBEC3 SBS2(+) De Novo χ_--^unb (×10⁴)")
        ax.bar(x_o + 0.19, s1_chi, 0.38, color="#1f77b4", alpha=0.88, label="SBS1 CpG (sym) Substrate-Coupled χ_--^unb (de novo ≡ 0)")
        ax.set_yscale("log"); ax.set_ylim(1.2e-4, 3.8); ax.set_xticks(x_o); ax.set_xticklabels([r["label"] for r in oc], rotation=50, ha="right", fontsize=7.5)
        ax.set_ylabel("U-Statistic Unbiased Chiral Susceptibility χ_--^unb × 10⁴ (log)")
        ax.set_title("B. Exact U-Statistic Unbiased Chiral Susceptibility (Zero 1/N Length Bias)")
        ax.legend(fontsize=8, loc="upper right", framealpha=0.9); ax.grid(True, axis="y", ls=":", alpha=0.35)

        fig.tight_layout(rect=[0, 0, 1, 0.945])
        p2 = fd / "fig_V2_onco_jacobian_loci.png"
        fig.savefig(p2, dpi=dpi, bbox_inches="tight"); plt.close(fig); saved.append(str(p2))
        print(f"  [FIG] {p2}")
    return saved


def write_virome_onco_markdown(payload, path):
    vr = payload.get("viruses") or []
    oc = payload.get("onco_loci") or []
    ins_rows = payload.get("proviral_insertions_V2") or []
    ref_info = payload.get("human_reference_source", "pooled_human_loci")
    st_dict = payload.get("self_test") or {}
    n_st = len([k for k in st_dict if k not in ("all_passed", "pilot_catalog_count")])
    L = [
        f"# Zero-RNG Toric Virome & Oncogenomics Audit (v{VERSION})", "",
        f"> **Methodology:** 100% deterministic, zero-RNG Klein $V_4$ & $V_4$-closed $496 \\oplus 3600$ screen engine (894/1024 dimensions inside $W_{{-+}}$) · {n_st}/{n_st} startup self-tests passed · Human Reference Arrow: `{ref_info}`.", "",
        "## 1. Theorem V1 — 33 Real NCBI RefSeq Viruses on the Codon Torus", "",
        "| Virus | Class | Mode | Len (bp) | GC % | ‖W_-+‖²_unb | z_iid(W_-+) | PR2-ok (unb) % | 3600-D unb % (raw) | CpG-free % | cos_raw | cos_unb (z≥3) | Best Human Chr (cos) | Top-1 Motif |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|",
    ]
    for r in vr:
        kv = r["klein_d3"]
        a11 = r.get("arrow_D11") or {}
        cu = r.get("cos_Wmp_unb_V18") or {}
        pr2 = kv.get("time_PR2_compliant_share_unb")
        if pr2 is None:
            pr2 = kv.get("time_PR2_compliant_share")
        pr2_s = f"{fnum(100.0*_fval(pr2, np.nan),2)}%" + ("*" if kv.get("W_mm_unb_censored") else "")
        ne_u = a11.get("pure_nonenv_arrow_share_unb", a11.get("pure_nonenv_arrow_share"))
        ne_r = a11.get("pure_nonenv_arrow_share")
        cu_s = fnum(cu.get("cos_unb"), 3, True) if cu.get("status") == "resolved" else "unresolved"
        t1 = (a11.get("top_motifs_all") or [None])[0]
        t1s = f"`{t1['motif_6mer']}` ({t1['flux_share_pct']:.1f}%)" if t1 else "n/a"
        bc = f"{r.get('best_human_chr','n/a')} ({fnum(r.get('best_human_chr_cos'),3,True)})" if r.get("best_human_chr") else "n/a"
        L.append(
            f"| **{r['label']}** | {r['vclass']} | {r['mode']} | {r['length_bp']:,} | {100*r['gc']:.1f}% | "
            f"{_fsci(kv.get('norm2_W_mp_unb'))} | {fnum(kv.get('z_iid_W_mp'),1)} | "
            f"{pr2_s} | {fnum(100.0*_fval(ne_u, np.nan),1)}% ({fnum(100.0*_fval(ne_r, np.nan),1)}%) | "
            f"{fnum(100.0*_fval(a11.get('cpg_free_arrow_share'), np.nan),1)}% | "
            f"{fnum(r.get('cos_Wmp_vs_human_grand'),3,True)} | {cu_s} | {bc} | {t1s} |"
        )

    if oc:
        L += [
            "", "## 2. Theorem V4 — 28 Real Human Cancer Driver, Fragile-Site & Control Loci (NCBI RefSeqGene)", "",
            "| Locus | Chr | Role | Len (bp) | GC % | ‖W_-+‖²_unb | PR2-ok (unb) % | Non-Env 3600-D % | CpG-free % | cos(W_-+, Human) | Pal d=6 Share % |",
            "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for r in oc:
            kv = r["klein_d3"]
            a11 = r.get("arrow_D11") or {}
            pr2 = kv.get("time_PR2_compliant_share_unb")
            if pr2 is None:
                pr2 = kv.get("time_PR2_compliant_share")
            pal6 = ((r.get("palindrome_spectrum") or {}).get("6") or {}).get("pal_Wmp_energy_share")
            L.append(
                f"| **{r['label']}** | {r['chr']} | {r['role']} | {r['length_bp']:,} | {100*r['gc']:.1f}% | "
                f"{_fsci(kv.get('norm2_W_mp_unb'))} | {fnum(100.0*_fval(pr2, np.nan),2)}% | "
                f"{fnum(100.0*_fval(a11.get('pure_nonenv_arrow_share_unb', a11.get('pure_nonenv_arrow_share')), np.nan),1)}% | "
                f"{fnum(100.0*_fval(a11.get('cpg_free_arrow_share'), np.nan),1)}% | "
                f"{fnum(r.get('cos_Wmp_vs_human_grand'),3,True)} | {fnum(100.0*_fval(pal6, np.nan),2)}% |"
            )

        L += [
            "", "## 3. Theorems V3 & V3.3 — Unbiased COSMIC Mutagenesis Jacobian & Exact De Novo vs. Substrate-Coupled Split", "",
            "| Locus | SBS1 Unb Chiral % | SBS1 De Novo χ_-- | SBS1 Sub χ_-- (×10⁴) | APOBEC-SBS2(+) Unb Chiral % | APOBEC De Novo χ_--^unb (×10⁴) | APOBEC De Novo % | ε_crit (1/4N) | UV-SBS7a(+) Unb Chiral % | Tobacco-SBS4(+) Unb Chiral % |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for r in oc:
            mj = r.get("mutagenesis_V3") or {}
            s1 = mj.get("SBS1_CpG_sym") or {}
            ap = mj.get("SBS2_APOBEC_plus") or {}
            uv = mj.get("SBS7a_UV_plus") or {}
            tb = mj.get("SBS4_Tobacco_plus") or {}
            L.append(
                f"| **{r['label']}** | {fnum(100.0*_fval(s1.get('J_chiral_share_of_directed_unb'), np.nan),2)}% | "
                f"{_fsci(s1.get('chi_mm_denovo_x1e4'))} | "
                f"{fnum(s1.get('chi_mm_unb_x1e4'),4)} | "
                f"{fnum(100.0*_fval(ap.get('J_chiral_share_of_directed_unb'), np.nan),2)}% | "
                f"{fnum(ap.get('chi_mm_denovo_unb_x1e4'),3)} | "
                f"{fnum(100.0*_fval(ap.get('denovo_chiral_share'), np.nan),1)}% | "
                f"{fnum(ap.get('eps_crit_4N'),3)} | "
                f"{fnum(100.0*_fval(uv.get('J_chiral_share_of_directed_unb'), np.nan),2)}% | "
                f"{fnum(100.0*_fval(tb.get('J_chiral_share_of_directed_unb'), np.nan),2)}% |"
            )

    if ins_rows:
        L += [
            "", "## 4. Theorem V2 — Deterministic In-Silico Concatenation Proviral Model into Real Human Onco-Loci", "",
            "| Virus -> Host Locus | Virus % | Host PR2-ok % | +Strand PR2-ok % | +Strand E/floor | HH E/floor | cos(ΔW_-- raw) | cos(ΔB_-- count-space) | cos(ΔB_-+ count-space) |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for row in ins_rows:
            L.append(
                f"| `{row['virus']}` -> `{row['host']}` | {100*row['virus_fraction']:.1f}% | "
                f"{fnum(100.0*_fval(row.get('host_PR2_ok_unb'), np.nan),2)}% | "
                f"{fnum(100.0*_fval(row.get('plus_PR2_ok_unb'), np.nan),2)}% | "
                f"{fnum(row.get('plus_W_mm_snr_over_4N'),1)}x | {fnum(row.get('head2head_W_mm_snr_over_4N'),1)}x | "
                f"{fnum(row.get('cos_dWmm_plus_vs_minus'),3,True)} | "
                f"{fnum(row.get('cos_dBmm_plus_vs_minus'),4,True)} | {fnum(row.get('cos_dBmp_plus_vs_minus'),4,True)} |"
            )
    L.append("")
    atomic_write(path, "\n".join(L) + "\n")




# =====================================================================
# 11. PAPER 5 EXPANSION (v2.1.0): T2T MICRO-NEIGHBORHOODS, 10-LAG & 1-KB TILES
#     [THM-V5] Exact T2T-CHM13v2.0 +/-50 kb Chromosomal Locus Micro-Scanner
#     [THM-V6] 10-Lag (d = 3..171 bp) Viral & Cancer-Locus Spectrum
#     [THM-V7] 1-kb Tiled COSMIC Mutagenesis Vulnerability Scanner
# =====================================================================
V21_LAGS = (3, 6, 9, 10, 11, 15, 30, 90, 147, 171)


def _load_human_euch_multilag_arrows(atlas_npz_path=None):
    """Loads grand-pooled human euchromatin W_-+ tensors across all 10 lags strictly from SCRIPT_DIR."""
    candidates = []
    if atlas_npz_path:
        candidates.append(Path(atlas_npz_path))
    for rel in (
        "human_t2t_euch_ref_v2_8.npz",
        "human_t2t_euch_ref_v2_6.npz",
        "atlas_b3_cache.npz",
        "results_human_t2t_atlas_v1_5_3/atlas_b3_cache.npz",
    ):
        candidates.append((SCRIPT_DIR / rel).resolve())
    for cand in candidates:
        if cand.exists():
            npz = np.load(cand)
            out = {}
            for d in V21_LAGS:
                B_sum = np.zeros((NC, NC), dtype=np.float64)
                for k in HUMAN_T2T_24:
                    key_l = f"{k}__euch_lag_{int(d)}" if int(d) != 3 else f"{k}__euch"
                    if key_l in npz.files:
                        B_sum += np.array(npz[key_l], dtype=np.float64)
                if float(B_sum.sum()) > 0.0:
                    kv = klein_v4_pr2_audit(B_sum)
                    if kv is not None:
                        out[int(d)] = kv["_W_mp"]
            npz.close()
            if out:
                return out
    return {}


def compute_multilag_profile_v6(nt, human_euch_lags=None, lags=V21_LAGS):
    """[THM-V6] Exact 10-lag (d=3..171 bp) toric spectrum on a viral or cancer-locus sequence."""
    L = int(nt.size)
    rows = []
    for d in lags:
        ld = int(d)
        if L <= ld + 30:
            continue
        B = count_triplet_pairs_lag(nt, d=ld)
        kv = klein_v4_pr2_audit(B)
        if kv is None or kv["N"] < 50:
            continue
        sc = two_sided_screen_full(B)
        a11 = thm_D11_arrow_screen_and_cpg_audit(B)
        Wmp = kv["_W_mp"]
        tot = float(np.sum(Wmp * Wmp))
        pal_e = float(np.sum(Wmp[np.arange(NC), MIRROR] ** 2))
        cos_h = None
        if human_euch_lags and ld in human_euch_lags:
            cos_h = _cos_mat_frob(Wmp, human_euch_lags[ld])
        rows.append({
            "lag_bp": ld,
            "N": kv["N"],
            "W_mp_unb": kv.get("norm2_W_mp_unb"),
            "W_mm_unb": kv.get("norm2_W_mm_unb"),
            "PR2_ok_unb": kv.get("time_PR2_compliant_share_unb"),
            "I_nonbr_MM": sc["I_nonbridge_MM"] if sc else None,
            "pure_nonenv_3600": a11["pure_nonenv_arrow_share"] if a11 else None,
            "pure_nonenv_3600_unb": a11["pure_nonenv_arrow_share_unb"] if a11 else None,
            "isotropy_index_raw": (
                a11["pure_nonenv_arrow_share"] / (447.0 / 512.0)
                if (a11 and a11.get("pure_nonenv_arrow_share") is not None) else None
            ),
            "cpg_free_share": a11["cpg_free_arrow_share"] if a11 else None,
            "pal_Wmp_share": (pal_e / tot) if tot > 1e-30 else None,
            "cos_vs_human_euch_lag": cos_h,
        })
    return rows


def scan_locus_vulnerability_tiles_v7(nt, tile_bp=1000):
    """[THM-V7] Non-overlapping 1-kb tile scan of COSMIC signature chiral injection across a locus."""
    L = int(nt.size)
    n_tiles = L // int(tile_bp)
    if n_tiles < 3:
        return None
    sigs = ("SBS1_CpG_sym", "SBS2_APOBEC_plus", "SBS13_APOBEC_plus", "SBS7a_UV_plus", "SBS4_Tobacco_plus")
    sig_sub = {sn: COSMIC_SIG_TABLES[sn] for sn in sigs}
    tiles = []
    for ti in range(n_tiles):
        s0 = ti * int(tile_bp)
        s1 = s0 + int(tile_bp)
        seg = nt[s0:s1]
        if int(np.count_nonzero(seg < 4)) < int(0.8 * tile_bp):
            continue
        mj = analyze_mutagenesis_signatures(seg, sig_tables=sig_sub, eps_eval=0.05, compute_denovo=False)
        row = {"tile_idx": ti, "start_kb": s0 / 1e3, "end_kb": s1 / 1e3}
        for sn in sigs:
            d = mj.get(sn) or {}
            sh_u = d.get("J_chiral_share_of_directed_unb")
            c_mm = d.get("chi_mm_unb_x1e4")
            c_mp = d.get("chi_mp_unb_x1e4")
            n_elig = int(round(_fval(d.get("eligible_site_frac"), 0.0) * _fval(d.get("N_interior"), 0.0)))
            u_clip = bool((c_mm is None or float(c_mm) <= 0.0) or (c_mp is None or float(c_mp) <= 0.0))
            deg_half = bool(sh_u is not None and abs(float(sh_u) - 0.5) < 1e-12)
            row[sn + "_chiral_share"] = sh_u
            row[sn + "_chi_mm_x1e4"] = c_mm
            row[sn + "_cos_Wmp"] = d.get("cos_Wmp_native_vs_J")
            row[sn + "_n_elig_sites"] = n_elig
            row[sn + "_u_clipped"] = u_clip
            row[sn + "_degenerate_exact_half"] = deg_half
        tiles.append(row)
    if not tiles:
        return None
    summary = {}
    for sn in sigs:
        vals = [float(t[sn + "_chi_mm_x1e4"]) for t in tiles if _finite(t.get(sn + "_chi_mm_x1e4"))]
        cvals = [float(t[sn + "_chiral_share"]) for t in tiles if _finite(t.get(sn + "_chiral_share"))]
        cvals_clean = [
            float(t[sn + "_chiral_share"]) for t in tiles
            if _finite(t.get(sn + "_chiral_share"))
            and (not t.get(sn + "_u_clipped"))
            and (not t.get(sn + "_degenerate_exact_half"))
        ]
        n_deg = sum(1 for t in tiles if t.get(sn + "_u_clipped") or t.get(sn + "_degenerate_exact_half"))
        if vals:
            bi = int(np.argmax(vals))
            summary[sn] = {
                "median_chi_mm_x1e4": float(np.median(vals)),
                "max_chi_mm_x1e4": float(vals[bi]),
                "peak_tile_kb": f"{tiles[bi]['start_kb']:.0f}-{tiles[bi]['end_kb']:.0f} kb",
                "median_chiral_share": float(np.median(cvals)) if cvals else None,
                "median_chiral_share_nondegenerate": float(np.median(cvals_clean)) if cvals_clean else None,
                "n_degenerate_or_clipped_tiles": int(n_deg),
            }
    return {"tile_bp": int(tile_bp), "n_tiles": len(tiles), "summary": summary, "tiles": tiles}


def locate_anchor_on_chromosome(nt_chr, nt_locus, anchor_len=60):
    """[THM-V5 v2.6] SNP-robust deterministic multi-scale (60 -> 36 -> 24 bp) anchor locator
    recovering exact 0-based start coordinate & strand of nt_locus inside nt_chr."""
    L_loc = int(nt_locus.size)
    L_chr = int(nt_chr.size)
    if L_loc < 60 or L_chr < 60:
        return None, None
    chr_bytes = nt_chr.tobytes()

    def _verify_candidate(start_est, strand):
        # Verify with 7 short 18-bp probes across the locus
        hits = 0
        checks = 0
        for frac in (0.10, 0.25, 0.40, 0.50, 0.65, 0.80, 0.90):
            off = int(frac * max(L_loc - 20, 1))
            probe = nt_locus[slice(off, off + 18)]
            if probe.size < 18 or int(probe.max()) >= 4:
                continue
            checks += 1
            if strand == "+":
                p0 = start_est + off
                p_ref = probe
            else:
                p0 = start_est + (L_loc - off - 18)
                p_ref = reverse_complement_nt(probe)
            lo = max(0, p0 - 4000)
            hi = min(L_chr, p0 + 4000)
            if hi > lo and chr_bytes.find(p_ref.tobytes(), lo, hi) >= 0:
                hits += 1
        return hits >= max(2, checks // 2)

    lens_try = [int(anchor_len)]
    for al in (60, 36, 24):
        if al not in lens_try and al < L_loc:
            lens_try.append(al)

    for alen in lens_try:
        max_off = L_loc - alen
        if max_off < 0:
            continue
        offsets = sorted(set(
            [0, min(200, max_off), min(1000, max_off)]
            + [int(max_off * j / 24) for j in range(1, 24)]
        ))
        for off in offsets:
            cand = nt_locus[slice(off, off + alen)]
            if cand.size < alen or int(cand.max()) >= 4:
                continue
            cb = cand.tobytes()
            pos = chr_bytes.find(cb)
            while pos >= 0:
                st_est = int(pos - off)
                if _verify_candidate(st_est, "+"):
                    return st_est, "+"
                pos = chr_bytes.find(cb, pos + 1)
            rcb = reverse_complement_nt(cand).tobytes()
            pos_rc = chr_bytes.find(rcb)
            while pos_rc >= 0:
                st_est = int(pos_rc - (L_loc - off - alen))
                if _verify_candidate(st_est, "-"):
                    return st_est, "-"
                pos_rc = chr_bytes.find(rcb, pos_rc + 1)
    return None, None


def audit_t2t_neighborhood_v5(chr_dir, onco_nt_cache, viral_nt_cache, human_ref, half_win_bp=50000, sub_bp=5000, auto_download=False):
    """[THM-V5] Audits exact +/-50 kb T2T-CHM13v2.0 chromosomal neighborhoods around
    all 28 human loci directly on NC_0609xx.1 FASTA files inside SCRIPT_DIR/chr_cache."""
    chr_dir = Path(chr_dir)
    if auto_download:
        chr_dir.mkdir(parents=True, exist_ok=True)
    if not chr_dir.exists():
        return []
    target_loci = list(ONCO_LOCI_28.keys())
    by_chr = {}
    for lname in target_loci:
        if lname in onco_nt_cache and lname in ONCO_LOCI_28:
            ch = ONCO_LOCI_28[lname]["chr"]
            by_chr.setdefault(ch, []).append(lname)
    if not by_chr:
        return []

    w_h_grand = human_ref.get("grand_W_mp") if human_ref else None
    w_h_unit = (w_h_grand / max(float(np.linalg.norm(w_h_grand, ord="fro")), TINY)) if w_h_grand is not None else None
    viral_Wmp = {}
    viral_Wmp_resH = {}
    for vn, vnt in viral_nt_cache.items():
        kv_v = klein_v4_pr2_audit(count_triplet_pairs_lag(vnt, d=3))
        if kv_v is not None:
            wv = kv_v["_W_mp"]
            viral_Wmp[vn] = wv
            if w_h_unit is not None:
                viral_Wmp_resH[vn] = wv - float(np.sum(wv * w_h_unit)) * w_h_unit

    out_rows = []
    for ch, lnames in by_chr.items():
        ckey = f"Human_T2T_{ch}"
        if ckey not in HUMAN_T2T_24:
            continue
        acc = HUMAN_T2T_24[ckey]["acc"]
        fp = chr_dir / f"{ckey}_{acc}.fna.gz"
        if not fp.exists():
            hits = list(chr_dir.glob(f"*{acc}*.fna*"))
            if hits:
                fp = hits[0]
        if not fp.exists() or fp.stat().st_size <= 150:
            if auto_download:
                fp = download_chromosome_fasta(ckey, HUMAN_T2T_24[ckey], chr_dir)
            else:
                continue
        print(f"  [THM-V5 T2T] Loading {ch} ({fp.name}) for exact neighborhood micro-scan of {','.join(lnames)}...")
        nt_chr, _ = load_chromosome_nt(fp)
        L_chr = int(nt_chr.size)
        chr_base_wmp = (human_ref.get("chr_W_mp") or {}).get(ch) if human_ref else None

        for lname in lnames:
            nt_loc = onco_nt_cache[lname]
            start_bp, strand = locate_anchor_on_chromosome(nt_chr, nt_loc, anchor_len=60)
            if start_bp is None:
                print(f"    [WARN THM-V5] Could not locate {lname} on {ch}")
                continue
            center_bp = max(0, min(L_chr - 1, start_bp + int(nt_loc.size) // 2))
            w0 = max(0, center_bp - int(half_win_bp))
            w1 = min(L_chr, center_bp + int(half_win_bp))
            seg_win = nt_chr[w0:w1]
            B3_w = count_triplet_pairs_lag(seg_win, d=3)
            kv_w = klein_v4_pr2_audit(B3_w)
            a11_w = thm_D11_arrow_screen_and_cpg_audit(B3_w)
            B6_w = count_triplet_pairs_lag(seg_win, d=6)
            kv_6 = klein_v4_pr2_audit(B6_w)
            pal6_e = float(np.sum(kv_6["_W_mp"][np.arange(NC), MIRROR] ** 2)) / max(float(np.sum(kv_6["_W_mp"] ** 2)), TINY)
            cnt4 = np.bincount(seg_win[seg_win < 4], minlength=4).astype(np.float64)
            gc_w = float((cnt4[1] + cnt4[2]) / max(float(cnt4.sum()), 1.0))

            n_sub = (w1 - w0) // int(sub_bp)
            sub_nonenv, sub_wmp = [], []
            for si in range(n_sub):
                sb = seg_win[si * sub_bp:(si + 1) * sub_bp]
                Bs = count_triplet_pairs_lag(sb, d=3)
                kvs = klein_v4_pr2_audit(Bs)
                a11s = thm_D11_arrow_screen_and_cpg_audit(Bs)
                if kvs and a11s:
                    sub_wmp.append(kvs.get("norm2_W_mp_unb"))
                    sub_nonenv.append(a11s.get("pure_nonenv_arrow_share"))

            v_cos = {vn: _cos_mat_frob(vwmp, kv_w["_W_mp"]) for vn, vwmp in viral_Wmp.items()}
            best_v = max(v_cos.items(), key=lambda kv_i: _fval(kv_i[1], -2.0)) if v_cos else (None, None)
            best_res = (None, None)
            cos_h_win = None
            if w_h_unit is not None:
                cos_h_win = _cos_mat_frob(kv_w["_W_mp"], w_h_grand)
                wl_resH = kv_w["_W_mp"] - float(np.sum(kv_w["_W_mp"] * w_h_unit)) * w_h_unit
                v_cos_res = {vn: _cos_mat_frob(vw_r, wl_resH) for vn, vw_r in viral_Wmp_resH.items()}
                if v_cos_res:
                    best_res = max(v_cos_res.items(), key=lambda kv_i: _fval(kv_i[1], -2.0))
            out_rows.append({
                "locus": lname,
                "chr": ch,
                "role": ONCO_LOCI_28[lname]["role"],
                "t2t_center_mb": center_bp / 1e6,
                "locus_strand": strand,
                "window_bp": int(w1 - w0),
                "window_gc": gc_w,
                "W_mp_unb": kv_w.get("norm2_W_mp_unb"),
                "PR2_ok_unb": kv_w.get("time_PR2_compliant_share_unb"),
                "pure_nonenv_3600": a11_w.get("pure_nonenv_arrow_share"),
                "pure_nonenv_3600_unb": a11_w.get("pure_nonenv_arrow_share_unb"),
                "cpg_free_share": a11_w.get("cpg_free_arrow_share"),
                "pal_d6_Wmp_share": pal6_e,
                "cos_vs_chr_euch": _cos_mat_frob(kv_w["_W_mp"], chr_base_wmp) if chr_base_wmp is not None else None,
                "cos_vs_human_grand": cos_h_win,
                "viral_cos_Wmp": v_cos,
                "top_matching_virus": best_v[0],
                "top_matching_virus_cos": best_v[1],
                "top_residual_virus": best_res[0],
                "top_residual_virus_cos": best_res[1],
                "sub5kb_nonenv_3600": sub_nonenv,
                "sub5kb_W_mp_unb": sub_wmp,
            })
            print(f"    -> {lname:6s} @ {ch}:{center_bp/1e6:.2f} Mb ({strand}) | GC={100*gc_w:.1f}% | "
                  f"NonEnv3600={fnum(100*_fval(a11_w.get('pure_nonenv_arrow_share'),0),1)}% | "
                  f"Pal6={100*pal6_e:.1f}% | top_vir={best_v[0]} ({fnum(best_v[1],3,True)})")
        del nt_chr
    return out_rows


def self_test_virome_onco_v21(verbose=False):
    """[THM-V5..V7] 3 additional deterministic Zero-RNG self-tests for v2.1.0."""
    ok = {}
    L = 4000
    _ix = np.arange(L, dtype=np.uint64) + np.uint64(101)
    _zz = (_ix ^ (_ix >> np.uint64(30))) * np.uint64(0xbf58476d1ce4e5b9)
    _zz = (_zz ^ (_zz >> np.uint64(27))) * np.uint64(0x94d049bb133111eb)
    nt = ((_zz ^ (_zz >> np.uint64(31))) & np.uint64(3)).astype(np.uint8)

    # 1. Exact T2T anchor locator on both + and - strands
    anc = nt[600:900]
    p_fwd, st_f = locate_anchor_on_chromosome(nt, anc, anchor_len=60)
    p_rev, st_r = locate_anchor_on_chromosome(nt, reverse_complement_nt(anc), anchor_len=60)
    assert p_fwd == 600 and st_f == "+", (p_fwd, st_f)
    assert p_rev == 600 and st_r == "-", (p_rev, st_r)
    ok["thmV5_t2t_exact_anchor_locator_ok"] = True

    # 2. Multi-lag spectrum consistency at d=3 with direct klein_v4_pr2_audit
    ml = compute_multilag_profile_v6(nt, lags=(3, 6, 15))
    kv3 = klein_v4_pr2_audit(count_triplet_pairs_lag(nt, d=3))
    assert len(ml) == 3 and abs(ml[0]["W_mp_unb"] - kv3["norm2_W_mp_unb"]) < 1e-18
    ok["thmV6_multilag_profile_exact_ok"] = True

    # 3. 1-kb tile mutagenesis scanner exactness
    tl = scan_locus_vulnerability_tiles_v7(nt, tile_bp=1000)
    assert tl is not None and tl["n_tiles"] == 4 and "SBS2_APOBEC_plus" in tl["summary"]
    ok["thmV7_tile_vulnerability_scanner_ok"] = True
    ok.update(self_test_virome_onco_v22(verbose=False))
    if verbose:
        print(json.dumps(_clean(ok), indent=2))
    return ok


def append_v21_markdown_and_figures(payload, out_dir, md_path, dpi=300):
    t2t_rows = payload.get("t2t_neighborhoods_V5") or []
    vr = payload.get("viruses") or []
    oc = payload.get("onco_loci") or []

    # Append Sections 5, 6, 7 to Markdown report
    extra_md = []
    if t2t_rows:
        extra_md += [
            "", "## 5. Theorem V5 — Exact T2T-CHM13v2.0 ±50 kb Chromosomal Neighborhoods Around Cancer & Viral Hotspots", "",
            "| Locus | T2T Coord | Strand | Win GC % | ‖W_-+‖²_unb | PR2-ok % | Non-Env 3600-D % | Pal d=6 % | cos(vs Chr Euch) | Top Matching Virus (cos) |",
            "|---|---|---:|---:|---:|---:|---:|---:|---:|---|",
        ]
        for r in t2t_rows:
            extra_md.append(
                f"| **{r['locus']}** | `{r['chr']}:{r['t2t_center_mb']:.2f} Mb` | `{r['locus_strand']}` | "
                f"{100*r['window_gc']:.1f}% | {_fsci(r.get('W_mp_unb'))} | "
                f"{fnum(100.0*_fval(r.get('PR2_ok_unb'), np.nan),2)}% | "
                f"{fnum(100.0*_fval(r.get('pure_nonenv_3600'), np.nan),1)}% | "
                f"{fnum(100.0*_fval(r.get('pal_d6_Wmp_share'), np.nan),1)}% | "
                f"{fnum(r.get('cos_vs_chr_euch'),3,True)} | "
                f"`{r.get('top_matching_virus')}` ({fnum(r.get('top_matching_virus_cos'),3,True)}) |"
            )

    # Section 6: Multi-lag summary (d=3 vs d=6 vs d=30 vs d=171) for key viruses & loci
    ml_vr = [r for r in vr if r.get("multilag_V6")]
    if ml_vr:
        extra_md += [
            "", "## 6. Theorem V6 — Multi-Lag Toric Decay & Human Impedance Across Lags (d = 3, 6, 30, 171 bp)", "",
            "| Virus | Class | d=3 ‖W_-+‖² (cos_Hum) | d=6 ‖W_-+‖² (Pal %) | d=30 ‖W_-+‖² (cos_Hum) | d=171 ‖W_-+‖² (cos_Hum) |",
            "|---|---|---|---|---|---|",
        ]
        for r in ml_vr:
            by_d = {row["lag_bp"]: row for row in r["multilag_V6"]}
            def _fmt_d(ld, show_pal=False):
                rw = by_d.get(ld)
                if not rw:
                    return "n/a"
                extra = f"Pal {100*_fval(rw.get('pal_Wmp_share'),0):.1f}%" if show_pal else f"cos {fnum(rw.get('cos_vs_human_euch_lag'),2,True)}"
                return f"{_fsci(rw.get('W_mp_unb'))} ({extra})"
            extra_md.append(
                f"| **{r['label']}** | {r['vclass']} | {_fmt_d(3)} | {_fmt_d(6, True)} | {_fmt_d(30)} | {_fmt_d(171)} |"
            )

    # Section 7: 1-kb Tile Vulnerability Hotspots across Cancer Loci
    tl_oc = [r for r in oc if r.get("tile_vulnerability_V7")]
    if tl_oc:
        extra_md += [
            "", "## 7. Theorem V7 — 1-kb Tiled COSMIC Mutagenesis Chiral-Susceptibility Hotspots (χ_-- × 10⁴)", "",
            "| Locus | Tiles (1 kb) | APOBEC SBS2(+) Med / Peak χ_--×10⁴ (Hotspot) | APOBEC SBS13(+) Peak χ_--×10⁴ (Hotspot) | UV SBS7a(+) Peak χ_--×10⁴ | Tobacco SBS4(+) Peak χ_--×10⁴ |",
            "|---|---:|---|---|---:|---:|",
        ]
        for r in tl_oc:
            tv = r["tile_vulnerability_V7"]
            sm = tv["summary"]
            s2 = sm.get("SBS2_APOBEC_plus") or {}
            s13 = sm.get("SBS13_APOBEC_plus") or {}
            uv = sm.get("SBS7a_UV_plus") or {}
            tb = sm.get("SBS4_Tobacco_plus") or {}
            extra_md.append(
                f"| **{r['label']}** | {tv['n_tiles']} | "
                f"{fnum(s2.get('median_chi_mm_x1e4'),2)} / **{fnum(s2.get('max_chi_mm_x1e4'),2)}** (`{s2.get('peak_tile_kb','n/a')}`) | "
                f"**{fnum(s13.get('max_chi_mm_x1e4'),2)}** (`{s13.get('peak_tile_kb','n/a')}`) | "
                f"{fnum(uv.get('max_chi_mm_x1e4'),2)} | {fnum(tb.get('max_chi_mm_x1e4'),2)} |"
            )

    if extra_md:
        old_md = Path(md_path).read_text(encoding="utf-8")
        atomic_write(md_path, old_md.rstrip() + "\n" + "\n".join(extra_md) + "\n")

    if not HAVE_MPL:
        return []
    fd = Path(out_dir) / "figures"
    fd.mkdir(parents=True, exist_ok=True)
    saved = []

    # FIG V3: Multi-Lag Viral vs Human Euchromatin Decay + T2T +/-50 kb Neighborhood Profiles
    if ml_vr:
        fig, axes = plt.subplots(1, 2, figsize=(18.5, 6.8))
        fig.suptitle(f"Theorems V5 & V6: Multi-Lag Toric Decay (d = 3..171 bp) & T2T-CHM13v2.0 ±50 kb Locus Micro-Scans (v{VERSION})",
                     fontsize=13.5, fontweight="bold", y=0.98)
        ax = axes[0]
        show_v = ["HIV1_HXB2", "HTLV1", "HPV16", "HPV6b", "HBV", "EBV_HHV4", "SARS_CoV_2", "Phage_T4"]
        by_lbl_v = {r["label"]: r for r in ml_vr}
        for vn in show_v:
            rv = by_lbl_v.get(vn)
            if not rv:
                continue
            ds = [row["lag_bp"] for row in rv["multilag_V6"]]
            ws = [max(_fval(row.get("W_mp_unb"), 1e-8), 1e-8) for row in rv["multilag_V6"]]
            ax.plot(range(len(ds)), ws, marker="o", ms=4.2, lw=1.6, label=vn)
        ax.set_yscale("log"); ax.set_ylim(2.5e-9, 3.5e-5)
        ax.set_xticks(range(len(V21_LAGS))); ax.set_xticklabels([str(d) for d in V21_LAGS])
        ax.set_xlabel("Lag d (bp)"); ax.set_ylabel("Unbiased PR2-Compliant Flux ||W_-+||²_unb (log)")
        ax.set_title("A. Multi-Lag Toric Decay Across Integrating vs. Control Viruses")
        ax.legend(fontsize=7.5, ncol=2, loc="lower left", framealpha=0.9); ax.grid(True, ls=":", alpha=0.35)

        ax = axes[1]
        if t2t_rows:
            focal_v5 = {"MYC", "TERT", "FHIT", "TP53", "BRCA1", "PIK3CA", "ERBB2", "APOE"}
            for tr_row in t2t_rows:
                if tr_row["locus"] not in focal_v5:
                    continue
                ys = tr_row.get("sub5kb_nonenv_3600") or []
                if ys:
                    xs_kb = np.linspace(-47.5, 47.5, len(ys))
                    ax.plot(xs_kb, [100.0 * _fval(v, np.nan) for v in ys], marker=".", lw=1.5,
                            label=f"{tr_row['locus']} ({tr_row['chr']})")
            hb3 = payload.get("human_reference_baselines") or {}
            h_ne3 = 100.0 * _fval(hb3.get("grand_nonenv_3600_unb", hb3.get("grand_nonenv_3600_raw")), 0.370)
            ax.axvline(0.0, color="k", ls="--", lw=1.0, label="Locus Center")
            ax.axhline(h_ne3, color="#17becf", ls=":", lw=1.1, label=f"Genome Euchromatin Mean ({h_ne3:.1f}%)")
            ax.set_ylim(33.0, 93.0)
            ax.set_xlabel("Distance from Locus Center on T2T-CHM13v2.0 (kb)")
            ax.set_ylabel("5-kb Sector Raw Non-Envelope 3600-D Share (%)")
            ax.set_title("B. Local 3600-D Syntax Field Around Human Cancer / Viral Hotspots (±50 kb)")
            ax.legend(fontsize=7.5, ncol=2, loc="lower right", bbox_to_anchor=(0.99, 0.08), framealpha=0.92); ax.grid(True, ls=":", alpha=0.35)
        else:
            ax.text(0.5, 0.5, "Pass --chr-dir to render T2T +/-50 kb tracks", ha="center", va="center", transform=ax.transAxes)

        fig.tight_layout(rect=[0, 0, 1, 0.945])
        p3 = fd / "fig_V3_t2t_neighborhood_multilag.png"
        fig.savefig(p3, dpi=dpi, bbox_inches="tight"); plt.close(fig); saved.append(str(p3))
        print(f"  [FIG] {p3}")

    # FIG V4: 1-kb Tiled APOBEC3 Chiral Vulnerability Profiles Across Key Cancer Drivers
    focus_genes = ["TP53", "PIK3CA", "BRAF", "TERT", "MYC", "BRCA1"]
    by_lbl_o = {r["label"]: r for r in tl_oc}
    avail_fg = [g for g in focus_genes if g in by_lbl_o]
    if avail_fg:
        n_p = len(avail_fg)
        fig, axes = plt.subplots(2, 3, figsize=(18.5, 8.2))
        fig.suptitle(f"Theorem V7: 1-kb Tiled APOBEC3 Chiral Susceptibility (χ_-- × 10⁴) Across Key Cancer Loci (v{VERSION})",
                     fontsize=13.5, fontweight="bold", y=0.98)
        for idx_g, ax in enumerate(axes.ravel()):
            if idx_g >= n_p:
                ax.axis("off")
                continue
            gn = avail_fg[idx_g]
            tiles = by_lbl_o[gn]["tile_vulnerability_V7"]["tiles"]
            xk = [t["start_kb"] + 0.5 for t in tiles]
            y_s2 = [_fval(t.get("SBS2_APOBEC_plus_chi_mm_x1e4"), 0.0) for t in tiles]
            y_s13 = [_fval(t.get("SBS13_APOBEC_plus_chi_mm_x1e4"), 0.0) for t in tiles]
            y_s1 = [_fval(t.get("SBS1_CpG_sym_chi_mm_x1e4"), 0.0) for t in tiles]
            ax.plot(xk, y_s2, color="#d62728", lw=1.4, label="APOBEC SBS2 (+)")
            ax.plot(xk, y_s13, color="#ff7f0e", lw=1.2, alpha=0.85, label="APOBEC SBS13 (+)")
            ax.plot(xk, y_s1, color="#1f77b4", lw=1.0, ls="--", label="SBS1 CpG (sym; sub-coupled only, de novo≡0)")
            ax.set_title(f"{gn} ({by_lbl_o[gn]['chr']})", fontsize=10.5, fontweight="bold")
            ax.set_xlabel("Locus Coordinate (kb)"); ax.set_ylabel("Chiral Susceptibility χ_-- × 10⁴")
            if idx_g == 0:
                ax.set_ylim(top=max(max(y_s2 + y_s13 + y_s1) * 1.24, 0.5))
                ax.legend(fontsize=7.5, loc="upper right", framealpha=0.92)
            ax.grid(True, ls=":", alpha=0.35)
        fig.tight_layout(rect=[0, 0, 1, 0.945])
        p4 = fd / "fig_V4_locus_tile_vulnerability.png"
        fig.savefig(p4, dpi=dpi, bbox_inches="tight"); plt.close(fig); saved.append(str(p4))
        print(f"  [FIG] {p4}")
    return saved




# =====================================================================
# 12. PAPER 5 EXPANSION (v2.2.0): DUPLEX SYNONYMOUS BARRIER & GROUP TESTS
#     [THM-V8] 6-Frame Duplex Stop Barrier & Synonymous Mutation Disruption
#     [THM-V9] Exact Nonparametric Mann-Whitney Group Discrimination
# =====================================================================
_STOP_TRIPLETS_IDX = frozenset({
    TRIPLET_IDX["TAA"], TRIPLET_IDX["TAG"], TRIPLET_IDX["TGA"]
})
_ANTISENSE_STOP_TRIPLETS_IDX = frozenset({
    TRIPLET_IDX["TTA"], TRIPLET_IDX["CTA"], TRIPLET_IDX["TCA"]
})

_GENETIC_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
_CODON_AA_ARR = np.array([_GENETIC_CODE[c] for c in TRIPLETS], dtype="<U1")


def _extract_longest_orf_nt(nt, min_codons=60, max_codons=2500):
    """Finds the longest valid ATG -> STOP open reading frame across all 3 forward phases."""
    L = int(nt.size)
    best_span = None
    best_len = 0
    atg_idx = TRIPLET_IDX["ATG"]
    for ph in (0, 1, 2):
        n_cod = (L - ph) // 3
        if n_cod < min_codons:
            continue
        sub = nt[ph:ph + 3 * n_cod].reshape(n_cod, 3)
        ok_c = (sub[:, 0] < 4) & (sub[:, 1] < 4) & (sub[:, 2] < 4)
        c_ids = np.where(
            ok_c,
            (sub[:, 0].astype(np.int32) << 4) | (sub[:, 1].astype(np.int32) << 2) | sub[:, 2].astype(np.int32),
            -1,
        )
        cur_start = -1
        for ci in range(n_cod):
            cid = int(c_ids[ci])
            if cid < 0:
                cur_start = -1
                continue
            if cur_start < 0:
                if cid == atg_idx:
                    cur_start = ci
            else:
                if cid in _STOP_TRIPLETS_IDX:
                    span_cod = ci - cur_start
                    if span_cod >= min_codons and span_cod > best_len:
                        best_len = span_cod
                        use_cod = min(span_cod, max_codons)
                        best_span = nt[ph + 3 * cur_start:ph + 3 * (cur_start + use_cod)]
                    cur_start = -1
    return best_span, int(best_len)


def audit_synonymous_duplex_barrier_v8(onco_nt_cache, onco_results):
    """[THM-V8] Audits 6-frame stop-codon barrier and enumerates EVERY single-nucleotide
    synonymous wobble substitution in the longest coding ORF of each human locus,
    measuring antisense (-0) shadow-stop creation/destruction & 3600-D non-envelope shift."""
    rows = []
    for r in onco_results:
        name = r["label"]
        nt = onco_nt_cache.get(name)
        if nt is None or int(nt.size) < 300:
            continue
        # 6-frame stop density across full locus
        rc = reverse_complement_nt(nt)
        stop_dens = []
        for strand_nt in (nt, rc):
            for ph in (0, 1, 2):
                nc_f = (int(strand_nt.size) - ph) // 3
                if nc_f < 10:
                    continue
                sub = strand_nt[ph:ph + 3 * nc_f].reshape(nc_f, 3)
                ok = (sub[:, 0] < 4) & (sub[:, 1] < 4) & (sub[:, 2] < 4)
                cids = (sub[ok, 0].astype(np.int32) << 4) | (sub[ok, 1].astype(np.int32) << 2) | sub[ok, 2].astype(np.int32)
                if cids.size > 0:
                    n_st = sum(int(np.count_nonzero(cids == st)) for st in _STOP_TRIPLETS_IDX)
                    stop_dens.append(float(n_st) / float(cids.size))
        mean_6frame_stop = float(np.mean(stop_dens)) if stop_dens else None

        orf_nt, full_orf_codons = _extract_longest_orf_nt(nt, min_codons=60, max_codons=2000)
        if orf_nt is None or int(orf_nt.size) < 180:
            continue
        B0 = count_triplet_pairs_lag(orf_nt, d=3)
        N0 = float(B0.sum())
        if N0 <= 10.0:
            continue
        kv0 = klein_v4_pr2_audit(B0)
        Wmp0_ne = kv0["_W_mp"] - pi_env_project(kv0["_W_mp"])

        n_cod = int(orf_nt.size) // 3
        n_syn_total = 0
        n_anti_stop_flip = 0
        sum_dWmp_ne2 = 0.0
        sum_dWmm2 = 0.0

        # Enumerate all synonymous 3rd-position wobble transitions across strictly interior codons
        # Wobble base is at p0+2; all 6 affected d=3 6-mers span [p0-3, p0+7] inclusive (11 nt).
        for ci in range(1, n_cod - 2):
            p0 = 3 * ci
            b1, b2, b3 = int(orf_nt[p0]), int(orf_nt[p0 + 1]), int(orf_nt[p0 + 2])
            if b1 >= 4 or b2 >= 4 or b3 >= 4:
                continue
            cid_old = (b1 << 4) | (b2 << 2) | b3
            aa_old = _CODON_AA_ARR[cid_old]
            if aa_old == "*":
                continue
            was_anti = cid_old in _ANTISENSE_STOP_TRIPLETS_IDX
            lo = p0 - 3
            hi = p0 + 8
            sub_old = orf_nt[lo:hi]
            B_sub_old = count_triplet_pairs_lag(sub_old, d=3)
            for b3_new in range(4):
                if b3_new == b3:
                    continue
                cid_new = (b1 << 4) | (b2 << 2) | b3_new
                if _CODON_AA_ARR[cid_new] != aa_old:
                    continue
                n_syn_total += 1
                is_anti = cid_new in _ANTISENSE_STOP_TRIPLETS_IDX
                if was_anti != is_anti:
                    n_anti_stop_flip += 1
                sub_new = sub_old.copy()
                sub_new[p0 + 2 - lo] = b3_new
                dB = count_triplet_pairs_lag(sub_new, d=3) - B_sub_old
                dM = dB / N0
                TP = dM.T
                RP = dM.T[np.ix_(MIRROR, MIRROR)]
                SP = dM[np.ix_(MIRROR, MIRROR)]
                dWmp = 0.25 * (dM - TP + RP - SP)
                dWmm = 0.25 * (dM - TP - RP + SP)
                dWmp_ne = dWmp - pi_env_project(dWmp)
                sum_dWmp_ne2 += float(np.sum(dWmp_ne * dWmp_ne))
                sum_dWmm2 += float(np.sum(dWmm * dWmm))

        if n_syn_total > 0:
            rows.append({
                "locus": name,
                "chr": r["chr"],
                "role": r["role"],
                "mean_6frame_stop_freq": mean_6frame_stop,
                "longest_orf_codons": full_orf_codons,
                "audited_codons": n_cod,
                "n_synonymous_wobble_hits": n_syn_total,
                "antisense_m0_stop_flip_share": float(n_anti_stop_flip) / float(n_syn_total),
                "mean_dWmp_3600_norm2_per_syn": sum_dWmp_ne2 / float(n_syn_total),
                "mean_dBmp_3600_count2_per_syn": (N0 * N0) * sum_dWmp_ne2 / float(n_syn_total),
                "barrier_ratio_to_705_256": ((N0 * N0) * sum_dWmp_ne2 / float(n_syn_total)) / (705.0 / 256.0),
                "mean_dWmm_norm2_per_syn": sum_dWmm2 / float(n_syn_total),
                "native_Wmp_3600_norm2": float(np.sum(Wmp0_ne * Wmp0_ne)),
            })
    return rows


def exact_mann_whitney_u(x, y):
    """[THM-V9] Exact deterministic Zero-RNG Mann-Whitney U test with rank-biserial r_rb
    and tie-corrected two-sided p-value (exact DP when n1+n2 <= 36 and tie-free)."""
    xa = np.asarray([float(v) for v in x if _finite(v)], dtype=np.float64)
    ya = np.asarray([float(v) for v in y if _finite(v)], dtype=np.float64)
    n1, n2 = int(xa.size), int(ya.size)
    if n1 < 2 or n2 < 2:
        return None
    all_v = np.concatenate([xa, ya])
    order = np.argsort(all_v, kind="mergesort")
    sv = all_v[order]
    n = n1 + n2
    ranks = np.empty(n, dtype=np.float64)
    tie_cubed_sum = 0.0
    i = 0
    while i < n:
        j = i + 1
        while j < n and sv[j] == sv[i]:
            j += 1
        t_len = j - i
        if t_len > 1:
            tie_cubed_sum += float(t_len ** 3 - t_len)
        avg_r = 0.5 * (i + 1 + j)
        ranks[order[i:j]] = avg_r
        i = j
    R1 = float(np.sum(ranks[:n1]))
    U1 = R1 - 0.5 * n1 * (n1 + 1)
    U2 = float(n1 * n2) - U1
    U_min = min(U1, U2)
    r_rb = (U1 - U2) / float(n1 * n2)  # >0 means x tends to exceed y

    if tie_cubed_sum == 0.0 and n <= 36:
        u_int = int(round(U_min))
        dp = np.zeros((n1 + 1, n1 * n2 + 1), dtype=np.float64)
        dp[0, 0] = 1.0
        for item_idx in range(1, n + 1):
            max_k = min(item_idx, n1)
            for k in range(max_k, 0, -1):
                shift = item_idx - k
                if shift <= n1 * n2:
                    dp[k, shift:] += dp[k - 1, :n1 * n2 + 1 - shift]
        tot_comb = float(math.comb(n, n1))
        tail_le = float(np.sum(dp[n1, :u_int + 1])) / tot_comb
        p_val = min(1.0, 2.0 * tail_le) if u_int * 2 < n1 * n2 else 1.0
        method = "exact_dp"
    else:
        mu_u = 0.5 * n1 * n2
        var_u = (n1 * n2 / 12.0) * ((n + 1.0) - tie_cubed_sum / max(float(n * (n - 1)), 1.0))
        if var_u <= 1e-30:
            p_val = 1.0
        else:
            z = (abs(U1 - mu_u) - 0.5) / math.sqrt(var_u)
            p_val = min(1.0, math.erfc(max(z, 0.0) / math.sqrt(2.0)))
        method = "tie_corrected_normal"

    p_min_disc = min(1.0, 2.0 / float(math.comb(n, n1)))
    return {
        "n1": n1, "n2": n2,
        "median1": float(np.median(xa)), "median2": float(np.median(ya)),
        "U1": float(U1), "U2": float(U2),
        "rank_biserial_r": float(r_rb),
        "p_twosided": float(p_val),
        "p_min_discrete": float(p_min_disc),
        "method": method,
    }


def compute_virome_onco_group_tests_v9(viral_results, onco_results, t2t_rows=None):
    """[THM-V9] Automated nonparametric group contrasts across the 33 viruses & 28 human loci."""
    tests = []
    if viral_results:
        int_classes = {"retrovirus_int", "dna_onco_small", "hepatitis_hcc", "herpes_large", "vector_dna"}
        v_int = [r for r in viral_results if r["vclass"] in int_classes]
        v_rna = [r for r in viral_results if r["vclass"] == "rna_acute_ctrl"]
        v_hr = [r for r in viral_results if r["vclass"] == "dna_onco_small"]
        v_lr = [r for r in viral_results if r["vclass"] == "dna_benign_ctrl"]
        v_hum = [r for r in viral_results if r["vclass"] != "phage_outgroup"]
        v_phg = [r for r in viral_results if r["vclass"] == "phage_outgroup"]

        contrasts_v = [
            ("Onco/Integrating Viruses vs Acute RNA Controls", v_int, v_rna),
            ("High-Risk DNA Oncoviruses vs Benign Controls", v_hr, v_lr),
            ("Human Viruses (n=30) vs Bacteriophages (n=3)", v_hum, v_phg),
        ]
        metrics_v = [
            ("cos_Wmp_vs_human_grand", lambda r: r.get("cos_Wmp_vs_human_grand")),
            ("cos_Wmp_unb_vs_human_grand", lambda r: (r.get("cos_Wmp_unb_V18") or {}).get("cos_unb")),
            ("cos_Wmp_nonenv3600_vs_human", lambda r: (r.get("beyond_dinuc_V11") or {}).get("cos_nonenv3600_vs_human")),
            ("cos_Wmp_cpgfree_vs_human", lambda r: (r.get("beyond_dinuc_V11") or {}).get("cos_cpg_free_vs_human")),
            ("cos_Wmp_beyond_M1_vs_human", lambda r: (r.get("beyond_dinuc_V11") or {}).get("cos_markov1_resid_vs_human")),
            ("cos_Wmp_beyond_M2_vs_human", lambda r: (r.get("beyond_dinuc_V11") or {}).get("cos_markov2_resid_vs_human")),
            ("pure_nonenv_3600_unb", lambda r: (r.get("arrow_D11") or {}).get("pure_nonenv_arrow_share_unb")),
            ("z_iid_W_mm", lambda r: r["klein_d3"].get("z_iid_W_mm")),
            ("cpg_free_share", lambda r: (r.get("arrow_D11") or {}).get("cpg_free_arrow_share")),
        ]
        for cname, g1, g2 in contrasts_v:
            for mname, fn in metrics_v:
                res = exact_mann_whitney_u([fn(r) for r in g1], [fn(r) for r in g2])
                if res is not None:
                    res["contrast"] = cname
                    res["metric"] = mname
                    if "Bacteriophages" in cname and mname == "cos_Wmp_vs_human_grand":
                        res["note"] = "Compositional in raw 1024-D; null beyond Markov-2"
                    tests.append(res)

    if onco_results:
        o_drv = [r for r in onco_results if not str(r["role"]).startswith("control")]
        o_ctl = [r for r in onco_results if str(r["role"]).startswith("control")]
        metrics_o = [
            ("APOBEC_SBS2_chi_mm_unb_x1e4", lambda r: ((r.get("mutagenesis_V3") or {}).get("SBS2_APOBEC_plus") or {}).get("chi_mm_unb_x1e4")),
            ("APOBEC_SBS2_denovo_chi_mm_x1e4", lambda r: ((r.get("mutagenesis_V3") or {}).get("SBS2_APOBEC_plus") or {}).get("chi_mm_denovo_unb_x1e4")),
            ("APOBEC_SBS2_unb_chiral_share", lambda r: ((r.get("mutagenesis_V3") or {}).get("SBS2_APOBEC_plus") or {}).get("J_chiral_share_of_directed_unb")),
            ("pure_nonenv_3600_unb", lambda r: (r.get("arrow_D11") or {}).get("pure_nonenv_arrow_share_unb")),
            ("pal_d6_Wmp_share", lambda r: ((r.get("palindrome_spectrum") or {}).get("6") or {}).get("pal_Wmp_energy_share")),
        ]
        for mname, fn in metrics_o:
            res = exact_mann_whitney_u([fn(r) for r in o_drv], [fn(r) for r in o_ctl])
            if res is not None:
                res["contrast"] = "Cancer Drivers/Fragile (n=24) vs Non-Cancer Controls (n=4)"
                res["metric"] = mname
                tests.append(res)
    if t2t_rows:
        t_drv = [r for r in t2t_rows if not str(r["role"]).startswith("control")]
        t_ctl = [r for r in t2t_rows if str(r["role"]).startswith("control")]
        for mname, fn in (
            ("T2T_100kb_NonEnv_3600_unb", lambda r: r.get("pure_nonenv_3600_unb", r.get("pure_nonenv_3600"))),
            ("T2T_100kb_Pal_d6_share", lambda r: r.get("pal_d6_Wmp_share")),
            ("T2T_100kb_cos_vs_chr_euch", lambda r: r.get("cos_vs_chr_euch")),
        ):
            res = exact_mann_whitney_u([fn(r) for r in t_drv], [fn(r) for r in t_ctl])
            if res is not None:
                res["contrast"] = "Length-Matched 100-kb T2T: Drivers vs Controls"
                res["metric"] = mname
                tests.append(res)
    return tests


def self_test_virome_onco_v22(verbose=False):
    """[THM-V8..V9] 4 additional deterministic Zero-RNG self-tests for v2.2.0 (60/60 total)."""
    ok = {}
    # 1. Exact Mann-Whitney U DP check on known sample [10, 12, 15] vs [1, 2, 3, 4] -> U1=12, U2=0, p=2/35
    mw = exact_mann_whitney_u([10.0, 12.0, 15.0], [1.0, 2.0, 3.0, 4.0])
    assert mw is not None and mw["U2"] == 0.0 and abs(mw["p_twosided"] - (2.0 / 35.0)) < 1e-12, mw
    ok["thmV9_exact_mann_whitney_dp_err"] = abs(mw["p_twosided"] - (2.0 / 35.0))

    # 2. Genetic code table completeness & stop codon integrity
    assert len(_GENETIC_CODE) == 64 and _CODON_AA_ARR.size == 64
    assert sum(1 for aa in _CODON_AA_ARR if aa == "*") == 3
    ok["thmV8_genetic_code_64_ok"] = True

    # 3. Synthetic ORF synonymous wobble audit: ATG + 80x(TTG/TTA Leu) + TAA
    # TTG (Leu) -> TTA (Leu) is synonymous AND creates antisense -0 stop (TTA <-> TAA)!
    orf_syn = np.tile(np.array([3, 3, 2], dtype=np.uint8), 120)  # TTG repeats (366 bp > 300)
    orf_full = np.concatenate([np.array([0, 3, 2], dtype=np.uint8), orf_syn, np.array([3, 0, 0], dtype=np.uint8)])
    v8_res = audit_synonymous_duplex_barrier_v8(
        {"SYN_LEU": orf_full},
        [{"label": "SYN_LEU", "chr": "chr17", "role": "oncogene"}],
    )
    assert len(v8_res) == 1 and v8_res[0]["antisense_m0_stop_flip_share"] > 0.9, v8_res
    ok["thmV8_synonymous_antisense_stop_flip_ok"] = True

    # 4. Catalog count cross-check
    assert len(VIRAL_CATALOG_33) == 33 and len(ONCO_LOCI_28) == 28
    ok["thmV9_catalog_counts_verified"] = True
    ok.update(self_test_virome_onco_v23(verbose=False))
    if verbose:
        print(json.dumps(_clean(ok), indent=2))
    return ok


def append_v22_markdown_and_figures(payload, out_dir, md_path, dpi=300):
    v8_rows = payload.get("synonymous_duplex_barrier_V8") or []
    v9_rows = payload.get("group_discrimination_V9") or []
    extra_md = []
    if v8_rows:
        extra_md += [
            "", "## 8. Theorem V8 — 6-Frame Duplex Stop Barrier & Synonymous Mutation Toric Disruption", "",
            "| Locus | Chr | Role | 6-Frame Stop % | Longest ORF (codons) | Synonymous Hits | Antisense (-0) Stop Flip % | Mean Δ‖B_-+^3600‖² / Syn (count²) |",
            "|---|---|---|---:|---:|---:|---:|---:|",
        ]
        for r in v8_rows:
            extra_md.append(
                f"| **{r['locus']}** | {r['chr']} | {r['role']} | "
                f"{fnum(100.0*_fval(r.get('mean_6frame_stop_freq'), np.nan),2)}% | "
                f"{r['longest_orf_codons']:,} | {r['n_synonymous_wobble_hits']:,} | "
                f"{fnum(100.0*_fval(r.get('antisense_m0_stop_flip_share'), np.nan),2)}% | "
                f"{fnum(r.get('mean_dBmp_3600_count2_per_syn'),3)} |"
            )

    if v9_rows:
        extra_md += [
            "", "## 9. Theorem V9 — Nonparametric Mann–Whitney Group Discrimination (Zero-RNG)", "",
            "| Contrast | Metric | n1 / n2 | Median 1 vs Median 2 | Rank-Biserial r | Two-Sided p |",
            "|---|---|---:|---|---:|---:|",
        ]
        for r in v9_rows:
            extra_md.append(
                f"| {r['contrast']} | `{r['metric']}` | {r['n1']} / {r['n2']} | "
                f"{fnum(r['median1'],4)} vs {fnum(r['median2'],4)} | "
                f"{fnum(r['rank_biserial_r'],3,True)} | {_fsci(r['p_twosided'])} |"
            )

    if extra_md:
        old_md = Path(md_path).read_text(encoding="utf-8")
        atomic_write(md_path, old_md.rstrip() + "\n" + "\n".join(extra_md) + "\n")

    if not HAVE_MPL:
        return []
    fd = Path(out_dir) / "figures"
    fd.mkdir(parents=True, exist_ok=True)
    saved = []
    vr = payload.get("viruses") or []
    if vr and v8_rows:
        fig, axes = plt.subplots(1, 2, figsize=(18.5, 6.8))
        fig.suptitle(f"Theorems V8 & V9: Viral Class Toric Separation & Synonymous Duplex-Barrier Audit (v{VERSION})",
                     fontsize=13.5, fontweight="bold", y=0.98)
        ax = axes[0]
        class_order = ["retrovirus_int", "dna_onco_small", "hepatitis_hcc", "herpes_large",
                       "vector_dna", "dna_benign_ctrl", "rna_acute_ctrl", "phage_outgroup"]
        data_boxes, box_lbls = [], []
        for cls in class_order:
            vals = [_fval(r.get("cos_Wmp_vs_human_grand"), np.nan) for r in vr if r["vclass"] == cls]
            vals = [v for v in vals if _finite(v)]
            if vals:
                data_boxes.append(vals)
                box_lbls.append(f"{cls}\n(n={len(vals)})")
        if data_boxes:
            ax.boxplot(data_boxes, patch_artist=True,
                       boxprops=dict(facecolor="#cfe2f3", edgecolor="#1f77b4"),
                       medianprops=dict(color="#d62728", lw=2.0))
            ax.set_xticklabels(box_lbls)
            ax.axhline(0.0, color="k", ls=":", lw=0.9)
            ax.tick_params(axis="x", rotation=25, labelsize=8)
            ax.set_ylabel("cos( W_-+^virus , W_-+^human_euch )")
            ax.set_title("A. Viral Class Alignment to Human Euchromatin Arrow")
            ax.grid(True, axis="y", ls=":", alpha=0.35)

        ax = axes[1]
        v8_sorted = sorted(v8_rows, key=lambda r: _fval(r.get("mean_dBmp_3600_count2_per_syn"), 0.0))[-16:]
        ys = np.arange(len(v8_sorted))
        vals_s = [_fval(r.get("mean_dBmp_3600_count2_per_syn"), 0.0) for r in v8_sorted]
        cols_s = ["#2ca02c" if str(r["role"]).startswith("control") else "#d62728" for r in v8_sorted]
        ax.barh(ys, vals_s, color=cols_s, alpha=0.88)
        ax.axvline(705.0 / 256.0, color="#1a1a1a", ls="--", lw=1.25, label="THM-V15 Generic Bound 705/256 (2.754)")
        ax.set_xlim(0.0, 2.95)
        ax.set_yticks(ys); ax.set_yticklabels([f"{r['locus']} ({r['chr']})" for r in v8_sorted], fontsize=8)
        ax.set_xlabel("Length-Invariant 3600-D Disruption ||ΔB_-+^3600||² per Synonymous Wobble Hit (count²)")
        ax.set_title("B. Synonymous Wobble Disruption in Coding ORFs (Red = Cancer, Green = Control)")
        ax.legend(fontsize=8.0, loc="lower right", framealpha=0.94); ax.grid(True, axis="x", ls=":", alpha=0.35)

        fig.tight_layout(rect=[0, 0, 1, 0.945])
        p5 = fd / "fig_V5_group_discrimination_and_synonymous.png"
        fig.savefig(p5, dpi=dpi, bbox_inches="tight"); plt.close(fig); saved.append(str(p5))
        print(f"  [FIG] {p5}")
    return saved




# =====================================================================
# 13. PAPER 5 EXPANSION (v2.6.0): GC-CONTROLLED APOBEC3 EVASION,
#     BEYOND-DINUCLEOTIDE MIMICRY & ALL-ORDERS FINITE-DOSE MUTAGENESIS
#     [THM-V10.2] GC- & Dinucleotide-Controlled APOBEC3 TCW O/E Evasion
#     [THM-V11]   Exact Markov-1 Dinucleotide Null & Beyond-Dinuc Mimicry
#     [THM-V12]   Exact All-Orders Finite-Dose E_eps[B] & Double-Hit Hessian
# =====================================================================
def count_dinuc_4x4(nt, chunk_size=4_000_000):
    """Exact 4x4 adjacent dinucleotide count matrix across valid ACGT bases (never crosses N)."""
    L = int(nt.size)
    D = np.zeros((4, 4), dtype=np.float64)
    if L < 2:
        return D
    for s0 in range(0, L - 1, int(chunk_size)):
        s1 = min(L - 1, s0 + int(chunk_size))
        a = nt[s0:s1]
        b = nt[s0 + 1:s1 + 1]
        ok = (a < 4) & (b < 4)
        if np.any(ok):
            code = (a[ok].astype(np.int32) << 2) | b[ok].astype(np.int32)
            D += np.bincount(code, minlength=16).reshape(4, 4)
    return D


def exact_markov1_from_dinuc_counts(D):
    """[THM-V11] Exact closed-form stationary Markov-1 64x64 triplet-pair probability table
    P_M1(u, v) = p1(a0) T(a1|a0) T(a2|a1) T(b0|a2) T(b1|b0) T(b2|b1) (Zero-RNG)."""
    D = np.asarray(D, dtype=np.float64)
    tot = float(D.sum())
    if tot <= 0.0:
        return np.full((NC, NC), 1.0 / float(NC * NC), dtype=np.float64)
    P2 = D / tot
    p1 = P2.sum(axis=1)
    T = P2 / np.maximum(p1[:, None], TINY)
    # Vectorized over all 64x64 triplet pairs
    a0, a1, a2 = C_N1, C_N2, C_N3
    p_left = p1[a0] * T[a0, a1] * T[a1, a2]          # (64,)
    p_right = T[a0, a1] * T[a1, a2]                  # (64,) = T(b1|b0)*T(b2|b1)
    P_m1 = p_left[:, None] * T[a2[:, None], a0[None, :]] * p_right[None, :]
    s = float(P_m1.sum())
    return (P_m1 / s) if s > 0.0 else P_m1


def exact_markov1_from_pair_table(B):
    """Extracts exact 4x4 dinucleotide counts from a 64x64 d=3 pair table B and builds P_M1."""
    B = np.asarray(B, dtype=np.float64)
    # Use suffix dinucleotide (a1,a2) and prefix dinucleotide (b0,b1) and junction (a2,b0)
    D = np.zeros((4, 4), dtype=np.float64)
    row_s = B.sum(axis=1)
    col_s = B.sum(axis=0)
    for x in range(NC):
        D[C_N1[x], C_N2[x]] += row_s[x] + col_s[x]
        D[C_N2[x], C_N3[x]] += row_s[x] + col_s[x]
    # Add junction dinucleotide (a2, b0)
    for a2_i in range(4):
        m_r = (C_N3 == a2_i)
        for b0_i in range(4):
            m_c = (C_N1 == b0_i)
            D[a2_i, b0_i] += float(B[np.ix_(m_r, m_c)].sum())
    return exact_markov1_from_dinuc_counts(D)


def exact_markov2_from_trinuc_counts(c64):
    """[THM-V11.2] Exact closed-form stationary Markov-2 (64-trinucleotide) 64x64 pair table:
        P_M2(a0 a1 a2, b0 b1 b2) = p3(a0 a1 a2) * T2(b0|a1 a2) * T2(b1|a2 b0) * T2(b2|b0 b1).
    Provably satisfies I_nonbridge(P_M2) == 0 identically (Zero-RNG)."""
    c = np.asarray(c64, dtype=np.float64)
    tot = float(c.sum())
    if tot <= 0.0:
        return np.full((NC, NC), 1.0 / float(NC * NC), dtype=np.float64)
    p3 = c / tot
    # T2[dinuc(a,b), c] = P(c | a, b)
    p3_mat = p3.reshape(16, 4)
    row_s = p3_mat.sum(axis=1, keepdims=True)
    T2 = np.where(row_s > 1e-30, p3_mat / np.maximum(row_s, TINY), 0.25)
    # For left u=(a0,a1,a2) and right v=(b0,b1,b2):
    # suf_u = 4*a1 + a2 = SUF16[u]; junc_uv = 4*a2 + b0; pre_v = 4*b0 + b1 = OUT16[v]
    junc_uv = (C_N3[:, None] << 2) | C_N1[None, :]
    P_m2 = (
        p3[:, None]
        * T2[SUF16[:, None], C_N1[None, :]]
        * T2[junc_uv, C_N2[None, :]]
        * T2[OUT16[None, :], C_N3[None, :]]
    )
    s = float(P_m2.sum())
    return (P_m2 / s) if s > 0.0 else P_m2


def exact_markov2_from_pair_table(B):
    """Fallback Markov-2 builder from 64x64 d=3 pair table B (averages row & col triplet marginals)."""
    B = np.asarray(B, dtype=np.float64)
    return exact_markov2_from_trinuc_counts(0.5 * (B.sum(axis=1) + B.sum(axis=0)))


def exact_spearman_rho(x, y):
    """Exact deterministic Zero-RNG tie-corrected Spearman rank correlation & two-sided t-approx p-value."""
    pairs = [(float(a), float(b)) for a, b in zip(x, y) if _finite(a) and _finite(b)]
    n = len(pairs)
    if n < 5:
        return None
    xa = np.array([p_[0] for p_ in pairs], dtype=np.float64)
    ya = np.array([p_[1] for p_ in pairs], dtype=np.float64)

    def _avg_ranks(v):
        ord_i = np.argsort(v, kind="mergesort")
        sv = v[ord_i]
        rk = np.empty(v.size, dtype=np.float64)
        i = 0
        while i < v.size:
            j = i + 1
            while j < v.size and sv[j] == sv[i]:
                j += 1
            rk[ord_i[i:j]] = 0.5 * (i + 1 + j)
            i = j
        return rk

    rx = _avg_ranks(xa)
    ry = _avg_ranks(ya)
    rx -= float(np.mean(rx))
    ry -= float(np.mean(ry))
    den = math.sqrt(float(np.dot(rx, rx)) * float(np.dot(ry, ry)))
    if den <= 1e-30:
        return {"n": n, "rho": 0.0, "p_twosided": 1.0}
    rho = max(-1.0, min(1.0, float(np.dot(rx, ry) / den)))
    if abs(rho) >= 1.0 - 1e-12:
        p_val = 0.0
    else:
        t_stat = abs(rho) * math.sqrt((n - 2.0) / max(1.0 - rho * rho, 1e-30))
        try:
            from scipy.stats import t as _sp_t
            p_val = float(2.0 * _sp_t.sf(t_stat, n - 2))
        except Exception:
            p_val = float(math.erfc(t_stat / math.sqrt(2.0)))
    return {"n": n, "rho": rho, "p_twosided": min(1.0, max(0.0, p_val))}


def residualized_group_test_v11(viral_results, group1_classes, group2_classes, metric_fn, include_cpg_oe=False):
    # [THM-V11.3 & V18] Exact Zero-RNG covariate-adjusted contrast controlling simultaneously
    # for viral GC%, log10(length_bp), and optionally ZAP CpG O/E via exact OLS residualization.
    sub = [r for r in viral_results if r["vclass"] in group1_classes or r["vclass"] in group2_classes]
    rows = []
    for r in sub:
        y = metric_fn(r)
        cpg = r.get("cpg_oe")
        if _finite(y) and ((not include_cpg_oe) or _finite(cpg)):
            rows.append((
                r["vclass"] in group1_classes,
                float(y),
                float(r["gc"]),
                math.log10(max(float(r["length_bp"]), 1.0)),
                float(cpg) if _finite(cpg) else 0.0,
            ))
    if len(rows) < 8:
        return None
    Y = np.array([t[1] for t in rows], dtype=np.float64)
    cols = [np.ones(len(rows)), [t[2] for t in rows], [t[3] for t in rows]]
    if include_cpg_oe:
        cols.append([t[4] for t in rows])
    X = np.column_stack(cols)
    beta, _, _, _ = np.linalg.lstsq(X, Y, rcond=None)
    resid = Y - X @ beta
    g1 = [float(resid[i]) for i, t in enumerate(rows) if t[0]]
    g2 = [float(resid[i]) for i, t in enumerate(rows) if not t[0]]
    mw = exact_mann_whitney_u(g1, g2)
    if mw is None:
        return None
    mw["beta_gc"] = float(beta[1])
    mw["beta_log10L"] = float(beta[2])
    if include_cpg_oe:
        mw["beta_cpg_oe"] = float(beta[3])
    mw["covariates"] = ["gc", "log10L", "cpg_oe"] if include_cpg_oe else ["gc", "log10L"]
    return mw


def compute_t2t_rank1_diagnostics_v18(viral_results, t2t_rows):
    # [THM-V18] Tests Rank-1 factorization cos(v, l) ~ cos(v, H) * cos(l, H) across all
    # 33 x 28 viral-neighborhood pairs and summarizes raw vs human-residual top winners.
    if not viral_results or not t2t_rows:
        return None
    v_h = {r["label"]: r.get("cos_Wmp_vs_human_grand") for r in viral_results}
    xs, ys = [], []
    raw_cnt, res_cnt = {}, {}
    for row in t2t_rows:
        tv = row.get("top_matching_virus")
        if tv:
            raw_cnt[tv] = raw_cnt.get(tv, 0) + 1
        rv = row.get("top_residual_virus")
        if rv:
            res_cnt[rv] = res_cnt.get(rv, 0) + 1
        clh = row.get("cos_vs_human_grand")
        vmap = row.get("viral_cos_Wmp") or {}
        if not _finite(clh):
            continue
        for vn, cvl in vmap.items():
            cvh = v_h.get(vn)
            if _finite(cvl) and _finite(cvh):
                xs.append(float(cvh) * float(clh))
                ys.append(float(cvl))
    if len(xs) < 10:
        return None
    X = np.column_stack([np.ones(len(xs)), xs])
    Y = np.asarray(ys, dtype=np.float64)
    beta, _, _, _ = np.linalg.lstsq(X, Y, rcond=None)
    pred = X @ beta
    ss_res = float(np.sum((Y - pred) ** 2))
    ss_tot = float(np.sum((Y - float(np.mean(Y))) ** 2))
    r2 = 1.0 - ss_res / max(ss_tot, 1e-30)
    return {
        "n_pairs": len(xs),
        "rank1_r2": float(r2),
        "rank1_slope": float(beta[1]),
        "rank1_intercept": float(beta[0]),
        "raw_top_virus_counts": dict(sorted(raw_cnt.items(), key=lambda kv: -kv[1])),
        "residual_top_virus_counts": dict(sorted(res_cnt.items(), key=lambda kv: -kv[1])),
    }


def klein_orbit_attribution_v13(W_a, W_b=None, top_k=10):
    """[THM-V13] Exact Zero-RNG Klein V_4 orbit attribution of Frobenius inner product <W_a, W_b>_F
    (or self-energy ||W_a||^2 when W_b is None). Partition of total inner product sums to 100.00%."""
    A = np.asarray(W_a, dtype=np.float64)
    B = A if W_b is None else np.asarray(W_b, dtype=np.float64)
    prod = A * B
    tot = float(np.sum(prod))
    if abs(tot) <= 1e-36:
        return {"total_inner": 0.0, "sum_orbit_shares_pct": 0.0, "top_orbits": []}
    seen = set()
    orbits = []
    for u in range(NC):
        for v in range(NC):
            if u == v:
                continue
            ru, rv = int(MIRROR[v]), int(MIRROR[u])
            canon = min((u, v), (v, u), (ru, rv), (rv, ru))
            if canon in seen:
                continue
            seen.add(canon)
            members = {(u, v), (v, u), (ru, rv), (rv, ru)}
            contrib = sum(float(prod[a_, b_]) for a_, b_ in members)
            if abs(contrib) <= 1e-36:
                continue
            # Pick positive-flux orientation in A
            best_m = max(members, key=lambda ab: float(A[ab[0], ab[1]]))
            mu, mv = best_m
            mru, mrv = int(MIRROR[mv]), int(MIRROR[mu])
            orbits.append({
                "motif_6mer": TRIPLETS[mu] + "-" + TRIPLETS[mv],
                "mirror_6mer": TRIPLETS[mru] + "-" + TRIPLETS[mrv],
                "orbit_size": len(members),
                "has_CpG": bool(_ANY_CG_PAIR_MASK[mu, mv]),
                "contrib": contrib,
                "share_pct": 100.0 * contrib / tot,
            })
    orbits.sort(key=lambda d: d["contrib"], reverse=True)
    sum_all_pct = sum(o["share_pct"] for o in orbits)
    return {
        "total_inner": tot,
        "n_orbits": len(orbits),
        "sum_orbit_shares_pct": sum_all_pct,
        "top_orbits": orbits[:int(top_k)],
    }


def viral_beyond_dinuc_mimicry_v11(nt, B3, human_ref):
    """[THM-V11 & THM-V11.2] Decomposes viral-to-human toric alignment across 6 strict controls:
    1) Full 1024-D W_-+; 2) 496-D bridge envelope; 3) 3600-D non-envelope;
    4) Strictly CpG-free; 5) Beyond-Markov-1 (dinucleotide) residual;
    6) Beyond-Markov-2 (full 64-trinucleotide / codon-usage) residual."""
    if human_ref is None:
        return None
    kv = klein_v4_pr2_audit(B3)
    W_v = kv["_W_mp"]
    W_h = human_ref["grand_W_mp"]

    W_v_env = pi_env_project(W_v)
    W_v_ne = W_v - W_v_env
    W_h_env = pi_env_project(W_h)
    W_h_ne = W_h - W_h_env

    no_cg = ~_ANY_CG_PAIR_MASK
    W_v_ncg = np.where(no_cg, W_v, 0.0)
    W_h_ncg = np.where(no_cg, W_h, 0.0)
    W_v_ne_ncg = np.where(no_cg, W_v_ne, 0.0)
    W_h_ne_ncg = np.where(no_cg, W_h_ne, 0.0)

    D_v = count_dinuc_4x4(nt)
    P_m1_v = exact_markov1_from_dinuc_counts(D_v)
    W_v_m1 = klein_v4_pr2_audit(P_m1_v)["_W_mp"]
    W_v_perp_m1 = W_v - W_v_m1
    W_h_perp_m1 = human_ref.get("grand_W_mp_perp_m1")

    T3_v = count_triplets_1d(nt)
    P_m2_v = exact_markov2_from_trinuc_counts(T3_v)
    W_v_m2 = klein_v4_pr2_audit(P_m2_v)["_W_mp"]
    W_v_perp_m2 = W_v - W_v_m2
    W_h_perp_m2 = human_ref.get("grand_W_mp_perp_m2")
    W_h_m2 = human_ref.get("grand_W_mp_m2")
    W_v_gs_m2 = W_v - (float(np.sum(W_v * W_v_m2)) / max(float(np.sum(W_v_m2 * W_v_m2)), TINY)) * W_v_m2
    W_h_gs_m2 = (
        W_h - (float(np.sum(W_h * W_h_m2)) / max(float(np.sum(W_h_m2 * W_h_m2)), TINY)) * W_h_m2
        if W_h_m2 is not None else None
    )

    return {
        "cos_total_1024_vs_human": _cos_mat_frob(W_v, W_h),
        "cos_env496_vs_human": _cos_mat_frob(W_v_env, W_h_env),
        "cos_nonenv3600_vs_human": _cos_mat_frob(W_v_ne, W_h_ne),
        "cos_cpg_free_vs_human": _cos_mat_frob(W_v_ncg, W_h_ncg),
        "cos_nonenv3600_cpg_free_vs_human": _cos_mat_frob(W_v_ne_ncg, W_h_ne_ncg),
        "cos_markov1_only_vs_human": _cos_mat_frob(W_v_m1, human_ref.get("grand_W_mp_m1")) if human_ref.get("grand_W_mp_m1") is not None else None,
        "cos_markov1_resid_vs_human": _cos_mat_frob(W_v_perp_m1, W_h_perp_m1) if W_h_perp_m1 is not None else None,
        "cos_markov2_only_vs_human": _cos_mat_frob(W_v_m2, W_h_m2) if W_h_m2 is not None else None,
        "cos_markov2_resid_vs_human": _cos_mat_frob(W_v_perp_m2, W_h_perp_m2) if W_h_perp_m2 is not None else None,
        "cos_markov2_gs_orth_vs_human": _cos_mat_frob(W_v_gs_m2, W_h_gs_m2) if W_h_gs_m2 is not None else None,
        "cos_self_vs_markov1": _cos_mat_frob(W_v, W_v_m1),
        "cos_self_vs_markov2": _cos_mat_frob(W_v, W_v_m2),
    }


def exact_finite_dose_polynomial_v12(nt, rate_tbl, chunk_size=500_000):
    """[THM-V12] Exact O(L) Zero-RNG all-orders polynomial expansion of independent-site
    mutagenesis expectation: E_eps[B] = sum_{m=0}^6 eps^m * C_m, valid for ALL eps in [0, 1].
    C_0 == B_0, C_1 == J_S (exact Jacobian), 2*C_2 == H_S (exact double-hit kataegis Hessian)."""
    L = int(nt.size)
    C_flat = np.zeros((7, NC * NC), dtype=np.float64)
    if L < 8:
        return C_flat.reshape(7, NC, NC)
    row_nz = np.count_nonzero(rate_tbl > 0.0, axis=1)
    assert int(row_nz.max()) <= 1, "rate_tbl must have <= 1 target base per context"
    target_b = np.argmax(rate_tbl, axis=1).astype(np.int32)
    target_w = rate_tbl[np.arange(64), target_b]

    max_i = L - 6
    shifts = np.array([10, 8, 6, 4, 2, 0], dtype=np.int32)
    clear_masks = ~(3 << shifts)

    for s0 in range(1, max_i, int(chunk_size)):
        s1 = min(max_i, s0 + int(chunk_size))
        if s1 <= s0:
            continue
        ok = (
            (nt[s0 - 1:s1 - 1] < 4) & (nt[s0:s1] < 4) & (nt[s0 + 1:s1 + 1] < 4)
            & (nt[s0 + 2:s1 + 2] < 4) & (nt[s0 + 3:s1 + 3] < 4)
            & (nt[s0 + 4:s1 + 4] < 4) & (nt[s0 + 5:s1 + 5] < 4)
            & (nt[s0 + 6:s1 + 6] < 4)
        )
        if not np.any(ok):
            continue
        idx = np.flatnonzero(ok) + s0
        b_arr = [nt[idx + p].astype(np.int32) for p in range(-1, 7)]
        pair_old = (
            (b_arr[1] << 10) | (b_arr[2] << 8) | (b_arr[3] << 6)
            | (b_arr[4] << 4) | (b_arr[5] << 2) | b_arr[6]
        )
        ws = np.empty((idx.size, 6), dtype=np.float64)
        bns = np.empty((idx.size, 6), dtype=np.int32)
        mask6 = np.zeros(idx.size, dtype=np.int32)
        for p in range(6):
            ctx = (b_arr[p] << 4) | (b_arr[p + 1] << 2) | b_arr[p + 2]
            wp = target_w[ctx]
            ws[:, p] = wp
            bns[:, p] = target_b[ctx]
            mask6 |= ((wp > 0.0).astype(np.int32) << p)

        # Fast path for windows with 0 mutable sites
        m0 = (mask6 == 0)
        if np.any(m0):
            C_flat[0] += np.bincount(pair_old[m0], minlength=NC * NC)

        # Group remaining windows by active 6-bit mask (1..63)
        nz_masks = np.unique(mask6[~m0])
        for mk in nz_masks:
            mk_i = int(mk)
            sel = np.flatnonzero(mask6 == mk_i)
            po_s = pair_old[sel]
            ws_s = ws[sel]
            bns_s = bns[sel]
            active = [p for p in range(6) if (mk_i >> p) & 1]
            k_act = len(active)
            for sub in range(1 << k_act):
                mut_pos = [active[j] for j in range(k_act) if (sub >> j) & 1]
                unm_pos = [active[j] for j in range(k_act) if not ((sub >> j) & 1)]
                deg_s = len(mut_pos)
                pair_s = po_s.copy()
                w_prod = np.ones(sel.size, dtype=np.float64)
                for p in mut_pos:
                    pair_s = (pair_s & clear_masks[p]) | (bns_s[:, p] << shifts[p])
                    w_prod *= ws_s[:, p]
                poly = np.zeros((sel.size, k_act + 1), dtype=np.float64)
                poly[:, deg_s] = w_prod
                for q in unm_pos:
                    wq = ws_s[:, q]
                    for d_i in range(k_act, 0, -1):
                        poly[:, d_i] -= wq * poly[:, d_i - 1]
                for m_deg in range(deg_s, k_act + 1):
                    col_w = poly[:, m_deg]
                    if np.any(col_w != 0.0):
                        C_flat[m_deg] += np.bincount(pair_s, weights=col_w, minlength=NC * NC)

    return C_flat.reshape(7, NC, NC)


def analyze_finite_dose_curve_v12(nt, sig_names=("SBS1_CpG_sym", "SBS2_APOBEC_plus", "SBS13_APOBEC_plus", "SBS7a_UV_plus", "SBS4_Tobacco_plus"),
                                  doses=(0.0, 0.01, 0.02, 0.05, 0.10, 0.25, 0.50, 1.00)):
    """[THM-V12] Evaluates exact all-orders finite-dose mutagenesis curves & double-hit
    kataegis Hessian curvature across doses eps in [0, 1] (Zero-RNG)."""
    out = {}
    for sname in sig_names:
        smeta = COSMIC_SIG_TABLES[sname]
        C = exact_finite_dose_polynomial_v12(nt, smeta["rates"])
        N = float(C[0].sum())
        if N <= 10.0:
            continue
        H_mat = 2.0 * C[2] / N
        H_mm = 0.25 * (H_mat - H_mat.T - H_mat.T[np.ix_(MIRROR, MIRROR)] + H_mat[np.ix_(MIRROR, MIRROR)])
        J_mat = C[1] / N
        J_mm = 0.25 * (J_mat - J_mat.T - J_mat.T[np.ix_(MIRROR, MIRROR)] + J_mat[np.ix_(MIRROR, MIRROR)])
        # Evaluate exact trajectory & linear-vs-exact deviation
        curve = []
        for eps in doses:
            e = float(eps)
            B_ex = C[0].copy()
            epow = e
            for m in range(1, 7):
                B_ex += epow * C[m]
                epow *= e
            B_ex = np.maximum(B_ex, 0.0)
            kv_ex = klein_v4_pr2_audit(B_ex)
            B_lin = np.maximum(C[0] + e * C[1], 0.0)
            kv_lin = klein_v4_pr2_audit(B_lin)
            nl_dev = abs(kv_ex["norm2_W_mm"] - kv_lin["norm2_W_mm"]) / max(kv_ex["norm2_W_mm"], 1e-30)
            curve.append({
                "eps": e,
                "mass_err": abs(float(B_ex.sum()) - N),
                "norm2_W_mp": kv_ex["norm2_W_mp"],
                "norm2_W_mm": kv_ex["norm2_W_mm"],
                "norm2_W_mp_unb": kv_ex["norm2_W_mp_unb"],
                "norm2_W_mm_unb": kv_ex["norm2_W_mm_unb"],
                "PR2_ok_expect": kv_ex["time_PR2_compliant_share"],
                "PR2_ok_unb": kv_ex["time_PR2_compliant_share_unb"],
                "nonlinear_Wmm_rel_dev": nl_dev,
            })
        out[sname] = {
            "N_interior": int(N),
            "hessian_W_mm_norm2_x1e4": 1e4 * float(np.sum(H_mm * H_mm)),
            "cos_J_mm_vs_H_mm": _cos_mat_frob(J_mm, H_mm),
            "curve": curve,
        }
    return out


def compute_viral_apobec_restriction_v10(viral_results, viral_nt_cache=None):
    """[THM-V10 & THM-V10.2] Summarizes APOBEC3A/B (SBS2/SBS13) & CpG (SBS1) susceptibility
    across all 33 viruses with GC-controlled (base O/E) and Markov-1 dinucleotide-controlled
    (dinuc O/E) TCW-motif evolutionary depletion metrics."""
    rows = []
    for r in viral_results:
        mj = r.get("mutagenesis_V3") or {}
        s2 = mj.get("SBS2_APOBEC_plus") or {}
        s13 = mj.get("SBS13_APOBEC_plus") or {}
        s1 = mj.get("SBS1_CpG_sym") or {}
        if not s2 or s2.get("skipped"):
            continue
        tcw_obs = _fval(s2.get("eligible_site_frac"), np.nan)
        tcw_oe_base = tcw_oe_dinuc = cpg_oe = None
        if viral_nt_cache and r["label"] in viral_nt_cache:
            vnt = viral_nt_cache[r["label"]]
            D4 = count_dinuc_4x4(vnt)
            s_d = float(D4.sum())
            if s_d > 10.0:
                P2 = D4 / s_d
                p1 = P2.sum(axis=1)
                exp_base = float(p1[3] * p1[1] * (p1[0] + p1[3]))
                exp_dinuc = float(P2[3, 1] * (P2[1, 0] + P2[1, 3]) / max(float(p1[1]), TINY))
                if exp_base > 1e-12 and _finite(tcw_obs):
                    tcw_oe_base = float(tcw_obs) / exp_base
                if exp_dinuc > 1e-12 and _finite(tcw_obs):
                    tcw_oe_dinuc = float(tcw_obs) / exp_dinuc
                exp_cpg = float(p1[1] * p1[2])
                if exp_cpg > 1e-12:
                    cpg_oe = float(P2[1, 2]) / exp_cpg
        rows.append({
            "virus": r["label"],
            "vclass": r["vclass"],
            "mode": r["mode"],
            "length_bp": r["length_bp"],
            "gc": r["gc"],
            "tcw_eligible_frac": s2.get("eligible_site_frac"),
            "tcw_oe_base": tcw_oe_base,
            "tcw_oe_dinuc": tcw_oe_dinuc,
            "cpg_oe": cpg_oe,
            "sbs2_chiral_share_unb": s2.get("J_chiral_share_of_directed_unb"),
            "sbs2_chi_mm_x1e4": s2.get("chi_mm_unb_x1e4"),
            "sbs2_chi_mm_denovo_x1e4": s2.get("chi_mm_denovo_unb_x1e4"),
            "sbs2_eps_crit": s2.get("eps_crit_4N"),
            "sbs2_cos_Wmp": s2.get("cos_Wmp_native_vs_J"),
            "sbs13_chi_mm_x1e4": s13.get("chi_mm_unb_x1e4"),
            "sbs1_cpg_eligible_frac": s1.get("eligible_site_frac"),
        })
    int_classes = {"retrovirus_int", "dna_onco_small", "hepatitis_hcc", "herpes_large", "vector_dna"}
    evade_classes = {"retrovirus_int", "dna_onco_small", "dna_benign_ctrl"}
    mw_tcw = exact_mann_whitney_u(
        [x["tcw_eligible_frac"] for x in rows if x["vclass"] in int_classes],
        [x["tcw_eligible_frac"] for x in rows if x["vclass"] == "rna_acute_ctrl"],
    )
    mw_oe_base = exact_mann_whitney_u(
        [x["tcw_oe_base"] for x in rows if x["vclass"] in int_classes],
        [x["tcw_oe_base"] for x in rows if x["vclass"] == "rna_acute_ctrl"],
    )
    mw_oe_dinuc = exact_mann_whitney_u(
        [x["tcw_oe_dinuc"] for x in rows if x["vclass"] in int_classes],
        [x["tcw_oe_dinuc"] for x in rows if x["vclass"] == "rna_acute_ctrl"],
    )
    mw_oe_evade = exact_mann_whitney_u(
        [x["tcw_oe_base"] for x in rows if x["vclass"] in evade_classes],
        [x["tcw_oe_base"] for x in rows if x["vclass"] == "rna_acute_ctrl"],
    )
    mw_snr = exact_mann_whitney_u(
        [x["sbs2_chi_mm_x1e4"] for x in rows if x["vclass"] in int_classes],
        [x["sbs2_chi_mm_x1e4"] for x in rows if x["vclass"] == "rna_acute_ctrl"],
    )
    return {
        "n_viruses": len(rows),
        "rows": rows,
        "mw_tcw_int_vs_rna": mw_tcw,
        "mw_tcw_oe_base_int_vs_rna": mw_oe_base,
        "mw_tcw_oe_dinuc_int_vs_rna": mw_oe_dinuc,
        "mw_tcw_oe_base_retro_pap_vs_rna": mw_oe_evade,
        "mw_sbs2_snr_int_vs_rna": mw_snr,
    }


def self_test_virome_onco_v23(verbose=False):
    """[THM-V10] 2 additional deterministic Zero-RNG self-tests for v2.3.0."""
    ok = {}
    tcw_half = np.tile(np.array([3, 1, 0, 3, 1, 3, 2, 0], dtype=np.uint8), 120)
    tcw_pr2 = np.concatenate([tcw_half, np.full(20, 255, dtype=np.uint8), reverse_complement_nt(tcw_half)])
    mj = analyze_mutagenesis_signatures(tcw_pr2, eps_eval=0.05)
    assert mj["SBS2_APOBEC_plus"]["eligible_site_frac"] > 0.10, mj["SBS2_APOBEC_plus"]
    assert mj["SBS2_APOBEC_plus"]["J_norm2_W_mm"] > 1e-3, mj["SBS2_APOBEC_plus"]
    assert mj["SBS2_APOBEC_plus"]["delta_W_mm_over_human_euch"] > 1e4, mj["SBS2_APOBEC_plus"]
    ok["thmV10_viral_apobec_tcw_restriction_ok"] = True

    g_seq = np.full(600, 2, dtype=np.uint8)
    mj_g = analyze_mutagenesis_signatures(g_seq, eps_eval=0.05)
    assert mj_g["SBS2_APOBEC_plus"]["eligible_site_frac"] == 0.0
    assert mj_g["SBS2_APOBEC_plus"]["J_norm2_W_mm"] == 0.0
    ok["thmV10_zero_tcw_exact_null_ok"] = True
    ok.update(self_test_virome_onco_v25(verbose=False))
    ok.update(self_test_virome_onco_v26(verbose=False))
    if verbose:
        print(json.dumps(_clean(ok), indent=2))
    return ok


def self_test_virome_onco_v26(verbose=False):
    """[THM-V3.3, V2, V5, V10.2, V11, V12] 7 additional exact Zero-RNG self-tests for v2.6.0 (72/72 total)."""
    ok = {}

    # 1. THM-V3.3: Symmetric SBS1 has IDENTICALLY ZERO de novo chiral injection even on PR2-broken DNA
    _ix = np.arange(3000, dtype=np.uint64) + np.uint64(409)
    _zz = (_ix ^ (_ix >> np.uint64(30))) * np.uint64(0xbf58476d1ce4e5b9)
    _zz = (_zz ^ (_zz >> np.uint64(27))) * np.uint64(0x94d049bb133111eb)
    _u = ((_zz ^ (_zz >> np.uint64(31))) % np.uint64(10)).astype(np.uint8)
    # Strongly G/C-skewed PR2-broken synthetic sequence (50% G, 20% C, 20% T, 10% A)
    s_asym = np.where(_u < 5, 2, np.where(_u < 7, 1, np.where(_u < 9, 3, 0))).astype(np.uint8)
    mj_as = analyze_mutagenesis_signatures(s_asym, {
        "SBS1_CpG_sym": COSMIC_SIG_TABLES["SBS1_CpG_sym"],
        "SBS2_APOBEC_plus": COSMIC_SIG_TABLES["SBS2_APOBEC_plus"],
    })
    s1_as = mj_as["SBS1_CpG_sym"]
    ap_as = mj_as["SBS2_APOBEC_plus"]
    assert s1_as["chi_mm_denovo_x1e4"] < 1e-26 and s1_as["chi_mm_sub_x1e4"] > 1e-4, s1_as
    assert ap_as["chi_mm_denovo_unb_x1e4"] > 0.01 and ap_as["denovo_chiral_share"] > 0.5, ap_as
    ok["thmV3_3_symmetric_denovo_zero_err"] = s1_as["chi_mm_denovo_x1e4"]

    # 2. THM-V3.3 identity: J_denovo + J_sub == J identically
    rc_as = reverse_complement_nt(s_asym)
    _, Jf, _ = exact_signature_jacobian(s_asym, COSMIC_SIG_TABLES["SBS2_APOBEC_plus"]["rates"])
    _, Jr, _ = exact_signature_jacobian(rc_as, COSMIC_SIG_TABLES["SBS2_APOBEC_plus"]["rates"])
    id_err = float(np.max(np.abs((0.5 * (Jf + Jr) + 0.5 * (Jf - Jr)) - Jf)))
    assert id_err < 1e-15, id_err
    ok["thmV3_3_denovo_plus_sub_identity_err"] = id_err

    # 3. THM-V2: Dilution-free count-space antipodal law even on PR2-broken host
    v_syn = np.tile(np.array([2, 2, 0, 1, 3, 2, 0, 0, 3, 1, 2], dtype=np.uint8), 100)
    ins_as = deterministic_proviral_insertion_audit(s_asym, v_syn)
    assert ins_as["cos_dBmm_plus_vs_minus"] < -0.99 and ins_as["cos_dBmp_plus_vs_minus"] > 0.99, ins_as
    ok["thmV2_count_space_antipodal_ok"] = True

    # 4. THM-V5: SNP-robust multi-scale locator succeeds even with a SNP every 50 bp
    anc_snp = s_asym[slice(500, 1500)].copy()
    for p_snp in range(25, anc_snp.size, 50):
        anc_snp[p_snp] = (anc_snp[p_snp] + 1) % 4
    p_f, st_f = locate_anchor_on_chromosome(s_asym, anc_snp, anchor_len=60)
    p_r, st_r = locate_anchor_on_chromosome(s_asym, reverse_complement_nt(anc_snp), anchor_len=60)
    assert p_f == 500 and st_f == "+" and p_r == 500 and st_r == "-", (p_f, st_f, p_r, st_r)
    ok["thmV5_snp_robust_locator_ok"] = True

    # 5. THM-V11: Stationary Markov-1 64x64 table satisfies I_nonbridge(P_M1) == 0 identically
    P_m1 = exact_markov1_from_dinuc_counts(count_dinuc_4x4(s_asym))
    sc_m1 = two_sided_screen_full(P_m1 * 1e9)
    assert abs(float(P_m1.sum()) - 1.0) < 1e-14 and sc_m1["I_nonbridge"] < 1e-14, sc_m1
    ok["thmV11_markov1_nonbridge_zero_err"] = sc_m1["I_nonbridge"]

    # 6. THM-V12: Exact all-orders polynomial E_eps[B] matches exhaustive 2^m enumeration on tiny sequence
    s_mini = np.array([0, 3, 1, 0, 3, 1, 3, 2, 3, 1, 0, 3, 1, 3, 0, 2, 3, 1, 0, 1], dtype=np.uint8)
    rt_ap = COSMIC_SIG_TABLES["SBS2_APOBEC_plus"]["rates"]
    C_mini = exact_finite_dose_polynomial_v12(s_mini, rt_ap)
    B0_m, J_m, _ = exact_signature_jacobian(s_mini, rt_ap)
    assert float(np.max(np.abs(C_mini[0] - B0_m))) == 0.0
    assert float(np.max(np.abs(C_mini[1] - J_m))) < 1e-14
    # Exhaustive 2^m verification over all mutable interior sites of s_mini at eps=0.37
    eps_t = 0.37
    mut_sites = []
    for t_i in range(1, s_mini.size - 1):
        ctx_i = (int(s_mini[t_i - 1]) << 4) | (int(s_mini[t_i]) << 2) | int(s_mini[t_i + 1])
        bn_i = int(np.argmax(rt_ap[ctx_i]))
        w_i = float(rt_ap[ctx_i, bn_i])
        if w_i > 0.0:
            mut_sites.append((t_i, bn_i, w_i))
    assert 3 <= len(mut_sites) <= 10, len(mut_sites)
    B_exh = np.zeros((NC, NC), dtype=np.float64)
    for mask_e in range(1 << len(mut_sites)):
        sm_e = s_mini.copy()
        prob = 1.0
        for j_e, (t_i, bn_i, w_i) in enumerate(mut_sites):
            if (mask_e >> j_e) & 1:
                sm_e[t_i] = bn_i
                prob *= eps_t * w_i
            else:
                prob *= (1.0 - eps_t * w_i)
        B_exh += prob * count_triplet_pairs_lag(sm_e[slice(1, sm_e.size - 1)], d=3)
    B_poly = sum((eps_t ** m_) * C_mini[m_] for m_ in range(7))
    v12_err = float(np.max(np.abs(B_poly - B_exh)))
    assert v12_err < 1e-13, v12_err
    # Check PR2 preservation at eps=1.0 for symmetric SBS1 on PR2-exact duplex
    s_pr2_m = np.concatenate([s_asym[slice(0, 600)], np.full(20, 255, dtype=np.uint8), reverse_complement_nt(s_asym[slice(0, 600)])])
    C_s1 = exact_finite_dose_polynomial_v12(s_pr2_m, COSMIC_SIG_TABLES["SBS1_CpG_sym"]["rates"])
    B_s1_1 = sum(C_s1[m_] for m_ in range(7))
    kv_s1_1 = klein_v4_pr2_audit(B_s1_1)
    assert kv_s1_1["norm2_W_mm"] < 1e-32, kv_s1_1
    ok["thmV12_all_orders_vs_exhaustive_err"] = v12_err

    # 7. THM-V10.2: TCW O/E is GC-invariant (~1.0) on both 70%-GC and 30%-GC random sequences
    _ix_gc = np.arange(20000, dtype=np.uint64) + np.uint64(9001)
    _z_gc = (_ix_gc ^ (_ix_gc >> np.uint64(30))) * np.uint64(0xbf58476d1ce4e5b9)
    _z_gc = (_z_gc ^ (_z_gc >> np.uint64(27))) * np.uint64(0x94d049bb133111eb)
    _r20 = ((_z_gc ^ (_z_gc >> np.uint64(31))) % np.uint64(20)).astype(np.uint8)
    s_hi_gc = np.where(_r20 < 7, 1, np.where(_r20 < 14, 2, np.where(_r20 < 17, 0, 3))).astype(np.uint8)
    s_lo_gc = np.where(_r20 < 3, 1, np.where(_r20 < 6, 2, np.where(_r20 < 13, 0, 3))).astype(np.uint8)
    oes = []
    for sg in (s_hi_gc, s_lo_gc):
        D_g = count_dinuc_4x4(sg)
        P2_g = D_g / float(D_g.sum())
        p1_g = P2_g.sum(axis=1)
        obs_g = analyze_mutagenesis_signatures(sg, {"SBS2_APOBEC_plus": COSMIC_SIG_TABLES["SBS2_APOBEC_plus"]}, compute_denovo=False)["SBS2_APOBEC_plus"]["eligible_site_frac"]
        oes.append(obs_g / float(p1_g[3] * p1_g[1] * (p1_g[0] + p1_g[3])))
    assert abs(oes[0] - 1.0) < 0.08 and abs(oes[1] - 1.0) < 0.08, oes
    ok["thmV10_2_tcw_oe_gc_invariance_ok"] = True

    # 8. THM-V3.3 exact identity: chi_mm_denovo_unb(s) == chi_mm_unb(s || gap || rc) to machine precision
    s_pr2_dn = np.concatenate([s_asym, np.full(20, 255, dtype=np.uint8), rc_as])
    mj_pr2_dn = analyze_mutagenesis_signatures(s_pr2_dn, {"SBS2_APOBEC_plus": COSMIC_SIG_TABLES["SBS2_APOBEC_plus"]}, compute_denovo=False)["SBS2_APOBEC_plus"]
    dn_eq_err = abs(ap_as["chi_mm_denovo_unb_x1e4"] - mj_pr2_dn["chi_mm_unb_x1e4"])
    assert dn_eq_err < 1e-12, (dn_eq_err, ap_as["chi_mm_denovo_unb_x1e4"], mj_pr2_dn["chi_mm_unb_x1e4"])
    ok["thmV3_3_denovo_unb_equals_pr2_duplex_err"] = dn_eq_err

    # 9. THM-V3.3 exact Parallelogram Law: 0.5*(||J(s)_--||^2 + ||J(rc)_--||^2) == ||J_dn_--||^2 + ||J_sub_--||^2
    mj_rc_as = analyze_mutagenesis_signatures(rc_as, {"SBS2_APOBEC_plus": COSMIC_SIG_TABLES["SBS2_APOBEC_plus"]}, compute_denovo=False)["SBS2_APOBEC_plus"]
    lhs_p = 0.5 * (ap_as["chi_mm_x1e4"] + mj_rc_as["chi_mm_x1e4"])
    rhs_p = ap_as["chi_mm_denovo_x1e4"] + ap_as["chi_mm_sub_x1e4"]
    par_err = abs(lhs_p - rhs_p)
    assert par_err < 1e-13, (par_err, lhs_p, rhs_p)
    ok["thmV3_3_parallelogram_law_err"] = par_err

    # 10. THM-V11.2: Stationary Markov-2 (64-trinucleotide) table satisfies I_nonbridge(P_M2) == 0 identically!
    P_m2 = exact_markov2_from_trinuc_counts(count_triplets_1d(s_asym))
    sc_m2 = two_sided_screen_full(P_m2 * 1e9)
    assert abs(float(P_m2.sum()) - 1.0) < 1e-14 and sc_m2["I_nonbridge"] < 1e-14, sc_m2
    ok["thmV11_2_markov2_nonbridge_zero_err"] = sc_m2["I_nonbridge"]

    # 11. THM-V12 exact conservation: sum(C_m) == 0 for every degree m = 1..6
    m_cons_err = max(abs(float(C_mini[m_].sum())) for m_ in range(1, 7))
    assert m_cons_err < 1e-13, m_cons_err
    ok["thmV12_mass_conservation_all_degrees_err"] = m_cons_err

    # 12. THM-V13 exact Klein-orbit attribution sums to 100.000%
    W_test_a = klein_v4_pr2_audit(count_triplet_pairs_lag(s_asym, d=3))["_W_mp"]
    W_test_b = klein_v4_pr2_audit(count_triplet_pairs_lag(s_hi_gc, d=3))["_W_mp"]
    att_t = klein_orbit_attribution_v13(W_test_a, W_test_b, top_k=5)
    att_err = abs(att_t["sum_orbit_shares_pct"] - 100.0)
    assert att_err < 1e-11, att_t
    ok["thmV13_orbit_attribution_sum_err"] = att_err
    ok.update(self_test_virome_onco_v27(verbose=False))

    if verbose:
        print(json.dumps(_clean(ok), indent=2))
    return ok


def self_test_virome_onco_v27(verbose=False):
    """[THM-V14..V17] 6 exact Zero-RNG algebraic proofs for v2.7.0 (83/83 total)."""
    ok = {}

    def _p_mm(M):
        return 0.25 * (M - M.T - M.T[np.ix_(MIRROR, MIRROR)] + M[np.ix_(MIRROR, MIRROR)])

    def _p_mp(M):
        return 0.25 * (M - M.T + M.T[np.ix_(MIRROR, MIRROR)] - M[np.ix_(MIRROR, MIRROR)])

    # 1. THM-V14a: All-orders (m=1..6) de novo chiral null for symmetric SBS1 on PR2-broken DNA
    _ix = np.arange(800, dtype=np.uint64) + np.uint64(31415)
    _zz = (_ix ^ (_ix >> np.uint64(30))) * np.uint64(0xbf58476d1ce4e5b9)
    _zz = (_zz ^ (_zz >> np.uint64(27))) * np.uint64(0x94d049bb133111eb)
    _u = ((_zz ^ (_zz >> np.uint64(31))) % np.uint64(10)).astype(np.uint8)
    s_asym = np.where(_u < 5, 2, np.where(_u < 7, 1, np.where(_u < 9, 3, 0))).astype(np.uint8)
    rc_as = reverse_complement_nt(s_asym)
    C_f = exact_finite_dose_polynomial_v12(s_asym, COSMIC_SIG_TABLES["SBS1_CpG_sym"]["rates"])
    C_r = exact_finite_dose_polynomial_v12(rc_as, COSMIC_SIG_TABLES["SBS1_CpG_sym"]["rates"])
    all_ord_sym_err = max(float(np.sum(_p_mm(C_f[m_] + C_r[m_]) ** 2)) for m_ in range(1, 7))
    assert all_ord_sym_err < 1e-26, all_ord_sym_err
    ok["thmV14_all_orders_symmetric_denovo_zero_err"] = all_ord_sym_err

    # 2. THM-V14b: Exact universal isolated-CpG double-hit angle == -sqrt(10/17)
    s_cpg = np.array([NT_LUT[ord(ch)] for ch in "AGACTGAATCGAGGCCATAA"], dtype=np.uint8)
    C_cpg = exact_finite_dose_polynomial_v12(s_cpg, COSMIC_SIG_TABLES["SBS1_CpG_sym"]["rates"])
    cos_cpg = _cos_mat_frob(_p_mm(C_cpg[1]), _p_mm(2.0 * C_cpg[2]))
    cpg_angle_err = abs(cos_cpg - (-math.sqrt(10.0 / 17.0)))
    assert cpg_angle_err < 1e-14, (cos_cpg, -math.sqrt(10.0 / 17.0), cpg_angle_err)
    ok["thmV14_isolated_cpg_angle_err"] = cpg_angle_err

    # 3. THM-V14c: Isolated TCW motif has identically ZERO APOBEC3 double-hit Hessian
    s_tcw = np.array([NT_LUT[ord(ch)] for ch in "AGAGGAGATCAGGAGGAGAA"], dtype=np.uint8)
    C_tcw = exact_finite_dose_polynomial_v12(s_tcw, COSMIC_SIG_TABLES["SBS2_APOBEC_plus"]["rates"])
    tcw_h_err = float(np.sum(C_tcw[2] ** 2))
    assert tcw_h_err == 0.0 and float(np.sum(_p_mm(C_tcw[1]) ** 2)) > 0.0, tcw_h_err
    ok["thmV14_isolated_tcw_zero_hessian_err"] = tcw_h_err

    # 4. THM-V15: Universal 705/256 single-base duplex barrier + exact local-vs-global ORF identity
    s_11 = np.array([2, 1, 0, 3, 2, 0, 1, 3, 1, 3, 3], dtype=np.uint8)  # GCATG[A->G]CTCTT
    s_11m = s_11.copy()
    s_11m[5] = 2
    dB_11 = count_triplet_pairs_lag(s_11m, d=3) - count_triplet_pairs_lag(s_11, d=3)
    dB_mp = _p_mp(dB_11)
    dB_3600 = dB_mp - pi_env_project(dB_mp)
    barrier_err = abs(float(np.sum(dB_3600 ** 2)) - (705.0 / 256.0))
    assert barrier_err < 1e-13, (float(np.sum(dB_3600 ** 2)), 705.0 / 256.0)
    # Local 11-nt slice == full-ORF difference identically
    orf_chk = np.concatenate([np.array([0, 3, 2], dtype=np.uint8),
                              np.tile(np.array([2, 1, 0, 3, 2, 0, 1, 3, 1], dtype=np.uint8), 30),
                              np.array([3, 0, 0], dtype=np.uint8)])
    p0_t = 3 * 12
    orf_mut = orf_chk.copy()
    orf_mut[p0_t + 2] = (int(orf_mut[p0_t + 2]) + 1) % 4
    dB_full = count_triplet_pairs_lag(orf_mut, d=3) - count_triplet_pairs_lag(orf_chk, d=3)
    dB_loc = count_triplet_pairs_lag(orf_mut[p0_t - 3:p0_t + 8], d=3) - count_triplet_pairs_lag(orf_chk[p0_t - 3:p0_t + 8], d=3)
    loc_glob_err = float(np.max(np.abs(dB_full - dB_loc)))
    assert loc_glob_err == 0.0, loc_glob_err
    ok["thmV15_universal_705_256_barrier_err"] = barrier_err
    ok["thmV15_local_11nt_vs_full_orf_err"] = loc_glob_err

    # 5. THM-V16: Exact unbiased U-statistic 496+3600 traces (130, 894) & Pythagorean sum == ||W_-+||^2_unb
    assert abs(float(_D_ENV_WMP_DIAG.sum()) - 130.0) < 1e-12
    assert abs(float(_D_NONENV_WMP_DIAG.sum()) - 894.0) < 1e-12
    B_as = count_triplet_pairs_lag(s_asym, d=3)
    kv_as = klein_v4_pr2_audit(B_as)
    a11_as = thm_D11_arrow_screen_and_cpg_audit(B_as)
    pyth_u_err = abs((a11_as["norm2_W_mp_env_496_unb"] + a11_as["norm2_W_mp_nonenv_3600_unb"]) - kv_as["norm2_W_mp_unb"])
    assert pyth_u_err < 1e-16, pyth_u_err
    ok["thmV16_unbiased_496_3600_pythag_err"] = pyth_u_err

    # 6. THM-V17: Exact duplex chiral interference law & Cauchy-Schwarz strand-selection bound
    _, Jf, _ = exact_signature_jacobian(s_asym, COSMIC_SIG_TABLES["SBS2_APOBEC_plus"]["rates"])
    _, Jr, _ = exact_signature_jacobian(rc_as, COSMIC_SIG_TABLES["SBS2_APOBEC_plus"]["rates"])
    Jmm_f = _p_mm(Jf)
    Jmm_r = _p_mm(Jr)
    Jmm_dn = 0.5 * (Jmm_f + Jmm_r)
    Jmm_sb = 0.5 * (Jmm_f - Jmm_r)
    lhs_i = float(np.sum(Jmm_f ** 2) - np.sum(Jmm_r ** 2))
    rhs_i = 4.0 * float(np.sum(Jmm_dn * Jmm_sb))
    interf_err = abs(lhs_i - rhs_i) / max(float(np.sum(Jmm_f ** 2) + np.sum(Jmm_r ** 2)), 1.0)
    f_dn = float(np.sum(Jmm_dn ** 2)) / max(float(np.sum(Jmm_dn ** 2) + np.sum(Jmm_sb ** 2)), 1e-30)
    assert interf_err < 1e-14 and abs(lhs_i) / (float(np.sum(Jmm_f ** 2) + np.sum(Jmm_r ** 2))) <= 2.0 * math.sqrt(f_dn * (1.0 - f_dn)) + 1e-12
    ok["thmV17_chiral_interference_identity_err"] = interf_err
    ok.update(self_test_virome_onco_v28(verbose=False))

    if verbose:
        print(json.dumps(_clean(ok), indent=2))
    return ok


def self_test_virome_onco_v28(verbose=False):
    # [THM-V18] 10 exact Zero-RNG algebraic & statistical self-tests for Patch V18 (94/94 total).
    from itertools import product as _iter_prod
    ok = {}

    # 1. Exact degenerate U-statistic null variance vs exhaustive enumeration on 4x4 Klein model (N=3)
    K_t = 4
    mir_t = np.array([3, 2, 1, 0])
    A_t = np.array([
        [0.10, 0.05, 0.08, 0.03],
        [0.05, 0.12, 0.04, 0.02],
        [0.08, 0.04, 0.12, 0.05],
        [0.03, 0.02, 0.05, 0.10],
    ], dtype=np.float64)
    P_t = 0.25 * (A_t + A_t.T + A_t.T[np.ix_(mir_t, mir_t)] + A_t[np.ix_(mir_t, mir_t)])
    P_t /= P_t.sum()
    pf = P_t.ravel()
    Em = np.eye(K_t * K_t).reshape(K_t * K_t, K_t, K_t)
    H_mp = np.zeros((K_t * K_t, K_t * K_t))
    H_mm = np.zeros((K_t * K_t, K_t * K_t))
    for i_ in range(K_t * K_t):
        M_i = Em[i_]
        Pmp_i = 0.25 * (M_i - M_i.T + M_i.T[np.ix_(mir_t, mir_t)] - M_i[np.ix_(mir_t, mir_t)])
        Pmm_i = 0.25 * (M_i - M_i.T - M_i.T[np.ix_(mir_t, mir_t)] + M_i[np.ix_(mir_t, mir_t)])
        H_mp[i_] = Pmp_i.ravel()
        H_mm[i_] = Pmm_i.ravel()
    Kmp = H_mp @ H_mp.T
    Kmm = H_mm @ H_mm.T
    var_mp_ex = var_mm_ex = 0.0
    for x1, x2, x3 in _iter_prod(range(K_t * K_t), repeat=3):
        pr = pf[x1] * pf[x2] * pf[x3]
        ump = (Kmp[x1, x2] + Kmp[x1, x3] + Kmp[x2, x3]) / 3.0
        umm = (Kmm[x1, x2] + Kmm[x1, x3] + Kmm[x2, x3]) / 3.0
        var_mp_ex += pr * ump * ump
        var_mm_ex += pr * umm * umm
    p2_t = float(np.sum(P_t * P_t))
    d2_t = float(np.sum(np.diag(P_t) ** 2))
    l2_t = float(np.sum(P_t[np.arange(K_t), mir_t] ** 2))
    cl_mp = (p2_t - d2_t + l2_t) / (2.0 * 3.0 * 2.0)
    cl_mm = (p2_t - d2_t - l2_t) / (2.0 * 3.0 * 2.0)
    var_err = max(abs(var_mp_ex - cl_mp), abs(var_mm_ex - cl_mm))
    assert var_err < 1e-15, (var_err, var_mp_ex, cl_mp, var_mm_ex, cl_mm)
    ok["thmV18_null_var_exact_vs_enum_err"] = var_err

    # 2. De-attenuated cosine gate (noise -> unresolved; clean signal -> cos_unb >= cos_raw)
    _ix = np.arange(8000, dtype=np.uint64) + np.uint64(271828)
    _zz = (_ix ^ (_ix >> np.uint64(30))) * np.uint64(0xbf58476d1ce4e5b9)
    _zz = (_zz ^ (_zz >> np.uint64(27))) * np.uint64(0x94d049bb133111eb)
    s_rnd = ((_zz ^ (_zz >> np.uint64(31))) & np.uint64(3)).astype(np.uint8)
    B_sym_noise = count_triplet_pairs_lag(s_rnd, d=3)
    B_sym_noise = 0.5 * (B_sym_noise + B_sym_noise.T)
    s_sig = np.tile(np.array([0, 0, 1, 2, 3, 1, 2, 0, 3, 2, 1], dtype=np.uint8), 500)
    B_sig = count_triplet_pairs_lag(s_sig, d=3)
    W_sig = klein_v4_pr2_audit(B_sig)["_W_mp"]
    cu_noise = _cos_wmp_unbiased_vs_ref(B_sym_noise, W_sig, z_gate=3.0)
    cu_sig = _cos_wmp_unbiased_vs_ref(B_sig, W_sig, z_gate=3.0)
    assert cu_noise["status"] == "unresolved" and cu_noise["cos_unb"] is None, cu_noise
    assert cu_sig["status"] == "resolved" and cu_sig["cos_unb"] >= cu_sig["cos_raw"], cu_sig
    ok["thmV18_cos_unb_recovers_pop_and_gates_noise_ok"] = True

    # 3. Exact Pythagorean cosine identity on 496 + 3600 orthogonal decomposition
    s_b = np.tile(np.array([2, 1, 0, 3, 2, 0, 1, 3, 1, 0, 2, 2, 3], dtype=np.uint8), 400)
    Wa = klein_v4_pr2_audit(B_sig)["_W_mp"]
    Wb = klein_v4_pr2_audit(count_triplet_pairs_lag(s_b, d=3))["_W_mp"]
    Wa_e, Wa_n = pi_env_project(Wa), Wa - pi_env_project(Wa)
    Wb_e, Wb_n = pi_env_project(Wb), Wb - pi_env_project(Wb)
    sa = float(np.sum(Wa_e ** 2) / np.sum(Wa ** 2))
    sb = float(np.sum(Wb_e ** 2) / np.sum(Wb ** 2))
    lhs_c = _cos_mat_frob(Wa, Wb)
    rhs_c = math.sqrt(sa * sb) * _cos_mat_frob(Wa_e, Wb_e) + math.sqrt((1.0 - sa) * (1.0 - sb)) * _cos_mat_frob(Wa_n, Wb_n)
    cos_pyth_err = abs(lhs_c - rhs_c)
    assert cos_pyth_err < 1e-14, (cos_pyth_err, lhs_c, rhs_c)
    ok["thmV18_cos_496_3600_pythagorean_identity_err"] = cos_pyth_err

    # 4. Exact Benjamini-Hochberg FDR check on known vector [0.001, 0.04, 0.03, 0.5]
    q_bh = benjamini_hochberg_fdr([0.001, 0.04, 0.03, 0.5])
    q_exp = [0.004, 0.05333333333333334, 0.05333333333333334, 0.5]
    bh_err = max(abs(a - b) for a, b in zip(q_bh, q_exp))
    assert bh_err < 1e-15, (bh_err, q_bh, q_exp)
    ok["thmV18_bh_fdr_exact_values_err"] = bh_err

    # 5. Exact 3-term concatenation energy identity in deterministic_proviral_insertion_audit
    ins_18 = deterministic_proviral_insertion_audit(s_rnd[:3000], s_sig[:700])
    e3_err = max(ins_18["concat_energy_3term_W_mm"]["decomp_err"], ins_18["concat_energy_3term_W_mp"]["decomp_err"])
    assert e3_err < 1e-11, ins_18
    ok["thmV18_concat_3term_energy_identity_err"] = e3_err

    # 6. Exact chiral excess over 1/2 identity: chiral_share == 1/2 - <A, RA> / (2||A||^2)
    mj_18 = analyze_mutagenesis_signatures(s_rnd[:3000], {
        "SBS2_APOBEC_plus": COSMIC_SIG_TABLES["SBS2_APOBEC_plus"],
        "SBS13_APOBEC_plus": COSMIC_SIG_TABLES["SBS13_APOBEC_plus"],
    }, compute_denovo=False)
    ch_id_err = max(
        abs(mj_18[sn]["J_chiral_share_of_directed"] - (0.5 + mj_18[sn]["J_chiral_excess_over_half"]))
        for sn in mj_18
    )
    assert ch_id_err < 1e-15, ch_id_err
    ok["thmV18_chiral_excess_over_half_identity_err"] = ch_id_err

    # 7. Three exact chiral regimes on uniform PR2-balanced duplex: SBS1==0, SBS2~0.5, self-mirror SBS13>0.60
    s_u_pr2 = np.concatenate([s_rnd, np.full(20, 255, dtype=np.uint8), reverse_complement_nt(s_rnd)])
    mj_3r = analyze_mutagenesis_signatures(s_u_pr2, {
        "SBS1_CpG_sym": COSMIC_SIG_TABLES["SBS1_CpG_sym"],
        "SBS2_APOBEC_plus": COSMIC_SIG_TABLES["SBS2_APOBEC_plus"],
        "SBS13_APOBEC_plus": COSMIC_SIG_TABLES["SBS13_APOBEC_plus"],
    }, compute_denovo=False)
    assert mj_3r["SBS1_CpG_sym"]["J_chiral_share_of_directed"] < 1e-25, mj_3r["SBS1_CpG_sym"]
    assert abs(mj_3r["SBS2_APOBEC_plus"]["J_chiral_share_of_directed"] - 0.5) < 0.08, mj_3r["SBS2_APOBEC_plus"]
    assert mj_3r["SBS13_APOBEC_plus"]["J_chiral_share_of_directed_unb"] > 0.60 and mj_3r["SBS13_APOBEC_plus"]["J_chiral_share_of_directed"] > 0.55 and mj_3r["SBS13_APOBEC_plus"]["J_R_overlap_cos"] < -0.10, mj_3r["SBS13_APOBEC_plus"]
    ok["thmV18_sbs13_self_mirror_exceeds_half_ok"] = True

    # 8. Unit-dose gate for eps_crit > 1.0 on short sparse sequence
    s_short_sparse = s_rnd[:1500]
    mj_sp = analyze_mutagenesis_signatures(s_short_sparse, {"SBS2_APOBEC_plus": COSMIC_SIG_TABLES["SBS2_APOBEC_plus"]}, compute_denovo=False)["SBS2_APOBEC_plus"]
    assert (mj_sp["eps_crit_4N"] is not None) and (mj_sp["eps_crit_4N"] > 1.0) and (mj_sp["eps_crit_status"] == "no_physical_transition"), mj_sp
    ok["thmV18_eps_crit_unit_dose_gate_ok"] = True

    # 9. 3-covariate residualization orthogonality check (gc, log10L, cpg_oe)
    syn_vr = []
    for k_ in range(16):
        gc_k = 0.35 + 0.04 * (k_ % 8)
        ln_k = 4000 + 5000 * ((k_ * 7) % 11)
        cpg_k = 0.15 + 0.06 * ((k_ * 5) % 13)
        grp_k = "herpes_large" if k_ < 8 else "rna_acute_ctrl"
        y_k = 0.3 * gc_k + 0.05 * math.log10(ln_k) - 0.1 * cpg_k + (0.20 if k_ < 8 else 0.0)
        syn_vr.append({"vclass": grp_k, "gc": gc_k, "length_bp": ln_k, "cpg_oe": cpg_k, "y": y_k})
    rw_3 = residualized_group_test_v11(syn_vr, {"herpes_large"}, {"rna_acute_ctrl"}, lambda r: r["y"], include_cpg_oe=True)
    assert rw_3 is not None and rw_3["p_twosided"] < 1e-3 and rw_3["covariates"] == ["gc", "log10L", "cpg_oe"], rw_3
    ok["thmV18_3cov_cpg_residualization_ok"] = True

    # 10. Zero-RNG design-effect audit (stride=6 + block jackknife)
    de_t = klein_design_effect_audit(s_sig, d=3)
    assert de_t is not None and de_t["n_blocks_jk"] >= 5 and de_t["N_disj6"] > 500 and de_t["z_iid_mp_disj6"] > 3.0, de_t
    ok["thmV18_design_effect_jackknife_positive_ok"] = True
    ok.update(self_test_virome_onco_v29(verbose=False))

    if verbose:
        print(json.dumps(_clean(ok), indent=2))
    return ok


def self_test_virome_onco_v29(verbose=False):
    """[THM-V19 & THM-V15b] 3 exact Zero-RNG algebraic self-tests for v2.9.0 (97/97 total)."""
    from fractions import Fraction
    ok = {}

    def _v4(M):
        MT = M.T
        MR = M.T[np.ix_(MIRROR, MIRROR)]
        MS = M[np.ix_(MIRROR, MIRROR)]
        return (
            0.25 * (M + MT + MR + MS),
            0.25 * (M + MT - MR - MS),
            0.25 * (M - MT + MR - MS),
            0.25 * (M - MT - MR + MS),
        )

    # 1. THM-V19a: All 8 Watson-Crick self-conjugate COSMIC channels satisfy R(J) == -J
    # identically on a genuinely non-reversible (D != D.T) stationary PR2 Markov-1 chain.
    t_arr = np.arange(4000, dtype=np.int64)
    s_nr = ((t_arr * 7 + (t_arr // 3) * 5 + (t_arr * t_arr) // 11) % 4).astype(np.uint8)
    rc_nr = (3 - s_nr[::-1]).astype(np.uint8)
    D_pr2 = np.zeros((4, 4), dtype=np.float64)
    for seq_i in (s_nr, rc_nr):
        b_i = np.roll(seq_i, -1)
        D_pr2 += np.bincount((seq_i.astype(np.int32) << 2) | b_i.astype(np.int32), minlength=16).reshape(4, 4)
    D_pr2 /= float(D_pr2.sum())
    assert float(np.max(np.abs(D_pr2 - D_pr2.T))) > 0.01
    p1_pr2 = D_pr2.sum(axis=1)
    T_pr2 = D_pr2 / p1_pr2[:, None]

    idx8 = np.arange(1 << 16, dtype=np.int32)
    b8 = [(idx8 >> (2 * (7 - k_))) & 3 for k_ in range(8)]
    po8 = (b8[1] << 10) | (b8[2] << 8) | (b8[3] << 6) | (b8[4] << 4) | (b8[5] << 2) | b8[6]
    w_m1 = p1_pr2[b8[0]].copy()
    for k_ in range(7):
        w_m1 *= T_pr2[b8[k_], b8[k_ + 1]]

    def _eval_j_enum(w8_vec, rules):
        jf = np.zeros(4096, dtype=np.float64)
        for p_ in range(6):
            sh_ = 2 * (5 - p_)
            cm_ = ~(3 << sh_)
            for (ca, cb, cc), bn in rules:
                m_ = (b8[p_] == ca) & (b8[p_ + 1] == cb) & (b8[p_ + 2] == cc)
                po_ = po8[m_]
                pn_ = (po_ & cm_) | (bn << sh_)
                wm_ = w8_vec[m_]
                jf += np.bincount(pn_, weights=wm_, minlength=4096)
                jf -= np.bincount(po_, weights=wm_, minlength=4096)
        return jf.reshape(NC, NC)

    max_anti_err = 0.0
    for cb in (1, 3):
        bn = 3 - cb
        for ca in range(4):
            cc = 3 - ca
            J_sc = _eval_j_enum(w_m1, [((ca, cb, cc), bn)])
            RJ_sc = J_sc.T[np.ix_(MIRROR, MIRROR)]
            max_anti_err = max(max_anti_err, float(np.max(np.abs(J_sc + RJ_sc))))
    assert max_anti_err < 1e-16, max_anti_err
    ok["thmV19_self_conjugate_R_antiinv_max_err"] = max_anti_err

    # 2. THM-V19b: Exact 5-level rational spectrum of all 96 COSMIC single-strand channels
    # on uniform 4^8 DNA + exact SBS2 (23/48), SBS13 (7/10), and T-even SBS4 null.
    w_unif = np.full(1 << 16, 1.0 / float(1 << 16), dtype=np.float64)
    level_counts = {}
    for cb in (1, 3):
        for bn in range(4):
            if bn == cb:
                continue
            for ca in range(4):
                for cc in range(4):
                    J_c = _eval_j_enum(w_unif, [((ca, cb, cc), bn)])
                    _, _, Wmp_c, Wmm_c = _v4(J_c)
                    nmp_c, nmm_c = float(np.sum(Wmp_c ** 2)), float(np.sum(Wmm_c ** 2))
                    fr = Fraction(nmm_c / (nmp_c + nmm_c)).limit_denominator(1000)
                    level_counts[fr] = level_counts.get(fr, 0) + 1
    expected_96 = {
        Fraction(5, 12): 16, Fraction(3, 7): 4, Fraction(11, 24): 48,
        Fraction(1, 2): 20, Fraction(1, 1): 8,
    }
    assert level_counts == expected_96, level_counts
    J_s2 = _eval_j_enum(w_unif, [((3, 1, 0), 3), ((3, 1, 3), 3)])
    J_s13 = _eval_j_enum(w_unif, [((3, 1, 0), 2), ((3, 1, 3), 2)])
    J_s4 = _eval_j_enum(w_unif, [((a_, 2, c_), 3) for a_ in range(4) for c_ in range(4)])
    _, _, mp2, mm2 = _v4(J_s2)
    _, _, mp13, mm13 = _v4(J_s13)
    _, _, mp4, mm4 = _v4(J_s4)
    assert abs(float(np.sum(mm2 ** 2)) / float(np.sum(mp2 ** 2) + np.sum(mm2 ** 2)) - 23.0 / 48.0) < 1e-14
    assert abs(float(np.sum(mm13 ** 2)) / float(np.sum(mp13 ** 2) + np.sum(mm13 ** 2)) - 7.0 / 10.0) < 1e-14
    assert float(np.sum(mp4 ** 2) + np.sum(mm4 ** 2)) < 1e-30
    ok["thmV19_cosmic96_five_rational_levels_ok"] = True

    # 3. THM-V15b: Full 4-quadrant exact rational quanta on 496 + 3600 for generic single-base hit
    s_11 = np.array([2, 1, 0, 3, 2, 0, 1, 3, 1, 3, 3], dtype=np.uint8)
    s_11m = s_11.copy()
    s_11m[5] = 2
    dB_11 = count_triplet_pairs_lag(s_11m, d=3) - count_triplet_pairs_lag(s_11, d=3)
    exp_q = [
        (125.0 / 512.0, 1411.0 / 512.0),
        (63.0 / 256.0, 705.0 / 256.0),
        (63.0 / 256.0, 705.0 / 256.0),
        (127.0 / 512.0, 1409.0 / 512.0),
    ]
    q_err = 0.0
    for Q_i, (ee, en) in zip(_v4(dB_11), exp_q):
        Qe_i = pi_env_project(Q_i)
        Qn_i = Q_i - Qe_i
        q_err = max(q_err, abs(float(np.sum(Qe_i ** 2)) - ee), abs(float(np.sum(Qn_i ** 2)) - en))
    assert q_err < 1e-13, q_err
    ok["thmV15_four_quadrant_rational_quanta_err"] = q_err
    ok.update(self_test_virome_onco_v291(verbose=False))

    if verbose:
        print(json.dumps(_clean(ok), indent=2))
    return ok


def self_test_virome_onco_v291(verbose=False):
    """[THM-V20] 3 exact Zero-RNG algebraic self-tests for v2.9.1 (100/100 total):
    (1) Per-orbit chirality rule across all 2016 T-odd basis vectors & chiral blindness of all 64 RC palindromes;
    (2) Isolated CpG exact 1/2 U-statistic share vs isolated TCW zero U-statistic;
    (3) Closed-form sparse-tile cosine law -6*sqrt(2/(17*N)) (exact) & -sqrt(6p/(2+5p))."""
    from fractions import Fraction
    ok = {}

    def _p_pm(M):
        return 0.25 * (M + M.T - M.T[np.ix_(MIRROR, MIRROR)] - M[np.ix_(MIRROR, MIRROR)])

    def _p_mp(M):
        return 0.25 * (M - M.T + M.T[np.ix_(MIRROR, MIRROR)] - M[np.ix_(MIRROR, MIRROR)])

    def _p_mm(M):
        return 0.25 * (M - M.T - M.T[np.ix_(MIRROR, MIRROR)] + M[np.ix_(MIRROR, MIRROR)])

    # 1. THM-V20a: Exact per-orbit rule on all 2016 T-odd basis vectors + 64 RC-palindromic 6-mers
    iu, iv = np.triu_indices(NC, k=1)
    pal = (iv == MIRROR[iu])
    assert int(pal.sum()) == 32
    E_mp_exp = 0.5 * (1.0 + pal.astype(np.float64))
    E_mm_exp = 0.5 * (1.0 - pal.astype(np.float64))
    assert abs(float(E_mp_exp.sum()) - 1024.0) < 1e-12
    assert abs(float(E_mm_exp.sum()) - 992.0) < 1e-12
    inv_s2 = 1.0 / math.sqrt(2.0)
    orb_err = 0.0
    for k_ in range(iu.size):
        u_, v_ = int(iu[k_]), int(iv[k_])
        e_ = np.zeros((NC, NC), dtype=np.float64)
        e_[u_, v_] = inv_s2
        e_[v_, u_] = -inv_s2
        orb_err = max(
            orb_err,
            abs(float(np.sum(_p_mp(e_) ** 2)) - E_mp_exp[k_]),
            abs(float(np.sum(_p_mm(e_) ** 2)) - E_mm_exp[k_]),
        )
    pal_err = 0.0
    for u_ in range(NC):
        v_ = int(MIRROR[u_])
        E_pal = np.zeros((NC, NC), dtype=np.float64)
        E_pal[u_, v_] = 1.0
        pal_err = max(pal_err, float(np.max(np.abs(_p_pm(E_pal)))), float(np.max(np.abs(_p_mm(E_pal)))))
    for site_ in ("GAATTC", "GGATCC", "AAGCTT", "CTCGAG", "GTCGAC", "CCGCGG"):
        assert int(MIRROR[TRIPLET_IDX[site_[:3]]]) == TRIPLET_IDX[site_[3:]]
    assert Fraction(992, 2016) == Fraction(31, 63)
    assert Fraction(894, 1024) == Fraction(447, 512)
    assert Fraction(32, 1024) == Fraction(1, 32)
    err_v20a = max(orb_err, pal_err)
    assert err_v20a < 1e-14, err_v20a
    ok["thmV20_palindrome_chiral_blindness_and_orbit_rule_err"] = err_v20a

    # 2. THM-V20b: Isolated CpG -> exact U-statistic cross-term 5/2 in both W_-+ & W_-- => share == 1/2;
    #    isolated TCW (single site) -> exact U-statistic u_mp == u_mm == 0.0
    s_cpg = np.array([NT_LUT[ord(ch)] for ch in "AGACTGAATCGAGGCCATAA"], dtype=np.uint8)
    B0_c, J_c, _ = exact_signature_jacobian(s_cpg, COSMIC_SIG_TABLES["SBS1_CpG_sym"]["rates"])
    dg_c = exact_signature_jacobian_diag(s_cpg, COSMIC_SIG_TABLES["SBS1_CpG_sym"]["rates"])
    N_c = float(B0_c.sum())
    nmp_c = float(np.sum(_p_mp(J_c) ** 2))
    nmm_c = float(np.sum(_p_mm(J_c) ** 2))
    ump_c = (nmp_c - dg_c["W_mp"]) / (N_c * (N_c - 1.0))
    umm_c = (nmm_c - dg_c["W_mm"]) / (N_c * (N_c - 1.0))
    sh_c = umm_c / (ump_c + umm_c)

    s_tcw = np.array([NT_LUT[ord(ch)] for ch in "AGAGGAGATCAGGAGGAGAA"], dtype=np.uint8)
    B0_t, J_t, _ = exact_signature_jacobian(s_tcw, COSMIC_SIG_TABLES["SBS2_APOBEC_plus"]["rates"])
    dg_t = exact_signature_jacobian_diag(s_tcw, COSMIC_SIG_TABLES["SBS2_APOBEC_plus"]["rates"])
    ump_t = float(np.sum(_p_mp(J_t) ** 2)) - dg_t["W_mp"]
    umm_t = float(np.sum(_p_mm(J_t) ** 2)) - dg_t["W_mm"]
    err_v20b = max(
        abs(sh_c - 0.5),
        abs((nmp_c - dg_c["W_mp"]) - 2.5),
        abs((nmm_c - dg_c["W_mm"]) - 2.5),
        abs(ump_t),
        abs(umm_t),
    )
    assert err_v20b < 1e-14, (err_v20b, sh_c, nmp_c - dg_c["W_mp"], nmm_c - dg_c["W_mm"], ump_t, umm_t)
    ok["thmV20_isolated_cpg_exact_half_ustat_err"] = err_v20b

    # 3. THM-V20c: Exact sparse-tile cosine law -6*sqrt(2/(17*N)) on isolated CpG
    #    and Binomial sparsity law -sqrt(6p/(2+5p)) across 16 deterministic SplitMix64 1-kb tiles
    cos_cpg_emp = _cos_mat_frob(_p_mp(B0_c), _p_mp(J_c))
    cos_cpg_ex = -6.0 * math.sqrt(2.0 / (17.0 * N_c))
    err_v20c = abs(cos_cpg_emp - cos_cpg_ex)
    assert err_v20c < 1e-14, (cos_cpg_emp, cos_cpg_ex, err_v20c)
    for sn_, p_th in (("SBS2_APOBEC_plus", 1.0 / 32.0), ("SBS7a_UV_plus", 1.0 / 8.0), ("SBS4_Tobacco_plus", 1.0 / 4.0)):
        rt_ = COSMIC_SIG_TABLES[sn_]["rates"]
        c_list = []
        for tk_ in range(16):
            ix_ = np.arange(1000, dtype=np.uint64) + np.uint64(100000 + tk_ * 1009)
            zz_ = (ix_ ^ (ix_ >> np.uint64(30))) * np.uint64(0xbf58476d1ce4e5b9)
            zz_ = (zz_ ^ (zz_ >> np.uint64(27))) * np.uint64(0x94d049bb133111eb)
            nt_k = ((zz_ ^ (zz_ >> np.uint64(31))) & np.uint64(3)).astype(np.uint8)
            B0_k, J_k, _ = exact_signature_jacobian(nt_k, rt_)
            c_list.append(_cos_mat_frob(_p_mp(B0_k), _p_mp(J_k)))
        th_c = -math.sqrt(6.0 * p_th / (2.0 + 5.0 * p_th))
        assert abs(float(np.mean(c_list)) - th_c) < 0.04, (sn_, float(np.mean(c_list)), th_c)
    ok["thmV20_sparse_tile_cos_formula_err"] = err_v20c

    if verbose:
        print(json.dumps(_clean(ok), indent=2))
    return ok



def append_v23_markdown_and_figures(payload, out_dir, md_path, dpi=300):
    v10 = payload.get("viral_apobec_restriction_V10") or {}
    rows = v10.get("rows") or []
    vr = payload.get("viruses") or []
    fd_v12 = payload.get("finite_dose_curves_V12") or {}
    extra_md = []

    if rows:
        extra_md += [
            "", "## 10. Theorems V10 & V10.2 — GC- & Dinucleotide-Controlled Viral APOBEC3 TCW Evasion (33 Viruses)", "",
            "| Virus | Class | GC % | TCW % | TCW O/E (Base) | TCW O/E (Dinuc) | CpG O/E | SBS2(+) Unb Chiral % | SBS2(+) De Novo χ_--^unb (×10⁴) |",
            "|---|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for r in rows:
            extra_md.append(
                f"| **{r['virus']}** | {r['vclass']} | {100*r['gc']:.1f}% | "
                f"{fnum(100.0*_fval(r.get('tcw_eligible_frac'), np.nan),2)}% | "
                f"**{fnum(r.get('tcw_oe_base'),3)}** | {fnum(r.get('tcw_oe_dinuc'),3)} | "
                f"{fnum(r.get('cpg_oe'),3)} | "
                f"{fnum(100.0*_fval(r.get('sbs2_chiral_share_unb'), np.nan),2)}% | "
                f"{fnum(r.get('sbs2_chi_mm_denovo_x1e4'),3)} |"
            )
        for lbl_mw, key_mw in (
            ("TCW O/E Base (Integrating/Onco vs Acute RNA)", "mw_tcw_oe_base_int_vs_rna"),
            ("TCW O/E Dinuc (Integrating/Onco vs Acute RNA)", "mw_tcw_oe_dinuc_int_vs_rna"),
            ("TCW O/E Base (Retroviruses+Papilloma/Polyoma vs Acute RNA)", "mw_tcw_oe_base_retro_pap_vs_rna"),
        ):
            mw = v10.get(key_mw)
            if mw:
                extra_md.append(
                    f"> **Mann–Whitney {lbl_mw}:** median `{mw['median1']:.3f}` vs `{mw['median2']:.3f}`, "
                    f"`r_rb = {mw['rank_biserial_r']:+.3f}`, `p = {_fsci(mw['p_twosided'])}`."
                )

    v11_vr = [r for r in vr if r.get("beyond_dinuc_V11")]
    if v11_vr:
        extra_md += [
            "", "## 11. Theorems V11 & V11.2 — Beyond-Dinucleotide (Markov-1) & Beyond-Trinucleotide (Markov-2) Viral Mimicry Proof", "",
            "| Virus | Class | GC % | cos Total (1024-D) | cos 3600-D | cos CpG-Free | cos 3600-D CpG-Free | cos Beyond-M1 (Dinuc) | cos Beyond-M2 (Trinuc) | cos(Self, M2) |",
            "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for r in sorted(v11_vr, key=lambda x: _fval(x["beyond_dinuc_V11"].get("cos_markov2_resid_vs_human"), -2.0), reverse=True):
            b = r["beyond_dinuc_V11"]
            extra_md.append(
                f"| **{r['label']}** | {r['vclass']} | {100*r['gc']:.1f}% | "
                f"{fnum(b.get('cos_total_1024_vs_human'),3,True)} | "
                f"{fnum(b.get('cos_nonenv3600_vs_human'),3,True)} | "
                f"{fnum(b.get('cos_cpg_free_vs_human'),3,True)} | "
                f"{fnum(b.get('cos_nonenv3600_cpg_free_vs_human'),3,True)} | "
                f"{fnum(b.get('cos_markov1_resid_vs_human'),3,True)} | "
                f"**{fnum(b.get('cos_markov2_resid_vs_human'),3,True)}** | "
                f"{fnum(b.get('cos_self_vs_markov2'),3,True)} |"
            )
        sp_map = payload.get("spearman_controls_V11") or {}
        for sk, slbl in (
            ("gc_vs_cos_total", "GC % vs Total 1024-D Alignment"),
            ("gc_vs_cos_beyond_m2", "GC % vs Beyond-Markov-2 Alignment"),
            ("tcw_oe_vs_apobec_denovo", "TCW O/E vs De Novo APOBEC3 χ_--^unb"),
        ):
            sp = sp_map.get(sk)
            if sp:
                extra_md.append(
                    f"> **Spearman Control ({slbl}, n={sp['n']}):** `rho_S = {sp['rho']:+.3f}`, `p = {_fsci(sp['p_twosided'])}`."
                )
        sp_l = sp_map.get("log10L_vs_cos_beyond_m2")
        if sp_l:
            extra_md.append(
                f"> **Spearman Control (log10 Length vs Beyond-Markov-2 Alignment, n={sp_l['n']}):** `rho_S = {sp_l['rho']:+.3f}`, `p = {_fsci(sp_l['p_twosided'])}`."
            )
        for ak, albl in (
            ("beyond_m2_int_vs_rna_adj_gc_L", "Beyond-M2 Excess Onco/Int vs Acute RNA (GC + log10L Residualized)"),
            ("beyond_m2_nucDNA_vs_cytoRNA_adj_gc_L", "Beyond-M2 Excess Nuclear dsDNA vs Acute RNA (GC + log10L Residualized)"),
            ("beyond_m2_gs_nucDNA_vs_cytoRNA_adj_gc_L", "Beyond-M2 Gram-Schmidt Nuclear dsDNA vs Acute RNA (GC + log10L Residualized)"),
        ):
            ar = sp_map.get(ak)
            if ar:
                extra_md.append(
                    f"> **Covariate-Adjusted Control ({albl}, n={ar['n1']}/{ar['n2']}):** "
                    f"resid median `{ar['median1']:+.4f}` vs `{ar['median2']:+.4f}`, "
                    f"`r_rb = {ar['rank_biserial_r']:+.3f}`, `p = {_fsci(ar['p_twosided'])}`."
                )

    if fd_v12:
        extra_md += [
            "", "## 12. Theorem V12 — Exact All-Orders Finite-Dose Mutagenesis Phase Transition & Double-Hit Kataegis Hessian", "",
            "| Target | Signature | Hessian ‖H_--‖²×10⁴ | cos(J_--, H_--) | PR2-ok (ε=0) | PR2-ok (ε=0.10) | PR2-ok (ε=0.50) | PR2-ok (ε=1.00) | Nonlinear Dev (ε=0.25) |",
            "|---|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for tname, sig_map in fd_v12.items():
            for sn in ("SBS1_CpG_sym", "SBS2_APOBEC_plus", "SBS13_APOBEC_plus"):
                sd = sig_map.get(sn)
                if not sd:
                    continue
                by_e = {round(pt["eps"], 3): pt for pt in sd["curve"]}
                def _pr(ev):
                    p = by_e.get(ev)
                    return f"{fnum(100.0*_fval(p.get('PR2_ok_expect', p.get('PR2_ok_unb')), np.nan),2)}%" if p else "n/a"
                nl25 = (by_e.get(0.25) or {}).get("nonlinear_Wmm_rel_dev")
                extra_md.append(
                    f"| **{tname}** | `{sn}` | {fnum(sd.get('hessian_W_mm_norm2_x1e4'),3)} | "
                    f"{fnum(sd.get('cos_J_mm_vs_H_mm'),3,True)} | {_pr(0.0)} | {_pr(0.1)} | {_pr(0.5)} | {_pr(1.0)} | "
                    f"{fnum(100.0*_fval(nl25, np.nan),1)}% |"
                )

    v13 = payload.get("orbit_attribution_V13") or {}
    if v13:
        extra_md += [
            "", "## 13. Theorem V13 — Exact Klein-Orbit 6-mer Motif Attribution (Viral Mimicry & APOBEC3 Kataegis)", "",
        ]
        for idx_a, (title_a, key_a) in enumerate((
            ("KSHV_HHV8 vs. Human Euchromatin (Total 1024-D W_-+)", "KSHV_HHV8_vs_human_total"),
            ("KSHV_HHV8 vs. Human Euchromatin (Beyond-Markov-2 Residual W_-+^⊥M2)", "KSHV_HHV8_vs_human_beyond_M2"),
            ("EBV_HHV4 vs. Human Euchromatin (Beyond-Markov-2 Residual W_-+^⊥M2)", "EBV_HHV4_vs_human_beyond_M2"),
            ("TP53 APOBEC3 SBS2(+) Chiral Kataegis Injection (J_-- Self-Energy)", "TP53_APOBEC_SBS2_Jmm"),
            ("PIK3CA APOBEC3 SBS2(+) Chiral Kataegis Injection (J_-- Self-Energy)", "PIK3CA_APOBEC_SBS2_Jmm"),
        )):
            att = v13.get(key_a)
            if not att or not att.get("top_orbits"):
                continue
            extra_md += [
                f"### 13{chr(65 + idx_a)}. {title_a}", "",
                "| Rank | 6-mer Motif | Mirror 6-mer | Orbit Size | Share of Total Inner Product % | CpG? |",
                "|---:|---|---|---:|---:|:---:|",
            ]
            for rk, ob in enumerate(att["top_orbits"][:8], 1):
                extra_md.append(
                    f"| {rk} | `{ob['motif_6mer']}` | `{ob['mirror_6mer']}` | {ob['orbit_size']} | "
                    f"{ob['share_pct']:+.2f}% | {'Yes' if ob['has_CpG'] else 'No'} |"
                )
            extra_md.append("")

    fdr18 = payload.get("fdr_ledger_V18") or {}
    r1_18 = payload.get("t2t_rank1_V18") or {}
    if fdr18.get("rows") or r1_18:
        extra_md += ["", "## 14. Theorem V18 — Global BH-FDR Ledger, ZAP 3-Covariate Control & T2T Rank-1 Diagnostics", ""]
        if r1_18:
            extra_md.append(
                f"> **T2T Rank-1 Factorization ({r1_18['n_pairs']} pairs):** "
                f"`cos(v, l) ~ cos(v, H)*cos(l, H)` explains **`R² = {r1_18['rank1_r2']:.4f}`** "
                f"(slope `{r1_18['rank1_slope']:.3f}`). Raw top winners: `{r1_18['raw_top_virus_counts']}`; "
                f"Human-residual `W^⊥H` top winners: `{r1_18['residual_top_virus_counts']}`."
            )
        if fdr18.get("rows"):
            extra_md += [
                "",
                f"> **Global Benjamini–Hochberg FDR (`m = {fdr18['m_total_tests']}` tests; `{fdr18['n_primary']}` primary, "
                f"`{fdr18['n_exploratory']}` exploratory; `{fdr18['n_bh_sig_q05']}` survive `q_BH ≤ 0.05`):**", "",
                "| Rank | Tier | Section | Test | Two-Sided p | p_min | q_BH |",
                "|---:|:---:|:---:|---|---:|---:|---:|",
            ]
            for rk_f, fr in enumerate(fdr18["rows"][:15], 1):
                extra_md.append(
                    f"| {rk_f} | `{fr['tier']}` | {fr['section']} | `{fr['test']}` | "
                    f"{_fsci(fr['p_twosided'])} | {_fsci(fr.get('p_min_discrete'))} | {fnum(fr['q_bh_all'], 4)} |"
                )

    if extra_md:
        old_md = Path(md_path).read_text(encoding="utf-8")
        atomic_write(md_path, old_md.rstrip() + "\n" + "\n".join(extra_md) + "\n")

    if not HAVE_MPL:
        return []
    fd = Path(out_dir) / "figures"
    fd.mkdir(parents=True, exist_ok=True)
    saved = []
    cls_col = {
        "retrovirus_int": "#d62728", "dna_onco_small": "#ff7f0e", "dna_benign_ctrl": "#2ca02c",
        "hepatitis_hcc": "#9467bd", "herpes_large": "#1f77b4", "vector_dna": "#17becf",
        "rna_acute_ctrl": "#7f7f7f", "phage_outgroup": "#8c564b",
    }

    # FIG V6 (Upgraded v2.6): GC-Controlled TCW O/E vs De Novo APOBEC3 Chiral Susceptibility
    if rows:
        fig, axes = plt.subplots(1, 2, figsize=(18.5, 7.2))
        fig.suptitle(f"Theorem V10.2: GC-Controlled Viral APOBEC3 TCW-Motif Evasion (O/E) & De Novo Restriction (v{VERSION})",
                     fontsize=13.5, fontweight="bold", y=0.98)
        ax = axes[0]
        ax.axvline(1.0, color="k", ls="--", lw=1.0, alpha=0.6, label="Neutral Composition (TCW O/E = 1.0)")
        rows_plot = sorted(rows, key=lambda r: _fval(r.get("tcw_oe_base"), 0.0))
        pts_v6 = []
        for r in rows_plot:
            oe = _fval(r.get("tcw_oe_base"), np.nan)
            chi = _fval(r.get("sbs2_chi_mm_denovo_x1e4"), np.nan)
            if _finite(oe) and _finite(chi) and chi > 0:
                ax.scatter(oe, chi, s=58, color=cls_col.get(r["vclass"], "#333"), alpha=0.88, edgecolors="k", lw=0.45, zorder=4)
                pts_v6.append((oe, chi, r["virus"]))
        v6_xlim, v6_ylim = (0.41, 1.58), (2.2e-2, 2.2e-1)
        ax.set_yscale("log"); ax.set_xlim(*v6_xlim); ax.set_ylim(*v6_ylim)
        pos_v6 = _layout_scatter_labels([p[0] for p in pts_v6], [p[1] for p in pts_v6], [p[2] for p in pts_v6],
                                        v6_xlim, v6_ylim, log_y=True, min_dx=0.080, min_dy=0.045, iters=240)
        for (px_, py_, lb_), (lx_, ly_, draw_ln) in zip(pts_v6, pos_v6):
            if draw_ln:
                ax.plot([px_, lx_], [py_, ly_], color="#888888", lw=0.55, alpha=0.65, zorder=2)
            ax.text(lx_, ly_, lb_, fontsize=6.5, va="center", zorder=5)
        ax.set_xlabel("GC-Controlled APOBEC3 Target Ratio TCW O/E (Observed / Mononucleotide Expected)")
        ax.set_ylabel("Unbiased De Novo APOBEC3 Chiral Susceptibility χ_--^de novo × 10⁴ (log)")
        ax.set_title("A. Composition-Normalized TCW Evasion vs. De Novo Chiral Response")
        ax.legend(fontsize=8, loc="upper left", framealpha=0.9); ax.grid(True, ls=":", alpha=0.35)

        ax = axes[1]
        rows_s = sorted(rows, key=lambda r: _fval(r.get("tcw_oe_base"), 0.0))
        ys = np.arange(len(rows_s))
        ax.barh(ys, [_fval(r.get("tcw_oe_base"), 0.0) for r in rows_s],
                color=[cls_col.get(r["vclass"], "#333") for r in rows_s], alpha=0.88)
        ax.axvline(1.0, color="k", ls="--", lw=1.1, label="Neutral O/E = 1.0")
        ax.set_yticks(ys); ax.set_yticklabels([f"{r['virus']} ({100*r['gc']:.0f}% GC)" for r in rows_s], fontsize=7.5)
        ax.set_xlabel("GC-Controlled APOBEC3 Target Ratio TCW O/E")
        ax.set_title("B. True Evolutionary TCW Depletion Ranking (GC-Normalized)")
        ax.legend(fontsize=8, loc="lower right", framealpha=0.9); ax.grid(True, axis="x", ls=":", alpha=0.35)

        fig.tight_layout(rect=[0, 0, 1, 0.945])
        p6 = fd / "fig_V6_viral_apobec_restriction.png"
        fig.savefig(p6, dpi=dpi, bbox_inches="tight"); plt.close(fig); saved.append(str(p6))
        print(f"  [FIG] {p6}")

    # FIG V7 (NEW v2.6): Theorem V11 Beyond-Dinucleotide Mimicry & Theorem V12 All-Orders Finite-Dose Curves
    if v11_vr and fd_v12:
        fig, axes = plt.subplots(1, 2, figsize=(18.5, 7.6))
        fig.suptitle(f"Theorems V11–V12: Beyond-Markov-1/2 Viral Mimicry & Exact All-Orders Kataegis Dose Curves (v{VERSION})",
                     fontsize=13.0, fontweight="bold", y=0.98)
        ax = axes[0]
        v11_s = sorted(v11_vr, key=lambda r: _fval(r["beyond_dinuc_V11"].get("cos_markov2_resid_vs_human"), -2.0))
        ys = np.arange(len(v11_s))
        c_tot = [_fval(r["beyond_dinuc_V11"].get("cos_total_1024_vs_human"), 0.0) for r in v11_s]
        c_m1r = [_fval(r["beyond_dinuc_V11"].get("cos_markov1_resid_vs_human"), 0.0) for r in v11_s]
        c_m2r = [_fval(r["beyond_dinuc_V11"].get("cos_markov2_resid_vs_human"), 0.0) for r in v11_s]
        ax.barh(ys + 0.26, c_tot, 0.26, color="#aec7e8", label="Total 1024-D cos(W_-+)")
        ax.barh(ys, c_m1r, 0.26, color="#ff7f0e", alpha=0.88, label="Beyond-Markov-1 (Dinuc) cos(W_-+^⊥M1)")
        ax.barh(ys - 0.26, c_m2r, 0.26, color="#d62728", alpha=0.92, label="Beyond-Markov-2 (64-Trinuc) cos(W_-+^⊥M2)")
        ax.axvline(0.0, color="k", lw=0.9)
        ax.set_yticks(ys); ax.set_yticklabels([r["label"] for r in v11_s], fontsize=7.5)
        ax.set_xlabel("Cosine Alignment to Human T2T Euchromatin")
        ax.set_title("A. THM-V11/V11.2: Viral Mimicry Beyond Dinucleotide (M1) & 64-Trinucleotide (M2) Nulls")
        ax.legend(fontsize=7.8, loc="lower right", framealpha=0.9); ax.grid(True, axis="x", ls=":", alpha=0.35)

        ax = axes[1]
        for tname, col_t in (("FHIT", "#1f77b4"), ("TP53", "#d62728"), ("PIK3CA", "#ff7f0e"), ("TERT", "#9467bd"), ("HPV16", "#2ca02c"), ("HBV", "#8c564b")):
            if tname not in fd_v12:
                continue
            sd_ap = fd_v12[tname].get("SBS2_APOBEC_plus")
            sd_s1 = fd_v12[tname].get("SBS1_CpG_sym")
            if sd_ap:
                ex = [pt["eps"] for pt in sd_ap["curve"]]
                yp = [100.0 * _fval(pt.get("PR2_ok_expect", pt.get("PR2_ok_unb")), np.nan) for pt in sd_ap["curve"]]
                ax.plot(ex, yp, marker="o", ms=4.0, lw=1.8, color=col_t, label=f"{tname}: APOBEC3 SBS2(+)")
            if sd_s1 and tname in ("FHIT", "TERT"):
                ex = [pt["eps"] for pt in sd_s1["curve"]]
                yp = [100.0 * _fval(pt.get("PR2_ok_expect", pt.get("PR2_ok_unb")), np.nan) for pt in sd_s1["curve"]]
                ax.plot(ex, yp, marker="s", ms=3.5, lw=1.3, ls="--", color=col_t, alpha=0.75, label=f"{tname}: SBS1 CpG (sym)")
        ax.set_ylim(37, 103); ax.set_xlim(-0.02, 1.03)
        ax.set_xlabel("Finite Mutagenesis Dose ε (Exact All-Orders Expectation E_ε[B])")
        ax.set_ylabel("Exact Expected PR2-Compliant Directed Flux Share (%)")
        ax.set_title("B. THM-V12: Exact All-Orders Chiral Phase Transition (APOBEC3 vs. SBS1)")
        ax.legend(fontsize=7.5, loc="lower left", ncol=2, framealpha=0.9); ax.grid(True, ls=":", alpha=0.35)

        fig.tight_layout(rect=[0, 0, 1, 0.945])
        p7 = fd / "fig_V7_beyond_dinuc_mimicry_and_dose.png"
        fig.savefig(p7, dpi=dpi, bbox_inches="tight"); plt.close(fig); saved.append(str(p7))
        print(f"  [FIG] {p7}")

    return saved


def run_virome_onco_audit(args):
    """Master driver for Paper 5 (--virome-onco-audit): 100% standalone inside SCRIPT_DIR."""
    import shutil
    t0 = time.time()
    st = self_test(verbose=False)
    out = Path(args.out_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    v_dir = Path(args.viral_dir).resolve()
    o_dir = Path(args.onco_dir).resolve()
    c_dir = Path(args.chr_dir).resolve()
    local_ref_npz = (SCRIPT_DIR / "human_t2t_euch_ref_v2_8.npz").resolve()

    if getattr(args, "fresh_cache", False):
        print("  [FRESH-CACHE] Clearing local viral_cache, onco_cache, chr_cache & reference npz inside SCRIPT_DIR...")
        for d_cl in (v_dir, o_dir, c_dir):
            if d_cl.exists() and str(d_cl).startswith(str(SCRIPT_DIR)):
                shutil.rmtree(d_cl)
        if local_ref_npz.exists():
            local_ref_npz.unlink()
    elif getattr(args, "fresh_virome_cache", False):
        print("  [FRESH-VIROME-CACHE] Clearing local viral_cache & onco_cache inside SCRIPT_DIR...")
        for d_cl in (v_dir, o_dir):
            if d_cl.exists() and str(d_cl).startswith(str(SCRIPT_DIR)):
                shutil.rmtree(d_cl)

    v_dir.mkdir(parents=True, exist_ok=True)
    o_dir.mkdir(parents=True, exist_ok=True)
    c_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 96)
    print(f"PAPER 5 ZERO-RNG VIROME & ONCOGENOMICS AUDIT (v{VERSION}) | 33 Viruses + 28 Human Loci | Standalone SCRIPT_DIR")
    print(f"  SCRIPT_DIR = {SCRIPT_DIR}")
    print("=" * 96)

    human_ref = _load_human_reference_arrows(getattr(args, "atlas_npz", ""), require_all_24=True)
    if (human_ref is None or human_ref.get("n_chr", 0) < 24) and args.auto_download:
        build_local_human_reference(c_dir, local_ref_npz, auto_download=True)
        human_ref = _load_human_reference_arrows(str(local_ref_npz))
    if human_ref:
        print(f"  [REF] Loaded {human_ref['n_chr']}-chr Human T2T Euchromatin Reference from {human_ref['source']}")
    else:
        print("  [REF] No local T2T reference found (pass --auto-download to build inside SCRIPT_DIR); will pool 28 human loci.")

    # 1. Audit 28 Human Cancer & Control Loci (first, so they can also serve as fallback reference & insertion hosts)
    onco_results, onco_nt_cache, onco_failures = [], {}, []
    B_human_fallback = np.zeros((NC, NC), dtype=np.float64)
    # Clean stale wrong FHIT_NG_008074.1 (SBF2) if present from v2.4
    stale_fhit = o_dir / "FHIT_NG_008074.1.fna.gz"
    if stale_fhit.exists():
        stale_fhit.unlink()
        print("  [CLEAN] Removed stale FHIT_NG_008074.1.fna.gz (was SBF2 on chr11)")
    for i, (name, meta) in enumerate(ONCO_LOCI_28.items(), 1):
        meta_use = dict(meta)
        try:
            if meta_use["acc"] == "AUTO_REFSEQGENE":
                cached_hits = sorted(o_dir.glob(f"{name}_NG_*.fna.gz"))
                if cached_hits:
                    fp = cached_hits[0]
                    meta_use["acc"] = fp.name.replace(f"{name}_", "").replace(".fna.gz", "")
                elif args.auto_download:
                    real_acc = _resolve_refseqgene_acc(name, expected_chr=meta_use["chr"])
                    print(f"  [NCBI RESOLVE] {name} ({meta_use['chr']}) -> {real_acc}")
                    meta_use["acc"] = real_acc
                    fp = _download_with_retry_v2(name, meta_use, o_dir)
                else:
                    print(f"  [skip onco {i:2d}/28] {name} not cached (use --auto-download)")
                    continue
            else:
                fp = o_dir / f"{name}_{meta_use['acc']}.fna.gz"
                if not fp.exists() or fp.stat().st_size <= 150:
                    if args.auto_download:
                        fp = _download_with_retry_v2(name, meta_use, o_dir)
                    else:
                        print(f"  [skip onco {i:2d}/28] {name} not cached (use --auto-download)")
                        continue
            nt, hdr = load_chromosome_nt(fp)
            hdr_ok = ((f"({name})" in hdr) or (f"({name}@)" in hdr)) and (meta_use["acc"] in hdr)
            assert hdr_ok, f"Accession/Header mismatch for {name} ({meta_use['acc']}): got {hdr!r}"
            nv = int(np.count_nonzero(nt < 4))
            if nv < 1000:
                continue
            B3 = count_triplet_pairs_lag(nt, d=3)
            B_human_fallback += B3
            kv = klein_v4_pr2_audit(B3)
            sc = two_sided_screen_full(B3)
            a11 = thm_D11_arrow_screen_and_cpg_audit(B3)
            mut_v3 = analyze_mutagenesis_signatures(nt, eps_eval=0.05)
            pal_sp = _spaced_palindrome_spectrum(nt)
            cnt4 = np.bincount(nt[nt < 4], minlength=4).astype(np.float64)
            gc = float((cnt4[1] + cnt4[2]) / max(float(cnt4.sum()), 1.0))
            onco_nt_cache[name] = nt
            rec = {
                "label": name, "acc": meta_use["acc"], "chr": meta_use["chr"], "role": meta_use["role"],
                "notes": meta["notes"], "length_bp": int(nt.size), "valid_bp": nv, "gc": gc,
                "fasta_sha256": sha256_file(fp),
                "klein_d3": strip_private(kv), "screen_d3": sc, "arrow_D11": a11,
                "palindrome_spectrum": pal_sp, "mutagenesis_V3": mut_v3,
                "_W_mp": kv["_W_mp"], "_B3": B3,
            }
            onco_results.append(rec)
            ap_ch = 100.0 * _fval((mut_v3.get("SBS2_APOBEC_plus") or {}).get("J_chiral_share_of_directed_unb"), 0.0)
            s1_ch = 100.0 * _fval((mut_v3.get("SBS1_CpG_sym") or {}).get("J_chiral_share_of_directed_unb"), 0.0)
            print(f"  [ONCO {i:2d}/28] {name:7s} ({meta['chr']:5s}, {nt.size/1e3:6.1f} kb) | "
                  f"PR2={fnum(100*_fval(kv.get('time_PR2_compliant_share_unb'),0),2)}% | "
                  f"NonEnv={fnum(100*_fval(a11.get('pure_nonenv_arrow_share'),0),1)}% | "
                  f"SBS1_ch={s1_ch:.2f}% vs APOBEC_ch={ap_ch:.1f}%")
        except Exception as e:
            print(f"  [ERROR onco] {name}: {e}")
            onco_failures.append({"label": name, "error": repr(e)})

    if human_ref is None and float(B_human_fallback.sum()) > 0.0:
        kv_fb = klein_v4_pr2_audit(B_human_fallback)
        P_m1_fb = exact_markov1_from_pair_table(B_human_fallback)
        kv_m1_fb = klein_v4_pr2_audit(P_m1_fb)
        P_m2_fb = exact_markov2_from_pair_table(B_human_fallback)
        kv_m2_fb = klein_v4_pr2_audit(P_m2_fb)
        human_ref = {"source": "pooled_28_RefSeqGene_loci", "n_chr": 0, "grand_W_mp": kv_fb["_W_mp"], "grand_W_mp_m1": kv_m1_fb["_W_mp"], "grand_W_mp_perp_m1": kv_fb["_W_mp"] - kv_m1_fb["_W_mp"], "grand_W_mp_m2": kv_m2_fb["_W_mp"], "grand_W_mp_perp_m2": kv_fb["_W_mp"] - kv_m2_fb["_W_mp"], "chr_W_mp": {}}

    for r in onco_results:
        wmp = r.pop("_W_mp")
        b3_o = r.pop("_B3")
        if human_ref is not None:
            r["cos_Wmp_vs_human_grand"] = _cos_mat_frob(wmp, human_ref["grand_W_mp"])
            r["cos_Wmp_unb_V18"] = _cos_wmp_unbiased_vs_ref(b3_o, human_ref["grand_W_mp"], z_gate=3.0)

    # 2. Audit 33 Real NCBI RefSeq Viruses
    viral_results, viral_nt_cache, viral_failures = [], {}, []
    for i, (name, meta) in enumerate(VIRAL_CATALOG_33.items(), 1):
        fp = v_dir / f"{name}_{meta['acc']}.fna.gz"
        try:
            if not fp.exists() or fp.stat().st_size <= 150:
                if args.auto_download:
                    fp = _download_with_retry_v2(name, meta, v_dir)
                else:
                    print(f"  [skip virus {i:2d}/33] {name} not cached (use --auto-download)")
                    continue
            nt, hdr = load_chromosome_nt(fp)
            assert meta["acc"] in hdr, f"Viral accession mismatch for {name}: expected {meta['acc']}, got {hdr!r}"
            nv = int(np.count_nonzero(nt < 4))
            if nv < 500:
                continue
            B3 = count_triplet_pairs_lag(nt, d=3)
            kv = klein_v4_pr2_audit(B3)
            sc = two_sided_screen_full(B3)
            a11 = thm_D11_arrow_screen_and_cpg_audit(B3)
            pal_sp = _spaced_palindrome_spectrum(nt)
            cnt4 = np.bincount(nt[nt < 4], minlength=4).astype(np.float64)
            gc = float((cnt4[1] + cnt4[2]) / max(float(cnt4.sum()), 1.0))

            cos_grand = _cos_mat_frob(kv["_W_mp"], human_ref["grand_W_mp"]) if human_ref else None
            cos_unb_v = _cos_wmp_unbiased_vs_ref(B3, human_ref["grand_W_mp"], z_gate=3.0) if human_ref else None
            D4_v = count_dinuc_4x4(nt)
            s_d4 = float(D4_v.sum())
            cpg_oe_v = None
            if s_d4 > 10.0:
                P2_v = D4_v / s_d4
                p1_v = P2_v.sum(axis=1)
                exp_cg_v = float(p1_v[1] * p1_v[2])
                if exp_cg_v > 1e-12:
                    cpg_oe_v = float(P2_v[1, 2]) / exp_cg_v
            deff_v = klein_design_effect_audit(nt, d=3)
            chr_cos = {}
            best_chr, best_cos = None, -2.0
            if human_ref and human_ref.get("chr_W_mp"):
                for clbl, cwmp in human_ref["chr_W_mp"].items():
                    cv = _cos_mat_frob(kv["_W_mp"], cwmp)
                    if cv is not None:
                        chr_cos[clbl] = cv
                        if cv > best_cos:
                            best_cos, best_chr = cv, clbl

            viral_nt_cache[name] = nt

            rec = {
                "label": name, "acc": meta["acc"], "vclass": meta["vclass"], "mode": meta["mode"],
                "notes": meta["notes"], "length_bp": int(nt.size), "valid_bp": nv, "gc": gc,
                "fasta_sha256": sha256_file(fp),
                "klein_d3": strip_private(kv), "screen_d3": sc, "arrow_D11": a11,
                "palindrome_spectrum": pal_sp,
                "cpg_oe": cpg_oe_v,
                "cos_Wmp_vs_human_grand": cos_grand,
                "cos_Wmp_unb_V18": cos_unb_v,
                "design_effect_V18": deff_v,
                "best_human_chr": best_chr, "best_human_chr_cos": best_cos if best_chr else None,
                "cos_Wmp_by_chr24": chr_cos,
            }
            viral_results.append(rec)
            pr2_v = 100.0 * _fval(kv.get("time_PR2_compliant_share_unb", kv.get("time_PR2_compliant_share")), 0.0)
            print(f"  [VIR  {i:2d}/33] {name:14s} ({nt.size:7,d} bp, GC={100*gc:4.1f}%) | "
                  f"PR2={pr2_v:6.2f}% | NonEnv={fnum(100*_fval(a11.get('pure_nonenv_arrow_share'),0),1)}% | "
                  f"cos_Hum={fnum(cos_grand,3,True)} | best={best_chr or 'n/a'}")
        except Exception as e:
            print(f"  [ERROR virus] {name}: {e}")
            viral_failures.append({"label": name, "error": repr(e)})

    # 3. Theorem V2: Deterministic Proviral Insertions into Classic Human Integration Hotspots
    ins_pairs = [
        ("HBV", "TERT"), ("HBV", "MYC"),
        ("HPV16", "MYC"), ("HPV16", "FHIT"), ("HPV18", "MYC"), ("HPV6b", "MYC"),
        ("HIV1_HXB2", "MYC"), ("HIV1_HXB2", "TP53"), ("HIV1_HXB2", "APOE"),
        ("HTLV1", "TP53"), ("MCPyV", "TP53"), ("AAV2", "APOE"),
    ]
    ins_rows = []
    for vn, hn in ins_pairs:
        if vn in viral_nt_cache and hn in onco_nt_cache:
            res_i = deterministic_proviral_insertion_audit(onco_nt_cache[hn], viral_nt_cache[vn])
            if res_i is not None:
                res_i["virus"] = vn
                res_i["host"] = hn
                ins_rows.append(res_i)
                print(f"  [THM-V2 INS] {vn:10s} -> {hn:6s} | +PR2={fnum(100*_fval(res_i['plus_PR2_ok_unb'],0),2)}% | "
                      f"W_-- SNR={fnum(res_i['plus_W_mm_snr_over_4N'],1)}x | HH cancel={fnum(res_i['chiral_cancellation_ratio_hh'],1)}x")

    # 4. Theorems V5, V6, V7: T2T Neighborhood Micro-Scan, 10-Lag Spectrum & 1-kb Vulnerability Tiles
    ref_npz_for_lags = getattr(args, "atlas_npz", "") or (str(local_ref_npz) if local_ref_npz.exists() else "")
    human_euch_lags = _load_human_euch_multilag_arrows(ref_npz_for_lags)
    for r in viral_results:
        vnt = viral_nt_cache.get(r["label"])
        if vnt is not None:
            r["multilag_V6"] = compute_multilag_profile_v6(vnt, human_euch_lags=human_euch_lags)
            r["mutagenesis_V3"] = analyze_mutagenesis_signatures(vnt, eps_eval=0.05)
            r["beyond_dinuc_V11"] = viral_beyond_dinuc_mimicry_v11(vnt, count_triplet_pairs_lag(vnt, d=3), human_ref)
    for r in onco_results:
        ont = onco_nt_cache.get(r["label"])
        if ont is not None:
            r["multilag_V6"] = compute_multilag_profile_v6(ont, human_euch_lags=human_euch_lags)
            r["tile_vulnerability_V7"] = scan_locus_vulnerability_tiles_v7(ont, tile_bp=1000)
    t2t_rows = audit_t2t_neighborhood_v5(c_dir, onco_nt_cache, viral_nt_cache, human_ref, auto_download=args.auto_download)
    syn_v8 = audit_synonymous_duplex_barrier_v8(onco_nt_cache, onco_results)
    grp_v9 = compute_virome_onco_group_tests_v9(viral_results, onco_results, t2t_rows=t2t_rows)
    vir_v10 = compute_viral_apobec_restriction_v10(viral_results, viral_nt_cache=viral_nt_cache)
    fd_v12 = {}
    for fname in ("TP53", "PIK3CA", "TERT", "HRAS", "MYC", "BRCA1", "FHIT"):
        if fname in onco_nt_cache:
            fd_v12[fname] = analyze_finite_dose_curve_v12(onco_nt_cache[fname])
    for vname in ("HPV16", "HIV1_HXB2", "MCPyV", "HBV"):
        if vname in viral_nt_cache:
            fd_v12[vname] = analyze_finite_dose_curve_v12(viral_nt_cache[vname])
    # Spearman GC/composition controls + THM-V13 exact Klein-orbit 6-mer attribution
    sp_v11 = {
        "gc_vs_cos_total": exact_spearman_rho([r["gc"] for r in viral_results], [r.get("cos_Wmp_vs_human_grand") for r in viral_results]),
        "gc_vs_cos_beyond_m2": exact_spearman_rho(
            [r["gc"] for r in viral_results],
            [(r.get("beyond_dinuc_V11") or {}).get("cos_markov2_resid_vs_human") for r in viral_results],
        ),
        "tcw_oe_vs_apobec_denovo": exact_spearman_rho(
            [x.get("tcw_oe_base") for x in vir_v10.get("rows", [])],
            [x.get("sbs2_chi_mm_denovo_x1e4") for x in vir_v10.get("rows", [])],
        ),
        "log10L_vs_cos_beyond_m2": exact_spearman_rho(
            [math.log10(max(float(r["length_bp"]), 1.0)) for r in viral_results],
            [(r.get("beyond_dinuc_V11") or {}).get("cos_markov2_resid_vs_human") for r in viral_results],
        ),
        "beyond_m2_int_vs_rna_adj_gc_L": residualized_group_test_v11(
            viral_results,
            {"retrovirus_int", "dna_onco_small", "hepatitis_hcc", "herpes_large", "vector_dna"},
            {"rna_acute_ctrl"},
            lambda r: (r.get("beyond_dinuc_V11") or {}).get("cos_markov2_resid_vs_human"),
        ),
        "beyond_m2_nucDNA_vs_cytoRNA_adj_gc_L": residualized_group_test_v11(
            viral_results,
            {"dna_onco_small", "dna_benign_ctrl", "herpes_large", "vector_dna"},
            {"rna_acute_ctrl"},
            lambda r: (r.get("beyond_dinuc_V11") or {}).get("cos_markov2_resid_vs_human"),
        ),
        "beyond_m2_gs_nucDNA_vs_cytoRNA_adj_gc_L": residualized_group_test_v11(
            viral_results,
            {"dna_onco_small", "dna_benign_ctrl", "herpes_large", "vector_dna"},
            {"rna_acute_ctrl"},
            lambda r: (r.get("beyond_dinuc_V11") or {}).get("cos_markov2_gs_orth_vs_human"),
        ),
        "beyond_m2_int_vs_rna_adj_gc_L_cpg": residualized_group_test_v11(
            viral_results,
            {"retrovirus_int", "dna_onco_small", "hepatitis_hcc", "herpes_large", "vector_dna"},
            {"rna_acute_ctrl"},
            lambda r: (r.get("beyond_dinuc_V11") or {}).get("cos_markov2_resid_vs_human"),
            include_cpg_oe=True,
        ),
        "beyond_m2_nucDNA_vs_cytoRNA_adj_gc_L_cpg": residualized_group_test_v11(
            viral_results,
            {"dna_onco_small", "dna_benign_ctrl", "herpes_large", "vector_dna"},
            {"rna_acute_ctrl"},
            lambda r: (r.get("beyond_dinuc_V11") or {}).get("cos_markov2_resid_vs_human"),
            include_cpg_oe=True,
        ),
        "beyond_m2_gs_nucDNA_vs_cytoRNA_adj_gc_L_cpg": residualized_group_test_v11(
            viral_results,
            {"dna_onco_small", "dna_benign_ctrl", "herpes_large", "vector_dna"},
            {"rna_acute_ctrl"},
            lambda r: (r.get("beyond_dinuc_V11") or {}).get("cos_markov2_gs_orth_vs_human"),
            include_cpg_oe=True,
        ),
    }
    # Global BH-FDR ledger across all V9 + V10 + V11 tests (explicit primary vs exploratory split)
    fdr_refs = []
    for t_r in grp_v9:
        if _finite(t_r.get("p_twosided")):
            fdr_refs.append(("V9", f"{t_r.get('contrast')} | {t_r.get('metric')}", "exploratory", t_r))
    for k10 in ("mw_tcw_int_vs_rna", "mw_tcw_oe_base_int_vs_rna", "mw_tcw_oe_dinuc_int_vs_rna",
                "mw_tcw_oe_base_retro_pap_vs_rna", "mw_sbs2_snr_int_vs_rna"):
        t10 = vir_v10.get(k10)
        if isinstance(t10, dict) and _finite(t10.get("p_twosided")):
            fdr_refs.append(("V10", k10, "exploratory", t10))
    primary_v11 = {"tcw_oe_vs_apobec_denovo", "beyond_m2_nucDNA_vs_cytoRNA_adj_gc_L"}
    for k11, t11 in sp_v11.items():
        if isinstance(t11, dict) and _finite(t11.get("p_twosided")):
            fdr_refs.append(("V11", k11, "primary" if k11 in primary_v11 else "exploratory", t11))
    q_all = benjamini_hochberg_fdr([ref[3]["p_twosided"] for ref in fdr_refs])
    fdr_rows = []
    for (sec_t, name_t, tier_t, obj_t), qv in zip(fdr_refs, q_all):
        obj_t["q_bh_all"] = float(qv)
        obj_t["tier"] = tier_t
        fdr_rows.append({
            "section": sec_t, "test": name_t, "tier": tier_t,
            "p_twosided": float(obj_t["p_twosided"]),
            "p_min_discrete": obj_t.get("p_min_discrete"),
            "q_bh_all": float(qv),
        })
    fdr_rows.sort(key=lambda d: d["p_twosided"])
    fdr_ledger = {
        "m_total_tests": len(fdr_rows),
        "n_primary": sum(1 for d in fdr_rows if d["tier"] == "primary"),
        "n_exploratory": sum(1 for d in fdr_rows if d["tier"] == "exploratory"),
        "n_bh_sig_q05": sum(1 for d in fdr_rows if d["q_bh_all"] <= 0.05),
        "rows": fdr_rows,
    }
    t2t_rank1 = compute_t2t_rank1_diagnostics_v18(viral_results, t2t_rows)
    att_v13 = {}
    if human_ref is not None:
        for vn in ("KSHV_HHV8", "EBV_HHV4"):
            if vn in viral_nt_cache:
                vnt = viral_nt_cache[vn]
                Wv = klein_v4_pr2_audit(count_triplet_pairs_lag(vnt, d=3))["_W_mp"]
                Wv_m2 = klein_v4_pr2_audit(exact_markov2_from_trinuc_counts(count_triplets_1d(vnt)))["_W_mp"]
                att_v13[f"{vn}_vs_human_total"] = klein_orbit_attribution_v13(Wv, human_ref["grand_W_mp"], top_k=10)
                if human_ref.get("grand_W_mp_perp_m2") is not None:
                    att_v13[f"{vn}_vs_human_beyond_M2"] = klein_orbit_attribution_v13(Wv - Wv_m2, human_ref["grand_W_mp_perp_m2"], top_k=10)
    for g_att in ("TP53", "PIK3CA"):
        if g_att in onco_nt_cache:
            B0_a, J_a, _ = exact_signature_jacobian(onco_nt_cache[g_att], COSMIC_SIG_TABLES["SBS2_APOBEC_plus"]["rates"])
            Na = max(float(B0_a.sum()), 1.0)
            dMa = J_a / Na
            Jmm_a = 0.25 * (dMa - dMa.T - dMa.T[np.ix_(MIRROR, MIRROR)] + dMa[np.ix_(MIRROR, MIRROR)])
            att_v13[f"{g_att}_APOBEC_SBS2_Jmm"] = klein_orbit_attribution_v13(Jmm_a, None, top_k=10)
    print(f"  [THM-V11..V13] Markov-1/2 mimicry (33 vir), all-orders dose curves ({len(fd_v12)} targets) & Klein-orbit attribution ({len(att_v13)} maps) computed.")

    payload = {
        "version": VERSION,
        "provenance": collect_provenance(),
        "self_test": st,
        "human_reference_source": human_ref["source"] if human_ref else None,
        "human_reference_baselines": {
            "grand_PR2_ok_unb": human_ref.get("grand_PR2_ok_unb") if human_ref else None,
            "grand_nonenv_3600_unb": human_ref.get("grand_nonenv_3600_unb") if human_ref else None,
            "grand_nonenv_3600_raw": human_ref.get("grand_nonenv_3600_raw") if human_ref else None,
        },
        "n_viruses_ok": len(viral_results), "n_viruses_failed": len(viral_failures), "viral_failures": viral_failures,
        "n_onco_ok": len(onco_results), "n_onco_failed": len(onco_failures), "onco_failures": onco_failures,
        "viruses": viral_results,
        "onco_loci": onco_results,
        "proviral_insertions_V2": ins_rows,
        "t2t_neighborhoods_V5": t2t_rows,
        "t2t_rank1_V18": t2t_rank1,
        "synonymous_duplex_barrier_V8": syn_v8,
        "group_discrimination_V9": grp_v9,
        "viral_apobec_restriction_V10": vir_v10,
        "finite_dose_curves_V12": fd_v12,
        "spearman_controls_V11": sp_v11,
        "fdr_ledger_V18": fdr_ledger,
        "orbit_attribution_V13": att_v13,
        "runtime_seconds": time.time() - t0,
    }
    jp = out / "virome_onco_audit_v2.json"
    mp = out / "virome_onco_audit_v2.md"
    atomic_write(jp, json.dumps(_clean(payload), indent=2))
    write_virome_onco_markdown(payload, mp)
    figs = [] if args.no_plots else make_virome_onco_figures(payload, out, dpi=args.dpi)
    figs += append_v21_markdown_and_figures(payload, out, mp, dpi=args.dpi)
    figs += append_v22_markdown_and_figures(payload, out, mp, dpi=args.dpi)
    figs += append_v23_markdown_and_figures(payload, out, mp, dpi=args.dpi)
    print("\n" + "=" * 96)
    print(f"VIROME & ONCO AUDIT DONE (v{VERSION}): {len(viral_results)}/33 viruses, {len(onco_results)}/28 loci, "
          f"{len(ins_rows)} insertions in {time.time()-t0:.1f} s")
    print(f"[JSON] {jp}\n[MD]   {mp}\n[FIGS] {len(figs)} figure(s)")
    print("=" * 96)
    return 0


def main(argv=None):
    import gc
    p = argparse.ArgumentParser(description=f"Human T2T-CHM13v2.0 Complete 24-Chr Atlas v{VERSION} (Zero-RNG)")
    p.add_argument("--chr-dir", default=str(SCRIPT_DIR / "chr_cache"))
    p.add_argument("--out-dir", default="")
    p.add_argument("--self-test", action="store_true")
    p.add_argument("--all-24", action="store_true")
    p.add_argument("--chromosomes", default="", help="Subset, e.g. chr4,chr19,chr21,chrY")
    p.add_argument("--auto-download", action="store_true")
    p.add_argument("--resume", action="store_true")
    p.add_argument("--sector-bp", type=int, default=100000)
    p.add_argument("--euch-h-thresh", type=float, default=0.0, help="<=0 uses parameter-free THM-D7 Otsu-BIC (default); >0 uses legacy fixed H_3mer threshold")
    p.add_argument("--micro-lags", default=",".join(map(str, MICRO_LAGS_DEFAULT)))
    p.add_argument("--replot-only", type=str, default="")
    p.add_argument("--euch-multilag-check", type=str, default="", help="Optional per-chr pure-euchromatin 10-lag check, e.g. chr2,chr19")
    p.add_argument("--skip-diskmaps", action="store_true", help="Skip regenerating the 24 individual diskmaps during --replot-only")
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument("--no-plots", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--virome-onco-audit", action="store_true", help="Run Paper 5 Zero-RNG 33-Virus & 28-Cancer-Locus Audit")
    p.add_argument("--viral-dir", default=str(SCRIPT_DIR / "viral_cache"))
    p.add_argument("--onco-dir", default=str(SCRIPT_DIR / "onco_cache"))
    p.add_argument("--atlas-npz", default="", help="Optional explicit path to a precomputed .npz reference")
    p.add_argument("--fresh-cache", action="store_true", help="Wipe and re-download ALL local caches (viral, onco, 24 T2T chr & ref npz) inside SCRIPT_DIR")
    p.add_argument("--fresh-virome-cache", action="store_true", help="Wipe and re-download only viral_cache & onco_cache inside SCRIPT_DIR")
    args = p.parse_args(argv)
    _explicit_out_dir = bool(args.out_dir)
    if not args.out_dir:
        args.out_dir = str(SCRIPT_DIR / ("results_virome_onco_v2_8" if args.virome_onco_audit else "results_human_t2t_atlas_v2_8"))

    if args.self_test:
        print("=" * 96 + f"\nSELF-TEST — HUMAN T2T-CHM13v2.0 24-CHR ATLAS v{VERSION} (Theorems D1-D19 + V1-V18, Zero-RNG)\n" + "=" * 96)
        self_test(verbose=True); print("ALL SELF-TESTS PASSED."); return 0
    if args.replot_only:
        t_rp = time.time()
        jp = Path(args.replot_only); payload = json.loads(jp.read_text(encoding="utf-8"))
        out_p = Path(args.out_dir).resolve() if _explicit_out_dir else jp.parent
        if "provenance_run" not in payload and "provenance" in payload:
            payload["provenance_run"] = payload["provenance"]
        payload["version"] = VERSION
        payload["provenance_replot"] = collect_provenance()
        payload["self_test"] = self_test(verbose=False)
        apply_v15_post_audit(payload)
        if args.euch_multilag_check:
            ec_res = run_euch_2d_multilag_control(payload, args.chr_dir, args.euch_multilag_check)
            if ec_res is not None:
                payload["euchromatin_multilag_control_D12"] = ec_res
        atomic_write(out_p / "human_t2t_atlas_v1.json", json.dumps(_clean(payload), indent=2))
        write_human24_markdown(payload, out_p / "human_t2t_atlas_v1.md")
        figs = [] if args.no_plots else make_human24_figures(payload, out_p, dpi=args.dpi, skip_diskmaps=args.skip_diskmaps)
        print(f"[v1.5-REPLOT] Self-test + post-audit + {len(figs)} figure(s) + Markdown updated in {time.time()-t_rp:.1f} s -> {out_p}")
        return 0

    if args.virome_onco_audit:
        return run_virome_onco_audit(args)

    st = self_test(verbose=False)
    prov = collect_provenance()
    out = Path(args.out_dir).resolve(); out.mkdir(parents=True, exist_ok=True); chr_dir = Path(args.chr_dir)
    if args.chromosomes:
        sel_keys = [raw if raw in HUMAN_T2T_24 else f"Human_T2T_{raw}"
                    for raw in (s.strip() for s in args.chromosomes.split(",") if s.strip())
                    if raw in HUMAN_T2T_24 or f"Human_T2T_{raw}" in HUMAN_T2T_24]
    else:
        sel_keys = list(HUMAN_T2T_24.keys())

    targets = []
    for k in sel_keys:
        meta = HUMAN_T2T_24[k]
        fp = chr_dir / f"{k}_{meta['acc']}.fna.gz"
        if not fp.exists() and chr_dir.exists():
            hits = list(chr_dir.glob(f"*{meta['acc']}*.fna*"))
            if hits: fp = hits[0]
        if args.dry_run:
            print(f"  [{'CACHED' if fp.exists() else 'NEED_DOWNLOAD':13s}] {k:16s} | {meta['acc']:12s} | {meta['arm_class']:14s} | {meta['notes']}")
            continue
        targets.append((k, meta, fp))
    if args.dry_run:
        return 0

    print("=" * 96 + f"\nHUMAN T2T-CHM13v2.0 24-CHR ATLAS v{VERSION} | Target={len(targets)} chr | sector={args.sector_bp//1000} kb | Zero-RNG\n" + "=" * 96)
    partial_path, b3_npz_path = out / "human_t2t_atlas_v1.partial.json", out / "atlas_b3_cache.npz"
    results, failures, done_keys, b3_store = [], [], set(), {}
    if args.resume and partial_path.exists():
        try:
            pjs = json.loads(partial_path.read_text(encoding="utf-8"))
            if pjs.get("version") != VERSION:
                print(f"[RESUME] Ignoring checkpoint from older version {pjs.get('version')} != {VERSION}")
            else:
                results = pjs.get("results", [])
                done_keys = {r["key"] for r in results}
                if b3_npz_path.exists():
                    npz = np.load(b3_npz_path)
                    for fn in npz.files: b3_store[fn] = np.array(npz[fn])
                    npz.close()
                print(f"[RESUME] Loaded {len(results)}/24 chromosome(s) from {partial_path.name}")
        except Exception as e:
            print(f"[RESUME] Could not load checkpoint: {e}")

    t_all = time.time()
    order_idx = {ck: ci for ci, ck in enumerate(HUMAN_T2T_24.keys())}
    for i, (k, meta, fp) in enumerate(targets, 1):
        if k in done_keys:
            print(f"--- [{i:2d}/{len(targets)}] {k} (resumed, skipping) ---")
            continue
        print(f"\n--- [{i:2d}/{len(targets)}] {k} ({meta['acc']}) ---")
        try:
            if not fp.exists():
                if args.auto_download: fp = download_chromosome_fasta(k, meta, chr_dir)
                else: print(f"[skip] {k} not in {chr_dir} (use --auto-download)"); continue
            r = analyze_chromosome(k, meta, fp, args)
            if r is not None:
                for suf, fk in (("all", "_B3_all"), ("euch", "_B3_euch"), ("sat", "_B3_sat"), ("strata", "_B_strata")):
                    mat = r.pop(fk, None)
                    if mat is not None: b3_store[f"{k}__{suf}"] = mat
                for lag_d, bmat in (r.pop("_B_lags", None) or {}).items():
                    if bmat is not None: b3_store[f"{k}__lag_{int(lag_d)}"] = bmat
                for lag_d, bmat in (r.pop("_B_euch_lags", None) or {}).items():
                    if bmat is not None: b3_store[f"{k}__euch_lag_{int(lag_d)}"] = bmat
                r.pop("_P3_global", None)
                tmp_npz = out / "atlas_b3_cache.tmp.npz"
                np.savez_compressed(tmp_npz, **b3_store); os.replace(tmp_npz, b3_npz_path)
                results.append(r); results.sort(key=lambda item: order_idx.get(item["key"], 999))
                atomic_write(partial_path, json.dumps(_clean({"version": VERSION, "results": results}), indent=2))
            gc.collect()
        except Exception as e:
            import traceback; print(f"[ERROR] {k}: {e}"); traceback.print_exc()
            failures.append({"key": k, "error": repr(e)})

    grand = grand_genome_audit(results, b3_store)
    xmat = cross_chromosome_matrix_human24(results, b3_store)
    payload = {
        "version": VERSION, "assembly": "Homo sapiens T2T-CHM13v2.0 (GCF_009914755.1)",
        "provenance": prov,
        "self_test": st, "config": vars(args), "n_chromosomes": len(results),
        "n_failed": len(failures), "failures": failures, "runtime_seconds": time.time() - t_all,
        "grand_genome": grand, "cross_chromosome": xmat, "results": results,
    }
    apply_v15_post_audit(payload)
    if args.euch_multilag_check:
        ec_res = run_euch_2d_multilag_control(payload, chr_dir, args.euch_multilag_check)
        if ec_res is not None:
            payload["euchromatin_multilag_control_D12"] = ec_res
    json_path, md_path = out / "human_t2t_atlas_v1.json", out / "human_t2t_atlas_v1.md"
    atomic_write(json_path, json.dumps(_clean(payload), indent=2))
    write_human24_markdown(payload, md_path)
    figs = [] if args.no_plots else make_human24_figures(payload, out, dpi=args.dpi, skip_diskmaps=args.skip_diskmaps)

    print("\n" + "=" * 96 + f"\nHUMAN T2T ATLAS DONE (v{VERSION}): {len(results)} OK, {len(failures)} failed in {time.time()-t_all:.1f} s")
    if grand.get("whole_genome"):
        eg = grand.get("euchromatin_grand") or {}
        ge = eg.get("klein_d3") or {}
        ga = eg.get("arrow_D11") or {}
        gs = (grand.get("satellite_grand") or {}).get("klein_d3") or {}
        pr2_eu = ge.get("time_PR2_compliant_share_unb")
        if pr2_eu is None: pr2_eu = ge.get("time_PR2_compliant_share")
        pr2_sa = gs.get("time_PR2_compliant_share_unb")
        if pr2_sa is None: pr2_sa = gs.get("time_PR2_compliant_share")
        n_epd = sum(1 for row in (grand.get("entropy_phase_diagram_D10") or []) if not row.get("skipped"))
        n_ml = len(grand.get("grand_multilag_spectrum_D12") or [])
        n_eml = len(grand.get("grand_euch_multilag_spectrum_D12") or [])
        print(f"[GRAND] Total={grand['total_gb']:.3f} Gb | Euch PR2-ok(unb)={fnum(100*(pr2_eu or 0),2)}% | "
              f"Sat PR2-ok(unb)={fnum(100*(pr2_sa or 0),2)}% | Euch NonEnv3600={fnum(100*(ga.get('pure_nonenv_arrow_share') or 0),1)}% | "
              f"CpG-free={fnum(100*(ga.get('cpg_free_arrow_share') or 0),1)}% | D10={n_epd}/8 | D12={n_ml}w/{n_eml}e lags")
    print(f"[JSON] {json_path}\n[MD]   {md_path}\n[FIGS] {len(figs)} figure(s) in {out / 'figures'}\n" + "=" * 96)
    return 0


if __name__ == "__main__":
    sys.exit(main())
