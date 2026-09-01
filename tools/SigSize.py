"""
SHUTTLE-NIST signature-size analysis: Shannon entropy lower bound and
rANS reservation sizing.

Signature structure (single merged rANS stream)
-----------------------------------------------
    sigma = seedC || rANS(HighBits(z^(0..lenS)) || hint)
                  || polyz0_lo_pack || L * polyz1_lo_pack

HighBits(z) and hint share a single rANS context. Their scales are
comparable, but their PMFs are not identical: HighBits(z) is a rounded
discrete Gaussian and hint follows a bucket-crossing distribution. The
single table is generated from their coefficient-count-weighted mixture.

Per-component Shannon entropy (bit/coef, continuous-Gaussian
approximation):

    H(z^(0))   = (1/2) * log_2(2*pi*e * (r/alpha_1)^2)              # post-CompressY
    H(z^(i))   = (1/2) * log_2(2*pi*e * r^2)         (i = 1..lenS)  # main slots
    H(hint)    = (1/2) * log_2(2*pi*e * (2*r/alpha_h)^2)
               = (1/2) * log_2(8*pi*e * (r/alpha_h)^2)

Physical constraint: a discrete random variable has non-negative
entropy. When sigma < 1/sqrt(2*pi*e) ~= 0.242 the continuous formula
goes negative -- an artifact of the formula breaking down, not a real
quantity. The exact discrete entropy stays >= 0; we take max(H, 0).

Total entropy (bytes, including seedC):

    B_sig_lower = (n*H(z^(0)) + n*lenS*H(z^(i)) + n*lenE*max(H(hint), 0)) / 8
                + seedC_B

Actual rANS-encoded bytes exceed the Shannon lower bound slightly due to:
  - probability quantization (t = 10 bit frequency rounding);
  - 4 B flush + 2 B length prefix per stream;
  - tail margin Delta (sized by Gauss CLT at p_rans* = 2^-30).
"""

import math
from typing import Dict, List, Tuple


# -------------------------------------------------------------------------
# SHUTTLE-NIST-{120, 180, 255} parameter sets.
# -------------------------------------------------------------------------
SHUTTLE_NIST_PARAMS: Dict[str, dict] = {
    "SHUTTLE-NIST-120": {
        "n":       256,
        "q":       32257,
        "lenS":    3,
        "lenE":    2,
        "eta":     1,
        "tau":     58,
        "r":       70,
        "alpha_1": 8,
        "alpha_h": 128,
        "T_sigma": 11,
        "seedC_B": 32,
    },
    "SHUTTLE-NIST-180": {
        "n":       256,
        "q":       15361,
        "lenS":    4,
        "lenE":    3,
        "eta":     1,
        "tau":     80,
        "r":       82,
        "alpha_1": 8,
        "alpha_h": 128,
        "T_sigma": 11,
        "seedC_B": 48,
    },
    "SHUTTLE-NIST-255": {
        "n":       256,
        "q":       7681,
        "lenS":    5,
        "lenE":    4,
        "eta":     1,
        "tau":     128,
        "r":       99,
        "alpha_1": 16,
        "alpha_h": 128,
        "T_sigma": 11,
        "seedC_B": 64,
    },
}


def gauss_entropy_continuous(sigma: float) -> float:
    """Continuous-Gaussian entropy (bit):

        H(N(0, sigma^2)) = (1/2) * log_2(2*pi*e * sigma^2)

    Note: when sigma < 1/sqrt(2*pi*e) ~= 0.242 the formula gives a
    negative value (an artifact of the formula breaking down for narrow
    sigmas). The true discrete entropy stays >= 0; see
    gauss_entropy_discrete.
    """
    if sigma <= 0:
        return float("-inf")
    return 0.5 * math.log2(2 * math.pi * math.e * sigma * sigma)


def gauss_entropy_discrete(sigma: float, max_k: int = 50) -> float:
    """Exact entropy of D_{Z, sigma} at center 0, computed by per-bin
    PMF summation:

        p(k) propto exp(-k^2 / (2 * sigma^2)),   k in [-max_k, max_k]

    With max_k large enough the tail contribution is negligible
    (max_k = 50 covers more than 5*sigma for sigma <= 30).

    Comparison with the continuous formula:
      - sigma >> 1: H_discrete ~ H_continuous, deviation O(e^{-pi*sigma^2}).
      - sigma << 1: H_discrete -> 0+ (mass concentrates at 0),
        whereas H_continuous -> -inf.
    """
    if sigma <= 0:
        return 0.0
    # unnormalized weights
    weights = [math.exp(-k * k / (2 * sigma * sigma)) for k in range(-max_k, max_k + 1)]
    Z = sum(weights)
    H = 0.0
    for w in weights:
        p = w / Z
        if p > 0:
            H -= p * math.log2(p)
    return H


def compute_shannon_lower_bound(p: dict) -> dict:
    """Strict Shannon entropy lower bound on sig size (bytes), with the
    per-component contributions.

    The returned dict contains three estimates of the hint contribution:
      - hint_B_continuous: continuous formula evaluated literally
                           (may be negative).
      - hint_B_physical:   continuous formula clipped at 0
                           (physical Shannon lower bound).
      - hint_B_discrete:   exact discrete entropy (closest to the real
                           discrete distribution).

    Primary result: total_B_physical.
    """
    n       = p["n"]
    lenS    = p["lenS"]
    lenE    = p["lenE"]
    r       = p["r"]
    alpha_1 = p["alpha_1"]
    alpha_h = p["alpha_h"]
    seedC_B = p["seedC_B"]

    # --- 1) Slot 0 z^(0) (post-CompressY, scale r/alpha_1) ---
    sigma_z0 = r / alpha_1
    H_z0     = gauss_entropy_continuous(sigma_z0)
    z0_B     = n * H_z0 / 8.0

    # --- 2) Main slots z^(i), i = 1..lenS (scale r) ---
    H_zi     = gauss_entropy_continuous(r)
    z_main_B = n * lenS * H_zi / 8.0

    # --- 3) Hint h (scale sigma_h = 2*r/alpha_h; the factor 2 comes from
    #               the 2*z_2 term inside MakeHint) ---
    sigma_h         = 2.0 * r / alpha_h
    H_h_continuous  = gauss_entropy_continuous(sigma_h)
    H_h_physical    = max(H_h_continuous, 0.0)
    H_h_discrete    = gauss_entropy_discrete(sigma_h)
    hint_B_continuous = n * lenE * H_h_continuous / 8.0
    hint_B_physical   = n * lenE * H_h_physical   / 8.0
    hint_B_discrete   = n * lenE * H_h_discrete   / 8.0

    return {
        "sigma_z0":          sigma_z0,
        "sigma_h":           sigma_h,
        "H_z0":              H_z0,
        "H_zi":              H_zi,
        "H_h_continuous":    H_h_continuous,
        "H_h_physical":      H_h_physical,
        "H_h_discrete":      H_h_discrete,
        "z0_B":              z0_B,
        "z_main_B":          z_main_B,
        "hint_B_continuous": hint_B_continuous,
        "hint_B_physical":   hint_B_physical,
        "hint_B_discrete":   hint_B_discrete,
        "seedC_B":           seedC_B,
        "total_B_continuous": z0_B + z_main_B + hint_B_continuous + seedC_B,
        "total_B_physical":   z0_B + z_main_B + hint_B_physical   + seedC_B,
        "total_B_discrete":   z0_B + z_main_B + hint_B_discrete   + seedC_B,
    }


# ============================================================================
# rANS stream reservation analysis.
# ============================================================================
#
# One rANS stream encodes n_b i.i.d. symbols. The output byte length
# L_bytes is a random variable; we want a reservation R such that
# P[L_bytes > R] <= p_ovf_star.
#
# Two models:
#   1. Gauss CLT (works for lambda >= 30 and a PMF that is not too skewed):
#        L_bytes ~ N(mu_b, sigma_b^2)
#        R = mu_b + Phi^{-1}(1 - p_star) * sigma_b
#   2. Poisson (for lambda < 30 -- in particular for very skewed PMFs
#      with lambda < 1):
#        N = number of nonzero symbols ~ Poisson(lambda), lambda = n_b * q,
#        where q = P[nonzero]
#        L_bits ~ const + ell(+/-1) * N
#        R = ceil((const + ell(+/-1) * N_max(p_star)) / 8)
#
# Switching threshold: use Poisson when lambda < 30.
# ============================================================================


def phi_inv(p: float) -> float:
    """Standard-normal inverse CDF (quantile function), Phi^{-1}(p).

    Used for tail-quantile queries: for a target tail probability q,
    take Phi^{-1}(1 - q).

    Implementation: a fast Acklam approximation gives the initial point,
    then one Newton step on Phi(x) = p polishes it.
    """
    if p <= 0.0:
        return float("-inf")
    if p >= 1.0:
        return float("inf")
    # Acklam approximation
    a = [-3.969683028665376e+01,  2.209460984245205e+02,
         -2.759285104469687e+02,  1.383577518672690e+02,
         -3.066479806614716e+01,  2.506628277459239e+00]
    b = [-5.447609879822406e+01,  1.615858368580409e+02,
         -1.556989798598866e+02,  6.680131188771972e+01,
         -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01,
         -2.400758277161838e+00, -2.549732539343734e+00,
          4.374664141464968e+00,  2.938163982698783e+00]
    d = [ 7.784695709041462e-03,  3.224671290700398e-01,
          2.445134137142996e+00,  3.754408661907416e+00]
    p_low, p_high = 0.02425, 1 - 0.02425
    if p < p_low:
        q = math.sqrt(-2 * math.log(p))
        x = (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) / \
            ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1)
    elif p <= p_high:
        q = p - 0.5
        r = q * q
        x = (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5]) * q / \
            (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1)
    else:
        q = math.sqrt(-2 * math.log(1 - p))
        x = -(((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) / \
             ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1)
    # one Newton refinement
    e = 0.5 * (1 + math.erf(x / math.sqrt(2))) - p
    pdf = math.exp(-x * x / 2) / math.sqrt(2 * math.pi)
    if pdf > 0:
        x -= e / pdf
    return x


def quantized_pmf(sigma: float, t: int, M_voc: int) -> Tuple[Dict[int, float], Dict[int, int]]:
    """Quantize the discrete Gaussian D_{Z, sigma} on the alphabet
    [-M_voc, M_voc] with prob_bits = t:
      - p[s] = true normalized probability (float), over s in [-M_voc, M_voc]
      - g[s] = quantized frequency (int, sum = 2^t, each >= 1).

    M_voc is the rANS vocabulary radius; it determines which symbols
    receive freq >= 1. Symbols outside [-M_voc, M_voc] are not in the
    vocabulary and do not enter the normalization.

    Suitable for the z-hi context (its PMF is genuinely discrete
    Gaussian). In the hint context with sigma_h < 1 this PMF severely
    underestimates the nonzero rate; use quantized_pmf_hint instead.
    """
    raw = {k: math.exp(-k * k / (2 * sigma * sigma)) for k in range(-M_voc, M_voc + 1)}
    Z = sum(raw.values())
    p = {k: w / Z for k, w in raw.items()}

    # quantize: g[s] = max(1, round(p[s] * 2^t)); absorb rounding slack
    # in the largest bucket.
    quanta = 2 ** t
    g = {k: max(1, round(prob * quanta)) for k, prob in p.items()}
    diff = quanta - sum(g.values())
    if diff != 0:
        kmax = max(g.keys(), key=lambda k: g[k])
        g[kmax] += diff
    return p, g


def quantize_pmf(p: Dict[int, float], t: int) -> Dict[int, int]:
    """Quantize an arbitrary PMF to positive rANS frequencies."""
    quanta = 2 ** t
    g = {k: max(1, round(prob * quanta)) for k, prob in p.items()}
    diff = quanta - sum(g.values())
    if diff != 0:
        kmax = max(g.keys(), key=lambda k: g[k])
        g[kmax] += diff
    return g


def round_half_up_ratio(x: int, step: int) -> int:
    """Integer round-half-up of x/step, matching the C implementation."""
    return math.floor((x + step // 2) / step)


def rounded_response_pmf(r: float, step: int, M_voc: int,
                         inner_step: int = 0,
                         T_sigma: float = 11.0) -> Dict[int, float]:
    """PMF of one response symbol after the implementation's rounding."""
    bound = int(T_sigma * r)
    weights = {x: math.exp(-x * x / (2 * r * r))
               for x in range(-bound, bound + 1)}
    Z = sum(weights.values())
    pmf = {s: 0.0 for s in range(-M_voc, M_voc + 1)}
    for x, weight in weights.items():
        y = x if inner_step == 0 else round_half_up_ratio(x, inner_step)
        symbol = round_half_up_ratio(y, step)
        if symbol in pmf:
            pmf[symbol] += weight / Z
    total = sum(pmf.values())
    return {s: probability / total for s, probability in pmf.items()}


def hint_pmf_correct(r: float, alpha_h: int, M_voc: int,
                     T_sigma: float = 11.0) -> Dict[int, float]:
    """True MakeHint output PMF (bucket-crossing model).

    Key formula:

       hint_k = round(u) - round(u - delta),  u ~ Uniform[0,1),
                                              delta = 2z/alpha_h

    Average over z ~ D_{Z, r} (truncated at 11*sigma). This PMF is the
    *correct* hint distribution and replaces the legacy discrete-Gaussian
    D_{Z, 2r/alpha_h} approximation (which is only good when
    sigma_h >> 1).
    """
    M_z = int(T_sigma * r)
    z_weights = {z: math.exp(-z * z / (2 * r * r))
                 for z in range(-M_z, M_z + 1)}
    Z = sum(z_weights.values())
    z_pmf = {z: w / Z for z, w in z_weights.items()}

    pmf: Dict[int, float] = {h: 0.0 for h in range(-M_voc, M_voc + 1)}
    for z, pz in z_pmf.items():
        delta = 2.0 * z / alpha_h
        k = int(math.floor(delta + 0.5))    # round half up
        f = delta - k
        if -M_voc <= k <= M_voc:
            pmf[k] += pz * (1.0 - abs(f))
        if abs(f) > 1e-15:
            sign = 1 if f > 0 else -1
            kp = k + sign
            if -M_voc <= kp <= M_voc:
                pmf[kp] += pz * abs(f)

    total = sum(pmf.values())
    return {h: p / total for h, p in pmf.items()}


GAUSS_REGIME_THRESH = 1.0   # regime split (sigma_h >= 1: Gaussian; < 1: Bernoulli)


def quantized_pmf_hint(r: float, alpha_h: int, t: int, M_voc: int
                       ) -> Tuple[Dict[int, float], Dict[int, int]]:
    """Like quantized_pmf, but for the hint context: regime-aware.

    sigma_h = 2r/alpha_h >= 1.0 -> discrete Gaussian D_{Z, sigma_h}.
    sigma_h < 1.0               -> bucket-crossing (Bernoulli-like) PMF.
    """
    sigma_h = 2.0 * r / alpha_h
    if sigma_h >= GAUSS_REGIME_THRESH:
        return quantized_pmf(sigma_h, t, M_voc)
    p = hint_pmf_correct(r, alpha_h, M_voc)
    quanta = 2 ** t
    g = {k: max(1, round(prob * quanta)) for k, prob in p.items()}
    diff = quanta - sum(g.values())
    if diff != 0:
        kmax = max(g.keys(), key=lambda k: g[k])
        g[kmax] += diff
    return p, g


def rans_rate_and_var(p: Dict[int, float], g: Dict[int, int], t: int) -> Tuple[float, float]:
    """Given a theoretical PMF p and quantized frequencies g, compute
    under rANS encoding:
      - E[ell] = sum_s p[s] * (-log2(g[s] / 2^t))
                 expected bits per symbol (includes the quantization loss)
      - Var[ell] = sum_s p[s] * (-log2(g[s] / 2^t))^2 - E[ell]^2

    These give mu_b and sigma_b in the Gauss CLT model.
    """
    quanta = 2 ** t
    E_ell, E_ell2 = 0.0, 0.0
    for s, prob in p.items():
        if prob > 0 and g[s] > 0:
            ell = -math.log2(g[s] / quanta)
            E_ell  += prob * ell
            E_ell2 += prob * ell * ell
    var_ell = E_ell2 - E_ell * E_ell
    return E_ell, max(var_ell, 0.0)  # numerical guard: var must not go negative


def gauss_reservation(n_b: int, E_ell: float, var_ell: float, p_star: float) -> Tuple[float, float, float]:
    """Gauss CLT model: R = mu_b - Phi^{-1}(p_star) * sigma_b.

    Returns (mu_b, sigma_b, R) in bytes:
      mu_b = n_b * E_ell / 8 + 4   (4 = flush bytes)
      sigma_b = sqrt(n_b * var_ell) / 8.
    """
    mu_b    = n_b * E_ell / 8.0 + 4.0
    sigma_b = math.sqrt(n_b * var_ell) / 8.0
    # Using the lower-tail quantile directly avoids losing p_star in the
    # subtraction 1-p_star for targets below machine epsilon.
    R       = mu_b - phi_inv(p_star) * sigma_b
    return mu_b, sigma_b, R


def poisson_tail_kmax(lam: float, p_star: float, k_search_max: int = 200) -> int:
    """Smallest k satisfying P[Poisson(lam) >= k] <= p_star.

    Returns N_max = k - 1, the largest permissible nonzero count.

    Computes the cumulative P[Poisson(lam) >= k] exactly via
    P[Poisson(lam) >= k] = 1 - sum_{j=0}^{k-1} lam^j e^{-lam} / j!.
    """
    # recurrence for P[N = j]: p_j = p_{j-1} * lam / j, p_0 = e^{-lam}
    log_p0 = -lam
    p0 = math.exp(log_p0)
    cdf = p0  # P[N <= 0]
    if 1 - cdf <= p_star:
        return 0  # k = 1 suffices, so N_max = 0
    pj = p0
    for j in range(1, k_search_max + 1):
        pj *= lam / j
        cdf += pj
        tail = max(1 - cdf, 0.0)
        if tail <= p_star:
            return j  # at k = j+1 the tail is <= p_star, hence N_max = j
    return k_search_max


def _M_voc_from_bound(sigma_dist: float, T_sigma: float = 11.0,
                      ceil_bound: bool = True) -> int:
    """Safe vocabulary radius via round-half-up tight bound.

    We use ceil(T_sigma * sigma) as a safe radius; the true tight bound
    comes from the l_infty bound on z and is typically a tight match
    for the Gaussian regime.
    """
    M = math.ceil(T_sigma * sigma_dist) if ceil_bound else round(T_sigma * sigma_dist)
    return max(M, 1)  # at least 1 (to allow +/-1 nonzero symbols)


def poisson_reservation(n_b: int, sigma_dist: float, t: int, p_star: float,
                        M_voc: int) -> Tuple[float, int, float, float]:
    """Poisson-model reservation (for very small sigma_dist / very small
    lambda, where the Gauss CLT does not apply).

    Steps:
      1. From the theoretical PMF take q = P[nonzero] = 1 - p(0);
         lambda = n_b * q.
      2. N_max = poisson_tail_kmax(lambda, p_star).
      3. Quantization: g(0), g(+/-1) determine ell(0), ell(+/-1).
      4. L_bits ~ const + ell(+/-1) * N, with const = 32 + n_b * ell(0).
      5. R = ceil((const + ell(+/-1) * N_max) / 8).

    Returns (lambda, N_max, const_bits, R_bytes).
    """
    p, g = quantized_pmf(sigma_dist, t, M_voc)
    quanta = 2 ** t
    p0 = p.get(0, 0.0)
    q = 1 - p0
    lam = n_b * q
    N_max = poisson_tail_kmax(lam, p_star)
    ell_0  = -math.log2(g[0]   / quanta) if g[0] > 0 else 0.0
    ell_pm = -math.log2(g[1]   / quanta) if g.get(1, 0) > 0 else float("inf")
    const_bits = 32 + n_b * ell_0
    R_bits = const_bits + ell_pm * N_max
    R_bytes = math.ceil(R_bits / 8)
    return lam, N_max, const_bits, R_bytes


# ----------------------------------------------------------------------------
# Automatic Gauss vs Poisson selection.
# ----------------------------------------------------------------------------
LAMBDA_THRESHOLD = 30.0  # Poisson when lambda < 30, else Gauss CLT


def reservation_auto(n_b: int, sigma_dist: float, t: int, p_star: float,
                     M_voc: int) -> Dict:
    """Pick Gauss or Poisson automatically (assumes the PMF is the
    discrete Gaussian).

    Used for the z-hi context. Rule: lambda < 30 -> Poisson, else Gauss.
    """
    p, g = quantized_pmf(sigma_dist, t, M_voc)
    p0 = p.get(0, 0.0)
    q = 1 - p0
    lam = n_b * q
    E_ell, var_ell = rans_rate_and_var(p, g, t)

    if lam < LAMBDA_THRESHOLD:
        lam_p, N_max, const_bits, R_p = poisson_reservation(n_b, sigma_dist, t, p_star, M_voc)
        return {
            "model":     "poisson",
            "lambda":    lam_p,
            "N_max":     N_max,
            "const_bits": const_bits,
            "mu_b":      n_b * E_ell / 8.0 + 4.0,
            "sigma_b":   math.sqrt(n_b * var_ell) / 8.0,
            "R_bytes":   R_p,
            "ell_pm1":   -math.log2(g.get(1, 1) / (2**t)) if g.get(1, 0) > 0 else float("nan"),
        }
    else:
        mu_b, sigma_b, R_g = gauss_reservation(n_b, E_ell, var_ell, p_star)
        return {
            "model":   "gauss",
            "lambda":  lam,
            "mu_b":    mu_b,
            "sigma_b": sigma_b,
            "R_bytes": math.ceil(R_g),
        }


def reservation_hint(n_b: int, r: float, alpha_h: int, t: int, p_star: float,
                     M_voc: int) -> Dict:
    """Regime-aware reservation for the hint context:

    sigma_h >= 1.0 -> discrete Gaussian D_{Z, sigma_h} (Gauss regime).
    sigma_h <  1.0 -> bucket-crossing Bernoulli-like PMF (narrow regime).
    """
    p, g = quantized_pmf_hint(r, alpha_h, t, M_voc)
    p0 = p.get(0, 0.0)
    lam = n_b * (1 - p0)
    E_ell, var_ell = rans_rate_and_var(p, g, t)
    mu_b, sigma_b, R_g = gauss_reservation(n_b, E_ell, var_ell, p_star)
    sigma_h = 2.0 * r / alpha_h
    model = "gauss-hint-gaussian" if sigma_h >= GAUSS_REGIME_THRESH \
                                   else "gauss-hint-bernoulli"
    return {
        "model":   model,
        "lambda":  lam,
        "mu_b":    mu_b,
        "sigma_b": sigma_b,
        "R_bytes": math.ceil(R_g),
    }


def compute_two_stream_reservation(p: dict, p_rans_star: float = 2**-20) -> dict:
    """Two-stream (z-hi + hint) reservation, with the per-stream budget
    split evenly: p_k = p_rans_star / 2 (union bound).

    The vocabulary radii M_voc come from the configured 11r design bounds.
    """
    n        = p["n"]
    lenS     = p["lenS"]
    lenE     = p["lenE"]
    r        = p["r"]
    alpha_h  = p["alpha_h"]
    alpha_r  = p.get("alpha_r", 64)
    T_sigma  = p.get("T_sigma", 11)
    eta      = p.get("eta", 1)
    tau      = p.get("tau", 30)
    t        = p.get("t", 10)

    By = T_sigma * r

    p_per = p_rans_star / 2

    # Stream 1: z-hi (n*(lenS+1) coefficients, scale r/alpha_r)
    n_b1     = n * (lenS + 1)
    sigma_zhi = r / alpha_r
    M_voc_zhi = math.ceil((By + tau * eta) / alpha_r)
    res1 = reservation_auto(n_b1, sigma_zhi, t, p_per, M_voc_zhi)

    # Stream 2: hint (n*lenE coefficients). Uses the bucket-crossing
    # PMF in place of the D_{Z, 2r/alpha_h} discrete-Gaussian
    # approximation (which severely underestimates the nonzero rate
    # when sigma_h << 1).
    n_b2     = n * lenE
    sigma_h  = 2 * r / alpha_h  # nominal, for reporting only
    M_voc_h  = math.floor(2 * (By + tau * eta) / alpha_h) + 1
    M_voc_h  = max(M_voc_h, 1)
    res2 = reservation_hint(n_b2, r, alpha_h, t, p_per, M_voc_h)

    return {
        "stream_zhi":  res1,
        "stream_hint": res2,
        "M_voc_zhi":   M_voc_zhi,
        "M_voc_h":     M_voc_h,
        "R_total":     res1["R_bytes"] + res2["R_bytes"],
        "p_rans_star": p_rans_star,
    }


def compute_merged_reservation(p: dict, p_rans_star: float = 2**-30) -> dict:
    """SHUTTLE-NIST single-stream design: one rANS stream encoding
    HighBits(z^(0..lenS)) || hint. The table and capacity calculation use
    the coefficient-count-weighted mixture of the rounded-response and
    bucket-crossing PMFs.
    """
    n        = p["n"]
    lenS     = p["lenS"]
    lenE     = p["lenE"]
    r        = p["r"]
    alpha_1  = p["alpha_1"]
    alpha_h  = p["alpha_h"]
    alpha_r  = p.get("alpha_r", 64)
    T_sigma  = p.get("T_sigma", 11)
    eta      = p.get("eta", 1)
    tau      = p.get("tau", 30)
    t        = p.get("t", 10)

    if alpha_h != 2 * alpha_r:
        # This utility currently covers the deployed parameter family, in
        # which alpha_h = 2*alpha_r keeps the two constituent scales close.
        return None

    By = T_sigma * r
    alpha_0p = alpha_r // alpha_1
    sigma = r / alpha_r

    # M_voc is the maximum of the three configured 11r-based design radii.
    # Out-of-vocabulary probability is accounted for separately.
    M_voc_zi = math.ceil((By + tau * eta) / alpha_r)
    M_voc_z0 = math.ceil((math.ceil(By / alpha_1) + 1) / alpha_0p)
    M_voc_h  = (2 * (By + tau * eta)) // alpha_h + 1
    M_voc = max(M_voc_zi, M_voc_z0, M_voc_h, 1)

    p_z0 = rounded_response_pmf(r, alpha_0p, M_voc,
                                inner_step=alpha_1,
                                T_sigma=T_sigma)
    p_zi = rounded_response_pmf(r, alpha_r, M_voc,
                                T_sigma=T_sigma)
    p_h = hint_pmf_correct(r, alpha_h, M_voc, T_sigma=T_sigma)

    block_count = 1 + lenS + lenE
    p_mix = {
        s: (p_z0[s] + lenS * p_zi[s] + lenE * p_h[s]) / block_count
        for s in p_z0
    }
    g_mix = quantize_pmf(p_mix, t)

    E_z0, V_z0 = rans_rate_and_var(p_z0, g_mix, t)
    E_zi, V_zi = rans_rate_and_var(p_zi, g_mix, t)
    E_h, V_h = rans_rate_and_var(p_h, g_mix, t)

    n_b = n * block_count
    mu_b = n * (E_z0 + lenS * E_zi + lenE * E_h) / 8.0 + 4.0
    sigma_b = math.sqrt(n * (V_z0 + lenS * V_zi + lenE * V_h)) / 8.0
    R_bytes = math.ceil(mu_b - phi_inv(p_rans_star) * sigma_b)
    res = {
        "model": "gauss-mixture",
        "mu_b": mu_b,
        "sigma_b": sigma_b,
        "R_bytes": R_bytes,
    }

    return {
        "stream":     res,
        "M_voc":      M_voc,
        "n_b":        n_b,
        "sigma":      sigma,
        "pmf":        p_mix,
        "freqs":      g_mix,
        "R_bytes":    res["R_bytes"],
        "p_rans_star": p_rans_star,
    }


def compute_entropy_decomposition(p: dict, p_rans_star: float = 2**-30) -> dict:
    """SHUTTLE-NIST single-stream byte breakdown:

      mu_b         merged rANS body (including the single 4 B flush)
      lo_z0_B      n * log2(alpha_0') / 8           (bit-pack)
      lo_zi_B      n * lenS * log2(alpha_r) / 8     (bit-pack)
      len_prefix   2 B (one stream)
      seedC        32/48/64 B
      Delta        R - mu_b (single-stream tail margin)
      actual sig   expected + Delta
    """
    n        = p["n"]
    lenS     = p["lenS"]
    alpha_1  = p["alpha_1"]
    alpha_r  = p.get("alpha_r", 64)
    alpha_0p = alpha_r // alpha_1
    seedC_B  = p["seedC_B"]

    out = compute_merged_reservation(p, p_rans_star)
    if out is None:
        return None
    rs = out["stream"]

    mu_b      = rs["mu_b"]
    lo_z0_B   = n * math.log2(alpha_0p) / 8
    lo_zi_B   = n * lenS * math.log2(alpha_r) / 8
    len_pref  = 2               # one merged stream * 2 B
    expected_sig = (mu_b + lo_z0_B + lo_zi_B
                    + len_pref + seedC_B)
    delta_total = rs["R_bytes"] - mu_b
    actual_sig  = expected_sig + delta_total

    return {
        "mu_b":         mu_b,
        "lo_z0_B":      lo_z0_B,
        "lo_zi_B":      lo_zi_B,
        "len_prefix_B": len_pref,
        "seedC_B":      seedC_B,
        "expected_sig": expected_sig,
        "delta_total":  delta_total,
        "actual_sig":   actual_sig,
        "R":            rs["R_bytes"],
        "M_voc":        out["M_voc"],
        "n_b":          out["n_b"],
        "sigma":        out["sigma"],
    }


def print_entropy_table(p: dict, p_rans_star: float = 2**-30) -> None:
    """Print the per-component byte breakdown (single-stream design).

    If the parameter set does not satisfy alpha_h = 2*alpha_r, skip.
    """
    e = compute_entropy_decomposition(p, p_rans_star)
    if e is None:
        print(f"\n  (skip merged-stream table: alpha_h={p['alpha_h']} != 2*alpha_r=128)")
        return
    star = int(math.log2(p_rans_star))
    print(f"\n  Signature byte breakdown (single-stream), p_rans^* = 2^{star}:")
    print(f"  {'Component':<54} {'Value (B)':>10}")
    print(f"  {'-'*54} {'-'*10}")
    print(f"  {'mu_b: merged rANS body (with 4 B flush, Gauss)':<54} {e['mu_b']:>10.2f}")
    print(f"  {'LowBits z^(0) bit-pack (n*log2(alpha_0)/8)':<54} {e['lo_z0_B']:>10.2f}")
    print(f"  {'LowBits z^(1..lenS) bit-pack (n*lenS*log2(alpha_r)/8)':<54} {e['lo_zi_B']:>10.2f}")
    print(f"  {'Length prefix (2 B)':<54} {e['len_prefix_B']:>10}")
    print(f"  {'seedC':<54} {e['seedC_B']:>10}")
    print(f"  {'-'*54} {'-'*10}")
    print(f"  {'Expected sig length (with quantization, no Delta)':<54} {e['expected_sig']:>10.2f}")
    print(f"  {'-'*54} {'-'*10}")
    print(f"  {'+ Delta (single-stream tail margin)':<54} {e['delta_total']:>10.2f}")
    print(f"  {'-'*54} {'-'*10}")
    print(f"  {'= Actual sig (single-stream, recommended sizing)':<54} {e['actual_sig']:>10.2f}")


def print_reservation_table(p: dict) -> None:
    """Print the single-stream reservation at three p_ovf_star targets."""
    out_sample = compute_merged_reservation(p, 2**-30)
    if out_sample is None:
        return
    print(f"\n  rANS reservation (single merged stream):")
    print(f"  vocab radius: M_voc = {out_sample['M_voc']},  n_b = {out_sample['n_b']},  "
          f"sigma = {out_sample['sigma']:.4f}")
    print(f"  {'p_rans_star':<14} {'model':<8} {'mu_b':>8} {'sigma_b':>9} {'R(B)':>6}")
    print(f"  {'-'*14} {'-'*8} {'-'*8} {'-'*9} {'-'*6}")
    for p_star in (2**-20, 2**-30, 2**-40):
        out = compute_merged_reservation(p, p_star)
        rs = out["stream"]
        print(f"  2^{int(math.log2(p_star)):<12} {rs['model']:<8} "
              f"{rs['mu_b']:>8.2f} {rs['sigma_b']:>9.2f} {rs['R_bytes']:>6}")


def print_one_mode(name: str, p: dict) -> None:
    res = compute_shannon_lower_bound(p)

    print(f"\n[{name}]")
    print(f"  parameters: n={p['n']}, lenS={p['lenS']}, lenE={p['lenE']}, "
          f"r={p['r']}, alpha_1={p['alpha_1']}, alpha_h={p['alpha_h']}")
    print(f"  sigma_z0 = r/alpha_1     = {res['sigma_z0']:.4f}")
    print(f"  sigma_h  = 2r/alpha_h    = {res['sigma_h']:.4f}")
    print()
    print(f"  Shannon entropy lower bound (no flush, no length prefix, no quantization; with seedC):")
    print(f"    H(z^(0))   = {res['H_z0']:+7.4f} bit/coef  -> "
          f"n*H/8        = {res['z0_B']:+8.2f} B")
    print(f"    H(z^(i))   = {res['H_zi']:+7.4f} bit/coef  -> "
          f"n*lenS*H/8   = {res['z_main_B']:+8.2f} B")
    print(f"    H(hint)    = {res['H_h_physical']:+7.4f} bit/coef  (clipped >=0; "
          f"continuous = {res['H_h_continuous']:+7.4f})")
    print(f"      discrete   {res['H_h_discrete']:+7.4f} bit/coef  -> "
          f"n*lenE*H/8   = {res['hint_B_discrete']:+8.2f} B")
    print(f"    seedC                                     = "
          f"{res['seedC_B']:+8.2f} B")
    print(f"    Shannon lower bound (clipped + seedC)     ="
          f" {res['total_B_physical']:+8.2f} B")
    # rANS-implementation byte breakdown (with flush, length prefix,
    # quantization and seedC included).
    print_entropy_table(p, p_rans_star=2**-30)

    # Reservation comparison across three p_rans^* targets.
    print_reservation_table(p)


def main() -> None:
    print("=" * 78)
    print("SHUTTLE-NIST signature size: Shannon lower bound + rANS reservation")
    print("=" * 78)
    print()

    for name, p in SHUTTLE_NIST_PARAMS.items():
        print_one_mode(name, p)

    print()
    print("=" * 78)
    print("SHUTTLE-NIST Shannon lower bound summary")
    print("=" * 78)
    print()
    print(f"  {'Parameter set':<22}{'Shannon lower (B)':>20}{'Discrete exact (B)':>20}")
    print(f"  {'-'*22}{'-'*20}{'-'*20}")
    for name, p in SHUTTLE_NIST_PARAMS.items():
        res = compute_shannon_lower_bound(p)
        print(f"  {name:<22}{res['total_B_physical']:20.2f}"
              f"{res['total_B_discrete']:20.2f}")


if __name__ == "__main__":
    main()
