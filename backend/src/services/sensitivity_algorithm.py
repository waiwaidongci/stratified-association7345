"""Fixed-margin swap sensitivity for sparse stratified tables.

A correction with k>0 moves one success from method B to method A and one
failure from method A to method B; k<0 reverses that. Every corrected table
therefore keeps its row totals (n) and its success/failure margins. At most
``swap_limit`` such pairs are allowed in each stratum.

Probabilities are handled with exact integer weights: the conditional pmf is
W(s)/D with D = prod C(n_i, r_i), so tail probabilities can be compared,
minimised and maximised without floating-point tie ambiguity. The common OR is
compared as an integer pair (numerator, denominator) for the same reason.
"""
import hashlib
import json
import math
from fractions import Fraction
from itertools import product

from scipy.special import log_ndtr

from .cmh_algorithm import validate_tables

MAX_COMBINATIONS = 200_000
MAX_SWAP_LIMIT = 500
MAX_SUPPORT = 1201
MAX_EXAMPLES = 12
METHOD_VERSION = 'sensitivity-fixed-margin-swap-v1'
INTERPRETATION = (
    '每层只在固定边际内互换成功/失败标签，不增加试片、不改变层总数或行列边际；'
    '枚举结果描述少量误判如何改变共同OR方向、CMH与均值距离双侧精确p，'
    '不是新的确证检验，也不改变原研究的冻结快照与结论。'
)
TAIL_RULE = (
    '对每个修正表，取S=Σa在固定边际条件分布中的位置，纳入所有 |s−E[S]| ≥ |S−E[S]| 的支持点概率；'
    '距离与并列用整数权重精确判定，不是Fisher概率排序双侧。'
)
NUMERIC = (
    '各层超几何分布以整数权重卷积一次，得到每个支持点的精确权重；'
    '尾概率为整数权重之和除以共同分母，极小值同时保留科学计数文本与log值。'
)


def _tables(tables):
    if not isinstance(tables, list):
        raise ValueError('分层表须为列表')
    return sorted(
        ({'code': t.get('code', str(i)), 'name': t.get('name', ''),
          'a': t['a'], 'b': t['b'], 'c': t['c'], 'd': t['d']}
         for i, t in enumerate(tables)),
        key=lambda t: t['code'],
    )


def plan_sensitivity(snapshot, swap_limit):
    """Validate the frozen tables and return the complete enumeration plan."""
    if isinstance(swap_limit, bool) or not isinstance(swap_limit, int) or not 0 <= swap_limit <= MAX_SWAP_LIMIT:
        raise ValueError(f'每层互换上限须为0..{MAX_SWAP_LIMIT}的整数')
    rows = validate_tables(_tables(snapshot['tables']))
    if not 2 <= len(rows) <= 8:
        raise ValueError('敏感性复核需要2..8个完整分层')
    per_layer = []
    combinations = 1
    for t in rows:
        # a+k, b-k, c-k, d+k stay non-negative and |k| <= swap_limit.
        lower = max(-swap_limit, -t['a'], -t['d'])
        upper = min(swap_limit, t['b'], t['c'])
        count = upper - lower + 1
        per_layer.append({'code': t['code'], 'a': t['a'], 'b': t['b'], 'c': t['c'], 'd': t['d'],
                          'lower': lower, 'upper': upper, 'count': count})
        combinations *= count
    if combinations > MAX_COMBINATIONS:
        raise ValueError(
            f'修正组合数{combinations}超过资源上限{MAX_COMBINATIONS}；'
            '请降低每层互换上限后重新冻结'
        )
    support_size = sum(min(t['row_a'], t['success_total']) - max(0, t['row_a'] - t['failure_total']) + 1
                       for t in rows) - len(rows) + 1
    if support_size > MAX_SUPPORT:
        raise ValueError(f'条件分布支持域{support_size}点超过资源上限{MAX_SUPPORT}；请减少层数或每层总数')
    if not any(t['variance_a'] > 0 for t in rows):
        raise ValueError('全部层无有效条件方差，CMH敏感性未定义')
    return {'rows': rows, 'per_layer': per_layer, 'combination_count': combinations,
            'combination_cap': MAX_COMBINATIONS, 'support_size': support_size}


def _integer_weights(rows):
    """Convolve per-stratum integer hypergeometric weights once."""
    weights = [1]
    offset = 0
    for t in rows:
        r, successes, failures = t['row_a'], t['success_total'], t['failure_total']
        lo, hi = max(0, r - failures), min(r, successes)
        layer = [math.comb(successes, a) * math.comb(failures, r - a) for a in range(lo, hi + 1)]
        combined = [0] * (len(weights) + len(layer) - 1)
        for i, w in enumerate(weights):
            if w:
                for j, v in enumerate(layer):
                    combined[i + j] += w * v
        weights = combined
        offset += lo
    return offset, weights


def _tail_numerators(rows, support, weights):
    """Integer tail numerator for every support point under the mean-distance rule."""
    scale = 1
    for t in rows:
        scale = math.lcm(scale, t['n'])
    mean_numerator = sum(t['row_a'] * t['success_total'] * (scale // t['n']) for t in rows)
    distances = [abs(int(s) * scale - mean_numerator) for s in support]
    order = sorted(range(len(support)), key=lambda i: distances[i], reverse=True)
    numerators = [0] * len(support)
    cumulative = 0
    i = 0
    while i < len(order):
        distance = distances[order[i]]
        j = i
        while j < len(order) and distances[order[j]] == distance:
            j += 1
        cumulative += sum(weights[k] for k in order[i:j])
        for k in range(i, j):
            numerators[order[k]] = cumulative
        i = j
    denominator = 1
    for t in rows:
        denominator *= math.comb(t['n'], t['row_a'])
    return numerators, denominator, mean_numerator, scale


def _direction(num, den):
    """Exact direction of the common OR from integer numerator/denominator."""
    if den == 0:
        return 'undefined' if num == 0 else 'higher'
    if num == den:
        return 'equal'
    return 'higher' if num > den else 'lower'


def _ratio(num, den):
    if den == 0:
        return {'value': None, 'status': 'unbounded' if num > 0 else 'undefined'}
    return {'value': num / den, 'status': 'finite'}


def _scientific(logp):
    """Render 10**logp without depending on the platform long-double range."""
    exponent = math.floor(logp / math.log(10))
    mantissa = math.exp(logp - exponent * math.log(10))
    if not 1 <= mantissa < 10:
        mantissa, exponent = mantissa / 10, exponent + 1
    return f'{mantissa:.8f}e{exponent:+03d}'


def _probability(numerator, denominator):
    numeric = numerator / denominator
    logp = math.log(numerator) - math.log(denominator) if numerator > 0 else None
    if logp is None:
        display = '0.e+00'
    else:
        display = _scientific(logp)
    return {
        'value': numeric if numeric > 0 else None,
        'log': logp,
        'display': display,
        'underflow': bool(numerator > 0 and numeric == 0),
    }


def _asymptotic(statistic):
    logp = float(math.log(2) + log_ndtr(-math.sqrt(statistic)))
    p = math.exp(logp)
    return {
        'value': p if p > 0 else None,
        'log': logp,
        'display': _scientific(logp),
        'underflow': p == 0,
    }


class Extreme:
    """Track a minimum or maximum with exact comparisons and bounded examples."""

    def __init__(self, high=False):
        self.high = high
        self.best = None
        self.ties = 0
        self.examples = []

    def offer(self, value, item):
        if self.best is None or (value > self.best if self.high else value < self.best):
            self.best, self.ties, self.examples = value, 1, [item]
        elif value == self.best:
            self.ties += 1
            if len(self.examples) < MAX_EXAMPLES:
                self.examples.append(item)
        else:
            return

    def report(self, value):
        return {'value': value, 'tie_count': self.ties,
                'examples': self.examples,
                'examples_truncated': self.ties > len(self.examples)}


def analyze_sensitivity(snapshot):
    alpha = snapshot['alpha']
    swap_limit = snapshot['swap_limit']
    plan = plan_sensitivity(snapshot, swap_limit)
    rows = plan['rows']
    if isinstance(alpha, bool) or not isinstance(alpha, (int, float)) or not math.isfinite(alpha) or not .001 <= alpha <= .1:
        raise ValueError('alpha须在.001..0.1之间')

    expected = math.fsum(t['expected_a'] for t in rows)
    variance = math.fsum(t['variance_a'] for t in rows)
    if variance <= 0:
        raise ValueError('全部层无有效条件方差，CMH敏感性未定义')

    support_offset, weights = _integer_weights(rows)
    support = list(range(support_offset, support_offset + len(weights)))
    if len(support) > MAX_SUPPORT:
        raise ValueError(f'卷积支持域超过{MAX_SUPPORT}个点的资源限额')
    numerators, denominator, mean_numerator, scale = _tail_numerators(rows, support, weights)
    if sum(weights) != denominator:
        raise ValueError('整数权重与固定边际总质量不一致')

    alpha_threshold = Fraction(alpha) * denominator
    # A p within float noise of the frozen alpha cannot be ordered exactly as the
    # original longdouble comparison did; the original run's flag stays authoritative.
    alpha_tie_tolerance = alpha_threshold * Fraction(1, 10**12)
    base_sum = sum(t['a'] for t in rows)
    base_numerator = numerators[base_sum - support_offset]

    # One common integer scale lets every corrected OR be compared exactly.
    or_scale = 1
    for t in rows:
        or_scale = math.lcm(or_scale, t['n'])
    base_num = sum(t['a'] * t['d'] * (or_scale // t['n']) for t in rows)
    base_den = sum(t['b'] * t['c'] * (or_scale // t['n']) for t in rows)
    baseline_direction = _direction(base_num, base_den)
    base_statistic = (base_sum - expected) ** 2 / variance

    baseline = {
        'observed_sum_a': base_sum,
        'adjusted_or': _ratio(base_num, base_den),
        'direction': baseline_direction,
        'cmh': base_statistic,
        'asymptotic_p': _asymptotic(base_statistic),
        'exact_p': _probability(base_numerator, denominator),
        'below_alpha': bool(base_numerator <= alpha_threshold),
        'alpha_tie': bool(abs(base_numerator - alpha_threshold) <= alpha_tie_tolerance),
        'mean_numerator': mean_numerator,
        'mean_denominator': scale,
        'mean_distance': abs(base_sum * scale - mean_numerator) / scale,
    }

    direction_counts = {'higher': 0, 'lower': 0, 'equal': 0, 'undefined': 0}
    flip_count = changed_count = below_alpha_count = 0
    flip_examples = []
    p_low, p_high = Extreme(), Extreme(high=True)
    cmh_low, cmh_high = Extreme(), Extreme(high=True)

    def flip_sign(direction):
        return (1 if direction == 'higher' else -1 if direction == 'lower' else 0) * \
               (1 if baseline_direction == 'higher' else -1 if baseline_direction == 'lower' else 0)

    ranges = [range(x['lower'], x['upper'] + 1) for x in plan['per_layer']]
    ordinal = 0
    for shifts in product(*ranges):
        ordinal += 1
        observed = base_sum + sum(shifts)
        tail_numerator = numerators[observed - support_offset]
        num = sum((t['a'] + k) * (t['d'] + k) * (or_scale // t['n']) for t, k in zip(rows, shifts))
        den = sum((t['b'] - k) * (t['c'] - k) * (or_scale // t['n']) for t, k in zip(rows, shifts))
        direction = _direction(num, den)
        statistic = (observed - expected) ** 2 / variance
        below_alpha = tail_numerator <= alpha_threshold
        flipped = flip_sign(direction) == -1

        direction_counts[direction] += 1
        if direction != baseline_direction:
            changed_count += 1
        if flipped:
            flip_count += 1
        if below_alpha:
            below_alpha_count += 1

        item = {
            'ordinal': ordinal,
            'observed_sum_a': observed,
            'strata': [{'code': t['code'], 'k': k, 'a': t['a'] + k, 'b': t['b'] - k,
                        'c': t['c'] - k, 'd': t['d'] + k}
                       for t, k in zip(rows, shifts)],
            'adjusted_or': _ratio(num, den),
            'direction': direction,
            'cmh': statistic,
            'exact_p': _probability(tail_numerator, denominator),
            'below_alpha': below_alpha,
        }
        if flipped and len(flip_examples) < MAX_EXAMPLES:
            flip_examples.append(item)
        p_low.offer(tail_numerator, item)
        p_high.offer(tail_numerator, item)
        cmh_low.offer(statistic, item)
        cmh_high.offer(statistic, item)

    canonical = {'version': METHOD_VERSION, 'snapshot': snapshot}
    fingerprint = hashlib.sha256(json.dumps(
        canonical, ensure_ascii=False, sort_keys=True, separators=(',', ':'),
        allow_nan=False).encode()).hexdigest()

    return {
        'method_version': METHOD_VERSION,
        'fingerprint': fingerprint,
        'alpha': alpha,
        'swap_limit': swap_limit,
        'source': snapshot.get('source', {}),
        'capacity': {
            'combination_count': plan['combination_count'],
            'combination_cap': plan['combination_cap'],
            'per_layer': plan['per_layer'],
            'support_lower': support_offset,
            'support_upper': support[-1],
            'support_size': len(support),
            'tail_denominator': str(denominator),
            'example_cap': MAX_EXAMPLES,
        },
        'baseline': baseline,
        'baseline_note': '基线由同一冻结四格以整数权重重算，与原研究的longdouble结果只可能在浮点末位有差别；原研究记录不被改写。',
        'direction': {
            'baseline': baseline_direction,
            'flip_defined': baseline_direction in ('higher', 'lower'),
            'counts': direction_counts,
            'flip_count': flip_count,
            'changed_count': changed_count,
            'examples': flip_examples,
            'examples_truncated': flip_count > len(flip_examples),
        },
        'exact_p': {
            'min': p_low.report(_probability(p_low.best, denominator)),
            'max': p_high.report(_probability(p_high.best, denominator)),
        },
        'cmh': {
            'min': cmh_low.report(cmh_low.best),
            'max': cmh_high.report(cmh_high.best),
        },
        'significance': {
            'baseline_below_alpha': baseline['below_alpha'],
            'below_alpha_count': below_alpha_count,
            'above_alpha_count': plan['combination_count'] - below_alpha_count,
        },
        'tail_rule': TAIL_RULE,
        'numeric': NUMERIC,
        'interpretation': INTERPRETATION,
    }
