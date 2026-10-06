"""UV-C: finite-family amplitude design with the frozen UV3 certificate engine.

The declared candidates are b_a=a*z^4*(1-z/4)^5, a in AMPLITUDES. All fixed
resources, time points, interval arithmetic and thermal tails are identical to
verify_uv_smooth_calibrated_intervals.py. Only its profile factory is injected
temporarily; the original binding is restored after every candidate. This
does not modify the frozen base source. Standard library only.
"""

import argparse
from decimal import Decimal
from fractions import Fraction as F
from math import comb
from pathlib import Path
import json

import verify_uv_smooth_calibrated_intervals as base


AMPLITUDES = ('0.3', '0.4', '0.5', '0.55', '0.6', '0.625', '0.64', '0.65')
TARGET_TAU = '0.001'


def profile_factory(amplitude):
    amplitude = F(amplitude)

    def polynomial_data(depth):
        power = [F(0)] * 10
        for j in range(6):
            power[4+j] = amplitude * 256 * (-1)**j * comb(5, j)
        derivatives = [power]
        for _ in range(4):
            previous = derivatives[-1]
            derivatives.append([F(j+1, 4)*previous[j+1]
                                for j in range(len(previous)-1)])
        padded = [list(p) + [F(0)]*(10-len(p)) for p in derivatives]
        lower = [F(17, 20)*padded[0][j]+padded[2][j]-2*padded[1][j]
                 for j in range(10)]
        lower[0] += F(17, 20)
        higher = [padded[0][j]-padded[2][j]+2*padded[1][j]
                  for j in range(10)]
        higher[0] += 1
        low_cells = base.subdivide(base.degree9_bernstein(lower), 6)
        high_cells = base.subdivide(base.degree9_bernstein(higher), 6)
        low_min = min(min(c) for c in low_cells)
        high_min = min(min(c) for c in high_cells)
        assert low_min > 0 and high_min > 0
        cell_coeffs = [base.subdivide(base.degree9_bernstein(p), depth)
                       for p in derivatives]
        intervals = []
        for index in range(2**depth):
            intervals.append([base.I(base.fi(min(cells[index])).lo,
                                     base.fi(max(cells[index])).hi)
                              for cells in cell_coeffs])
        return intervals, {
            'amplitude_exact': str(amplitude),
            'b_power_coefficients_in_x': [str(x) for x in power],
            'endpoint_potential_lower': '0.15',
            'endpoint_potential_upper': '2',
            'positivity_dyadic_cells': 64,
            'lower_polynomial_min_Bernstein_coefficient_exact': str(low_min),
            'upper_polynomial_min_Bernstein_coefficient_exact': str(high_min),
            'path_bound_reason': 'V_s(z) is monotone in s at each z; 0<=s<=1',
        }
    return polynomial_data


def evaluate_candidate(amplitude, depth):
    original_factory = base.polynomial_data
    try:
        base.polynomial_data = profile_factory(amplitude)
        return base.benchmark(depth)
    finally:
        base.polynomial_data = original_factory


def run_design(depth):
    results, comparison = [], []
    for amplitude in AMPLITUDES:
        certificate = evaluate_candidate(amplitude, depth)
        ramp = next(r for r in certificate['continuously_calibrated_ramps']
                    if r['tau'] == TARGET_TAU)
        readout = ramp['bounded_readout']
        comparison.append({
            'amplitude': amplitude,
            'dissipated_work_interval': ramp['dissipated_work'],
            'off_label_probability_interval': ramp['off_label_probability'],
            'second_moment_upper': ramp['second_work_moment_upper'],
            'accepted_bin_probability_interval': readout['reported_bin_probability'],
            '95percent_detection_preparations': readout['trials_for_at_least_one_event_95percent'],
            '10percent_relative_frequency_95percent_preparations': readout['trials_for_10percent_relative_frequency_error_95percent_Chebyshev'],
        })
        results.append({'amplitude': amplitude, 'certificate': certificate})
        print(json.dumps(comparison[-1]), flush=True)
    selected = max(comparison,
                   key=lambda row: Decimal(row['accepted_bin_probability_interval'][0]))
    baseline = next(row for row in comparison if row['amplitude'] == '0.6')
    p_selected = base.I(selected['accepted_bin_probability_interval'][0])
    p_baseline = base.I(baseline['accepted_bin_probability_interval'][0])
    reduction = 1-base.I(selected['95percent_detection_preparations'])/baseline['95percent_detection_preparations']
    return {
        'scope': 'certified guarantee optimization over a declared eight-element amplitude family; not global or actual-probability optimality',
        'resources': {
            'r': '0.1', 'beta': '1', 'tau': TARGET_TAU,
            'outer_length': '1.5', 'scaled_profile_support': ['0', '4'],
            'path_potential_bounds': ['0.15', '2'],
            'calibration': 'h_s=1 for every s',
            'initial_report_gate': ['7', '10'],
            'final_report_band': ['20', '2000'],
            'energy_error_max_each': '0.5',
            'unit_efficiency': 'initial ground and final true energy [gamma,1999.5]',
            'instrument': 'ideal TPM plus classical reporting, or nondemolition diagonal Kraus operators',
        },
        'candidate_amplitudes': list(AMPLITUDES),
        'integration_depth': depth,
        'objective': 'maximize the certified lower accepted-bin probability, equivalently minimize its sufficient 95percent detection budget',
        'comparison': comparison,
        'baseline_amplitude': '0.6',
        'selected_amplitude': selected['amplitude'],
        'selected_to_baseline_certified_lower_probability_ratio': (p_selected/p_baseline).pair(),
        'sufficient_detection_budget_fractional_reduction': reduction.pair(),
        'all_candidate_certificates': results,
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--depth', type=int, default=12)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if not 8 <= args.depth <= 16:
        raise SystemExit('Use an integration depth between 8 and 16.')
    result = run_design(args.depth)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({key: result[key] for key in (
        'selected_amplitude',
        'selected_to_baseline_certified_lower_probability_ratio',
        'sufficient_detection_budget_fractional_reduction')}, indent=2))
