"""Independent UV-C saved-record QA and fresh selected-candidate rebuild."""

from decimal import Decimal as D
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json

import verify_uv_bounded_design_intervals as design


def main():
    directory = Path(__file__).resolve().parent
    source = directory/'PHY_UV_BOUNDED_DESIGN_CERTIFICATE_2026_10_06.json'
    saved = json.loads(source.read_text(encoding='utf-8'))
    assert saved['candidate_amplitudes'] == list(design.AMPLITUDES)
    assert saved['integration_depth'] == 12
    assert len(saved['comparison']) == len(saved['all_candidate_certificates']) == 8
    rows = saved['comparison']
    for row, entry in zip(rows, saved['all_candidate_certificates']):
        assert row['amplitude'] == entry['amplitude']
        certificate = entry['certificate']
        assert (certificate['r'], certificate['beta'], certificate['outer_length']) == ('0.1', '1', '1.5')
        assert certificate['integration']['dyadic_cells'] == 4096
        profile = certificate['profile']
        assert F(profile['amplitude_exact']) == F(row['amplitude'])
        assert profile['endpoint_potential_lower'] == '0.15'
        assert profile['endpoint_potential_upper'] == '2'
        assert F(profile['lower_polynomial_min_Bernstein_coefficient_exact']) > 0
        assert F(profile['upper_polynomial_min_Bernstein_coefficient_exact']) > 0
        ramp = next(r for r in certificate['continuously_calibrated_ramps'] if r['tau'] == '0.001')
        readout = ramp['bounded_readout']
        assert readout['reported_initial_energy_gate'] == ['7', '10']
        assert readout['reported_final_energy_band'] == ['20', '2000']
        assert readout['max_absolute_error_each_energy'] == '0.5'
        assert row['accepted_bin_probability_interval'] == readout['reported_bin_probability']
        assert D(row['accepted_bin_probability_interval'][0]) > 0
        assert D(row['dissipated_work_interval'][0]) > 0
        assert row['95percent_detection_preparations'] == readout['trials_for_at_least_one_event_95percent']
        p = design.base.I(row['accepted_bin_probability_interval'][0])
        trials = row['95percent_detection_preparations']
        threshold = design.base.I('0.05').ln()
        assert (trials*(1-p).ln()).hi < threshold.lo
        assert ((trials-1)*(1-p).ln()).lo > threshold.hi
    selected = max(rows, key=lambda row: D(row['accepted_bin_probability_interval'][0]))
    assert saved['selected_amplitude'] == selected['amplitude'] == '0.65'
    baseline = next(row for row in rows if row['amplitude'] == '0.6')
    # This is an optimization of the certified guarantee; actual probability
    # superiority over the baseline is NOT implied by the overlapping ranges.
    assert D(selected['accepted_bin_probability_interval'][0]) < D(baseline['accepted_bin_probability_interval'][1])
    original_factory = design.base.polynomial_data
    fresh = design.evaluate_candidate(selected['amplitude'], 12)
    assert design.base.polynomial_data is original_factory
    selected_saved = next(entry['certificate'] for entry in saved['all_candidate_certificates']
                          if entry['amplitude'] == selected['amplitude'])
    assert fresh == selected_saved
    original = json.loads((directory/'PHY_UV_SMOOTH_CALIBRATED_CERTIFICATE_2026_10_06.json').read_text(encoding='utf-8'))['certificate']
    baseline_saved = next(entry['certificate'] for entry in saved['all_candidate_certificates'] if entry['amplitude'] == '0.6')
    baseline_copy = json.loads(json.dumps(baseline_saved))
    del baseline_copy['profile']['amplitude_exact']
    assert baseline_copy == original
    files = ('verify_uv_bounded_design_intervals.py', 'verify_uv_smooth_calibrated_intervals.py',
             'verify_uv_step_wall_intervals.py', source.name)
    result = {
        'status': 'PASS',
        'all_eight_saved_candidates_common_resource_checks': 'PASS',
        'all_eight_95percent_budget_ceilings_rechecked': 'PASS',
        'baseline_matches_original_smooth_certificate_except_amplitude_metadata': True,
        'fresh_selected_candidate_exactly_matches_saved': True,
        'source_profile_factory_restored': True,
        'selected_amplitude': selected['amplitude'],
        'integration_depth': 12,
        'actual_probability_superiority_over_baseline_certified': False,
        'reason': 'selected lower probability bound is below baseline upper bound',
        'sha256': {name: hashlib.sha256((directory/name).read_bytes()).hexdigest() for name in files},
        'fresh_selected_certificate': fresh,
    }
    target = directory/'PHY_UV_BOUNDED_DESIGN_INDEPENDENT_QA_2026_10_06.json'
    target.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k: v for k, v in result.items() if k != 'fresh_selected_certificate'}, indent=2))


if __name__ == '__main__':
    main()
