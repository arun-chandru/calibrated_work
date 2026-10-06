"""Full-Gibbs smooth calibrated-wall certificate; standard library only.

Uses the independently audited Decimal interval primitives and exact step-wall
initial eigenvalue brackets from verify_uv_step_wall_intervals.py. The smooth
final wall is not spatially discretized as an eigenproblem. Exact rational
Bernstein bounds and interval rectangle integration enclose its matrix
elements; Rayleigh/Temple enclose its ground energy. All omitted Gibbs states
are bounded explicitly. See PHY_UV_SMOOTH_CALIBRATED_CERTIFICATE_2026_10_06.md.
"""

from decimal import Decimal as D, ROUND_CEILING
from fractions import Fraction as F
from math import comb, ceil
from pathlib import Path
import argparse
import json

from verify_uv_step_wall_intervals import I, PI, eigenpair


R = I('0.1')
BETA = I(1)
COLLAR = I('0.5')
VLOW = I('0.15')
LENGTH = I(4)


def fi(value):
    """Outward interval for an exact rational."""
    value = F(value)
    return I(value.numerator) / value.denominator


def mag(value):
    return max(abs(value.lo), abs(value.hi))


def upper(value):
    return I(value.hi)


def nonnegative(value):
    assert value.hi >= 0
    return I(max(D(0), value.lo), value.hi)


def degree9_bernstein(power):
    power = list(power) + [F(0)] * (10 - len(power))
    return [sum((power[j] * F(comb(k, j), comb(9, j))
                 for j in range(k + 1)), F(0)) for k in range(10)]


def split_bernstein(coefficients):
    row = list(coefficients)
    left, right = [row[0]], [row[-1]]
    while len(row) > 1:
        row = [(row[j] + row[j + 1]) / 2 for j in range(len(row) - 1)]
        left.append(row[0])
        right.append(row[-1])
    return left, list(reversed(right))


def subdivide(coefficients, depth):
    cells = [coefficients]
    for _ in range(depth):
        cells = [part for cell in cells for part in split_bernstein(cell)]
    return cells


def polynomial_data(depth):
    # x=z/4; b=(3/5)*4^4*x^4*(1-x)^5. All derivatives below are wrt z.
    power = [F(0)] * 10
    for j in range(6):
        power[4+j] = F(768, 5) * (-1)**j * comb(5, j)
    derivatives = [power]
    for _ in range(4):
        previous = derivatives[-1]
        derivatives.append([F(j+1, 4)*previous[j+1]
                            for j in range(len(previous)-1)])
    padded = [list(p) + [F(0)]*(10-len(p)) for p in derivatives]
    # (1+b)(V1-v)=1-v+(1-v)b+b''-2b'.
    lower = [(1-F(3, 20))*padded[0][j]+padded[2][j]-2*padded[1][j]
             for j in range(10)]
    lower[0] += 1-F(3, 20)
    # (1+b)(2-V1)=1+b-b''+2b'.
    higher = [padded[0][j]-padded[2][j]+2*padded[1][j]
              for j in range(10)]
    higher[0] += 1
    low_cells = subdivide(degree9_bernstein(lower), 6)
    high_cells = subdivide(degree9_bernstein(higher), 6)
    low_min = min(min(c) for c in low_cells)
    high_min = min(min(c) for c in high_cells)
    assert low_min > 0 and high_min > 0
    cell_coeffs = [subdivide(degree9_bernstein(p), depth) for p in derivatives]
    intervals = []
    for index in range(2**depth):
        intervals.append([I(fi(min(cells[index])).lo,
                            fi(max(cells[index])).hi) for cells in cell_coeffs])
    return intervals, {
        'b_power_coefficients_in_x': [str(x) for x in power],
        'endpoint_potential_lower': '0.15',
        'endpoint_potential_upper': '2',
        'positivity_dyadic_cells': 64,
        'lower_polynomial_min_Bernstein_coefficient_exact': str(low_min),
        'upper_polynomial_min_Bernstein_coefficient_exact': str(high_min),
        'path_bound_reason': 'V_s(z) is monotone in s at each z; 0<=s<=1',
    }


def geometric_gaussian_tail(c, start):
    return (-c*start**2).exp() / (1-(-c*(2*start+1)).exp())


def weighted_gaussian_tail(c, start):
    """Enclose sum_{n>=start} n^2 exp(-c n^2) from above."""
    x = (-c*(2*start+1)).exp()
    return (-c*start**2).exp() * (
        start**2/(1-x) + 2*start*x/(1-x)**2 + x*(1+x)/(1-x)**3)


def ceil_upper(value):
    return int(value.hi.to_integral_value(rounding=ROUND_CEILING))


def rectangle_integrals(ground, cells):
    count = len(cells)
    dz = LENGTH/count
    integrals = {key: I(0) for key in
                 ('mean', 'second', 'ramp_second', 'norm_change',
                  'overlap_change', 'residual_cross', 'residual_square',
                  'fourth', 'asymptotic_energy', 'asymptotic_second')}
    norms = {key: D(0) for key in ('F', 'Fprime', 'Fsecond', 'g')}
    for j, derivs in enumerate(cells):
        b, bp, bpp, bppp, bpppp = derivs
        b = nonnegative(b)  # b>=0 analytically.
        denominator = 1+b
        g = bpp-2*bp
        gp = bppp-2*bpp
        gpp = bpppp-2*bppp
        diff = g/denominator
        diffp = gp/denominator-g*bp/denominator**2
        diffpp = (gpp/denominator-2*gp*bp/denominator**2
                  -g*bpp/denominator**2+2*g*bp**2/denominator**3)
        for key, value in zip(norms, (diff, diffp, diffpp, g)):
            norms[key] = max(norms[key], mag(value))
        z = I((j*dz).lo, ((j+1)*dz).hi)
        phase = ground['alpha']*(COLLAR-R*z)
        ep, em = phase.exp(), (-phase).exp()
        psi = ground['tail']*(ep-em)/2
        psip = -ground['tail']*ground['alpha']*(ep+em)/2
        # Stable exact expression for psip+psi/r; avoid interval cancellation
        # between two O(psi/r) terms whose difference is much smaller.
        correction = ground['tail']*((1/R-ground['alpha'])*ep
                                     -(1/R+ground['alpha'])*em)/2
        psi2 = psi**2
        residual = -2*bp/R*correction
        fourth_vector = (-diffpp*psi/R**4-2*diffp*psip/R**3
                         +diff**2*psi/R**4)
        physical_measure = R*dz
        values = {
            # Integration by parts, b'(0)=b'(4)=0, gives this identical mean.
            'mean': psi2*bp**2/(R**2*denominator**2)
                    -2*psi*correction*bp/(R*denominator),
            'second': psi2*diff**2/R**4,
            'ramp_second': psi2*g**2/R**4,
            'norm_change': psi2*(2*b+b**2),
            'overlap_change': psi2*b,
            'residual_cross': psi*denominator*residual,
            'residual_square': residual**2,
            'fourth': fourth_vector**2,
        }
        for key, value in values.items():
            integrals[key] = integrals[key]+physical_measure*value
        exponential = (-2*z).exp()
        integrals['asymptotic_energy'] += dz*exponential*bp**2/denominator**2
        integrals['asymptotic_second'] += dz*exponential*diff**2
    return integrals, norms


def benchmark(depth):
    cells, profile = polynomial_data(depth)
    initial0 = eigenpair(1, R, I(0))
    initial1 = eigenpair(2, R, I(0))
    integrals, norms = rectangle_integrals(initial0, cells)
    e0 = initial0['e']
    weight = (-BETA*e0).exp()
    weight2 = (-BETA*initial1['e']).exp()
    c_outer = BETA*PI**2/I('1.5')**2
    # Exact root for label 2, then every initial label >=3 by min-max.
    initial_tail = weight2+geometric_gaussian_tail(c_outer, 3)
    initial_energy_tail = (initial1['e']*weight2
                           +PI**2*weighted_gaussian_tail(c_outer, 3))
    z0 = weight+I(0, initial_tail.hi)
    initial_excited_probability = upper(initial_tail/z0)
    initial_energy = (weight*e0+I(0, initial_energy_tail.hi))/z0

    dmax = I(norms['F'])/R**2
    dpmax = I(norms['Fprime'])/R**3
    dppmax = I(norms['Fsecond'])/R**4
    gmax = I(norms['g'])/R**2
    mean_tail = upper(dmax*initial_tail)
    mean = (weight*integrals['mean']+I(-mean_tail.hi, mean_tail.hi))/z0
    second = (weight*integrals['second']+I(0, (dmax**2*initial_tail).hi))/z0
    mg2 = (weight*integrals['ramp_second']
           +I(0, (gmax**2*initial_tail).hi))/z0
    fourth_tail = 3*((dppmax**2+dmax**4)*initial_tail
                     +4*dpmax**2*initial_energy_tail)
    fourth = (weight*integrals['fourth']+I(0, fourth_tail.hi))/z0
    variance = nonnegative(second-mean**2)

    # Direct-sum Neumann interface lower comparison. Its first eigenvalue is
    # the interior DN ground level; the next is min(interior DN2, exterior ND1).
    dn2 = (I('1.5')*PI)**2
    nd1 = VLOW/R**2+(PI/(2*COLLAR))**2
    gap_threshold = I(min(dn2.lo, nd1.lo), min(dn2.hi, nd1.hi))
    final_heat_tail = ((-BETA*dn2).exp()/(1-(-4*BETA*PI**2).exp())
                       +(-BETA*nd1).exp()
                       /(1-(-2*BETA*(PI/COLLAR)**2).exp()))
    norm2 = 1+integrals['norm_change']
    shift = integrals['residual_cross']/norm2
    rayleigh = e0+shift
    residual_variance = nonnegative(integrals['residual_square']/norm2-shift**2)
    assert rayleigh.hi < gap_threshold.lo
    temple = rayleigh-residual_variance/(gap_threshold-rayleigh)
    e1 = I(temple.lo, rayleigh.hi)
    z1 = (-BETA*e1).exp()+I(0, final_heat_tail.hi)
    free = -(z1.ln()-z0.ln())/BETA
    dissipated = mean-free
    assert dissipated.lo > 0

    # Complete ground off-label lower bound from Cauchy applied to H1-e1.
    excitation_mean = e0+integrals['mean']-e1
    excitation_second = (integrals['second']
                         -2*(e1-e0)*integrals['mean']+(e1-e0)**2)
    assert excitation_mean.lo > 0 and excitation_second.lo > 0
    ground_off_lower = (I(excitation_mean.lo)**2/upper(excitation_second)).lo
    # A trial-vector angle upper bound uses only the verified residual and gap.
    trial_overlap = (1+integrals['overlap_change'])/norm2.sqrt()
    trial_off_upper = max(D(0), (1-trial_overlap**2).hi)
    residual_angle = residual_variance.sqrt()/(gap_threshold-rayleigh)
    ground_off_upper = min(D(1), ((I(trial_off_upper).sqrt()
                                  +upper(trial_overlap)*residual_angle)**2).hi)
    off_lower = (weight*I(ground_off_lower)/z0).lo
    off_upper = min(D(1), ((weight*I(ground_off_upper)
                           +I(0, initial_tail.hi))/z0).hi)
    off = I(off_lower, off_upper)
    ratio = R**2*variance/dissipated
    asymptotic_ratio = (integrals['asymptotic_second']
                        /integrals['asymptotic_energy'])

    result = {
        'r': '0.1', 'beta': '1', 'outer_length': '1.5',
        'profile': profile,
        'integration': {'method': 'exact rational Bernstein ranges + interval rectangles',
                        'dyadic_cells': len(cells),
                        'ground_integrals': {k: v.pair() for k, v in integrals.items()}},
        'initial_ground_energy': e0.pair(),
        'initial_second_energy': initial1['e'].pair(),
        'initial_partition': z0.pair(),
        'initial_excited_heat_upper': str(initial_tail.hi),
        'initial_excited_probability_upper': str(initial_excited_probability.hi),
        'initial_excited_energy_numerator_upper': str(initial_energy_tail.hi),
        'initial_internal_energy': initial_energy.pair(),
        'dimensionless_profile_norm_upper': {k: str(v) for k, v in norms.items()},
        'physical_difference_norm_upper': str(dmax.hi),
        'physical_difference_first_derivative_norm_upper': str(dpmax.hi),
        'physical_difference_second_derivative_norm_upper': str(dppmax.hi),
        'physical_ramp_envelope_norm_upper': str(gmax.hi),
        'final_second_eigenvalue_lower': str(gap_threshold.lo),
        'final_excited_heat_upper': str(final_heat_tail.hi),
        'trial_norm_squared': norm2.pair(),
        'trial_rayleigh': rayleigh.pair(),
        'trial_residual_variance': residual_variance.pair(),
        'final_ground_energy': e1.pair(),
        'final_ground_relative_shift': ((e1-e0)/e0).pair(),
        'final_partition': z1.pair(),
        'free_energy_change': free.pair(),
        'mean_work': mean.pair(),
        'second_work_moment': second.pair(),
        'work_variance': variance.pair(),
        'fourth_work_moment': fourth.pair(),
        'dissipated_work': dissipated.pair(),
        'off_label_probability': off.pair(),
        'finite_r_pressure_ratio': ratio.pair(),
        'asymptotic_pressure_ratio': asymptotic_ratio.pair(),
        'ramp_envelope_second_moment': mg2.pair(),
    }

    ramps = []
    for taustr in ('0.0005', '0.001', '0.002'):
        tau = I(taustr)
        work_error = mg2*(tau/3+gmax*tau**2/20)
        norm_error = tau*mg2.sqrt()/2
        lower_root = max(D(0), (off.sqrt()-norm_error).lo)
        upper_root = (off.sqrt()+norm_error).hi
        probability = I((I(lower_root)**2).lo,
                        min(D(1), (I(upper_root)**2).hi))
        diss = dissipated+I(-work_error.hi, work_error.hi)
        moment_upper = upper(mg2*(1+gmax*tau/6)**2)
        assert probability.lo > 0 and diss.lo > 0
        # UV-A: finite-band energy readout, error <=0.5 at each resolved energy.
        # Accept initial report in [7,10], final report in [20,2000].
        # The first gate selects the initial ground label exactly.
        low_true = gap_threshold-e0
        assert (e0-I('0.5')).lo > 7 and (e0+I('0.5')).hi < 10
        assert (initial1['e']-I('0.5')).lo > 10
        assert gap_threshold.lo > D('20.5')
        assert (e1+I('0.5')).hi < 20
        high_true = I('1999.5')-e0
        excluded_upper_tail = upper(moment_upper/high_true**2)
        bin_lower = (probability-initial_excited_probability-excluded_upper_tail).lo
        bin_upper = probability.hi
        assert bin_lower > 0
        detection_trials = ceil_upper(I('0.05').ln()/(1-I(bin_lower)).ln())
        frequency_trials = ceil_upper(1/(I('0.05')*I('0.1')**2*I(bin_lower)))
        mean_trials = ceil_upper(moment_upper/(I('0.05')*(I('0.1')*I(diss.lo))**2))
        ramps.append({
            'tau': taustr, 'tau_over_r_squared': (tau/R**2).pair(),
            'dissipated_work': diss.pair(),
            'mean_work_bridge_error_upper': str(work_error.hi),
            'off_label_probability': probability.pair(),
            'second_work_moment_upper': str(moment_upper.hi),
            'mean_work_10percent_relative_95percent_Chebyshev_trials': mean_trials,
            'bounded_readout': {
                'reported_initial_energy_gate': ['7', '10'],
                'reported_final_energy_band': ['20', '2000'],
                'positive_reported_work_range_for_accepted_events': ['10', '1993'],
                'max_absolute_error_each_energy': '0.5',
                'max_absolute_error_work': '1',
                'instrument_contract': 'ideal TPM plus classical reporting, or nondemolition diagonal Kraus operators',
                'guaranteed_numerical_report_unit_efficiency': 'initial ground energy and final true energies [gamma,1999.5]; no in-band dropout',
                'overflow_flag_above_true_final_energy': '2000.5',
                'true_initial_ground_off_work_lower': str(low_true.lo),
                'true_final_energy_upper_for_guaranteed_inclusion': '1999.5',
                'true_work_threshold_used_for_tail_bound': high_true.pair(),
                'discarded_high_work_probability_upper': str(excluded_upper_tail.hi),
                'reported_bin_probability': [str(bin_lower), str(bin_upper)],
                'trials_for_at_least_one_event_95percent': detection_trials,
                'trials_for_10percent_relative_frequency_error_95percent_Chebyshev': frequency_trials,
            },
        })
    result['continuously_calibrated_ramps'] = ramps

    # Sudden ratio sampling: iid reported work with |report-W|<=1e-4.
    readout_error = I('0.0001')
    statistical_mean_error = I('0.001')
    statistical_second_error = I('0.5')
    confidence_failure = I('0.05')
    reported_second = (second.sqrt()+readout_error)**2
    reported_fourth = (fourth.sqrt().sqrt()+readout_error)**4
    n_mean = ceil_upper(2*reported_second/(confidence_failure*statistical_mean_error**2))
    n_second = ceil_upper(2*reported_fourth/(confidence_failure*statistical_second_error**2))
    raw_second_bias = 2*readout_error*second.sqrt()+readout_error**2
    eps_m = statistical_mean_error+readout_error
    eps_s = statistical_second_error+raw_second_bias
    possible_mean = mean+I(-eps_m.hi, eps_m.hi)
    possible_second = second+I(-eps_s.hi, eps_s.hi)
    free_mid = (free.lo+free.hi)/2
    possible_diss = possible_mean-I(free_mid)
    possible_variance = possible_second-possible_mean**2
    assert possible_diss.lo > 0 and possible_variance.lo > 0
    estimator_range = R**2*possible_variance/possible_diss
    result['sudden_ratio_readout'] = {
        'max_absolute_error_each_energy': '0.00005',
        'max_absolute_error_work': str(readout_error.hi),
        'statistical_mean_absolute_error': str(statistical_mean_error.hi),
        'statistical_raw_second_moment_absolute_error': str(statistical_second_error.hi),
        'mean_reporting_bias_upper': str(readout_error.hi),
        'raw_second_reporting_bias_upper': str(raw_second_bias.hi),
        'independent_preparations_sufficient': max(n_mean, n_second),
        'mean_constraint_trials': n_mean,
        'raw_second_constraint_trials': n_second,
        'simultaneous_confidence_at_least': '0.95',
        'free_energy_subtraction_value': str(free_mid),
        'free_energy_subtraction_systematic_error_upper': str(mag(free-I(free_mid))),
        'ratio_estimator_guaranteed_range_on_confidence_event': estimator_range.pair(),
        'scope': 'sudden TPM only; finite-time ratio is not inferred',
    }
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--depth', type=int, default=12,
                        help='2^depth rigorous integration cells; default 4096')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if not 8 <= args.depth <= 16:
        raise SystemExit('Use an integration depth between 8 and 16.')
    result = benchmark(args.depth)
    document = {'method': 'Decimal outward intervals, exact Bernstein bounds, Temple and full Gibbs tails',
                'precision': 58, 'pi_enclosure': PI.pair(), 'certificate': result}
    payload = json.dumps(document, indent=2)
    if args.output:
        args.output.write_text(payload+'\n', encoding='utf-8')
    print(payload)
