"""Certified finite-r pressure ratios for the SAME smooth wall and six states.

Standard library only; positive canonical N=2 configuration sums retain all
tails. Two dressed final-state residuals give Temple/Kato--Temple intervals.
No spectrum mesh, alternating partition recursion, or finite final-state TPM
cutoff occurs. See PHY_UV_PREPARATION_COMPARISON_2026_10_06.md.
"""

from decimal import Decimal as D
from pathlib import Path
import argparse
import json

from verify_uv_step_wall_intervals import I, PI, eigenpair
from verify_uv_smooth_calibrated_intervals import (
    polynomial_data, nonnegative, mag, geometric_gaussian_tail,
)

R = I('0.025')
COLLAR = I('0.5')
LENGTH = I(4)
OUTER = I('1.5')
VLOW = I('0.15')


def integrate(initial, cells):
    """Exact Bernstein ranges and outward rectangle quadrature."""
    keys = ('mean', 'second', 'norm_change', 'residual_cross', 'residual_square')
    records = [{key: I(0) for key in keys} for _ in initial]
    cross, energy, second = I(0), I(0), I(0)
    dmax = D(0)
    dz = LENGTH/len(cells)
    for j, derivs in enumerate(cells):
        b, bp, bpp = derivs[:3]
        b = nonnegative(b)
        den = 1+b
        diff = (bpp-2*bp)/den
        dmax = max(dmax, mag(diff))
        z = I((j*dz).lo, ((j+1)*dz).hi)
        measure = R*dz
        vectors = []
        for p, record in zip(initial, records):
            phase = p['alpha']*(COLLAR-R*z)
            ep, em = phase.exp(), (-phase).exp()
            psi = p['tail']*(ep-em)/2
            correction = p['tail']*((1/R-p['alpha'])*ep
                                    -(1/R+p['alpha'])*em)/2
            residual = -2*bp/R*correction
            vectors.append((psi, correction))
            values = (
                psi**2*bp**2/(R**2*den**2)-2*psi*correction*bp/(R*den),
                psi**2*diff**2/R**4,
                psi**2*(2*b+b**2),
                psi*den*residual,
                residual**2,
            )
            for key, value in zip(keys, values):
                record[key] += measure*value
        p, cp = vectors[0]
        q, cq = vectors[1]
        cross += measure*(p*q*bp**2/(R**2*den**2)
                          -(p*cq+cp*q)*bp/(R*den))
        exponential = (-2*z).exp()
        energy += dz*exponential*bp**2/den**2
        second += dz*exponential*diff**2
    return records, cross, I(dmax)/R**2, second/energy


def enclose_final(initial, records, comparison):
    """First two final levels, including an isolated excited-state interval."""
    final, diagnostics = [], []
    for n, (p, record) in enumerate(zip(initial, records)):
        norm2 = 1+record['norm_change']
        shift = record['residual_cross']/norm2
        rho = p['e']+shift
        eps2 = nonnegative(record['residual_square']/norm2-shift**2)
        upper_separator = I(comparison[n+1]['e'].lo)
        assert rho.hi < upper_separator.lo
        lower = rho-eps2/(upper_separator-rho)
        if n == 0:
            # Ground-state Rayleigh principle and Temple inequality.
            upper = rho
            lower_separator = None
        else:
            # Min-max supplies e_2 >= comparison[1] > e_1^upper and
            # e_2 <= (2*pi)^2 < comparison[2] <= e_3.
            lower_separator = I(final[0].hi)
            assert comparison[1]['e'].lo > lower_separator.hi
            assert ((2*PI)**2).hi < upper_separator.lo
            assert rho.lo > lower_separator.hi
            assert eps2.hi < ((rho-lower_separator)*(upper_separator-rho)).lo
            upper = rho+eps2/(rho-lower_separator)
        enclosure = I(lower.lo, upper.hi)
        final.append(enclosure)
        diagnostics.append({
            'level_label': n+1,
            'trial_norm_squared': norm2.pair(),
            'rayleigh': rho.pair(),
            'residual_variance': eps2.pair(),
            'lower_separator': None if lower_separator is None else lower_separator.pair(),
            'upper_separator': upper_separator.pair(),
            'energy_enclosure': enclosure.pair(),
        })
    return final, diagnostics


def low_weights(energies, beta, statistics):
    one = [(-beta*e).exp() for e in energies]
    if statistics == 'one_particle':
        return one, one
    if statistics == 'bose_N2':
        return one, [one[0]**2, one[0]*one[1], one[1]**2]
    assert statistics == 'fermi_N2'
    return one, [one[0]*one[1]]


def tail_one(beta, higher_energies):
    # Labels 3 and 4 explicitly, every label >=5 by the free outer box.
    return (sum(((-beta*e).exp() for e in higher_energies), I(0))
            +geometric_gaussian_tail(beta*PI**2/OUTER**2, 5))


def partition(energies, higher, beta, statistics):
    one, weights = low_weights(energies, beta, statistics)
    tail = tail_one(beta, higher)
    if statistics != 'one_particle':
        # Sum over low/high unordered pairs plus high/high pairs. T^2
        # majorizes every allowed high/high Bose or Fermi configuration.
        tail = sum(one, I(0))*tail+tail**2
    z = sum(weights, I(0))+I(0, tail.hi)
    return z, weights, tail


def preparation(beta_string, statistics, initial, final, comparison, records,
                cross, dmax, asymptotic):
    beta = I(beta_string)
    e0 = [p['e'] for p in initial[:2]]
    z0, weights, tail0 = partition(e0, [p['e'] for p in initial[2:]],
                                  beta, statistics)
    z1, _, tail1 = partition(final, [p['e'] for p in comparison[2:]],
                            beta, statistics)
    d0, d1 = (record['mean'] for record in records)
    q0, q1 = (record['second'] for record in records)
    if statistics == 'one_particle':
        means, seconds, n = [d0, d1], [q0, q1], 1
    elif statistics == 'bose_N2':
        means = [2*d0, d0+d1, 2*d1]
        seconds = [2*q0+2*d0**2,
                   q0+q1+2*d0*d1+2*cross**2,
                   2*q1+2*d1**2]
        n = 2
    else:
        means, seconds = [d0+d1], [q0+q1+2*d0*d1-2*cross**2]
        n = 2
    mean_error = (n*dmax*tail0).hi
    second_error = ((n*dmax)**2*tail0).hi
    mean = (sum((w*m for w, m in zip(weights, means)), I(0))
            +I(-mean_error, mean_error))/z0
    raw_second = (sum((w*q for w, q in zip(weights, seconds)), I(0))
                  +I(0, second_error))/z0
    variance = raw_second-mean**2
    free = -(z1.ln()-z0.ln())/beta
    dissipated = mean-free
    assert variance.lo > 0 and dissipated.lo > 0
    ratio = R**2*variance/dissipated
    relative = ratio/asymptotic-1
    return {
        'statistics': statistics, 'particle_number': n, 'beta': beta_string,
        'initial_partition': z0.pair(), 'final_partition': z1.pair(),
        'initial_omitted_configuration_probability_upper': str((tail0/z0).hi),
        'final_omitted_configuration_probability_upper': str((tail1/z1).hi),
        'mean_work': mean.pair(), 'second_work_moment': raw_second.pair(),
        'work_variance': variance.pair(), 'free_energy_change': free.pair(),
        'dissipated_work': dissipated.pair(),
        'finite_r_pressure_ratio': ratio.pair(),
        'ratio_minus_common_limit': (ratio-asymptotic).pair(),
        'relative_difference_from_common_limit': relative.pair(),
        'absolute_relative_difference_upper': str(mag(relative)),
    }


def benchmark(depth):
    cells, profile = polynomial_data(depth)
    initial = [eigenpair(n, R, I(0)) for n in range(1, 5)]
    comparison = [eigenpair(n, R/VLOW.sqrt(), I(0)) for n in range(1, 5)]
    records, cross, dmax, asymptotic = integrate(initial[:2], cells)
    final, diagnostics = enclose_final(initial[:2], records, comparison)
    preparations = [preparation(beta, stat, initial, final, comparison, records,
                                cross, dmax, asymptotic)
                    for beta in ('0.5', '1')
                    for stat in ('one_particle', 'bose_N2', 'fermi_N2')]
    ratio_ranges = [I(*p['finite_r_pressure_ratio']) for p in preparations]
    common_envelope = I(min(p.lo for p in ratio_ranges), max(p.hi for p in ratio_ranges))
    worst = max(D(p['absolute_relative_difference_upper']) for p in preparations)
    return {
        'method': '58-digit outward Decimal intervals; exact Bernstein bounds; Temple/Kato-Temple; complete positive canonical tails',
        'r': str(R.lo), 'outer_length': str(OUTER.lo), 'profile': profile,
        'integration_dyadic_cells': len(cells),
        'initial_first_four_energies': [p['e'].pair() for p in initial],
        'lower_comparison_first_four_energies': [p['e'].pair() for p in comparison],
        'final_first_two_certification': diagnostics,
        'matrix_integrals_first_two_modes': [{k: v.pair() for k, v in record.items()} for record in records],
        'off_diagonal_D01': cross.pair(), 'physical_difference_norm_upper': str(dmax.hi),
        'common_asymptotic_pressure_ratio': asymptotic.pair(),
        'preparations': preparations,
        'all_preparations_finite_r_ratio_envelope': common_envelope.pair(),
        'all_preparations_relative_difference_from_common_limit_upper': str(worst),
        'pairwise_ratio_difference_upper': str((I(common_envelope.hi)-I(common_envelope.lo)).hi),
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--depth', type=int, default=12)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if not 8 <= args.depth <= 16:
        raise SystemExit('Use depth between 8 and 16.')
    result = benchmark(args.depth)
    payload = json.dumps(result, indent=2)
    if args.output:
        args.output.write_text(payload+'\n', encoding='utf-8')
    print(payload)
