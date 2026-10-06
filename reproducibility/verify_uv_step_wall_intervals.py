"""Outward-rounded finite-r TPM benchmark for matched step walls.

Standard library only. Decimal arithmetic is outward rounded; elementary
sin/cos use Taylor's explicit remainder, pi uses Machin's identity, and the
full Gibbs tail is bounded by min-max and a geometric majorant. No spatial
mesh, truncated final-state sum, or numerical integration is used.
See PHY_UV_FINITE_PARAMETER_BENCHMARK_2026_10_06.md for mathematical scope.
"""
from decimal import Decimal as D, localcontext, getcontext, ROUND_FLOOR, ROUND_CEILING
import math
import json
import argparse
from pathlib import Path

PREC = 58
getcontext().prec = PREC


def rounded(fn, direction):
    with localcontext() as c:
        c.prec = PREC
        c.rounding = direction
        return +fn()


class I:
    def __init__(self, lo, hi=None):
        if isinstance(lo, I):
            self.lo, self.hi = lo.lo, lo.hi
        else:
            self.lo, self.hi = D(lo), D(lo if hi is None else hi)
        assert self.lo <= self.hi

    def __add__(self, other):
        o = I(other)
        return I(rounded(lambda: self.lo+o.lo, ROUND_FLOOR),
                 rounded(lambda: self.hi+o.hi, ROUND_CEILING))
    __radd__ = __add__

    def __neg__(self):
        return I(self.hi.copy_negate(), self.lo.copy_negate())

    def __sub__(self, other):
        return self + (-I(other))

    def __rsub__(self, other):
        return I(other) + (-self)

    def __mul__(self, other):
        o = I(other)
        pairs = [(a,b) for a in (self.lo,self.hi) for b in (o.lo,o.hi)]
        return I(min(rounded(lambda: a*b, ROUND_FLOOR) for a,b in pairs),
                 max(rounded(lambda: a*b, ROUND_CEILING) for a,b in pairs))
    __rmul__ = __mul__

    def __truediv__(self, other):
        o = I(other)
        assert not o.lo <= 0 <= o.hi
        inv = I(rounded(lambda: D(1)/o.hi, ROUND_FLOOR),
                rounded(lambda: D(1)/o.lo, ROUND_CEILING))
        return self*inv

    def __rtruediv__(self, other):
        return I(other)/self

    def __pow__(self, n):
        assert isinstance(n,int) and n >= 0
        if n == 0:
            return I(1)
        # Repeated products, not Decimal.power (whose C implementation has a
        # weaker general rounding guarantee). Endpoint integer powers use
        # a sign-aware outward multiplication chain.
        def endpoint(x):
            v=I(1)
            for _ in range(n):
                v=v*I(x)
            return v
        low_end, high_end=endpoint(self.lo),endpoint(self.hi)
        a,b,c,d=low_end.lo,high_end.lo,low_end.hi,high_end.hi
        low = D(0) if n % 2 == 0 and self.lo <= 0 <= self.hi else min(a,b)
        return I(low, max(c,d))

    def exp(self):
        # Decimal exp/ln/sqrt are correctly rounded; one ulp outward suffices.
        return I(self.lo.exp().next_minus(), self.hi.exp().next_plus())

    def ln(self):
        assert self.lo > 0
        return I(self.lo.ln().next_minus(), self.hi.ln().next_plus())

    def sqrt(self):
        assert self.lo >= 0
        return I(max(D(0),self.lo.sqrt().next_minus()), self.hi.sqrt().next_plus())

    def sinh(self):
        return (self.exp()-(-self).exp())/2

    def cosh(self):
        return (self.exp()+(-self).exp())/2

    def trig(self, cosine=False):
        # Lagrange remainder after the indicated Taylor polynomial.
        x2 = self*self
        term = I(1) if cosine else self
        total = term
        for j in range(1,81):
            den = (2*j-1)*(2*j) if cosine else (2*j)*(2*j+1)
            term = -term*x2/den
            total = total+term
        degree = 160 if cosine else 161
        mag = I(max(abs(self.lo),abs(self.hi)))
        error = (mag**(degree+1)/I(math.factorial(degree+1))).hi
        return total + I(-error,error)

    def sin(self):
        return self.trig()

    def cos(self):
        return self.trig(True)

    def pair(self):
        return [str(self.lo),str(self.hi)]


def atan_small(x):
    x = I(x)
    term, total = x, x
    for n in range(1,100):
        term = -term*x*x
        total = total+term/(2*n+1)
    error = (x**201/201).hi
    return total+I(-error,error)


PI = 16*atan_small(I(1)/5)-4*atan_small(I(1)/239)
R = I('1.5')
THETA = I('0.5')


def characteristic(k, length, height):
    alpha = (height-k*k).sqrt()
    t = alpha*(R-length)
    return k*(k*length).cos()/(k*length).sin()+alpha*t.cosh()/t.sinh()


def eigenpair(n,r,a):
    length = 1+a*r
    height = 1/((1-a)*r)**2
    # The nth sub-barrier root lies in this half period. Use verified signs
    # at interior rational points; monotonicity of k cot(kL)+alpha coth(alpha D)
    # identifies its eigenvalue label.
    left = ((I(n)-I('0.499'))*PI/length).hi
    right = ((I(n)-I('0.001'))*PI/length).lo
    assert right*right < height.lo
    assert characteristic(I(left),length,height).lo > 0
    assert characteristic(I(right),length,height).hi < 0
    for _ in range(150):
        mid = (left+right)/2
        f = characteristic(I(mid),length,height)
        if f.lo > 0:
            left = mid
        elif f.hi < 0:
            right = mid
        else:
            assert right-left < D('1e-32')
            break
    assert right-left < D('1e-32')
    k = I(left,right)
    alpha = (height-k*k).sqrt()
    span = R-length
    sh = (alpha*span).sinh()
    inner = length/2-(2*k*length).sin()/(4*k)
    wall = (2*alpha*span).sinh()/(4*alpha*sh**2)-span/(2*sh**2)
    amp = 1/(inner+(k*length).sin()**2*wall).sqrt()
    return dict(k=k,e=k*k,alpha=alpha,length=length,span=span,
                amp=amp,tail=amp*(k*length).sin()/sh)


def wall_integral(p,start,end):
    aa = p['alpha']
    # p.tail multiplies sinh(alpha*(R-y)).
    return p['tail']**2*((2*aa*(R-start)).sinh()/(4*aa)
                         -(2*aa*(R-end)).sinh()/(4*aa)-(end-start)/2)


def overlap(p,q):
    k,l,aa,bb = p['k'],q['k'],p['alpha'],q['alpha']
    inside = ((k-l).sin()/(2*(k-l))-(k+l).sin()/(2*(k+l))) * p['amp']*q['amp']
    def primitive(y):
        z = aa*(R-y)
        return (-aa*z.cosh()*(l*y).sin()-l*z.sinh()*(l*y).cos())/(aa*aa+l*l)
    middle = p['tail']*q['amp']*(primitive(q['length'])-primitive(I(1)))
    span = R-q['length']
    outer = p['tail']*q['tail']*((((aa+bb)*span).sinh()/(2*(aa+bb)))
                               -(((aa-bb)*span).sinh()/(2*(aa-bb))))
    return inside+middle+outer


def benchmark(rstr):
    r,a = I(rstr),I('0.5')
    count=3
    initial=[eigenpair(n,r,I(0)) for n in range(1,count+1)]
    final=[eigenpair(n,r,a) for n in range(1,count+1)]
    thermal0=[(-THETA*p['e']).exp() for p in initial]
    thermal1=[(-THETA*p['e']).exp() for p in final]
    heat_c=THETA*PI**2/R**2
    tail=(-heat_c*(count+1)**2).exp()/(1-(-heat_c*(2*count+3)).exp())
    z0=sum(thermal0,I(0))+I(0,tail.hi)
    z1=sum(thermal1,I(0))+I(0,tail.hi)
    dv=3/r**2
    mean=I(0)
    second=I(0)
    prob=I(0)
    diagonal_mean=I(0)
    for p,q,weight in zip(initial,final,thermal0):
        first=wall_integral(p,I(1),1+a*r)
        last=wall_integral(p,1+a*r,R)
        mean=mean+weight*(-first+3*last)/r**2
        second=second+weight*(first+9*last)/r**4
        overlap2=overlap(p,q)**2
        assert 0 < overlap2.lo <= overlap2.hi < 1
        prob=prob+weight*(1-overlap2)
        diagonal_mean=diagonal_mean+weight*overlap2*(q['e']-p['e'])
    mean=(mean+I(-(dv*tail).hi,(dv*tail).hi))/z0
    second=(second+I(0,(dv**2*tail).hi))/z0
    prob=(prob+I(0,tail.hi))/z0
    diagonal_mean=(diagonal_mean+I(-(dv*tail).hi,(dv*tail).hi))/z0
    free=-(z1.ln()-z0.ln())/THETA
    dissipated=mean-free
    variance=second-mean**2
    conditional=r**2*(mean-diagonal_mean)/prob
    out={
        'r':rstr,'theta':'0.5','outer_length':'1.5','computed_initial_levels':count,
        'ground_energy':initial[0]['e'].pair(),
        'ground_relative_shift':((final[0]['e']-initial[0]['e'])/initial[0]['e']).pair(),
        'free_energy_change':free.pair(),
        'mean_work':mean.pair(),'dissipated_work':dissipated.pair(),
        'off_label_probability':prob.pair(),'work_variance':variance.pair(),
        'conditional_scaled_off_label_work_mean':conditional.pair(),
        'initial_thermal_tail_probability_upper':str((tail/z0).hi),
        'mean_work_over_r':(mean/r).pair(),
        'probability_over_r_cubed':(prob/r**3).pair(),
        'r_times_variance':(r*variance).pair(),
    }
    assert dissipated.lo > 0 and prob.lo > 0
    return out


def ramp_bounds(record, taustr):
    r,tau=I(record['r']),I(taustr)
    mean=I(*record['mean_work'])
    variance=I(*record['work_variance'])
    sudden=I(*record['dissipated_work'])
    prob=I(*record['off_label_probability'])
    m2=variance+mean**2
    d=3/r**2
    error=m2*(tau/3+d*tau**2/20)
    norm_error=tau*m2.sqrt()/2
    lower=(prob.sqrt()-norm_error).lo
    upper=(prob.sqrt()+norm_error).hi
    p_lower=(I(max(D(0),lower))**2).lo
    p_upper=min(D(1),(I(upper)**2).hi)
    return dict(r=record['r'],tau=taustr,
                protocol='H(s)=H0+s(H1-H0); endpoint calibration only',
                second_sudden_moment=m2.pair(),
                mean_work_error_upper=str(error.hi),
                dissipated_work=(sudden+I(-error.hi,error.hi)).pair(),
                off_label_probability=[str(p_lower),str(p_upper)])


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--precision',type=int,default=58)
    parser.add_argument('--ramp-only',type=Path,
                        help='Read the saved interval certificate and verify finite-time bridge bounds.')
    args=parser.parse_args()
    if args.precision != 58:
        raise SystemExit('Precision is fixed at58 for the documented remainder budget.')
    if args.ramp_only:
        record=json.loads(args.ramp_only.read_text(encoding='utf-8'))['results'][0]
        assert record['r']=='0.1'
        print(json.dumps({'finite_r_linear_ramp_bounds':[
            ramp_bounds(record,t) for t in ('0.001','0.002','0.003')]},indent=2))
        raise SystemExit(0)
    results=[benchmark(r) for r in ('0.1','0.05','0.025')]
    print(json.dumps({'method':'Decimal outward intervals + exact roots + complete Gibbs tail',
                      'pi_enclosure':PI.pair(),'results':results},indent=2))
