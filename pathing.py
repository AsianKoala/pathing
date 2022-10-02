from geometry import *

class DifferentiablePoint():
    def __init__(self, zero=0, first=0,  second=0, third=0):
        self.first = first
        self.zero = zero
        self.second = second
        self.third = third


class Quintic:
    def __init__(self, start: DifferentiablePoint, end: DifferentiablePoint) -> None:
        self.end = end
        self.start = start
        coeffs = [0]*6
        coeffs[5] = start.zero
        coeffs[4] = start.first
        coeffs[3] = start.second / 2.0
        coeffs[2] = -(20 * start.zero + 12 * start.first + 3 * start.second - 
                20 * end.zero + 8 * end.first - end.second) / 2.0
        coeffs[1] = (30 * start.zero + 16 * start.first + 3 * start.second
                - 30 * end.zero + 14 * end.first - 2 * end.second) / 2.0
        coeffs[0] = -(12 * start.zero + 6 * start.first + start.second
                - 12 * end.zero + 6 * end.first - end.second) / 2.0
        self.coeffs = coeffs

    def get(self, t: float) -> DifferentiablePoint:
        coeffs = self.coeffs
        return DifferentiablePoint(
            coeffs[0] * pow(t, 5) + coeffs[1] * pow(t, 4) + coeffs[2] * pow(t, 3) + coeffs[3] * pow(t, 2) + coeffs[4] * t + coeffs[5],
            5 * coeffs[0] * pow(t, 4) + 4 * coeffs[1] * pow(t, 3) + 3 * coeffs[2] * pow(t, 2) + 2 * coeffs[3] * t + coeffs[4],
            20 * coeffs[0] * pow(t, 3) + 12 * coeffs[1] * pow(t, 2) + 6 * coeffs[2] * t + 2 * coeffs[3])

    def __str__(self): 
        return "{}t^5 + {}t^4 + {}t^3 + {}t^2 + {}t + {}".format(round(self.coeffs[0]), round(self.coeffs[1]), 
            round(self.coeffs[2]), round(self.coeffs[3]), round(self.coeffs[4]), round(self.coeffs[5]))

# To actually create paths, we can't have splines parametrized on [0,1] with t,
# instead splines have to be parametrized according to arc length s.
# This problem is obvious when combining splines, cause t is completely
# arbitrary in different parts of a path (e.g.) for the first spline,
# t=1 might correspond to 10 inches of arc length, while t=1 at the
# another spline segment might correspond to 20 inches.
# Parametrizing is pretty simple itself (second week of multi var calc)
# 1. s(t) = int 0->t |r'(u)| du
# 2. t(s) = inverse function of s(t)
# 3. plug t(s) into r(t) (our spline)
# 4. r(t(s)) now is our give parametrized by arc length
# from here on out r((t(s)) will just be referred to as r(s) since its already reparamed
# Ryan (rbrott) talks about this in section 5 of his Quintic Splines for FTC paper
# r(s) obviously has different derivatives now, but they are fairly simple to find
# by just chain ruling stuff
# d/ds r(s) = 
# = d/ds (r(t(s)))
# = r'(t(s)) * t'(s)
# r'(t(s)) is just r'(t)
# = r'(t) * t'(s)
# t(s) is a pain to compute analytically, but since its an inverse to s(t),
# t'(s) = 1 / s'(t)
# so d/ds r(s) = r'(t) * (1 / s'(t))
# s'(t) = d/dt int 0->t |r'(u)| du
# = d/dt (|R'(t)| - |R'(0)|)
# = |r'(t)|
# so d/dt r(s) = r'(t) * (1 / |r'(t)|)
# this is cool and all, but we need a way to reparametrize from s->t still
# see https://github.com/GrappleRobotics/Pathfinder/blob/master/Pathfinder/src/include/grpl/pf/path/arc_parameterizer.h
class Spline:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.length = 0
        tParams = [(0,1)]
        # see https://github.com/GrappleRobotics/Pathfinder/blob/master/Pathfinder/src/include/grpl/pf/path/arc_parameterizer.h
        # on how to parametrize a spline


    def get(self, t, n) -> Vector:
        xt = self.x.get(t)
        yt = self.y.get(t)
        if n == 1: return Vector(xt.first, yt.first)
        elif n == 2: return Vector(xt.second, yt.second)
        elif n == 3: return Vector(xt.third, yt.third)
        else: return Vector(xt.zero, yt.zero)

    def start(self): return self.get(0.0)
    def end(self): return self.get(self.length)

if __name__ == "__main__":
    polynomial = Quintic(
            DifferentiablePoint(0, 0, 0),
            DifferentiablePoint(16, 0, 0)
            )

    print(polynomial.get(0.5).zero)

    




