from geometry import *
from collections import deque
from typing import List

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
#
# this arc class is used to reparametrize, also pulled from above link ^
# template <typename output_iterator_t>
# size_t parameterize(spline<2> &spline, output_iterator_t &&curve_begin, const size_t max_curve_count,
#                     double t_lo = 0, double t_hi = 1) {
#   _has_overrun = false;
#   if (max_curve_count <= 0) {
#     _has_overrun = true;
#     return 0;
#   }
#
#   double t_mid = (t_hi + t_lo) / 2.0;
#   double k_lo  = spline.curvature(t_lo);
#   double k_hi  = spline.curvature(t_hi);
#
#   augmented_arc2d arc{spline.position(t_lo), spline.position(t_mid), spline.position(t_hi), k_lo, k_hi};
#
#   bool subdivide = (fabs(k_hi - k_lo) > _max_delta_curvature) || (arc.length() > _max_arc_length);
#
#   if (subdivide) {
#     output_iterator_t head = curve_begin;
#
#     size_t len = parameterize(spline, head, max_curve_count, t_lo, t_mid);
#     len += parameterize(spline, head, max_curve_count - len, t_mid, t_hi);
#     return len;
#   } else {
#     arc.set_curvature(k_lo, k_hi);
#     *(curve_begin++) = arc;
#     return 1;
#   }
# }

class Arc:
    def __init__(self, start: Vector, mid: Vector, end: Vector):
        coeff_matrix = [
                [2 * (start.x - end.x), 2 * (start.y - end.y)],
                [2 * [start.x - mid.x, 2 * (start.y - mid.y)]]]
        coeff_det = coeff_matrix[0][0] * coeff_matrix[1][1] - coeff_matrix[0][1] * coeff_matrix[1][0]

        if coeff_det == 0:
            self.curvature = 0
            self.ref = start
            self.delta = end.minus(start)
            self.length = self.delta.norm()
            self.angle_offset = atan2(self.delta.y, self.delta.x)
        else:
            sNN, mNN, eNN = start.norm * start.norm, mid.norm * mid.norm, end.norm * end.norm
            rVec = [sNN - eNN, sNN - mNN]
            inverse1 = Vector(coeff_matrix[1][1] / coeff_det, -coeff_matrix[0][1] / coeff_det)
            inverse2 = Vector(-coeff_matrix[1][0] / coeff_det, coeff_matrix[0][0] / coeff_det)
            self.ref = Vector(inverse1.dot(rVec), inverse2.dot(rVec))
            self.angle_offset = atan2(start.minus(self.ref).y, start.minus(self.ref).x)
            angle1 = atan2(end.minus(self.ref).y, end.minus(self.ref).x)
            self.curvature = 1.0 / (start - self.ref).norm()
            self.length = abs(angle1 - self.angle_offset) / self.curvature
            if angle1 < self.angle_offset: self.curvature *= -1
        self.curvature_set = False

    def set_curvature(self, start_k, end_k):
        self.curvature = start_k
        self.dk_ds = (end_k - start_k) / self.length
        self.curvature_set = True

    def get_correct_curvature(self, s):
        if self.curvature_set:
            return self.curvature + s * self.dk_ds

    # s is distance along curve
    def linearlyInterpolate(self, s) -> Vector:
        return self.ref + self.delta * (s / self.length)

    def get(self, s) -> Vector:
        if self.curvature != 0:
            return self.ref + Vector.polar(1.0 / self.curvature, self.angle_offset + (s * self.curvature))
        else: 
            return self.linearlyInterpolate(s)

    def set_t(self, tStart, tEnd):
        self.tStart = tStart
        self.tEnd = tEnd
        self.dt = tEnd - tStart

    def interpolateSAlongT(self, s):
        return self.tStart + self.dt * (s / self.length())

class Spline:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.length = 0
        tParams = deque[(0,1)]
        self.arcs: List[Arc] = []
        its = 0
        while not len(tParams) == 0:
            curr = tParams.popleft()
            midT = (curr[0] + curr[1]) / 2

            startV = self.get(curr[0])
            endV = self.get(curr[1])
            midV = self.get(midT)

            klo = self.getK(curr[0])
            khi = self.getK(curr[1])
            dk = abs(khi = klo)
            arc = Arc(startV, midV, endV)
            subdivide = dk > 0.01 or arc.length > 1.0 
            if subdivide:
                tParams.append((midT, curr[1]))
                tParams.append((curr[0], midT))
            else:
                arc.set_curvature(klo, khi)
                arc.set_t(curr[0], curr[1])
                self.arcs.append(arc)
                self.length += arc.length
                its+=1

            if its > 1000:
                raise Exception("we fucked up")

    def get(self, t, n) -> Vector:
        xt = self.x.get(t)
        yt = self.y.get(t)
        if n == 1: return Vector(xt.first, yt.first)
        elif n == 2: return Vector(xt.second, yt.second)
        elif n == 3: return Vector(xt.third, yt.third)
        else: return Vector(xt.zero, yt.zero)

    # from multi: k = (a x v) / |v|^3
    def getK(self, t) -> float: return self.get(t, 2).cross(self.get(t, 1)) / pow(self.get(t, 1).norm(), 3)

    # now that we have our spline parametrized into arcs,
    # we can find the corresponding t with s by iterating across
    # our arc array and finding what iteration it is at
    # ok this is a shitty name but like... is it really?
    def t(self, s): 
        if s < 0: return 0
        if s > self.length: return 1
        arcLengthSum = 0
        intdt = 0
        its = 0
        while arcLengthSum < s:
            workingarc = self.arcs[its]
            if arcLengthSum + workingarc.length > s:
                return intdt + workingarc.interpolateSAlongT(s)
            arcLengthSum += workingarc.length
            intdt += workingarc.dt
            its += 1
        raise Exception("i think ur pretty fucking bad a coding neil")

    # s'(t) = d/dt int 0->t |r'(u)| du
    # = |r'(t)|
    # expand |r'(t)|
    # s'(t) = |r'(t)|
    # = norm(x'(t), y'(t))
    # = sqrt(x'(t)^2 + y'(t)^2)
    # take another deriv
    # s''(t) = (1 / sqrt(x'(t)^2 + y'(t)^2)) * (2 * x'(t) * x''(t) + 2 * y'(t) * y''(t))
    # ok lets try to simplify that xdddd
    # 2 * x'(t) * x''(t) + 2 * y'(t) * y''(t) can factor 2 out
    # = 2 * (x'(t) * x''(t) + y'(t) * y''(t)) last part of this is just tDeriv dot tDeriv2
    # = 2 * tDeriv dot tDeriv2
    # and that denom is just s'(t)
    # s''(t) = (2 * tDeriv dot tDeriv2) / sDeriv(t)
    # i dont want to take another fucking derivative 
    def sDeriv(self, t, n):
        if n == 1: return self.get(t, 1).norm()
        if n == 2: (2 * self.get(t, 1).dot(self.get(t, 2))) / self.sDeriv(t)
        raise Exception("im lazy and didn't implement more derivatives")

    def deriv(self, s, n):
        t = self.t(s)
        if n == 1: return self.get(t, 1).norm()
        if n == 2: return self.get(t, 2) * pow(self.sDeriv(t, 1), 2) + self.get(t, 1) * self.sDeriv(t, 2)
        raise Exception("im lazy and didn't implement more derivatives")

    def start(self): return self.get(0.0)
    def end(self): return self.get(self.length)

    
class SplineWithAngle(Spline):
    def __init__(self, x, y):
        super().__init__(x, y)

    def angle(self, s, n=0):
        if n == 0: return self.deriv(s, 1).angle()
        if n == 1: return self.deriv(s).cross(self.deriv(s, 2))
        raise Exception("didn't implement more angle derivs")

if __name__ == "__main__":
    polynomial = Quintic(
            DifferentiablePoint(0, 0, 0),
            DifferentiablePoint(16, 0, 0)
            )

    print(polynomial.get(0.5).zero)

    




