from geometry import *
from collections import deque
from typing import List
from math import radians
import scipy.linalg as la

class DifferentiablePoint():
    def __init__(self, zero=0, first=0,  second=0):
        self.first = first
        self.zero = zero
        self.second = second

class DifferentiablePoint2d:
    def __init__(self, zero: Vector, first: Vector):
        self.x = DifferentiablePoint(zero.x, first.x)
        self.y = DifferentiablePoint(zero.y, first.y)

class Quintic:
    def __init__(self, start: DifferentiablePoint, end: DifferentiablePoint) -> None:
        self.end = end
        self.start = start
        COEFF_MATRIX = np.array([
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
            [0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 2.0, 0.0, 0.0],
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            [5.0, 4.0, 3.0, 2.0, 1.0, 0.0],
            [20.0, 12.0, 6.0, 2.0, 0.0, 0.0]
        ])

        target = np.array([
            [start.zero],
            [start.first],
            [0.0],
            [end.zero],
            [end.first],
            [0.0]
        ])

        solved = la.solve(COEFF_MATRIX, target)
        self.coeffs = [float(x) for x in solved]

    def get(self, t: float) -> DifferentiablePoint:
        coeffs = self.coeffs
        return DifferentiablePoint(
            coeffs[0] * pow(t, 5) + coeffs[1] * pow(t, 4) + coeffs[2] * pow(t, 3) + coeffs[3] * pow(t, 2) + coeffs[4] * t + coeffs[5],
            5 * coeffs[0] * pow(t, 4) + 4 * coeffs[1] * pow(t, 3) + 3 * coeffs[2] * pow(t, 2) + 2 * coeffs[3] * t + coeffs[4],
            20 * coeffs[0] * pow(t, 3) + 12 * coeffs[1] * pow(t, 2) + 6 * coeffs[2] * t + 2 * coeffs[3])

    def __str__(self): 
        return "{}t^5 + {}t^4 + {}t^3 + {}t^2 + {}t + {}".format(round(self.coeffs[0]), round(self.coeffs[1]), 
            round(self.coeffs[2]), round(self.coeffs[3]), round(self.coeffs[4]), round(self.coeffs[5]))

class Arc:
    def __init__(self, start: Vector, mid: Vector, end: Vector):
        coeff_matrix = [
                [2 * (start.x - end.x), 2 * (start.y - end.y)],
                [2 * (start.x - mid.x), 2 * (start.y - mid.y)]]
        coeff_det = coeff_matrix[0][0] * coeff_matrix[1][1] - coeff_matrix[0][1] * coeff_matrix[1][0]

        if coeff_det == 0:
            self.curvature = 0
            self.ref = start
            self.delta = end.minus(start)
            self.length = self.delta.norm()
            self.angle_offset = atan2(self.delta.y, self.delta.x)
        else:
            sNN, mNN, eNN = start.norm() * start.norm(), mid.norm() * mid.norm(), end.norm() * end.norm()
            # rVec = [sNN - eNN, sNN - mNN]
            rVec = Vector(sNN -  eNN, sNN - mNN)
            inverse1 = Vector(coeff_matrix[1][1] / coeff_det, -coeff_matrix[0][1] / coeff_det)
            inverse2 = Vector(-coeff_matrix[1][0] / coeff_det, coeff_matrix[0][0] / coeff_det)
            self.ref = Vector(inverse1.dot(rVec), inverse2.dot(rVec))
            self.angle_offset = atan2(start.minus(self.ref).y, start.minus(self.ref).x)
            angle1 = atan2(end.minus(self.ref).y, end.minus(self.ref).x)
            self.curvature = 1.0 / (start.minus(self.ref)).norm()
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
            return self.ref.plus(Vector.fromPolar(1.0 / self.curvature, self.angle_offset + (s * self.curvature)))
        else: 
            return self.linearlyInterpolate(s)

    def set_t(self, tStart, tEnd):
        self.tStart = tStart
        self.tEnd = tEnd
        self.dt = tEnd - tStart

    # start + dt *  (s + 2 * length)
    #               -------------
    #                   length
    def interpolateSAlongT(self, s):
        return self.tStart + self.dt * (s / self.length)

class Spline:
    def __init__(self, x: Quintic, y: Quintic) -> None:
        self.x = x
        self.y = y
        tPairs = deque()
        tPairs.append((0.0, 1.0))
        its = 0
        self.arcs = []
        self.length = 0.0
        while len(tPairs) != 0:
            curr = tPairs.popleft()
            midT = (curr[0] + curr[1]) / 2.0

            startV = self.rt(curr[0])
            endV = self.rt(curr[1])
            midV = self.rt(midT)

            startK = self.rt(curr[0], 2).cross(self.rt(curr[1], 1)) / pow(self.rt(curr[0], 1).norm(), 3)
            endK = self.rt(curr[1], 2).cross(self.rt(curr[1], 1)) / pow(self.rt(curr[1], 1).norm(), 3)
            arc = Arc(startV, midV, endV)
            
            subdivide = abs(endK - startK) > 0.01 or arc.length > 0.25
            if subdivide:
                tPairs.append((midT, curr[1]))
                tPairs.append((curr[0], midT))
            else:
                arc.set_curvature(startK, endK)
                arc.set_t(curr[0], curr[1])
                self.arcs.append(arc)
                self.length += arc.length
                its += 1

        self.arcs.sort(key=lambda arc: arc.tStart)

    def invArc(self, s):
        if s <= 0.0: return 0.0
        if s >= self.length: return 1.0
        acc = 0.0
        for arc in self.arcs:
            if acc + arc.length > s:
                return arc.interpolateSAlongT(acc - s)
            acc += arc.length

    def rt(self, t, n = 0) -> Vector:
        xt = self.x.get(t)
        yt = self.y.get(t)
        if n == 1: return Vector(xt.first, yt.first)
        elif n == 2: return Vector(xt.second, yt.second)
        else: return Vector(xt.zero, yt.zero)

    def dsdt(self, t, n = 1):
        if n == 1: return self.rt(t, 1).norm()
        elif n == 2: return (2 * self.rt(t, 1).dot(self.rt(t, 2))) / self.dsdt(t)

    def dtds(self, t, n = 1):
        if n == 1: return 1.0 / self.dsdt(t)
        elif n == 2: return -self.dsdt(t, 2) / pow(self.dsdt(t), 3)

    def rs(self, s, n = 0) -> Vector:
        t = self.invArc(s)
        if n == 0: return self.rt(t)
        elif n == 1: return self.rt(t, 1).normalized()
        elif n == 2: return self.rt(t, 2) * pow(self.dtds(t), 2) + self.rt(t, 1) * self.dtds(t, 2)

    def get(self, s, n = 0) -> Pose:
        if n == 0: 
            v = self.rs(s)
            return Pose(v.x, v.y, self.rs(s, 1).angle())
        elif n == 1: 
            v = self.rs(s, 1)
            return Pose(v.x, v.y, self.rs(s, 1).cross(self.rs(s, 2)))
        elif n == 2: 
            v = self.rs(s, 2)
            return Pose(v.x, v.y, 0.0)

class Path:
    def __init__(self, *poses: Pose) -> None:
        curr = poses[0]
        self.curve_segments = []
        self.length = 0.0
        for target in poses[1:]:
            cv = curr.vec
            tv = target.vec
            r = cv.dist(tv)
            s = DifferentiablePoint2d(cv, Vector.fromPolar(r, curr.heading))
            e = DifferentiablePoint2d(tv, Vector.fromPolar(r, target.heading))
            spline = Spline(Quintic(s.x, e.x), Quintic(s.y, e.y))
            self.curve_segments.append(spline)
            self.length += spline.length
            curr = target

    def get(self, s, n = 0) -> Pose:
        if s <= 0.0: return self.curve_segments[0].get(0.0, n)
        if s >= self.length: return self.curve_segments[-1].get(self.curve_segments[-1].size, n)
        acc = 0.0
        for spline in self.curve_segments:
            if acc + spline.length > s: return spline.get(s - acc, n)
            acc += spline.length

    def clamp(self, x, a, b):
        if x < a: return a
        if x > b: return b
        return x

    def project(self, p: Vector, pGuess):
        acc = pGuess
        for _ in range(10):
            acc = self.clamp(acc + (p - self.get(acc).vec).dot(self.get(acc, 1).vec), 0.0, self.length)
        return acc

if __name__ == "__main__":
    start = Pose()
    path =  Path(start, Pose(20, 20, radians(90)))
    print(path.get(path.length / 2.0).vec)

