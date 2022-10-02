from math import atan2, hypot, sqrt, pi
import numpy as np

def epsilonEquals(this, other, epsilon=0.0000001): return abs(this - other) < epsilon
def angleNorm(angle):
    a = angle
    tau = 2 * pi
    while a > pi: a -= tau
    while a < -pi: a += tau
    return a

class Vector:
    def __init__(self, x=0.0, y=0.0):
        self.x = float(x)
        self.y = float(y)

    @classmethod
    def polar(cls, r, theta): return Vector(float(r * np.cos(theta)), float(r * np.sin(theta)))

    def angle(self): return float(atan2(self.y, self.x))

    def dot(self, other): return other.x * self.x + other.y * self.y

    def cross(self, other) -> float:
        ret = self.x * other.y - self.y * other.x
        return ret 

    def angleBetween(self, other): return np.arccos((self.dot(other) / (self.norm() * other.norm())))

    def norm(self): return sqrt(self.x ** 2 + self.y ** 2)

    def normalized(self):
        len = self.norm()
        return Vector(self.x / len, self.y / len)

    def projectOnto(self, other): return other.times(self.dot(other)).div(other.dot(other))

    def times(self, d): return Vector(float(self.x * d), float(self.y * d))

    def div(self, d): return self.times(1.0/d)

    def rotate(self, angle):
        return Vector(
            float(self.x * np.cos(angle) - self.y * np.sin(angle)),
            float(self.x * np.sin(angle) + self.y * np.cos(angle))
        )

    def dist(self, other): return hypot(self.x-other.x, self.y-other.y)

    def plus(self, other): return Vector(self.x + other.x, self.y + other.y)
    
    def neg(self): return Vector(-self.x, -self.y)

    def minus(self, other): return self.plus(other.neg())

    def __str__(self): return "{:.5g}, {:.5g}".format(self.x, self.y)

class Pose:
    def __init__(self, x=0, y=0, heading=0):
        self.x = float(x)
        self.y = float(y)
        self.vec = Vector(self.x, self.y)
        self.heading = float(heading)

    def times(self, s): return Pose(self.vec.times(s), self.heading * s)

