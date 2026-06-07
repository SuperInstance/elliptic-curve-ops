//! # Elliptic Curve Operations
//!
//! A pure-Rust library for elliptic curve arithmetic over finite fields.
//! Supports the short Weierstrass form y² = x³ + ax + b over prime fields,
//! with point addition, scalar multiplication, order computation, and
//! point compression/decompression.
//!
//! ## Quick Start
//!
//! ```
//! use elliptic_curve_ops::{EllipticCurve, Point};
//!
//! // y² = x³ + 2x + 3 over GF(97)
//! let curve = EllipticCurve::new(2, 3, 97);
//!
//! let p = Point::new(3, 6);
//! let q = Point::new(1, 43);
//!
//! let r = curve.add(&p, &q);
//! assert!(curve.is_on_curve(&r));
//! ```

use std::fmt;

/// A point on an elliptic curve, or the point at infinity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Point {
    pub x: i64,
    pub y: i64,
    pub is_infinity: bool,
}

impl Point {
    /// Create a new affine point.
    pub fn new(x: i64, y: i64) -> Self {
        Point { x, y, is_infinity: false }
    }

    /// The point at infinity (additive identity).
    pub fn infinity() -> Self {
        Point { x: 0, y: 0, is_infinity: true }
    }

    /// Negate a point: -(x, y) = (x, -y).
    pub fn negate(&self, p: i64) -> Self {
        if self.is_infinity {
            return Self::infinity();
        }
        Point::new(self.x, (-self.y).rem_euclid(p))
    }
}

impl fmt::Display for Point {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_infinity {
            write!(f, "O (∞)")
        } else {
            write!(f, "({}, {})", self.x, self.y)
        }
    }
}

/// An elliptic curve in short Weierstrass form: y² = x³ + ax + b over GF(p).
#[derive(Debug, Clone, Copy)]
pub struct EllipticCurve {
    pub a: i64,
    pub b: i64,
    pub p: i64,
}

impl EllipticCurve {
    /// Create a new curve y² = x³ + ax + b over GF(p).
    ///
    /// # Panics
    /// Panics if the discriminant 4a³ + 27b² ≡ 0 (mod p) (singular curve).
    pub fn new(a: i64, b: i64, p: i64) -> Self {
        let disc = (4 * a.pow(3) + 27 * b.pow(2)).rem_euclid(p);
        assert!(disc != 0, "Curve is singular: discriminant is zero mod {}", p);
        EllipticCurve { a, b, p }
    }

    /// Check if a point lies on the curve.
    pub fn is_on_curve(&self, point: &Point) -> bool {
        if point.is_infinity {
            return true;
        }
        let lhs = (point.y.pow(2)).rem_euclid(self.p);
        let rhs = (point.x.pow(3) + self.a * point.x + self.b).rem_euclid(self.p);
        lhs == rhs
    }

    /// Point addition using the chord-tangent construction.
    ///
    /// If P = Q, uses the tangent line (doubling).
    /// If P = -Q, returns the point at infinity.
    /// If one point is infinity, returns the other.
    pub fn add(&self, p1: &Point, p2: &Point) -> Point {
        if p1.is_infinity {
            return *p2;
        }
        if p2.is_infinity {
            return *p1;
        }
        if p1.x == p2.x && p1.y != p2.y {
            return Point::infinity();
        }
        if *p1 == *p2 && p1.y == 0 {
            return Point::infinity();
        }

        let lambda = if p1 == p2 {
            // Tangent: λ = (3x₁² + a) / (2y₁)
            let num = (3 * p1.x.pow(2) + self.a).rem_euclid(self.p);
            let den = (2 * p1.y).rem_euclid(self.p);
            (num * mod_inverse(den, self.p)).rem_euclid(self.p)
        } else {
            // Chord: λ = (y₂ - y₁) / (x₂ - x₁)
            let num = (p2.y - p1.y).rem_euclid(self.p);
            let den = (p2.x - p1.x).rem_euclid(self.p);
            (num * mod_inverse(den, self.p)).rem_euclid(self.p)
        };

        let x3 = (lambda.pow(2) - p1.x - p2.x).rem_euclid(self.p);
        let y3 = (lambda * (p1.x - x3) - p1.y).rem_euclid(self.p);

        Point::new(x3, y3)
    }

    /// Scalar multiplication using the double-and-add algorithm.
    pub fn scalar_mul(&self, point: &Point, k: i64) -> Point {
        if k == 0 || point.is_infinity {
            return Point::infinity();
        }
        let k = k.rem_euclid(self.p); // reduce k
        if k == 0 {
            return Point::infinity();
        }

        let mut result = Point::infinity();
        let mut addend = *point;
        let mut k = k;

        while k > 0 {
            if k & 1 == 1 {
                result = self.add(&result, &addend);
            }
            addend = self.add(&addend, &addend);
            k >>= 1;
        }

        result
    }

    /// Compute the order of a point (smallest n > 0 such that nP = O).
    pub fn point_order(&self, point: &Point) -> i64 {
        if point.is_infinity {
            return 1;
        }
        let curve_order = self.curve_order();
        // Order divides curve order; find smallest divisor
        let mut n = curve_order;
        // Try dividing by small primes
        let mut d = 2;
        while d * d <= n {
            while n % d == 0 && self.scalar_mul(point, n / d).is_infinity {
                n /= d;
            }
            d += 1;
        }
        n
    }

    /// Compute the order of the curve (number of rational points) via enumeration.
    ///
    /// Note: O(√p) Schoof's algorithm is not implemented; this uses direct enumeration.
    pub fn curve_order(&self) -> i64 {
        let mut count = 1i64; // point at infinity
        for x in 0..self.p {
            let rhs = (x.pow(3) + self.a * x + self.b).rem_euclid(self.p);
            if rhs == 0 {
                count += 1; // (x, 0)
            } else {
                // Check if rhs is a quadratic residue
                let legendre = mod_pow(rhs, (self.p - 1) / 2, self.p);
                if legendre == 1 {
                    count += 2; // (x, y) and (x, -y)
                }
            }
        }
        count
    }

    /// Compress a point to a single x-coordinate and a sign bit.
    ///
    /// Returns `(x, is_y_odd)`.
    pub fn compress(&self, point: &Point) -> Option<(i64, bool)> {
        if point.is_infinity {
            return None;
        }
        let is_odd = point.y % 2 != 0;
        Some((point.x, is_odd))
    }

    /// Decompress a point from an x-coordinate and sign bit.
    pub fn decompress(&self, x: i64, is_y_odd: bool) -> Option<Point> {
        let rhs = (x.pow(3) + self.a * x + self.b).rem_euclid(self.p);
        let y = mod_sqrt(rhs, self.p)?;
        let y_odd = y % 2 != 0;
        let y_final = if y_odd == is_y_odd { y } else { (self.p - y).rem_euclid(self.p) };
        let point = Point::new(x, y_final);
        debug_assert!(self.is_on_curve(&point));
        Some(point)
    }

    /// Find all points on the curve (enumeration).
    pub fn all_points(&self) -> Vec<Point> {
        let mut points = vec![Point::infinity()];
        for x in 0..self.p {
            let rhs = (x.pow(3) + self.a * x + self.b).rem_euclid(self.p);
            if rhs == 0 {
                points.push(Point::new(x, 0));
            } else {
                let legendre = mod_pow(rhs, (self.p - 1) / 2, self.p);
                if legendre == 1 {
                    let y = mod_sqrt(rhs, self.p).unwrap();
                    points.push(Point::new(x, y));
                    points.push(Point::new(x, (self.p - y).rem_euclid(self.p)));
                }
            }
        }
        points
    }
}

/// Modular exponentiation via repeated squaring.
pub fn mod_pow(base: i64, exp: i64, modulus: i64) -> i64 {
    if modulus == 1 { return 0; }
    let mut result = 1i64;
    let mut base = base.rem_euclid(modulus);
    let mut exp = exp;
    while exp > 0 {
        if exp & 1 == 1 {
            result = (result * base) % modulus;
        }
        exp >>= 1;
        base = (base * base) % modulus;
    }
    result
}

/// Modular inverse via extended Euclidean algorithm.
pub fn mod_inverse(a: i64, m: i64) -> i64 {
    let (_, x, _) = extended_gcd(a.rem_euclid(m), m);
    x.rem_euclid(m)
}

/// Extended GCD: returns (gcd, x, y) such that a*x + b*y = gcd.
pub fn extended_gcd(a: i64, b: i64) -> (i64, i64, i64) {
    if a == 0 {
        return (b, 0, 1);
    }
    let (g, x, y) = extended_gcd(b % a, a);
    (g, y - (b / a) * x, x)
}

/// Modular square root via Tonelli-Shanks algorithm.
pub fn mod_sqrt(n: i64, p: i64) -> Option<i64> {
    if p == 2 { return Some(n % 2); }
    let n = n.rem_euclid(p);
    if n == 0 { return Some(0); }

    // Check QR
    if mod_pow(n, (p - 1) / 2, p) != 1 {
        return None;
    }

    // Factor out powers of 2 from p-1: p-1 = Q * 2^S
    let mut s = 0i64;
    let mut q = p - 1;
    while q % 2 == 0 {
        q /= 2;
        s += 1;
    }

    if s == 1 {
        let r = mod_pow(n, (p + 1) / 4, p);
        return Some(r);
    }

    // Find a non-residue z
    let mut z = 2i64;
    while mod_pow(z, (p - 1) / 2, p) != p - 1 {
        z += 1;
    }

    let mut m = s;
    let mut c = mod_pow(z, q, p);
    let mut t = mod_pow(n, q, p);
    let mut r = mod_pow(n, (q + 1) / 2, p);

    loop {
        if t == 1 { return Some(r); }
        let mut i = 1;
        let mut t2i = (t * t) % p;
        while t2i != 1 {
            t2i = (t2i * t2i) % p;
            i += 1;
        }
        let b = mod_pow(c, mod_pow(2, m - i - 1, p - 1), p);
        m = i;
        c = (b * b) % p;
        t = (t * c) % p;
        r = (r * b) % p;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // y² = x³ + 2x + 3 over GF(97)
    fn test_curve() -> EllipticCurve {
        EllipticCurve::new(2, 3, 97)
    }

    #[test]
    fn test_curve_creation() {
        let curve = test_curve();
        assert_eq!(curve.a, 2);
        assert_eq!(curve.b, 3);
        assert_eq!(curve.p, 97);
    }

    #[test]
    #[should_panic]
    fn test_singular_curve_panics() {
        // 4*0³ + 27*0² = 0 mod p → singular
        EllipticCurve::new(0, 0, 97);
    }

    #[test]
    fn test_point_on_curve() {
        let curve = test_curve();
        // (3, 6): 36 mod 97 = 36; 27+6+3=36 mod 97 ✓
        let p = Point::new(3, 6);
        assert!(curve.is_on_curve(&p));
    }

    #[test]
    fn test_point_not_on_curve() {
        let curve = test_curve();
        let p = Point::new(1, 1);
        assert!(!curve.is_on_curve(&p));
    }

    #[test]
    fn test_infinity_on_curve() {
        let curve = test_curve();
        assert!(curve.is_on_curve(&Point::infinity()));
    }

    #[test]
    fn test_point_addition_distinct() {
        let curve = test_curve();
        let p = Point::new(3, 6);
        let q = Point::new(1, 43);
        let r = curve.add(&p, &q);
        assert!(curve.is_on_curve(&r));
        assert!(!r.is_infinity);
    }

    #[test]
    fn test_point_doubling() {
        let curve = test_curve();
        let p = Point::new(3, 6);
        let doubled = curve.add(&p, &p);
        assert!(curve.is_on_curve(&doubled));
    }

    #[test]
    fn test_add_inverse_gives_infinity() {
        let curve = test_curve();
        let p = Point::new(3, 6);
        let neg_p = p.negate(97);
        let r = curve.add(&p, &neg_p);
        assert!(r.is_infinity);
    }

    #[test]
    fn test_add_identity() {
        let curve = test_curve();
        let p = Point::new(3, 6);
        let r = curve.add(&p, &Point::infinity());
        assert_eq!(r, p);
    }

    #[test]
    fn test_scalar_multiplication() {
        let curve = test_curve();
        let p = Point::new(3, 6);
        let double = curve.add(&p, &p);
        let scalar_double = curve.scalar_mul(&p, 2);
        assert_eq!(double, scalar_double);
    }

    #[test]
    fn test_scalar_mul_zero() {
        let curve = test_curve();
        let p = Point::new(3, 6);
        let r = curve.scalar_mul(&p, 0);
        assert!(r.is_infinity);
    }

    #[test]
    fn test_scalar_mul_identity() {
        let curve = test_curve();
        let p = Point::new(3, 6);
        let r = curve.scalar_mul(&p, 1);
        assert_eq!(r, p);
    }

    #[test]
    fn test_curve_order() {
        let curve = test_curve();
        let order = curve.curve_order();
        assert!(order > 0);
        assert!(order <= 2 * curve.p + 1); // Hasse bound
    }

    #[test]
    fn test_hasse_bound() {
        let curve = test_curve();
        let n = curve.curve_order();
        let p = curve.p as f64;
        let lower = (p + 1.0 - 2.0 * p.sqrt()).ceil() as i64;
        let upper = (p + 1.0 + 2.0 * p.sqrt()).floor() as i64;
        assert!(n >= lower && n <= upper, "Order {} violates Hasse bound [{}, {}]", n, lower, upper);
    }

    #[test]
    fn test_point_order() {
        let curve = test_curve();
        let p = Point::new(3, 6);
        let order = curve.point_order(&p);
        assert!(order > 0);
        let result = curve.scalar_mul(&p, order);
        assert!(result.is_infinity);
    }

    #[test]
    fn test_compression_decompression() {
        let curve = test_curve();
        let p = Point::new(3, 6);
        let (x, is_odd) = curve.compress(&p).unwrap();
        let decompressed = curve.decompress(x, is_odd).unwrap();
        assert_eq!(decompressed, p);
    }

    #[test]
    fn test_compress_infinity_returns_none() {
        let curve = test_curve();
        assert!(curve.compress(&Point::infinity()).is_none());
    }

    #[test]
    fn test_all_points_count_matches_order() {
        let curve = test_curve();
        let points = curve.all_points();
        assert_eq!(points.len() as i64, curve.curve_order());
    }

    #[test]
    fn test_negate() {
        let curve = test_curve();
        let p = Point::new(3, 6);
        let neg = p.negate(curve.p);
        assert!(curve.is_on_curve(&neg));
        let sum = curve.add(&p, &neg);
        assert!(sum.is_infinity);
    }

    #[test]
    fn test_mod_inverse() {
        assert_eq!(mod_inverse(3, 7), 5);
        assert_eq!((3 * mod_inverse(3, 7)) % 7, 1);
    }

    #[test]
    fn test_mod_pow() {
        assert_eq!(mod_pow(2, 10, 1000), 1024 % 1000);
        assert_eq!(mod_pow(3, 4, 5), 81 % 5);
    }
}
