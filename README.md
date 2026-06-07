# elliptic-curve-ops

A pure-Rust library for **elliptic curve arithmetic** over finite fields, providing
point operations on curves in short Weierstrass form **y² = x³ + ax + b** over
prime fields GF(p).

## Why This Matters

Elliptic curves are the backbone of modern public-key cryptography. Algorithms like
ECDSA (digital signatures), ECDH (key exchange), and EdDSA all rely on the algebraic
structure of elliptic curve groups. Understanding the underlying arithmetic — point
addition, scalar multiplication, and the group order — is essential for anyone working
in cryptography, number theory, or blockchain systems.

This library provides a clean, well-documented, and fully tested implementation of
these core operations, making it suitable for both educational exploration and as a
foundation for cryptographic protocol implementations.

## Features

- **Point addition** using the geometric chord-tangent construction
- **Point doubling** via the tangent line formula
- **Scalar multiplication** with the efficient double-and-add algorithm
- **Curve order computation** (number of rational points) with Hasse bound verification
- **Point order** computation for individual points
- **Point compression/decompression** (x-coordinate + sign bit)
- **Modular arithmetic helpers**: mod_pow, mod_inverse, mod_sqrt (Tonelli-Shanks)
- **Point enumeration** over small fields

## Mathematical Background

### Short Weierstrass Form

An elliptic curve over a prime field GF(p) (with p > 3) is defined by:

```
y² = x³ + ax + b    where 4a³ + 27b² ≢ 0 (mod p)
```

The discriminant condition ensures the curve is non-singular (no cusps or self-intersections).

### Group Law

Points on the curve form an **abelian group** under the following operations:

1. **Identity**: The point at infinity O is the additive identity: P + O = P
2. **Inverse**: For P = (x, y), the inverse is -P = (x, -y mod p)
3. **Addition (chord rule)**: Given two distinct points P₁, P₂, draw a line through them.
   The third intersection with the curve, reflected over the x-axis, gives P₁ + P₂.

The slope λ is computed as:

```
λ = (y₂ - y₁) · (x₂ - x₁)⁻¹  mod p
x₃ = λ² - x₁ - x₂            mod p
y₃ = λ(x₁ - x₃) - y₁         mod p
```

4. **Doubling (tangent rule)**: When P₁ = P₂, use the tangent line:

```
λ = (3x₁² + a) · (2y₁)⁻¹  mod p
```

### Scalar Multiplication

Scalar multiplication kP = P + P + ... + P (k times) is computed using the
**double-and-add** algorithm, requiring only O(log k) point additions.

### Hasse's Theorem

For a curve over GF(p), the number of rational points N satisfies:

```
|N - (p + 1)| ≤ 2√p
```

### Point Compression

Since y only appears as y² in the curve equation, a point can be compressed to
its x-coordinate plus a single bit indicating whether y is odd or even. Decompression
requires computing a modular square root, which is done via the Tonelli-Shanks algorithm.

## Usage

Add to your `Cargo.toml`:

```toml
[dependencies]
elliptic-curve-ops = "0.1.0"
```

### Basic Point Operations

```rust
use elliptic_curve_ops::{EllipticCurve, Point};

// Create curve y² = x³ + 2x + 3 over GF(97)
let curve = EllipticCurve::new(2, 3, 97);

// Define points
let p = Point::new(3, 6);
let q = Point::new(10, 28);

// Verify they're on the curve
assert!(curve.is_on_curve(&p));
assert!(curve.is_on_curve(&q));

// Point addition
let sum = curve.add(&p, &q);
println!("P + Q = {}", sum);

// Point doubling
let doubled = curve.add(&p, &p);
println!("2P = {}", doubled);
```

### Scalar Multiplication

```rust
use elliptic_curve_ops::{EllipticCurve, Point};

let curve = EllipticCurve::new(2, 3, 97);
let p = Point::new(3, 6);

// Compute 5P using double-and-add
let result = curve.scalar_mul(&p, 5);
println!("5P = {}", result);

// Verify: 5P = P + P + P + P + P
let mut check = p;
for _ in 0..4 {
    check = curve.add(&check, &p);
}
assert_eq!(result, check);
```

### Curve Order and Point Order

```rust
use elliptic_curve_ops::{EllipticCurve, Point};

let curve = EllipticCurve::new(2, 3, 97);
let p = Point::new(3, 6);

// Number of rational points on the curve
let order = curve.curve_order();
println!("Curve order: {}", order);

// Order of specific point
let pt_order = curve.point_order(&p);
println!("Point order: {}", pt_order);

// Verify: pt_order * P should be infinity
let result = curve.scalar_mul(&p, pt_order);
assert!(result.is_infinity);
```

### Point Compression

```rust
use elliptic_curve_ops::{EllipticCurve, Point};

let curve = EllipticCurve::new(2, 3, 97);
let p = Point::new(3, 6);

// Compress to x-coordinate + sign bit
let (x, is_odd) = curve.compress(&p).unwrap();
println!("Compressed: (x={}, y_is_odd={})", x, is_odd);

// Decompress back
let recovered = curve.decompress(x, is_odd).unwrap();
assert_eq!(recovered, p);
```

## API Reference

| Type/Function | Description |
|---|---|
| `Point` | Affine point (x, y) or infinity |
| `Point::infinity()` | The point at infinity |
| `Point::negate(p)` | Negate point modulo p |
| `EllipticCurve` | Curve y² = x³ + ax + b over GF(p) |
| `curve.add(p1, p2)` | Point addition (chord-tangent) |
| `curve.scalar_mul(p, k)` | Scalar multiplication (double-and-add) |
| `curve.curve_order()` | Number of rational points |
| `curve.point_order(p)` | Order of a specific point |
| `curve.compress(p)` | Compress point to (x, sign_bit) |
| `curve.decompress(x, sign)` | Recover point from compressed form |
| `curve.all_points()` | Enumerate all curve points |
| `mod_pow(base, exp, m)` | Fast modular exponentiation |
| `mod_inverse(a, m)` | Modular multiplicative inverse |
| `mod_sqrt(n, p)` | Modular square root (Tonelli-Shanks) |

## Testing

```bash
cargo test
```

All operations are verified against known results, including Hasse bound checks,
round-trip compression/decompression, and group law properties.

## License

MIT License. See [LICENSE](LICENSE) for details.
