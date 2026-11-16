//! The few 3-vector operations mesh frames need.

pub type Vec3 = [f64; 3];

/// MCNP's `VEC` for a cylindrical FMESH whose card gives none (MCNP 6.2 and
/// 6.3: "no VEC keyword on fmesh card; default 1 0 0 is used"; with an `AXS`
/// along x that is a fatal error, as the two are parallel).
pub const DEFAULT_VEC: Vec3 = [1.0, 0.0, 0.0];

pub fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub fn add(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

pub fn sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

pub fn scale(a: Vec3, s: f64) -> Vec3 {
    [a[0] * s, a[1] * s, a[2] * s]
}

pub fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// `a` scaled to length 1; None when it is zero or not finite.
pub fn normalize(a: Vec3) -> Option<Vec3> {
    let n = dot(a, a).sqrt();
    (n > 0.0 && n.is_finite()).then(|| scale(a, 1.0 / n))
}

/// `v` without its part along the unit vector `axis`, normalized: a
/// cylinder's `VEC` as MCNP uses it. None when `v` is parallel to `axis`.
pub fn perpendicular_unit(v: Vec3, axis: Vec3) -> Option<Vec3> {
    let p = sub(v, scale(axis, dot(v, axis)));
    // Parallel up to rounding: what is left is noise.
    if dot(p, p).sqrt() <= 1e-12 * dot(v, v).sqrt() {
        return None;
    }
    normalize(p)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn perpendicular_part() {
        assert_eq!(
            perpendicular_unit(DEFAULT_VEC, [0.0, 0.0, 1.0]),
            Some([1.0, 0.0, 0.0])
        );
        assert_eq!(perpendicular_unit(DEFAULT_VEC, [1.0, 0.0, 0.0]), None);
        let a = normalize([1.0, 1.0, 0.0]).unwrap();
        let v = perpendicular_unit(DEFAULT_VEC, a).unwrap();
        assert!(dot(v, a).abs() < 1e-15);
        assert!((v[0] - 0.5f64.sqrt()).abs() < 1e-15 && (v[1] + 0.5f64.sqrt()).abs() < 1e-15);
        assert_eq!(normalize([0.0; 3]), None);
    }
}
