//! NTRU-style experiments over finite group algebras.
//!
//! Classical NTRU works in `Z[x] / (x^N - 1)`, which is the group algebra
//! `Z[C_N]`. This crate replaces the cyclic group by an explicit finite group
//! and runs the same left-sided NTRU-style construction over `Z[G]`.

use std::collections::BTreeMap;
use std::error::Error;
use std::fmt;

/// Errors returned by the research implementation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum NtruError {
    /// A modulus was not usable.
    InvalidModulus,
    /// A group or coefficient vector had the wrong shape.
    Shape,
    /// A modular inverse does not exist.
    NotInvertible,
    /// Parameters are inconsistent.
    InvalidParameters(String),
    /// Optional FFT backend rejected the input.
    Backend(String),
}

impl fmt::Display for NtruError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidModulus => write!(f, "invalid modulus"),
            Self::Shape => write!(f, "shape mismatch"),
            Self::NotInvertible => write!(f, "element is not invertible"),
            Self::InvalidParameters(message) => write!(f, "invalid parameters: {message}"),
            Self::Backend(message) => write!(f, "backend error: {message}"),
        }
    }
}

impl Error for NtruError {}

/// A deterministic pseudo-random generator used to keep examples reproducible
/// without adding a dependency.
#[derive(Clone, Debug)]
pub struct Lcg {
    state: u64,
}

impl Lcg {
    /// Construct a generator from a seed.
    pub fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    /// Return the next `u64`.
    pub fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.state
    }

    /// Return a value in `0..upper`.
    pub fn gen_range(&mut self, upper: usize) -> usize {
        assert!(upper > 0);
        (self.next_u64() as usize) % upper
    }
}

/// Group families currently implemented by the crate.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GroupKind {
    /// The cyclic group `C_n`.
    Cyclic,
    /// The dihedral group with rotation order `n`, of order `2n`.
    Dihedral,
    /// The symmetric group `S_n`.
    Symmetric,
}

/// An explicit finite group with a fixed coefficient ordering.
#[derive(Clone, Debug)]
pub struct FiniteGroup {
    kind: GroupKind,
    n: usize,
    name: String,
    elements: Vec<Vec<usize>>,
}

impl FiniteGroup {
    /// Construct the cyclic group `C_n`.
    pub fn cyclic(n: usize) -> Result<Self, NtruError> {
        if n == 0 {
            return Err(NtruError::InvalidParameters("n must be positive".into()));
        }
        Ok(Self {
            kind: GroupKind::Cyclic,
            n,
            name: format!("C{n}"),
            elements: (0..n).map(|i| vec![i]).collect(),
        })
    }

    /// Construct the dihedral group with rotation order `n`.
    pub fn dihedral(n: usize) -> Result<Self, NtruError> {
        if n < 3 {
            return Err(NtruError::InvalidParameters(
                "dihedral rotation order must be at least 3".into(),
            ));
        }
        Ok(Self {
            kind: GroupKind::Dihedral,
            n,
            name: format!("D{n} (order {})", 2 * n),
            elements: (0..2 * n).map(|i| vec![i]).collect(),
        })
    }

    /// Construct the symmetric group `S_n`.
    pub fn symmetric(n: usize) -> Result<Self, NtruError> {
        if n == 0 {
            return Err(NtruError::InvalidParameters("n must be positive".into()));
        }
        let elements = all_permutations(n);
        Ok(Self {
            kind: GroupKind::Symmetric,
            n,
            name: format!("S{n}"),
            elements,
        })
    }

    /// Return the group family.
    pub fn kind(&self) -> &GroupKind {
        &self.kind
    }

    /// Return the group parameter.
    pub fn n(&self) -> usize {
        self.n
    }

    /// Return the display name.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Return the group order.
    pub fn order(&self) -> usize {
        self.elements.len()
    }

    /// Return the coefficient ordering.
    pub fn elements(&self) -> &[Vec<usize>] {
        &self.elements
    }

    /// Return the product index for two basis indices.
    pub fn multiply_index(&self, lhs: usize, rhs: usize) -> usize {
        match self.kind {
            GroupKind::Cyclic => (lhs + rhs) % self.n,
            GroupKind::Dihedral => self.multiply_dihedral(lhs, rhs),
            GroupKind::Symmetric => self.multiply_symmetric(lhs, rhs),
        }
    }

    fn multiply_dihedral(&self, lhs: usize, rhs: usize) -> usize {
        let n = self.n;
        let lhs_reflection = lhs >= n;
        let rhs_reflection = rhs >= n;
        let lhs_k = if lhs_reflection { lhs - n } else { lhs };
        let rhs_k = if rhs_reflection { rhs - n } else { rhs };

        match (lhs_reflection, rhs_reflection) {
            (false, false) => (lhs_k + rhs_k) % n,
            (false, true) => n + ((rhs_k + n - lhs_k) % n),
            (true, false) => n + ((lhs_k + rhs_k) % n),
            (true, true) => (rhs_k + n - lhs_k) % n,
        }
    }

    fn multiply_symmetric(&self, lhs: usize, rhs: usize) -> usize {
        let lhs_images = &self.elements[lhs];
        let rhs_images = &self.elements[rhs];
        let product: Vec<_> = rhs_images.iter().map(|image| lhs_images[*image]).collect();
        self.elements
            .iter()
            .position(|candidate| *candidate == product)
            .expect("product permutation is in S_n")
    }
}

fn all_permutations(n: usize) -> Vec<Vec<usize>> {
    let mut current: Vec<_> = (0..n).collect();
    let mut out = vec![current.clone()];
    while next_permutation(&mut current) {
        out.push(current.clone());
    }
    out
}

fn next_permutation(values: &mut [usize]) -> bool {
    if values.len() < 2 {
        return false;
    }
    let mut pivot = values.len() - 2;
    while values[pivot] >= values[pivot + 1] {
        if pivot == 0 {
            values.reverse();
            return false;
        }
        pivot -= 1;
    }
    let mut successor = values.len() - 1;
    while values[successor] <= values[pivot] {
        successor -= 1;
    }
    values.swap(pivot, successor);
    values[pivot + 1..].reverse();
    true
}

/// Backend call counters from a run.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct BackendStats {
    counts: BTreeMap<String, usize>,
}

impl BackendStats {
    fn record(&mut self, name: &str) {
        *self.counts.entry(name.to_string()).or_insert(0) += 1;
    }

    /// Return all backend counters.
    pub fn counts(&self) -> &BTreeMap<String, usize> {
        &self.counts
    }

    /// Format counters as a compact comma-separated string.
    pub fn format(&self) -> String {
        self.counts
            .iter()
            .map(|(name, count)| format!("{name}={count}"))
            .collect::<Vec<_>>()
            .join(", ")
    }
}

/// A private key for the NTRU-style construction.
#[derive(Clone, Debug)]
pub struct PrivateKey {
    /// Small private element.
    pub f: Vec<i64>,
    /// Inverse of `f` modulo `p`.
    pub f_p: Vec<u64>,
    /// Inverse of `f` modulo `q`.
    pub f_q: Vec<u64>,
    /// Number of ternary samples tried.
    pub attempts: usize,
}

/// Result of one trial.
#[derive(Clone, Debug)]
pub struct TrialResult {
    /// Whether decryption recovered the plaintext.
    pub success: bool,
    /// Whether the centered lift matched `pgr + fm` over the integers.
    pub no_wrap: bool,
    /// Private-key search attempts.
    pub key_attempts: usize,
}

/// Summary for repeated trials.
#[derive(Clone, Debug)]
pub struct TrialSummary {
    /// Requested trial count.
    pub requested_trials: usize,
    /// Completed trial count.
    pub completed_trials: usize,
    /// Successful decryptions.
    pub successes: usize,
    /// No-wrap checks passed.
    pub no_wraps: usize,
    /// Key generation failures.
    pub keygen_failures: usize,
    /// Average key search attempts.
    pub avg_key_attempts: Option<f64>,
    /// Maximum key search attempts.
    pub max_key_attempts: Option<usize>,
    /// Backend counters.
    pub backend_stats: BackendStats,
}

/// NTRU-style experiment state over a group algebra.
#[derive(Clone, Debug)]
pub struct GroupAlgebraNtru {
    group: FiniteGroup,
    p: u64,
    q: u64,
    d: usize,
    stats: BackendStats,
}

impl GroupAlgebraNtru {
    /// Construct a scheme instance.
    pub fn new(group: FiniteGroup, p: u64, q: u64, d: usize) -> Result<Self, NtruError> {
        if !is_prime(p) {
            return Err(NtruError::InvalidParameters("p must be prime".into()));
        }
        if q <= 1 {
            return Err(NtruError::InvalidParameters(
                "q must be greater than 1".into(),
            ));
        }
        if gcd(p, q) != 1 {
            return Err(NtruError::InvalidParameters(
                "p and q must be coprime".into(),
            ));
        }
        if group.order() < 2 * d + 1 {
            return Err(NtruError::InvalidParameters(
                "group is too small for the requested ternary density".into(),
            ));
        }
        Ok(Self {
            group,
            p,
            q,
            d,
            stats: BackendStats::default(),
        })
    }

    /// Return the group.
    pub fn group(&self) -> &FiniteGroup {
        &self.group
    }

    /// Return `p`.
    pub fn p(&self) -> u64 {
        self.p
    }

    /// Return `q`.
    pub fn q(&self) -> u64 {
        self.q
    }

    /// Return `d`.
    pub fn d(&self) -> usize {
        self.d
    }

    /// Return backend counters.
    pub fn backend_stats(&self) -> &BackendStats {
        &self.stats
    }

    /// Return the multiplicative identity as residues modulo `modulus`.
    pub fn one_mod(&self, modulus: u64) -> Vec<u64> {
        let mut out = vec![0; self.group.order()];
        out[0] = 1 % modulus;
        out
    }

    /// Sample an element in `T(d1,d2)`.
    pub fn random_ternary(
        &self,
        rng: &mut Lcg,
        d1: usize,
        d2: usize,
    ) -> Result<Vec<i64>, NtruError> {
        let support = d1 + d2;
        if support > self.group.order() {
            return Err(NtruError::InvalidParameters(
                "ternary support too large".into(),
            ));
        }
        let mut positions: Vec<_> = (0..self.group.order()).collect();
        for i in 0..support {
            let j = i + rng.gen_range(self.group.order() - i);
            positions.swap(i, j);
        }
        let mut out = vec![0; self.group.order()];
        for &idx in &positions[..d1] {
            out[idx] = 1;
        }
        for &idx in &positions[d1..support] {
            out[idx] = -1;
        }
        Ok(out)
    }

    /// Sample a centered plaintext modulo `p`.
    pub fn random_plaintext(&self, rng: &mut Lcg) -> Vec<i64> {
        let half_width = ((self.p - 1) / 2) as i64;
        (0..self.group.order())
            .map(|_| rng.gen_range((2 * half_width + 1) as usize) as i64 - half_width)
            .collect()
    }

    /// Multiply two modular coefficient vectors.
    pub fn multiply_mod(
        &mut self,
        lhs: &[u64],
        rhs: &[u64],
        modulus: u64,
    ) -> Result<Vec<u64>, NtruError> {
        if lhs.len() != self.group.order() || rhs.len() != self.group.order() || modulus == 0 {
            return Err(NtruError::Shape);
        }

        if let Some(product) = self.try_fft_multiply(lhs, rhs, modulus)? {
            return Ok(product);
        }

        self.stats.record("naive:multiply");
        self.naive_multiply_mod(lhs, rhs, modulus)
    }

    /// Invert a modular coefficient vector.
    pub fn inverse_mod(&mut self, element: &[u64], modulus: u64) -> Result<Vec<u64>, NtruError> {
        if element.len() != self.group.order() || modulus == 0 {
            return Err(NtruError::Shape);
        }

        match self.group.kind() {
            GroupKind::Cyclic => {
                let inverse = cyclic_inverse_mod(element, modulus)?;
                self.stats.record("cyclic:invert");
                self.verify_inverse(element, &inverse, modulus)?;
                Ok(inverse)
            }
            GroupKind::Dihedral => {
                match dihedral_inverse_via_cyclic(element, self.group.n(), modulus) {
                    Ok(inverse) => {
                        self.stats.record("dihedral-cyclic:invert");
                        self.verify_inverse(element, &inverse, modulus)?;
                        Ok(inverse)
                    }
                    Err(error) => {
                        if is_prime(modulus) && gcd(2 * self.group.n() as u64, modulus) == 1 {
                            return Err(error);
                        }
                        let inverse = self.linear_solve_inverse(element, modulus)?;
                        self.stats.record("linear-solve:invert");
                        self.verify_inverse(element, &inverse, modulus)?;
                        Ok(inverse)
                    }
                }
            }
            GroupKind::Symmetric => {
                let inverse = self.linear_solve_inverse(element, modulus)?;
                self.stats.record("linear-solve:invert");
                self.verify_inverse(element, &inverse, modulus)?;
                Ok(inverse)
            }
        }
    }

    /// Try to invert, returning `None` for noninvertible elements.
    pub fn try_inverse_mod(&mut self, element: &[u64], modulus: u64) -> Option<Vec<u64>> {
        self.inverse_mod(element, modulus).ok()
    }

    /// Generate an NTRU-style private key.
    pub fn random_private_key(
        &mut self,
        rng: &mut Lcg,
        max_tries: usize,
    ) -> Result<PrivateKey, NtruError> {
        for attempts in 1..=max_tries {
            let f = self.random_ternary(rng, self.d + 1, self.d)?;
            let f_mod_p = signed_to_mod(&f, self.p);
            let f_p = match self.try_inverse_mod(&f_mod_p, self.p) {
                Some(inverse) => inverse,
                None => continue,
            };
            let f_mod_q = signed_to_mod(&f, self.q);
            let f_q = match self.try_inverse_mod(&f_mod_q, self.q) {
                Some(inverse) => inverse,
                None => continue,
            };
            return Ok(PrivateKey {
                f,
                f_p,
                f_q,
                attempts,
            });
        }
        Err(NtruError::NotInvertible)
    }

    /// Compute the public key `h = F_q g mod q`.
    pub fn public_key(&mut self, f_q: &[u64], g: &[i64]) -> Result<Vec<u64>, NtruError> {
        self.multiply_mod(f_q, &signed_to_mod(g, self.q), self.q)
    }

    /// Encrypt `m` as `e = p h r + m mod q`.
    pub fn encrypt(&mut self, m: &[i64], h: &[u64], r: &[i64]) -> Result<Vec<u64>, NtruError> {
        let hr = self.multiply_mod(h, &signed_to_mod(r, self.q), self.q)?;
        let m_mod_q = signed_to_mod(m, self.q);
        Ok(hr
            .iter()
            .zip(m_mod_q.iter())
            .map(|(value, message)| {
                add_mod(mul_mod(self.p % self.q, *value, self.q), *message, self.q)
            })
            .collect())
    }

    /// Decrypt by computing `F_p * center_lift(f e) mod p`.
    pub fn decrypt(&mut self, e: &[u64], f: &[i64], f_p: &[u64]) -> Result<Vec<i64>, NtruError> {
        let a_q = self.multiply_mod(&signed_to_mod(f, self.q), e, self.q)?;
        let a_lift = center_lift_vector(&a_q, self.q);
        let b_p = self.multiply_mod(f_p, &signed_to_mod(&a_lift, self.p), self.p)?;
        Ok(center_lift_vector(&b_p, self.p))
    }

    /// Run one randomized trial.
    pub fn trial(&mut self, rng: &mut Lcg, max_key_tries: usize) -> Result<TrialResult, NtruError> {
        let private = self.random_private_key(rng, max_key_tries)?;
        let g = self.random_ternary(rng, self.d, self.d)?;
        let r = self.random_ternary(rng, self.d, self.d)?;
        let m = self.random_plaintext(rng);
        let h = self.public_key(&private.f_q, &g)?;
        let e = self.encrypt(&m, &h, &r)?;
        let decrypted = self.decrypt(&e, &private.f, &private.f_p)?;

        let a_q = self.multiply_mod(&signed_to_mod(&private.f, self.q), &e, self.q)?;
        let a_lift = center_lift_vector(&a_q, self.q);
        let expected = add_i64_vectors(
            &scalar_mul_i64(self.p as i64, &self.multiply_integer(&g, &r)),
            &self.multiply_integer(&private.f, &m),
        );

        Ok(TrialResult {
            success: decrypted == m,
            no_wrap: a_lift == expected,
            key_attempts: private.attempts,
        })
    }

    /// Run repeated randomized trials.
    pub fn run_trials(&mut self, trials: usize, max_key_tries: usize, seed: u64) -> TrialSummary {
        let mut rng = Lcg::new(seed);
        let mut results = Vec::new();
        let mut keygen_failures = 0;
        for _ in 0..trials {
            match self.trial(&mut rng, max_key_tries) {
                Ok(result) => results.push(result),
                Err(NtruError::NotInvertible) => keygen_failures += 1,
                Err(_) => keygen_failures += 1,
            }
        }
        let successes = results.iter().filter(|result| result.success).count();
        let no_wraps = results.iter().filter(|result| result.no_wrap).count();
        let attempts: Vec<_> = results.iter().map(|result| result.key_attempts).collect();
        let avg_key_attempts = if attempts.is_empty() {
            None
        } else {
            Some(attempts.iter().sum::<usize>() as f64 / attempts.len() as f64)
        };
        let max_key_attempts = attempts.into_iter().max();
        TrialSummary {
            requested_trials: trials,
            completed_trials: results.len(),
            successes,
            no_wraps,
            keygen_failures,
            avg_key_attempts,
            max_key_attempts,
            backend_stats: self.stats.clone(),
        }
    }

    /// Multiply over the integers without modular reduction.
    pub fn multiply_integer(&self, lhs: &[i64], rhs: &[i64]) -> Vec<i64> {
        let mut out = vec![0; self.group.order()];
        for (lhs_index, lhs_value) in lhs.iter().enumerate() {
            if *lhs_value == 0 {
                continue;
            }
            for (rhs_index, rhs_value) in rhs.iter().enumerate() {
                if *rhs_value == 0 {
                    continue;
                }
                let product_index = self.group.multiply_index(lhs_index, rhs_index);
                out[product_index] += lhs_value * rhs_value;
            }
        }
        out
    }

    fn naive_multiply_mod(
        &self,
        lhs: &[u64],
        rhs: &[u64],
        modulus: u64,
    ) -> Result<Vec<u64>, NtruError> {
        let mut out = vec![0; self.group.order()];
        for (lhs_index, lhs_value) in lhs.iter().enumerate() {
            let lhs_value = *lhs_value % modulus;
            if lhs_value == 0 {
                continue;
            }
            for (rhs_index, rhs_value) in rhs.iter().enumerate() {
                let rhs_value = *rhs_value % modulus;
                if rhs_value == 0 {
                    continue;
                }
                let product_index = self.group.multiply_index(lhs_index, rhs_index);
                out[product_index] = add_mod(
                    out[product_index],
                    mul_mod(lhs_value, rhs_value, modulus),
                    modulus,
                );
            }
        }
        Ok(out)
    }

    fn try_fft_multiply(
        &mut self,
        lhs: &[u64],
        rhs: &[u64],
        modulus: u64,
    ) -> Result<Option<Vec<u64>>, NtruError> {
        match self.group.kind() {
            GroupKind::Dihedral => self.try_dihedral_fft_multiply(lhs, rhs, modulus),
            GroupKind::Symmetric => self.try_symmetric_fft_multiply(lhs, rhs, modulus),
            GroupKind::Cyclic => Ok(None),
        }
    }

    #[cfg(feature = "fft")]
    fn try_dihedral_fft_multiply(
        &mut self,
        lhs: &[u64],
        rhs: &[u64],
        modulus: u64,
    ) -> Result<Option<Vec<u64>>, NtruError> {
        if self.group.n() < 3
            || !self.group.n().is_power_of_two()
            || gcd(2 * self.group.n() as u64, modulus) != 1
        {
            return Ok(None);
        }
        if is_prime(modulus) && !(modulus - 1).is_multiple_of(self.group.n() as u64) {
            return Ok(None);
        }
        let omega = match fft_dihedral::root_of_unity(self.group.n(), modulus) {
            Ok(root) => root,
            Err(_) => return Ok(None),
        };
        let n = self.group.n();
        match fft_dihedral::dihedral_multiply_fft(
            &lhs[..n],
            &lhs[n..],
            &rhs[..n],
            &rhs[n..],
            modulus,
            omega,
        ) {
            Ok((rotations, reflections)) => {
                self.stats.record("fft-dihedral:multiply");
                Ok(Some(rotations.into_iter().chain(reflections).collect()))
            }
            Err(_) => Ok(None),
        }
    }

    #[cfg(not(feature = "fft"))]
    fn try_dihedral_fft_multiply(
        &mut self,
        _lhs: &[u64],
        _rhs: &[u64],
        _modulus: u64,
    ) -> Result<Option<Vec<u64>>, NtruError> {
        Ok(None)
    }

    #[cfg(feature = "fft")]
    fn try_symmetric_fft_multiply(
        &mut self,
        lhs: &[u64],
        rhs: &[u64],
        modulus: u64,
    ) -> Result<Option<Vec<u64>>, NtruError> {
        if !is_prime(modulus) || modulus <= self.group.n() as u64 {
            return Ok(None);
        }
        let plan = match fft_symmetric::SymmetricFft::new(self.group.n(), modulus) {
            Ok(plan) => plan,
            Err(_) => return Ok(None),
        };
        match plan.multiply(lhs, rhs) {
            Ok(product) => {
                self.stats.record("fft-symmetric:multiply");
                Ok(Some(product))
            }
            Err(_) => Ok(None),
        }
    }

    #[cfg(not(feature = "fft"))]
    fn try_symmetric_fft_multiply(
        &mut self,
        _lhs: &[u64],
        _rhs: &[u64],
        _modulus: u64,
    ) -> Result<Option<Vec<u64>>, NtruError> {
        Ok(None)
    }

    fn linear_solve_inverse(&self, element: &[u64], modulus: u64) -> Result<Vec<u64>, NtruError> {
        let matrix = self.left_multiplication_matrix(element, modulus);
        solve_modular_linear_system(&matrix, &self.one_mod(modulus), modulus)
    }

    fn left_multiplication_matrix(&self, element: &[u64], modulus: u64) -> Vec<Vec<u64>> {
        let mut matrix = vec![vec![0; self.group.order()]; self.group.order()];
        for (column, _) in self.group.elements().iter().enumerate() {
            for (index, coeff) in element.iter().enumerate() {
                let coeff = *coeff % modulus;
                if coeff == 0 {
                    continue;
                }
                let row = self.group.multiply_index(index, column);
                matrix[row][column] = add_mod(matrix[row][column], coeff, modulus);
            }
        }
        matrix
    }

    fn verify_inverse(
        &mut self,
        element: &[u64],
        inverse: &[u64],
        modulus: u64,
    ) -> Result<(), NtruError> {
        let identity = self.one_mod(modulus);
        if self.multiply_mod(element, inverse, modulus)? != identity {
            return Err(NtruError::NotInvertible);
        }
        if self.multiply_mod(inverse, element, modulus)? != identity {
            return Err(NtruError::NotInvertible);
        }
        Ok(())
    }
}

/// Return whether `value` is prime.
pub fn is_prime(value: u64) -> bool {
    if value < 2 {
        return false;
    }
    const SMALL: [u64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
    for prime in SMALL {
        if value == prime {
            return true;
        }
        if value.is_multiple_of(prime) {
            return false;
        }
    }
    let mut d = value - 1;
    let mut s = 0;
    while d.is_multiple_of(2) {
        s += 1;
        d /= 2;
    }
    for base in SMALL {
        if base >= value {
            continue;
        }
        let mut x = pow_mod(base, d, value);
        if x == 1 || x == value - 1 {
            continue;
        }
        let mut probably_prime = false;
        for _ in 0..s - 1 {
            x = mul_mod(x, x, value);
            if x == value - 1 {
                probably_prime = true;
                break;
            }
        }
        if !probably_prime {
            return false;
        }
    }
    true
}

/// Greatest common divisor.
pub fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let r = a % b;
        a = b;
        b = r;
    }
    a
}

/// Modular addition.
pub fn add_mod(lhs: u64, rhs: u64, modulus: u64) -> u64 {
    ((lhs as u128 + rhs as u128) % modulus as u128) as u64
}

/// Modular subtraction.
pub fn sub_mod(lhs: u64, rhs: u64, modulus: u64) -> u64 {
    ((lhs as u128 + modulus as u128 - rhs as u128) % modulus as u128) as u64
}

/// Modular multiplication.
pub fn mul_mod(lhs: u64, rhs: u64, modulus: u64) -> u64 {
    ((lhs as u128 * rhs as u128) % modulus as u128) as u64
}

/// Modular exponentiation.
pub fn pow_mod(mut base: u64, mut exp: u64, modulus: u64) -> u64 {
    let mut acc = 1 % modulus;
    base %= modulus;
    while exp > 0 {
        if exp & 1 == 1 {
            acc = mul_mod(acc, base, modulus);
        }
        base = mul_mod(base, base, modulus);
        exp >>= 1;
    }
    acc
}

fn extended_gcd(a: i128, b: i128) -> (i128, i128, i128) {
    if b == 0 {
        (a, 1, 0)
    } else {
        let (g, x, y) = extended_gcd(b, a % b);
        (g, y, x - (a / b) * y)
    }
}

/// Multiplicative inverse modulo `modulus`.
pub fn inv_mod(value: u64, modulus: u64) -> Option<u64> {
    let (g, x, _) = extended_gcd(value as i128, modulus as i128);
    if g != 1 {
        None
    } else {
        Some(x.rem_euclid(modulus as i128) as u64)
    }
}

/// Center lift one residue.
pub fn center_lift(value: u64, modulus: u64) -> i64 {
    let residue = value % modulus;
    if residue > modulus / 2 {
        residue as i64 - modulus as i64
    } else {
        residue as i64
    }
}

/// Center lift a vector.
pub fn center_lift_vector(values: &[u64], modulus: u64) -> Vec<i64> {
    values
        .iter()
        .map(|value| center_lift(*value, modulus))
        .collect()
}

/// Convert signed coefficients to canonical residues.
pub fn signed_to_mod(values: &[i64], modulus: u64) -> Vec<u64> {
    values
        .iter()
        .map(|value| (*value as i128).rem_euclid(modulus as i128) as u64)
        .collect()
}

fn add_i64_vectors(lhs: &[i64], rhs: &[i64]) -> Vec<i64> {
    lhs.iter().zip(rhs.iter()).map(|(a, b)| a + b).collect()
}

fn scalar_mul_i64(scalar: i64, values: &[i64]) -> Vec<i64> {
    values.iter().map(|value| scalar * value).collect()
}

/// Solve a modular linear system over `Z/modulus Z` using unit pivots.
pub fn solve_modular_linear_system(
    matrix: &[Vec<u64>],
    rhs: &[u64],
    modulus: u64,
) -> Result<Vec<u64>, NtruError> {
    let size = rhs.len();
    if matrix.len() != size || matrix.iter().any(|row| row.len() != size) || modulus <= 1 {
        return Err(NtruError::Shape);
    }
    let mut rows: Vec<Vec<u64>> = matrix
        .iter()
        .zip(rhs.iter())
        .map(|(row, rhs)| {
            let mut augmented: Vec<_> = row.iter().map(|value| value % modulus).collect();
            augmented.push(*rhs % modulus);
            augmented
        })
        .collect();

    for col in 0..size {
        let pivot_row = (col..size).find(|row| gcd(rows[*row][col], modulus) == 1);
        let Some(pivot_row) = pivot_row else {
            return Err(NtruError::NotInvertible);
        };
        if pivot_row != col {
            rows.swap(col, pivot_row);
        }
        let pivot_inv = inv_mod(rows[col][col], modulus).ok_or(NtruError::NotInvertible)?;
        for entry in &mut rows[col] {
            *entry = mul_mod(*entry, pivot_inv, modulus);
        }
        let pivot = rows[col].clone();
        for (row_index, row) in rows.iter_mut().enumerate() {
            if row_index == col {
                continue;
            }
            let factor = row[col];
            if factor == 0 {
                continue;
            }
            for (idx, entry) in row.iter_mut().enumerate().take(size + 1) {
                *entry = sub_mod(*entry, mul_mod(factor, pivot[idx], modulus), modulus);
            }
        }
    }

    Ok(rows.into_iter().map(|row| row[size] % modulus).collect())
}

/// Multiply in `(Z/mZ)[x] / (x^n - 1)`.
pub fn cyclic_convolution_mod(
    lhs: &[u64],
    rhs: &[u64],
    modulus: u64,
) -> Result<Vec<u64>, NtruError> {
    if lhs.len() != rhs.len() || modulus == 0 {
        return Err(NtruError::Shape);
    }
    let n = lhs.len();
    let mut out = vec![0; n];
    for (i, lhs_value) in lhs.iter().enumerate() {
        let lhs_value = *lhs_value % modulus;
        if lhs_value == 0 {
            continue;
        }
        for (j, rhs_value) in rhs.iter().enumerate() {
            let rhs_value = *rhs_value % modulus;
            if rhs_value == 0 {
                continue;
            }
            out[(i + j) % n] = add_mod(
                out[(i + j) % n],
                mul_mod(lhs_value, rhs_value, modulus),
                modulus,
            );
        }
    }
    Ok(out)
}

/// Apply `x -> x^-1` in `(Z/mZ)[x] / (x^n - 1)`.
pub fn cyclic_involution(values: &[u64], modulus: u64) -> Vec<u64> {
    let n = values.len();
    let mut out = vec![0; n];
    for (index, value) in values.iter().enumerate() {
        out[(n - index) % n] = *value % modulus;
    }
    out
}

/// Invert an element of `(Z/mZ)[x] / (x^n - 1)`.
pub fn cyclic_inverse_mod(values: &[u64], modulus: u64) -> Result<Vec<u64>, NtruError> {
    if values.is_empty() || modulus <= 1 {
        return Err(NtruError::Shape);
    }
    let normalized: Vec<_> = values.iter().map(|value| value % modulus).collect();
    if is_prime(modulus) {
        cyclic_inverse_prime_field(&normalized, modulus)
    } else {
        cyclic_inverse_linear_solve(&normalized, modulus)
    }
}

fn cyclic_inverse_linear_solve(values: &[u64], modulus: u64) -> Result<Vec<u64>, NtruError> {
    let n = values.len();
    let mut matrix = vec![vec![0; n]; n];
    for row in 0..n {
        for col in 0..n {
            matrix[row][col] = values[(row + n - col) % n] % modulus;
        }
    }
    let mut rhs = vec![0; n];
    rhs[0] = 1 % modulus;
    solve_modular_linear_system(&matrix, &rhs, modulus)
}

fn cyclic_inverse_prime_field(values: &[u64], modulus: u64) -> Result<Vec<u64>, NtruError> {
    let n = values.len();
    let mut modulus_poly = vec![0; n + 1];
    modulus_poly[0] = (modulus - 1) % modulus;
    modulus_poly[n] = 1;
    let inverse = polynomial_inverse_mod(values, &modulus_poly, modulus)?;
    Ok(polynomial_mod_xn_minus_1(&inverse, n, modulus))
}

fn polynomial_inverse_mod(
    values: &[u64],
    modulus_poly: &[u64],
    modulus: u64,
) -> Result<Vec<u64>, NtruError> {
    let mut old_r = trim_polynomial(modulus_poly, modulus);
    let mut r = trim_polynomial(values, modulus);
    let mut old_t = Vec::new();
    let mut t = vec![1];
    while !r.is_empty() {
        let (quotient, remainder) = polynomial_divmod(&old_r, &r, modulus)?;
        old_r = r;
        r = remainder;
        let next_t = polynomial_sub(&old_t, &polynomial_mul(&quotient, &t, modulus), modulus);
        old_t = t;
        t = next_t;
    }
    if old_r.len() != 1 || old_r[0] == 0 {
        return Err(NtruError::NotInvertible);
    }
    let scale = inv_mod(old_r[0], modulus).ok_or(NtruError::NotInvertible)?;
    Ok(old_t
        .into_iter()
        .map(|coeff| mul_mod(scale, coeff, modulus))
        .collect())
}

fn polynomial_divmod(
    numerator: &[u64],
    denominator: &[u64],
    modulus: u64,
) -> Result<(Vec<u64>, Vec<u64>), NtruError> {
    let mut remainder = trim_polynomial(numerator, modulus);
    let denominator = trim_polynomial(denominator, modulus);
    if denominator.is_empty() {
        return Err(NtruError::Shape);
    }
    if remainder.len() < denominator.len() {
        return Ok((Vec::new(), remainder));
    }
    let mut quotient = vec![0; remainder.len() - denominator.len() + 1];
    let denominator_lead_inv =
        inv_mod(*denominator.last().unwrap(), modulus).ok_or(NtruError::NotInvertible)?;
    while !remainder.is_empty() && remainder.len() >= denominator.len() {
        let shift = remainder.len() - denominator.len();
        let coeff = mul_mod(*remainder.last().unwrap(), denominator_lead_inv, modulus);
        quotient[shift] = coeff;
        for (index, denominator_coeff) in denominator.iter().enumerate() {
            remainder[shift + index] = sub_mod(
                remainder[shift + index],
                mul_mod(coeff, *denominator_coeff, modulus),
                modulus,
            );
        }
        remainder = trim_polynomial(&remainder, modulus);
    }
    Ok((trim_polynomial(&quotient, modulus), remainder))
}

fn polynomial_mul(lhs: &[u64], rhs: &[u64], modulus: u64) -> Vec<u64> {
    if lhs.is_empty() || rhs.is_empty() {
        return Vec::new();
    }
    let mut out = vec![0; lhs.len() + rhs.len() - 1];
    for (i, lhs_value) in lhs.iter().enumerate() {
        let lhs_value = *lhs_value % modulus;
        if lhs_value == 0 {
            continue;
        }
        for (j, rhs_value) in rhs.iter().enumerate() {
            let rhs_value = *rhs_value % modulus;
            if rhs_value == 0 {
                continue;
            }
            out[i + j] = add_mod(out[i + j], mul_mod(lhs_value, rhs_value, modulus), modulus);
        }
    }
    trim_polynomial(&out, modulus)
}

fn polynomial_sub(lhs: &[u64], rhs: &[u64], modulus: u64) -> Vec<u64> {
    let size = lhs.len().max(rhs.len());
    let mut out = vec![0; size];
    for (index, entry) in out.iter_mut().enumerate() {
        let lhs_value = lhs.get(index).copied().unwrap_or(0);
        let rhs_value = rhs.get(index).copied().unwrap_or(0);
        *entry = sub_mod(lhs_value, rhs_value, modulus);
    }
    trim_polynomial(&out, modulus)
}

fn polynomial_mod_xn_minus_1(values: &[u64], n: usize, modulus: u64) -> Vec<u64> {
    let mut out = vec![0; n];
    for (index, coeff) in values.iter().enumerate() {
        out[index % n] = add_mod(out[index % n], *coeff % modulus, modulus);
    }
    out
}

fn trim_polynomial(values: &[u64], modulus: u64) -> Vec<u64> {
    let mut out: Vec<_> = values.iter().map(|value| value % modulus).collect();
    while out.last() == Some(&0) {
        out.pop();
    }
    out
}

/// Invert in a dihedral group algebra via a cyclic-ring inverse.
pub fn dihedral_inverse_via_cyclic(
    element: &[u64],
    n: usize,
    modulus: u64,
) -> Result<Vec<u64>, NtruError> {
    if element.len() != 2 * n || n < 3 || modulus <= 1 {
        return Err(NtruError::Shape);
    }
    let rotations: Vec<_> = element[..n].iter().map(|value| value % modulus).collect();
    let reflections: Vec<_> = element[n..].iter().map(|value| value % modulus).collect();
    let rotation_bar = cyclic_involution(&rotations, modulus);
    let reflection_bar = cyclic_involution(&reflections, modulus);
    let rotation_norm = cyclic_convolution_mod(&rotations, &rotation_bar, modulus)?;
    let reflection_norm = cyclic_convolution_mod(&reflection_bar, &reflections, modulus)?;
    let determinant: Vec<_> = rotation_norm
        .iter()
        .zip(reflection_norm.iter())
        .map(|(lhs, rhs)| sub_mod(*lhs, *rhs, modulus))
        .collect();
    let determinant_inverse = cyclic_inverse_mod(&determinant, modulus)?;
    let inverse_rotations = cyclic_convolution_mod(&rotation_bar, &determinant_inverse, modulus)?;
    let neg_reflections: Vec<_> = reflections
        .iter()
        .map(|value| if *value == 0 { 0 } else { modulus - *value })
        .collect();
    let inverse_reflections =
        cyclic_convolution_mod(&neg_reflections, &determinant_inverse, modulus)?;
    Ok(inverse_rotations
        .into_iter()
        .chain(inverse_reflections)
        .collect())
}

/// Check that `Z[C_N]` multiplication is cyclic polynomial multiplication.
pub fn verify_cyclic_group_algebra_model(n: usize) -> Result<bool, NtruError> {
    let group = FiniteGroup::cyclic(n)?;
    let scheme = GroupAlgebraNtru::new(group, 3, 41, 2)?;
    for i in 0..n {
        for j in 0..n {
            if scheme.group.multiply_index(i, j) != (i + j) % n {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cyclic_model_matches_quotient_ring() {
        assert!(verify_cyclic_group_algebra_model(7).unwrap());
    }

    #[test]
    fn dihedral_cyclic_inverse_handles_small_fields() {
        let group = FiniteGroup::dihedral(16).unwrap();
        let mut fast = GroupAlgebraNtru::new(group.clone(), 3, 97, 3).unwrap();
        let dense = GroupAlgebraNtru::new(group, 3, 97, 3).unwrap();
        let mut rng = Lcg::new(20260506);
        for _ in 0..100 {
            let element = fast.random_ternary(&mut rng, 4, 3).unwrap();
            let element_mod = signed_to_mod(&element, 3);
            let dense_inverse = dense.linear_solve_inverse(&element_mod, 3);
            if let Ok(expected) = dense_inverse {
                let inverse = fast.inverse_mod(&element_mod, 3).unwrap();
                assert_eq!(inverse, expected);
                assert!(fast
                    .backend_stats()
                    .counts()
                    .contains_key("dihedral-cyclic:invert"));
                return;
            }
        }
        panic!("could not find an invertible sample");
    }

    #[test]
    fn symmetric_s3_round_trips() {
        let group = FiniteGroup::symmetric(3).unwrap();
        let mut scheme = GroupAlgebraNtru::new(group, 5, 67, 2).unwrap();
        let summary = scheme.run_trials(3, 50, 20260506);
        assert_eq!(summary.completed_trials, 3);
        assert_eq!(summary.successes, 3);
        assert_eq!(summary.no_wraps, 3);
    }

    #[test]
    fn dihedral_trials_round_trip() {
        let group = FiniteGroup::dihedral(8).unwrap();
        let mut scheme = GroupAlgebraNtru::new(group, 3, 97, 2).unwrap();
        let summary = scheme.run_trials(2, 100, 20260506);
        assert_eq!(summary.completed_trials, 2);
        assert_eq!(summary.successes, 2);
    }
}
