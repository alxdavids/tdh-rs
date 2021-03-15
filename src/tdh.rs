use rand_core::{RngCore, OsRng};
use bit_vec::BitVec;
use curve25519_dalek::ristretto::RistrettoPoint;
use curve25519_dalek::constants::RISTRETTO_BASEPOINT_POINT as G;
use curve25519_dalek::scalar::Scalar;

const RISTRETTO_BYTE_LEN: usize = 32;

fn random_element() -> RistrettoPoint {
  let mut rng = OsRng;
  RistrettoPoint::random(&mut rng)
}

fn random_scalar() -> Scalar {
  let mut buf = [0u8; RISTRETTO_BYTE_LEN];
  random_bytes(&mut buf);
  Scalar::from_bytes_mod_order(buf)
}

pub fn random_bytes(out: &mut [u8; RISTRETTO_BYTE_LEN]){
  let mut rng = OsRng;
  rng.fill_bytes(out);
}

#[derive(Debug)]
pub struct Matrix {
  pub entries: [Vec<RistrettoPoint>; 2],
  pub n: u64,
}
impl Matrix {
  fn new(n: u64) -> Self {
    let inner = [Vec::with_capacity(n as usize), Vec::with_capacity(n as usize)];
    Self {
      entries: inner,
      n: n,
    }
  }
}

pub struct Hint(Vec<u8>);

pub struct HashKey {
  A: Matrix,
  n: u64,
}

pub struct EncodeKey {
  U: RistrettoPoint,
  B: Matrix,
  n: u64,
}

pub struct Trapdoor {
  s: Scalar,
  t: Scalar,
}

pub struct TrapdoorHash;
impl TrapdoorHash {
  pub fn sample(n: u64) -> HashKey {
    let mut A = Matrix::new(n);
    for b in 0..2 {
      for _ in 0..n {
        let r = random_element();
        A.entries[b].push(r)
      }
    }
    return HashKey { A: A, n: n };
  }

  pub fn generate(hk: &HashKey, idx: u64) -> (EncodeKey, Trapdoor) {
    let s = random_scalar();
    let t = random_scalar();
    let u = s*G;
    let mut B = Matrix::new(hk.n);
    for b in 0..2 {
      for j in 0..hk.n {
        let mut val = s*hk.A.entries[b][j as usize];
        if j == idx && b == 1 {
          val = val + t*G;
        }
        B.entries[b].push(val);
      }
    }
    (EncodeKey { U: u, B: B, n: hk.n}, Trapdoor { s: s, t: t })
  }

  pub fn hash(hk: &HashKey, x: &[u8], rho: &[u8; RISTRETTO_BYTE_LEN]) -> RistrettoPoint {
    let bv = BitVec::from_bytes(x);
    assert_eq!(bv.len(), hk.n as usize);
    let r = Scalar::from_bytes_mod_order(*rho);
    let mut h = r*G;
    for j in 0..hk.n {
      h = match bv.get(j as usize) {
        Some(false) => h + hk.A.entries[0][j as usize],
        Some(true) => h + hk.A.entries[1][j as usize],
        _ => panic!("Invalid response"),
      };
    }
    h
  }

  pub fn e(ek: &EncodeKey, x: &[u8], rho: &[u8; RISTRETTO_BYTE_LEN]) -> RistrettoPoint {
    let bv = BitVec::from_bytes(x);
    assert_eq!(bv.len(), ek.n as usize);
    let r = Scalar::from_bytes_mod_order(*rho);
    let mut hint = r*ek.U;
    for j in 0..ek.n {
      hint = match bv.get(j as usize) {
        Some(false) => hint + ek.B.entries[0][j as usize],
        Some(true) => hint + ek.B.entries[1][j as usize],
        _ => panic!("Invalid response"),
      };
    }
    hint
  }

  pub fn d(td: &Trapdoor, h: &RistrettoPoint) -> [RistrettoPoint; 2] {
    let e0 = td.s * h;
    let e1 = td.s * h + td.t*G;
    return [e0, e1]
  }
}

#[cfg(test)]
mod tests {
  #![allow(soft_unstable)]
  extern crate test;
  use test::Bencher;
  use bit_vec::BitVec;
  use rand::distributions::{Distribution, Uniform};
  use super::{TrapdoorHash, random_bytes};

  #[test]
  fn sanity_check() {
    let n = 16;
    let hk = TrapdoorHash::sample(n);
    for i in 0..n {
      let (ek, td) = TrapdoorHash::generate(&hk, i);
      let x = [1u8; 2];
      let bv = BitVec::from_bytes(&x);
      let mut rho = [0u8; 32];
      random_bytes(&mut rho);
      let h = TrapdoorHash::hash(&hk, &x, &rho);
      let e = TrapdoorHash::e(&ek, &x, &rho);
      let [e0, e1] = TrapdoorHash::d(&td, &h);
      match bv.get(i as usize) {
        Some(false) => assert_eq!(e, e0),
        Some(true) => assert_eq!(e, e1),
        _ => panic!("Unexpected failure"),
      }
    }
  }

  #[bench]
  fn benchmark_sample_16(b: &mut Bencher) {
    let n = 16;
    benchmark_sample(b, n);
  }

  #[bench]
  fn benchmark_sample_128(b: &mut Bencher) {
    let n = 32;
    benchmark_sample(b, n)
  }

  #[bench]
  fn benchmark_sample_256(b: &mut Bencher) {
    let n = 256;
    benchmark_sample(b, n)
  }
  
  #[bench]
  fn benchmark_sample_1024(b: &mut Bencher) {
    let n = 1024;
    benchmark_sample(b, n)
  }
  
  #[bench]
  fn benchmark_generate_16(b: &mut Bencher) {
    let n = 16;
    benchmark_generate(b, n);
  }

  #[bench]
  fn benchmark_generate_128(b: &mut Bencher) {
    let n = 32;
    benchmark_generate(b, n)
  }

  #[bench]
  fn benchmark_generate_256(b: &mut Bencher) {
    let n = 256;
    benchmark_generate(b, n)
  }
  
  #[bench]
  fn benchmark_generate_1024(b: &mut Bencher) {
    let n = 1024;
    benchmark_generate(b, n)
  }
  
  #[bench]
  fn benchmark_hash_16(b: &mut Bencher) {
    let n = 16;
    benchmark_hash(b, n);
  }

  #[bench]
  fn benchmark_hash_128(b: &mut Bencher) {
    let n = 32;
    benchmark_hash(b, n)
  }

  #[bench]
  fn benchmark_hash_256(b: &mut Bencher) {
    let n = 256;
    benchmark_hash(b, n)
  }
  
  #[bench]
  fn benchmark_hash_1024(b: &mut Bencher) {
    let n = 1024;
    benchmark_hash(b, n)
  }
  
  #[bench]
  fn benchmark_e_16(b: &mut Bencher) {
    let n = 16;
    benchmark_e(b, n);
  }

  #[bench]
  fn benchmark_e_128(b: &mut Bencher) {
    let n = 32;
    benchmark_e(b, n)
  }

  #[bench]
  fn benchmark_e_256(b: &mut Bencher) {
    let n = 256;
    benchmark_e(b, n)
  }
  
  #[bench]
  fn benchmark_e_1024(b: &mut Bencher) {
    let n = 1024;
    benchmark_e(b, n)
  }
  
  #[bench]
  fn benchmark_d_16(b: &mut Bencher) {
    let n = 16;
    benchmark_d(b, n);
  }

  #[bench]
  fn benchmark_d_128(b: &mut Bencher) {
    let n = 32;
    benchmark_d(b, n)
  }

  #[bench]
  fn benchmark_d_256(b: &mut Bencher) {
    let n = 256;
    benchmark_d(b, n)
  }
  
  #[bench]
  fn benchmark_d_1024(b: &mut Bencher) {
    let n = 1024;
    benchmark_d(b, n)
  }

  fn benchmark_sample(b: &mut Bencher, n: u64) {
    b.iter(|| {
      TrapdoorHash::sample(n)
    });
  }

  fn benchmark_generate(b: &mut Bencher, n: u64) {
    let hk = TrapdoorHash::sample(n);
    let dist = Uniform::from(0..n);
    let idx = dist.sample(&mut rand::thread_rng());
    b.iter(|| {
      TrapdoorHash::generate(&hk, idx)
    });
  }

  fn benchmark_hash(b: &mut Bencher, n: u64) {
    let hk = TrapdoorHash::sample(n);
    let mut x = Vec::new();
    random_bytes_arbitrary_length(&mut x, (n/8) as usize);
    let rho = [0u8; 32];
    b.iter(|| {
      TrapdoorHash::hash(&hk, &x, &rho)
    });
  }

  fn benchmark_e(b: &mut Bencher, n: u64) {
    let hk = TrapdoorHash::sample(n);
    let mut x = Vec::new();
    random_bytes_arbitrary_length(&mut x, (n/8) as usize);
    let dist = Uniform::from(0..n);
    let idx = dist.sample(&mut rand::thread_rng());
    let (ek, _) = TrapdoorHash::generate(&hk, idx);
    let rho = [0u8; 32];
    b.iter(|| {
      TrapdoorHash::e(&ek, &x, &rho)
    });
  }

  fn benchmark_d(b: &mut Bencher, n: u64) {
    let hk = TrapdoorHash::sample(n);
    let mut x = Vec::new();
    random_bytes_arbitrary_length(&mut x, (n/8) as usize);
    let dist = Uniform::from(0..n);
    let idx = dist.sample(&mut rand::thread_rng());
    let (_, td) = TrapdoorHash::generate(&hk, idx);
    let rho = [0u8; 32];
    let h = TrapdoorHash::hash(&hk, &x, &rho);
    b.iter(|| {
      TrapdoorHash::d(&td, &h)
    });
  }

  pub fn random_bytes_arbitrary_length(out: &mut Vec<u8>, byte_len: usize){
    let mut rng = rand::thread_rng();
    let dist = Uniform::from(0..255);
    let mut bytes = Vec::new();
    for _i in 0..byte_len {
      let byte = dist.sample(&mut rng);
      bytes.push(byte);
    }
    out.clear();
    out.extend_from_slice(&bytes);
  }
}