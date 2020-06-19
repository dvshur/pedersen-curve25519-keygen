mod polynom;
mod pow;

use curve25519_dalek::{constants, edwards::EdwardsPoint, scalar::Scalar};
use digest::Digest;
use ed25519_dalek::Sha512;
use polynom::Polynom;
use pow::Pow;
use rand::{CryptoRng, RngCore};
use std::convert::TryInto;

// Secret sharing parameters
// Parties number
const N: usize = 5;
// Threshold
const T: usize = 3;

const PREFIX: [u8; 32] = [255u8; 32];

#[allow(non_snake_case)]
fn main() {
    // Parties indexes
    let xs: [Scalar; N] = [
        Scalar::from(1u8),
        Scalar::from(2u8),
        Scalar::from(3u8),
        Scalar::from(4u8),
        Scalar::from(5u8),
    ];

    // can be arbitrary point in a subgroup, but equal for all participants
    let PEDERSEN_BASE_POINT =
        &Scalar::from_bytes_mod_order(PREFIX) * &constants::ED25519_BASEPOINT_TABLE;

    let mut csprng = rand::thread_rng();

    // generate (x_i, h_i) — a key pair for each participant
    let key_pairs = {
        let mut arr = [(Scalar::default(), EdwardsPoint::default()); N];
        for i in 0..N {
            arr[i] = generate_key_pair(&mut csprng);
        }
        arr
    };

    // generate blinders
    let r_arr = {
        let mut arr = [Scalar::default(); N];
        for i in 0..N {
            arr[i] = Scalar::random(&mut csprng);
        }
        arr
    };

    // pedersen commits
    let commits: [EdwardsPoint; N] = {
        let mut arr = [EdwardsPoint::default(); N];
        for i in 0..N {
            arr[i] = &key_pairs[i].1 + &(&r_arr[i] * &PEDERSEN_BASE_POINT);
        }
        arr
    };

    // unblind commits
    let commits_unblinded = {
        let mut arr = [EdwardsPoint::default(); N];
        for i in 0..N {
            arr[i] = &commits[i] - &(&r_arr[i] * &PEDERSEN_BASE_POINT);
        }
        arr
    };

    // resulting public key
    let public_key: EdwardsPoint = commits_unblinded.iter().sum();

    // random polynoms for each participants
    let polynoms = {
        let mut arr = Vec::with_capacity(N);
        for i in 0..N {
            arr.push(Polynom::random(&mut csprng, &key_pairs[i].0, T - 1));
        }
        arr
    };

    // F_i_j — exponents of polynom coefficients of each party
    let F_arr = {
        let mut arr = [[EdwardsPoint::default(); T]; N];
        // F_i_0 for each party equals h_i, a public key
        for i in 0..N {
            arr[i][0] = key_pairs[i].1
        }
        for i in 0..N {
            for j in 1..(T - 1) {
                arr[i][j] = &polynoms[i].coeffs[j] * &constants::ED25519_BASEPOINT_TABLE;
            }
        }
        arr
    };

    // s_i_j — polynom values at participants indices
    // should be encrypted in production
    let s_arr = {
        let mut arr = [[Scalar::default(); N]; N];
        for i in 0..N {
            for j in 0..N {
                arr[i][j] = polynoms[i].at(&xs[j]);
            }
        }
        arr
    };

    // each P_i verifies that shares received from all P_j
    // are consistent with previously published F_i_j values
    for i in 0..N {
        for j in 0..T {
            if i == j {
                continue;
            }

            let j_index = xs[j];

            let sum: EdwardsPoint = (0..T)
                .map(|i| j_index.pow(i as u64))
                .zip(F_arr[i].iter())
                .map(|(j_pow, F_i_j)| &j_pow * F_i_j)
                .sum();

            assert_eq!(sum, &s_arr[i][j] * &constants::ED25519_BASEPOINT_TABLE);
        }
    }
    // if the check above fails, P_i broadcasts that an error has been found,
    // publishes s_i_j and the signature and aborts

    // each participant reconstructs his shares
    let shares = {
        let mut arr = [Scalar::zero(); N];
        for i in 0..N {
            for j in 0..N {
                arr[i] += s_arr[j][i];
            }
        }
        arr
    };

    // check correctness
    let reconstruct_participants: &[Scalar; T] = xs[0..T].try_into().unwrap();
    let reconstruct_shares: &[Scalar; T] = shares[0..T].try_into().unwrap();

    let secret = shamir_reconstruct(reconstruct_participants, reconstruct_shares);

    assert_eq!(&secret * &constants::ED25519_BASEPOINT_TABLE, public_key);
}

fn lagrange_coeffs_at_zero(xs: &[Scalar; T]) -> [Scalar; T] {
    let mut cs = [Scalar::one(); T];

    for i in 0..T {
        for j in 0..T {
            if i != j {
                cs[i] *= xs[j] * (xs[j] - xs[i]).invert();
            }
        }
    }

    cs
}

fn generate_key_pair<T>(csprng: &mut T) -> (Scalar, EdwardsPoint)
where
    T: CryptoRng + RngCore,
{
    let mut bytes = [0u8; 32];
    csprng.fill_bytes(&mut bytes);

    bytes = Sha512::digest(&bytes)[00..32].try_into().unwrap();

    // do a conversion as per RFC
    bytes[0] &= 248;
    bytes[31] &= 127;
    bytes[31] |= 64;

    let sk = Scalar::from_bits(bytes);

    (sk, &sk * &constants::ED25519_BASEPOINT_TABLE)
}

fn shamir_reconstruct(xs: &[Scalar; T], shares: &[Scalar; T]) -> Scalar {
    let lagrange_coeffs = lagrange_coeffs_at_zero(xs);

    let mut res = Scalar::zero();
    for i in 0..T {
        res += lagrange_coeffs[i] * shares[i];
    }

    res
}
