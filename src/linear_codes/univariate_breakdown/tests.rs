#[cfg(test)]
mod tests {
    use crate::linear_codes::utils::calculate_t;
    use ark_bls12_377::Fq;

    #[test]
    fn test_calculate_t_with_good_parameters() {
        assert!(calculate_t::<Fq>(128, (4, 1), 2_usize.pow(32)).unwrap() < 200);
        assert!(calculate_t::<Fq>(256, (4, 1), 2_usize.pow(32)).unwrap() < 400);
    }
}