test_that("rho_opt works", {
    x1 <- c(1.0, 2.0, 2.5, 3.0, 4.0)
    x2 <- c(-1.0, -2.0, -2.5, -3.0, -4.0)
    expected_output <- c(0.1538462, 0.6153846, 0.9072139, 1.0, 1.0)
    expect_equal(round(rho_opt(x = x1, c = 1.0), 6), round(expected_output, 6))
    expect_equal(round(rho_opt(x = x2, c = 1.0), 6), round(expected_output, 6))
})

test_that("psi_opt works", {
    x1 <- c(1.0, 2.0, 2.5, 3.0, 4.0)
    x2 <- c(-1.0, -2.0, -2.5, -3.0, -4.0)
    expected_output_1 <- c(0.3076923, 0.6153846, 0.4421154, -2.186285e-15, 0.0)
    expected_output_2 <- c(-0.3076923, -0.6153846, -0.4421154, 2.186285e-15, 0.0)
    expect_equal(round(psi_opt(x = x1, c = 1.0), 6), round(expected_output_1, 6))
    expect_equal(round(psi_opt(x = x2, c = 1.0), 6), round(expected_output_2, 6))
})

test_that("derpsi_opt works", {
    x1 <- c(1.0, 2.0, 2.5, 3.0, 4.0)
    x2 <- c(-1.0, -2.0, -2.5, -3.0, -4.0)
    expected_output <- c(0.3076923, 0.3076923, -0.965462, 0.0, 0.0)
    expect_equal(round(derpsi_opt(x = x1, c = 1.0), 6), round(expected_output, 6))
    expect_equal(round(derpsi_opt(x = x2, c = 1.0), 6), round(expected_output, 6))
})