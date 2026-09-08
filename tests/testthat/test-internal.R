test_that("EMPTY is <EMPTY>", {
    expect_equal(EMPTY, "<EMPTY>")
})

test_that("UNKNOWN is <UNKNOWN>", {
    expect_equal(UNKNOWN, "<UNKNOWN>")
})

test_that("explode() works", {
    expect_equal(explode("A,B,C"), c("A", "B", "C"))
    expect_equal(explode("A, B, C"), c("A", "B", "C"))
    expect_equal(explode("A, B, C", ", "), c("A", "B", "C"))
    expect_equal(explode("A, B, C", ","), c("A", "B", "C"))
    expect_equal(explode("A; B; C", ";"), c("A", "B", "C"))
})

test_that("n_ending() works", {
    expect_equal(n_ending("a+", "+"), 1)
    expect_equal(n_ending("a++", "+"), 2)
    expect_equal(n_ending("a+b+++", "+"), 3)
    expect_equal(n_ending("a-b--", "-"), 2)
})

test_that("revert_cell_name() works", {
    expect_equal(revert_cell_name("CD4"), "CD4")
    expect_equal(revert_cell_name("CD4.1"), "CD4.1")
    expect_equal(revert_cell_name("CD4..1"), "CD4")
    expect_equal(revert_cell_name("C.D.4..1"), "C.D.4")
})

test_that("do_call() works with non-named args", {
    expect_equal(do_call(sum, list(1, 2, 3)), 6)
})

test_that("do_call() works with named args", {
    expect_equal(do_call(sum, list(x = 1, y = 2, z = 3)), 6)
})
