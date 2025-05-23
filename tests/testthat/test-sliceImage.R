require(testthat)
require(Cardinal)

context("sliceImage")

test_that("sliceImage", {

	path <- CardinalIO::exampleImzMLFile("continuous")
	s <- readImzML(path)
	s1 <- s
	s2 <- s
	runNames(s1) <- "run0"
	runNames(s2) <- "run1"
	s <- cbind(s1, s2)

	mz <- mz(s)[1:2]
	rs1 <- sliceImage(s, i=1)
	rs2 <- sliceImage(s, i=1:2)
	rs3 <- sliceImage(s, mz=mz[1L])
	rs4 <- sliceImage(s, mz=mz)
	rs5 <- sliceImage(s, mz=mz, run=1)

	expect_equivalent(dim(rs1), c(3,3,2))
	expect_equivalent(dim(rs2), c(3,3,2,2))
	expect_equivalent(dim(rs3), c(3,3,2))
	expect_equivalent(dim(rs4), c(3,3,2,2))
	expect_equal(rs1, rs2[,,,1L])
	expect_equal(rs1, rs3)
	expect_equal(rs2, rs4)
	expect_equal(rs3, rs4[,,,1L])
	expect_equal(rs5, rs4[,,1L,])

})
