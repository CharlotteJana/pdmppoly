context("momentList")
#t1 tests für validate_momentList

mList <- structure(list(rawMomentOrders = rbind(1:3, diag(3), 4:6),
                        rawMoments = list("A", "m1", "m2", "m3", "B"),
                        centralMomentOrders = expand.grid(list(0:2, 0:2, 0:2)),
                        centralMoments = as.list(letters, "lastElement")),
                   class = "momentList")

test_that("cov works", {
  skip("besser cov testen als cov.momentList")
  
  cov <- cov.momentList(mList)
  expect_equal(cov, list(list("c", "e", "k"),
                         list("e", "g", "m"),
                         list("k", "m", "s")))
})

test_that("mean works", {
  mean <- mean(mList)
  expect_equal(mean, list("m1", "m2", "m3"))
})