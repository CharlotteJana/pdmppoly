context("momentList")
#t1 tests für validate_momentList

test_that("mean works", {
  
  mList <- structure(list(rawMomentOrders = rbind(1:3, diag(3), 4:6),
                          rawMoments = list("A", "m1", "m2", "m3", "B"),
                          centralMomentOrders = expand.grid(list(0:2, 0:2, 0:2)),
                          centralMoments = as.list(letters, "lastElement")),
                     class = "momentList")
  
  mean <- mean(mList)
  expect_equal(mean, list("m1", "m2", "m3"))
})

test_that("cov works for n > 1", {
  skip("besser cov testen als cov.momentList")
  
  mList <- structure(list(rawMomentOrders = rbind(1:3, diag(3), 4:6),
                          rawMoments = list("A", "m1", "m2", "m3", "B"),
                          centralMomentOrders = expand.grid(list(0:2, 0:2, 0:2)),
                          centralMoments = as.list(letters, "lastElement")),
                     class = "momentList")
  
  cov <- cov.momentList(mList)
  expect_equal(cov, list(list("c", "e", "k"),
                         list("e", "g", "m"),
                         list("k", "m", "s")))
})

test_that("cov uses transformMoment appropriately", {
  skip("not implemented yet")
  
  mList <- momentList(rawMomentOrders = expand.grid(list(0:2, 0:2)),
                      rawMoments = as.list(letters[1:9]))
  
  cov <- cov.momentList(mList)
  expect_equal(cov, list(list("c", "e", "k"),
                         list("e", "g", "m"),
                         list("k", "m", "s")))
})

test_that("cov works for n = 1", {
  skip("besser cov testen als cov.momentList")
  
  mList <- structure(list(rawMomentOrders = cbind(c(3, 1, 5)),
                          rawMoments = list("A", "m1", "B"),
                          centralMomentOrders = cbind(c(0, 1, 2)),
                          centralMoments = list(1, 0, "C")),
                     class = "momentList")
  
  cov <- cov.momentList(mList)
  expect_equal(cov, list(list("C")))
})
