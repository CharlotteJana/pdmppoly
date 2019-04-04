context("transformMoment")

test_that("transformMoment works for n = 3, type = 'raw' and all central moments given", {
  
  mList <- structure(list(rawMomentOrders = diag(3),
                          rawMoments = list("m1", "m2", "m3"),
                          centralMomentOrders = expand.grid(list(0:1, 0:1, 0:2)),
                          centralMoments = as.list(letters[1:12])),
                     class = "momentList")
  
  expect_warning( # central moments of order 0 should be 1 and central moments of order 1 should be 0
    res1 <- transformMoment(order = c(1,1,2),
                            type = 'raw',
                            closure = "",
                            momentList = mList)
  )
  
  momentRaw <- quote(1 * ("m1"^1 * "m2"^1 * "m3"^2) * ("a") + 
                     1 * ("m1"^0 * "m2"^1 * "m3"^2) * ("b") + 
                     1 * ("m1"^1 * "m2"^0 * "m3"^2) * ("c") + 
                     1 * ("m1"^0 * "m2"^0 * "m3"^2) * ("d") + 
                     2 * ("m1"^1 * "m2"^1 * "m3"^1) * ("e") + 
                     2 * ("m1"^0 * "m2"^1 * "m3"^1) * ("f") + 
                     2 * ("m1"^1 * "m2"^0 * "m3"^1) * ("g") + 
                     2 * ("m1"^0 * "m2"^0 * "m3"^1) * ("h") + 
                     1 * ("m1"^1 * "m2"^1 * "m3"^0) * ("i") + 
                     1 * ("m1"^0 * "m2"^1 * "m3"^0) * ("j") + 
                     1 * ("m1"^1 * "m2"^0 * "m3"^0) * ("k") + 
                     1 * ("m1"^0 * "m2"^0 * "m3"^0) * ("l"))
  
  res2 <- structure(list(rawMomentOrders = rbind(diag(3), c(1, 1, 2)),
                         rawMoments = list("m1", "m2", "m3", momentRaw),
                         centralMomentOrders = mList$centralMomentOrders,
                         centralMoments = mList$centralMoments),
                    class = "momentList")
  
  expect_equal(res1, res2, check.attributes = FALSE)
})


test_that("transformMoment works for n = 3, type = 'central' and all raw moments given", {
  
  mList <- structure(list(rawMomentOrders = expand.grid(list(0:1, 0:1, 0:2)),
                          rawMoments = as.list(letters[1:12]),
                          centralMomentOrders = NULL,
                          centralMoments = NULL),
                     class = "momentList")
  
  res1 <- transformMoment(order = c(1,1,2),
                          type = 'central',
                          closure = "",
                          momentList = mList)
  
  momentCentr <- quote(1 * 1 * ("b"^1 * "c"^1 * "e"^2) * ("a") + 
                       -1 * 1 * ("b"^0 * "c"^1 * "e"^2) * ("b") + 
                       -1 * 1 * ("b"^1 * "c"^0 * "e"^2) * ("c") + 
                       1 * 1 * ("b"^0 * "c"^0 * "e"^2) * ("d") + 
                       -1 * 2 * ("b"^1 * "c"^1 * "e"^1) * ("e") + 
                       1 * 2 * ("b"^0 * "c"^1 * "e"^1) * ("f") + 
                       1 * 2 * ("b"^1 * "c"^0 * "e"^1) * ("g") + 
                       -1 * 2 * ("b"^0 * "c"^0 * "e"^1) * ("h") + 
                       1 * 1 * ("b"^1 * "c"^1 * "e"^0) * ("i") + 
                       -1 * 1 * ("b"^0 * "c"^1 * "e"^0) * ("j") + 
                       -1 * 1 * ("b"^1 * "c"^0 * "e"^0) * ("k") + 
                       1 * 1 * ("b"^0 * "c"^0 * "e"^0) * ("l"))
  
  res2 <- structure(list(rawMomentOrders = mList$rawMomentOrders,
                         rawMoments = mList$rawMoments,
                         centralMomentOrders = t(c(1, 1, 2)),
                         centralMoments = list(momentCentr)),
                    class = "momentList")
  
  expect_equal(res1, res2, check.attributes = FALSE)
})

test_that("transformMoment works for n = 2 and some central moments missing", {
  
  cOrders <- expand.grid(list(0:1, 0:2))
  cM <- list(0, 1, 1, "a", "b")
  
  rOrders <- expand.grid(list(0:2, 0:2))
  rOrders <- rOrders[rowSums(rOrders) > 0 & rowSums(rOrders) <= 2, ]
  
  mList <- structure(list(rawMomentOrders = rOrders,
                          rawMoments = list("m1", "A", "m2", "B", "C"),
                          centralMomentOrders = c,
                          centralMoments = cM),
                     class = "momentList")
  
  res1 <- transformMoment(order = c(1,2),
                          type = 'raw',
                          closure = "zero",
                          momentList = mList)
  
  res2 <- structure(list(rawMomentOrders = rbind(diag(3), c(1, 1, 2)),
                         rawMoments = list("m1", "m2", "m3", momentRaw),
                         centralMomentOrders = mList$centralMomentOrders,
                         centralMoments = mList$centralMoments),
                    class = "momentList")
  
  expect_equal(res1, res2, check.attributes = FALSE)
})

