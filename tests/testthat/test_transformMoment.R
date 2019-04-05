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

test_that("transformMoment works with only the last moment(s) missing", {
  
  cOrders <- expand.grid(list(0:1, 0:2))
  cOrders <- cOrders[-6, ]
  cMoments <- list(0, 1, 1, "a", "b")
  
  rOrders <- expand.grid(list(0:2, 0:2))
  rOrders <- rOrders[rowSums(rOrders) > 0 & rowSums(rOrders) <= 2, ]
  rMoments <- list("m1", "A", "m2", "B", "C")
  
  mList <- structure(list(rawMomentOrders = rOrders,
                          rawMoments = rMoments,
                          centralMomentOrders = cOrders,
                          centralMoments = cMoments),
                     class = "momentList")
  
  #------- (1,2) missing in rOrders and cOrders, type = 'raw' ----------
  
  res1 <- transformMoment(order = c(1,2),
                          type = 'raw',
                          closure = "zero",
                          momentList = mList)
  
  momentRaw <- quote(1 * ("m1"^1 * "m2"^2) * (0) + 
                     1 * ("m1"^0 * "m2"^2) * (1) + 
                     2 * ("m1"^1 * "m2"^1) * (1) + 
                     2 * ("m1"^0 * "m2"^1) * ("a") + 
                     1 * ("m1"^1 * "m2"^0) * ("b") + 
                     1 * ("m1"^0 * "m2"^0) * (0))
  
  res2 <- structure(list(rawMomentOrders = rbind(rOrders, c(1, 2)),
                         rawMoments = append(rMoments, momentRaw),
                         centralMomentOrders = rbind(cOrders, c(1,2)),
                         centralMoments = append(cMoments, 0)),
                    class = "momentList")
  
  expect_equal(res1, res2, check.attributes = FALSE)
  
  #------- (1,2) missing in rOrders and cOrders, type = 'central' ----------
  
  res1 <- transformMoment(order = c(1,2),
                          type = 'central',
                          closure = "zero",
                          momentList = mList)
  
  res2 <- structure(list(rawMomentOrders = rOrders,
                         rawMoments = rMoments,
                         centralMomentOrders = rbind(cOrders, c(1,2)),
                         centralMoments = append(cMoments, 0)),
                    class = "momentList")
  
  expect_equal(res1, res2, check.attributes = FALSE)
  
  #------- (1,2) only missing in rOrders, type = 'raw' ----------
  
  cOrders <- expand.grid(list(0:1, 0:2))
  cMoments <- list(0, 1, 1, "a", "b", "c")
  mList$centralMomentOrders <- cOrders
  mList$centralMoments <- cMoments
  
  res1 <- transformMoment(order = c(1,2),
                          type = 'raw',
                          closure = "zero",
                          momentList = mList)
  
  momentRaw <- quote(1 * ("m1"^1 * "m2"^2) * (0) + 
                       1 * ("m1"^0 * "m2"^2) * (1) + 
                       2 * ("m1"^1 * "m2"^1) * (1) + 
                       2 * ("m1"^0 * "m2"^1) * ("a") + 
                       1 * ("m1"^1 * "m2"^0) * ("b") + 
                       1 * ("m1"^0 * "m2"^0) * ("c"))
  
  res2 <- structure(list(rawMomentOrders = rbind(rOrders, c(1, 2)),
                         rawMoments = append(rMoments, momentRaw),
                         centralMomentOrders = cOrders,
                         centralMoments = cMoments),
                    class = "momentList")
  
  expect_equal(res1, res2, check.attributes = FALSE)
  
  #------- (1,2) only missing in rOrders, type = 'central' ----------
  
  res1 <- transformMoment(order = c(1,2),
                          type = 'central',
                          closure = "zero",
                          momentList = mList)
  
  expect_equal(res1, mList, check.attributes = FALSE)
  
  #------- (1,2) only missing in cOrders, type = 'raw' ----------
  
  cOrders <- cOrders[-6, ]
  cMoments <- cMoments[-6]
  rOrders <- rbind(rOrders, c(1, 2))
  rMoments <- append(rMoments, "D")
  mList$centralMomentOrders <- cOrders
  mList$centralMoments <- cMoments
  mList$rawMomentOrders <- rOrders
  mList$rawMoments <- rMoments
  
  res1 <- transformMoment(order = c(1,2),
                          type = 'raw',
                          closure = "zero",
                          momentList = mList)
  
  expect_equal(res1, mList, check.attributes = FALSE)

  #------- (1,2) only missing in cOrders, type = 'central' ----------
  
  res1 <- transformMoment(order = c(1,2),
                          type = 'central',
                          closure = "zero",
                          momentList = mList)
  
  momentCentr <- quote(-1 * ("m1"^1 * "m2"^2) * ("0") + 
                       1 * ("m1"^0 * "m2"^2) * ("m1") + 
                       2 * ("m1"^1 * "m2"^1) * ("m2") + 
                       2 * ("m1"^0 * "m2"^1) * ("B") + 
                       -1 * ("m1"^1 * "m2"^0) * ("C") + 
                       1 * ("m1"^0 * "m2"^0) * ("D"))
  
  res2 <- structure(list(rawMomentOrders = rOrders,
                         rawMoments = rMoments,
                         centralMomentOrders = rbind(cOrders, c(1,2)),
                         centralMoments = append(cMoments, momentCentr)),
                    class = "momentList")
                    
  expect_equal(res1, res2, check.attributes = FALSE)
  
  #------- (1,2) not missing at all ----------
  
  cOrders <- expand.grid(list(0:1, 0:2))
  cMoments <- list(0, 1, 1, "a", "b", "c")
  mList$centralMomentOrders <- cOrders
  mList$centralMoments <- cMoments
  
  res1 <- transformMoment(order = c(1,2),
                          type = 'raw',
                          closure = "zero",
                          momentList = mList)
  
  res2 <- transformMoment(order = c(1,2),
                          type = 'central',
                          closure = "zero",
                          momentList = mList)
  
  expect_equal(res1, mList, check.attributes = FALSE)
  expect_equal(res2, mList, check.attributes = FALSE)
})

test_that("transformMoment works with a missing moment in the middle and at the end", {
  skip("not implemented yet")
})

test_that("transformMoment does not need central moments of order 0 and 1 as input", {
  skip("not implemented yet")
})