# Regression tests for the cleaning helpers; requires only base R.
# Usage: Rscript tests/test_read_and_clean.R
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) != 1L) stop("Run these tests with Rscript.")
repo_root <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)))

run_tests <- function() {
  scripts <- c("identify gene reaction associations.R", "Fig2.R", "Fig3.R",
               "FigS2.R", "FigS3.R", "FigS4.R", "FigS5.R")
  # Extract the actual helper without loading plotting packages or running analysis.
  helpers <- lapply(scripts, function(script) {
    expressions <- parse(file.path(repo_root, "Code", "R", script))
    assignments <- Filter(function(expr) {
      is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
        identical(expr[[2L]], as.name("read_and_clean"))
    }, as.list(expressions))
    stopifnot(length(assignments) == 1L)
    environment <- new.env(parent = globalenv())
    eval(assignments[[1L]], envir = environment)
    environment$read_and_clean
  })
  fixture <- tempfile(fileext = ".csv")
  on.exit(unlink(fixture), add = TRUE)
  clean_fixture <- function(helper, data) {
    write.table(data, fixture, sep = ",", row.names = FALSE, na = "NA")
    helper(fixture)
  }

  for (i in seq_along(helpers)) {
    helper <- helpers[[i]]
    # Distinguish a missing row from a missing column; retain zeros and partial NA.
    mixed <- data.frame(ID = c("missing", "partial", "zero", "negative", "infinite"),
                        rxn_a = c(NA, 1, 0, -2, Inf),
                        rxn_b = c(NA, NA, 2, 3, -Inf),
                        empty = rep(NA_real_, 5))
    result <- clean_fixture(helper, mixed)
    stopifnot(is.data.frame(result),
              identical(result$ID, c("partial", "zero", "negative")),
              identical(names(result), c("ID", "rxn_a", "rxn_b")),
              identical(as.numeric(result$rxn_a), c(1, 0, 0)),
              is.na(result$rxn_b[[1L]]),
              identical(as.numeric(result$rxn_b[-1L]), c(2, 3)))

    single <- clean_fixture(helper, data.frame(ID = "valid", rxn = 0, empty = NA_real_))
    stopifnot(is.data.frame(single), identical(dim(single), c(1L, 2L)),
              identical(names(single), c("ID", "rxn")), single$rxn == 0)

    empty <- clean_fixture(helper, data.frame(ID = c("a", "b"),
                                             rxn = c(NA_real_, NA_real_)))
    stopifnot(is.data.frame(empty), identical(dim(empty), c(0L, 1L)),
              identical(names(empty), "ID"))
    id_only <- clean_fixture(helper, data.frame(ID = c("a", "b")))
    stopifnot(is.data.frame(id_only), identical(dim(id_only), c(0L, 1L)),
              identical(names(id_only), "ID"))
    no_rows <- clean_fixture(helper, data.frame(ID = character(), rxn = numeric()))
    stopifnot(is.data.frame(no_rows), identical(dim(no_rows), c(0L, 1L)),
              identical(names(no_rows), "ID"))

    # Check retained IDs, reaction names and values against an independent
    # matrix-based calculation on both archived inference tables.
    for (filename in c("Max_flux_screen_8.csv", "Max_flux_screen_8_Re.csv")) {
      path <- file.path(repo_root, "Results", "screens", filename)
      raw <- read.table(path, header = TRUE, sep = ",")
      values <- as.matrix(raw[-1L])
      values[is.infinite(values)] <- NA
      values[!is.na(values) & values < 0] <- 0
      keep_rows <- apply(values, 1L, function(row) any(!is.na(row)))
      kept <- values[keep_rows, , drop = FALSE]
      keep_cols <- apply(kept, 2L, function(column) any(!is.na(column)))
      expected <- kept[, keep_cols, drop = FALSE]
      result <- helper(path)
      stopifnot(is.data.frame(result),
                identical(result[[1L]], raw[[1L]][keep_rows]),
                identical(names(result), c(names(raw)[1L], colnames(expected))),
                isTRUE(all.equal(unname(as.matrix(result[-1L])), unname(expected),
                                 check.attributes = FALSE, tolerance = 0)))
      if (i == 1L) {
        cat(sprintf("%s: %d -> %d rows; %d -> %d reaction columns\n",
                    filename, nrow(raw), nrow(result), ncol(raw) - 1L, ncol(result) - 1L))
      }
    }
    cat("PASS:", scripts[[i]], "\n")
  }
  # Exercise the actual downstream expressions after all-NA columns are removed.
  collect_calls <- function(expr, predicate) {
    found <- if (is.call(expr) && predicate(expr)) list(expr) else list()
    if (is.call(expr) || is.expression(expr) || is.pairlist(expr)) {
      children <- as.list(expr)
      for (i in seq_along(children)) {
        if (identical(children[[i]], quote(expr = ))) next
        found <- c(found, collect_calls(children[[i]], predicate))
      }
    }
    found
  }
  for (script in c("Fig3.R", "FigS2.R", "FigS3.R", "FigS4.R", "FigS5.R")) {
    expressions <- parse(file.path(repo_root, "Code", "R", script))
    counts <- collect_calls(expressions, function(expr) {
      identical(expr[[1L]], as.name("<-")) &&
        identical(expr[[2L]], as.name("tmp_count"))
    })
    stopifnot(length(counts) == 2L)
    for (count in counts) {
      environment <- new.env(parent = globalenv())
      environment$current_data <- data.frame(ID = c("a", "b", "c"), rxn = c(0, NA, 2))
      environment$rxn <- "rxn"
      for (row in 1:3) {
        environment$rowindex <- row
        stopifnot(eval(count, environment) == c(1L, 0L, 0L)[row])
      }
      environment$rxn <- "removed"
      stopifnot(eval(count, environment) == 0L)
    }
  }
  expressions <- parse(file.path(repo_root, "Code", "R", scripts[[1L]]))
  blocked <- collect_calls(expressions, function(expr) {
    identical(expr[[1L]], as.name("map2_lgl"))
  })
  stopifnot(length(blocked) == 3L)
  for (call in blocked) {
    body <- call[[4L]][[2L]]
    environment <- new.env(parent = globalenv())
    environment$T1auto <- environment$T1hetero <- environment$T1mixo <-
      data.frame(EnzymeID = "a", rxn = 0)
    environment$.x <- "a"
    environment$.y <- "rxn"
    stopifnot(isTRUE(eval(body, environment)))
    environment$.y <- "removed"
    stopifnot(identical(eval(body, environment), FALSE))
    environment$.y <- "rxn"
    environment$.x <- "absent"
    stopifnot(identical(eval(body, environment), FALSE))
  }
  cat("PASS: downstream zero-flux counts and missing-column guards\n")
  cat("All cleaning regression tests passed.\n")
}
run_tests()
