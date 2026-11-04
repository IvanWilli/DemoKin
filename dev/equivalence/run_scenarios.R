# Run every scenario against one implementation and save outputs + timings.
#
# Usage: Rscript run_scenarios.R <package_dir> <out_dir> [scenario_regex] [fast]
#
# The package is loaded from source with pkgload::load_all() so the exact
# source tree is what gets timed, not whatever happens to be installed.
# scenario_regex: "" (or omitted) runs every scenario.
# fast: TRUE runs every scenario with `fast = TRUE` (the C++ path); FALSE or
# omitted runs the default R code. scenarios.R itself never mentions `fast`:
# the switch is made by shims defined below in the global environment, which
# the scenario closures resolve before the package exports.

args <- commandArgs(trailingOnly = TRUE)
pkg_dir <- args[1]
out_dir <- args[2]
pattern <- if (length(args) >= 3 && nzchar(args[3])) args[3] else NULL
fast <- length(args) >= 4 && toupper(args[4]) == "TRUE"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(pkgload::load_all(pkg_dir, export_all = FALSE, quiet = TRUE))

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
source(file.path(dirname(script_file), "scenarios.R"))

if (fast) {
  # every exported function that has a `fast` argument is shadowed by a shim
  # that forces fast = TRUE; rename_kin(), output_period_cohort_combination()
  # and plot_diagram() have no such argument and are left alone
  ns <- asNamespace("DemoKin")
  fast_exports <- Filter(function(fn) "fast" %in% names(formals(get(fn, ns))), getNamespaceExports(ns))
  for (fn in fast_exports) {
    assign(fn, local({ f <- get(fn, ns); function(...) f(..., fast = TRUE) }), envir = globalenv())
  }
  message("fast = TRUE for: ", paste(sort(fast_exports), collapse = ", "))
}

inp <- build_inputs()
scen <- build_scenarios(inp)
if (!is.null(pattern)) scen <- scen[grepl(pattern, names(scen))]

timings <- data.frame(scenario = names(scen), seconds = NA_real_, status = NA_character_, stringsAsFactors = FALSE)
for (i in seq_along(scen)) {
  nm <- names(scen)[i]
  cat(sprintf("[%s] %-28s ", format(Sys.time(), "%H:%M:%S"), nm))
  gc()
  t0 <- proc.time()[["elapsed"]]
  res <- tryCatch(scen[[nm]](), error = function(e) structure(list(message = conditionMessage(e)), class = "scenario_error"))
  el <- proc.time()[["elapsed"]] - t0
  timings$seconds[i] <- el
  timings$status[i] <- if (inherits(res, "scenario_error")) paste("ERROR:", res$message) else "ok"
  saveRDS(res, file.path(out_dir, paste0(nm, ".rds")))
  cat(sprintf("%8.2f s  %s\n", el, timings$status[i]))
  rm(res)
}
# a subset run (scenario_regex) updates its rows and keeps the others
tf <- file.path(out_dir, "_timings.csv")
if (!is.null(pattern) && file.exists(tf)) {
  old <- read.csv(tf, stringsAsFactors = FALSE)
  timings <- rbind(old[!old$scenario %in% timings$scenario, ], timings)
}
write.csv(timings, tf, row.names = FALSE)
cat("done\n")
