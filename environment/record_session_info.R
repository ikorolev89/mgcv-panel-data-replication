#!/usr/bin/env Rscript
packages <- c("mgcv", "ggplot2", "plm", "sandwich", "lmtest")
for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Install required package: ", p)
}
dir.create("environment", showWarnings = FALSE)
versions <- data.frame(package = packages, version = vapply(packages, function(p) as.character(packageVersion(p)), character(1)))
write.csv(versions, "environment/package-versions.csv", row.names = FALSE)
for (p in packages) suppressPackageStartupMessages(library(p, character.only = TRUE))
writeLines(capture.output(sessionInfo()), "environment/session-info.txt")
ip <- installed.packages()
deps <- unique(c(packages, unlist(tools::package_dependencies(packages, db = ip, recursive = TRUE))))
write.csv(data.frame(package = deps, version = ip[deps, "Version"]), "environment/dependency-versions.csv", row.names = FALSE)
print(versions, row.names = FALSE)
