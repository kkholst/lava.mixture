'.onLoad' <- function(lib, pkg="lava.mixture") {
  desc <- utils::packageDescription(pkg)
  cat("Loading '", desc$Package, "' package...\n", sep="")
  cat("Version    : ", desc$Version, "\n", sep="")
}
