'.onAttach' <- function(lib, pkg="lava.mixture") {
  desc <- utils::packageDescription(pkg)
  packageStartupMessage("\nLoading '", desc$Package, "' package...\n",
                          "Version    : ", desc$Version, "\n",
                          "Overview: help(package=", desc$Package, ")\n"); 
}
