# make my own quick package
library(outsider.devtools)
# construct a skeleton file structure for the module
module_path <- module_skeleton(program_name = 'echo', flpth = getwd())
# check the file structure
module_check(flpth = module_path)
# look-up key identifying names: R package name, Docker image name
module_identities(flpth = module_path)
# build the R package and Docker image
module_build(flpth = module_path, tag = 'latest')
# test the module
module_test(flpth = module_path)
