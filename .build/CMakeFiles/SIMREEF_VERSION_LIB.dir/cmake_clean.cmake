file(REMOVE_RECURSE
  "libSIMREEF_VERSION_LIB.pdb"
  "libSIMREEF_VERSION_LIB.a"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/SIMREEF_VERSION_LIB.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
