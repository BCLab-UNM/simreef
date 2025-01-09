file(REMOVE_RECURSE
  "libUPCXX_UTILS_VERSION_LIB.pdb"
  "libUPCXX_UTILS_VERSION_LIB.a"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/UPCXX_UTILS_VERSION_LIB.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
