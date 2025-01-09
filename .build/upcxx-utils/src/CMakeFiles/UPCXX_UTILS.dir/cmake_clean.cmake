file(REMOVE_RECURSE
  "libUPCXX_UTILS.pdb"
  "libUPCXX_UTILS.a"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/UPCXX_UTILS.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
