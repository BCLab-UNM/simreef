# SIMReeF

SIMReeF (Spatial Interaction Model of Reef Foraging) is a high-performance agent-based model (ABM) capable of explicitly simulating interactions between vast numbers of fish and algae within a 3D spatial grid. 

## To compile the code
```cd simforager```
```./build_hopper.sh release```

## To compile the C++ code

```module load boost```

```g++ --std=c++17 -I $BOOST_INC vonmises.cpp -o vonmises```
