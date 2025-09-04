# 🕸️ FoodWeb

This repository provides tools analyze community food webs with respect to:  
- **Stability**  
- **Sensitive and influential species**  
- **Impacts of perturbations**  

---

## 🔄 Data Transformation

**`extractcommunites.jl`**  
- Extracts food web and species information from excel files into a list of `Taxa` (see [`types.jl`](./types.jl)).  
- Defines `Community` objects (see [`types.jl`](./types.jl)) that can be analyzed.  
- Saves these communities to a file (default: `communities.jld2`).  

---

## 📊 Application

**`analysis.jl`**  
- Loads a community from `communities.jld2`.  
- Analyzes food web stability using functions from [`backgroundFunctions.jl`](./backgroundFunctions.jl).  
- Provides the option to:  
  - Detect **sensitive and influential taxa**, with visualization (.png files)
  - Detect the **impact of a given perturbation** on taxa, with visualization (.png files)

---

## 🚀 Getting Started

1. Clone this repository.  
2. Install required Julia packages (`Pkg.add("JLD2")`, etc.).  
3. Run `extractcommunites.jl` to build and save communities.  
4. Run `analysis.jl` to perform stability and sensitivity analyses. 
