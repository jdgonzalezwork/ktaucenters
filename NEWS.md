# ktaucenters 1.0.0
This release focuses on accelerating clusters estimation by writing code in C++.

## Breaking changes

* The functions `ktaucenters_aux` and `denpoints` have been removed from the API.

* The funtion `ROBINDEN` have been renamed to `robinden`.
  
## New features

* Added `ktaucentersfast` with improved computational performance and an API similar to `kmeans` (by Douglas Carmona).
