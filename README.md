# anfoamRec
Application for 3D reconstruction of foam

## Installation
Dependencies:
vtk-python

The code depends on third-party applications, all must be callable from python/shell, place their executables for example to /usr/local/bin:
- `neper` for tessellation (http://neper.sourceforge.net/)
- `se_api` for modifying foam structure, converting geo to Surface Evolver format, analyzing tools (https://github.com/kolarji2/SE_api)
- `meshconv` for converting between stl and ply format (http://www.patrickmin.com/meshconv/)
- `binvox` for generating voxel structure from ply format (http://www.patrickmin.com/binvox/)
- `vox_fill` for filling empty holes in foam structure in voxel format
- `evolver` for structure relaxation (http://facstaff.susqu.edu/brakke/evolver/evolver.html)
PackingGeneration - for dense packing available on https://github.com/VasiliBaranov/packing-generation

## Usage
$python main.py -c examples/example_sp1.json -gra
$python main.py -c examples/example_el1.json -gra
