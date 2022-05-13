# Back projection

This program reads a file containing the neutron double scattering events, and reconstructs an image of the neutron source. Also available in [Github](https://github.com/fm140905/NPRE554-Back-projection-in-neutron-imaging.git)

Require [ROOT](https://root.cern/install/) installation.

## Compile
- Install `root`(v6.20 or above), `cmake` and `git`.
- Open a new terminal and run 
    ```bash
    git clone https://github.com/fm140905/NPRE554-Back-projection-in-neutron-imaging.git
    ```
- `cd` into the repo that you have just cloned.
- Run 
  ```bash
  mkdir -p build
  cd build
  cmake ..
  cmake --build .
  ```
  to compile the project.

## Test
- In the main folder run `bin/main`. You should see a status bar and an image should pop up when it reaches 100%.
## Use
- `coincidences_s1.txt` and `coincidences_s2.txt` contain the list of coincidence events. Make a new file of the same format from your simulation/experiment data. 
- Change the input parameters in `input/input.txt` accordingly.
- In the main folder run `bin/main`.