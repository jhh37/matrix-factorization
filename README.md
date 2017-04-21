# Secrets of Matrix Factorization (ICCV 2015)
This is a repository for the matrix factorization
project published in ICCV 2015.

## Contributors
- Je Hyeong Hong (jhh37@outlook.com)
- Andrew Fitzgibbon (awf@microsoft.com)

## A how-to-pull guide
Clone the repository using the following code
```
git clone --recursive https://github.com/jhh37/matrix-factorization
```

## Code
### List of algorithms
Below is a list of algorithms included.
  1. [JE] Our implementation of DRW (Damped Ruhe and Wedin) algorithms
  2. [JE] Our Ceres implementation of joint optimisation algorithms
  3. [CH] Chen's LM_S, LM_M series (our modification)
  4. [CO] Chen's LM_S, LM_M series (original: http://sist.sysu.edu.cn/~chenpei/)
  5. [DB] Del Bue's BALM (our implementation but not necessarily fast)
  6. [NB] Boumal's RTRMC (original: https://perso.uclouvain.be/nicolas.boumal/RTRMC/)
  7. [PG] Gotardo's CSF (original: http://www2.ece.ohio-state.edu/~gotardop/)
  8. [RC] Cabral's ALM series (our implementation but not necessarily fast)
  9. [TO] Okatani's Damped Wiberg (DW) (original: http://www.vision.is.tohoku.ac.jp/us/download/)

### A how-to-run demo
You may run the script [ run_demo.m ] inside the [
Code ] folder.

### Running experiments
Before starting, there are some issues you need to
watch out for:
  
  1. Prior to running any experiment, you must
have run [ Code/setup.m ], which adds relevant
paths and compile required mex files.
  2. All commands must run in the [ Code ]
directory.
  3. To run our Ceres-based algorithms (assuming
Ceres-solver from http://ceres-solver.org/ is
already installed), you need to compile the source
(which is pulled as a submodule) present in [
Code/Algorithms/JE/lrmf_ceres ].
    1. The CMake file is given.
    2. It may be convenient to compile in the [
Code/Algorithms/JE/lrmf_ceres/build ] directory
since this is git-ignored.
    3. Place the compiled executable inside the [
Code/Algorithms/JE ] folder with the name [
lrmf_ceres_exec ] or [ lrmf_ceres_exec.exe ].
    4. Reading sparse matrices is currently not
supported by [ LRMF-Ceres ].


## Video Spotlight
Our 1-min video spotlight can be found at
https://youtu.be/Kk0NIeNnc5M.

## Acknowledgement
- The work was supported by Microsoft and Toshiba
Research Europe. 
- We thank Roberto Cipolla, Christopher Zach,
Bamdev Mishra and the anonymous reviewers for
their comments.
- The conference travel was funded by the British
Machine Vision Association (BMVA), Cambridge
University Engineering Department (Rex Moir Fund),
Christ's College (University of Cambridge) and
Cambridge Philosophical Society.
