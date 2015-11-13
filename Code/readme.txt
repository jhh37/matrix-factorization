# LIST OF ALGORITHMS
Below is a list of algorithms included in the Algorithm folder.

1. [JE] Our implementation of DRW (Damped Ruhe and Wedin) algorithms
2. [JE] Our Ceres implementation of joint optimisation algorithms
3. [CH] Chen’s LM_S, LM_M series (our modification)
4. [CO] Chen’s LM_S, LM_M series (original: http://sist.sysu.edu.cn/~chenpei/)
5. [DB] Del Bue’s BALM (our implementation but not necessarily fast)
6. [NB] Boumal’s RTRMC (original: https://perso.uclouvain.be/nicolas.boumal/RTRMC/)
7. [PG] Gotardo’s CSF (original: http://www2.ece.ohio-state.edu/~gotardop/)
8. [RC] Cabrals’ ALM series (our implementation but not necessarily fast)
9. [TO] Okatani’s Damped Wiberg (DW) (original: http://www.vision.is.tohoku.ac.jp/us/download/)

# HOW TO RUN
For now, you can run the script [ run_demo.m ] and view its source to observe the code structure.

# RUNNING EXPERIMENTS
Before starting, there are some issues you may need to watch out for:

1. All commands must run in the [ Code ] directory, where this readme.txt file
is present.

2. To run our Ceres-implemented algorithms (assuming Ceres-solver from http://ceres-solver.org/ is already installed), you need to compile the source (inserted as a submodule) present in [ Algorithms/JE/ceres_lrmc ]. Both CMake and Visual Studio project files are given, and you may like to build inside [ Algorithms/JE/ceres_lrmc/build ] since it is git-ignored. Then place the executable inside the folder [ Algorithms/JE ] with the name [ je_ceres_lrmc ] or [ je_ceres_lrmc.exe ].

3. Our Ceres implementation on NET (Netflix_2k) is not yet supported as reading sparse matrices is currently not supported by [ Ceres-LRMC ].

4. Executing run_demo.m will call for setup.m which should setup all the
external libraries and compile mex files.
