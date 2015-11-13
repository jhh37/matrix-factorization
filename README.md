# Secrets of Matrix Factorization (ICCV '15)
Je Hyeong Hong (jhh37@cam.ac.uk) and Andrew Fitzgibbon (awf@microsoft.com)

# Abstract
Matrix factorization (or low-rank matrix completion) with missing data is a key computation in many computer vision and machine learning tasks, and is also related to a broader class of nonlinear optimization problems such as bundle adjustment. The problem has received much atten- tion recently, with renewed interest in variable-projection approaches, yielding dramatic improvements in reliability and speed. However, on a wide class of problems, no one approach dominates, and because the various approaches have been derived in a multitude of different ways, it has been difficult to unify them. This paper provides a uni- fied derivation of a number of recent approaches, so that similarities and differences are easily observed. We also present a simple meta-algorithm which wraps any existing algorithm, yielding 100% success rate on many standard datasets. Given 100% success, the focus of evaluation must turn to speed, as 100% success is trivially achieved if we do not care about speed. Again our unification allows a num- ber of generic improvements applicable to all members of the family to be isolated, yielding a unified algorithm that outperforms our re-implementation of existing algorithms, which in some cases already outperform the original au- thors’ publicly available codes.

# Code
Please refer to the readme text file inside the Code folder.

# Video
Our 1-minute video spotlight can be found at https://youtu.be/Kk0NIeNnc5M.
