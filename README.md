# Secrets of Matrix Factorization (ICCV 2015)

# Contributors
- Je Hyeong Hong (jhh37@cam.ac.uk)
- Andrew Fitzgibbon (awf@microsoft.com)

# Abstract
Matrix factorization (or low-rank matrix
completion) with missing data is a key computation
in many computer vision and machine learning
tasks, and is also related to a broader class of
nonlinear optimization problems such as bundle
adjustment. The problem has received much
attention recently, with renewed interest in
variable-projection approaches, yielding dramatic
improvements in reliability and speed. However, on
a wide class of problems, no one approach
dominates, and because the various approaches have
been derived in a multitude of different ways, it
has been difficult to unify them. This paper
provides a unified derivation of a number of
recent approaches, so that similarities and
differences are easily observed. We also present a
simple meta-algorithm which wraps any existing
algorithm, yielding 100% success rate on many
standard datasets. Given 100% success, the focus
of evaluation must turn to speed, as 100% success
is trivially achieved if we do not care about
speed. Again our unification allows a number of
generic improvements applicable to all members of
the family to be isolated, yielding a unified
algorithm that outperforms our re-implementation
of existing algorithms, which in some cases
already outperform the original authors’
publicly available codes.

# Code
Please refer to [ readme.txt ] inside the [ Code ]
directory.

# Video Spotlight
Our 1-min video spotlight can be found at
https://youtu.be/Kk0NIeNnc5M.

# Acknowledgment
- The work was supported by Microsoft and Toshiba
Research Europe. We thank Roberto Cipolla,
Christopher Zach, Bamdev Mishra and the anonymous
reviewers for their invaluable comments.
- The conference travel was funded by the British
Machine Vision Association (BMVA), Cambridge
University Engineering Department (Rex Moir Fund),
Christ's College (University of Cambridge) and
Cambridge Philosophical Society.
