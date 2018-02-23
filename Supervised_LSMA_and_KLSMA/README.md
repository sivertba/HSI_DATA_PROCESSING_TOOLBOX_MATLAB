#  Supervised LSMA and KLSMA

Linear spectral mixture analysis (LSMA) is a mathematical theory that models a data samples as
linear mixtures of a finite number of basic spectral constituents with appropriate weights from
which data samples can be solved by finding these weights via a linear inverse problem. Linear
spectral unmixing is one of its applications. 

It assumes that there are p basic material substances
\{m_j\}^p_{j=1}
that can be used to represent data sample vectors in linear forms with their corresponding
abundance fractions 
\{a_j\}^p_{j=1}
that are unknown parameters. The spectral unmixing is then performed
by finding best estimates of
\{\hat{a}_j\}^p_{j=1}
, and referred to as unmixed abundance
fractions. In real applications, two physical constraints must be imposed on the used linear
mixing model, which are abundance sum-to-one constraint (ASC), 
\sum^p_{j=1} a_j = 1 
, and abundance
nonnegativity constraint (ANC), a_j >= 0 for all 1 <= j <= p. 
Three LSMA-based techniques developed in Chang (2003a) have been widely used for spectral unmixing. 

These are unconstrained orthogonal subspace projection (OSP)/least squares OSP, partially ANC-constrained method, Least
Squares nonnegativity-constrained least squares (NCLS) and a fully abundance-constrained
method, Constrained Least Squares (FCLS). Since OSP/LSOSP, NCLS, and FCLS are linear techniques,
they may have difficulty solving linear nonseparable problems. To address this issue, these
three techniques are further extended to their kernel-based counterparts, called Kernel OSP/LSOSP
(KOSP/KLSOSP), Kernel NCKS (KNCLS), and Kernel FCLS (KFCLS)