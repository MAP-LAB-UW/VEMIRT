# VEMIRT: A Package for High-Dimensional IRT Models

VEMIRT is created to assist researchers to conduct exploratory and
confirmatory multidimensional item response theory (MIRT) analysis and
cooresponding item differential functioning (DIF) analysis. The core
computation engine of VEMIRT is a family of Gaussian Variational EM
algorithms that are considerably more efficient than currently available
algorithms in other software packages, especially when the number of
latent factors exceeds four.

## Identifying the number of factors

[`pa_poly`](https://MAP-LAB-UW.github.io/VEMIRT/reference/pa_poly.md)
identifies the number of factors via parallel analysis.

## Exploratory factor analysis

- [`E2PL_gvem_rot`](https://MAP-LAB-UW.github.io/VEMIRT/reference/E2PL_gvem_rot.md)
  conducts M2PL Analysis with post-hoc rotation (Promax & CF-Quartimax)

- [`E2PL_gvem_lasso`](https://MAP-LAB-UW.github.io/VEMIRT/reference/E2PL_gvem_lasso.md)
  conducts M2PL Analysis with Lasso penalty

- [`E2PL_gvem_adaptlasso`](https://MAP-LAB-UW.github.io/VEMIRT/reference/E2PL_gvem_adaptlasso.md)
  conducts M2PL Analysis with adaptive Lasso penalty

- [`E2PL_iw`](https://MAP-LAB-UW.github.io/VEMIRT/reference/C2PL_iw.md)
  conducts importance sampling to correct bias for M2PL analysis

- [`E3PL_sgvem_rot`](https://MAP-LAB-UW.github.io/VEMIRT/reference/E3PL_sgvem_rot.md)
  conducts stochastic GVEM to further improve the computational
  effficiency for exploratory M3PL analysis

- [`E3PL_sgvem_lasso`](https://MAP-LAB-UW.github.io/VEMIRT/reference/E3PL_sgvem_lasso.md)
  conducts M3PL Analysis with Lasso penalty

- [`E3PL_sgvem_adaptlasso`](https://MAP-LAB-UW.github.io/VEMIRT/reference/E3PL_sgvem_adaptlasso.md)
  conducts M3PL Analysis with adaptive Lasso penalty

- [`MGRM_gvem`](https://MAP-LAB-UW.github.io/VEMIRT/reference/MGRM_gvem.md)
  conducts GVEM for the multidimensional graded response model with
  post-hoc rotation

- [`MGPCM_gvem`](https://MAP-LAB-UW.github.io/VEMIRT/reference/MGPCM_gvem.md)
  conducts GVEM for the multidimensional partial credit model with
  post-hoc rotation

## Confirmatory factor analysis

- [`C2PL_gvem`](https://MAP-LAB-UW.github.io/VEMIRT/reference/C2PL_gvem.md)
  conducts GVEM for confirmatory M2PL analysis

- [`C2PL_bs`](https://MAP-LAB-UW.github.io/VEMIRT/reference/C2PL_bs.md)
  conducts bootstrap sampling to correct bias and produce standard
  errors for confirmatory M2PL analysis

- [`C2PL_iw`](https://MAP-LAB-UW.github.io/VEMIRT/reference/C2PL_iw.md)
  conducts importance sampling to correct bias for M2PL analysis

- [`C2PL_iw2`](https://MAP-LAB-UW.github.io/VEMIRT/reference/C2PL_iw2.md)
  conducts IW-GVEM for confirmatory M2PL analysis (alternative
  implementation to
  [`C2PL_iw`](https://MAP-LAB-UW.github.io/VEMIRT/reference/C2PL_iw.md))

- [`C3PL_sgvem`](https://MAP-LAB-UW.github.io/VEMIRT/reference/C3PL_sgvem.md)
  conducts stochastic GVEM for confirmatory M3PL analysis

- [`MGRM_gvem`](https://MAP-LAB-UW.github.io/VEMIRT/reference/MGRM_gvem.md)
  conducts GVEM for the multidimensional graded response model

- [`MGPCM_gvem`](https://MAP-LAB-UW.github.io/VEMIRT/reference/MGPCM_gvem.md)
  conducts GVEM for the multidimensional partial credit model

## Differential item functioning analysis

- [`D1PL_em`](https://MAP-LAB-UW.github.io/VEMIRT/reference/D1PL_em.md)
  conducts DIF analysis for M1PL models using EM algorithms

- [`D1PL_gvem`](https://MAP-LAB-UW.github.io/VEMIRT/reference/D1PL_gvem.md)
  conducts DIF analysis for M1PL models using GVEM algorithms

- [`D2PL_em`](https://MAP-LAB-UW.github.io/VEMIRT/reference/D2PL_em.md)
  conducts DIF analysis for M2PL models using EM algorithms

- [`D2PL_pair_em`](https://MAP-LAB-UW.github.io/VEMIRT/reference/D2PL_pair_em.md)
  conducts DIF analysis for 2PL models using EM algorithms with group
  pairwise truncated \\L_1\\ penalty

- [`D2PL_gvem`](https://MAP-LAB-UW.github.io/VEMIRT/reference/D2PL_gvem.md)
  conducts DIF analysis for M2PL models using GVEM algorithms

- [`D2PL_lrt`](https://MAP-LAB-UW.github.io/VEMIRT/reference/D2PL_lrt.md)
  conducts DIF analysis for M2PL models using the likelihood ratio test

## Shiny apps for VEMIRT

- [`shinyVEMIRT`](https://MAP-LAB-UW.github.io/VEMIRT/reference/shinyVEMIRT.md)
  Run the shiny app for VEMIRT

- [`DIFdashboard`](https://MAP-LAB-UW.github.io/VEMIRT/reference/DIFdashboard.md)
  Run the shiny app for DIF Dashboard

## See also

Useful links:

- <https://MAP-LAB-UW.github.io/VEMIRT>

- <https://github.com/MAP-LAB-UW/VEMIRT>

## Author

**Maintainer**: Weicong Lyu <weiconglyu@um.edu.mo>
([ORCID](https://orcid.org/0000-0002-8159-2161))

Authors:

- Yijun Cheng <chengxb@uw.edu>
  ([ORCID](https://orcid.org/0000-0002-0671-9193))

- Jiaying Xiao <jxiao6@uw.edu>
  ([ORCID](https://orcid.org/0000-0001-9513-6477))

- He Ren <heren@uw.edu> ([ORCID](https://orcid.org/0000-0003-1208-7089))

- Ruoyi Zhu <zhux0445@uw.edu>
  ([ORCID](https://orcid.org/0000-0002-8159-2161))

- Gongjun Xu <gongjun@umich.edu>
  ([ORCID](https://orcid.org/0000-0003-4023-5413))

- Chun Wang <wang4066@uw.edu>
  ([ORCID](https://orcid.org/0000-0003-2695-9781))
