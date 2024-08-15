# Recast-CMS-SUS-19-013

This repository contains the complete high accuracy recreation of the <a href="https://link.springer.com/article/10.1007%2FJHEP09%282020%29149">CMS-SUS-19-013</a> search.

We start with the baseline selection with requirements 1-7 from the original publication. These selections are applied online, immediately after the monte-carlo generation on the computing cluster to severely shrink the required storage space.

All analysis steps that do not concern the reweighting to real detector data are followed. The result are the ingredients that are needed for the statistical analysis.

We verify our simulated data by reproducing the statistically indistinguishable exclusion bounds. The events are then used in the publication on <a href="https://journals.aps.org/prd/abstract/10.1103/PhysRevD.109.096031">resonant anomaly detection</a>.

Dependencies are: ROOT, MadGraph, Delphes, ExRootAnalysis, Boost and GSL
