# MTF-Lib
This repository contains resources for multiple target filtering (MTF). It includes different types of filtering algorithms such as:
1. Traditional multi-target filters such as Global Nearest Neighbour (GNN), Joint Probabilistic Data Association Filter (JPDAF), Multiple Hypothesis Tracking (MHT), etc. Some implementations, particularly GNN and JPDA,  are available in https://github.com/USNavalResearchLaboratory/TrackerComponentLibrary For instance, look into AssignmentAlgorithms/singleScanUpdate.m and SampleCode/BasicTrackingExamples/demo2DDataAssociation.m. For MHT implementation, look into http://rehg.org/mht/ and https://ingemarcox.cs.ucl.ac.uk/?page_id=9
2. Random Finite Set (RFS)-based multi-target filtering algorithms such as Probability Hypothesis Density (PHD) filter, Cardinalized Probability Hypothesis Density (CPHD) filter, Cardinality Balanced Multi-Bernoulli (CB-MB) filter, Labeled Multi-Bernoulli (LMB) filter, Generalized Labeled Multi-Bernoulli (GLMB) filter, etc. For this part, look into http://ba-tuong.vo-au.com/codes.html.
3. Stochastic populations-based filter such as Hypothesized and Independent Stochastic Population (HISP) filter https://jeremiehoussineau.com/software/, etc.
4. Multiple target, multiple type filtering algorithm such as [N-type PHD filter](https://github.com/nathanlem1/MTF-Lib/tree/master/N-type-PHD-Filter/N-type-GM-PHD-Filter), etc
5. My PHD filter implementations in both [Matlab](https://github.com/nathanlem1/MTF-Lib/tree/master/GM-PHD-Filter/GM-PHD-Filter-Matlab) and [Python](https://github.com/nathanlem1/MTF-Lib/tree/master/GM-PHD-Filter/GM-PHD-Filter-Python) are included.

Different implementation schemes such as Kalman filter (KF), extended Kalman filter (EKF), uncented Kalman filter (UKF) and  sequential Monte Carlo (SMC) or particle filter (PF) of these multi-target filtering algorithms are included.  

Single target filtering in different implementation schemes for a quick start such as KF, EKF, UKF and  PF are also included. Particularly, KF is implemented in both Matlab and Python.

It is also worth looking at the following links for some related information on https://github.com/USNavalResearchLaboratory/TrackerComponentLibrary, https://github.com/chisyliu/Sensor-Fusion-and-Nonlinear-Filtering-SSY345 and https://github.com/rlabbe/Kalman-and-Bayesian-Filters-in-Python.
