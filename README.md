# MTF-Lib
Multiple Target Filtering Library

I am planning to prepare a kind of multi-target filtering library which includes many types of multi-target filtering algorithms such as:
1. Traditional multi-target filters such as Global Nearest Neighbour (GNN), Joint Probabilistic Data Association Filter (JPDAF), Multiple Hypothesis Tracking (MHT), etc. Some implementations, particularly GNN and JPDA,  are available in https://github.com/USNavalResearchLaboratory/TrackerComponentLibrary For instance, look into AssignmentAlgorithms/singleScanUpdate.m and SampleCode/BasicTrackingExamples/demo2DDataAssociation.m. For MHT implementation, look into http://rehg.org/mht/ and https://ingemarcox.cs.ucl.ac.uk/?page_id=9
2. Random Finite Set (RFS)-based multi-target filtering algorithms such as Probability Hypothesis Density (PHD) filter, Cardinalized Probability Hypothesis Density (CPHD) filter, Cardinality Balanced Multi-Bernoulli (CB-MB) filter, Labeled Multi-Bernoulli (LMB) filter, Generalized Labeled Multi-Bernoulli (GLMB) filter, etc. This part is taken from http://ba-tuong.vo-au.com/codes.html.
3. Stochastic populations-based filter such as Hypothesized and Independent Stochastic Population (HISP) filter https://jeremiehoussineau.com/software/, etc.
4. Multiple target, multiple type filtering algorithm such as N-type PHD filter, etc

Different implementation schemes such as Kalman filter (KF), extended Kalman filter (EKF), uncented Kalman filter (UKF) and  sequential Monte Carlo (SMC) or particle filter (PF) of these multi-target filtering algorithms will be included. 

Single target filtering in different implementation schemes for a quick start such as KF, EKF, UKF and  PF will also be included.

The codes will be uploaded soon!
