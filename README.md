# DRESIM
Dynamic Resilient Infrastructure Simulation and Management (http://tempdemoproject.s3-website-us-east-1.amazonaws.com/)
A demo for the cascading failure and resilient repair. 
See the User Manual file for the interface description, code explanation, and the usage instruction. 

The file SingleRepair allows one node repair at a time. The file SimultaneousRepair allows the repair of multiple nodes at the same time. 

Yalmip (https://yalmip.github.io/) is required to run the files. Proper LP solvers such as Gurobi (https://www.gurobi.com/) are required. 

Factored graph and variable elimination are applied to obtain the suboptimal decision for large-scale interdependent infrastructure network. 

An overview video of the project can be found at https://www.youtube.com/watch?v=J9ugLhDMSmo

References: 
1. A Factored MDP Approach to Optimal Mechanism Design for Resilient Large-Scale Interdependent Critical Infrastructures (https://ieeexplore.ieee.org/document/8064531)
2. DISTRIBUTED AND OPTIMAL RESILIENT PLANNING OF LARGE-SCALE INTERDEPENDENT CRITICAL INFRASTRUCTURES 
(https://ieeexplore.ieee.org/abstract/document/8632399)
