# Energy Efficiency Fairness for Multi-Pair Wireless-Powered Relaying System

This repo is cloned from [this Code Ocean capsule](https://codeocean.com/capsule/4717808/tree/v1), which contains the code for the following scientific paper:

Kien-Giang Nguyen, Quang-Doanh Vu,  Le-Nam Tran, and Markku Juntti, "Energy Efficiency Fairness for Multi-Pair Wireless-Powered Relaying Systems," IEEE J. Sel. Areas Commun., vol. 37, no. 2, pp. 357-373, Feb. 2019.

## Instructions
The "**main_EEFairness_EHRelay_OW.m**" makes use of [Yalmip](https://yalmip.github.io/) as a parser and [MOSEK](https://www.mosek.com/) as the internal convex conic solver for speed. If you run the code on the Code Ocean platform, it automatically provides these packages (see environment/postInstall). You can of course download the code and run it on your local machine, but you need to install these packages yourself. You can substitute MOSEK with any convex conic solver such as SDPT3, Sedumi, Gurobi, and Cplex, to name a few.
