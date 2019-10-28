**Study: Stout, A., A. Van Stelten-Carlson, H. Marquis, M. Ballou, B. Reilly, G. Loneragan, K. Nightingale, and R. Ivanek. Public health impact of foodborne exposure to naturally occurring virulence-attenuated Listeria monocytogenes: inference from mouse and mathematical models. Interface Focus (under review)**
*(link to the article will be added upon article acceptance for publication)* 

This study tests hypotheses that (i) natural exposure to virulence-attenuated (vA) strains of L. monocytogenes through food can confer protective immunity against listeriosis attributable to fully- virulent (fV) strains and (ii) current food safety measures to minimize exposure to both of L. monocytogenes strains may have adverse population level outcomes.  To test these hypotheses, we evaluated the host response to L. monocytogenes in a mouse infection model and through mathematical modeling in a human population. Here we include the 5 elements representing the experimental data and mathematical model underlying the study:
* A: experimental data described in the manuscript,
* B and C: modeling code described in the manuscript and supplementary materials, and 
* D and E: model predictions described in the manuscript and supplementary materials. 

**The 5 elements are explained below:**

**A. Mouse Data.xlsb**
Data from the vaccine challenge experiments described in Table 3, in which L. monocytogenes recovery from internal organs was evaluated. The results of analysis of this data is depicted in Figure 2. 

**B. Two_strain_Lm_model_with_intermediate_immune_state.R** 
A system of ordinary differential equations describing a two-strain model of L. monocytogenes exposures and cross-protection to predict the annual number of listeriosis cases evaluated at the steady state. The model structure is presented in Figure 1. 

**C. MonteCarloSimulation_Lmonocytogenes.R**
The model coded in “Two_strain_Lm_model_with_intermediate_immune_state.R” with Monte Carlo simulation for parameters representing: (1) Probability of colonization by either strain for an individual in the susceptible compartment (kappa_SC); (2) Probability of colonization by either strain for an individual in the protected compartment (kappa_PC), where kappa_PC can be smaller or larger than kappa_SC as defined by the variable parameter “Factor x”; and (3) Rate of immunity loss. The model predicts the annual number of listeriosis cases evaluated at the steady state in a total of 10,000 iterations. 

**D. Model Outputs.xlsb**
Predictions from the model coded in “Two_strain_Lm_model_with_intermediate_immune_state.R” evaluated for different numbers of exposures to L. monocytogenes annually, including a constant ratio of virulence attenuated (VA) and fully virulent (FV) strains and evaluated under the parameter sets estimated under the “wide” and “narrow” calibration windows, as well under a special case of a simplified model without immune boosting. The model predictions are shown in Figure 3.  
The dataset also shows predictions from the same mathematical model evaluated for different numbers of exposures to L. monocytogenes annually, but assuming policies where either (1) virulence attenuated strains are targeted while fully virulent strains are held constant, or (2) fully virulent strains are targeted while virulence attenuated strains are held constant. All scenarios are evaluated for parameter sets estimated under the “wide” and “narrow” calibration window. The model predictions are shown in Figure 4.   

**E. MonteCarlo with factor range 0.5 to 1.5.xlsb**
Predictions of the annual number of listeriosis cases evaluated at the steady state in the Monte Carlo model coded in “MonteCarloSimulation_Lmonocytogenes.R” run with 10,000 iterations. The model predictions are shown in Figure S5 in the Supplementary Materials.  
