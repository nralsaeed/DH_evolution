# DH_evolution
Box model and parameter study for the Mars water and DH evolution from 3Ga to present

DH_modules.py contains all the different versions of the model as callable functions and the parameter space routine 
              that can be imported to notebooks etc. as well as some other useful functions 

              included are:
              GEL(x)
              rGEL(x)
              box_model
              loop_box_model
              var_outg
              var_loss
              var_box_model
              loop_var_box_model
              parameter_study
              parameter_study_var
              find_init

DH_plots_generation.ipynb contains all the plotting routines and data reduction 

var_R32_results_average_MVN.csv   contains the results of the parameter study with R=0.32
                                  Columns : 
                                  Initial: initial amount of water at 3.3Ga in m GEL
                                  Remainder: current amount of water today in m GEL
                                  Outgassed: total amount of water outgassed over the past 3 by in m GEL
                                  DH: final D/H ratio in units of VSMOW
                                  Escape: rate of escape of H2O in grams/year averaged out over the 3 by
                                  Rows: 
                                  succesful runs of the model

var_R32_results_average_MVN.csv   contains the results of the parameter study with R=0
                                  Columns : 
                                  Initial: initial amount of water at 3.3Ga in m GEL
                                  Remainder: current amount of water today in m GEL
                                  Outgassed: total amount of water outgassed over the past 3 by in m GEL
                                  DH: final D/H ratio in units of VSMOW
                                  Escape: array of rates of escape of H2O in grams/year
                                  Rows: 
                                  succesful runs of the model
