'''
This file contains all the different versions of the model as callable functions and the parameter space routine 
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

'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tnrange, tqdm_notebook
from multiprocessing import Pool

#Defining constants and useful conversions that will be used in the code

#The terrestrial DH ratio 
DH_before = 1.56e-4                   #terrestial value of D/H (Vienna Standard Mean Ocean Water) 

#function that converts H2O amount from m GEL to grams
def GEL(x):
    ''' Calculates the amount of H2O in grams that is in a global layer equivalent given in meters 
    Inputs: 
    x = GEL in meters
    Outputs:
    amount of water in grams
    '''
    SA = 1.448e14 #m^2 surface area of Mars
    den = 1e6 #g/m^3 density of water
    water = x*den*SA  #height * density * Surface Area
    return water

#function that converts H2O amount from grams to m GEL
def rGEL(y):
    '''Calculates the global equivalent of water in meters given amount of water in grams
    Reverse of GEL(x)
    Inputs: 
    y = amount of water in grams
    Outputs:
    GEL in meters'''
    SA = 1.448e14 #m^2 surface area of Mars
    den = 1e6 #g/m^3 density of water
    gel = y/(den*SA) # weight/(density*SurfaceArea)
    return gel

def box_model(Iw,X,V,R):
    '''The box model 
    Takes in inputs for boundary conditions (initial water in the past, water today, 
    outgassing rate and fractionation factor) and calculates an H2O escape rate to satisfy 
    the boundary conditions. Calculates the time dependant HDO escape rate and calculates 
    the amount of water in each box at each timestep. Outputs a dataframe with all the values 
    in each timestep as well as the final enrichment today and the calculated escape rate of water. 
    Inputs:
    Iw: The initial reservoir of water in meters GEL
    X: Water budget today in meters GEL
    V: total amount of water outgassed in the past tf (3.3 billion) years
    R: the fractionation factor of D escape compared to H escape 
    	(efficiency of D escape relative to H which ranges from 1e-5 to 0.32)
    Outputs:
    DF: Dataframe that gives the following information at each timestep:
        'Space'   :   amount of H2O in space that has escaped
        'Exch_Res':   amount of H2O remaining in the exchangeable reservoir
        'HDO'     :   amount of HDO remaining in the exchangeable reservoir
        'DH'      :   the D/H ratio in the exchangeable reservoir
        't'       :   time (in years) elapsed since the beginning of the model run (from tf (3 bn) yrs ago)
    enrichment: The final calculated enrichment at the end of the box model run ('today')
    H2O_loss: The average rate of escape of H2O (in grams per year)
    '''
    
    #start with defining the starting time in Ga
    tf = 3.3e9 #total time in years (we start at 3bn yrs ago)
    #amount of timesteps the model has to run through (the more the better but the more computation time needed)
    n = 100000 #total steps 
    #calculatin the time increment
    dt = tf/n
    
    #initial conditions
    
    #The terrestrial DH ratio we will be measuring enrichment against
    DH_before = 1.56e-4                   #terrestial value of D/H 
    
    #Initial enrichment (as a multiple of the terrestial value) tf (3 billion) years ago 
    	#(acquired from Mahaffy et al 2015, Gale crater sediment)
    d = 3
    
    ##from inputs 
    #(converting to grams of H2O):
    remaining_H2O = GEL(X) #X meter global layer equivalent of water we see today
    initial_H2O = GEL(Iw)
    #calculating the initial amount of HDO given intial inputs:
    initial_HDO = initial_H2O*2*(d*DH_before) #D/H ratio is 3*SMOW in the Hesperian so d = 3 in input
    #calculating the constant rate of outgassing from input:
    H2O_supply =  GEL(V)/tf  #H2O grams supplied by outgassing per year 
    #calculating the corresponding constant rate of outgassing of HDO
    HDO_supply = H2O_supply*2.0*DH_before #HDO grams supplied by outgassing (same as H2O but with terrestrial D/H ratio)
    #calculating the constant rate of escape of H2O needed to satisfy the input conditions
    H2O_loss = ((initial_H2O-remaining_H2O)/tf)+H2O_supply   #H2O grams lost to space per year 
                                                                #(amnt needed to balance out 10 m remaining global layer)

    #creating empty arrays to store variables in the loop
    x = np.zeros([n+1]) #space reservoir
    y = np.zeros_like(x) #h2o in exch reservoir
    z = np.zeros_like(x) #hdo in exch reservoir
    DH = np.zeros_like(x) #d/h in exch reservoir

    #create the time array
    t = np.linspace(0, tf, n+1)
    
    #assigning the initial conditions to the first entry in the arrays 
    x[0] = 0.0 #initial amount in escape box
    y[0] = initial_H2O #initial amount of water in exchangeable reservoir at t=0
    z[0] = initial_HDO #initial amount of deuterated water in exchangeable reservoir
    DH[0] = d*DH_before #initial D/H ratio of the reservoir
    
    #The Box Model loop
    for i in range(1,n+1):  #start at 1 since the 0 entries were the initial conditions
        Escaped_box = x[i-1]+ (H2O_loss*dt)     #how much water has escaped to space in each timestep 
                                                #[water in previous timestep + loss*time increment]
        Atmosphere_box = y[i-1] + (dt*H2O_supply) - (dt*H2O_loss) #how much water is in the exchangeable reservoir 
                                                                  #[water in previous timestep + supply*timeinc - loss*timeinc]
        HDO_loss = (R*z[i-1]*H2O_loss)/(y[i-1]) #how much HDO is lost in each time step 
                                                #(refer to notes for derivation of this equation, related to the fractionating nature of escape)
        HDO = z[i-1] + (dt*(HDO_supply - HDO_loss)) #how much HDO is left in the exchangeable reservoir
                                                    #[HDO in previous timestep + supply*timeinc - loss*timeinc]
        DH_change = HDO/(2.0*Atmosphere_box) #DH ratio in each time step [HDO/2*H2O in exch reservoir]
        #saving these values to their corresponding arrays
        x[i] = Escaped_box
        y[i] = Atmosphere_box
        z[i] = HDO
        DH[i] = DH_change
        #end of box model
    
    #Creating the pandas dataframe of the values output by the boxmodel:
    df1=[]
    for i in range(0,n+1):
        df1.append([x[i],y[i],z[i],DH[i], t[i]])
    DF = pd.DataFrame(df1)
    DF.columns = ['Space','Exch_Res','HDO','DH','t']
    
    #calculating the enrichment at the last timestep (today):
    enrichment = DH[n]/DH_before
    
    return DF, enrichment, H2O_loss

# Same as box_model function but doesn't create/save a dataframe and has less total steps, meant for parameter space study since DF isn't needed
def loop_box_model(Iw,X,V,R):
    '''The box model 
    Takes in inputs for boundary conditions (initial water in the past, water today, 
    outgassing rate and fractionation factor) and calculates an H2O escape rate to satisfy 
    the boundary conditions. Calculates the time dependant HDO escape rate and calculates 
    the amount of water in each box at each timestep. Outputs a dataframe with all the values 
    in each timestep as well as the final enrichment today and the calculated escape rate of water. 
    Inputs:
    Iw: The initial reservoir of water in meters GEL
    X: Water budget today in meters GEL
    V: total amount of water outgassed in the past tf (3 billion) years
    R: the fractionation factor of D escape compared to H escape 
    	(efficiency of D escape relative to H which ranges from 1e-5 to 0.32)
    Outputs:
    DF: Dataframe that gives the following information at each timestep:
        'Space'   :   amount of H2O in space that has escaped
        'Exch_Res':   amount of H2O remaining in the exchangeable reservoir
        'HDO'     :   amount of HDO remaining in the exchangeable reservoir
        'DH'      :   the D/H ratio in the exchangeable reservoir
        't'       :   time (in years) elapsed since the beginning of the model run (from tf (3 bn) yrs ago)
    enrichment: The final calculated enrichment at the end of the box model run ('today')
    H2O_loss: The average rate of escape of H2O (in grams per year)
    '''
    
    #start with defining the starting time in Ga
    tf = 3.3e9 #total time in years (we start at 3bn yrs ago)
    #amount of timesteps the model has to run through (the more the better but the more computation time needed)
    n = 10000 #total steps (usually 100000 but can lower to make parameter space study faster)
    #calculatin the time increment
    dt = tf/n
    
    #initial conditions
    
    #The terrestrial DH ratio we will be measuring enrichment against
    DH_before = 1.56e-4                   #terrestial value of D/H 
    
    #Initial enrichment (as a multiple of the terrestial value) tf (3 billion) years ago 
    	#(acquired from Mahaffy et al 2015, Gale crater sediment)
    d = 3
    
    ##from inputs 
    #(converting to grams of H2O):
    remaining_H2O = GEL(X) #X meter global layer equivalent of water we see today
    initial_H2O = GEL(Iw)
    #calculating the initial amount of HDO given intial inputs:
    initial_HDO = initial_H2O*2*(d*DH_before) #D/H ratio is 3*SMOW in the Hesperian so d = 3 in input
    #calculating the constant rate of outgassing from input:
    H2O_supply =  GEL(V)/tf  #H2O grams supplied by outgassing per year 
    #calculating the corresponding constant rate of outgassing of HDO
    HDO_supply = H2O_supply*2.0*DH_before #HDO grams supplied by outgassing (same as H2O but with terrestrial D/H ratio)
    #calculating the constant rate of escape of H2O needed to satisfy the input conditions
    H2O_loss = ((initial_H2O-remaining_H2O)/tf)+H2O_supply   #H2O grams lost to space per year 
                                                                #(amnt needed to balance out 10 m remaining global layer)

    #creating empty arrays to store variables in the loop
    x = np.zeros([n+1]) #space reservoir
    y = np.zeros_like(x) #h2o in exch reservoir
    z = np.zeros_like(x) #hdo in exch reservoir
    DH = np.zeros_like(x) #d/h in exch reservoir

    #create the time array
    t = np.linspace(0, tf, n+1)
    
    #assigning the initial conditions to the first entry in the arrays 
    x[0] = 0.0 #initial amount in escape box
    y[0] = initial_H2O #initial amount of water in exchangeable reservoir at t=0
    z[0] = initial_HDO #initial amount of deuterated water in exchangeable reservoir
    DH[0] = d*DH_before #initial D/H ratio of the reservoir
    
    #The Box Model loop
    for i in range(1,n+1):  #start at 1 since the 0 entries were the initial conditions
        Escaped_box = x[i-1]+ (H2O_loss*dt)     #how much water has escaped to space in each timestep 
                                                #[water in previous timestep + loss*time increment]
        Atmosphere_box = y[i-1] + (dt*H2O_supply) - (dt*H2O_loss) #how much water is in the exchangeable reservoir 
                                                                  #[water in previous timestep + supply*timeinc - loss*timeinc]
        HDO_loss = (R*z[i-1]*H2O_loss)/(y[i-1]) #how much HDO is lost in each time step 
                                                #(refer to notes for derivation of this equation, related to the fractionating nature of escape)
        HDO = z[i-1] + (dt*(HDO_supply - HDO_loss)) #how much HDO is left in the exchangeable reservoir
                                                    #[HDO in previous timestep + supply*timeinc - loss*timeinc]
        DH_change = HDO/(2.0*Atmosphere_box) #DH ratio in each time step [HDO/2*H2O in exch reservoir]
        #saving these values to their corresponding arrays
        x[i] = Escaped_box
        y[i] = Atmosphere_box
        z[i] = HDO
        DH[i] = DH_change
        #end of box model
    
    #calculating the enrichment at the last timestep (today):
    enrichment = DH[n]/DH_before
    
    return enrichment, H2O_loss

#NON-CONSTANT O|UTGASSING AND LOSS RATES:
#decreasing outgassing rate following cooling curve of planet (Jakosky & Shock 1998)
#decreasing loss rate such that final loss rate today matches MAVEN + total loss balances out initial and final water + outgassed

def var_outg(t,tf,V):
    ''' Function that creates an array of outgassing rates that are monotonically decreasing with time, 
    following the cooling curve of the planet
    Assumption: outgassing rate reaches 0 in 2bn years from present for Mars
    input: time array [t] at some discretization, total time [tf], total outgassed water in m GEL [V]
    output: array of outgassingg rates'''
    m = ((GEL(V))/((tf*5e9)-(((tf)**2)/2)))
    out = -m*t+(5e9*m)
    return out

def var_loss(t,tf,Iw,X,V):
	''' Function that creates an array of escape rates that are monotonically decreasing with time, 
    such that the escape rate today matches the MAVEN escape rates for Hydrogen.
    Assumption: escape rate MAVEN measures is representative of the current epoch
    input: time array [t] at some discretization, total time [tf], Initial water budget m GEL [Iw]
    , Total water budget today m GEL [X], total outgassed water in m GEL [V]
    output: array of escape rates'''
    mvn_loss = 8700*3.154e+7
    total_loss = GEL(Iw)+GEL(V)-GEL(X)
    m = (total_loss-(mvn_loss*tf))/(((tf)**2)-(((tf)**2)/2))
    loss_rate = -m*t + (mvn_loss + (m*tf))
    return loss_rate

#variable outgassing box model, here the outgassing is decreasing in time following the cooling of the planet. 
def var_box_model(Iw,X,V,R):
    '''The box model 
    Takes in inputs for boundary conditions (initial water in the past, water today, 
    outgassing rate and fractionation factor) and calculates an H2O escape rate to satisfy 
    the boundary conditions. Calculates the time dependant HDO escape rate and calculates 
    the amount of water in each box at each timestep. Outputs a dataframe with all the values 
    in each timestep as well as the final enrichment today and the calculated escape rate of water. 
    Inputs:
    Iw: The initial reservoir of water in meters GEL
    X: Water budget today in meters GEL
    V: total amount of water outgassed in the past tf (3.3 billion) years
    R: the fractionation factor of D escape compared to H escape 
    	(efficiency of D escape relative to H which ranges from 1e-5 to 0.32)
    Outputs:
    DF: Dataframe that gives the following information at each timestep:
        'Space'   :   amount of H2O in space that has escaped
        'Exch_Res':   amount of H2O remaining in the exchangeable reservoir
        'HDO'     :   amount of HDO remaining in the exchangeable reservoir
        'DH'      :   the D/H ratio in the exchangeable reservoir
        't'       :   time (in years) elapsed since the beginning of the model run (from tf (3 bn) yrs ago)
    enrichment: The final calculated enrichment at the end of the box model run ('today')
    H2O_loss: The average rate of escape of H2O (in grams per year)
    '''
    
    #start with defining the starting time in Ga
    tf = 3.3e9 #total time in years (we start at 3bn yrs ago)
    #amount of timesteps the model has to run through (the more the better but the more computation time needed)
    n = 100000 #total steps 
    #calculatin the time increment
    dt = tf/n
    
    #create the time array
    t = np.linspace(0, tf, n+1)

    #initial conditions
    
    #The terrestrial DH ratio we will be measuring enrichment against
    DH_before = 1.56e-4                   #terrestial value of D/H 
    
    #Initial enrichment (as a multiple of the terrestial value) tf (3 billion) years ago 
    	#(acquired from Mahaffy et al 2015, Gale crater sediment)
    d = 3
    
    ##from inputs 
    #(converting to grams of H2O):
    remaining_H2O = GEL(X) #X meter global layer equivalent of water we see today
    initial_H2O = GEL(Iw)
    #calculating the initial amount of HDO given intial inputs:
    initial_HDO = initial_H2O*2*(d*DH_before) #D/H ratio is 3*SMOW in the Hesperian so d = 3 in input
    #calculating the constant rate of outgassing from input:
    H2O_supply =  var_outg(t,tf,V)  #H2O grams supplied by outgassing per year 
    #calculating the corresponding constant rate of outgassing of HDO
    HDO_supply = H2O_supply*2.0*DH_before #HDO grams supplied by outgassing (same as H2O but with terrestrial D/H ratio)
    #calculating the constant rate of escape of H2O needed to satisfy the input conditions
    H2O_loss = var_loss(t,tf,Iw,X,V)  #H2O grams lost to space per year 
                                                                #(amnt needed to balance out 10 m remaining global layer)

    #creating empty arrays to store variables in the loop
    x = np.zeros([n+1]) #space reservoir
    y = np.zeros_like(x) #h2o in exch reservoir
    z = np.zeros_like(x) #hdo in exch reservoir
    DH = np.zeros_like(x) #d/h in exch reservoir

    #create the time array
    t = np.linspace(0, tf, n+1)
    
    #assigning the initial conditions to the first entry in the arrays 
    x[0] = 0.0 #initial amount in escape box
    y[0] = initial_H2O #initial amount of water in exchangeable reservoir at t=0
    z[0] = initial_HDO #initial amount of deuterated water in exchangeable reservoir
    DH[0] = d*DH_before #initial D/H ratio of the reservoir
    
    #The Box Model loop
    for i in range(1,n+1):  #start at 1 since the 0 entries were the initial conditions
        Escaped_box = x[i-1]+ (H2O_loss[i-1]*dt)     #how much water has escaped to space in each timestep 
                                                #[water in previous timestep + loss*time increment]
        Atmosphere_box = y[i-1] + (dt*H2O_supply[i-1]) - (dt*H2O_loss[i-1]) #how much water is in the exchangeable reservoir 
                                                                  #[water in previous timestep + supply*timeinc - loss*timeinc]
        HDO_loss = (R*z[i-1]*H2O_loss[i-1])/(y[i-1]) #how much HDO is lost in each time step 
                                                #(refer to notes for derivation of this equation, related to the fractionating nature of escape)
        HDO = z[i-1] + (dt*(HDO_supply[i-1] - HDO_loss)) #how much HDO is left in the exchangeable reservoir
                                                    #[HDO in previous timestep + supply*timeinc - loss*timeinc]
        DH_change = HDO/(2.0*Atmosphere_box) #DH ratio in each time step [HDO/2*H2O in exch reservoir]
        #saving these values to their corresponding arrays
        x[i] = Escaped_box
        y[i] = Atmosphere_box
        z[i] = HDO
        DH[i] = DH_change
        #end of box model
    
    #Creating the pandas dataframe of the values output by the boxmodel:
    df1=[]
    for i in range(0,n+1):
        df1.append([x[i],y[i],z[i],DH[i], t[i]])
    DF = pd.DataFrame(df1)
    DF.columns = ['Space','Exch_Res','HDO','DH','t']
    
    #calculating the enrichment at the last timestep (today):
    enrichment = DH[n]/DH_before
    
    return DF, enrichment, H2O_loss

#variable outgassing box model, here the outgassing is decreasing in time following the cooling of the planet. 
def loop_var_box_model(Iw,X,V,R):
    '''The box model 
    Takes in inputs for boundary conditions (initial water in the past, water today, 
    outgassing rate and fractionation factor) and calculates an H2O escape rate to satisfy 
    the boundary conditions. Calculates the time dependant HDO escape rate and calculates 
    the amount of water in each box at each timestep. Outputs a dataframe with all the values 
    in each timestep as well as the final enrichment today and the calculated escape rate of water. 
    Inputs:
    Iw: The initial reservoir of water in meters GEL
    X: Water budget today in meters GEL
    V: total amount of water outgassed in the past tf (3.3 billion) years
    R: the fractionation factor of D escape compared to H escape 
    	(efficiency of D escape relative to H which ranges from 1e-5 to 0.32)
    Outputs:
    DF: Dataframe that gives the following information at each timestep:
        'Space'   :   amount of H2O in space that has escaped
        'Exch_Res':   amount of H2O remaining in the exchangeable reservoir
        'HDO'     :   amount of HDO remaining in the exchangeable reservoir
        'DH'      :   the D/H ratio in the exchangeable reservoir
        't'       :   time (in years) elapsed since the beginning of the model run (from tf (3 bn) yrs ago)
    enrichment: The final calculated enrichment at the end of the box model run ('today')
    H2O_loss: The average rate of escape of H2O (in grams per year)
    '''
    
    #start with defining the starting time in Ga
    tf = 3.3e9 #total time in years (we start at 3bn yrs ago)
    #amount of timesteps the model has to run through (the more the better but the more computation time needed)
    n = 10000 #total steps 
    #calculatin the time increment
    dt = tf/n
    
    #create the time array
    t = np.linspace(0, tf, n+1)

    #initial conditions
    
    #The terrestrial DH ratio we will be measuring enrichment against
    DH_before = 1.56e-4                   #terrestial value of D/H 
    
    #Initial enrichment (as a multiple of the terrestial value) tf (3 billion) years ago 
    	#(acquired from Mahaffy et al 2015, Gale crater sediment)
    d = 3
    
    ##from inputs 
    #(converting to grams of H2O):
    remaining_H2O = GEL(X) #X meter global layer equivalent of water we see today
    initial_H2O = GEL(Iw)
    #calculating the initial amount of HDO given intial inputs:
    initial_HDO = initial_H2O*2*(d*DH_before) #D/H ratio is 3*SMOW in the Hesperian so d = 3 in input
    #calculating the constant rate of outgassing from input:
    H2O_supply =  var_outg(t,tf,V)  #H2O grams supplied by outgassing per year 
    #calculating the corresponding constant rate of outgassing of HDO
    HDO_supply = H2O_supply*2.0*DH_before #HDO grams supplied by outgassing (same as H2O but with terrestrial D/H ratio)
    #calculating the constant rate of escape of H2O needed to satisfy the input conditions
    H2O_loss = var_loss(t,tf,Iw,X,V)  #H2O grams lost to space per year 
                                                                #(amnt needed to balance out 10 m remaining global layer)

    #creating empty arrays to store variables in the loop
    x = np.zeros([n+1]) #space reservoir
    y = np.zeros_like(x) #h2o in exch reservoir
    z = np.zeros_like(x) #hdo in exch reservoir
    DH = np.zeros_like(x) #d/h in exch reservoir

    #create the time array
    t = np.linspace(0, tf, n+1)
    
    #assigning the initial conditions to the first entry in the arrays 
    x[0] = 0.0 #initial amount in escape box
    y[0] = initial_H2O #initial amount of water in exchangeable reservoir at t=0
    z[0] = initial_HDO #initial amount of deuterated water in exchangeable reservoir
    DH[0] = d*DH_before #initial D/H ratio of the reservoir
    
    #The Box Model loop
    for i in range(1,n+1):  #start at 1 since the 0 entries were the initial conditions
        Escaped_box = x[i-1]+ (H2O_loss[i-1]*dt)     #how much water has escaped to space in each timestep 
                                                #[water in previous timestep + loss*time increment]
        Atmosphere_box = y[i-1] + (dt*H2O_supply[i-1]) - (dt*H2O_loss[i-1]) #how much water is in the exchangeable reservoir 
                                                                  #[water in previous timestep + supply*timeinc - loss*timeinc]
        HDO_loss = (R*z[i-1]*H2O_loss[i-1])/(y[i-1]) #how much HDO is lost in each time step 
                                                #(refer to notes for derivation of this equation, related to the fractionating nature of escape)
        HDO = z[i-1] + (dt*(HDO_supply[i-1] - HDO_loss)) #how much HDO is left in the exchangeable reservoir
                                                    #[HDO in previous timestep + supply*timeinc - loss*timeinc]
        DH_change = HDO/(2.0*Atmosphere_box) #DH ratio in each time step [HDO/2*H2O in exch reservoir]
        #saving these values to their corresponding arrays
        x[i] = Escaped_box
        y[i] = Atmosphere_box
        z[i] = HDO
        DH[i] = DH_change
        #end of box model
    
    #calculating the enrichment at the last timestep (today):
    enrichment = DH[n]/DH_before
    
    return enrichment, H2O_loss

#PARAMETER SPACE STUDY SECTION

def parameter_study(a,b,c,d,e,f,R, filename):
    '''
    
    Creates an array of all possible combinations for the ranges 
    provided in the inputs for the three different variables and runs the model for each combination
    Inputs:
    a: start value for Initial Water
    b: end value for Initial Water
    c: start value for Remaining Water
    d: end value for Remaining Water
    e: start value for Outgassing
    f: end value for Outgassing
    R: fractionation factor
    filename: name of file to save dataframe of results to (must be in the form 'filename.csv')
    Output: 
    
    Output:
    A dataframe with all the succesful runs of the model,
    dataframe contains columns of the Initial amount, Remaining amount, Outgassed amount, DH ratio and Escape values listed
    '''
    
    #initial H2O total number of points
    n1 = b-a
    #remaining H2O total number of points
    n2 = d-c
    #outgassed H2O total number of points
    n3 = f-e

    # creating an array for each variable individually spaced by one unit m GEL
    init =  np.linspace(a, b, n1+1)
    rem  =  np.linspace(c, d, n2+1)
    outg =  np.linspace(e, f, n3+1)

    #Combining the three arrays to make one large array that contains all possible 
    #combinations of the three variable in the provided ranges
    vals = []
    for i in range(0,n1+1):
        for j in range(0,n2+1):
            for k in range(0,n3+1):
                vals.append([init[i], rem[j], outg[k],R])
                
    
    #create an empty array to store the results of the parameter study in            
    result = []

    #looping through the parameter space as listed in vals and running the model for each entry in vals
    for i in tnrange(len(vals)):  #tnrange shows a live progress-bar in the notebook
        enrichment, a = loop_box_model(*vals[i])
        # checking result of model against success criteria
        # if enrichment matches current atmosphere, append entry to result.. otherwise discard 
        if 5 <= enrichment <= 6:
            result.append([vals[i][0], vals[i][1], vals[i][2], enrichment, a])

    #convert to pandas dataframe
    Dataframe_result= pd.DataFrame(result)
    Dataframe_result.columns = ['Initial','Remainder','Outgassed','DH','Escape']

    #save dataframe to csv file externally
    Dataframe_result.to_csv(filename)
    return Dataframe_result

def parameter_study_var(a,b,c,d,e,f,R, filename):
    '''
    
    Creates an array of all possible combinations for the ranges 
    provided in the inputs for the three different variables and runs the model for each combination
    Inputs:
    a: start value for Initial Water
    b: end value for Initial Water
    c: start value for Remaining Water
    d: end value for Remaining Water
    e: start value for Outgassing
    f: end value for Outgassing
    R: fractionation factor
    filename: name of file to save dataframe of results to (must be in the form 'filename.csv')
    Output: 
    
    Output:
    A dataframe with all the succesful runs of the model,
    dataframe contains columns of the Initial amount, Remaining amount, Outgassed amount, DH ratio and Escape values listed
    '''
    
    #initial H2O total number of points
    n1 = b-a
    #remaining H2O total number of points
    n2 = d-c
    #outgassed H2O total number of points
    n3 = f-e

    # creating an array for each variable individually spaced by one unit m GEL
    init =  np.linspace(a, b, n1+1)
    rem  =  np.linspace(c, d, n2+1)
    outg =  np.linspace(e, f, n3+1)

    #Combining the three arrays to make one large array that contains all possible 
    #combinations of the three variable in the provided ranges
    vals = []
    for i in range(0,n1+1):
        for j in range(0,n2+1):
            for k in range(0,n3+1):
                vals.append([init[i], rem[j], outg[k],R])
                
    
    #create an empty array to store the results of the parameter study in            
    result = []

    #looping through the parameter space as listed in vals and running the model for each entry in vals
    for i in tnrange(len(vals)):  #tnrange shows a live progress-bar in the notebook
        enrichment, a = loop_var_box_model(*vals[i])
        # checking result of model against success criteria
        # if enrichment matches current atmosphere, append entry to result.. otherwise discard 
        if 5 <= enrichment <= 6:
            result.append([vals[i][0], vals[i][1], vals[i][2], enrichment, a])

    #convert to pandas dataframe
    Dataframe_result= pd.DataFrame(result)
    Dataframe_result.columns = ['Initial','Remainder','Outgassed','DH','Escape']

    #save dataframe to csv file externally
    Dataframe_result.to_csv(filename)
    return Dataframe_result

def find_init(a,b,R,X,V,e):
    '''Loops through different values for initial water and outputs the amount of initial water needed for given enrichment
    given constant outgassing and remining water
    Input: 
    a: beginning value for range of initial water
    b: end value for range of initial water
    X: Remaining (current) water budget m GEL
    V: total amount of water outgassed in the past tf  years
    e: current day enrichment we want to achieve
    Output:
    init[locmin]: initial water needed to recreate enrichment
    DH_init[locmin]: the corresponding enrichment (should be equal or close to the input e)
    DH_loss[locmin]: the corresponding water loss'''
    nn = b-a
    init = np.linspace(a, b, nn+1)
    DH_init = np.zeros_like(init)
    DH_loss = np.zeros_like(init)
    check = True
    mult = 1
    while check == True:
        for i in range(nn+1):
            DH_init[i] = lvbm(init[i],X,V,R)[0]
            DH_loss[i] = lvbm(init[i],X,V,R)[1][0]
        if abs(DH_init - e).min() >= 0.001:
            locmin =  abs(DH_init - e).argmin()
            mult = 0.1*mult
            init = np.linspace(init[locmin]+(init[locmin]*mult), init[locmin]-(init[locmin]*mult), nn+1)
            check == True
        else:
            locmin =  abs(DH_init - e).argmin()
            check == False
            break
    return init[locmin], DH_init[locmin], DH_loss[locmin]


