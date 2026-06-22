Team and Contributors

D-FENSE

Americo Cunha Jr (team leader) (LNCC - UERJ)

Paulo Esquef (LNCC)

Emanuelle Paixão (LNCC)


2. Repository Structure

Main
|_>  Aggredated Data

     (Data used to train the model)
     
|_>  LNCC_SurgeModel26_Dengue  (model source and results)

     |_> validation1    (tranining for validation 1) 
     
         |_> matlab    (matlab scritps to generate the results)
         
         |_> spreadsheets  (CSV files with the results; one for each state)
         
         |_> plots   (PDF files with plots of the results; one for each state)
         
     |_> validation2  (see validation 1)
     
         |_> matlab
         
         |_> spreadsheets
         
         |_> plots
         
     |_> validation3 (see validation 1)
     
         |_> matlab
         
         |_> spreadsheets
         
         |_> plots
         
     |_> validation4  (see validation 1)
     
         |_> matlab
         
         |_> spreadsheets
         
         |_> plots
         
         
3. Libraries and Dependencies

Proprietary functions of Matlab toolboxes:
table.m  % matlab function for managing table creation
readtable.m % matlab function for reading a CSV file 
writetable.m % matlab function for creating a CSV file from a table
buffer.m  % matlab function for data organization into columns (block segmentation)
circshift.m % matlab function used to align the incidence curves by their maxima, via circular shifts. 
lsqcurvefit.m % matlab function for curve fitting: finds parameters of an incidence model to an average incidence curve obtained from the smoothed time-series of dengue cases.  
lognfit.m  % matlab function used to estimate the parameters of a log normal distribution from a set of data
lognrnd.m % matlab function used to generate values from a log normal distribution
prctile.m % matlab function for calculating the percentiles of a set of numbers

Extra functions:
ssa_modPE.m  % Singular Spectrum Analysis method used for lowpass filtering the data (included in folder matlab)
alignDC.m  % function that aligns the smoothed incidence (surge) curves (52 epiweeks) so that their maxima are aligned to the middle of the window (sample 26).
 

4. Data and Variables
For each state, the input time series aggregates the sum of all dengue cases, over a epiweek.
Data is organized in 27 CSV files, one per state, stored in folder 'Aggregated Data'.
Only the time series of dengue cases has been used. No climate data has been used.
  
5. Model Training

Data processing: an average incidence curve is obtained by first lowpass filtering the time series of dengue cases via the Singular Spectrum Analysis method. Then, respecting the cutoff date for each validation (1 to 4), the time series is organized in blocks that span from epiweek (EW) 41 of a given year to EW 40 of the subsequent year, a total of 52 EWs. Then, the set of incidence curves are aligned in time so that the maxima of each curve gets centered at sample 26 (middle of the block). This is carried out by ancilary function alignDC.m. Finally, the arithmetic mean of the aligned incidence curves is calculated to produce a typical incidence curve (data driven). 

Model training:
From the typical incidence curve, the parameters of an incidence model are estimated via function lsqcurvefit.m. From the model we generate a modeled incidence curve. 
Then, we obtain a gain (g), so that g times the modeled incidence curve matches one of the observed incidence curve (smoothed and aligned). We do that for all observed incidence curves. From the set of values of g we estimate the parameters of a log normal distribution. From that distribution, we generate 10000 values of g. The set of 10000 forecasts for the next season is produced by g * modeled incidence curve. From this set, for a given time instant, we calculate the necesary percentiles (via function prctile.m) to produce the median prediction and the specified prediction intervals.


CODE (example for validation 1) : 
dcases=dcases(1:ind_EW25_2022);  % crops time-series at EW 25 of 2022 (including)

% lowpass filtering
L=26; % window length for the SSA filter
nsv=6; % number of selected eigenvalues (ordered)
[dcasesf]=ssa_modPE(dcases,L,nsv); % filtered time series 

dcasesfs=dcasesf(41:ind_EW25_2022);  % selects filtered data from EW 41/2010 to EW 25/2022

SS=52;  % seasonality of 52 Epidemic Weeks (EW)

DC=buffer(dcasesfs,SS); % organizes vector dcases in a matrix with 52 rows
% each column of DC has cases for 52 EWs. 

[DCalign, ind_max]=alignDC(DC); % aligns surges so that they all peak
% in the middle of the window between EW 41 year to EW 40 subsequent year

typ_DC=mean(DCalign')'; % typical surge waveform in 52 EWs:
% from EW 41 to EW 40 (subsequent year)  

% Surge Model Estimation (from typ_DC)

% Initial parameters of the logistic model
L=120000;  k=0.3; n0=26;  % time-shift
P=[L,k,n0];  % initial parameter vector
n=0:length(typ_DC)-1; n=n(:);  % time basis

% Surge Model Estimation (nonlinear) 
options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective');
fun = @(P,n) (P(1).*P(2).*exp(P(2).*(n-(P(3)))))./(1+exp(P(2).*(n-(P(3)))).^2); % surge model
P_est = lsqcurvefit(fun,P,n,typ_DC,[500;0.1;24],[370000;0.5;28],options); % find estimated model parameters P_est
P_est(3)=round(P_est(3)); % rounds estimated n0, since it is supposed to be an integer 

Model_surge=fun(P_est,n);  % Synthesizes a modeled incidence curve from estimated model 

% Now, we calculate a set of gains g so tha an observed surge is given by g * Model_surge

ns=size(DCalign,2); % number of observed surges (number of columns of DCalign)

g=zeros(ns,1);  % initializes with zeros vector to store the set of gains
x=Model_surge;   % surge template: attributes Model_surge to x (for notation clarity)  
for kk=1:ns
    a=DCalign(:,kk);   % observed surge at column kk
    g(kk)=a'*x./(x'*x);   % calculates amplitude gain for each observed surge at column kk
end

[param]=lognfit(g);  % estimates log nomal distribution from g

mg=param(1);  % mean of value of g (from 2010 to 2022)
sigma=param(2); % standard deviation of value of g (from 2010 to 2022)

MC=10000;  % number of montecarlo runs for forecast surges in 2023

g_MC=lognrnd(mg,sigma,MC,1);  % set of randomly generated gains

for kk=1:MC
    forecast_cases_v1(:,kk)=g_MC(kk)*Model_surge; % set of MC realizations of the surge forecast
end

set_prctile=[2.5 5 10 25 50 75 90 95 97.5]; % 2.5 to 97.5% percentiles
PP=prctile(forecast_cases_v1',set_prctile); % calculates the percentiles and stores in PP


5. Data Usage Restriction

We respected the cutoff date of EW 25 year for each validation. In the temporal organization of the data into blocks of 52 EWs (from EW41 of year to EW40 subsequent year), via command DC=buffer(dcasesfs,52), the last block ends at EW 25 and dengue cases from EW 26 to EW 40 are set to zero. This is not a problem for the subsequent time alignment of the surges in each block because the peak of the last the surge happens at about EW 16, i.e., before EW 25. 


6. Predictive Uncertainty How are your prediction intervals computed?
We have carried out a Monte Carlo simulation with 10000 runs, for each state. Out of the 10000 forecast realizations, for a given time instant, we calculate the specified percentiles (via prctile.m), to produce the specified prediction intervals.    


8. References
  None. 
