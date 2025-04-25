### Define auxiliary functions for fitting
import scipy.optimize as opt
import numpy as np
def true_p(dv,a,b):
  #return 1/(1+np.exp((dv-a)/b))
  lo = log_odds(dv,a,b)
  l = log_p(lo)
  l[l>20]= 20 #Saturate to avoid numerical issues
  return np.exp(l)
def obs_p(dv,a,b,c):
  p = (1-c)*true_p(dv,a,b) + (c)*.5
  return p
def log_odds(dv,a,b,c=0): #log p/(1-p)
  if c == 0:
    return -(dv-a)/b
  else:
    p = obs_p(dv,a,b,c)
    return np.log(p/(1-p)) #Check if there is a computationally better approach
def log_p(log_odds):
  #compute log(p) from log-odds, i.e. log (p/(1-p)), avoiding numerical issues
  th = 20 #Threshold avoid numerical issues
  out = np.zeros(len(log_odds))
  out[log_odds > th] = 0 #Can we be more precise here? Do we need to be?
  out[log_odds < -th] = log_odds[log_odds < -th]
  out[np.abs(log_odds)<=th] = -np.log(1+np.exp(-log_odds[np.abs(log_odds)<=th]))
  return out
def exp_RT(t_nd, noise_std, a, b, dv):
  #p2 = true_p(dv,a,b)
  tol = 1e-8
  tol2 = 20
  #factor = (2*p2-1)/(np.log(1-p2) - np.log(p2))
  #p2[np.abs(dv-a)<tol]=.5
  #dv[np.abs(dv-a)<tol] = a+1e-5 #Just to avoid numerical issues
  #factor = (2*p2-1)*b/(dv-a)
  #factor[np.abs(p2-.5)<tol] = -.5 #Just to avoid numerical issues
  #factor[np.abs(p2-1)<tol] = -0.0000000000000001 #54
  #factor[np.abs(p2)<tol] = -0.00000000000000001 #54 #Should never happen, just in case
  lo = log_odds(dv,a,b)
  lp = log_p(lo)
  factor = np.zeros(len(lo)) #factor = (1-2*p2)/lo
  factor[np.abs(lo)<tol] = -.5
  factor[np.abs(lo)>tol2] = -1/np.abs(lo[np.abs(lo)>tol2]) #This will be roughly 0
  aux = true_p(dv[factor==0], a, b) #Only computing probabilities in a sensible range to avoid exp overflow
  factor[factor==0] = (1-2*aux)/lo[factor==0]
  t = t_nd - factor * (2/noise_std**2)
  return t
def neg_log_likelihood(params,v,dv,nreps=None):
  if nreps is None:
    nreps = np.ones(len(v))
  #nreps is the number of repetitions of each dv, used to compute the likelihood. When a different number of reps is used per dv, it is necessary to pass this parameter
  a,b,c = params
  #p2 = obs_p(dv,a,b,c)
  #return -np.sum(nreps*(v* np.log(p2) + (1-v)*np.log(1-p2)))
  lo = log_odds(dv,a,b,c)
  l = log_p(lo)
  return -np.sum(nreps*((v-1)*lo+l))
def fit_params(v, dv, bounds = ((-100, 100), (-200, -0.01), (0,.5)), params_init = [0, -100, .05], method = 'logl'):
  #Optimize setting a very low tolerance for stopping criteria
  if method == 'logl':
    result = opt.minimize(neg_log_likelihood, params_init, args=(v,dv), bounds=bounds, tol=1e-8)
  elif method == 'lsq':
    cost = lambda params: np.sum((obs_p(dv, *params) - v)**2)
    result = opt.minimize(cost, params_init, bounds=bounds, tol=1e-8)
  a,b,c = result.x
  if not result.success:
    print(result)
  return a,b,c
def time_RMS(params, dv, t):
  t_nd, noise_std, a, b = params
  err = np.sum((t - exp_RT(t_nd, noise_std, a, b, dv))**2)
  return err
def fit_RT(t,dv, bounds = ((-3,3), (0.001, 1), (-100, 100), (-200, -0.01))):
  #Set params_init to middle point of bounds:
  params_init = [np.mean(b) for b in bounds]
  #If any of the bounds was inf, set to the other bound (if finite), else, set to 0
  params_init[bounds[0]==-np.inf] = bounds[1][bounds[0]==-np.inf]
  params_init[bounds[1]==np.inf] = bounds[0][bounds[1]==np.inf]
  params_init[bounds[0]==-np.inf and bounds[1]==np.inf] = 0
  #result = opt.minimize(time_RMS, params_init, args=(dv, t),bounds=bounds, tol=1e-8)
  #Change optimizer to a least squares one
  #First define a lambda function as t-exp_RT(params,dv)
  err = lambda params, dv, t: t-exp_RT(*params, dv)
  #Reformat bounds as a tuple of two np arrays, one for lb and one for ub
  bounds = (np.array([b[0] for b in bounds]), np.array([b[1] for b in bounds]))
  result = opt.least_squares(err, params_init, args=(dv, t), bounds=bounds, xtol=1e-8, ftol=1e-8)
  t_nd, noise_std, a, b = result.x
  #If any optimal parameter is at the boundary or close enough (1%?), print a warning
  if np.any(np.abs(result.x/bounds[0] -1 )<1e-2) or np.any(np.abs(result.x/bounds[1]-1)<1e-2):
    print('Warning: at least one parameter is at the boundary of the optimization region')
    for i in range(len(result.x)):
      if result.x[i] == bounds[0][i] or result.x[i] == bounds[1][i]:
        print('Parameter', i, 'is at the boundary')
  if not result.success:
    print(result)
  return t_nd, noise_std, a, b
def neg_log_RTpdf(params, dv, t):
  #TODO: everything here was assuming that the distance between barriers was =1, but in my parametrization it is =2. I need to fix this
  t_nd, noise_std, a, b = params
  noise_std = noise_std/np.sqrt(2) #This is a hack to fix the fact that I was assuming that the distance between barriers was 1, but in my parametrization it is 2
  #p = true_p(dv,a,b) #=1/(1+np.exp((dv-a)/b)) = 1/(1+np.exp(-v/s^2))
  lo = log_odds(dv,a,b) #=np.log(p/(1-p))
  lp = log_p(lo) #=np.log(p)
  p = np.exp(lp) #I could just use lp and lo and avoid this step
  sum = 0
  #Sum over odd k from 1 to 10
  ts2 = (t-t_nd)*noise_std**2/2
  ts2[ts2<0] = 0 #This is a hack to avoid numerical issues when t<t_nd
  N = 100
  for k in range(1, N, 2):
    #sum += k*np.sin(k*np.pi/2) * np.exp(-(t-t_nd)*(k**2*np.pi**2+lo**2)*noise_std**2/2) 
    aux = 1 if k%4==1 else -1
    sum += k* aux * np.exp(-ts2*(k**2*np.pi**2))  #NOTE: the exponential can be split into two factors, and one of them can be taken out of the sum, so we only need to compute the sum over the other term and add an appropriate term to get the log likelihood
  logsum = np.log(sum) - ts2*lo**2
  #Now numerically fix ill-conditioned situations
  th = 6 / N #The convergence of the above series in N terms can only be assumed if t-t_nd is at least 6/N, though this probably depends on noise_std
  logsum[t<(t_nd+th)] = -1/(t[t<(t_nd+th)]-t_nd) #TODO: compute the limit of the log-sum when t->t_nd and use that expression
  logsum[t<=t_nd] = -1000 - 100*(t_nd-t) #This is a hack to avoid numerical issues when t<=t_nd, is roughly equivalent to not resolving timesteps below 1e-3 but maintianing differentiability
  #pdf = np.log(np.pi) + 2*np.log(noise_std) - np.log(p*(1-p))/2  +np.log(sum) #log(p*(1-p)) could be replaced by -log odds +2*log(p)
  pdf = np.log(np.pi) + 2*np.log(noise_std) + lo/2 - lp + logsum
  return -np.sum(pdf)
def fit_RT2(t, dv, bounds = ((-2,2), (0.001, 1), (-100, 100), (-200, -0.01)), params_init = [2, 0.2, 0, -100]):
  result = opt.minimize(neg_log_RTpdf, params_init, args=(dv,t), bounds=bounds, tol=1e-8)
  t_nd, noise_std, a, b = result.x
  if not result.success:
    print(result)
  return t_nd, noise_std, a, b

#Load data in detailed form
import pandas as pd
import numpy as np
def load_data():  
  # Load the data matrix as a Pandas DataFrame
  data = pd.read_csv('../data/Dataset.csv')

  ### Do some visualization sanity checks
  # Print the DataFrame
  print(data.head())

  # Group the DataFrame by subID and pertSize
  individual_data = data.groupby(['subID', 'pertSize'])
  grouped_data = data.groupby(['pertSize'])

  # Compute the average value of the leftResponse field for each group
  avg_left_response_by_subject = individual_data['leftResponse'].mean()
  avg_left_response_aggregated = grouped_data['leftResponse'].mean()
  test = grouped_data['leftResponse'].count()

  #Plot, in a single plot, every subject's average leftResponse values for each pertSize
  import matplotlib.pyplot as plt
  for subID in data['subID'].unique():
      plt.plot(avg_left_response_by_subject[subID], label='Subject ' + str(subID))
  #Now, add a new figure and plot the log-odds of the leftResponse values for each pertSize
  plt.figure()
  # for subID in data['subID'].unique():
  #     p = avg_left_response_by_subject[subID]
  #     #Get log-odds:
  #     lo = np.log(p/(1-p))
  #     #PLot, using markers, not lines:
  #     plt.plot(lo, label='Subject ' + str(subID), marker='o')

  # # Print the average leftResponse values
  # print(avg_left_response_by_subject)
  # print(avg_left_response_aggregated)
  # detailed_dv = data['pertSize'].unique()
  # print(data.groupby(['subID'])['leftResponse'].count())
  # print(data.groupby(['pertSize'])['noResponse'].mean())
  # print(data.groupby(['pertSize'])['leftResponse'].mean())
  # print(data.groupby(['pertSize'])['rightResponse'].mean())
  # print(data.groupby(['pertSize'])['prevSize'].mean())
  return data

    