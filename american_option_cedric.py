# -*- coding: utf-8 -*-
"""
Éditeur de Spyder 

Ceci est un script temporaire.
"""
import numpy as np
import keras
from keras.layers import Dense, Activation
from keras.models import Sequential
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
sc = StandardScaler()
#Process of the underlyings X_i_j_t
n = 2000
T = 48
K = 90
m = 12
r = 0.05
initial_stocks = 100*np.ones(m)
vol_matrix = np.matrix[[0.3024,0.1354,0.0722,0.1367,0.1647],[0.1354,0.2270,0.0613,0.1264,0.1610],[0.0722,0.0613,0.0717,0.0884,0.2937,0.1394],[0.1641,0.1610,0.0699,0.1394,0.2535]] 
def sample_X_t(t):
    X_t = np.zeros(m)
    for i in range(m):
        log_trend = 0
        for j in range(m):
            log_trend += vol_matrix[i][j]*np.sqrt(t)*np.random.normal(0,1) - 0.5*vol_matrix[i][j]**2*t
        X_t[i] = initial_stocks[i]*np.exp(r*t)*np.exp(log_trend)
    return X_t
def sample_X(T,n):
    X_l_t = np.zeros((T+1,T))
    for t in range(T+1):
        for l in range(T):
            X_l_t[t][l] = np.zeros(n)
            for i in range(n):
                X_l_t[t][l][i] = sample_X_t(t)
    return X_l_t
def pay_off_put_american(X_t,t,K):
    s = np.sum(X_t)
    return np.exp(-r*t)*max(K-s/len(X_t),0)
#estimation des qt et choix du k optimal pour déterminer tau étoile
def qtau_estimate():
    sample = sample_X(T,n)
    X_t_plus_1 = sample[t+1][t+1]
    f_t_plus_1 = pay_off_put_american(X_t_plus_1,t+1,K)
    q_t_plus_1 = q_t_estimates(t+1)
    y_t = np.maximum(f_t_plus_1,q_t_plus_1)
    X_t = sample[t][t]
    def q_t_estimates(t,k,Xt):
        if (t == T):
            return  0
        else:        
# Splitting the dataset into the Training set and Test set
            X_train, X_test, y_train, y_test = train_test_split(X_t, y_t, test_size = 0.5, random_state = 0)
            X_train = sc.fit_transform(X_train)
            X_test = sc.transform(X_test)
# Initialising the ANN
            model = Sequential()
# Adding the input layer and the first hidden layer
            model.add(Dense(32, activation = 'sigmoid', input_dim = 5))
# Adding the second hidden layer
            model.add(Dense(units = k, activation = 'sigmoid'))
# Adding the third hidden layer
            model.add(Dense(units = k, activation = 'sigmoid'))
# Adding the output layer
            model.add(Dense(units = 1))
        #model.add(Dense(1))
# Compiling the ANN
            model.compile(optimizer = 'adam', loss = 'mean_squared_error')
# Fitting the ANN to the Training set
            model.fit(X_train, y_train, batch_size = 10, epochs = 100)
            y_pred = model.predict(X_test)
        #plt.plot(y_test, color = 'red', label = 'Real data')
        #plt.plot(y_pred, color = 'blue', label = 'Predicted data')
        #plt.title('Prediction')
        #plt.legend()
        #plt.show()"
            return y_pred
    new_samples = np.zeros(n)
    for i in range(n):
        new_samples[i] = (X_t[i],q_t_estimates(t,k,Xt))
        
        

        
    
