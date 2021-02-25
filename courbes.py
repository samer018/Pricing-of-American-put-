# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 18:04:46 2019

@author: 33695
"""
import numpy as np
import matplotlib.pyplot as plt
spot_price = 100

# Long call
strike_price_long_call = 90
premium_long_call = 0

# Stock price range at expiration of the put
sT = np.arange(70,120,1)
def call_payoff(sT, strike_price, premium):
    return np.where(sT > strike_price, sT - strike_price, 0) - premium

payoff_long_call = call_payoff (sT, strike_price_long_call, premium_long_call)
# Plot
fig, ax = plt.subplots()
ax.spines['top'].set_visible(False) # Top border removed
ax.spines['right'].set_visible(False) # Right border removed
ax.spines['bottom'].set_position('zero') # Sets the X-axis in the center
ax.plot(sT,payoff_long_call,label='Long Call',color='g')
plt.xlabel('Stock Price')
plt.ylabel('Profit and loss')
plt.legend()
#plt.show()
stock_payoff = (sT - spot_price)*-1.0

fig, ax = plt.subplots()
ax.spines['top'].set_visible(False) # Top border removed
ax.spines['right'].set_visible(False) # Right border removed
ax.spines['bottom'].set_position('zero') # Sets the X-axis in the center
ax.plot(sT,stock_payoff,label='Stock Payoff',color='b')
plt.xlabel('Stock Price')
plt.ylabel('Profit and loss')
plt.legend()
#plt.show()
payoff_synthetic_long_put = payoff_long_call + stock_payoff
# Plot
fig, ax = plt.subplots()
ax.spines['top'].set_visible(False) # Top border removed
ax.spines['right'].set_visible(False) # Right border removed
ax.spines['bottom'].set_position('zero') # Sets the X-axis in the center
ax.plot(sT,payoff_synthetic_long_put,label=' Long Put')
plt.xlabel('Stock Price')
plt.ylabel('Profit and loss')
plt.legend()
plt.show()

from american_option_cedric import q_tau_estimate
