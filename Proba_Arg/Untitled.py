#!/usr/bin/env python
# coding: utf-8

# In[34]:


import numpy
import matplotlib.pyplot as plt


# In[45]:


kmax = numpy.arange(2, 20, 0.05)


# In[46]:


kmin = numpy.arange(2, 20, 0.05)


# In[47]:


def complexity(kmin, kmax):
    #eturn 2**(2*kmin-1) - kmax*kmax
    return 2**(n*(kmin+1)) + kmin - 1 - kmax**(2*n)


# In[48]:


# plot the data
fig = plt.figure()

zero_points = [(km, kM) for km in kmin for kM in kmax if complexity(km, kM) >= 0 and km <= kM]
zero_x = [x for (x, _) in zero_points]
zero_y = [y for (_, y) in zero_points]

zero_points2 = [(km, kM) for km in kmin for kM in kmax if complexity(km, kM) < 0 and km <= kM]
zero_x2 = [x for (x, _) in zero_points2]
zero_y2 = [y for (_, y) in zero_points2]

plt.plot(zero_x, zero_y, color="green",label="Upper bound Fast_AG < Lower bound Constellation")
plt.plot(zero_x2, zero_y2, color="red",label="Upper bound Fast_AG > Lower bound Constellation")
plt.plot(kmin, kmax)
plt.legend()
plt.ylabel("kmax")
plt.xlabel("kmin")

#print(zero_points)

# display the plot
plt.show()


# In[ ]:




