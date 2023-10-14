#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import api.chemistry as chem


# In[ ]:


chem.Molecule('NaCl').mass


# In[ ]:


chem.Molecule('NaCl').solubility


# In[ ]:


chem.Molecule('NaCl').bond


# In[ ]:


chem.Equation(['Na(OH)', 'HCl']).reaction


# In[ ]:


pred = chem.Equation(['NaCl', 'HCl']).predict()
pred.short


# In[ ]:


re, pr = chem.Equation(['Na(OH)', 'HCl']).balance(manual=['NaCl', 'H2O'])
print(f'{re.short} -> {pr.short}')


# In[ ]:


re, pr = chem.Equation(['Na(OH)', 'HCl']).balance()
print(f'{re.short} -> {pr.short}')

