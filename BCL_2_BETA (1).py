#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Importing necessary libraries
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from chembl_webresource_client.new_client import new_client


# In[2]:


# Initializing new client
activity = new_client.activity


# In[3]:


# CHEMBL ID for Apoptosis Regulator Bcl-2 in Homo sapiens
selected_target = 'CHEMBL4860'


# In[4]:


# Geting activity data for the selected target with IC50 values
res = activity.filter(target_chembl_id=selected_target).filter(standard_type="IC50")



# In[5]:


# Converting activity data to DataFrame
df = pd.DataFrame.from_dict(res)


# In[6]:


# Printing the first few rows of the DataFrame
print(df.head(3))


# In[7]:


df.to_csv('Bcl_2_data_raw.csv', index=False)


# In[8]:


# Preprocessing the data
df2 = df[df.standard_value.notna()]
df2['standard_value'] = pd.to_numeric(df2['standard_value'], errors='coerce')
df2.dropna(subset=['standard_value'], inplace=True)
df2['pIC50'] = -np.log10(df2['standard_value'] * (10**-9))


# In[9]:


# Categorize the bioactivity
active_cutoff = np.percentile(df2['pIC50'], 75)
inactive_cutoff = np.percentile(df2['pIC50'], 25)


# In[10]:


def categorize_activity(value, active_threshold, inactive_threshold):
    if value >= active_threshold:
        return "active"
    elif value <= inactive_threshold:
        return "inactive"
    else:
        return "intermediate"

df2['bioactivity_class'] = df2['pIC50'].apply(categorize_activity, args=(active_cutoff, inactive_cutoff))


# In[11]:


# Visualize the distribution of pIC50 values
sns.histplot(df2['pIC50'], bins=30, kde=False)
plt.xlabel('pIC50')
plt.ylabel('Frequency')
plt.title('Distribution of pIC50 values for Bcl-2 Target')
plt.show()


# In[12]:


# Print the head of the updated DataFrame
print(df2[['molecule_chembl_id', 'canonical_smiles', 'standard_value', 'pIC50', 'bioactivity_class']].head())


# In[13]:


# Filter to select only compounds with high potency (pIC50 > 8)
high_potency_compounds = df2[df2['pIC50'] > 8]

# Print the first few rows of the DataFrame with high potency compounds
print(high_potency_compounds.head())

# Optionally, you can save this filtered data to a new CSV file
high_potency_compounds.to_csv('high_potency_compounds_bcl.csv', index=False)


# In[14]:


plt.figure(figsize=(10, 5))
sns.boxplot(df2['pIC50'])
plt.title('Boxplot of pIC50 Values')
plt.show()


# In[15]:


plt.figure(figsize=(10, 5))
sns.histplot(df2['pIC50'], kde=True)
plt.title('Histogram of pIC50 Values')
plt.show()


# In[16]:


high_potency_compounds = df2[df2['pIC50'] > 8]


# In[17]:


# Sort the high_potency_compounds DataFrame based on pIC50 values in descending order
sorted_high_potency_compounds = high_potency_compounds.sort_values(by='pIC50', ascending=False)




# In[18]:


# Select the top 5 compounds
top_8_compounds = sorted_high_potency_compounds.head(8)

# Print the names (CHEMBL IDs) of the top 5 compounds
print("Top 8 Potential Drugs for Bcl-2 based on pIC50:")
print(top_8_compounds[['molecule_chembl_id', 'pIC50']])


# In[ ]:





# In[ ]:




