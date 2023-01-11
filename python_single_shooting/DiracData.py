import pandas as pd
import scipy.constants as constants
import sys

data = pd.read_excel(r'C:\Users\alexa\Documents\Shooting_Dirac\python_single_shooting\data.xlsx')
data.index = data['Unnamed: 0']
data.pop('Unnamed: 0')

pardf= pd.read_excel(r'C:\Users\alexa\Documents\Shooting_Dirac\python_single_shooting\parameters.xlsx')

p1 = {pardf.columns[x]:pardf[pardf.Scenario == 1].values[0][x] for x in range(1,8)}
p2 = {pardf.columns[x]:pardf[pardf.Scenario == 2].values[0][x] for x in range(1,7)}
p3 = {pardf.columns[x]:pardf[pardf.Scenario == 3].values[0][x] for x in range(1,8)}

meq = "mass energy equivalent in MeV"
mn = constants.physical_constants["neutron " + meq][0] #energy in MeV
mp = constants.physical_constants["proton " + meq][0] #energy in MeV
