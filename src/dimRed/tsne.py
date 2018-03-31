import pandas as pd
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from numpy.random import rand 

df = rand(100,200)
x = df
model = TSNE(n_components=2, random_state=0)
model.fit_transform(x)

res = model.fit_transform(x)

pd_res = pd.DataFrame(res[0:, 0:,])

plt.scatter(pd_res[50:,0] , pd_res[50:,1])

