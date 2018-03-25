import numpy as np
import matplotlib.pyplot as plt
from matplotlib import offsetbox
from sklearn import manifold, datasets, decomposition, discriminant_analysis
#import StringIO

'''
with open('../Single_cell_data/sample_data.txt') as f:
  rows,cols=np.fromfile(f, dtype=int, count=2, sep=" ")
  data = np.fromfile(f, dtype=int, count=cols*rows, sep=" ").reshape((rows,cols))
  #547,96
'''

'''
with open('../../Single_cell_data/sample_data.txt','r') as f:
    temp = f.readlines()
X = np.array(temp)
'''


#import numpy 
#import StringIO
#foo = X
#fn = StringIO.StringIO(foo) #make a file object from the string
#data = np.loadtxt(fn) #use loadtxt with default settings.
#print(data)

import numpy as np
data = np.loadtxt('../../Single_cell_data/sample_data.txt',delimiter='\t')
print(data)
print(type(data))
print(data.shape)
#------

#digits = datasets.load_digits()
X = data
y = np.r_[1:97]
#y = digits.target
n_samples, n_features = X.shape

'''
def embedding_plot(X, title):
    x_min, x_max = np.min(X, axis=0), np.max(X, axis=0)
    X = (X - x_min) / (x_max - x_min)

    plt.figure()
    ax = plt.subplot(aspect='equal')
    sc = ax.scatter(X[:,0], X[:,1], lw=0, s=40, c=y/10.)

    shown_images = np.array([[1., 1.]])
    for i in range(X.shape[0]):
        if np.min(np.sum((X[i] - shown_images) ** 2, axis=1)) < 1e-2: continue
        shown_images = np.r_[shown_images, [X[i]]]
        ax.add_artist(offsetbox.AnnotationBbox(offsetbox.OffsetImage(digits.images[i], cmap=plt.cm.gray_r), X[i]))

    plt.xticks([]), plt.yticks([])
    plt.title(title)



X_pca = decomposition.PCA(n_components=2).fit_transform(X)
#embedding_plot(X_pca, "PCA")

X_lda = discriminant_analysis.LinearDiscriminantAnalysis(n_components=2).fit_transform(X, y)
#embedding_plot(X_lda, "LDA")

X_tsne = manifold.TSNE(n_components=2, init='pca').fit_transform(X)
#embedding_plot(X_tsne,"t-SNE")

#plt.show()

'''

# other attempt
'''
print(__doc__)

import matplotlib.pyplot as plt

from sklearn import datasets
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

iris = datasets.load_iris()

#X = iris.data
#y = iris.target
target_names = iris.target_names

pca = PCA(n_components=2)
X_r = pca.fit(X).transform(X)

lda = LinearDiscriminantAnalysis(n_components=2)
X_r2 = lda.fit(X, y).transform(X)

# Percentage of variance explained for each components
print('explained variance ratio (first two components): %s'
      % str(pca.explained_variance_ratio_))
'''

'''
plt.figure()
colors = ['navy', 'turquoise', 'darkorange']
lw = 2

for color, i, target_name in zip(colors, [0, 1, 2], target_names):
    plt.scatter(X_r[y == i, 0], X_r[y == i, 1], color=color, alpha=.8, lw=lw,
                label=target_name)
plt.legend(loc='best', shadow=False, scatterpoints=1)
plt.title('PCA of IRIS dataset')

plt.figure()
for color, i, target_name in zip(colors, [0, 1, 2], target_names):
    plt.scatter(X_r2[y == i, 0], X_r2[y == i, 1], alpha=.8, color=color,
                label=target_name)
plt.legend(loc='best', shadow=False, scatterpoints=1)
plt.title('LDA of IRIS dataset')

plt.show()
'''




