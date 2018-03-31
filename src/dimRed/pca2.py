import numpy as np
import matplotlib.pyplot as plt
from matplotlib import offsetbox
from sklearn import manifold, datasets, decomposition, discriminant_analysis
 
# f = open('../Single_cell_data/sample_data.txt','r')
# vocab_temp = f.read().split()
# f.close()
# col = len(vocab_temp)
# print("Training column size:")
# print(col)

# #dataset = list()

# row = run('cat '+'/Users/ya/Documents/10percent/X_true.txt'+" | wc -l").split()[0]
# print("Training row size:")
# print(row)
# matrix_tmp = np.zeros((int(row),col), dtype=np.int64)
# print("Train Matrix size:")
# print(matrix_tmp.size)
#         # label_tmp.ndim must be equal to 1
# label_tmp = np.zeros((int(row)), dtype=np.int64)
# f = open('/Users/ya/Documents/10percent/X_true.txt','r')
# count = 0
# for line in f:
#     line_tmp = line.split()
#     #print(line_tmp)
#     for word in line_tmp[0:]:
#         if word not in vocab_temp:
#             continue
#         matrix_tmp[count][vocab_temp.index(word)] = 1
#     count = count + 1
# f.close()
# print("Train matrix is:\n ")
# print(matrix_tmp)
# print(label_tmp)
# print(len(label_tmp))
# print("No. of topics in train:")
# print(len(set(label_tmp)))
# print("Train Label size:")
# print(len(label_tmp))

# def embedding_plot(X, title):
#     x_min, x_max = np.min(X, axis=0), np.max(X, axis=0)
#     X = (X - x_min) / (x_max - x_min)
 
#     plt.figure()
#     ax = plt.subplot(aspect='equal')
#     sc = ax.scatter(X[:,0], X[:,1], lw=0, s=40, c=y/10.)
 
#     shown_images = np.array([[1., 1.]])
#     for i in range(X.shape[0]):
#         if np.min(np.sum((X[i] - shown_images) ** 2, axis=1)) < 1e-2: continue
#         shown_images = np.r_[shown_images, [X[i]]]
#         ax.add_artist(offsetbox.AnnotationBbox(offsetbox.OffsetImage(digits.images[i], cmap=plt.cm.gray_r), X[i]))
 
#     plt.xticks([]), plt.yticks([])
#     plt.title(title)

# import codecs
# from sklearn.decomposition import TruncatedSVD
# from sklearn.feature_extraction.text import CountVectorizer
# with codecs.open('../Single_cell_data/sample_data.txt', 'r', encoding='utf-8') as inf:
#     lines = inf.readlines()
# vectorizer = CountVectorizer(binary=True)
# X_train = vectorizer.fit_transform(lines)
# perform_pca = False
# if perform_pca:
#     n_components = 100
#     pca = TruncatedSVD(n_components)
#     X_train = pca.fit_transform(X_train)

# #X_pca = decomposition.PCA(n_components=2).fit_transform(X)
# #embedding_plot(X_train, "PCA")
# plt.plot(X_train)
# plt.show()



from sklearn.decomposition import PCA

with open('../Single_cell_data/sample_data.txt','r') as f:
	temp = f.readlines()
X = np.array(temp)
pca = PCA(n_components=2)
pca.fit(X)
PCA(copy=True, iterated_power='auto', n_components=2, random_state=None, svd_solver='auto', tol=0.0, whiten=False)
print(pca.explained_variance_ratio_)  
print(pca.singular_values_)  

print(temp)
