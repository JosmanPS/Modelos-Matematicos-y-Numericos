import numpy as np
import pandas as pd
import svm


data = pd.read_table('wdbc.data', sep=',', header=None)
X = np.array(data.ix[:, 2:])
y = data.ix[:, 1]
y = np.array([1.0 if label == 'B' else -1.0 for label in y]).ravel()
n_samples, n_features = X.shape

c = 0.1
trainer = svm.SVMTrainer(svm.Kernel.linear(), c)
predictor = trainer.train(X, y)

test = X[0].reshape(1, n_features)
predictions = [predictor.predict(sample.reshape(1, n_features))
               for sample in X]
