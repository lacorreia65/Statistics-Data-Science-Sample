# -*- coding: utf-8 -*-
"""
Created on Sat Jun 22 22:28:28 2019

@author: LuisAlvaro
"""

# Análise de Sentimentos
#
# Patricia Sayuri e Luis Correia
#

### ---------------------------------------------
import mglearn
import numpy as np
import pandas as pd
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.model_selection import train_test_split

train_data = pd.read_csv("data\\train_val_set.csv")
validate_data = pd.read_csv("data\\test_set.csv")

# Aplit Train database in 02 parts (80% - 20%) to validate model selection
df_train,df_test = train_test_split(train_data, test_size = 0.20, stratify=train_data.label, random_state=10)

text_train, y_train = df_train['review'], df_train['label']
text_test, y_test = df_test['review'], df_test['label']

#print("type of text_train: {}".format(type(text_train)))
#print("length of text_train: {}".format(len(text_train)))
#print("text_train[1]:\n{} \nLabel: {}".format(text_train[1]),format(y_train[1]))
#print("text_test[0]:\n{}".format(text_test[0]))

vect = CountVectorizer().fit(text_train)
X_train = vect.transform(text_train)
print("X_train:\n{}".format(repr(X_train)))

# Print Bag-Of-Words for X_Train 
print("bag_of_words: {}".format(repr(X_train)))
print("Dense representation of bag_of_words:\n{}".format(X_train.toarray()))
print("Vocabulary size: {}".format(len(vect.vocabulary_)))
print("Vocabulary content:\n {}".format(vect.vocabulary_))

feature_names = vect.get_feature_names()
print("Number of features: {}".format(len(feature_names)))
print("First 20 features:\n{}".format(feature_names[:20]))
print("Features 100 to 150:\n{}".format(feature_names[100:150]))
print("Every 10th feature:\n{}".format(feature_names[::1000]))

# ---- Building a Classifier - LogisticRegression -----
from sklearn.model_selection import cross_val_score
from sklearn.linear_model import LogisticRegression
scores = cross_val_score(LogisticRegression(class_weight='balanced', C = 0.1), 
                         X_train, y_train,cv=5)
print("Mean cross-validation accuracy: {:.4f}".format(np.mean(scores)))

# Balancing Classes
from sklearn.model_selection import GridSearchCV
param_grid = {'C': [0.001, 0.01, 0.1, 1, 10]}
grid = GridSearchCV(LogisticRegression(class_weight='balanced'), param_grid, cv=5)
grid.fit(X_train, y_train)
print("Best cross-validation score: {:.4f}".format(grid.best_score_))
print("Best parameters: ", grid.best_params_)

# Test Classifier - Try #1 - Logistic Regression
X_test = vect.transform(text_test)
print("*** Test Classifier try #1 (no manipulation) {:.4f}".format(grid.score(X_test, y_test)))

# Using SVM as Classifier
from sklearn import svm
svc = svm.SVC(gamma="scale")
grid = GridSearchCV(svc, param_grid, cv=5)
grid.fit(X_train, y_train)
print("Best cross-validation score: {:.4f}".format(grid.best_score_))
print("Best parameters: ", grid.best_params_)

# Test Classifier - Try #1 - SVM
X_test = vect.transform(text_test)
print("*** Test Classifier try #1 (no manipulation) {:.4f}".format(grid.score(X_test, y_test)))

vect = CountVectorizer(min_df=3).fit(text_train)
X_train = vect.transform(text_train)
print("X_train with min_df: {}".format(repr(X_train)))

feature_names = vect.get_feature_names()
print("First 50 features:\n{}".format(feature_names[:50]))
print("Features 20010 to 20030:\n{}".format(feature_names[20010:20030]))
print("Every 700th feature:\n{}".format(feature_names[::700]))

grid = GridSearchCV(LogisticRegression(), param_grid, cv=5)
grid.fit(X_train, y_train)
print("Best cross-validation score: {:.4f}".format(grid.best_score_))    

## Removing Stop Words
#from sklearn.feature_extraction.text import ENGLISH_STOP_WORDS
#print("Number of stop words: {}".format(len(ENGLISH_STOP_WORDS)))
#print("Every 10th stopword:\n{}".format(list(ENGLISH_STOP_WORDS)[::10]))

# Specifying stop_words="english" uses the built-in list = Sufixo 'NS-NonStop(words)
# We could also augment it and pass our own.
vectNS = CountVectorizer(min_df=5, stop_words="english").fit(text_train)
X_trainNS = vectNS.transform(text_train)
print("X_trainNS with stop words:\n{}".format(repr(X_trainNS)))

gridNS = GridSearchCV(LogisticRegression(class_weight='balanced'), param_grid, cv=5)
gridNS.fit(X_trainNS, y_train)
print("Best cross-validation score (no StopWords): {:.4f}".format(gridNS.best_score_))

# Test Classifier - Try #2
X_testNS = vectNS.transform(text_test)
print("*** Test Classifier try #2 (no StopWords) {:.4f}".format(gridNS.score(X_testNS, y_test)))


#--> TO DO 1 - Try ano ther approach, discarding frequently appearing words, by setting the 
#    max_df option of CountVectorizer to see how it influences the number of features and 
#    the performance

## Rescaling the Data with tf–idf
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.pipeline import make_pipeline
pipe = make_pipeline(TfidfVectorizer(min_df=5, norm=None),LogisticRegression(class_weight='balanced'))
param_grid = {'logisticregression__C': [0.001, 0.01, 0.1, 1, 10]}
grid = GridSearchCV(pipe, param_grid, cv=5)
grid.fit(text_train, y_train)
print("Best cross-validation score: {:.4f}".format(grid.best_score_))

vectorizer = grid.best_estimator_.named_steps["tfidfvectorizer"]
# transform the training dataset
X_train = vectorizer.transform(text_train)
# find maximum value for each of the features over the dataset
max_value = X_train.max(axis=0).toarray().ravel()
sorted_by_tfidf = max_value.argsort()
# get feature names
feature_names = np.array(vectorizer.get_feature_names())
print("Features with lowest tfidf:\n{}".format(feature_names[sorted_by_tfidf[:20]]))
print("Features with highest tfidf: \n{}".format(feature_names[sorted_by_tfidf[-20:]]))

# Finding the words that have low inverse document frequency 
sorted_by_idf = np.argsort(vectorizer.idf_)
print("Features with lowest idf:\n{}".format(feature_names[sorted_by_idf[:100]]))

# Investigating the Model's coefficient -> 25 largest and 25 smallest coefficients
mglearn.tools.visualize_coefficients(
        grid.best_estimator_.named_steps["logisticregression"].coef_,
        feature_names, n_top_features=40)

# Test Classifier - Try #3  *** WITH ERROR ***
X_testTF = vectorizer.transform(text_test)
print("*** Test Classifier try #3 (TF-IDF Adjust) {:.4f}".format(grid.score(X_testTF, y_test)))

## Try #4 - Bag-of-Words with More Than One Word (n-Grams)
pipe = make_pipeline(TfidfVectorizer(min_df=5), LogisticRegression(class_weight='balanced'))
# running the grid search takes a long time because of the
# relatively large grid and the inclusion of trigrams
param_grid = {"logisticregression__C": [0.1, 1, 10],  # [0.001, 0.01, 0.1, 1, 10, 100]
              "tfidfvectorizer__ngram_range": [(1, 3)]} # [(1, 1), (1, 2), (1, 3)]
gridNG = GridSearchCV(pipe, param_grid, cv=5)
gridNG.fit(text_train, y_train)
print("Best cross-validation score: {:.4f}".format(gridNG.best_score_))
print("Best parameters:\n{}".format(gridNG.best_params_))

# extract scores from grid_search
import matplotlib.pyplot as plt
scores = gridNG.cv_results_['mean_test_score'].reshape(-1, 3).T
# visualize heat map
heatmap = mglearn.tools.heatmap(
        scores, xlabel="C", ylabel="ngram_range", cmap="viridis", fmt="%.3f",
        xticklabels=param_grid['logisticregression__C'],
        yticklabels=param_grid['tfidfvectorizer__ngram_range'])
plt.colorbar(heatmap)

# Test Classifier - Try #4  *** WITH ERROR ***
X_testNG = vectorizer.transform(text_test)
print("*** Test Classifier try #4 (n-gram) {:.4f}".format(gridNG.score(X_testNG, y_test)))

# extract feature names and coefficients
vect = gridNG.best_estimator_.named_steps['tfidfvectorizer']
feature_names = np.array(vect.get_feature_names())
coef = gridNG.best_estimator_.named_steps['logisticregression'].coef_
mglearn.tools.visualize_coefficients(coef, feature_names, n_top_features=40)

# find 3-gram features
mask = np.array([len(feature.split(" ")) for feature in feature_names]) == 3
# visualize only 3-gram features
mglearn.tools.visualize_coefficients(coef.ravel()[mask],feature_names[mask], n_top_features=40)


### Testing new Tokenization, Stemming etc
import spacy
import nltk

# load spacy's English-language models
en_nlp = spacy.load('en')
# instantiate nltk's Porter stemmer
stemmer = nltk.stem.PorterStemmer()
# define function to compare lemmatization in spacy with stemming in nltk
def compare_normalization(doc):
    # tokenize document in spacy
    doc_spacy = en_nlp(doc)
    # print lemmas found by spacy
    print("Lemmatization:")
    print([token.lemma_ for token in doc_spacy])
    # print tokens found by Porter stemmer
    print("Stemming:")
    print([stemmer.stem(token.norm_.lower()) for token in doc_spacy])
