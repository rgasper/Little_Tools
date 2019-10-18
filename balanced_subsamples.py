# -*- coding: utf-8 -*-
"""
balanced_subsamples.py
    :author: Raymond Gasper 
    :created: 18 October 2019
"""
import logging
import configs
import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype
log = logging.getLogger(configs.logfile)

def balanced_subsample_np(X, y, sample_size=None, random_seed=None, imbalance_factor=1):
    """ return a balanced data set by sampling all classes with sample_size 
        current version is developed on assumption that the positive
        class is the minority.
    from https://stackoverflow.com/a/38563000 with modification

    Parameters:
    ===========
    X: {numpy.ndarrray, 2d}
    y: {numpy.ndarray, 1d}
    """
    assert isinstance(X, np.ndarray) and isinstance(y, np.ndarray)
    assert len(X.shape) == 2
    assert len(y.shape) == 1

    uniq_labels = np.unique(y)
    log.debug('target(y) vector has {} distinct labels'.format(uniq_labels.shape[0]))
    uniq_counts = {label: sum(y == label) for label in uniq_labels}
    
    min_elems = min([val for val in uniq_counts.itervalues()])
    log.debug('class with the least elements has {} elements'.format(min_elems))
    if sample_size is None:
        sample_size = int(min_elems*imbalance_factor)
    log.debug('sampling {} elements per class'.format(sample_size))

    # just re-seeds the generator "randomly" if fed None
    np.random.seed(random_seed)

    # find observation index of each class labels
    groupby_labels = {}
    for ii, label in enumerate(uniq_labels):
        obs_idx = [idx for idx, val in enumerate(y) if val == label]
        groupby_labels[label] = obs_idx

    # oversampling on observations of each label
    balanced_copy_idx = []
    for gb_label, gb_idx in groupby_labels.iteritems():
        over_sample_idx = np.random.choice(gb_idx, size=sample_size, replace=True).tolist()
        # if sample size is bigger than the number of available indices, it will repeat entries
        over_sample_idx = list(np.unique(over_sample_idx))
        balanced_copy_idx += over_sample_idx
    np.random.shuffle(balanced_copy_idx)

    return (X[balanced_copy_idx, :], y[balanced_copy_idx], balanced_copy_idx)


def balanced_subsample_pd(df, y_col, sample_size=None, random_seed=None, imbalance_factor=1):
    """ Warning! Can only take in numeric target data. Must convert to numeric labels before use.
        return a balanced data set by sampling all classes with sample_size 
        current version is developed on assumption that the positive
        class is the minority.
    from https://stackoverflow.com/a/38563000 with modification

    Parameters:
    ===========
    df: {pd.Dataframe}
    y_col: {column name : str}

    NOTE(gasperr) sample size doesn't really work as intended right now,
                    it will not balance at all if you supply it
    """
    assert isinstance(df, pd.core.frame.DataFrame)

    # just re-seeds the generator "randomly" if fed None
    np.random.seed(random_seed)
    class_elems = []
    min_elems = None

    uniq_ys = np.unique(df[y_col].values)

    # take observation elements of each label
    for yi in uniq_ys:
        elems = df.loc[df[y_col] == yi]
        class_elems.append((yi, elems))
        if min_elems == None or elems.shape[0] < min_elems:
            min_elems = elems.shape[0]

    if sample_size is None:
        sample_size = int(min_elems*imbalance_factor)

    # sampling on observations of each label
    dfs = []
    for label,df_this_label in class_elems:
        if df_this_label.shape[0] <= sample_size:
            elems = df_this_label.sample(n=df_this_label.shape[0], replace=False)
        else:
            elems = df_this_label.sample(n=sample_size, replace=False)
        dfs.append(elems)

    return pd.concat(dfs)


def test_balanced_sample_makers():
    ''' a little test function to make sure things are working. '''
    # load datasets
    # 50 samples of three different iris flower types, grouped in the dataset rows by iris type
    # sklearn data comes pre-split into features and target
    from sklearn.datasets import load_iris
    iris_pd = pd.read_csv('https://raw.githubusercontent.com/mwaskom/seaborn-data/master/iris.csv')
    iris_skl = load_iris()
    iris_np_X = iris_skl.data
    iris_np_y = iris_skl.target

    # remove entries from the first iris type in both datasets
    keep = 25
    z = 50-keep
    drop_inds = np.arange(0,z)
    iris_pd = iris_pd.drop(drop_inds)
    iris_np_X = np.delete(iris_np_X, list(drop_inds), axis=0)
    iris_np_y = np.delete(iris_np_y, list(drop_inds), axis=0)

    # split into feature, target for pandas. convert target to numeric
    features, target = [u'sepal_length', u'sepal_width', u'petal_length', u'petal_width'], u'species'
    iris_pd_feats, iris_pd_targ = iris_pd[features], iris_pd[target]
    iris_pd[target] = iris_pd[target].replace({'setosa':0, 'versicolor':1, 'virginica':2}) 

    # make some balanced datasets
    # and evaluate whether the test succeeded. BEWARE! Test assertions use assumptions.
    try:
        iris_pd_bal = balanced_subsample_pd(iris_pd, target)
        assert all([val == keep for val in iris_pd_bal[u'species'].value_counts().values])
        log.info("Test of dataset balancing functions succeeded! Works for pandas!")
    except:
        log.exception("SAD. dataset balancing for pandas failed.")
    try:
        iris_np_X_bal, iris_np_y_bal, iris_np_bal_inds  = balanced_subsample_np(iris_np_X, iris_np_y)
        assert iris_np_X_bal.shape[0] == iris_np_y_bal.shape[0]
        assert iris_np_X_bal.shape[0] == keep*3
        log.info("Test of dataset balancing functions succeeded! Works for numpy!")
    except:
        log.exception("SAD. dataset balancing for numpy failed.")


if __name__ == '__main__':
    test_balanced_sample_makers()
