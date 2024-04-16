#%%
import pandas as pd
from sklearn.model_selection import train_test_split
# %%
df_X = pd.read_csv('Train_call.txt', sep='\t')
# %%
df_X.head()
# %%
df_X.shape
# %%
df_X.describe()
# %%
df_X['Nclone'].hist()
# %%
df_y = pd.read_csv('Train_clinical.txt', sep='\t')
# %%
df_y.head()
# %%
df_y.Subgroup.unique()
# %%
# we gonna remap the Subgroup to a number
mapper = {
    'HER2+': 0,
    'HR+': 1,
    'Triple Neg': 2
}
# %%
df_y['y'] = df_y.apply(lambda x: mapper[x['Subgroup']], axis=1)

# %%
df_y.head()
# %%
df_X.T
# %%
gene_mapper = {
    'ERBB2': [17, 35104766, 35138441],
    'PGR': [11, 100414313, 100506465],
    'ESR1': [6, 152053324, 152466099],
    'ESR2': [14, 63763506, 63875021],
}
# %%
# select only df rows which fall into ranges or overlap
def region_selector(row):
    # checks if this chromosome has a gene
    gene_match = [
        gene for gene in gene_mapper if gene_mapper[gene][0] == row['Chromosome']
    ] # TODO: change if two genes on the same chromosome
    if gene_match:
        gene_match = gene_match[0]
        gene_range = gene_mapper[gene_match]
        if (row['Start'] <= gene_range[2] and row['End'] >= gene_range[1]) or \
        (row['Start'] >= gene_range[1] and row['End'] <= gene_range[2]) or \
        (row['Start'] <= gene_range[1] and row['End'] >= gene_range[1]):
            row['gene'] = gene_match
            return row
    else:
        return None
# %%
# df_XX = df_X.apply(region_selector, axis=1).dropna()
selected_rows = [
    region_selector(row) for _, row in df_X.iterrows() if region_selector(row) is not None
]
#%%
X_train_selected = pd.concat(selected_rows, axis=1, ignore_index=False).drop(['Chromosome', 'Start', 'End', 'Nclone', 'gene'])
X_train_selected = X_train_selected.astype('int64')
# %%
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
import matplotlib.pyplot as plt
import xgboost as xgb
#%%
cv = StratifiedKFold(n_splits=3, random_state=43, shuffle=True)
#%%
X_train_full = df_X.T.drop(['Chromosome', 'Start', 'End', 'Nclone'])#, 'gene'])
#%%
clf = RandomForestClassifier()
#%%
def perform_crossval_xgb(clf, X_train: pd.DataFrame, y_train: pd.DataFrame, cv: StratifiedKFold):
    scores = []
    # clf = RandomForestClassifier()
    clf = xgb.XGBClassifier(tree_method="hist", early_stopping_rounds=3)
    for train_index, test_index in cv.split(X_train, y_train):
        X_train_cv, X_test_cv = X_train.iloc[train_index], X_train.iloc[test_index]
        y_train_cv, y_test_cv = y_train.iloc[train_index], y_train.iloc[test_index]
        clf.fit(X_train_cv, y_train_cv, eval_set=[(X_test_cv, y_test_cv)], verbose=False)
        scores.append(clf.score(X_test_cv, y_test_cv))
    return scores
#%%
# XGbitionist
clf = xgb.XGBClassifier(tree_method="hist", early_stopping_rounds=3)
scores_full_xgb = perform_crossval_xgb(clf, X_train_full, df_y['y'], cv)
#%%
clf = xgb.XGBClassifier(tree_method="hist", early_stopping_rounds=3)
scores_selected_xgb = perform_crossval_xgb(clf, X_train_selected, df_y['y'], cv)
# %%
clf = RandomForestClassifier()
clf = LogisticRegression()
scores_full = cross_val_score(
    clf, X_train_full, df_y['y'], cv=cv
)    
# %%
clf = RandomForestClassifier()
clf = LogisticRegression()
scores_selected = cross_val_score(
    clf, X_train_selected, df_y['y'], cv=cv
)

# %%
import seaborn as sns
# %%
# plot boxplot with all points 
sns.boxplot(data=[scores_full, scores_selected], showmeans=True)
# plot labels
plt.xticks([0, 1], ['All data', 'Selected data'])
# %%
sns.set_style("whitegrid")
# plot all 4 experiments
sns.boxplot(data=[scores_full, scores_selected, scores_full_xgb, scores_selected_xgb], showmeans=True)
# plot labels
plt.xticks([0, 1, 2, 3], ['All data RF', 'Selected data RF', 'All data XGB', 'Selected data XGB'])
#%%
# plot same but with different color lines connecting points from the same cv round
sns.boxplot(data=[scores_full, scores_selected, scores_full_xgb, scores_selected_xgb], showmeans=True)
# plot labels
plt.xticks([0, 1, 2, 3], ['All data RF', 'Selected data RF', 'All data XGB', 'Selected data XGB'])
# plot lines
sns.swarmplot(data=[scores_full, scores_selected, scores_full_xgb, scores_selected_xgb], color='black')
# connect lines
for i in range(3):
    plt.plot([0, 1, 2, 3], [scores_full[i], scores_selected[i], scores_full_xgb[i], scores_selected_xgb[i]], alpha=0.5)
# %%
# plot scores as scatter plot of different color with x-axis being cv-round
sns.scatterplot(x=range(1, 1+len(scores_full)), y=scores_full, color='blue')
sns.scatterplot(x=range(1, 1+len(scores_selected)), y=scores_selected, color='red')
sns.scatterplot(x=range(1, 1+len(scores_selected_xgb)), y=scores_selected_xgb, color='orange', marker='x', linewidth=1.2)
sns.scatterplot(x=range(1, 1+len(scores_full_xgb)), y=scores_full_xgb, color='green', marker='|', linewidth=1.2)
plt.legend(['All data RF', 'Selected data RF', 'All data XGB', 'Selected data XGB'])
# %%
# show most important features for xgbosst
xgb.plot_importance(clf)
# %%
# print top 10 important features
importance = clf.feature_importances_
features = X_train_full.columns
importance_dict = dict(zip(features, importance))
sorted_importance = sorted(importance_dict.items(), key=lambda x: x[1], reverse=True)
sorted_importance[:50]
# %%
df_X.T.loc[:, 2748]
# %%
