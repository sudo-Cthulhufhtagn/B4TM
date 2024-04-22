# B4THDWN

This repository contatins solutions for B4TM course project of team 22.

- [B4THDWN](#b4thdwn)
  - [Research question](#research-question)
    - [Short version](#short-version)
    - [Elaborate explanation](#elaborate-explanation)
  - [Data sources](#data-sources)
  - [Task devision](#task-devision)
  - [Methodology](#methodology)
  - [Running code instructions](#running-code-instructions)
  - [Positions of genes](#positions-of-genes)

## Research question
### Short version
Is there a differece in prediction accuracy and feature importance if array CGH data is used only with locations of receptors in question (HER2, ER, PR). 

If there is a difference, find the locations and map them to genes which affect 
### Elaborate explanation
The general question in so called [instructions file](https://canvas.vu.nl/courses/75132/files/7501567?wrap=1) (which contains no strict instructions but rather vague requirements), is to "play with the data", which it is more of an engineering task, but ok... We are going into more details to understand whether the locations of the genes of interest are enough to precisely predict type of canser, or multi-region dependencies help increase accuracy(because region of the gene is sometimes smaller than what aCGH captures, but multi-site context might help improve metric and even find related genes).

## Data sources
* [train features](https://canvas.vu.nl/courses/75132/files/7501496/preview)
* [train labels](https://canvas.vu.nl/courses/75132/files/7501475/preview)
* [location to gene mapping](https://canvas.vu.nl/courses/75132/files/7501551?wrap=1)

## Task devision

Steps:
 1. Fit the model on all features
 2. Fit the model on features(locations) of genes of interest
 3. Perform t-test(or non-parametric) on 5-fold cross validation between points 1 and 2.
      
If p-value is <0.05, then:

  1. Do feature explanation indepentently for each predicted class, find regions which are important for whole dataset and identify genes there and check their funciotn.

> Somebody play with crossval, feature selection, hyperparameter tuning to achieve higher accuracy performance.

## Methodology
Linear model, random forest and xgboost will be utilized to fit both full data and selected regions, to test whether null hypothesis is valid. Following that, feature importance will be performed to identify regions which also affect improved accuracy on full data(if observed), and gene mappings should be performed to identify which genes might be affected and affecting the condition.

## Running code instructions
Create environment(with  either conda or venv/virtualenv) or proceed with default python.

1. Install requirements:
```
pip install -r requirements.txt
```

## Positions of genes

| Receptor | Corresponding gene(s) | Chr | Start | Stop | EnsID | Inside regions |
|---|---|---|---|---|---|---|
|HER2|ERBB2|17|35104766|35138441|ENSG00000141736|1|
|PR|PGR|11|100414313|100506465|ENSG00000082175|1|
|ER|ESR1|6|152053324|152466099|ENSG00000091831|3|
|ER|ESR2|14|63763506|63875021|ENSG00000140009|1|