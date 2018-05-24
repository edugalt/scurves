#References:

This repository contains codes and data used in the paper: 

Extracting information from S-curves of language change
Fakhteh Ghanbarnejad, Martin Gerlach, José M. Miotto, Eduardo G. Altmann
J. R. Soc. Interface 11, 20141044 (2014) 
http://rsif.royalsocietypublishing.org/content/11/101/20141044
http://arxiv.org/abs/1406.4498

The material of this repository has been originally published as:

https://figshare.com/articles/Timeseries_for_examples_of_language_change_obtained_from_the_googlengram_corpus/1172265

#Data

DATA files (*.csv) in the directory data/

Data of the timeseries used in the paper "Extracting information from S-curves of language change" by F. Ghanbarnejad, M. Gerlach, J. M. Miotto, and E. G. Altmann, [arxiv:1406.4498]. 

Each file contains the timeseries of one example of language change. The file name is given as '<class>_<corpus>_<word>.csv', where:

<class>: refers to the general class of examples of language change ('OrthographicReform1901', 'OrthographicReform1996', 'RussianNames', 'VerbsRegularization');
<corpus>: refers to the language of the googlengram-corpus considered ('de': German, 'en': English);
<word>: refers to the particular case in a class of change, e.g. 'dass' was considered for the change of 'OrthographicReform1901' as well as 'OrthographicReform1996'  in German.

Each file is formatted in the following way:
- the first row contains the header, i.e. 'year' \tab 'rho' \tab 'sigma_s' \newline
- each following row contains a tuple of the timeseries (t,rho,sigma_s), i.e. t \tab rho \tab sigma_s \newline, where

t:       year
rho:     the fraction of times one variant is used over the others
sigma_s: sampling error associated to 'rho' from sampling a finite number of words.

The timeseries are obtained from the raw data of the googlengram-corpus, see the Supplementary Material (Section I: Data) of the paper [arxiv:1406.4498] for details.


#Code

File code.py

The example code shows how to fit the mixed model, Eq. (28) of the Supplementary Material of the paper , to one timeseries (e.g., 'VerbsRegularization_en_abide.csv', i.e. the timeseries of the word 'abide' as an example of the class of verb regularization in the English googlengram corpus) and obtain the parameters of the best fit. It requires numpy, scipy, and matplotlib. 

Open the file Notebook.ipynb (or Notebook.html) for further details. The program will:
- read the timeseries (t,rho,sigma_s),
- numerically find the best set of parameters of the mixed model describing the data,
- plot the timeseries (t,rho) of the data and the timeseries (t_fit, rho_fit) of the mixed model with the obtained parameters. 

