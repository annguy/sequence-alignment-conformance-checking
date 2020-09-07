Sequence Alignment and Visualization 
==============================

useful alignment approaches in conformance checking

Project Organization
------------

    ├── LICENSE
    ├── Makefile           <- Makefile with commands like `make data` or `make train`
    ├── README.md          <- The top-level README for developers using this project.
    ├── data
    │   ├── external       <- Data from third party sources.
    │   ├── interim        <- Intermediate data that has been transformed.
    │   ├── processed      <- The final, canonical data sets for modeling.
    │   └── raw            <- The original, immutable data dump.
    │       ├── "CCC19 - Log CSV.csv"  
    │       ├── "CCC19_Log_normative_simulated_v1"
    │       └── "CCC19_Log_normative_simulated_v3"
    ├── docs               <- A default Sphinx project; see sphinx-doc.org for details
    │
    ├── notebooks          <- Jupyter notebooks. Alignment and visualization demo
    │   └── 1.0-wenyu-Sequence_Alignment_and_Visualization.ipynb
    ├── references         <- Data dictionaries, manuals, and all other explanatory materials.
    │
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment
    │
    ├── envtironment.yml   <- The requirements file for reproducing the analysis environment
    │
    ├── setup.py           <- makes project pip installable (pip install -e .) so src can be imported
    ├── src                <- Source code for use in this project.
    │   ├── __init__.py    <- Makes src a Python module
    │   │
    │   ├── preprocess           <- Scripts to preprocess event log
    │   │   └── data_process.py
    │   │
    │   ├── multiple_alignment   <- Scripts of different alignment approaches               
    │   │   ├── pairwise_alignment.py
    │   │   └── multiple_alignment.py
    │   │
    │   └── visualization  <- Scripts to create alignment visualizations
    │       └── visualization.py
    │
    └── tox.ini            <- tox file with settings for running tox; see tox.readthedocs.io


--------



<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>
