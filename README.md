# Sequence Alignment and Visualization 


Useful alignment approaches in conformance checking sequence, inspired by approaches in bioinformatics.
This project was conducted at the Machine Learning and Data Analytics Lab (CS 14), Department of Computer Science Friedrich-Alexander-University Erlangen-Nuremberg (FAU).
Data is provided by the 2019 Conformance Checking Challenge (CCC19).


## Citation and Contact
You find a PDF of the paper
[https://arxiv.org/abs/2010.11719](https://arxiv.org/abs/2010.11719).

If you use our work, please also cite the paper:
```
 @article{Nguyen_Zhang_Schwinn_Eskofier_2020, 
 title={Conformance Checking for a Medical Training Process Using Petri net Simulation and Sequence Alignment}, 
 url={http://arxiv.org/abs/2010.11719}, note={00000 
arXiv: 2010.11719}, journal={arXiv:2010.11719 [cs]}, 
author={Nguyen, An and Zhang, Wenyu and Schwinn, Leo and Eskofier, Bjoern}, 
year={2020}
}

```
If you would like to get in touch, please contact [an.nguyen@fau.de](mailto:an.nguyen@fau.de).


## Getting Start 
### Installing
A step by step series of examples that tell you how to get a development env running
Say what the step will be
```
git clone https://github.com/annguy/sequence-alignment-conformance-checking.git
```
Create new environment using the environment.yml 
```
conda env create -f environment.yml
```
Install src in the new environment 
```
pip install -e
```
### Examples
1. Pairwise Alignment

For example, seq1 =  "HEAHEE", seq2 = "PAHE"
```
from pairwise_alignment import Needleman_Wunsch
Needleman_Wunsch(seq1, seq2, gap_open_penalty=-2, gap_extend_penalty=-2, match_reward=1, mismatch_penalty=-2)
```
```
(['H', 'E', 'A', 'H', 'E', 'E'], ['P', '-', 'A', 'H', 'E', '-'], -3.0, 0.5)
```
The structure of result is (align1, align2, score, identity)

2. Visualization
```
from visualization import view_alignment
test = ['a-----bcdefghijkno-p-qrst--uvxyz012',
 'bflngcbcde--hilfaolpfqrstvwuvwyz012']

ID =['54%',
 '1_Pre']
 
p = view_alignment(test, ID,plot_width=1000,fontsize="16pt",text_font_size="19pt",height_adjust=58)
show(p)
```
![Image text](https://github.com/annguy/sequence-alignment-conformance-checking/blob/master/reports/figures/test.png)  

3. Save as "png" or "svg"
```
from bokeh.io import export_png
export_png(p,filename="test.png",height=100, width=1200)
```
```
p.output_backend = "svg"
export_svgs(p, filename="test.svg")
```


## Project Organization


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
    └── src                <- Source code for use in this project.
       ├── __init__.py    <- Makes src a Python module
       │
       ├── preprocess           <- Scripts to preprocess event log
       │   └── data_process.py
       │
       ├── multiple_alignment   <- Scripts of different alignment approaches               
       │   ├── pairwise_alignment.py
       │   └── multiple_alignment.py
       │
       └── visualization  <- Scripts to create alignment visualizations
            └── visualization.py


<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>

## Contributor 
[An Nguyen](https://www.mad.tf.fau.de/person/an-nguyen/)

Wenyu Zhang  


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details



