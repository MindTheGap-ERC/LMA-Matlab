# LMA-Matlab

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12805858.svg)](https://doi.org/10.5281/zenodo.12805858)

Matlab implementation of the diagenetic reactive-transport model by [L'Heureux (2018)](https://doi.org/10.1155/2018/4968315)

## Authors

__Niklas Hohmann__  
Utrecht University  
email: n.h.hohmann [at] uu.nl  
Web page: [uu.nl/staff/NHohmann](https://www.uu.nl/staff/NHHohmann)  
ORCID: [0000-0003-1559-1838](https://orcid.org/0000-0003-1559-1838)

__Hanno Spreeuw__  
Netherlands eScience Center  
email: h.spreeuw [at] esciencecenter.nl  
Web page: [www.esciencecenter.nl/team/dr-hanno-spreeuw/](https://www.esciencecenter.nl/team/dr-hanno-spreeuw)  
ORCID: [0000-0002-5057-0322](https://orcid.org/0000-0002-5057-0322)

## Usage

Run

```Matlab
run("ScenarioA.m")
```

in the Matlab console. Make sure the file is on your search path. If you want to modify the scenario, go to "scenarioA.m" and modify the paramters.

## Repository structure

* _.gitignore_ : untracked files
* _CITATION.cff_ : citation information
* _LICENSE_ : Apache 2.0 license text
* _LMA_solve.m_ : Solver for the PDE system
* _README.md_ : Readme file
* _ScenarioA.m_ : Parameters for scenario A from L'Heureux (2018)

## Citation

To cite this software please use

* Hohmann, N., & Spreeuw, H. (2024). Matlab implementation if the reactive-transport diagenetic model by L'Heureux (2018) (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.12805858

or see the information in the CITATION.cff file.

## Copyright

Copyright 2023 Netherlands eScience Center and Utrecht University

## License

Apache 2.0 License, see LICENSE file for license text.

## Funding information

Funded by the European Union (ERC, MindTheGap, StG project no 101041077). Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Research Council. Neither the European Union nor the granting authority can be held responsible for them.
![European Union and European Research Council logos](https://erc.europa.eu/sites/default/files/2023-06/LOGO_ERC-FLAG_FP.png)
