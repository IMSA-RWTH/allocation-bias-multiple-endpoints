# multEndp

The function <b> 'multEndp'</b> is an intuitive tool for assessing the impact of allocation bias on testing decisions in multiple endpoint clinical trials. This function aims to provide better scientific justification for the selection of an randomization procedure (RP) in multiple endpoint clinical trials.

## Motivation


Randomization is one of the most important design features of a clinical trial to prevent allocation bias. Selection of a RP on scientific arguments, i.e. quantification of the potential to mitigate allocation bias, increases the validity of a clinical trial. In clinical practice, the choice of a RP  is usually not based on scientific arguments. The function <b> 'multEndp'</b> can base the choice of randomization procedures in the planning phase of a clinical trial on scientific rationales.


## Usage

```
library(randomizeR)

multEndp(N=32,m=2,sigma=c(1,1), randomization='CR', eta=c(0.1024,0.1024), r_number=100000,procedure='Sidak')

#' Output:
#'    N m  eta.1  eta.2 randomization  mean_FWER go_nogo p
#'   32 2 0.1024 0.1024            CR 0.05018712 0.55387 0
```

## License
This project is free and open source software, licensed under GPL3.

