# ETBD
ETBD simulation 


## Setup
Requirements: R 3.0 or greater with the following packages installed and the following scripts sourced:

```R
package.vector = c('ape','caper','phyltools','plyr','picante','geiger','MASS','sads','stringr')
install.packages(package.vector)
source('code/run_sim.r')
source('code/analyze_sim.r')
source('code/senc_sim_fun.r')
source('code/senc_analysis_fun.r')
source('code/supplemental_analysis_functions.r')
```



## Run simulations
The parameters governing a given simulation are specified here. These parameters include:

<table>
  <tr>
    <td>t</td>
    <td>number of time steps</td>
  </tr>
  <tr>
    <td>pallo</td>
    <td>probability of allopatric spectiaion, only works if 'probleave' is more than 0</td>
  </tr>
  <tr>
    <td>psymp</td>
    <td>probability of sympatric specation</td>
  </tr>
  <tr>
    <td>SAD</td>
    <td>if TRUE indicates a fishers log-series species abundnace distribution</td>
  </tr>
  <tr>
    <td>GEO</td>
    <td>if TRUE indicates a geometric species abundnace distribution</td>
  </tr>
  <tr>
    <td>LOGNORM</td>
    <td>if TRUE will indicate a log-normal species abundnace distribution</td>
  </tr>
  <tr>
    <td>watchgrow</td>
    <td>if TRUE allows animation of tree growth</td>
  </tr>
  <tr>
    <td>SADmarg</td>
    <td>the margin of error used in the SAD forcing mechanism</td>
  </tr>
  <tr>
    <td>siteN</td>
    <td>number of sites (needs to be the same number os the length of JmaxV) or a square if isGrid is TRUE)</td>
  </tr>
  <tr>
    <td>probleave</td>
    <td>probability that a species will move to an adjacent site</td>
  </tr>
  <tr>
    <td>JmaxV</td>
    <td>vector of productivity zones (Jmax gradient) must be same length as numbe of sites</td>
  </tr>
   <tr>
    <td>NegExpEx</td>
    <td>if TRUE extinction is modeled as a negaitive exponential function</td>
  </tr>
    <tr>
    <td>exparm</td>
    <td>if NegExpEx is TRUE this is the negative exponential parameter for extinciton</td>
  </tr>
  <tr>
   <td>split</td>
    <td>a small portion (10%) of the parent species will split off and migrate, default = TRUE</td>
  </tr>
  <tr>
     <td>bud</td>
    <td>if TRUE a small portion (10%) of the parent species will emerge leaving the origianl population of the parent species the same</td>
  </tr>
    <td>isGrid</td>
    <td>if TRUE the space configuration turns from linear to a grid conected with chebyshev distance</td>
  </tr>
</table>








##simulation output 


The output files produced are listed below.

<table>
  <tr>
    <td>SIMULATION OUTPUT:</td>
    <td>a list of results</td>
  </tr>
  <tr>
    <td></td>
    <td>tree:   final tree</td>
  </tr>
  <tr>
    <td></td>
    <td>trees:   collection of every tree for time step</td>
  </tr>
  <tr>
    <td></td>
    <td>matrix_list:   final matrix of species locations and abundanecs</td>
  </tr>
  <tr>
    <td></td>
    <td>mig:   collection of matrix lists for every time step</td>  
  </tr>
  <tr>
    <td></td>
    <td>extincts:   list of all species that have gone extinct.</td>  
  </tr>
  <tr>
    <td></td>
    <td>allop:   matrix list of all species that allopatrically speciated</td>  
  </tr>
  <tr>
    <td></td>
    <td>symp:   matrix list of all species that have sympatrically speciated </td>  
  </tr>
  <tr>
    <td></td>
    <td>allotrip:   list of every allopatrically speciating species per time step</td>  
  </tr>
  <tr>
    <td></td>
    <td>symptrip:   list of every sympatrically speciating species per time step</td>  
  </tr>
  <tr>
    <td></td>
    <td>exsp:    matrix lists of extinct species per site </td>
  </tr>
</table>




















