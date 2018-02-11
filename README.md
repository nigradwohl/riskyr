---
title: "riskyr readme"
author: "Hansjörg Neth, SPDS, uni.kn"
date: "2017 12 20"
output: html_document
---


# The `riskyr` package

A toolbox for rendering risk literacy more transparent. 


## Motivation

> Solving a problem simply means representing it<br>
> so as to make the solution transparent. (H.A. Simon)[^1]

[^1]: Simon, H.A. (1996). _The Sciences of the Artificial_ (3rd ed.). The MIT Press, Cambridge, MA. (p. 132).


The issues addressed by `riskyr` are less of a _computational_ than of a _representational_ nature (i.e., concerning the representation in and translation between different formats of information).  Whereas people tend to find it difficult to understand and compute information expressed in terms of _probabilities_, the same information is often easy to understand and compute when expressed in terms of _frequencies_. But rather than just expressing probabilistic information in terms of frequencies, `riskyr` allows translating between formats and illustrates their relationships in a variety of transparent and interactive visualizations.

The basic assumptions and goals driving the development of `riskyr` are the following:

1. Effective training in risk literacy requires simple tools and transparent representations. 

2. We provide a set of (computational and representational) tools that allow various calculations, translations between formats, 
and many alternative views on the interplay between probabilities and frequencies. The functions and representations in the `riskyr` toolbox are not isolated, but complement, support and explain each other.

3. To facilitate active and explorative learning, all visualizations can be explored interactively, 
providing immediate feedback on the effect of changes in parameters.  
By providing many customization options, users can design and shape representations of risk-related information to suit their goals and needs.

<!-- riskyr logo: -->  
<a href = "https://github.com/hneth/riskyr">
<img src = "./inst/pix/riskyr_cube.png" alt = "riskyr" style = "width: 175px; float: right; border:25;"/>
</a>
<!-- ![riskyr](./inst/pix/riskyr_cube.png) --> 
<!-- knitr::include_graphics("./inst/pix/riskyr_cube.png") -->


## Rationale

We begin with some basic variables:

-   the prevalence value `prev` of some condition
-   the sensitivity value `sens` of some decision or diagnostic test (or probability of treatment success)
-   the specificity value `spec` of this decision or diagnostic test (or probability of side effects)
-   a population size `N`

and provide a variety of _perspectives_ on (and representations of) the _consequences_ of and _interplay_ between these variables:

1.  a _data table_ of individual cases  
2.  an _icon array_ (with a population vs. sample view, sorted or randomized)  
3.  a _tree of natural frequencies_  
4.  a _2x2 confusion/contingency table_ (corresponding to the leaves of the frequency tree)  
5.  a _mosaic plot_ that illustrates the proportions of the confusion table  
6.  a curve of _predictive values_ (PPV and NPV) as a function of `prev`  
7.  a plane of _predictive values_ (PPV and NPV) as a function of `sens` and `spec` for a given `prev`  
    <!-- 8. fact boxes (with additional details on benefits and harms of tests or treatments)  -->

A _library of scenarios_ illustrates example cases with known data from the literature.


## Features

### Ontology 

The `riskyr` universe describes the interplay between a total of 10 probabilities (3 of which are essential) 
and 11 frequencies (4 of which are essential). 

### Alternative perspectives

Classification results can be viewed from multiple perspectives, 
which correspond to alternative ways of splitting a population of `N` indivduals into subsets:

1. By _condition_: `TRUE` vs. `FALSE`, then by _decision_: `hi`, `mi`, `fa`, `cr`

2. By _decision_: `pos` vs. `neg`, then by _condition_: `hi`, `mi`, `fa`, `cr`

3. By _correspondence_ between condition and decision: Various metrics of _accuracy_.


### Translating between representational formats

A given _scenario_ is represented both in terms of probabilities and in terms of frequencies. 

A set of _conversion functions_ allow switching back and forth between both formats (i.e., compute frequencies from probabilities and probabilities from frequencies). 


## Package and Application

Our primary objective is to collect and develop a set of basic risk literacy tools in R. To maximise impact, we split our efforts into two complementary projects:

1. The `riskyr` package renders risk literacy more transparent by providing a set of risk-related tools and corresponding representations.

2. The corresponding R Shiny app `riskyrApp` allows using the `riskyr` toolbox in an interactive fashion without any coding.

The combination of package and application facilitates risk communication and supports instruction and training efforts in promoting risk literacy.


## ToDo

Things currently being implemented:

- enriched data set of example scenarios 
- unified options for plotting functions
- additional curves and planes


## About

`riskyr` originated out of a series of lectures and workshops on risk literacy in spring/summer 2017. 
The current version (`riskyr` 0.0.0.916, as of Feb. 8, 2018) is still under development. 
Its primary designers and developers are 
[Hansjörg Neth](https://www.spds.uni-konstanz.de/hans-neth), 
[Felix Gaisbauer](https://www.spds.uni-konstanz.de/felix-gaisbauer), and 
[Nico Gradwohl](https://www.spds.uni-konstanz.de/nico-gradwohl), 
who are researchers at the department of 
[Social Psychology and Decision Sciences](https://www.spds.uni-konstanz.de) at the 
[University of Konstanz](https://www.uni-konstanz.de/en/), Germany. 

The `riskyr` package is open source software written in [R](https://www.r-project.org/) and released under the 
[GPL 2](https://tldrlegal.com/license/gnu-general-public-license-v2) | 
[GPL 3](https://tldrlegal.com/license/gnu-general-public-license-v3-(gpl-3)) licenses. 

Please email at <contact.riskyr@gmail.com>  in case you want to use, adapt, or share this software.


### Contact

We appreciate your feedback, comments, or questions. 

- Please report any `riskyr`-related issues at <https://github.com/hneth/riskyr/issues>.

- For general inquiries, please email us at <contact.riskyr@gmail.com>. 

<!-- uni.kn logo: -->  
<!-- ![](./inst/pix/uniKn_logo.png) --> 
<a href="http://www.uni-konstanz.de">
<img src = "./inst/pix/uniKn_logo.png" alt = "uni.kn.logo" style = "width: 350px; float: right; border:20;"/>
</a>

### Reference

To cite `riskyr` in derivations and publications use:

- Neth, H., Gaisbauer, F., Gradwohl, N., & Gaissmaier, W. (2018).  
  `riskyr`: A toolbox for rendering risk literacy more transparent.  
  Social Psychology and Decision Sciences, University of Konstanz, Germany.  
  Computer software (R package version 0.0.0.916, Feb. 8, 2018).  
  Retrieved from <https://github.com/hneth/riskyr>.  

A BibTeX entry for LaTeX users is: 

    @Manual{,
      title = {riskyr: A toolbox for rendering risk literacy more transparent},
      author = {Hansjörg Neth and Felix Gaisbauer and Nico Gradwohl and Wolfgang Gaissmaier},
      year = {2018},
      organization = {Social Psychology and Decision Sciences, University of Konstanz},
      address = {Konstanz, Germany},
      note = {R package (version 0.0.0.916, Feb. 8, 2018)},
      url = {https://github.com/hneth/riskyr},
      }    
    
Calling `citation("riskyr")` in the package also displays this information.

<!-- eof -->
