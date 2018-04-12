# Waddle -- Waddington's epigenetic landscape

The _Waddington_ or _epigenetic landscape_ provides a theoretical framework for cell development and differentiation. A figure and citation from Waddington's original paper [1]:

![Part of the epigenetic landscape](pres/ball.PNG?raw=true "Workflow")

_"The path followed by the ball […] corresponds to the developmental history of a particular part of the egg.
There is first an alternative, towards the right or the left. Along the former path, a second alternative is offered; along the path to the left, the main channel continues leftwards, but there is an alternative path which, however, can only be reached over a threshold"_ (Waddington, 1940).


<!---This project aims to develop a set of _Julia_ tools and methods to 
explore the structure of such landscapes for a given model. We use single-cell data as well as simulated data. --->

We provide different tools and methods in Julia to 
(1) simulate, explore and analyse the landscape for a given model and 
(2) visualise the landscape certain genes or gene combinations (after applying dimensionality reduction) given a single-cell data set.

<!--- <img src="https://latex.codecogs.com/svg.latex?\Large&space;x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" title="\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" /> --->

<!---![\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}](https://latex.codecogs.com/svg.latex?x%3D%5Cfrac%7B-b%5Cpm%5Csqrt%7Bb%5E2-4ac%7D%7D%7B2a%7D)--->


## Workflow

![Workflow](pres/workflow.png?raw=true "Workflow")


### Simulations
#### Single-cell data
#### GPU
#### Optimisation: MAP

### Stability Analysis

### Data

### Dimensionality Reduction
Different dimensionality reduction methods are used, ia linear methods such as PCA, SVD, PPCA as well as non-linear methods such as kernelPCA, MDS, ICA and tSNE.

## Contributing
_TODO_
Please read [CONTRIBUTING.md] for details on our code of conduct, and the process for submitting pull requests to us.

<!---## Versioning
_TODO_ --->

## Authors

* **Felicia Burtscher** -- [burfel](https://github.com/burfel)
* **Lucas Ducrot** -- [lucasducrot](https://github.com/lucasducrot)
* **Madeleine Hall**  -- [mgh17](https://github.com/mgh17)
* **Luis Torada** --[lt2216](https://github.com/lt2216)

See also the list of [contributors](https://github.com/waddle-project/contributors) who participated in this project.

<!--- ## License
_TODO_ --->

## Acknowledgments
* [JuliaDiffEq](https://github.com/JuliaDiffEq) and its contributors
* [JuliaPlots](https://github.com/JuliaPlots) and its contributors
<!---* Hat tip to anyone whose code was used --->
* supervisor Prof Michael PH Stumpf for advice and guidance 
* members of the Theoretical Systems Biology Group at IC London, especially Rowan Brackston and Ivan Croydon Veleslavov, for their help and expertise in both the topic and Julia 
* Suhail A Islam for fixing technical issues 
* Prof Michael Sternberg and others involved in conducting and overseeing the MSc in Bioinformatics and Theoretical Systems Biology at IC London

## References

[1] Waddington CH: _Organisers and Genes._ 1940, Cambridge, UK: Cambridge University Press
