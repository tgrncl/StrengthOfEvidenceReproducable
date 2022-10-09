## Strength of Evidence
Reproducible R Script for Binder
Paper: Why Most Results of Socio-Technical Security User Studies Are False -- And what to do about it [1]

Thomas Groß
School of Computing 
Newcastle University
Newcastle upon Tyne
United Kingdom


## Summary

This Binder project analyses a dataset with effect sizes from scientific papers in security and privacy in the years 2006-2016.
It computes measures of strength of evidence for the statistical tests involved.
It spawns a Shiny server and client application as a User Interface.


## How to use

### Open RStudio

The Binder version of the software will load from
   https://mybinder.org/v2/gh/tgrncl/StrengthOfEvidenceReproducable.git/HEAD?urlpath=rstudio
   
### Open the Shiny app

Load and execute **app.R**

Note that the app will pre-compute a simulation of power and strength-of-evidence
values to improve responsiveness of the Shiny app. This will take a few minutes 
during the start-up of the Shiny app.

The Shiny app will open in a new browser window. Ensure that you permit the
pop-up window.

### Interact with the visualisation

You can choose the type of graph being displayed.
You can further set parameters for effect size thresholds, bias and prior.


## Aim

The aim of the software is to evaluate the strength of evidence of cyber security user studies.

The strength of evidence is measured in difference ways:
 - *How likely is a reported result true in reality, given the observations of the experiment?*
   (This is measured in the Positive Predictive Value (PPV) and the False Positive Rate (FPR).)
 - *What is the strength of evidence of the result itself independent of any prior knowledge/probabilities?*
   (We can measure this with a likelihood ratio of the evidence for and against the effect reported.)
 - *What prior knowledge/probability would we needed to have had to be confident in the reported result?*
   (We can evaluate this with what is called the Reverse Bayesian Prior. Essentially, this approach fixes an assumed likelihood of the result being true in reality and traces back how likely we must have assumed the result before the experiment to reach that conclusion.)


## Method

### Sample

The sample for this work is based on a Systematic Literature Review (SLR) conducted by Kovila Coopamootoo and Thomas Groß in 2017 [2].
The specification of the SLR and the retrieved sample incl. the literature list was published with a 2019 STAST paper [3] of Groß.
We include that specification in this repository as **SLR_spec_sample.pdf**.

Coopamootoo and Groß investigated cyber security user studies from relevant venues in the years 2006-2016.
In their work, they evaluated the user studies with respect to qualitative completeness of reporting.

### Analysis Approach

1. The  analysis is based on the given SLR sample of papers. For each paper, sample sizes are extracted.
2. Test statistics of statistical inferences were extracted in two ways (i) with the statistical evaluation tool **statcheck** and (ii) by coding descriptive statistics by hand.
3. From these test statistics, standard effect sizes are computed.
4. We also simulate statistical power, and strength of evidence measures, based on parameters effect size threshold, bias and prior.
5. Based on these analyses, we enable the user to evaluate what the strength of evidence of the field may be based on the users' own assumptions on different parameters. 

### Visualization

The Shiny app in this repository makes it possible to check the computations and explore different parameter constellations.


## License
Copyright (c) 2022 Thomas Groß, tgrncl at github.com

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.



## References

[1] T. Groß. Why Most Results of Socio-Technical Security User Studies Are False -- And what to do about it. In proceedings of the 12th Workshop on Socio-Technical Aspects in Security (STAST), 2022.
[2] K. Coopamootoo and T. Groß. Systematic evaluation for evidence-based methods in cyber security. Technical Report TR-1528, Newcastle University, 2017
[3] T. Groß. Fidelity of Statistical Reporting in 10 Years of Cyber Security User Studies. In proceedings of the 9th International Workshop on Socio-Technical Aspects in Security (STAST), 2019.