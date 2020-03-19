# Quantifying the non-reported new daily COVID-2019 cases by region in Spain 

#### Webpage: [underreported.github.io]{underreported.github.io} 

<p align="justify"> The present outbreak of  COVID-19 disease, caused by the SARS-CoV-2 virus, has put the planet in quarantine. On January 30, 2020, the World Health Organization (WHO) declared the COVID-19 outbreak a "public health emergency of international concern", and then a pandemic on March 11.</p>

<p align="justify"> Spain has become the fifth country worldwide with more infected cases, officially registering over 13 thousand cases in a short time. Although many critical and severe measures have been considered from the authorities to lessen the impact of the outbreak and help flatten the curve, they rely on numbers that could be unreliable and therefore misrepresent the implications of such pandemic. </p>

Counts in Spain due to the protocols used for testing, mainly include individuals with severe symptoms. The authorities have juste announced a new protocol with rapid tests to be implementend in a few days [elpais.com](https://elpais.com/sociedad/2020-03-18/el-numero-de-personas-contagiadas-por-coronavirus-crece-hasta-las-13716-un-18-mas-que-hace-un-dia.html).

<p align="justify"> Given the nature of our data, we can guess that the estimated number of cases that we are finding are in fact potentially severe cases, and presumably the size of the infected population (asymptomatic) is even higher.</p>

<p align="justify"> Accordingly, the current analysis aims to update the situation concerning COVID-19 daily, and particularly quantify the potential under-reporting in the official registered cases by region in Spain. Results herein can help to have a more realistic picture of the pandemic at a real time as well as to more accurately estimate essential measures such as the basic reproduction number or the fatality rate that are used for practitioners and politicians to make decisions.</p>

The data for the analysis have been  extracted from [eldiario.es](https://www.eldiario.es/sociedad/Consulta-evolucion-coronavirus-expansion-Espana_0_1005099739.html#mapaccaa), where official data are gathered.

<p align="justify"> Notice that this analysis con be easily reproduced for other countries. </p>

[Figures 1 (a)](https://github.com/underreported/COVID19_UR/blob/master/daily-plots/17-03-2020/Figure1a.pdf) and [(b)](https://github.com/underreported/COVID19_UR/blob/master/daily-plots/17-03-2020/Figure1b.pdf) and Table 1 provide a description and evolution of daily counts by region. 

<table class="table table-striped" style="width: auto !important; ">
<caption>Table 1: Summary of the daily COVID-19 cases from 27-02-20 to 17-03-2020 by region in Spain</caption>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Andalucia </th>
   <th style="text-align:right;"> Aragon </th>
   <th style="text-align:right;"> Asturias </th>
   <th style="text-align:right;"> Canarias </th>
   <th style="text-align:right;"> Cantabria </th>
   <th style="text-align:right;"> Castilla Leon </th>
   <th style="text-align:right;"> Catalunya </th>
   <th style="text-align:right;"> Extremadura </th>
   <th style="text-align:right;"> Galicia </th>
   <th style="text-align:right;"> La Rioja </th>
   <th style="text-align:right;"> Navarra </th>
   <th style="text-align:right;"> Pais Vasco </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> minimum </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> mean </td>
   <td style="text-align:right;"> 34.15 </td>
   <td style="text-align:right;"> 10.40 </td>
   <td style="text-align:right;"> 9.65 </td>
   <td style="text-align:right;"> 7.40 </td>
   <td style="text-align:right;"> 2.90 </td>
   <td style="text-align:right;"> 21.55 </td>
   <td style="text-align:right;"> 69.70 </td>
   <td style="text-align:right;"> 7.65 </td>
   <td style="text-align:right;"> 14.60 </td>
   <td style="text-align:right;"> 17.75 </td>
   <td style="text-align:right;"> 15.65 </td>
   <td style="text-align:right;"> 52.45 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> median </td>
   <td style="text-align:right;"> 7.00 </td>
   <td style="text-align:right;"> 3.50 </td>
   <td style="text-align:right;"> 2.00 </td>
   <td style="text-align:right;"> 4.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 6.50 </td>
   <td style="text-align:right;"> 16.00 </td>
   <td style="text-align:right;"> 1.00 </td>
   <td style="text-align:right;"> 1.50 </td>
   <td style="text-align:right;"> 14.00 </td>
   <td style="text-align:right;"> 0.50 </td>
   <td style="text-align:right;"> 23.50 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> maximum </td>
   <td style="text-align:right;"> 168.00 </td>
   <td style="text-align:right;"> 67.00 </td>
   <td style="text-align:right;"> 45.00 </td>
   <td style="text-align:right;"> 29.00 </td>
   <td style="text-align:right;"> 20.00 </td>
   <td style="text-align:right;"> 97.00 </td>
   <td style="text-align:right;"> 491.00 </td>
   <td style="text-align:right;"> 42.00 </td>
   <td style="text-align:right;"> 80.00 </td>
   <td style="text-align:right;"> 53.00 </td>
   <td style="text-align:right;"> 91.00 </td>
   <td style="text-align:right;"> 213.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> standard deviation </td>
   <td style="text-align:right;"> 51.49 </td>
   <td style="text-align:right;"> 16.47 </td>
   <td style="text-align:right;"> 13.82 </td>
   <td style="text-align:right;"> 8.66 </td>
   <td style="text-align:right;"> 5.41 </td>
   <td style="text-align:right;"> 29.84 </td>
   <td style="text-align:right;"> 120.95 </td>
   <td style="text-align:right;"> 12.44 </td>
   <td style="text-align:right;"> 23.65 </td>
   <td style="text-align:right;"> 17.80 </td>
   <td style="text-align:right;"> 24.42 </td>
   <td style="text-align:right;"> 67.91 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> dispersion index </td>
   <td style="text-align:right;"> 77.63 </td>
   <td style="text-align:right;"> 26.09 </td>
   <td style="text-align:right;"> 19.80 </td>
   <td style="text-align:right;"> 10.13 </td>
   <td style="text-align:right;"> 10.09 </td>
   <td style="text-align:right;"> 41.32 </td>
   <td style="text-align:right;"> 209.89 </td>
   <td style="text-align:right;"> 20.24 </td>
   <td style="text-align:right;"> 38.32 </td>
   <td style="text-align:right;"> 17.85 </td>
   <td style="text-align:right;"> 38.11 </td>
   <td style="text-align:right;"> 87.93 </td>
  </tr>
</tbody>
</table>

<p align="justify">  If the under-reporting is ignored, the daily counts can be appropriately modeled following: <img src="https://render.githubusercontent.com/render/math?math=exp(\alpha_0 + \alpha_1t)"> , since the number of daily COVID-19 cases overtime properly growths exponentially according to Figure 1.</p>
 
<p align="justify"> However, if we consider that the official number of daily cases does not reflect the total number of cases (e.g., a proportion of the cases is not observed, and thus the data are misreported), the model above does not make any sense, and therefore a more appropriate alternative should be considered. </p>

We shall base all the subsequent analysis in a model  introduced by [Fernández-Fontelo et al. (2016)](https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.7026). 

<p align="justify"> In that model, two different processes are considered: <img src="https://render.githubusercontent.com/render/math?math=X_n">  which is the true process but unobserved (latent), and <img src="https://render.githubusercontent.com/render/math?math=Y_n">  which is observed and potentially under-reported. In this application, the latent process is assumed to be Poisson distributed with time-dependent rate, <img src="https://render.githubusercontent.com/render/math?math=\lambda_t=exp(\beta_0 + \beta_1t)"> . The observed process will always be lower or equal than the latent process (due to the under-reporting) in such a way that <img src="https://render.githubusercontent.com/render/math?math=Y_n">  will be equal than <img src="https://render.githubusercontent.com/render/math?math=X_n">  (non under-reporting) with probability <img src="https://render.githubusercontent.com/render/math?math=1-\omega">; or <img src="https://render.githubusercontent.com/render/math?math=Y_n"> is <img src="https://render.githubusercontent.com/render/math?math=q \circ X_n">  with probability <img src="https://render.githubusercontent.com/render/math?math=\omega"> . Parameters <img src="https://render.githubusercontent.com/render/math?math=\omega">  and <img src="https://render.githubusercontent.com/render/math?math=q">  quantify the overall frequency and intensity of the phenomenon, which roughly speaking describe respectively the number of times the observed counts are not equal to the real ones, and the distance between the real and observed processes. </p>

<p align="justify"> Table 2 shows the estimates of the models parameters by region. They can be interpreted as follows quickly. For instance, for Catalunya, the overall frequency and intensity of under-reporting are roughly 0.55 and 0.37, respectively. This means that 55% of the counts throughout the period are not entirely reported and that averagely the proximity between the real and observed processes is 0.37 (being 1 when two processes are identical). </p>

(TO DO: TABLES AND FIGURES 17-03)

<p align="justify"> Using the Viterbi algorithm, the model also enables reconstructing the most likely sequence of real COVID-19 cases throughout the study. This allows us to have an estimated time series of truly daily cases and evaluate the impact of under-reporting over measures such as the basic reproduction number. Figure 2 shows the observed and reconstructed series over time by region. </p>

Table 3 shows the percentages of means counts that are not covered by the official registers. Thus, the highest the rate, the lower is the coverage, and therefore the severe is the impact of the under-reporting. 

(TO DO: TABLES AND FIGURES 17-03)

It is instructive to see what the difference would be on epidemic spread by fitting an epidemic model to the reconstructed series of counts and the observed counts recorded by public agencies. We fit the classic SIR (Susceptible-Infectious-Recovered) model. Table 4 shows the basic reproduction rate by using the reconstructed series (RE) and the observed (RR). 

## Updates
This document is daily updated. Last update 17-03-2020

## Authors
Amanda Fernández Fontelo, Alejandra Cabaña, Argimiro Arratia, David Moriña and Pere Puig
