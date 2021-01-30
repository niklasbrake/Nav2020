## What is this?
In this repository is the code associated with the following paper: 

**Brake N, Mancino M, Yan Y, Shimomura S, Kubo Y, Bowie D, and Khadra A (2021). "Intrinsic Gating Behavior of Voltage-Gated Sodium Channels Predetermines Regulation by Auxiliary β-subunits."** *****Submitted.*****

There are four sections to this repository:
* **Modelling**. This contains the functions used to fit the parameters of a Markov model of a voltage-gated ion channel to electrophysiology data.
* **Data Processing**. This contains the functions used to extract data from the .abf files produced from the pClamp software. These functions extract relevant features for fitting.
* **GUI**. This is a user interface that I used for facilitating structural changes to the Markov model. With it, you can draw states and connections between states, specify the functions dictating the rate constants, and automatically generates an \*.m file that encodes the system of ODEs described by the model.
* **Figures**. These are the functions that were used to generate the figures from the paper. Data for the representative traces are also included. Note these data files are about 300 MB in total. Summary data files were generated using the functions in the Data Processing folder. See the methods section of the manuscript for more details.

## Acknowledgements
I completed the work for this project in 2018-2019 during my PhD under the supervision of [Dr. Anmar Khadra](http://www.medicine.mcgill.ca/physio/khadralab/) and in collaboration with [Dr. Derek Bowie](http://www.medicine.mcgill.ca/pharma/dbowielab/) at McGill Univeristy. The electrophysiology data used for fitting was collected by Adamo Mancino, who at the time was a master's student in the Bowie lab. Figure 3 of the paper presents voltage-clamp fluormetry data. This data was collected by Adamo Mancino and Yuhao Yan while visting the lab of Dr. Yoshihiro Kubo at the National Institute for Physiological Sciences in Japan.

## License
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

To reference this code, please cite the journal article mentioned above.
