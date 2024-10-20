This is the README file for the FastEEC project available at 
    
  https://github.com/abudhraj/FastEEC

which provides a fast evaluation of higher-point energy correlators for both integer 
as well as non-integer values of N. Our code uses the fastjet package for jet 
reclustering and resolving branches, exploiting that for correlations at a given angle, 
radiation whose separation is much smaller can be treated as one. For integer N values, 
our default option is to cluster with Cambridge/Aachen and use a fixed resolution scale f, 
but we also found that the kT algorithm with resolution scale f = k_T,min^2 performs 
particularly well. These are provided in the files: `eec_fast.cc` (for C/A clustering) 
and `eec_fast_kt.cc` (for k_T clustering). For the case of non-integer N values, we use the
Cambridge/Aachen to re-cluster the jet but now with a fixed number of subjets `nsub` in order
to optimize the working efficiency of our code. *Note that for non-integer N, this is sometimes 
referred to as ν-point energy correlator or ν-correlator (in short).* The code is provided in the 
file: `eec_fast_nu_point.cc`. We now provide a new implementation of the N-point projected energy 
correlators applicable to any (integer or non-integer) N values that we recently proposed. This new
implementation achieves a computational time of M^2 log(M), for any N value. The code is provided
in the file: `eec_angles.cc`. 

Additionally, we provide the option to use higher powers of transverse momenta to 
compute the projected energy correlator for integer N. Once again, two flavors of jet 
clustering, i.e., Cambridge/Aachen and k_T are provided. In this case, as the higher power 
reduces the soft sensitivity, a reasonable description for weight kappa = 2 can even 
be obtained by using resolution scale as small as f=8, corresponding to a relative 
error of about 0.5-1.5% for various values of N. We achieve a substantial gain in 
computation time, up to as high as four orders of magnitude, for N=7. The files that 
implement higher powers of transverse momentum weights are: `eec_fast_weight.cc` 
(for C/A) and `eec_fast_kt_weight.cc` (for k_T).

**version 0.3: includes the implementation of projected N-point energy correlators for any (integer or non-integer) N, using the new proposed parameterization**

If you use our package in research, please include a citation to

  S. Alipour-fard, A. Budhraja, J. Thaler and W. Waalewijn, New Angles on Energy Correlators, arXiv: 2410.XXXX,

  A. Budhraja, H. Chen and W. Waalewijn, ν-point energy correletors with FastEEC: 
  small x physics from LHC jets, arXiv: 2409.12235.

and
  
  A. Budhraja, W. Waalewijn, Fast Evaluation of N-point Energy Correlators, 
  arXiv: 2406.08577. 

# Instructions for users

You can download our code freely from github. Unpack it, enter the resulting 
directory and run

```
  make
```

As the code uses FastJet libraries, make sure the directory paths are predefined 
before executing make. Alternatively, one can also compile a specific eec_fast.cc by 

```
  make eec_fast
```

Note that, if you have compiled FastJet with a C++11 compiler, you need to make the 
appropriate change in the Makefile so it also uses a C++11 compiler (if your 
compiler needs options to get that behaviour, e.g. -std=c++11, then include those 
same options in the CXXFLAGS).

When the code is executed from the command line the following inputs are needed

```
./eec_fast input_file events N f minbin nbins output_file
```
The above command line parameters require
> The input_file from which events should be read <br/>
> events > 0 is the number of events <br/>
> 1 < N < 9 specifies which point correlator to compute <br/>
> resolution factor f > 1 <br/>
> minbin < 0 specifies the minimum bin value as log10(R_L) needed for sampling output <br/>
> number of bins > 1 <br/>
> The output_file in which data should be stored <br/>

This message is also displayed any time incorrect number of command line inputs are 
entered by the user. The output is collected in the user-specified number of bins. 
This output is normalized such that the sum of entries of all the bins equals 1 for 
each event, due to momentum conservation.

For the higher power energy correlators, the user can also provide the value of 
kappa by which the transverse momenta are weighted as a command line input as

```
./eec_fast_weight input_file events N f kappa minbin nbins output_file
```
with
> weight kappa > 0 <br/>

In this case the output is no longer normalized because the observable is not 
collinear safe.

The same number of command line inputs are also used for the k_T versions of the 
code both with weight = 1 (`eec_fast_kt.cc`) and weight &ne; 1 
(`eec_fast_kt_weight.cc`) except for the difference that 
> resolution factor f' > 0 <br/>

is now allowed. Note that in this case the resolution factor is f = f' k_T,min^2, 
where the command-line parameter is f'.

For the code for non-integer ν-point correlators, the command line input replaces 
the jet resolution parameter with the number of subjets as
```
./eec_fast_nu_point input_file events N nsub minbin nbins output_file
```
with
> 17 > number of subjets (nsub) > 1 <br/>

The same code can be used to compute the ν-point energy correlator for integer values of
ν as well, though to reach the same level of accuracy it is slower.

For the new implementation of the N-point correlators, we do not use any approximations, 
hence the command line inputs only require the following 
```
./eec_angles input_file events N minbin nbins output_file
```
The same code can be used for any N and provides a drastic computational speed up, without
any approximations.

# Instructions for the input_file

The input_file consists of a line for each particle, with four numbers in each line.
These four numbers are the `jet_number, pT, eta and phi` of the particle, and must 
be provided in this order.

We include a sample input file generated using publically available "MIT Open Data",
which utilizes the reprocessed data on jets from the CMS 2011A Open Data. This 
data.txt file consists of 100 000 jets with transverse momenta 500 < pT < 550 GeV 
and rapidity |eta| < 1.9. 

# Instructions to read the output_file
 
The first line of the output file consists of events, nbins, minbin, maxbin, where 
maxbin = log10(1) = 0. The next line contains the histogram values, and the number 
of entries equals nbins. 

We have provided a set of output files for the specific case of N=4 point correlator 
with different files corresponding to: 
- `test_4pt.out`: that contains the result generated by the full method (using the 
available energy correlators package). 

The other files correspond to the different approximations we propose with: 
- `test_4pt_8.out`, `test_4pt_32.out` and `test_4pt_64.out` specifying the use 
of use constant resolution factors of 8, 32 and 64 with C/A clustering.
- The version with k_T clustering and f=k_T,min^2 is called `test_4pt_1_kt.out`.

There is also a small Mathematica notebook enclosed with this release called 
`eec_analysis.nb`. The notebook illustrates explicitly how the output files can be 
read. We also provide instructions to reproduce the plot shown in the paper for N=4 
case corresponding to both the normalized distribution as well as the relative error 
of the different approximations when compared to full.

# Advanced guide for specific implementations

The code assumes that the input consists of a single jet and that reclustering it 
with a jet radius R = 1.5 returns exactly one jet. The jet radius in the code does 
not need to match that used to obtain the jets used as input, in fact it should be 
sufficiently larger than it. If needed, the jet radius in the code can be modified 
by changing the following line in `eec_compute.h` file (for weight =1) or 
`eec_higher_weight.h` (for weight &ne; 1) or `eec_nu_point.h` (for ν-point correlators)

```
const double R = 1.5;
```
Indeed we have applied it to entire events for e+e- collisions, which also requires 
modifying the distance measure, clustering the whole event into one jet.

Although any number of bins can be specified, the code internally fixes a large
dimension in which the histogrammed data is stored. This value is set to 1000 in the 
current implementation, corresponding to a maximum of 1000 bins. For users who would 
want to utilize more bins, the dimension specified in the `eec_compute.h` file (for 
weight =1) or `eec_higher_weight.h` (for weight &ne; 1) or `eec_nu_point.h` (for ν-point 
correlators) must be changed in the line

```
std::vector<double> res(1000);
```
to the desired value. The rest of the code remains unchanged. 

Similarly, if one is interested in using wider jets with a very high resolution, 
the number of subjets will increase. The code constructs a two-dimensional array to 
compute the distance between all subjets and using an array of 200 X 200 that may 
need to be increased when there are very many subjets. This can be done by modifying 
the line of the code below

```
//Precalculate all bin positions
double bindist[200][200];
```
contained in the header files `eec_compute.h` file (for weight =1) or 
`eec_higher_weight.h` (for weight &ne; 1).
