# InpaintingRevisited
Supplemental material to the paper Audio Inpainting: Revisited and Reweighted.
by Ondřej Mokrý and Pavel Rajmic. See a short guide to the files contained
in the repository provided below. In case of any questions, feel free to
contact the authors.

This repository contains the source for the web page https://ondrejmokry.github.io/InpaintingRevisited/.

Due to space limitations, there are only selected examples
of the results in the paper. To provide a comprehensive
comparison of the method, supplemental results are presented in the folder `figures`.
All the present figures are generated with scripts and data
available in the folder `matlab/data`. 

**Subfolders of the matlab folder:**
  * `data` -- The folder contains results of all the performed experiments,
              together with scripts that plot the used figures; both data
              (`*.mat` files) and the scripts (`plot_*.m`) are divided into
              subfolders based on the respective method.
  * `Janssen` -- Functions for audio inpainting using AR-based method by
              Janssen et al. [1].
  * `OMP` -- Functions for audio inpainting using the OMP algorithm [2, 3].
  * `PemoQ` -- Functions for evaluation of the perceived quality of restored
              audio [4]. This folder is not part of the package due to copyright
              issues. However, it needs to be included to reproduce the presented
              results.
  * `reweighted l1 relaxation` -- Implementation of the l1-relaxation-based
              methods presented in the paper and supporting functions.
  * `SPAIN` -- Functions for audio inpainting using the SPAIN algorithm [5].
  * `wavs` -- This folder is created when the script `demo.m` is executed.
              The script then saves the restored signals as WAV files to
              this folder, sorted by signals.

**Separate files:**
  * `demo.m` -- The user is asked to choose signal, gap parameters (number,
              length) and algorithms from the preset options; the script
              then performs audio with such settings (if necessary, the
              parameters of the methods need to be changed in the m-file).
  * `EBU_SQAM.mat` -- Collection of signals used for the evaluation [6].
  * `global_test_...` -- Implementation of all the tests performed for the
              experimental part of the paper, except for the overall evaluation.
              These 4 scripts are the source codes which generated the files in
              the subfolder `data`.
  * `inpainting_comparison_full.m`,  `inpainting_comparison_half.m
              -- Implementation of the methods used for overall evaluation;
              the scripta reproduce all the data used in the paper, thus it
              takes a substantial amount of time to run them.

All the scripts were tested in MATLAB R2019a and using LTFAT [7, 8]
(version 2.3.1), which is not contained in this package.

For the plotting scripts to be functional, MATLAB version R2018b or higher
is necessary, because the function `sgtitle()` is used multiple times.

Furthermore, to work correctly, the interpreter for legend and titles
should be set to latex, use for example

  `set(groot, 'defaultLegendInterpreter','latex');`
  
  `set(groot, 'defaultTextInterpreter','latex');`

**References:**

[1] A. J. E. M. Janssen, R. N. J. Veldhuis, and L. B. Vries, “Adaptive
    interpolation of discrete-time signals that can be modeled as
    autoregressive processes,” IEEE Trans. Acoustics, Speech and Signal
    Processing, 1986.

[2] A. Adler, V. Emiya, M. Jafari, M. Elad, R. Gribonval, and M. Plumbley,
    “Audio Inpainting,” IEEE Trans. Audio, Speech, and Language Processing,
    2012.

[3] T. Blumensath, “Sparsify toolbox,” online, URL:
    https://www.southampton.ac.uk/engineering/about/staff/tb1m08.page#software

[4] R. Huber and B. Kollmeier, “PEMO-Q—A new method for objective
    audio quality assessment using a model of auditory perception,” IEEE
    Trans. Audio Speech Language Processing, 2006.

[5] O. Mokry, P. Zaviska, P. Rajmic, and V. Vesely, “Introducing SPAIN
    (SParse Audio INpainter),” in 2019 27th European Signal Processing
    IEEE/ACM TRANSACTIONS ON AUDIO, SPEECH, AND LANGUAGE PROCESSING 20
    Conference (EUSIPCO). IEEE, 2019.

[6] EBU SQAM CD: Sound quality assessment material recordings for
    subjective tests. online. URL: https://tech.ebu.ch/publications/sqamcd.
    [Online]. Available: https://tech.ebu.ch/publications/sqamcd

[7] Z. Prusa, P. L. Sondergaard, N. Holighaus, C. Wiesmeyr, and P. Balazs,
    “The Large Time-Frequency Analysis Toolbox 2.0,” in Sound, Music, and
    Motion. Springer International Publishing, 2014.

[8] P. L. Sondergaard. (2013) LTFAT webpage. URL:
    http://ltfat.sourceforge.net
