# uboone_noise
`rmsNoiseFreq.C` is online analysis code. It computes the ADC pedestal, metrics for it, 
and average FFT magnitude plots for given event processed by wire-cell display. It takes
a file name as input and creates a text file with pedestal metrics and a root file with 
the RMS of each channel with some identifiers. 

`getDeadChannels.C` prints the number of dead channels for each event listed in `eventsFileName`,
which is a text file with a list of all events that are to be processed. Easily make such a file
with `ls [files] >> [eventsFileName]`

`treeMaker.C` contains several functions. `makeFFTMetrics` is self explainatory and takes as 
input a text file containing root files for events to be processed and the number of events 
listed in the text file.`chirpID` takes an input file with a list of root files for events 
to be processed. It adds a brach `chirping_flag` to each root file in the input and marks it 
as true for a channel if that channel is dead but alive in another event. `addTime` is a 
function that has as yet been unused, but it extracts time information from celltree files 
and adds it to `Trun` in root files from analysis code. `makeMetricTree` compiles text files
containing ADC pedestal metrics into a root file. 

`plotter.C` contains methods for plotting text metrics (after using `makeMetricTree`), number
of dead channels, and FFT metrics. 

`plotRMS.C` plots the ADC pedestal and FFT spectra for each event root file listed in an input
text file. 
